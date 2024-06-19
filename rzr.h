#ifndef RZR_H

/*
TODO:
 - allow out-of-memory errors without crashing?
FIXME:
 - broken demo cases: difference between rounded boxes is broken, pattern is
   not working
 - A circle precalcs/stores a number of values equal to its subpixel radius
   regardless of image size.
*/

#include <assert.h>
#include <math.h>

struct rzr_tx {
	float basis0_x, basis0_y;
	// basis1_x=-basis0_y, basis1_y = basis0_x
	float origin_x, origin_y;
};

static inline float rzr_tx_get_scale(struct rzr_tx* tx)
{
	const float dx=tx->basis0_x, dy=tx->basis0_y;
	return sqrtf(dx*dx + dy*dy);
}

enum rzr_op_code {
	RZROP_PICK = 0,
	RZROP_SWAP,
	RZROP_DROP,

	// these might be useful?
	//RZROP_ZERO,
	//RZROP_ONE,

	RZROP_POLY,
	RZROP_VERTEX,

	RZROP_CIRCLE,

	RZROP_UNION,
	RZROP_INTERSECTION,
	RZROP_DIFFERENCE,

	RZROP_END,
};

struct rzr_op {
	enum rzr_op_code code;
	union {
		struct { int stack_index;  } pick;
		struct { int cx,cy,radius; } circle;
		struct { int n_vertices;   } poly;
		struct { int x,y;          } vertex;
	};
};

struct rzr {
	int width, height;
	int supersampling_factor;

	int tx_stack_cap;
	int tx_stack_height;
	struct rzr_tx* tx_stack;

	int prg_cap;
	int prg_length;
	struct rzr_op* prg;

	int in_poly, poly_op_index;

	// (FYI) use these to see how much memory rzr is using
	int tx_stack_height_max;  // maximum seen value of rzr.tx_stack_height
	int prg_length_max;       // maximum seen value of rzr.prg_length
	size_t scratch_alloc_max; // maximum seen allocation cursor in scratch allocations (see rzr_render())

};

#define RZR_MAX_SUPERSAMPLING_FACTOR (32)
#ifndef RZR_DEFAULT_SUPERSAMPLING_FACTOR
#define RZR_DEFAULT_SUPERSAMPLING_FACTOR (8)
#endif

static inline int rzr_get_supersampling_factor(struct rzr* rzr)
{
	return
		  rzr->supersampling_factor <= 0                           ? RZR_DEFAULT_SUPERSAMPLING_FACTOR
		: rzr->supersampling_factor > RZR_MAX_SUPERSAMPLING_FACTOR ? RZR_MAX_SUPERSAMPLING_FACTOR
		: rzr->supersampling_factor;
}

static inline void rzr_init(struct rzr* rzr, int tx_stack_cap, struct rzr_tx* tx_stack, int prg_cap, struct rzr_op* prg, int width, int height, float pixels_per_unit, int supersampling_factor)
{
	memset(rzr, 0, sizeof *rzr);
	rzr->tx_stack_cap = tx_stack_cap;
	rzr->tx_stack = tx_stack;
	rzr->prg_cap = prg_cap;
	rzr->prg = prg;
	rzr->supersampling_factor = supersampling_factor;

	const int ssf = rzr_get_supersampling_factor(rzr);
	const float subpixels_per_unit = pixels_per_unit * (float)ssf;

	rzr->tx_stack_height = 1;
	rzr->tx_stack_height_max = rzr->tx_stack_height;
	struct rzr_tx* tx0 = &rzr->tx_stack[rzr->tx_stack_height-1];
	tx0->basis0_x = subpixels_per_unit;
	tx0->basis0_y = 0.0f;
	tx0->origin_x = (float)(width*ssf)  * 0.5f;
	tx0->origin_y = (float)(height*ssf) * 0.5f;

	rzr->width = width;
	rzr->height = height;
}

static inline struct rzr_tx* rzr_get_current_tx(struct rzr* rzr)
{
	assert(rzr->tx_stack_height > 0);
	return &rzr->tx_stack[rzr->tx_stack_height-1];
}

static inline void rzr_tx_save(struct rzr* rzr)
{
	assert(rzr->tx_stack_height < rzr->tx_stack_cap);
	rzr->tx_stack[rzr->tx_stack_height++] = *rzr_get_current_tx(rzr);
	if (rzr->tx_stack_height > rzr->tx_stack_height_max) {
		rzr->tx_stack_height_max = rzr->tx_stack_height;
	}
}

static inline void rzr_tx_restore(struct rzr* rzr)
{
	assert(rzr->tx_stack_height > 1);
	rzr->tx_stack_height--;
}

static inline void rzr_tx_translate(struct rzr* rzr, float x, float y)
{
	struct rzr_tx* tx = rzr_get_current_tx(rzr);
	tx->origin_x += (x*tx->basis0_x - y*tx->basis0_y);
	tx->origin_y += (x*tx->basis0_y + y*tx->basis0_x);
}

static inline float rzr__deg2rad(float deg) { return deg * 0.017453292519943295f; }

static inline void rzr_tx_rotate(struct rzr* rzr, float degrees)
{
	struct rzr_tx* tx = rzr_get_current_tx(rzr);
	const float radians = rzr__deg2rad(degrees);
	const float u = cosf(radians);
	const float v = sinf(radians);
	const float bx = u*tx->basis0_x - v*tx->basis0_y;
	const float by = u*tx->basis0_y + v*tx->basis0_x;
	tx->basis0_x = bx;
	tx->basis0_y = by;
}

static inline void rzr_tx_scale(struct rzr* rzr, float scalar)
{
	struct rzr_tx* tx = rzr_get_current_tx(rzr);
	tx->basis0_x *= scalar;
	tx->basis0_y *= scalar;
}

static inline struct rzr_op* rzr_op(struct rzr* rzr, enum rzr_op_code code)
{
	assert(rzr->prg_length < rzr->prg_cap);
	struct rzr_op* op = &rzr->prg[rzr->prg_length++];
	if (rzr->prg_length > rzr->prg_length_max) rzr->prg_length_max = rzr->prg_length;
	memset(op, 0, sizeof *op);
	assert(0 <= code && code < RZROP_END);
	assert(!rzr->in_poly || code == RZROP_VERTEX);
	op->code = code;
	return op;
}

static inline int rzr_float_to_int(float f)
{
	return (int)roundf(f);
}

static inline void rzr_pick(struct rzr* rzr, int stack_index)
{
	struct rzr_op* op = rzr_op(rzr, RZROP_PICK);
	op->pick.stack_index = stack_index;
}

static inline void rzr_dup(struct rzr* rzr) { rzr_pick(rzr, -1); }
static inline void rzr_swap(struct rzr* rzr) { rzr_op(rzr, RZROP_SWAP); }
static inline void rzr_drop(struct rzr* rzr) { rzr_op(rzr, RZROP_DROP); }

static inline void rzr_union(struct rzr* rzr)        { rzr_op(rzr, RZROP_UNION); }
static inline void rzr_intersection(struct rzr* rzr) { rzr_op(rzr, RZROP_INTERSECTION); }
static inline void rzr_difference(struct rzr* rzr)   { rzr_op(rzr, RZROP_DIFFERENCE); }

static inline void rzr_circle(struct rzr* rzr, float radius)
{
	struct rzr_op* op = rzr_op(rzr, RZROP_CIRCLE);
	struct rzr_tx* tx = rzr_get_current_tx(rzr);
	op->circle.radius = rzr_float_to_int(radius * rzr_tx_get_scale(tx));
	op->circle.cx = rzr_float_to_int(tx->origin_x);
	op->circle.cy = rzr_float_to_int(tx->origin_y);
}

static inline void rzr_begin_poly(struct rzr* rzr)
{
	assert(!rzr->in_poly);
	rzr->poly_op_index = rzr->prg_length;
	struct rzr_op* op = rzr_op(rzr, RZROP_POLY);
	op->poly.n_vertices = -1;
	rzr->in_poly = 1;
}

static inline void rzr_end_poly(struct rzr* rzr)
{
	assert(rzr->in_poly);
	const int n_vertices = rzr->prg_length - (rzr->poly_op_index + 1);
	assert((n_vertices >= 3) && "a polygon must have at least 3 vertices");
	rzr->prg[rzr->poly_op_index].poly.n_vertices = n_vertices;
	rzr->in_poly = 0;
}

static inline void rzr_vertex(struct rzr* rzr, float x, float y)
{
	assert(rzr->in_poly);
	struct rzr_op* op = rzr_op(rzr, RZROP_VERTEX);
	struct rzr_tx* tx = rzr_get_current_tx(rzr);
	op->vertex.x = rzr_float_to_int(tx->origin_x + x*tx->basis0_x - y*tx->basis0_y);
	op->vertex.y = rzr_float_to_int(tx->origin_y + x*tx->basis0_y + y*tx->basis0_x);
}

static inline void rzr_star(struct rzr* rzr, int n, float outer_radius, float inner_radius)
{
	const int nn = n+n;
	rzr_begin_poly(rzr);
	float phi = 0.0f;
	const float inc = 6.283185307179586f / (float)nn;
	for (int i = 0; i < nn; i++) {
		const float d = i&1 ? outer_radius : inner_radius;
		rzr_vertex(rzr, d*cosf(phi), d*sinf(phi));
		phi += inc;
	}
	rzr_end_poly(rzr);
}

static inline void rzr_pattern(struct rzr* rzr, float* ws)
{
	int n=0;
	float ww = 0.0f;
	float* p = ws;
	while (*p > 0) { ww+=*p ; p++; n++; }
	assert(((n&1) == 0) && "ws count must be even");
	assert((n >= 2) && "empty pattern");
	assert(ww > 0.0f);

	float x = 0.0f;
	p = ws;
	int ns = 0;
	const float y0 = -1.0f;
	const float y1 =  1.0f;
	for (int i = 0; i < n; i++) {
		const float w = *(p++);
		if ((i&1) == 0) {
			rzr_begin_poly(rzr);
			rzr_vertex(rzr, x   , y0);
			rzr_vertex(rzr, x+w , y0);
			rzr_vertex(rzr, x+w , y1);
			rzr_vertex(rzr, x   , y1);
			rzr_end_poly(rzr);
			if (ns++) rzr_union(rzr);
		}
		x += w;
	}
}

static inline void rzr_line(struct rzr* rzr, float width)
{
	assert(!"TODO");
}

static inline void rzr_split(struct rzr* rzr)
{
	assert(!"TODO");
}

static inline void rzr_box(struct rzr* rzr, float width, float height)
{
	rzr_begin_poly(rzr);
	rzr_vertex(rzr, -width, -height);
	rzr_vertex(rzr,  width, -height);
	rzr_vertex(rzr,  width,  height);
	rzr_vertex(rzr, -width,  height);
	rzr_end_poly(rzr);
}

static inline void rzr_rounded_box(struct rzr* rzr, float width, float height, float radius)
{
	if (radius <= 0.0f) {
		rzr_box(rzr, width, height);
		return;
	}
	rzr_box(rzr, width, height);
	const float rh = radius*0.5f;
	for (int pass = 0; pass < 2; pass++) {
		const float d = (pass==0) ? rh : (pass==1) ? radius : 0.0f;
		for (int corner = 0; corner < 4; corner++) {
			rzr_tx_save(rzr);
			switch (corner) {
			case 0: rzr_tx_translate(rzr, -width+d, -height+d); break;
			case 1: rzr_tx_translate(rzr,  width-d, -height+d); break;
			case 2: rzr_tx_translate(rzr,  width-d,  height-d); break;
			case 3: rzr_tx_translate(rzr, -width+d,  height-d); break;
			default: assert(!"err");
			}
			if (pass == 0) {
				rzr_box(rzr, rh, rh);
				rzr_difference(rzr);
			} else if (pass == 1) {
				rzr_circle(rzr, radius);
				rzr_union(rzr);
			} else {
				assert(!"unreachable");
			}
			rzr_tx_restore(rzr);
		}
	}
}

static inline void rzr_arc(struct rzr* rzr, float aperture_degrees, float radius, float width)
{
	assert(!"TODO");
}

static inline void rzr_segment(struct rzr* rzr, float x0, float y0, float x1, float y1, float r)
{
	assert(!"TODO");
}

static inline void rzr_isosceles_triangle(struct rzr* rzr, float w, float h)
{
	rzr_begin_poly(rzr);
	rzr_vertex(rzr,  0, 0);
	rzr_vertex(rzr,  w, h);
	rzr_vertex(rzr, -w, h);
	rzr_end_poly(rzr);
}

static inline void rzr_isosceles_trapezoid(struct rzr* rzr, float r1, float r2, float h)
{
	if (r1 <= 0.0f) {
		rzr_isosceles_triangle(rzr, r2, h);
		return;
	}
	rzr_begin_poly(rzr);
	rzr_vertex(rzr, -r1, 0);
	rzr_vertex(rzr,  r1, 0);
	rzr_vertex(rzr,  r2, h);
	rzr_vertex(rzr, -r2, h);
	rzr_end_poly(rzr);
}


void rzr_render(struct rzr*, size_t scratch_cap, void* scratch, int stride, uint8_t* pixels);

#ifndef RZR_NO_SHORT_NAMES

#ifndef RZR_INSTANCE
#define RZR_INSTANCE (rzr)
#endif

#define Save()                       rzr_tx_save(RZR_INSTANCE)
#define Restore()                    rzr_tx_restore(RZR_INSTANCE)
#define Translate(x,y)               rzr_tx_translate(RZR_INSTANCE,x,y)
#define Rotate(degrees)              rzr_tx_rotate(RZR_INSTANCE,degrees)
#define Scale(scalar)                rzr_tx_scale(RZR_INSTANCE,scalar)

#define Dup()                        rzr_dup(RZR_INSTANCE)
#define Pick(i)                      rzr_pick(RZR_INSTANCE,i)
#define Swap()                       rzr_swap(RZR_INSTANCE)
#define Drop()                       rzr_drop(RZR_INSTANCE)

#define Union()                      rzr_union(RZR_INSTANCE)
#define Intersection()               rzr_intersection(RZR_INSTANCE)
#define Difference()                 rzr_difference(RZR_INSTANCE)

#define Circle(r)                    rzr_circle(RZR_INSTANCE,r)

#define BeginPoly()                  rzr_begin_poly(RZR_INSTANCE)
#define Vertex(x,y)                  rzr_vertex(RZR_INSTANCE,x,y)
#define EndPoly()                    rzr_end_poly(RZR_INSTANCE)

#define Star(n,o,i)                  rzr_star(RZR_INSTANCE,n,o,i)
#define Pattern(...)                 rzr_pattern(RZR_INSTANCE,(float[]) { __VA_ARGS__, 0 })
#define Line(w)                      rzr_line(RZR_INSTANCE,w)
#define Split(w)                     rzr_split(RZR_INSTANCE)
#define Arc(a,r,w)                   rzr_arc(RZR_INSTANCE,a,r,w)
#define Box(w,h)                     rzr_box(RZR_INSTANCE,w,h)
#define RoundedBox(w,h,r)            rzr_rounded_box(RZR_INSTANCE,w,h,r)
#define Segment(x0,y0,x1,y1,r)       rzr_segment(RZR_INSTANCE,x0,y0,x1,y1,r)
#define IsoscelesTriangle(w,h)       rzr_isosceles_triangle(RZR_INSTANCE,w,h)
#define IsoscelesTrapezoid(r1,r2,h)  rzr_isosceles_trapezoid(RZR_INSTANCE,r1,r2,h)

//#define Simplify(n,e)                              rzr_simplify(RZR_INSTANCE,n,e)
//#define CubicBezier(x0,y0,c0x,c0y,c1x,c1y,x1,y1)   rzr_cubic_bezier(RZR_INSTANCE,x0,y0,c0x,c0y,c1x,c1y,x1,y1)
// TODO I think curve and simplification should be something separate, maybe
// combined... but something that emits the proper Vertex() calls

#endif

#define RZR_H
#endif

/*
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2024 Anders Kaare Straadt
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
*/

#ifndef RZR_H

/*
TODO/FIXME:
 - Allow out-of-memory errors without crashing?
 - Improve robustness of rzr__xline()? It powers "Split", "Line" and "Pattern".
   it probably has numerical problems with near-vertical lines (prevcx~=cx),
   and lines crossing close to corners
 - A circle precalcs/stores a number of values equal to its subpixel radius
   regardless of image size.
 - `post_order_stack` and `xops` could use proper memory management; see
   comments.
 - The polygon xspan sorting/finding is fragile; see comment; I've also seen
   `assert(xe->side == (i2&1));` fail in the wild.
 - rzr_rounded_box() doesn't handle when radius is larger than box.
*/

#include <assert.h>
#include <stdint.h>
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
	RZROP_PICK = 1,
	RZROP_SWAP,
	RZROP_DROP,

	RZROP_ZERO,
	RZROP_ONE,

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

#define RZR_MXN_MAX_SIZES (10)
struct rzr {
	int bitmap_width,  bitmap_height;
	int virtual_width, virtual_height;
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

	int n_widths;
	int widths[RZR_MXN_MAX_SIZES];
	int n_heights;
	int heights[RZR_MXN_MAX_SIZES];

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

static inline void rzr__init_common(struct rzr* rzr, int tx_stack_cap, struct rzr_tx* tx_stack, int prg_cap, struct rzr_op* prg, float pixels_per_unit, int supersampling_factor, int virtual_width, int virtual_height, int bitmap_width, int bitmap_height)
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
	tx0->origin_x = (float)(virtual_width*ssf)  * 0.5f;
	tx0->origin_y = (float)(virtual_height*ssf) * 0.5f;

	rzr->bitmap_width = bitmap_width;
	rzr->bitmap_height = bitmap_height;
	rzr->virtual_width = virtual_width;
	rzr->virtual_height = virtual_height;

	rzr->n_widths = 1;
	rzr->n_heights = 1;
	rzr->widths[0]  = virtual_width;;
	rzr->heights[0] = virtual_height;
}

static inline void rzr_init(struct rzr* rzr, int tx_stack_cap, struct rzr_tx* tx_stack, int prg_cap, struct rzr_op* prg, int width, int height, float pixels_per_unit, int supersampling_factor)
{
	rzr__init_common(rzr, tx_stack_cap, tx_stack, prg_cap, prg, pixels_per_unit, supersampling_factor, width, height, width, height);
}

// rzr_init_MxN() configures rzr to render a M x N grid. `widths` and `heights`
// are arrays defining a "size matrix". Both arrays are 0-terminated, and
// negative values means "stretch".  Stretch columns/rows are rendered 1 pixel
// and 0 units wide/high.
// Example: widths  = [10, -1, 20, 0] (or: [10, "stretch", 20, "end"])
//          heights = [15, -1, 25, 0] (or: [15, "stretch", 25, "end"])
//          pixels_per_units = 10
//    +---+-+---+
//    |a  |b|  c|
//    |   | |   |
//    |   | |   |
//    +---+-+---+
//    |d  |e|  f|
//    +---+-+---+
//    |   | |   |
//    |   | |   |
//    |g  |h|  i|
//    +---+-+---+
//     pixel pos    pixel dim     units pos    units dim
// a     0 , 0       10 x 15    -1.5 , -2.0    1.0 x 1.5
// b    10 , 0        1 x 15     0.0 , -2.0      0 x 1.5
// c    11 , 0       20 x 15     0.0 , -2.0    2.0 x 1.5
// d     0 , 15      10 x 1     -1.5 ,  0.0    1.0 x 0
// e    10 , 15       1 x 1      0.0 ,  0.0      0 x 0
// f    11 , 15      20 x 1      0.0 ,  0.0    2.0 x 0
// g     0 , 16      10 x 25    -1.5 ,  0.0    1.0 x 2.5
// h    10 , 16       1 x 25     0.0 ,  0.0      0 x 2.5
// i    11 , 16      20 x 25     0.0 ,  0.0    2.0 x 2.5
// Total pixel width   = 10+1+20 = 31     unit width  = 1.0+2.0 = 3.0
// Total pixel heights = 15+1+25 = 41     unit height = 1.5+2.5 = 4.0
static inline void rzr_init_MxN(struct rzr* rzr, int tx_stack_cap, struct rzr_tx* tx_stack, int prg_cap, struct rzr_op* prg, const int* widths, const int* heights, float pixels_per_unit, int supersampling_factor)
{
	int virtual_width=0, virtual_height=0;
	int bitmap_width=0, bitmap_height=0;
	for (int axis = 0; axis < 2; axis++) {
		const int* sizes = (axis==0) ? widths : (axis==1) ? heights : NULL;
		int n=0;
		int bitmap_size=0, virtual_size=0;
		int n_consec_stretch = 0;
		for (const int* p = sizes; *p != 0; n++, p++) {
			const int size = *p;
			assert(size != 0);
			if (size < 0) {
				n_consec_stretch++;
				assert((n_consec_stretch < 2) && "consecutive ''stretch'' sizes are non-sensical and not supported");
				bitmap_size += 1;
			} else {
				n_consec_stretch = 0;
				bitmap_size += size;
				virtual_size += size;
			}
		}
		assert(0 < n && n <= RZR_MXN_MAX_SIZES);
		if (axis == 0) {
			bitmap_width=bitmap_size;
			virtual_width=virtual_size;
		} else if (axis==1) {
			bitmap_height=bitmap_size;
			virtual_height=virtual_size;
		} else {
			assert(!"unreachable");
		}
	}
	rzr__init_common(rzr, tx_stack_cap, tx_stack, prg_cap, prg, pixels_per_unit, supersampling_factor, virtual_width, virtual_height, bitmap_width, bitmap_height);
	for (int axis = 0; axis < 2; axis++) {
		const int* sizes = (axis==0) ? widths : (axis==1) ? heights : NULL;
		int* wp    = (axis==0) ? rzr->widths : (axis==1) ? rzr->heights : NULL;
		const int* p = sizes;
		while (*p != 0) *(wp++) = *(p++);
		const int n = p-sizes;
		if (axis == 0) {
			rzr->n_widths = n;
		} else if (axis == 1) {
			rzr->n_heights = n;
		} else {
			assert(!"unreachable");
		}
	}
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

static inline void rzr_invert_tx(struct rzr_tx* tx)
{
	const float I = 1.0f / (tx->basis0_x*tx->basis0_x + tx->basis0_y*tx->basis0_y);
	const float bx =  I * tx->basis0_x;
	const float by = -I * tx->basis0_y;
	const float ox = I * ((-tx->basis0_y * tx->origin_y) - (tx->basis0_x * tx->origin_x));
	const float oy = I * (( tx->basis0_y * tx->origin_x) - (tx->basis0_x * tx->origin_y));
	tx->basis0_x = bx;
	tx->basis0_y = by;
	tx->origin_x = ox;
	tx->origin_y = oy;
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

static inline void rzr_zero(struct rzr* rzr)
{
	rzr_op(rzr, RZROP_ZERO);
}

static inline void rzr_one(struct rzr* rzr)
{
	rzr_op(rzr, RZROP_ONE);
}

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
	struct rzr_op* op = &rzr->prg[rzr->poly_op_index];
	op->poly.n_vertices = n_vertices;

	int iprev = n_vertices-1;
	int64_t area = 0;
	for (int i = 0; i < n_vertices; i++) {
		struct rzr_op* op = &rzr->prg[rzr->poly_op_index+1+i];
		struct rzr_op* prevop = &rzr->prg[rzr->poly_op_index+1+iprev];
		assert(op->code == RZROP_VERTEX);
		assert(prevop->code == RZROP_VERTEX);
		iprev = i;
		area += (op->vertex.x - prevop->vertex.x) * (op->vertex.y + prevop->vertex.y);
	}
	rzr->in_poly = 0;
	if (area >= 0) {
		// invalid polygon; replace with zero
		rzr->prg_length = rzr->poly_op_index;
		rzr_zero(rzr);
	}
}

static inline void rzr_map_point(struct rzr_tx* tx, float x, float y, float* out_u, float* out_v)
{
	if (out_u) *out_u = tx->origin_x + x*tx->basis0_x - y*tx->basis0_y;
	if (out_v) *out_v = tx->origin_y + x*tx->basis0_y + y*tx->basis0_x;
}

static inline void rzr_vertex(struct rzr* rzr, float x, float y)
{
	assert(rzr->in_poly);
	struct rzr_op* op = rzr_op(rzr, RZROP_VERTEX);
	struct rzr_tx* tx = rzr_get_current_tx(rzr);
	float u,v;
	rzr_map_point(tx, x, y, &u, &v);
	op->vertex.x = rzr_float_to_int(u);
	op->vertex.y = rzr_float_to_int(v);
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

struct rzr__xlines {
	struct rzr* rzr;
	float  corner_coords[8];
	float  xmin, xmax;
};

static void rzr__xlines_init(struct rzr__xlines* rl, struct rzr* rzr)
{
	memset(rl, 0, sizeof *rl);
	rl->rzr = rzr;
	struct rzr_tx tx = *rzr_get_current_tx(rzr);
	struct rzr_tx invtx = tx;
	rzr_invert_tx(&invtx);
	const int ss = rzr_get_supersampling_factor(rzr);
	const float w = rzr->virtual_width*ss;
	const float h = rzr->virtual_height*ss;
	float* cc = rl->corner_coords;
	rzr_map_point(&invtx, 0, 0, &cc[0], &cc[1]); cc+=2;
	rzr_map_point(&invtx, w, 0, &cc[0], &cc[1]); cc+=2;
	rzr_map_point(&invtx, w, h, &cc[0], &cc[1]); cc+=2;
	rzr_map_point(&invtx, 0, h, &cc[0], &cc[1]); cc+=2;
	for (int i = 0; i < 8; i += 2) {
		const float x = rl->corner_coords[i];
		if (i == 0 || x < rl->xmin) rl->xmin = x;
		if (i == 0 || x > rl->xmax) rl->xmax = x;
	}
}

static int rzr__xline(struct rzr__xlines* rl, float x0, float x1)
{
	assert(x0 < x1);
	if (rl->xmax < x0 || rl->xmin > x1) return 0;
	int cprev = 3;
	const float* cc = rl->corner_coords;
	struct rzr* rzr = rl->rzr;
	rzr_begin_poly(rzr);
	for (int c = 0; c < 4; c++) {
		const float cx = cc[(c<<1)];
		const float cy = cc[(c<<1)+1];
		const float prevcx = cc[(cprev<<1)];
		const float prevcy = cc[(cprev<<1)+1];
		cprev = c;
		const float cxmin = cx < prevcx ? cx : prevcx;
		const float cxmax = cx > prevcx ? cx : prevcx;
		const int previnside = x0 <= prevcx && prevcx <= x1;
		if (previnside) {
			rzr_vertex(rzr, prevcx, prevcy);
		}
		if (cx == prevcx) continue;
		const int x0vis = cxmin <= x0 && x0 <= cxmax;
		const int x1vis = cxmin <= x1 && x1 <= cxmax;
		if (!x0vis && !x1vis) continue;
		const float y0 = !x0vis ? 0.0f : prevcy + ((x0-prevcx)/(cx-prevcx)) * (cy-prevcy);
		const float y1 = !x1vis ? 0.0f : prevcy + ((x1-prevcx)/(cx-prevcx)) * (cy-prevcy);
		if (x0vis && x1vis) {
			if (cx > prevcx) {
				if (!previnside) {
					rzr_vertex(rzr, x0, y0);
				}
				rzr_vertex(rzr, x1, y1);
			} else {
				if (!previnside) {
					rzr_vertex(rzr, x1, y1);
				}
				rzr_vertex(rzr, x0, y0);
			}
		} else if (x0vis && !x1vis) {
			rzr_vertex(rzr, x0, y0);
		} else if (x1vis && !x0vis) {
			rzr_vertex(rzr, x1, y1);
		} else {
			assert(!"unreachable");
		}
	}
	rzr_end_poly(rzr);
	return 1;
}

static inline void rzr_pattern(struct rzr* rzr, float* ws)
{
	int nw=0;
	float totw = 0.0f;
	float* p = ws;
	while (*p > 0) { totw+=*p ; p++; nw++; }
	assert(((nw&1) == 0) && "ws count must be even");
	assert((nw >= 2) && "empty pattern");
	assert(totw > 0.0f);
	const float itotw = 1.0f/totw;

	struct rzr__xlines rl;
	rzr__xlines_init(&rl, rzr);
	float i0 = floorf(rl.xmin * itotw);
	int ne = 0;
	for (;;) {
		const float x0 = i0*totw;
		if (x0 > rl.xmax) break;
		float px0 = x0;
		for (int i = 0; i < nw; i++) {
			float px1 = px0 + ws[i];
			if (((i&1)==0) && rzr__xline(&rl, px0, px1)) {
				if (ne > 0) rzr_union(rzr);
				ne++;
			}
			px0 = px1;
		}
		const int i0p = i0;
		i0 += 1.0f;
		if (!(i0 > i0p)) break; // paranoia
	}
	if (ne == 0) rzr_zero(rzr);
}


static inline void rzr_line(struct rzr* rzr, float width)
{
	const float wh = 0.5f*width;
	struct rzr__xlines rl;
	rzr__xlines_init(&rl, rzr);
	if (!rzr__xline(&rl, -wh, wh)) rzr_zero(rzr);
}

static inline void rzr_split(struct rzr* rzr)
{
	struct rzr__xlines rl;
	rzr__xlines_init(&rl, rzr);
	if (rl.xmax <= 0.0f || !rzr__xline(&rl, 0.0f, rl.xmax)) rzr_zero(rzr);
}

static inline void rzr_box(struct rzr* rzr, float width, float height)
{
	rzr_begin_poly(rzr);
	const float wh = 0.5f*width;
	const float hh = 0.5f*height;
	rzr_vertex(rzr, -wh, -hh);
	rzr_vertex(rzr,  wh, -hh);
	rzr_vertex(rzr,  wh,  hh);
	rzr_vertex(rzr, -wh,  hh);
	rzr_end_poly(rzr);
}

static inline void rzr_rounded_box(struct rzr* rzr, float width, float height, float radius)
{
	if (radius <= 0.0f) {
		rzr_box(rzr, width, height);
		return;
	}

	// FIXME: handle width < radius*2 and height < radius*2

	const float wh = 0.5f*width;
	const float hh = 0.5f*height;

	rzr_begin_poly(rzr);
	rzr_vertex(rzr , -wh        , -hh+radius );
	rzr_vertex(rzr , -wh+radius , -hh        );
	rzr_vertex(rzr ,  wh-radius , -hh        );
	rzr_vertex(rzr ,  wh        , -hh+radius );
	rzr_vertex(rzr ,  wh        ,  hh-radius );
	rzr_vertex(rzr ,  wh-radius ,  hh        );
	rzr_vertex(rzr , -wh+radius ,  hh        );
	rzr_vertex(rzr , -wh        ,  hh-radius );
	rzr_end_poly(rzr);

	for (int corner = 0; corner < 4; corner++) {
		rzr_tx_save(rzr);
		switch (corner) {
		case 0: rzr_tx_translate(rzr, -wh+radius, -hh+radius); break;
		case 1: rzr_tx_translate(rzr,  wh-radius, -hh+radius); break;
		case 2: rzr_tx_translate(rzr,  wh-radius,  hh-radius); break;
		case 3: rzr_tx_translate(rzr, -wh+radius,  hh-radius); break;
		default: assert(!"unreachable");
		}
		rzr_circle(rzr, radius);
		rzr_union(rzr);
		rzr_tx_restore(rzr);
	}
}

static inline void rzr_arc(struct rzr* rzr, float aperture_degrees, float radius, float width)
{
	if (aperture_degrees < 0.0f) aperture_degrees = 0.0f;
	if (aperture_degrees >= 180.0f) aperture_degrees = 180.0f;
	rzr_tx_save(rzr);
	rzr_tx_rotate(rzr, -(aperture_degrees));
	rzr_split(rzr);
	rzr_tx_restore(rzr);
	rzr_tx_save(rzr);
	rzr_tx_rotate(rzr, (aperture_degrees)+180.0f);
	rzr_split(rzr);
	rzr_tx_restore(rzr);
	if (aperture_degrees < 90.0f) {
		rzr_intersection(rzr);
	} else {
		rzr_union(rzr);
	}
	const float wh = 0.5f*width;
	rzr_circle(rzr, radius+wh);
	rzr_circle(rzr, radius-wh);
	rzr_difference(rzr);
	rzr_intersection(rzr);

	rzr_tx_save(rzr);
	rzr_tx_rotate(rzr, -aperture_degrees);
	rzr_tx_translate(rzr, 0, -radius);
	rzr_circle(rzr, wh);
	rzr_union(rzr);
	rzr_tx_restore(rzr);

	rzr_tx_save(rzr);
	rzr_tx_rotate(rzr, aperture_degrees);
	rzr_tx_translate(rzr, 0, -radius);
	rzr_circle(rzr, wh);
	rzr_union(rzr);
	rzr_tx_restore(rzr);
}

static inline void rzr_capsule(struct rzr* rzr, float x0, float y0, float x1, float y1, float r0, float r1)
{
	const float dx = x1-x0;
	const float dy = y1-y0;
	if (dx == 0.0f && dy == 0.0f) {
		rzr_circle(rzr, r0>r1?r0:r1);
		return;
	}
	const float a = -atan2(dy, dx);
	const float b = (r0==r1) ? 0.0f : asinf((r1-r0) / sqrtf(dx*dx + dy*dy));
	const float s0 = sinf(a-b);
	const float s1 = sinf(a+b);
	const float c0 = cosf(a-b);
	const float c1 = cosf(a+b);
	rzr_begin_poly(rzr);
	rzr_vertex(rzr, x1 + r1*s0, y1 + r1*c0);
	rzr_vertex(rzr, x0 + r0*s0, y0 + r0*c0);
	rzr_vertex(rzr, x0 - r0*s1, y0 - r0*c1);
	rzr_vertex(rzr, x1 - r1*s1, y1 - r1*c1);
	rzr_end_poly(rzr);

	rzr_tx_save(rzr);
	rzr_tx_translate(rzr, x0, y0);
	rzr_circle(rzr, r0);
	rzr_tx_restore(rzr);
	rzr_union(rzr);

	rzr_tx_save(rzr);
	rzr_tx_translate(rzr, x1, y1);
	rzr_circle(rzr, r1);
	rzr_tx_restore(rzr);
	rzr_union(rzr);
}

static inline void rzr_segment(struct rzr* rzr, float x0, float y0, float x1, float y1, float r)
{
	rzr_capsule(rzr, x0, y0, x1, y1, r, r);
}

static inline void rzr_triangle(struct rzr* rzr, float w, float h)
{
	rzr_begin_poly(rzr);
	rzr_vertex(rzr,  0, 0);
	rzr_vertex(rzr,  w, h);
	rzr_vertex(rzr, -w, h);
	rzr_end_poly(rzr);
}

static inline void rzr_trapezoid(struct rzr* rzr, float r1, float r2, float h)
{
	if (r1 <= 0.0f) {
		rzr_triangle(rzr, r2, h);
		return;
	}
	rzr_begin_poly(rzr);
	rzr_vertex(rzr, -r1, 0);
	rzr_vertex(rzr,  r1, 0);
	rzr_vertex(rzr,  r2, h);
	rzr_vertex(rzr, -r2, h);
	rzr_end_poly(rzr);
}

void rzr_render(struct rzr*, size_t scratch_cap, void* scratch_stor, int stride, uint8_t* pixels);

void rzr_subpixel_queries(struct rzr*, int n, int* subpixel_coord_pairs, int* out_inside);
int  rzr_subpixel_query(struct rzr*, int subx, int suby);
void rzr_queries(struct rzr*, int n, int* coord_pairs, int* out_inside);
int  rzr_query(struct rzr*, int x, int y);

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

#define Zero()                       rzr_zero(RZR_INSTANCE)
#define One()                        rzr_one(RZR_INSTANCE)

#define Circle(r)                    rzr_circle(RZR_INSTANCE,r)

#define BeginPoly()                  rzr_begin_poly(RZR_INSTANCE)
#define Vertex(x,y)                  rzr_vertex(RZR_INSTANCE,x,y)
#define EndPoly()                    rzr_end_poly(RZR_INSTANCE)

#define Star(n,o,i)                  rzr_star(RZR_INSTANCE,n,o,i)
#define Pattern(...)                 rzr_pattern(RZR_INSTANCE,(float[]) { __VA_ARGS__, 0 })
#define Line(w)                      rzr_line(RZR_INSTANCE,w)
#define Split()                      rzr_split(RZR_INSTANCE)
#define Arc(a,r,w)                   rzr_arc(RZR_INSTANCE,a,r,w)
#define Box(w,h)                     rzr_box(RZR_INSTANCE,w,h)
#define RoundedBox(w,h,r)            rzr_rounded_box(RZR_INSTANCE,w,h,r)
#define Segment(x0,y0,x1,y1,r)       rzr_segment(RZR_INSTANCE,x0,y0,x1,y1,r)
#define Capsule(x0,y0,x1,y1,r0,r1)   rzr_capsule(RZR_INSTANCE,x0,y0,x1,y1,r0,r1)
#define Triangle(w,h)                rzr_triangle(RZR_INSTANCE,w,h)
#define Trapezoid(r1,r2,h)           rzr_trapezoid(RZR_INSTANCE,r1,r2,h)

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

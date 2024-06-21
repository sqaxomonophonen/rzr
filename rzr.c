#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <stdio.h> // XXX remove

#include "rzr.h"

#define ARRAY_LENGTH(xs) (sizeof(xs)/sizeof((xs)[0]))

#define SPANLIST_MAX_LENGTH (1<<14)

struct span {
	int w,x;
};

// returns 1 if [a0;a1] overlaps [b0;b1]; both endpoints are inclusive, so
// [a0;a1] contains both a0 and a1 but not a0-1 and a1+1
static inline int ispans_overlap(int a0, int a1, int b0, int b1)
{
	return a0 <= b1 && b0 <= a1;
}

static inline int spans_overlap(struct span a, struct span b)
{
	const int a0 = a.x;
	const int a1 = (a0 + (int)a.w) - 1;
	const int b0 = b.x;
	const int b1 = (b0 + (int)b.w) - 1;
	return ispans_overlap(a0, a1, b0, b1);
}

static inline int merged_spans(struct span* a, const struct span* b)
{
	assert((a->x <= b->x) && "bad input");
	struct span a1 = *a;
	a1.w++;
	if (!spans_overlap(a1, *b)) return 0;
	const int nw = (b->x - a->x) + b->w;
	if (nw > a->w) a->w = nw;
	return 1;
}

static inline void add_span(struct span* spanlist, int spanlist_cap, int* n, struct span* span)
{
	assert(span != NULL);
	if (*n > 0 && merged_spans(&spanlist[(*n)-1], span)) return;
	assert(*n < spanlist_cap);
	spanlist[(*n)++] = *span;
}

static inline void add_span_imm(struct span* spanlist, int spanlist_cap, int* n, int x, int w)
{
	assert(0 <= x && x < (1<<16));
	assert(1 <= w && w < (1<<16));
	struct span span = {.x=x,.w=w};
	add_span(spanlist, spanlist_cap, n, &span);
}

static int spanlist_union(struct span* out_spanlist, int out_spanlist_cap, int na, struct span* aa, int nb, struct span* bb)
{
	int ia = 0;
	int ib = 0;
	int nc = 0;
	for (;;) {
		struct span* a = ia < na ? &aa[ia] : NULL;
		struct span* b = ib < nb ? &bb[ib] : NULL;
		struct span* add = NULL;
		if (a && b) {
			if (a->x < b->x) {
				add = a;
				ia++;
			} else {
				add = b;
				ib++;
			}
		} else {
			assert(!a || !b);
			break;
		}
		assert(add != NULL);
		add_span(out_spanlist, out_spanlist_cap, &nc, add);
	}

	if (ia < na) {
		assert(ib == nb);
		for (; ia < na; ia++) add_span(out_spanlist, out_spanlist_cap, &nc, &aa[ia]);
	} else if (ib < nb) {
		assert(ia == na);
		for (; ib < nb; ib++) add_span(out_spanlist, out_spanlist_cap, &nc, &bb[ib]);
	}

	assert(ia == na);
	assert(ib == nb);
	return nc;
}

static int spanlist_intersection(struct span* out_spanlist, int out_spanlist_cap, int na, struct span* aa, int nb, struct span* bb)
{
	int ia = 0;
	int ib = 0;
	int nc = 0;
	while (ia < na && ib < nb) {
		struct span* a = &aa[ia];
		struct span* b = &bb[ib];
		if (spans_overlap(*a, *b)) {
			assert(nc < out_spanlist_cap);
			struct span* c0 = nc > 0 ? &out_spanlist[nc-1] : NULL;
			struct span* c = &out_spanlist[nc++];
			const int x0 = a->x > b->x ? a->x : b->x;
			c->x = x0;
			const int ax1 = (int)a->x + (int)a->w;
			const int bx1 = (int)b->x + (int)b->w;
			c->w = ax1<bx1 ? ax1-x0 : bx1-x0;
			if (c0 != NULL && merged_spans(c0, c)) nc--;
			if (ax1 < bx1) {
				ia++;
			} else {
				ib++;
			}
		} else if (a->x < b->x) {
			ia++;
		} else {
			ib++;
		}
	}
	assert((ia == na) || (ib == nb));
	return nc;
}

#if 1
// NOTE there are alternative implementations below
static int spanlist_difference(struct span* out_spanlist, int out_spanlist_cap, int na, struct span* aa, int nb, struct span* bb)
{
	int nc = 0;
	int ib = 0;
	int ia = 0;
	while ((ia < na) && (ib < nb)) {
		struct span* a = &aa[ia];
		struct span* b = &bb[ib];
		const int ax1 = (int)a->x + (int)a->w;
		int bx1 = (int)b->x + (int)b->w;
		if (b->x >= ax1) {
			add_span(out_spanlist, out_spanlist_cap, &nc, a);
			ia++;
			continue;
		} else if (a->x >= bx1) {
			ib++;
			continue;
		} else {
			struct span rm = *a;
			for (;;) {
				const int rmx0 = (int)rm.x;
				const int rmx1 = rmx0 + (int)rm.w;
				if (b->x >= rmx1) {
					add_span(out_spanlist, out_spanlist_cap, &nc, &rm);
					break;
				}
				if (rmx0 < b->x) {
					assert(spans_overlap(rm, *b));
					add_span_imm(out_spanlist, out_spanlist_cap, &nc, rmx0, b->x - rmx0);
				}
				bx1 = (int)b->x + (int)b->w;
				if (rmx1 > bx1) {
					assert(spans_overlap(rm, *b));
					rm.x = bx1;
					rm.w = rmx1-bx1;
					ib++;
					if (ib == nb) {
						add_span(out_spanlist, out_spanlist_cap, &nc, &rm);
						break;
					}
					assert(ib <= nb);
					b = &bb[ib];
					continue;
				}
				break;
			}
			ia++;
		}
	}
	for (; ia < na; ia++) {
		add_span(out_spanlist, out_spanlist_cap, &nc, &aa[ia]);
	}
	return nc;
}
#endif

#if 0
static int spanlist_invert(struct span* out_spanlist, int out_spanlist_cap, int x0, int x1, int na, struct span* aa)
{
	int nb = 0;
	int x = x0;
	for (int i = 0; i <= na; i++) {
		int w=0,nx=-1;
		if (i < na) {
			struct span a = aa[i];
			w = a.x - x;
			nx = a.x + a.w;
		} else {
			assert(i == na);
			w = x1 - x;
		}
		if (w > 0) {
			add_span_imm(out_spanlist, out_spanlist_cap, &nb, x, w);
		}
		x = nx;
	}
	return nb;
}

// DIFFERENCE(a,b) == INTERSECTION(a,INVERT(b))
// NOTE: this implementation is "beautiful" because it reuses
// spanlist_intersection(), but it's also "ugly" because it uses local storage
// ("bbb") which hard to dimension (too low and it might be a problem for valid
// input; too high and it might bust the stack).
static int spanlist_difference(struct span* out_spanlist, int out_spanlist_cap, int na, struct span* aa, int nb, struct span* bb)
{
	if (na == 0) return 0;
	assert(na > 0);
	struct span bbb[SPANLIST_MAX_LENGTH+1]; // XXX fragile
	const int x0 = aa[0].x;
	const int x1 = aa[na-1].x + aa[na-1].w;
	const int nbbb = spanlist_invert(bbb, ARRAY_LENGTH(bbb), x0, x1, nb, bb);
	return spanlist_intersection(out_spanlist, out_spanlist_cap, na, aa, nbbb, bbb);
}
#endif


#if 0
// FIXME doesn't work
static int spanlist_difference(struct span* out_spanlist, int out_spanlist_cap, int na, struct span* aa, int nb, struct span* bb)
{
	if (na == 0) return 0;
	assert(na > 0);
	int ia = 0;
	int ib = 0;
	int nc = 0;
	int bx = aa[0].x;
	int nxbx = 0;
	while (ia < na && ib <= nb) {
		struct span* a = &aa[ia];
		struct span b0;
		if (ib < nb) {
			b0 = bb[ib];
		} else {
			b0.x = aa[na-1].x + aa[na-1].w;
			b0.w = 0;
		}
		struct span b = { .x=bx,  .w=b0.x-bx };
		nxbx = b0.x + b0.w;
		if (spans_overlap(*a, b)) {
			assert(nc < out_spanlist_cap);
			struct span* c0 = nc > 0 ? &out_spanlist[nc-1] : NULL;
			struct span* c = &out_spanlist[nc++];
			const int x0 = a->x > b.x ? a->x : b.x;
			c->x = x0;
			const int ax1 = (int)a->x + (int)a->w;
			const int bx1 = (int)b.x + (int)b.w;
			c->w = ax1<bx1 ? ax1-x0 : bx1-x0;
			if (c0 != NULL && merged_spans(c0, c)) nc--;
			if (ax1 < bx1) {
				ia++;
			} else {
				ib++;
				bx = nxbx;
			}
		} else if (a->x < b.x) {
			ia++;
		} else {
			ib++;
			bx = nxbx;
		}
	}
	assert((ia == na) || (ib == nb+1));
	return nc;
}
#endif

static int spanlist_binop(enum rzr_op_code opcode, struct span* out_spanlist, int out_spanlist_cap, int na, struct span* aa, int nb, struct span* bb)
{
	switch (opcode) {
	case RZROP_UNION:                 return spanlist_union(out_spanlist, out_spanlist_cap, na, aa, nb, bb);
	case RZROP_INTERSECTION:   return spanlist_intersection(out_spanlist, out_spanlist_cap, na, aa, nb, bb);
	case RZROP_DIFFERENCE:       return spanlist_difference(out_spanlist, out_spanlist_cap, na, aa, nb, bb);
	default: assert(!"invalid or unhandled spanlist binop case");
	}
}

struct scratch {
	size_t cap;
	void* base;
	size_t cursor, cursor_max;
	size_t tail_alloc_size;
	size_t tail_alloc_alignment;
	int tail_alloc_n;
};

static inline struct scratch scratch_save(struct scratch* s)
{
	assert(s->tail_alloc_size == 0);
	assert(s->tail_alloc_alignment == 0);
	assert(s->tail_alloc_n == 0);
	return *s;
}

static inline void scratch_restore(struct scratch* s, const struct scratch* saved)
{
	assert(s->tail_alloc_size == 0);
	assert(s->tail_alloc_alignment == 0);
	assert(s->tail_alloc_n == 0);
	const size_t cursor_max = saved->cursor_max;
	*s = *saved;
	s->cursor_max = cursor_max;
}

static struct scratch setup_scratch(size_t cap, void* base)
{
	assert(base != NULL);
	return (struct scratch) {
		.cap    = cap,
		.base   = base,
	};
}

static inline void* scratch_allocn_ex(struct scratch* scratch, int is_tail_alloc, int do_clear, size_t alignment, size_t element_size, int n)
{
	if (!is_tail_alloc) {
		assert(scratch->tail_alloc_size == 0);
		assert(scratch->tail_alloc_alignment == 0);
	} else {
		assert(element_size == 0);
		assert(alignment == 0);
		element_size = scratch->tail_alloc_size;
		alignment = scratch->tail_alloc_alignment;
		scratch->tail_alloc_n += n;
	}
	assert(((alignment & (alignment-1)) == 0) && "expected power-of-two alignment");
	assert(((element_size & (alignment-1)) == 0) && "expected aligned element_size");
	if (n == 0) return NULL;
	const size_t cursor0 = ((scratch->cursor + alignment-1) & ~(alignment-1));
	const size_t size = n*element_size;
	const size_t cursor1 = cursor0 + size;
	assert((cursor1 < scratch->cap) && "scratch allocation overflow");
	void* ptr = scratch->base + cursor0;
	if (do_clear) memset(ptr, 0, size);
	scratch->cursor = cursor1;
	if (scratch->cursor > scratch->cursor_max) scratch->cursor_max = scratch->cursor;
	return ptr;
}

#define SCRATCH_ALLOCN(N, TYPE) scratch_allocn_ex(SCRATCHP, 0, 1, _Alignof(TYPE), sizeof(TYPE), N)

#define SCRATCH_AT() (SCRATCHP->base+SCRATCHP->cursor)

#define SCRATCH_TAIL_ALLOC_BEGIN(TYPE) \
	( \
	assert(SCRATCHP->tail_alloc_size == 0), \
	assert(SCRATCHP->tail_alloc_alignment == 0), \
	SCRATCHP->tail_alloc_size = sizeof(TYPE), \
	SCRATCHP->tail_alloc_alignment = _Alignof(TYPE), \
	SCRATCHP->tail_alloc_n = 0, \
	SCRATCH_AT() \
	)

#define SCRATCH_TAIL_ALLOC_END() \
	( \
	assert(SCRATCHP->tail_alloc_size > 0), \
	assert(SCRATCHP->tail_alloc_alignment > 0), \
	SCRATCHP->tail_alloc_size = 0, \
	SCRATCHP->tail_alloc_alignment = 0, \
	SCRATCHP->tail_alloc_n = 0 \
	)

#define SCRATCH_TAIL_OFFSET() (SCRATCHP->tail_alloc_n)
#define SCRATCH_TAIL_REMAINING() (assert(SCRATCHP->tail_alloc_size>0), (SCRATCHP->cap - SCRATCHP->cursor) / SCRATCHP->tail_alloc_size)
#define SCRATCH_TAIL_ALLOC(N) scratch_allocn_ex(SCRATCHP, 1, 1, 0, 0, N)
#define SCRATCH_TAIL_YANK(N)  scratch_allocn_ex(SCRATCHP, 1, 0, 0, 0, N) // doesn't clear allocation; use together with SCRATCH_TAIL_REMAINING() to write ahead

static inline int get_next_skip(int n_sizes, const int* sizes, int current_value)
{
	int acc = 0;
	for (int i = 0; i < n_sizes; i++) {
		const int size = sizes[i];
		if (size < 0) {
			if (acc > current_value) {
				return acc;
			}
		} else {
			assert(size > 0);
			acc += size;
		}
	}
	return -1;
}

static inline int get_next_yskip(struct rzr* rzr, int current_y)
{
	return get_next_skip(rzr->n_heights, rzr->heights, current_y);
}


struct spanline {
	int n;
	struct span* spans;
};

static void subrender(uint8_t* out_pixels, int n_widths, const int* widths, int ss, struct spanline* sub_scanlines)
{
	int cursors[RZR_MAX_SUPERSAMPLING_FACTOR];
	assert(0 < ss && ss <= RZR_MAX_SUPERSAMPLING_FACTOR);
	memset(cursors, 0, ss*sizeof(cursors[0]));
	const int ssss = ss*ss;
	int x = 0;
	int virtual_width = 0;
	for (int i = 0; i < n_widths; i++) {
		const int w = widths[i];
		if (w < 0) continue;
		assert(w > 0);
		virtual_width += w;
	}
	int next_x_skip = get_next_skip(n_widths, widths, -1);
	int x_skip = 0;
	while (x < virtual_width) {
		const int xx0 = x*ss;
		const int xx1 = xx0+ss-1;
		int v_pixel=0, v_tail=0;
		int rlx = -1;
		for (int i = 0; i < ss; i++) {
			struct spanline* line = &sub_scanlines[i];
			for (;;) {
				int c = cursors[i];
				if (c >= line->n) {
					// no more spans
					if (rlx == -1) rlx = virtual_width;
					break;
				}
				struct span* span = &line->spans[c];
				const int sxx0 = span->x;
				assert(sxx0 >= 0);
				assert(span->w > 0);
				const int sxx1 = sxx0 + span->w - 1;
				assert(sxx1 < ss*virtual_width);
				if (sxx1 < xx0) {
					// current span lies before current
					// pixel cell; advance cursor
					cursors[i]++;
					continue;
				}
				if (sxx0 > xx1) {
					// current span lies after current
					// pixel cell; update run-length
					// target.
					const int sx0 = sxx0 / ss;
					if (rlx == -1 || sx0 < rlx) rlx = sx0;
					break;
				}
				assert(ispans_overlap(xx0, xx1, sxx0, sxx1));
				// current span overlaps with current pixel
				// cell. crop span to calculate pixel
				// contribution.
				const int oxx0 = xx0 > sxx0 ? xx0 : sxx0;
				const int oxx1 = xx1 < sxx1 ? xx1 : sxx1;
				assert(oxx1 >= oxx0);
				const int v = oxx1 - oxx0 + 1;
				v_pixel += v;
				if (sxx1 <= xx1) {
					// current span ends inside current
					// pixel cell; continue with next span,
					// same sub scanline.
					cursors[i]++;
					continue;
				} else {
					// current span exits current pixel
					// cell; calculate tail contribution;
					// update run-length target.
					v_tail += ss;
					const int sx1 = sxx1 / ss;
					if (rlx == -1 || sx1 < rlx) rlx = sx1;
					break;
				}
			}
		}

		if (x == next_x_skip) {
			next_x_skip = get_next_skip(n_widths, widths, x);
			x_skip++;
		}
		out_pixels[x+x_skip] = (v_pixel * 255) / ssss;
		x++;
		if (x == next_x_skip) {
			next_x_skip = get_next_skip(n_widths, widths, x);
			x_skip++;
		}
		if (rlx >= 0) {
			int n = rlx-x;
			if (x+n > virtual_width) n = virtual_width-x;
			if (n > 0) {
				const int v = (v_tail * 255) / ssss;
				while (n > 0) {
					int nn = n;
					int skip=0;
					if (next_x_skip >= 0 && x+n >= next_x_skip) {
						nn = next_x_skip-x;
						skip=1;
					}
					assert(nn > 0);
					memset(out_pixels+x+x_skip, v, nn);
					x += nn;
					if (skip) {
						x_skip++;
						next_x_skip = get_next_skip(n_widths, widths, x);
					}
					n -= nn;
				}
			}
		}
	}
}

static int int_compar(const void* va, const void* vb)
{
	return *(const int*)va - *(const int*)vb;
}

static inline int lineseg_eval_x_at_y(int x0, int y0, int x1, int y1, int y)
{
	assert(y1 > y0);
	assert(y0 <= y && y <= y1);
	if (y == y0) return x0;
	if (y == y1) return x1;
	// XXX probably slightly wrong... is it midpoint? also beware that
	// (x1-x0) may be negative, so there are integer division shennanigans?
	const int x = x0 + ((x1-x0) * (y-y0))/(y1-y0);
	return x;
}

struct poly_xedge {
	int y0;
	int x, side, xstep, numerator, denominator;
};

struct poly_xspan {
	struct poly_xedge edge[2];
};

static int poly_xedge_compar(const void* va, const void* vb)
{
	const struct poly_xedge* a = va;
	const struct poly_xedge* b = vb;
	const int dx = a->x - b->x;
	if (dx != 0) return dx;
	const int dside = a->side - b->side;
	return dside;
}

int rzr_subpixel_query(struct rzr* rzr, int subx, int suby)
{
	const int n_ops = rzr->prg_length;
	const size_t stack_cap = 1<<16;
	int stack_height = 0;
	uint8_t stack[stack_cap];
	#define PUSH(v) (assert(stack_height < stack_cap), stack[stack_height++]=(v))
	#define POP()   (assert(stack_height > 0), stack[--stack_height])
	for (int pc = 0; pc < n_ops; pc++) {
		struct rzr_op* op = &rzr->prg[pc];
		const int opcode = op->code;
		switch (opcode) {
		case RZROP_CIRCLE: {
			const int64_t dx = subx - op->circle.cx;
			const int64_t dy = suby - op->circle.cy;
			const int64_t d2 = dx*dx + dy*dy;
			const int64_t r = op->circle.radius;
			const int v = d2 <= r*r;
			PUSH(v);
		}	break;
		case RZROP_POLY: {
			const int n = op->poly.n_vertices;
			int iprev = n-1;
			int acc = 0;
			for (int i = 0; i < n; i++) {
				struct rzr_op* voprev = &rzr->prg[pc+1+iprev];
				iprev = i;
				struct rzr_op* vop = &rzr->prg[pc+1+i];
				assert(voprev->code == RZROP_VERTEX);
				assert(vop->code == RZROP_VERTEX);
				const int p0x = voprev->vertex.x;
				const int p0y = voprev->vertex.y;
				const int p1x = vop->vertex.x;
				const int p1y = vop->vertex.y;
				if (p0y == p1y) continue;
				const int y0 = p0y < p1y ? p0y : p1y;
				const int y1 = p0y > p1y ? p0y : p1y;
				assert(y0 <= y1);
				if (!(y0 <= suby && suby < y1)) continue;
				if (p1y < p0y) {
					// left edge; increment accumulator if
					// point lies on or after edge
					if (subx >= lineseg_eval_x_at_y(p1x, p1y, p0x, p0y, suby)) acc++;
				} else if (p1y > p0y) {
					// right edge; decrement accumulator if
					// point lies after edge
					if (subx > lineseg_eval_x_at_y(p0x, p0y, p1x, p1y, suby)) acc--;
				} else {
					assert(!"unreachable");
				}
			}
			const int v = acc>0;
			PUSH(v);
			pc += n;
		}	break;
		case RZROP_UNION: {
			const int b = POP();
			const int a = POP();
			PUSH(a || b);
		}	break;
		case RZROP_DIFFERENCE: {
			const int b = POP();
			const int a = POP();
			PUSH(a && !b);
		}	break;
		case RZROP_INTERSECTION: {
			const int b = POP();
			const int a = POP();
			PUSH(a && b);
		}	break;
		default: assert(!"unhandled case");
		}
	}
	const int r = POP();
	#undef POP
	#undef PUSH
	return r;
}

void rzr_subpixel_queries(struct rzr* rzr, int n, int* subpixel_coord_pairs, int* out_inside)
{
	int* p = out_inside;
	for (int i = 0; i < n; i++) {
		const int subx = subpixel_coord_pairs[i<<1];
		const int suby = subpixel_coord_pairs[(i<<1)+1];
		const int is_inside = rzr_subpixel_query(rzr, subx, suby);
		if (p != NULL) *(p++) = is_inside;
	}
}

int rzr_query(struct rzr* rzr, int x, int y)
{
	return rzr_subpixel_query(rzr, x*rzr->supersampling_factor, y*rzr->supersampling_factor);
}

void rzr_queries(struct rzr* rzr, int n, int* coord_pairs, int* out_inside)
{
	int* p = out_inside;
	for (int i = 0; i < n; i++) {
		const int x = coord_pairs[i<<1];
		const int y = coord_pairs[(i<<1)+1];
		const int is_inside = rzr_query(rzr, x, y);
		if (p != NULL) *(p++) = is_inside;
	}
}

static void render(struct rzr* rzr, struct scratch* SCRATCHP, int stride, uint8_t* pixels)
{
	const int NIL = -1;

	const int n_ops = rzr->prg_length;

	int* pc2ri = SCRATCH_ALLOCN(n_ops, int);
	int* pc2ysli = SCRATCH_ALLOCN(n_ops, int);
	int n_circle = 0;
	int circle_radius_sum = 0;
	int faux_stack_height = 0;
	int max_faux_stack_height = 0;
	int n_yspan_lists = 0;
	int n_polys = 0;
	int n_remaining_vertices = 0;
	int max_n_vertices_in_poly = 0;
	int max_poly_yspan_count = 0;
	int n_vertices_total = 0;
	// first program pass:
	//  - count stuff for allocation purposes
	//  - assign indices to stuff
	for (int pc = 0; pc < n_ops; pc++) {
		struct rzr_op* op = &rzr->prg[pc];
		int resource_index = NIL;
		int yspan_list_index = NIL;
		int stack_delta = 0;

		if (n_remaining_vertices > 0) {
			assert(op->code == RZROP_VERTEX);
			n_remaining_vertices--;
			continue;
		}

		switch (op->code) {

		case RZROP_CIRCLE: {
			assert(n_remaining_vertices == 0);
			resource_index = n_circle++;
			yspan_list_index = n_yspan_lists++;
			stack_delta = 1; // 1 PUSH
			circle_radius_sum += op->circle.radius;
		}	break;

		case RZROP_POLY: {
			resource_index = n_polys++;
			yspan_list_index = n_yspan_lists++;
			stack_delta = 1; // 1 PUSH
			assert(n_remaining_vertices == 0);
			const int n = op->poly.n_vertices;
			assert(n >= 3);
			n_remaining_vertices = n;
			max_poly_yspan_count += (n-1);
			n_vertices_total += n;
			if (n > max_n_vertices_in_poly) max_n_vertices_in_poly = n;
		}	break;

		case RZROP_VERTEX: {
			assert(!"NOT EXPECTED: error in RZROP_POLY handler?");
		}	break;

		case RZROP_UNION:
		case RZROP_INTERSECTION:
		case RZROP_DIFFERENCE: {
			assert(n_remaining_vertices == 0);
			stack_delta = (-2+1); // 2 POPs, 1 PUSH
			yspan_list_index = n_yspan_lists++;
		}	break;

		default: assert(!"unhandled case");
		}

		faux_stack_height += stack_delta;
		if (faux_stack_height > max_faux_stack_height) max_faux_stack_height = faux_stack_height;

		pc2ri[pc] = resource_index;
		pc2ysli[pc] = yspan_list_index;
	}
	assert(n_remaining_vertices == 0);

	int** circles_xs = SCRATCH_ALLOCN(n_circle, int*);
	int* xs_storage = SCRATCH_ALLOCN(circle_radius_sum, int);

	struct vertex { int x,y; };
	struct vertex* vertex_stor = SCRATCH_ALLOCN(n_vertices_total, struct vertex);
	int vertex_cursor = 0;

	struct poly {
		int  yspan_offset,  yspan_count;
		int vertex_offset, vertex_count;
	};
	struct poly* polys = SCRATCH_ALLOCN(n_polys, struct poly);

	struct poly_yspan {
		int y0,y1;
		int xspan_offset, xspan_count;
	};
	struct poly_yspan* poly_yspans = SCRATCH_ALLOCN(max_poly_yspan_count, struct poly_yspan);
	int poly_yspan_cursor = 0;

	int* vertex_y_stor = SCRATCH_ALLOCN(max_n_vertices_in_poly, int);

	const int stack_cap = max_faux_stack_height;
	int stack_height = 0;
	int* stack = stack_cap == 0 ? NULL : SCRATCH_ALLOCN(stack_cap, int);
	#define PUSH(v) (assert(stack_height < stack_cap), stack[stack_height++]=(v))
	#define POP()   (assert(stack_height > 0), stack[--stack_height])

	struct yspan_list { int offset,n; };
	struct yspan {
		int y0,y1;
		int pc;
		int ysi_a, ysi_b;
	};

	struct yspan_list* yspan_lists = n_yspan_lists == 0 ? NULL : SCRATCH_ALLOCN(n_yspan_lists, struct yspan_list);
	struct yspan* yspan_stor = SCRATCH_TAIL_ALLOC_BEGIN(struct yspan);

	#define NEW_YSPAN_LIST() \
	{ \
		yspan_lists[yspan_list_index].offset = SCRATCH_TAIL_OFFSET(); \
		yspan_lists[yspan_list_index].n = 0; \
		yspan_list_index++; \
	}

	#define PUSH_YSPAN(Y0,Y1,YSLI_A,YSLI_B) \
	{ \
		struct yspan* yspan = SCRATCH_TAIL_ALLOC(1); \
		yspan->pc = pc; \
		yspan->y0 = Y0; \
		yspan->y1 = Y1; \
		yspan->ysi_a = YSLI_A; \
		yspan->ysi_b = YSLI_B; \
		yspan_lists[yspan_list_index-1].n++; \
		assert(yspan_lists[yspan_list_index-1].n < 2 || yspan->y0 > (yspan-1)->y1); \
		/* printf("push yspan(pc=%d,[%d;%d],a=%d,b=%d)\n", yspan->pc, yspan->y0, yspan->y1, yspan->ysi_a, yspan->ysi_b); */ \
	}

	// second program pass:
	//  - do precalc for shapes
	//  - create and link y-spans into an AST-like representation by
	//    running the program (stack operations only)
	int yspan_list_index = 0;
	for (int pc = 0; pc < n_ops; pc++) {
		struct rzr_op* op = &rzr->prg[pc];
		const int ri = pc2ri[pc];

		const int opcode = op->code;
		switch (opcode) {

		case RZROP_POLY: {
			PUSH(pc);

			struct poly* poly = &polys[ri];
			poly->yspan_offset = poly_yspan_cursor;

			const int n = op->poly.n_vertices;
			assert((pc+n+1) <= n_ops);
			poly->vertex_offset = vertex_cursor;
			poly->vertex_count = n;
			for (int i = 0; i < n; i++) {
				struct rzr_op* vop = &rzr->prg[pc+i+1];
				assert(vop->code == RZROP_VERTEX);
				vertex_y_stor[i] = vop->vertex.y;
				assert(vertex_cursor < n_vertices_total);
				struct vertex* v = &vertex_stor[vertex_cursor++];
				v->x = vop->vertex.x;
				v->y = vop->vertex.y;
			}
			qsort(vertex_y_stor, n, sizeof(vertex_y_stor[0]), int_compar);
			int n_spans = 0;
			NEW_YSPAN_LIST();
			for (int i = 1; i < n; i++) {
				const int y0 = vertex_y_stor[i-1];
				const int y1 = vertex_y_stor[i];
				if (y1 > y0) {
					n_spans++;
					PUSH_YSPAN(y0, y1-1, NIL, NIL);
					assert(poly_yspan_cursor < max_poly_yspan_count);
					struct poly_yspan* poly_yspan = &poly_yspans[poly_yspan_cursor++];
					poly_yspan->y0 = y0;
					poly_yspan->y1 = y1-1;
				}
			}
			poly->yspan_count = n_spans;

			pc += n; // skip RZROP_VERTEX's
		}	break;

		case RZROP_VERTEX: {
			assert(!"NOT EXPECTED: error in RZROP_POLY handler?");
		}	break;

		case RZROP_CIRCLE: {
			PUSH(pc);

			const int radius = op->circle.radius;
			assert((radius > 0) && "XXX can't handle zero-radius circle here (do it elsewhere?)");

			NEW_YSPAN_LIST();
			PUSH_YSPAN(op->circle.cy-(radius-1), op->circle.cy+(radius-1), NIL, NIL);

			// "cache" quarter circle
			assert(xs_storage != NULL);
			int* xs = xs_storage;
			xs_storage += radius;

			assert(0 <= ri && ri < n_circle);
			assert(circles_xs != NULL);
			circles_xs[ri] = xs;

			// Thank you Casey Muratori for this lovely code :-) ( https://www.computerenhance.com/p/efficient-dda-circle-outlines )
			const int r2 = radius+radius;
			int y = 0;
			int x = radius;
			int dy = -2;
			int dx = r2+r2-4;
			int d = r2-1;

			int* xp = xs;
			int* xpend = xp + radius;
			while (y <= x) {
				assert(xs <= xp && xp < xpend);
				assert(SCRATCHP->base <= (void*)xp && (void*)xp < (SCRATCHP->base+SCRATCHP->cap));
				#if 0
				printf("WRITE/Y xs[%zd] = %d\n", xp-xs, x);
				#endif
				*(xp++) = x;
				d += dy;
				dy -= 4;
				y++;
				if (d < 0) {
					d += dx;
					dx -= 4;
					x--;
					assert(0 <= x && x < radius);
					xs[x] = y;
					assert(x < radius);
					#if 0
					printf("WRITE/X xs[%d] = %d\n", x, y);
					#endif
				}
			}
			#if 0
			for (int y = 0; y < radius; y++) {
				const int cx = xs[y];
				printf("%.3d,%.3d ", cx, y);
				for (int x = 0; x < radius; x++) {
					printf("%s", cx > x ? "[]" : "..");
				}
				printf("\n");
			}
			#endif

		}	break;

		case RZROP_UNION:
		case RZROP_INTERSECTION:
		case RZROP_DIFFERENCE: {
			const int ysli_b = pc2ysli[POP()];
			const int ysli_a = pc2ysli[POP()];
			PUSH(pc);

			assert(0 <= ysli_a && ysli_a < yspan_list_index);
			assert(0 <= ysli_b && ysli_b < yspan_list_index);

			NEW_YSPAN_LIST();

			struct yspan_list* ysl_a = &yspan_lists[ysli_a];
			struct yspan_list* ysl_b = &yspan_lists[ysli_b];

			int ca=0, cb=0;
			int crop_A_at=0, crop_B_at=0;
			int do_crop_A=0, do_crop_B=0;
			while (ca < ysl_a->n && cb < ysl_b->n) {
				const int ysi_a = ysl_a->offset+ca;
				const int ysi_b = ysl_b->offset+cb;

				struct yspan ysa = yspan_stor[ysi_a];
				if (do_crop_A) ysa.y0 = crop_A_at;

				struct yspan ysb = yspan_stor[ysi_b];
				if (do_crop_B) ysb.y0 = crop_B_at;

				struct ospan {
					enum { A, B, AB } input;
					int y0, y1;
				};
				const int max_ospans = 2;
				struct ospan ospans[max_ospans];
				int n_ospans = 0;

				#define PUSH_OSPAN(INPUT,Y0,Y1) \
				{ \
					assert(n_ospans < max_ospans); \
					struct ospan* ospan = &ospans[n_ospans++]; \
					ospan->input = INPUT; \
					ospan->y0 = Y0; \
					ospan->y1 = Y1; \
					assert(n_ospans < 2 || ospan->y0 > (ospan-1)->y1); \
					assert(A <= ospan->input && ospan->input <= AB); \
					assert(ospan->y1 >= ospan->y0); \
				}
				const int a0 = ysa.y0;
				const int a1 = ysa.y1;
				const int b0 = ysb.y0;
				const int b1 = ysb.y1;
				do_crop_A = do_crop_B = 0;
				if (ysa.y1 < ysb.y0) {
					// AA|..
					// ..|BB
					PUSH_OSPAN(A, a0, a1);
					ca++;
				} else if (ysb.y1 < ysa.y0) {
					// ..|AA
					// BB|..
					PUSH_OSPAN(B, b0, b1);
					cb++;
				} else if (ysa.y0 < ysb.y0 && ysa.y1 > ysb.y1) {
					// AAA|A
					// .BB|.
					PUSH_OSPAN(A,  a0,   b0-1);
					PUSH_OSPAN(AB, b0,   b1);
					crop_A_at =    b1+1;
					do_crop_A = 1;
					cb++;
				} else if (ysb.y0 < ysa.y0 && ysb.y1 > ysa.y1) {
					// .AA|.
					// BBB|B
					PUSH_OSPAN(B,  b0,   a0-1);
					PUSH_OSPAN(AB, a0,   a1);
					crop_B_at =    a1+1;
					do_crop_B = 1;
					ca++;
				} else if (ysa.y0 < ysb.y0 && ysa.y1 < ysb.y1) {
					// AAA|.
					// .BB|B
					PUSH_OSPAN(A,  a0,   b0-1);
					PUSH_OSPAN(AB, b0,   a1);
					crop_B_at =    a1+1;
					do_crop_B = 1;
					ca++;
				} else if (ysb.y0 < ysa.y0 && ysb.y1 < ysa.y1) {
					// .AA|A
					// BBB|.
					PUSH_OSPAN(B,  b0,   a0-1);
					PUSH_OSPAN(AB, a0,   b1);
					crop_A_at =    b1+1;
					do_crop_A = 1;
					cb++;
				} else if (ysa.y0 == ysb.y0 && ysa.y1 > ysb.y1) {
					// AA|AA
					// BB|..
					PUSH_OSPAN(AB, b0,   b1);
					crop_A_at =    b1+1;
					do_crop_A = 1;
					cb++;
				} else if (ysa.y0 == ysb.y0 && ysb.y1 > ysa.y1) {
					// AA|..
					// BB|BB
					PUSH_OSPAN(AB, a0,   a1);
					crop_B_at =    a1+1;
					do_crop_B = 1;
					ca++;
				} else if (ysa.y0 < ysb.y0 && ysa.y1 == ysb.y1) {
					// AAAA|
					// ..BB|
					PUSH_OSPAN(A,  a0, b0-1);
					PUSH_OSPAN(AB, b0, b1);
					ca++;
					cb++;
				} else if (ysb.y0 < ysa.y0 && ysa.y1 == ysb.y1) {
					// ..AA|
					// BBBB|
					PUSH_OSPAN(B,  b0, a0-1);
					PUSH_OSPAN(AB, a0, a1);
					ca++;
					cb++;
				} else if (ysa.y0 == ysb.y0 && ysa.y1 == ysb.y1) {
					// AAAA|
					// BBBB|
					PUSH_OSPAN(AB, a0, a1);
					ca++;
					cb++;
				} else {
					assert(!"unhandled span-vs-span case");
				}

				#undef PUSH_OSPAN

				// assert that we "agree" on overlaps
				int expect_overlap =  0;
				for (int i = 0; i < n_ospans; i++) {
					if (ospans[i].input == AB) {
						expect_overlap = 1;
						break;
					}
				}
				assert(expect_overlap == ispans_overlap(a0,a1,b0,b1));

				for (int i = 0; i < n_ospans; i++) {
					struct ospan* ospan = &ospans[i];
					int a=0,b=0,discard=0;;
					if (ospan->input == AB) {
						a=b=1; // handle AB operation
					} else if (ospan->input == A && (opcode == RZROP_UNION || opcode == RZROP_DIFFERENCE)) {
						a=1; // pass A as-is
					} else if (ospan->input == B && opcode == RZROP_UNION) {
						b=1; // pass B as-is
					} else if (ospan->input != AB && opcode == RZROP_INTERSECTION) {
						discard = 1; // A*nil=nil and nil*B=nil
					} else if (ospan->input == B && opcode == RZROP_DIFFERENCE) {
						discard = 1; // nil-b=nil
					} else {
						assert(!"unhandled ospan case");
					}
					if (discard) continue;
					PUSH_YSPAN(ospan->y0, ospan->y1, (a?ysi_a:NIL), (b?ysi_b:NIL));
				}
			}

			if (opcode == RZROP_UNION || opcode == RZROP_DIFFERENCE) while (ca < ysl_a->n) {
				const int ysi_a = ysl_a->offset+ca;
				struct yspan ysa = yspan_stor[ysi_a];
				if (do_crop_A) {
					ysa.y0 = crop_A_at;
					do_crop_A = 0;
				}
				PUSH_YSPAN(ysa.y0, ysa.y1, ysi_a, NIL);
				ca++;
			}

			if (opcode == RZROP_UNION) while (cb < ysl_b->n) {
				const int ysi_b = ysl_b->offset+cb;
				struct yspan ysb = yspan_stor[ysi_b];
				if (do_crop_B) {
					ysb.y0 = crop_B_at;
					do_crop_B = 0;
				}
				PUSH_YSPAN(ysb.y0, ysb.y1, NIL, ysi_b);
				cb++;
			}

			assert( (opcode != RZROP_UNION)        || (ca == ysl_a->n && cb == ysl_b->n));
			assert( (opcode != RZROP_DIFFERENCE)   || (ca == ysl_a->n));
			assert( (opcode != RZROP_INTERSECTION) || (ca == ysl_a->n || cb == ysl_b->n));

		}	break;

		default: assert(!"unhandled case");

		}
	}
	SCRATCH_TAIL_ALLOC_END();
	assert(vertex_cursor == n_vertices_total);
	assert(yspan_list_index == n_yspan_lists);

	// setup stuff for polygon span rendering
	struct poly_xedge* xedges_stor = SCRATCH_TAIL_ALLOC_BEGIN(struct poly_xedge);
	int n_xspans = 0;
	for (int i0 = 0; i0 < n_polys; i0++) {
		struct poly* poly = &polys[i0];
		const int vertex_count  = poly->vertex_count;
		const int vertex_offset = poly->vertex_offset;
		const int yspan_count   = poly->yspan_count;
		const int yspan_offset  = poly->yspan_offset;
		for (int i1 = 0; i1 < yspan_count; i1++) {
			int i2prev = vertex_count-1;
			struct poly_yspan* yspan = &poly_yspans[yspan_offset + i1];
			assert(yspan->y1 >= yspan->y0);
			const int tail_offset = SCRATCH_TAIL_OFFSET();
			assert(((tail_offset&1) == 0) && "expected even offset since allocations must come in pairs");
			yspan->xspan_offset = tail_offset >> 1;
			int n_left=0, n_right=0;
			struct poly_xedge* xedges_begin = SCRATCH_AT();
			for (int i2 = 0; i2 < vertex_count; i2++) {
				struct vertex* pv = &vertex_stor[vertex_offset + i2prev];
				struct vertex* v  = &vertex_stor[vertex_offset + i2];
				i2prev = i2;
				if (pv->y == v->y) continue; // ignore y=k line
				const int ey0 = pv->y < v->y ? pv->y : v->y;
				const int ey1 = pv->y > v->y ? pv->y : v->y;
				assert(ey1 > ey0);
				if (!ispans_overlap(yspan->y0, yspan->y1, ey0, ey1-1)) continue; // edge not in current y-span
				//printf("yspan [%d;%d] contains edge {%d,%d}-{%d,%d}\n", yspan->y0, yspan->y1, pv->x, pv->y, v->x, v->y);
				assert(ey0 <= yspan->y0);
				assert(ey1 >= yspan->y1);
				struct poly_xedge* xedge = SCRATCH_TAIL_ALLOC(1);
				int x0,x1;
				if (v->y < pv->y) {
					// left edge
					n_left++;
					xedge->side = 0;
					x0 = lineseg_eval_x_at_y(v->x, v->y, pv->x, pv->y, yspan->y0);
					x1 = lineseg_eval_x_at_y(v->x, v->y, pv->x, pv->y, yspan->y1);
				} else if (v->y > pv->y) {
					// right edge
					n_right++;
					xedge->side = 1;
					x0 = lineseg_eval_x_at_y(pv->x, pv->y, v->x, v->y, yspan->y0);
					x1 = lineseg_eval_x_at_y(pv->x, pv->y, v->x, v->y, yspan->y1);
				} else {
					assert(!"unreachable");
				}
				assert(xedge != NULL);
				xedge->y0 = yspan->y0;
				xedge->x = x0;
				xedge->denominator = yspan->y1 - yspan->y0 + 1;

				if (x1 > x0) {
					xedge->numerator = x1-x0;
					xedge->xstep = 1;
				} else {
					xedge->numerator = x0-x1;
					xedge->xstep = -1;
				}
				assert(xedge->denominator > 0);
			}
			assert(n_left == n_right);
			assert(n_left > 0 && n_right > 0);
			struct poly_xedge* xedges_end = SCRATCH_AT();
			const int n_xedges = xedges_end-xedges_begin;
			assert(((n_xedges&1) == 0) && "edges must come in pairs to form spans");
			assert(n_xedges >= 2);
			// XXX this feels fragile. can I design it better? I
			// want to group left/right edges into spans, so
			// qsort() is mostly right, but it's hard to make qsort
			// also treat the "side" field correctly. e.g. this is
			// the correct order:
			//   x=0   side=0 (left)
			//   x=100 side=1 (right)
			//   x=100 side=0 (left)
			//   x=200 side=1 (right)
			// but qsort, which orders by [x,side], does this:
			//   x=0   side=0 (left)
			//   x=100 side=0 (left)
			//   x=100 side=1 (right)
			//   x=200 side=1 (right)
			// and the following for-loop fixes that.
			qsort(xedges_begin, n_xedges, sizeof xedges_begin[0], poly_xedge_compar);
			for (int i2 = 1; i2 < n_xedges; i2+=2) {
				if (xedges_begin[i2].side == 1) continue;
				if (xedges_begin[i2].x == xedges_begin[i2+1].x && xedges_begin[i2+1].side == 1) {
					struct poly_xedge tmp = xedges_begin[i2];
					xedges_begin[i2] = xedges_begin[i2+1];
					xedges_begin[i2+1] = tmp;
				}
			}
			for (int i2 = 0; i2 < n_xedges; i2++) {
				struct poly_xedge* xe = &xedges_begin[i2];
				assert(xe->side == (i2&1));
			}
			yspan->xspan_count = n_xedges >> 1;
			n_xspans += yspan->xspan_count;
		}
	}
	SCRATCH_TAIL_ALLOC_END();

	//////////////////////////////////////
	// setup polygon render states
	struct xedge_state {
		int x,xstep;
		int numerator;
		int denominator;
		int accumulator;
	};

	struct xspan_state {
		int y;
		struct xedge_state edge_states[2];
	};

	struct xspan_state* xspan_states_stor = SCRATCH_ALLOCN(n_xspans, struct xspan_state);
	for (int i0 = 0; i0 < n_xspans; i0++) {
		const struct poly_xedge* src_edges = &xedges_stor[i0<<1];
		struct xspan_state* dst_span = &xspan_states_stor[i0];
		dst_span->y = src_edges->y0;
		assert(src_edges[1].y0 == dst_span->y);
		for (int i1 = 0; i1 < 2; i1++) {
			const struct poly_xedge* src = &src_edges[i1];
			assert(src->side == i1);
			struct xedge_state* dst = &dst_span->edge_states[i1];
			dst->x = src->x;
			dst->xstep = src->xstep;
			dst->numerator = src->numerator;
			assert(src->denominator > 0);
			dst->denominator = src->denominator;
			dst->accumulator = src->denominator;
		}
	}

	if (stack_height < 1) return;

	//////////////////////////
	// begin actual render

	const int top = POP();
	struct yspan_list* ysl = &yspan_lists[pc2ysli[top]];

	const int post_order_stack_cap = 1<<10; // XXX do a proper estimate / upper bound thing?
	int* post_order_stack = SCRATCH_ALLOCN(post_order_stack_cap, int);

	struct xop {
		int opcode;
		union {
			struct {
				int cx,cy,radius;
				int resource_index;
			} circle;
			struct {
				int yspans_i0, yspans_i1;
			} poly;
		};
	};

	const int xop_cap = 1<<10; // XXX do a proper upper bound thing?
	struct xop* xops = SCRATCH_ALLOCN(xop_cap, struct xop);
	#if 0
	printf("sizeof(struct xop)=%zd\n", sizeof(struct xop));
	//assert(!"STOP");
	#endif

	const int supsamp = rzr_get_supersampling_factor(rzr);
	assert(supsamp > 0);
	struct spanline* sub_scanlines = SCRATCH_ALLOCN(supsamp, struct spanline);
	//memset(sub_scanlines, 0, sizeof(sub_scanlines[0])*supsamp);

	const int virutal_width_subpix = rzr->virtual_width*supsamp;

	int y_render_min=-1, y_render_max=-1;
	const struct scratch SCRATCH_SNAPSHOT = scratch_save(SCRATCHP);
	struct span* span_stor = SCRATCH_TAIL_ALLOC_BEGIN(struct span);
	uint8_t* ppix = NULL;
	const int yy1 = rzr->virtual_height * supsamp - 1;
	int next_yy_skip = get_next_skip(rzr->n_heights, rzr->heights, -1);
	int yy_skip = 0;
	for (int i = 0; i < ysl->n; i++) {
		const int iroot = ysl->offset+i;
		struct yspan* yroot = &yspan_stor[iroot];
		if (i > 0) assert(yroot->y0 > yspan_stor[iroot-1].y1);
		assert(yroot->pc == top);
		if (yroot->y1 < 0 || yroot->y0 > yy1) continue;
		const int root_y0 = yroot->y0 < 0 ? 0 : yroot->y0;
		const int root_y1 = yroot->y1 > yy1 ? yy1 : yroot->y1;

		int post_order_stack_height = 0;
		int cursor = iroot;
		int n_xop = 0;

		// post-order tree traversal to create compact program
		// in `xops`
		do {
			while (cursor >= 0) {
				const int b = yspan_stor[cursor].ysi_b;
				if (b >= 0) {
					assert(post_order_stack_height < post_order_stack_cap);
					post_order_stack[post_order_stack_height++] = b;
				}
				assert(post_order_stack_height < post_order_stack_cap);
				post_order_stack[post_order_stack_height++] = cursor;
				cursor = yspan_stor[cursor].ysi_a;
			}
			cursor = post_order_stack[--post_order_stack_height];
			const int b = yspan_stor[cursor].ysi_b;
			if (b >= 0 && (post_order_stack_height == 0 ? NIL : post_order_stack[post_order_stack_height-1]) == b) {
				post_order_stack[post_order_stack_height-1] = cursor;
				cursor = b;
			} else {
				assert(cursor >= 0);
				struct yspan v = yspan_stor[cursor];
				const int pc = v.pc;
				assert(v.y0 <= root_y0 && v.y1 >= root_y1);
				struct rzr_op* op = &rzr->prg[pc];

				struct xop* xop = NULL;

				int passthru_binop=0, binop=0;

				#define SHAPE_XOP() \
					assert(v.ysi_a == NIL); \
					assert(v.ysi_b == NIL); \
					assert(n_xop < xop_cap); \
					xop = &xops[n_xop++]; \
					memset(xop, 0, sizeof *xop);

				switch (op->code) {

				case RZROP_CIRCLE: {
					SHAPE_XOP();
					xop->circle.cx = op->circle.cx;
					xop->circle.cy = op->circle.cy;
					xop->circle.radius = op->circle.radius;
					xop->circle.resource_index = pc2ri[pc];
				}	break;

				case RZROP_POLY: {
					SHAPE_XOP();
					struct poly* poly = &polys[pc2ri[pc]];
					const int n = poly->yspan_count;
					const int yspan_offset = poly->yspan_offset;
					struct poly_yspan* yp = &poly_yspans[yspan_offset];
					int i0=-1, i1=-1;
					for (int i = 0; i < n; i++, yp++) {
						const int o = ispans_overlap(root_y0, root_y1, yp->y0, yp->y1);
						if (o) {
							if (i0 == -1) i0 = i;
							i1 = i;
						} else if (!o && i0 >= 0) {
							break;
						}
					}
					assert(i0 >= 0);
					assert(i1 >= 0);
					xop->poly.yspans_i0 = yspan_offset + i0;
					xop->poly.yspans_i1 = yspan_offset + i1;
				}	break;

				case RZROP_VERTEX: {
					assert(!"NOT EXPECTED: error in RZROP_POLY handler?");
				}	break;

				case RZROP_UNION: {
					assert((v.ysi_a >= 0 || v.ysi_b >= 0) && "bad binop/union prep");
					if ((v.ysi_a == NIL) || (v.ysi_b == NIL)) {
						passthru_binop = 1;
					} else {
						assert(v.ysi_a >= 0 && v.ysi_b >= 0);
						binop = 1;
					}
				}	break;

				case RZROP_DIFFERENCE: {
					assert((v.ysi_a >= 0) && "bad binop/difference prep");
					if (v.ysi_b == NIL) {
						passthru_binop = 1;
					} else {
						binop = 1;
					}
				}	break;

				case RZROP_INTERSECTION: {
					assert((v.ysi_a >= 0 && v.ysi_b >= 0) && "bad binop/intersection prep");
					binop = 1;
				}	break;

				default: assert(!"unhandled opcode");
				}

				if (binop) {
					assert(xop == NULL);
					assert(n_xop < xop_cap);
					xop = &xops[n_xop++];
					memset(xop, 0, sizeof *xop);
					assert(n_xop >= 2);
				} else if (passthru_binop) {
					assert(xop == NULL);
				} else {
					assert(xop != NULL);
				}

				if (xop != NULL) {
					xop->opcode = op->code;
				}

				cursor = -1;
			}
		} while (post_order_stack_height > 0);

		// execute xop-program for each sub-scanline in [root_y0;root_y1]
		for (int y = root_y0; y <= root_y1; y++) {
			const int sub_scanline_i = (y % supsamp);
			const int yy = (y / supsamp);
			if (yy == next_yy_skip) {
				next_yy_skip = get_next_skip(rzr->n_heights, rzr->heights, yy);
				yy_skip++;
			}
			assert(0 <= yy && yy < rzr->virtual_height);
			if (y_render_min == -1 || yy < y_render_min) y_render_min = yy;
			if (y_render_max == -1 || yy > y_render_max) y_render_max = yy;

			struct rlist { int offset,count; };
			const int rstack_cap = 1<<13;
			int rstack_height = 0;
			struct rlist rstack[rstack_cap];

			#define NSPAN(N,SPAN0) \
			{ \
				const int _n = N; \
				struct span* _p = SPAN0; \
				for (int _i = 0; _i < _n; _i++) { \
					struct span _s = *(_p++); \
					if (_s.x < 0) { \
						_s.w += _s.x; \
						_s.x = 0; \
					} \
					if (_s.x+_s.w > virutal_width_subpix) _s.w = virutal_width_subpix-_s.x; \
					if (_s.w <= 0) continue; \
					struct span* _dst = SCRATCH_TAIL_ALLOC(1); \
					*_dst = _s; \
				} \
			}

			for (int xpc = 0; xpc < n_xop; xpc++) {
				struct xop* xop = &xops[xpc];
				const int offset0 = SCRATCH_TAIL_OFFSET();
				switch (xop->opcode) {

				case RZROP_CIRCLE: {
					const int cx = xop->circle.cx;
					const int cy = xop->circle.cy;
					const int radius = xop->circle.radius;
					const int ri = xop->circle.resource_index;
					int ay = y-cy;
					if (ay < 0) ay = -ay;
					int* xs = circles_xs[ri];
					assert(0 <= ay && ay < radius);
					const int dx = xs[ay];
					struct span span = {
						.x = cx-dx+1,
						.w = dx+dx-1,
					};
					NSPAN(1, &span);
				}	break;

				case RZROP_POLY: {
					const int yspans_i0 = xop->poly.yspans_i0;
					const int yspans_i1 = xop->poly.yspans_i1;
					struct poly_yspan* ps = NULL;
					for (int i = yspans_i0; i <= yspans_i1; i++) {
						struct poly_yspan* tps = &poly_yspans[i];
						if (!(tps->y0 <= y && y <= tps->y1)) continue;
						ps = tps;
						break;
					}
					assert(ps != NULL);
					const int xspan_offset = ps->xspan_offset;
					const int xspan_count = ps->xspan_count;

					for (int si = 0; si < xspan_count; si++) {
						struct xspan_state* xs = &xspan_states_stor[xspan_offset+si];
						while (xs->y <= y) {
							struct span span = {
								.x = xs->edge_states[0].x,
								.w = xs->edge_states[1].x - xs->edge_states[0].x,
							};
							if (span.w > 0 && xs->y == y) {
								NSPAN(1, &span);
							}
							for (int ei = 0; ei < 2; ei++) {
								struct xedge_state* es = &xs->edge_states[ei];
								assert(es->numerator >= 0);
								es->accumulator += es->numerator;
								const int n = es->accumulator / es->denominator;
								es->x += es->xstep * n;
								es->accumulator -= es->denominator * n;
							}
							xs->y++;
						}
					}

				}	break;

				case RZROP_UNION:
				case RZROP_DIFFERENCE:
				case RZROP_INTERSECTION: {
					assert(rstack_height >= 2);
					struct rlist rb = rstack[--rstack_height];
					struct rlist ra = rstack[--rstack_height];
					struct span* aa = span_stor + ra.offset;
					struct span* bb = span_stor + rb.offset;
					const int nr = spanlist_binop(xop->opcode, SCRATCH_AT(), SCRATCH_TAIL_REMAINING(), ra.count, aa, rb.count, bb);
					SCRATCH_TAIL_YANK(nr);
				}	break;

				default: assert(!"unhandled case");
				}

				const int nspans = SCRATCH_TAIL_OFFSET() - offset0;
				for (int i = 0; i < nspans; i++) {
					struct span s = span_stor[offset0+i];
					assert(s.x >= 0);
					assert(s.w > 0);
					assert(s.x+s.w <= virutal_width_subpix);
				}
				assert(rstack_height < rstack_cap);
				struct rlist* rl = &rstack[rstack_height++];
				rl->offset = offset0;
				rl->count = nspans;
			}

			assert(rstack_height > 0);
			struct rlist rx = rstack[--rstack_height];

			assert(0 <= sub_scanline_i && sub_scanline_i < supsamp);
			struct spanline* ssc = &sub_scanlines[sub_scanline_i];
			ssc->spans = span_stor + rx.offset;
			ssc->n = rx.count;

			ppix = pixels + (yy + yy_skip) * stride;
			if (sub_scanline_i == (supsamp-1)) {
				subrender(ppix, rzr->n_widths, rzr->widths, supsamp, sub_scanlines);
				// reset scratch
				SCRATCH_TAIL_ALLOC_END();
				scratch_restore(SCRATCHP, &SCRATCH_SNAPSHOT);
				span_stor = SCRATCH_TAIL_ALLOC_BEGIN(struct span);
				memset(sub_scanlines, 0, sizeof(sub_scanlines[0])*supsamp);
			}
		}
	}
	if (SCRATCH_TAIL_OFFSET() > 0) {
		subrender(ppix, rzr->n_widths, rzr->widths, supsamp, sub_scanlines);
	}
	SCRATCH_TAIL_ALLOC_END();

	if (y_render_min > 0) {
		uint8_t* p = pixels;
		int next_y_skip = get_next_yskip(rzr, -1);
		for (int y = 0; y < y_render_min; y++) {
			if (y == next_y_skip) {
				next_y_skip = get_next_yskip(rzr, y);
				p += stride;
			}
			memset(p, 0, rzr->bitmap_width);
			p += stride;
		}
	}

	if (y_render_max >= 0) {
		const int h = rzr->virtual_height;
		const int y0 = y_render_max+1;
		int next_y_skip = get_next_yskip(rzr, y0-1);
		uint8_t* p = pixels + y0*stride;
		for (int y = y0; y < h; y++) {
			if (y == next_y_skip) {
				next_y_skip = get_next_yskip(rzr, y);
				p += stride;
			}
			memset(p, 0, rzr->bitmap_width);
			p += stride;
		}
	}

	//printf("allocated/max:    %zd\n", SCRATCH.cursor_max);

	// render "stretch regions" using point queries (only applicable to MxN
	// renders)
	for (int axis = 0; axis < 2; axis++) {
		int vcursor = 0;
		int bcursor = 0;
		const int n_usizes = (axis==0) ? rzr->n_widths : (axis==1) ? rzr->n_heights : 0;
		const int* usizes  = (axis==0) ? rzr->widths   : (axis==1) ? rzr->heights   : NULL;
		for (int i0 = 0; i0 < n_usizes; i0++) {
			const int size = usizes[i0];
			if (size < 0) {
				const int subpos = vcursor * rzr->supersampling_factor;
				int n,xx,yy,dxx,dyy;
				uint8_t* p;
				int pinc;
				int n_vsizes;
				int* vsizes;
				if (axis == 0) {
					p = pixels + bcursor;
					n = rzr->virtual_height;
					xx = subpos;
					dxx = 0;
					yy = 0;
					dyy = 1;
					pinc = stride;
					n_vsizes = rzr->n_heights;
					vsizes = rzr->heights;
				} else if (axis == 1) {
					p = pixels + bcursor*stride;
					n = rzr->virtual_width;
					xx = 0;
					dxx = 1;
					yy = subpos;
					dyy = 0;
					pinc = 1;
					n_vsizes = rzr->n_widths;
					vsizes = rzr->widths;
				} else {
					assert(!"unreachable");
				}

				int next_vskip = get_next_skip(n_vsizes, vsizes, -1);
				for (int i1 = 0; i1 < n; i1++) {
					if (i1 == next_vskip) {
						*(p) = rzr_subpixel_query(rzr, xx, yy) ? 255 : 0;
						next_vskip = get_next_skip(n_vsizes, vsizes, i1);
						assert(next_vskip > i1 || next_vskip == -1);
						p += pinc;
					}
					int acc = 0;
					for (int i2 = 0; i2 < supsamp; i2++) {
						acc += rzr_subpixel_query(rzr, xx, yy);
						xx += dxx;
						yy += dyy;
					}
					*p = (acc * 255) / supsamp;
					p += pinc;
				}
				bcursor++;
			} else {
				assert(size > 0);
				vcursor += size;
				bcursor += size;
			}
		}
	}

	#undef PUSH_YSPAN
	#undef NEW_YSPAN_LIST
	#undef POP
	#undef PUSH
}

void rzr_render(struct rzr* rzr, size_t scratch_cap, void* scratch_stor, int stride, uint8_t* pixels)
{
	struct scratch s = setup_scratch(scratch_cap, scratch_stor);
	render(rzr, &s, stride, pixels);
	if (s.cursor_max > rzr->scratch_alloc_max) rzr->scratch_alloc_max = s.cursor_max;
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
#ifdef UNIT_TEST
// cc -DUNIT_TEST -Wall -O0 -g rzr.c -o test_rzr -lm && ./test_rzr

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

static void assert_valid_spanlist(int n, struct span* spanlist)
{
	for (int i = 0; i < n; i++) {
		struct span* span = &spanlist[i];
		assert((span->w > 0) && "spanlist must not contain empty w=0 spans");
		if (i > 0) {
			struct span* lastspan = &spanlist[i-1];
			assert((span->x >= (lastspan->x + lastspan->w)) && "spanlist must be ordered and non-overlapping");
			assert(!spans_overlap(*span, *lastspan));
		}
	}
}

static void assert_clean_spanlist(int n, struct span* spanlist)
{
	assert_valid_spanlist(n, spanlist);
	for (int i = 1; i < n; i++) {
		struct span* s0 = &spanlist[i-1];
		struct span* s1 = &spanlist[i];
		assert((s1->x > ((int)s0->x + (int)s0->w)) && "there's no gap between s0 and s1; these spans should be merged for the spanlist to be \"clean\"");
	}
}

static int parse_spanlist(struct span* spanlist, int spanlist_cap, const char* str)
{
	enum {
		OUTSIDE_SPAN=0,
		INSIDE_SPAN,
	} state = OUTSIDE_SPAN;
	int n = 0;
	struct span* span = NULL;
	for (const char* p = str; *p; p++) {
		const int x = p-str;
		if (state == OUTSIDE_SPAN && *p == ' ') {
			assert(span == NULL);
		} else if (state == OUTSIDE_SPAN && *p == '[') {
			assert(n < spanlist_cap);
			span = &spanlist[n++];
			span->x = x;
			state = INSIDE_SPAN;
		} else if (state == INSIDE_SPAN && *p == '-') {
			assert(span != NULL);
		} else if (state == INSIDE_SPAN && *p == ']') {
			span->w = (x+1) - span->x;
			state = OUTSIDE_SPAN;
			span = NULL;
		} else {
			assert(!"unhandled char/state");
		}
	}
	assert((state == OUTSIDE_SPAN) && "span was not terminated");
	assert(span == NULL);
	return n;
}

static void unparse_spanlist(char* str, int str_cap, int n_spans, struct span* spanlist)
{
	assert_valid_spanlist(n_spans, spanlist);
	memset(str, ' ', str_cap);
	char* end = str;
	for (int i = 0; i < n_spans; i++) {
		struct span* span = &spanlist[i];
		for (int i = 0; i < span->w; i++) {
			const char ch = (i == 0) ? '[' : (i == (span->w-1)) ? ']' : '-';
			const int x = span->x + i;
			assert((x+1) < str_cap);
			end = str+x+1;
			str[x] = ch;
		}
	}
	*end = 0;
}

static void test_spanlist_parser_and_unparser(void)
{
	struct span spanlist[1<<10];
	char unparsed[1<<10];
	{
		const int n = parse_spanlist(spanlist, ARRAY_LENGTH(spanlist), "");
		assert(n == 0);
		unparse_spanlist(unparsed, ARRAY_LENGTH(unparsed), n, spanlist);
		assert(strcmp("", unparsed) == 0);
	} {
		const int n = parse_spanlist(spanlist, ARRAY_LENGTH(spanlist), "              ");
		assert(n == 0);
		unparse_spanlist(unparsed, ARRAY_LENGTH(unparsed), n, spanlist);
		assert(strcmp("", unparsed) == 0);
	} {
		const int n = parse_spanlist(spanlist, ARRAY_LENGTH(spanlist), "[---]");
		//                                                              01234
		assert(n == 1);
		assert(spanlist[0].x == 0); assert(spanlist[0].w == 5);
		unparse_spanlist(unparsed, ARRAY_LENGTH(unparsed), n, spanlist);
		assert(strcmp("[---]", unparsed) == 0);
	} {
		const int n = parse_spanlist(spanlist, ARRAY_LENGTH(spanlist), "  [---]  ");
		//                                                              0123456
		assert(n == 1);
		assert(spanlist[0].x == 2); assert(spanlist[0].w == 5);
		unparse_spanlist(unparsed, ARRAY_LENGTH(unparsed), n, spanlist);
		assert(strcmp("  [---]", unparsed) == 0);
	} {
		const int n = parse_spanlist(spanlist, ARRAY_LENGTH(spanlist), "  [---]  [--] [-]");
		//                                                              01234567890123456
		assert(n == 3);
		assert(spanlist[0].x == 2);  assert(spanlist[0].w == 5);
		assert(spanlist[1].x == 9);  assert(spanlist[1].w == 4);
		assert(spanlist[2].x == 14); assert(spanlist[2].w == 3);
		unparse_spanlist(unparsed, ARRAY_LENGTH(unparsed), n, spanlist);
		assert(strcmp("  [---]  [--] [-]", unparsed) == 0);
	} {
		const int n = parse_spanlist(spanlist, ARRAY_LENGTH(spanlist), " [][] []");
		//                                                              01234567
		assert(n == 3);
		assert(spanlist[0].x == 1); assert(spanlist[0].w == 2);
		assert(spanlist[1].x == 3); assert(spanlist[1].w == 2);
		assert(spanlist[2].x == 6); assert(spanlist[2].w == 2);
		unparse_spanlist(unparsed, ARRAY_LENGTH(unparsed), n, spanlist);
		assert(strcmp(" [][] []", unparsed) == 0);
	}
}

static void test_span_overlaps(void)
{
	#define TEST_OVERLAP(expected_p,str_xs,str_ys) \
	{ \
		struct span xs[1<<10]; \
		struct span ys[1<<10]; \
		const int n_xs = parse_spanlist(xs, ARRAY_LENGTH(xs), str_xs); \
		assert((n_xs == 1) && "bad test: string should contain exactly one span"); \
		const int n_ys = parse_spanlist(ys, ARRAY_LENGTH(ys), str_ys); \
		assert((n_ys == 1) && "bad test: string should contain exactly one span"); \
		const int p = spans_overlap(xs[0], ys[0]); \
		if (p != expected_p) { \
			fprintf(stderr, \
				"Expected the following two spans to %s:\n" \
				" `%s`\n" \
				" `%s`\n" \
				"but they %s.\n", \
				expected_p ? "overlap" : "NOT overlap", \
				str_xs, str_ys, \
				p ? "DID overlap" : "did NOT overlap"); \
			abort(); \
		} \
	}

	TEST_OVERLAP(0,
	"[--]",
	"     [---]");

	TEST_OVERLAP(0,
	"      [--]",
	"[---]");

	TEST_OVERLAP(0,
	"[--]",
	"    [---]");

	TEST_OVERLAP(0,
	"     [--]",
	"[---]");

	TEST_OVERLAP(1,
	"[--]",
	"   [---]");

	TEST_OVERLAP(1,
	"      [--]",
	"  [---]");

	TEST_OVERLAP(1,
	"[-----]",
	"   [-----]");

	TEST_OVERLAP(1,
	"       [-----]",
	"   [-----]");

	TEST_OVERLAP(0,
	"          [-----]",
	"   [-----]");

	TEST_OVERLAP(1,
	"     [--]",
	"  [--------]");

	TEST_OVERLAP(1,
	"  [--------]",
	"     [--]");

	#undef TEST_OVERLAP
}

static inline int is_same_string_ignoring_trailing_spaces(const char* a, const char* b)
{
	size_t na = strlen(a);
	size_t nb = strlen(b);
	while (na > 0 && a[na-1] == ' ') na--;
	while (nb > 0 && b[nb-1] == ' ') nb--;
	if (na != nb) return 0;
	if (na == 0 && nb == 0) return 1;
	assert(na == nb);
	return memcmp(a, b, na) == 0;
}

// test a commutative binary operator; it asserts that A@B=C and B@A=C where
// "@" is an "operator" defined by "operator_fn"
#define TEST_COMMUTATIVE_BINOP(opname,operator_fn,str_aa,str_bb,str_expected_cc) \
{ \
	for (int pass = 0; pass < 2; pass++) { \
		struct span xs[1<<10]; \
		struct span ys[1<<10]; \
		struct span actual_cc[1<<10]; \
		const char* str_xs = (pass == 0) ? str_aa : (pass == 1) ? str_bb : NULL; \
		const char* str_ys = (pass == 0) ? str_bb : (pass == 1) ? str_aa : NULL; \
		const int n_xs = parse_spanlist(xs, ARRAY_LENGTH(xs), str_xs); \
		const int n_ys = parse_spanlist(ys, ARRAY_LENGTH(ys), str_ys); \
		const int n_actual_cc = operator_fn(actual_cc, ARRAY_LENGTH(actual_cc), n_xs, xs, n_ys, ys); \
		char str_actual_cc[1<<10]; \
		unparse_spanlist(str_actual_cc, sizeof(str_actual_cc), n_actual_cc, actual_cc ); \
		if (!is_same_string_ignoring_trailing_spaces(str_expected_cc, str_actual_cc)) { \
			fprintf(stderr, \
				"Expected the " opname " (pass=%d) between the following two inputs:\n" \
				" `%s`\n" \
				" `%s`\n" \
				"to yield:\n" \
				" `%s`\n" \
				"but instead we got:\n" \
				" `%s`\n", \
				pass, \
				str_xs, \
				str_ys, \
				str_expected_cc, \
				str_actual_cc); \
			abort(); \
		} \
		assert_clean_spanlist(n_actual_cc, actual_cc); \
	} \
}
#define TEST_UNION(str_aa,str_bb,str_expected_cc)        TEST_COMMUTATIVE_BINOP("UNION",spanlist_union,str_aa,str_bb,str_expected_cc)
#define TEST_INTERSECTION(str_aa,str_bb,str_expected_cc) TEST_COMMUTATIVE_BINOP("INTERSECTION",spanlist_intersection,str_aa,str_bb,str_expected_cc)

#define TEST_DIFFERENCE(str_xs,str_ys,str_expected_zs) \
{ \
	struct span xs[1<<10]; \
	struct span ys[1<<10]; \
	struct span actual_zs[1<<10]; \
	const int n_xs = parse_spanlist(xs, ARRAY_LENGTH(xs), str_xs); \
	const int n_ys = parse_spanlist(ys, ARRAY_LENGTH(ys), str_ys); \
	const int n_actual_zs = spanlist_difference(actual_zs, ARRAY_LENGTH(actual_zs), n_xs, xs, n_ys, ys); \
	char str_actual_zs[1<<10]; \
	unparse_spanlist(str_actual_zs, sizeof(str_actual_zs), n_actual_zs, actual_zs ); \
	if (!is_same_string_ignoring_trailing_spaces(str_expected_zs, str_actual_zs)) { \
		fprintf(stderr, \
			"Expected the DIFFERENCE between the following two inputs:\n" \
			" `%s`\n" \
			" `%s`\n" \
			"to yield:\n" \
			" `%s`\n" \
			"but instead we got:\n" \
			" `%s`\n", \
			str_xs, \
			str_ys, \
			str_expected_zs, \
			str_actual_zs); \
		abort(); \
	} \
	assert_clean_spanlist(n_actual_zs, actual_zs); \
}

static void test_binop_union(void)
{
	TEST_UNION("", "", "");

	TEST_UNION(
	" [--]                  [--]          [-----]     [--]  [----]    [--]   ",
	"      [--] [-]  [----]       [----]      [----]  [--]  [----]  [------] ",
	" [--] [--] [-]  [----] [--]  [----]  [--------]  [--]  [----]  [------] ");

	TEST_UNION(
	"  [------]   [-----]    [---]     [-]      [-]  [-] ",
	"    [--]       [-]    [----]   [---]   [---]        " ,
	"  [------]   [-----]  [-----]  [----]  [-----]  [-] ");

	// "auto-merge" tests
	TEST_UNION(
	" [---][----]  [---][----]  [---][----]",
	"              [---][----]  [------][-]",
	" [---------]  [---------]  [---------]");

	// sieve/merge tests
	TEST_UNION(
	"    [-]   [-]   [-]   [-]  [-]   [----]   [-]        [-]    [-]    [-] [-] [-][-]   [-]   [-]",
	" [-]   [-]   [-]   [-]        [-]      [-]   [-]  [-]   [-]   [-]   [-]    [-]   [-]   [-][-]",
	" [----------------------]  [-------------------]  [-------] [---]  [-----] [----------------]");

	TEST_UNION(
	" [-]    [-]   [-]     [-]       [----------------]          ",
	"      [------------------]  [-]   [-]   [-]   [-]   [-] [-] ",
	" [-]  [------------------]  [-] [----------------]  [-] [-] ");
}

static void test_binop_intersection(void)
{
	TEST_INTERSECTION("", "", "");

	TEST_INTERSECTION(
	"  [-----]     [-----]        [-----] [---]  [---] [---] [---] [---][---]       [---] [---]  [-] [-] [-] [-] ",
	"      [----]          [----]     [----]        [---] [---]               [][]  [---]            [-----]     ",
	"      [-]                        [-] []        [] [] [] []                     [---]            [-] [-]     ");

	TEST_INTERSECTION(
	"  [---]  [-]  [---]  [---] [---] [---]  [--][--][--] [-------------------] ",
	"        [---]           [---]             [------]      [-]   [-]   [-]    ",
	"         [-]            [] []             [------]      [-]   [-]   [-]    ");
}

static void test_binop_difference(void)
{
	TEST_DIFFERENCE("", "", "");

	TEST_DIFFERENCE(
	" [--]  [--]           [--] [--] [--]        [--]      ",
	"            [--] [--] [--] [--]     [--][--]    [--]  ",
	" [--]  [--]                     [--]        [--]      ");

	TEST_DIFFERENCE(
	"  [------]         [------] [---] [---] [-] [---]    [--]  [--]      [---------]        [--]    ",
	"      [------] [------]        [-------------]     [------------] [---]       [---] [---]  [--] ",
	"  [--]                 [--] [-]               [-]                      [-----]           []     ");

	TEST_DIFFERENCE(
	"   [--]  [---]   [-------------]     [-]   [-]   [-]     [-]   [-]   [-]    [---------]   ",
	" [--------]         [-]   [-]     [-------------------]  [-------------] [---]   []  [---]",
	"           [-]   [-]   [-]   [-]                                              [-]  []     ");

	TEST_DIFFERENCE(
	" [--]    [--][--]     ",
	"     [--]        [--] ",
	" [--]    [------]     ");
}

static void b2ss(char* dst, int n, int value)
{
	char* p = dst;
	for (int i = 0; i < n; i++) {
		char v = ((value >> i) & 1) ? '-' : ' ';
		*(p++) = v;
		*(p++) = v;
	}
	*p = 0;
	char pp = ' ';
	for (p = dst; *p; p++) {
		if (*p == '-' && pp == ' ') {
			*p = '[';
		}
		if (*p == '-' && (*(p+1) == 0 || *(p+1) == ' ')) {
			*p = ']';

		}
		pp = *p;
	}

	while (p > dst && *(p-1) == ' ') {
		p--;
		*p = 0;
	}
}

static inline double tim(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return 1e-9*(double)ts.tv_nsec + (double)ts.tv_sec;
}

static void test_binop_generic(void)
{
	const int L = 7;
	const int N = 1<<L;
	const double t0 = tim();
	for (int i0 = 0; i0 < N; i0++) {
		for (int i1 = 0; i1 < N; i1++) {
			char ss0[256], ss1[256], ss2[256];
			b2ss(ss0, L, i0);
			b2ss(ss1, L, i1);
			{
				b2ss(ss2, L, (i0 | i1));
				TEST_UNION(ss0, ss1, ss2);
				//printf("UNION of:\n`%s`\n`%s`\nis:\n`%s`\n", ss0, ss1, ss2);
			} {
				b2ss(ss2, L, (i0 & i1));
				TEST_INTERSECTION(ss0, ss1, ss2);
				//printf("INTERSECTION of:\n`%s`\n`%s`\nis:\n`%s`\n", ss0, ss1, ss2);
			} {
				b2ss(ss2, L, (i0 & ~i1));
				TEST_DIFFERENCE(ss0, ss1, ss2);
				//printf("DIFFERENCE of:\n`%s`\n`%s`\nis:\n`%s`\n", ss0, ss1, ss2);
			}
		}
	}

	const double dt = tim() - t0;
	const double n = (N*N*3);
	printf("binop generic: %f per second (n=%.0f, dt=%f)\n", n/dt, n, dt);
}

static struct rzr* getrz(float pixels_per_unit, int width, int height, int supersampling_factor)
{
	static struct rzr_tx tx[256];
	static struct rzr_op prg[256];
	static struct rzr rzr;
	rzr_init(&rzr, ARRAY_LENGTH(tx), tx, ARRAY_LENGTH(prg), prg, width, height, pixels_per_unit, supersampling_factor);
	return &rzr;
}

static struct rzr* getrz3x3(float pixels_per_unit, int width, int height, int supersampling_factor)
{
	static struct rzr_tx tx[256];
	static struct rzr_op prg[256];
	static struct rzr rzr;
	const int widths[]  = {width,  -1, width,  0};
	const int heights[] = {height, -1, height, 0};
	rzr_init_MxN(&rzr, ARRAY_LENGTH(tx), tx, ARRAY_LENGTH(prg), prg, widths, heights, pixels_per_unit, supersampling_factor);
	return &rzr;
}

static void test_point_queries(void)
{
	const int width = 512;
	const int height = 512;
	const int ss = 16;

	{
		struct rzr* rzr = getrz(width/2, width, height, ss);
		Circle(0.5);
		assert(rzr_query(rzr, width/2, height/2));
		assert(rzr_query(rzr, width/2 + 50, height/2 - 50));
		assert(!rzr_query(rzr, 0, 0));
		assert(!rzr_query(rzr, width-1, height-1));
		assert(!rzr_query(rzr, -100, -1000));
		assert(!rzr_query(rzr, width+5000, -500));
		assert(!rzr_query(rzr, width, 0));
		assert(!rzr_query(rzr, 0, height));

		assert(rzr_subpixel_query(rzr, (ss*width)/2, (ss*height)/2));
		assert(rzr_subpixel_query(rzr, (ss*width)/2 + 50, (ss*height)/2 - 100));
		assert(!rzr_subpixel_query(rzr, 0, 0));
		assert(!rzr_subpixel_query(rzr, width*ss, height*ss));
	}

	{
		struct rzr* rzr = getrz(width/2, width, height, ss);
		BeginPoly();
		Vertex(0, -1);
		Vertex(1, 1);
		Vertex(-1, 1);
		EndPoly();
		assert(!rzr_query(rzr, 0, 0));
		assert(!rzr_query(rzr, width, 0));
		assert(!rzr_query(rzr, 0, height/2));
		assert(!rzr_query(rzr, width, height/2));
		assert(rzr_query(rzr, width/2, height/2));
	}

	const int N = 10000;
	for (int i = 0; i < N; i++) {
		struct rzr* rzr = getrz(width/2, width, height, ss);
		Save();
		Rotate((float)(i*360) / (float)N);
		Star(6 + (i%18), 0.8f, 0.5f);
		Restore();
		assert(!rzr_query(rzr, 0, 0));
		assert(!rzr_query(rzr, width-1, 0));
		assert(rzr_query(rzr, width/2, height/2));
		assert(!rzr_query(rzr, width-1, height-1));
		assert(!rzr_query(rzr, 0, height-1));
		assert(!rzr_query(rzr, width-1, height-1));
	}
}


static void test_subrender(void)
{
	#define RTEST(STR_L0,STR_L1,STR_L2,STR_EXPECTED_PIXELS) \
	{ \
		const char* str_l0 = STR_L0; \
		const char* str_l1 = STR_L1; \
		const char* str_l2 = STR_L2; \
		struct span l0s[1<<10]; \
		struct span l1s[1<<10]; \
		struct span l2s[1<<10]; \
		const int n_l0s = parse_spanlist(l0s, ARRAY_LENGTH(l0s), str_l0); \
		const int n_l1s = parse_spanlist(l1s, ARRAY_LENGTH(l1s), str_l1); \
		const int n_l2s = parse_spanlist(l2s, ARRAY_LENGTH(l2s), str_l2); \
		const char* str_expected_pixels = STR_EXPECTED_PIXELS; \
		const int nx = strlen(str_expected_pixels); \
		assert(((nx%3) == 0) && "bad string for expected pixels; must be multiple-of-3"); \
		const int width = nx/3; \
		uint8_t expected_pixels[1<<10]; \
		const char* rp = str_expected_pixels; \
		for (int i0 = 0; i0 < width; i0++) { \
			const char ch = *(rp++); \
			assert(('0' <= ch && ch <= '9') && "bad string for expected pixels; expected digit"); \
			expected_pixels[i0] = (((int)ch - (int)'0') * 255) / 9; \
			for (int i1 = 0; i1 < 2; i1++) assert((*(rp++) == '.') && "bad string for expected pixels; must be padded with dots"); \
		} \
		uint8_t pixels[1<<10]; \
		memset(pixels, 0xfe, width); \
		struct spanline sub_scanlines[3] = { \
			{.n = n_l0s, .spans = l0s}, \
			{.n = n_l1s, .spans = l1s}, \
			{.n = n_l2s, .spans = l2s}, \
		}; \
		subrender(pixels, 1, &width, 3, sub_scanlines); \
		if (memcmp(pixels, expected_pixels, width) != 0) { \
			fprintf(stderr, "in subrender() of:\n%s\n%s\n%s\nexpected pixels:\n%s\n", str_l0, str_l1, str_l2, str_expected_pixels); \
			for (int i0 = 0; i0 < width; i0++) fprintf(stderr, "%.2x.", expected_pixels[i0]); \
			fprintf(stderr, "\nbut got:\n"); \
			for (int i0 = 0; i0 < width; i0++) fprintf(stderr, "%c..", '0' + (((int)(pixels[i0])*9+128)/255)); \
			fprintf(stderr, "\n"); \
			for (int i0 = 0; i0 < width; i0++) fprintf(stderr, "%.2x.", pixels[i0]); \
			fprintf(stderr, "\n"); \
			abort(); \
		} \
	}

	RTEST(
		"   [-------]     [-------]     [-------]     [-]     [-]     [-]     ",
		"   [-------]     [-------]     [-------]     [-]     [-]     [-]     ",
		"   [-------]     [-------]     [-------]     [-]     [-]     [-]     ",
		"0..9..9..9..0..3..9..9..6..0..6..9..9..3..0..9..0..3..6..0..6..3..0..");

	RTEST(
		"     [-------]     [-------]     [-------]     [-]     [-]     [-]   ",
		"    [-------]     [-------]     [-------]     [-]     [-]     [-]    ",
		"   [-------]     [-------]     [-------]     [-]     [-]     [-]     ",
		"0..6..9..9..3..1..8..9..8..1..3..9..9..6..0..6..3..1..7..1..3..6..0..");

	RTEST(
		"       [-------------]        [---------]               [------]  ",
		"     [-------------]       [---------]              [--------]    ",
		"   [-------------]      [---------]            [---------]        ",
		"0..4..8..9..9..9..5..1..3..6..9..8..5..2..0..1..3..5..7..7..5..1.."
	);

	RTEST(
		"            []  [] []  [-] []        [---] []  []                    ",
		"          [-]  []    [-]    []        []   [--] []                   ",
		"            [-------------]          []   [-]  [-]                   ",
		"0..0..0..2..6..7..5..7..5..4..0..0..5..4..7..4..5..0..0..0..0..0..0.."
	);

	#undef RTEST
}

static void bitmap_ascii_dump(int width, int height, uint8_t* pixels)
{
	printf("\n");
	for (int y = 0; y < height; y++) {
		{
			uint8_t* p = pixels + y*width;
			for (int x = 0; x < width; x++) {
				const int v = *(p++);
				printf("%.2x", v);
			}
		}
		printf("   ");
		{
			uint8_t* p = pixels + y*width;
			for (int x = 0; x < width; x++) {
				const int v = *(p++);
				const char c = " .:ioVM@"[v>>5];
				printf("%c%c", c, c);
			}
		}
		printf("\n");
	}
}

static void test_rzr(void)
{
	const int width = 32;
	const int height = 32;
	struct rzr* rzr = getrz(width/2, width, height, 16);
	Circle(1.0f);
	Circle(0.9f);
	Difference();
	Star(8, 0.9f, 0.3f);
	Union();
	uint8_t scratch[1<<18];
	const size_t scratch_sz = sizeof(scratch);
	uint8_t pixels[width*height];
	memset(pixels, 0xfe, sizeof pixels);
	rzr_render(rzr, scratch_sz, scratch, width, pixels);
	bitmap_ascii_dump(width, height, pixels);
}

static void test_3x3(void)
{
	const int tile_size = 5;
	struct rzr* rzr = getrz3x3(tile_size, tile_size, tile_size, 16);
	Circle(1.0f);
	Circle(0.45f);
	Difference();
	uint8_t scratch[1<<18];
	const size_t scratch_sz = sizeof(scratch);
	uint8_t pixels[rzr->bitmap_width * rzr->bitmap_height];
	memset(pixels, 0xfe, sizeof pixels);
	rzr_render(rzr, scratch_sz, scratch, rzr->bitmap_width, pixels);

	const char* msg = "(-: STREEEETCH :-)";

	const int nch = strlen(msg);
	const int inner_width = (nch+1)/2 + 2;
	const int inner_height = 3;
	const int total_width = tile_size + inner_width + tile_size;
	const int total_height = tile_size + inner_height + tile_size;

	#if 0
	bitmap_ascii_dump(rzr->bitmap_width, rzr->bitmap_height, pixels);
	#endif

	//const int msg_x0 = (total_width/2 - (nch+1)/2)*2;
	const int msg_x0 = total_width/2 - nch/4;
	const int msg_y = total_height/2;
	printf("\n");
	for (int y = 0; y < total_height; y++) {
		for (int x = 0; x < total_width; x++) {
			//uint8_t* p = pixels + y*width;
			const int qx =
				  x < tile_size               ? x
				: x < (total_width-tile_size) ? tile_size
				: x - total_width + tile_size*2+1;
			const int qy =
				  y < tile_size                ? y
				: y < (total_height-tile_size) ? tile_size
				: y - total_height + tile_size*2+1;
			assert(0 <= qx && qx < rzr->bitmap_width);
			assert(0 <= qy && qy < rzr->bitmap_height);
			if (y == msg_y && (msg_x0 <= x && x < (msg_x0+(nch+1)/2))) {
				for (int i = 0; i < 2; i++) {
					const int ii = 2*(x-msg_x0)+i;
					printf("%c", ii < nch ? msg[ii] : ' ');
				}
			} else {
				uint8_t* p = pixels + qx + qy*rzr->bitmap_width;
				const int v = *p;
				const char c = " .:ioVM@"[v>>5];
				printf("%c%c", c, c);
			}
		}
		printf("\n");
	}
}

int main(int argc, char** argv)
{
	test_spanlist_parser_and_unparser();
	test_span_overlaps();
	test_binop_union();
	test_binop_intersection();
	test_binop_difference();
	test_binop_generic();
	test_point_queries();
	test_subrender();
	test_rzr();
	test_3x3();
	printf("TESTS: OK\n");
	return EXIT_SUCCESS;
}

#endif//UNIT_TEST
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
#ifdef DEMOS
// cc -DDEMOS -Wall -O0 -g rzr.c -o demo_rzr -lm && ./demo_rzr
//   or
// cc -DDEMOS -Wall -O2 rzr.c -o demo_rzr -lm && ./demo_rzr

// Choose one of the following:
#define TEST_PERFORMANCE // uses static compile-time allocations; timings are not affected by malloc()
//#define TEST_MEMERROR // deliberately uses fresh malloc()s for each render to make valgrind better at finding problems

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define MSF_GIF_IMPL
#include "msf_gif.h"

static inline double tim(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return 1e-9*(double)ts.tv_nsec + (double)ts.tv_sec;
}

#define N_TX (1<<10)
#define N_PRG (1<<10)
#define SCRATCH_SZ (1<<24)

static struct {
	struct rzr rzr;
	#if defined(TEST_MEMERROR)
	struct rzr_tx* tx;
	struct rzr_op* prg;
	uint8_t* scratch;
	#elif defined(TEST_PERFORMANCE)
	struct rzr_tx tx[N_TX];
	struct rzr_op prg[N_PRG];
	uint8_t scratch[SCRATCH_SZ];
	#else
	#error "TEST_MEMERROR or TEST_PERFORMANCE must be defined"
	#endif

	double t0;

	// anim
	int anim_n_frames, anim_n_frames_remaining;
	int anim_width, anim_height;
	float anim_pixels_per_unit;
	int anim_supersampling_factor;
	uint8_t* anim_data;
	uint8_t* anim_p;

	// tile
	int tile_width, tile_height, tile_stride;
	int tile_n_columns, tile_n_rows;
	int tile_cursor;
	uint8_t* tile_bitmap;
} g;

static void begin_mem(void)
{
	#if defined(TEST_MEMERROR)
	assert(g.tx == NULL);
	assert(g.prg == NULL);
	assert(g.scratch == NULL);
	// deliberate use of malloc() instead of calloc(); calloc() clears
	// allocations which makes it less likely that valgrind can detect
	// usage of uninitialized values
	g.tx = malloc(N_TX * sizeof(g.tx[0]));
	g.prg = malloc(N_PRG * sizeof(g.prg[0]));
	g.scratch = malloc(SCRATCH_SZ);
	#elif defined(TEST_PERFORMANCE)
	// memory is static
	#else
	#error "TEST_MEMERROR or TEST_PERFORMANCE must be defined"
	#endif
}

static void end_mem(void)
{
	#if defined(TEST_MEMERROR)
	assert(g.scratch != NULL); free(g.scratch); g.scratch = NULL;
	assert(g.prg     != NULL); free(g.prg);     g.prg     = NULL;
	assert(g.tx      != NULL); free(g.tx);      g.tx      = NULL;
	#elif defined(TEST_PERFORMANCE)
	// memory is static
	#else
	#error "TEST_MEMERROR or TEST_PERFORMANCE must be defined"
	#endif
}

static void post_report(struct rzr* rzr)
{
	#if 0
	printf("render used tx:%d prg:%d scratch:%zd\n", rzr->tx_stack_height_max, rzr->prg_length_max, rzr->scratch_alloc_max);
	#endif
}

static void begin_tiles(int width, int height, int n_columns, int n_rows)
{
	g.tile_width = width;
	g.tile_height = height;
	g.tile_stride = width * n_columns;
	g.tile_n_columns = n_columns;
	g.tile_n_rows = n_rows;
	g.tile_cursor = 0;
	g.tile_bitmap = calloc(width*n_columns*height*n_rows, 1);
	g.t0 = tim();
}

static struct rzr* begin_tile(float pixels_per_unit, int supersampling_factor)
{
	rzr_init(&g.rzr, N_TX, g.tx, N_PRG, g.prg, g.tile_width, g.tile_height, pixels_per_unit, supersampling_factor);
	return &g.rzr;
}

static void end_tile(void)
{
	const int column = g.tile_cursor % g.tile_n_columns;
	const int row    = g.tile_cursor / g.tile_n_columns;
	assert(0 <= column && column < g.tile_n_columns);
	assert(0 <= row && row < g.tile_n_rows);
	uint8_t* p =
		g.tile_bitmap
		+ (column * g.tile_width)
		+ (row * g.tile_height) * g.tile_stride;
	rzr_render(&g.rzr, SCRATCH_SZ, g.scratch, g.tile_stride, p);
	g.tile_cursor++;
}

static void end_tiles(const char* path)
{
	const double dt = tim() - g.t0;
	if (path != NULL) {
		assert(stbi_write_png(path, g.tile_stride, g.tile_height*g.tile_n_rows, 1, g.tile_bitmap, g.tile_stride));
		printf("Wrote %s; ", path);
	}
	printf("render took %.5fs\n", dt);
	free(g.tile_bitmap);
}

static struct rzr* begin_render(int width, int height, float pixels_per_unit, int supersampling_factor)
{
	begin_mem();
	rzr_init(&g.rzr, N_TX, g.tx, N_PRG, g.prg, width, height, pixels_per_unit, supersampling_factor);
	return &g.rzr;
}

static void end_render(const char* path)
{
	const int w = g.rzr.bitmap_width;
	const int h = g.rzr.bitmap_height;
	uint8_t* pixels = malloc(w*h);
	const double t0 = tim();
	rzr_render(&g.rzr, SCRATCH_SZ, g.scratch, w, pixels);
	post_report(&g.rzr);
	const double dt = tim() - t0;
	if (path != NULL) {
		assert(stbi_write_png(path, w, h, 1, pixels, w));
		printf("Wrote %s; ", path);
	}
	printf("render took %.5fs\n", dt);
	end_mem();
}

static void begin_anim(int n_frames, int width, int height, float pixels_per_unit, int supersampling_factor)
{
	g.anim_n_frames = n_frames;
	g.anim_n_frames_remaining = g.anim_n_frames;
	g.anim_width = width;
	g.anim_height = height;
	g.anim_pixels_per_unit = pixels_per_unit;
	g.anim_supersampling_factor =  supersampling_factor;
	g.anim_data = malloc(n_frames*width*height);
	g.anim_p = g.anim_data;
	g.t0 = tim();
}

static struct rzr* begin_frame(void)
{
	assert(g.anim_n_frames_remaining > 0);
	begin_mem();
	rzr_init(&g.rzr, N_TX, g.tx, N_PRG, g.prg, g.anim_width, g.anim_height, g.anim_pixels_per_unit, g.anim_supersampling_factor);
	return &g.rzr;
}

static void end_frame(void)
{
	assert(g.anim_n_frames_remaining > 0);
	g.anim_n_frames_remaining--;
	rzr_render(&g.rzr, SCRATCH_SZ, g.scratch, g.anim_width, g.anim_p);
	post_report(&g.rzr);
	#if 0
	{
		// store individual frames as png
		char ppath[1<<10];
		snprintf(ppath, sizeof ppath, "_frame_%.4d.png", g.anim_n_frames - g.anim_n_frames_remaining);
		assert(stbi_write_png(ppath, g.anim_width, g.anim_height, 1, g.anim_p, g.anim_width));
		printf("Wrote %s\n", ppath);
	}
	#endif
	end_mem();
	g.anim_p += g.anim_width * g.anim_height;
}

static void end_anim(const char* path)
{
	assert(g.anim_n_frames_remaining == 0);
	const double dt = tim() - g.t0;

	uint8_t* rgba = malloc(4*g.anim_width*g.anim_height);
	MsfGifState gs = {};
	msf_gif_begin(&gs, g.anim_width, g.anim_height);
	for (int i = 0; i < g.anim_n_frames; i++) {
		uint8_t* data = g.anim_data + i*g.anim_width*g.anim_height;
		uint8_t* rp = data;
		uint8_t* wp = rgba;
		for (int y = 0; y < g.anim_height; y++) {
			for (int x = 0; x < g.anim_width; x++) {
				uint8_t v = *(rp++);
				*(wp++) = v;
				*(wp++) = v;
				*(wp++) = v;
				*(wp++) = 255;
			}
		}
		msf_gif_frame(&gs, rgba, 2, 8, g.anim_width * 4);
	}
	MsfGifResult result = msf_gif_end(&gs);
	if (result.data) {
		FILE* fp = fopen(path, "wb");
		assert(fp != NULL);
		assert(fwrite(result.data, result.dataSize, 1, fp) == 1);
		assert(fclose(fp) == 0);
	}
	msf_gif_free(result);
	free(rgba);
	free(g.anim_data);

	printf("animation of %d frames took %.4fs; wrote %s\n", g.anim_n_frames, dt, path);
}

int main(int argc, char** argv)
{
	#if 1
	{
		const int S = 512;
		struct rzr* rzr = begin_render(S, S, S/2, 16);

		const int N = 30;

		const float inc = 3.0f * 6.283185307179586f / (float)N;
		for (int i = 0; i < N; i++) {
			const float R0 = 1.0f - (float)i/(float)N;
			const float R1 = R0 - 0.5f/(float)N;

			Circle(R0);
			if (i > 0) Union();
			Save();
			const float t = R0-R1;
			const float tx = t*cosf(inc*(float)i);
			const float ty = t*sinf(inc*(float)i);
			Translate(tx, ty);
			Circle(R1);
			Restore();
			Difference();

		}
		end_render("_rzrdemo_sharp.png");
	}
	#endif

	#if 1
	{
		const int S = 128;
		begin_tiles(S, S, 4, 2);
		{
			struct rzr* rzr = begin_tile(S/2, 16);
			Star(40, 1.0f, 0.2f);
			end_tile();
		}
		{
			struct rzr* rzr = begin_tile(S/2, 16);
			Star(10, 1.0f, 0.5f);
			Circle(0.4f);
			Difference();
			end_tile();
		}
		{
			struct rzr* rzr = begin_tile(S/2, 16);
			const int N = 10;
			for (int i = 0; i < N; i++) {
				Circle(1.0f - (float)i/(float)N);
				if ((i&1) == 0 && i >= 2) Union();
				if ((i&1) == 1) Difference();
			}
			end_tile();
		}
		{
			struct rzr* rzr = begin_tile(S/2, 16);
			const int N = 40;
			for (int i = 0; i < N; i++) {
				Circle(1.0f - (float)i/(float)N);
				if ((i&1) == 0 && i >= 2) Union();
				if ((i&1) == 1) Difference();
			}
			Star(40, 0.9f, 0.4f);
			Intersection();
			end_tile();
		}
		{
			struct rzr* rzr = begin_tile(S/2, 16);
			Circle(1.0);
			Save();
			Translate(-0.1,0.2);
			Rotate(33);
			IsoscelesTriangle(0.2, 0.5);
			Difference();
			Restore();
			Save();
			Rotate(160);
			IsoscelesTrapezoid(0.1,0.4,0.8);
			Difference();
			Restore();
			Save();
			Translate(0.5,0.2);
			Rotate(22);
			Box(0.3, 0.1);
			Difference();
			Restore();
			end_tile();
		}
		{
			struct rzr* rzr = begin_tile(S/2, 16);
			Save();
			Rotate(11);
			RoundedBox(0.8,0.6,0.3);
			Rotate(-11);
			RoundedBox(0.4,0.4,0.15);
			Difference();
			const float d = 0.03f;
			RoundedBox(0.4-d,0.4-d,0.15-d);
			Union();
			RoundedBox(0.4-2*d,0.4-2*d,0.15-2*d);
			Difference();
			Restore();
			end_tile();
		}
		#if 0
		{
			// TODO/FIXME
			struct rzr* rzr = begin_tile(S/2, 16);
			Save();
			//Scale(1.1);
			Rotate(10.0);
			Pattern(0.02, 0.03, 0.07, 0.1);
			Restore();
			end_tile();
		}
		#endif
		end_tiles("_rzrdemo_zoo.png");
	}
	#endif

	#if 1
	{
		const int NFRAM = 200;
		const int S = 128;
		begin_anim(NFRAM, S, S, S/2, 16);
		for (int i0 = 0; i0 < NFRAM; i0++) {
			const float d = (float)i0 / (float)NFRAM;
			struct rzr* rzr = begin_frame();
			const int N = 3;
			const float t0 = d * 6.283185307179586f;

			BeginPoly();
			for (int i1 = 0; i1 < N; i1++) {
				const float t1 = (((float)i1 + d) * 6.283185307179586f) / (float)N;
				const float x = cosf(t1);
				const float y = sinf(t1);
				Vertex(x,y);
			}
			EndPoly();

			const float r = 0.4f + sinf(t0)*0.35f;

			Circle(r+0.1f);
			Difference();

			Circle(r+0.05f+0.01f);
			Union();
			Circle(r+0.05f-0.01f);
			Difference();

			Circle(r);
			Union();

			BeginPoly();
			const int N2 = 5;
			for (int i1 = 0; i1 < N2; i1++) {
				const float phi = (((float)i1 - d*2.0f) * 6.283185307179586f) / (float)N2;
				const float x = r * cosf(phi);
				const float y = r * sinf(phi);
				Vertex(x,y);
			}
			EndPoly();
			Difference();

			end_frame();
		}
		end_anim("_rzrdemo_trirot.gif");
	}
	#endif

	return EXIT_SUCCESS;
}

#endif//DEMOS
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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

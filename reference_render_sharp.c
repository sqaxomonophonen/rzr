/*
A small program used to verify the "visual correctness" of rzr.c by rendering
an image equivalent to _rzrdemo_sharp.png (demo in rzr.c), but using a much
smaller and slower brute-force implementation that is easier to reason about.

Run with:
$ cc -Wall -O3 reference_render_sharp.c -o reference_render_sharp -lm && ./reference_render_sharp
Then compare _rzref_sharp.png with _rzrdemo_sharp.png (see rzr.c on how to
render the latter).
*/

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

static inline int circle_test(double cx, double cy, double r, double px, double py)
{
	const double x = px-cx;
	const double y = py-cy;
	const double d2 = x*x + y*y;
	const double r2 = r*r;
	return d2<r2;
}

int main(int argc, char** argv)
{
	const int S = 512;
	const int SS = 16;
	uint8_t* pixels = malloc(S*S);
	uint8_t* p = pixels;
	const double ns = 255.0 / (double)(SS*SS);
	const double scale0 = 2.0 / (double)S;
	const double scale1 = 2.0 / (double)(S*SS);
	printf("Rendering... might take a while!\n");
	for (int y = 0; y < S; y++) {
		for (int x = 0; x < S; x++) {
			int n = 0;
			#if 0
			if (y == 0) printf("pixel\n");
			#endif
			for (int sy = 0; sy < SS; sy++) {
				for (int sx = 0; sx < SS; sx++) {
					const double px = -1.0 + (double)x*scale0 + (double)sx*scale1;
					const double py = -1.0 + (double)y*scale0 + (double)sy*scale1;
					#if 0
					if (y == 0 && sy == 0) printf(" px=%f\n", px);
					#endif
					int p = 0;
					const int N = 30;
					const double inc = 3.0 * 6.283185307179586 / (double)N;
					for (int i = 0; i < N; i++) {
						const double r0 = 1.0 - (double)i/(double)N;
						const double r1 = r0 - 0.5/(double)N;
						p |= circle_test(0.0, 0.0, r0, px, py);
						const double t = r0-r1;
						p &= !circle_test(t*cos(inc*(double)i), t*sin(inc*(double)i), r1, px, py);
					}
					n += p;
				}
			}
			*(p++) = round((double)n * ns);
		}
	}

	const char* path = "_rzref_sharp.png";
	assert(stbi_write_png(path, S, S, 1, pixels, S));
	printf("Wrote %s\n", path);
	return 0;
}

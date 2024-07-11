// $ ./rzr_doodle.py doodle_demo.c
// If everything works, you should see the result of the rzr_doodle() render
// function below in a new window, and it should refresh and update the image
// when you change this file and save it. Try e.g. changing "N".

#include "rzr.h"

#if 0
// For testing `--ignore-unresolved-symbols` which adds
// `-Wl,--unresolved-symbols=ignore-in-object-files` as linker argument via
// `cc`. This is dangerous (the `undefined_symbol()` crashes), but also useful
// if your source file isn't "clean" and has cascading dependencies.
void undefined_symbol(void);
void exported_symbol(void)
{
	undefined_symbol();
}
#endif

void rzr_doodle(struct rzr* rzr)
{
	const int N = 10;
	for (int i = 0; i < N; i++) {
		Circle((float)(N-i) / (float)N);
		if (i == 0) continue;
		if ((i&1)) {
			Difference();
		} else {
			Union();
		}
	}
}

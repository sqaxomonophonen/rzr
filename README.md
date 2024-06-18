## Features
 - Small codebase, easy to drop-in: just copy `rzr.c` and `rzr.h` into your project.
 - Monochrome-only, 8-bit grayscale output
 - Stack-based drawing:
   - Shapes (pushes 1): circles, polygons (concave), stars, ...
   - Boolean operations (pops 2, pushes 1): union, difference and intersection.
 - Transforms: translate, rotate, scale
 - Reasonably fast? Merges subpixel-scanline span lists into `memset()`'able pixel spans.
 - Doesn't `malloc()`; you provide the memory it uses. Depends on libm/`math.h`.

## Examples

![](./_rzrdemo_zoo.png)
![](./_rzrdemo_trirot.gif)
![](./_rzrdemo_sharp.png)

(Generated with `cc -DDEMOS -Wall -O0 -g rzr.c -o demo_rzr -lm && ./demo_rzr`)

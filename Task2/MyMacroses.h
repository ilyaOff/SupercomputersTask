#pragma once

#define CENTR(w, i, j, M, N) (w[i + j * (M)])
#define LEFT(w, i, j, M, N) (w[i - 1 + j * (M)])
#define RIGHT(w, i, j, M, N) (w[i + 1 + j * (M)])
#define TOP(w, i, j, M, N) ( w[i + (j+1) * (M)])
#define BOTTOM(w, i, j, M, N) ( w[i  + (j-1) * (M)])


#define MainFunctionParallel2(res,f, w, i, j, M, N, a, b, h1, h2) {\
double center1, left1, right1, top1, down1;\
center1 = CENTR(w, i, j, M, N);\
left1 = LEFT(w, i, j, M, N);\
right1 = RIGHT(w, i, j, M, N);\
top1 = TOP(w, i, j, M, N);\
down1 = BOTTOM(w, i, j, M, N);\
double dx = ((RIGHT(a, i, j, M, N)) * (right1 - center1) - (CENTR(a, i, j, M, N)) * (center1 - left1)) / (h1 * h1);\
double dy = ((TOP(b, i, j, M, N)) * (top1 - center1) -(CENTR(b, i, j, M, N)) * (center1 - down1)) / (h2 * h2);\
res = f -(dx + dy);};

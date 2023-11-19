#pragma once

#define CENTR(w, i, indexJ, M) (w[i + indexJ ])
#define LEFT(w, i, j, M) (w[i - 1 + indexJ])
#define RIGHT(w, i, j, M) (w[i + 1 + indexJ])
#define TOP(w, i, j, M) ( w[i + indexJ + M])
#define BOTTOM(w, i, j, M) ( w[i + indexJ - M])


#define MainFunctionParallel2(res,f, w, i, indexJ, M, N, a, b, h1, h2) {\
double center1, left1, right1, top1, down1;\
center1 = CENTR(w, i, indexJ, M);\
left1 = LEFT(w, i, indexJ, M);\
right1 = RIGHT(w, i, indexJ, M);\
top1 = TOP(w, i, indexJ, M);\
down1 = BOTTOM(w, i, indexJ, M);\
double dx = ((RIGHT(a, i, indexJ, M)) * (right1 - center1) - (CENTR(a, i, indexJ, M)) * (center1 - left1)) / (h1 * h1);\
double dy = ((TOP(b, i, indexJ, M)) * (top1 - center1) -(CENTR(b, i, indexJ, M)) * (center1 - down1)) / (h2 * h2);\
res = f -(dx + dy);};

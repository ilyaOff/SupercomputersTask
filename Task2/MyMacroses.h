#pragma once

#define CENTR(w, i, j, M, N) (w[i][j])
#define LEFT(w, i, j, M, N) (i - 1 == 0 ? 0 : w[i - 1][j])
#define RIGHT(w, i, j, M, N) (i + 1 == M ? 0 : w[i + 1][j])
#define TOP(w, i, j, M, N) (j + 1 == N ? 0 : w[i][j + 1])
#define BOTTOM(w, i, j, M, N) (j - 1 == 0 ? 0 : w[i][j - 1])


#define MainFunctionParallel(w, k) ( -(k) * (w))

//#define MainFunctionParallel2(res) {\
//double center1, left1, right1, top1, down1;\
//center1 = 1;\
//left1 = i < 15 ? 2 : 5;\
// res = -( w[i][j] + center1 + left1);};

#define MainFunctionParallel2(res,f, w, i, j, M, N, a, b, h1, h2) {\
	double center1, left1, right1, top1, down1;\
	center1 = w[i][j];\
	left1 = i - 1 == 0 ? 0 : w[i - 1][j];\
	right1 = i + 1 == M ? 0 : w[i + 1][j];\
	top1 = j + 1 == N ? 0 : w[i][j + 1];\
	down1 = j - 1 == 0 ? 0 : w[i][j - 1];\
	double dx = (a[i + 1][j] * (right1 - center1) - a[i][j] * (center1 - left1)) / (h1 * h1);\
	double dy = (b[i][j + 1] * (top1 - center1) - b[i][j] * (center1 - down1)) / (h2 * h2);\
	res = f -(dx + dy);};

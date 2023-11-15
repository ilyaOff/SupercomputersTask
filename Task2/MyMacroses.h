#pragma once

#define CENTR(w, i, j, M, N) (w[i][j])
#define LEFT(w, i, j, M, N) (i - 1 == 0 ? 0 : w[i - 1][j])
#define RIGHT(w, i, j, M, N) (i + 1 == M ? 0 : w[i + 1][j])
#define TOP(w, i, j, M, N) (j + 1 == N ? 0 : w[i][j + 1])
#define BOTTOM(w, i, j, M, N) (j - 1 == 0 ? 0 : w[i][j - 1])
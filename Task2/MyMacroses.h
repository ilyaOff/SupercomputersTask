#pragma once

#define MainFunctionParallel2(res,f, w, i, j, M, N, a, b, h1, h2) {\
	double center1, left1, right1, top1, down1;\
	center1 = w[i][j];\
	left1 =  w[i - 1][j];\
	right1 =  w[i + 1][j];\
	top1 =  w[i][j + 1];\
	down1 =  w[i][j - 1];\
	double dx = (a[i + 1][j] * (right1 - center1) - a[i][j] * (center1 - left1)) / (h1 * h1);\
	double dy = (b[i][j + 1] * (top1 - center1) - b[i][j] * (center1 - down1)) / (h2 * h2);\
	res = f -(dx + dy);};

#define Write2FileWithStep(name, rank, step, w, sizeX, sizeY) {\
		std::ostringstream oss;\
		oss << name << step << "_" << rank << ".txt"; \
		ofstream fout(oss.str());\
		SaveResults(w, sizeY, sizeX, fout);\
		fout.close();\
};

#define Write2File(name, rank, w, sizeX, sizeY) {\
		std::ostringstream oss;\
		oss << name << rank << ".txt"; \
		ofstream fout(oss.str());\
		SaveResults(w, sizeY, sizeX, fout);\
		fout.close();\
};
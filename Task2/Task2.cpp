#include <iostream>
#include <string> 
#include <sstream>
#include<fstream>
#include<omp.h>
#include "Point.h"
#include "Task2.h"
//#define WRITEFILE
#define SHOWINFO
#define SHOWDELTAGRAPHIC

using namespace std;
#define MainFunction2(res,f, w, i, j, M, N, a, b, h1, h2) {\
	double center1, left1, right1, top1, down1;\
	center1 = w[i][j];\
	left1 = i - 1 == 0 ? 0 : w[i - 1][j];\
	right1 = i + 1 == M ? 0 : w[i + 1][j];\
	top1 = j + 1 == N ? 0 : w[i][j + 1];\
	down1 = j - 1 == 0 ? 0 : w[i][j - 1];\
	double dx = (a[i + 1][j] * (right1 - center1) - a[i][j] * (center1 - left1)) / (h1 * h1);\
	double dy = (b[i][j + 1] * (top1 - center1) - b[i][j] * (center1 - down1)) / (h2 * h2);\
	res = f -(dx + dy);};

//Параметры области
const Point A = Point(0.0, 0.0);
const Point B = Point(3.0, 0.0);
const Point C = Point(2.0, 3.0);
const Point D = Point(0.0, 3.0);
const Point P0 = Point(0.0, 0.0);
const Point P1 = Point(3.0, 3.0);

const double DELTA = 0.000001;
double epsilon = 0.01;
const int KMAX = 100000000;

int N = 0;
int M = 0;
int TracingPeriod = 10000;


int main(int argc, char **argv)
{
	ofstream log("f/Log.txt");
	if (argc == 1)
		log << "no arguments!" << endl;
	else
		ReadParameters(argc, argv);

	log << "M = " << M << " N = " << N << endl;
	if (M <= 0 || N <= 0)
	{
		log << "invalid parametres (M, N)";
		log.close();
		return -1;
	}

	#ifdef  SHOWDELTAGRAPHIC
	ofstream deltaLog("f/DeltaLog.txt");
	#endif //  SHOWDELTAGRAPHIC

	double start = omp_get_wtime();

	double h1 = (P1.X - P0.X) / (M);
	double h2 = (P1.Y - P0.Y) / (N);
	epsilon = h1;
	if (h2 > h1)
		epsilon = h2;
	epsilon = epsilon * epsilon;

	int sizeX = M + 1;
	int sizeY = N + 1;
	double **w = new double *[sizeX];
	double **r = new double *[sizeX];
	double **a = new double *[sizeX];
	double **b = new double *[sizeX];
	double **F = new double *[sizeX];

	for (int i = 0; i < sizeX; ++i)
	{
		w[i] = new double[sizeY];
		r[i] = new double[sizeY];
		a[i] = new double[sizeY];
		b[i] = new double[sizeY];
		F[i] = new double[sizeY];
	}

	for (int i = 0; i < sizeX; ++i)
	{
		for (int j = 0; j < sizeY; ++j)
		{
			w[i][j] = 0;
			r[i][j] = 0;
			double x = P0.X + i * h1;
			double y = P0.Y + j * h2;
			a[i][j] = CalculateA(x, y, h1, h2);
			b[i][j] = CalculateB(x, y, h1, h2);
			F[i][j] = CalculateF(x, y, h1, h2);
		}
	}
	w[sizeX / 2][sizeY / 2] = 1;

	double tau = 0.0;
	double rA = 0.0;
	double tauNumerator = 0.0, tauDenominator = 0.0;
	double deltaSqr2 = 0.0, deltaSqr1 = 0.0, deltaSqr = 0.0;
	int k = 1, stopEquals = 2 * TracingPeriod;

	#ifdef SHOWINFO
	{
		ofstream fout("f/F.txt");
		SaveResults(F, sizeX, sizeY, fout);
		fout.close();
	}
	{
		ofstream fout("f/A.txt");
		SaveResults(a, sizeX, sizeY, fout);
		fout.close();
	}
	{
		ofstream fout("f/B.txt");
		SaveResults(b, sizeX, sizeY, fout);
		fout.close();
	}
	#endif // SHOWINFO

	for (; k < KMAX; ++k)
	{
		//посчитать невязку r
		for (int i = 1; i < M; ++i)
		{
			for (int j = 1; j < N; ++j)
			{
				MainFunction2(r[i][j], -F[i][j], w, i, j, M, N, a, b, h1, h2);
				//r[i][j] = MainFunction(w, i, j, M, N, a, b, h1, h2) - F[i][j];
			}
		}
		//посчитать итерационный параметр
		tauNumerator = 0.0, tauDenominator = 0.0;
		for (int i = 1; i < M; ++i)
		{
			for (int j = 1; j < N; ++j)
			{
				MainFunction2(rA, 0, r, i, j, M, N, a, b, h1, h2);
				//rA = MainFunction(r, i, j, M, N, a, b, h1, h2);
				tauNumerator += rA * r[i][j];
				tauDenominator += rA * rA;
			}
		}
		tau = tauNumerator / tauDenominator;
		deltaSqr = 0.0;
		//посчитать w(k+1)
		//посчитать точность
		for (int i = 1; i < M; ++i)
		{
			for (int j = 1; j < N; ++j)
			{
				double step = tau * r[i][j];
				w[i][j] -= step;

				deltaSqr += step * step;
			}
		}
		#if defined SHOWINFO || defined SHOWDELTAGRAPHIC
		#if defined SHOWINFO
		{
			if (k % TracingPeriod == 0)
			{
				cout << k << endl;
				log << k << ")";
				log << " delta^2 = " << deltaSqr;
				log << " delta^2(k-1) = " << deltaSqr1;
				log << " delta^2(k-2) = " << deltaSqr2;
				log << " tau = " << tau;
				/*log << " tauNumerator = " << tauNumerator;
				log << " tauDenominator = " << tauDenominator;*/
				log << endl;

				#ifdef WRITEFILE
				std::ostringstream oss;
				oss << "f/result" << k << ".txt";
				ofstream fout(oss.str());
				SaveResults(w, N, M, fout);
				fout.close();
				#endif
			}
		}
		#endif // SHOWINFO

		#ifdef SHOWDELTAGRAPHIC
		{
			deltaLog << deltaSqr << endl;
		}
		#endif // SHOWDELTAGRAPHIC
		#endif // OR DEFINED

		if (deltaSqr < DELTA * DELTA)
			break;
		deltaSqr2 = deltaSqr1;
		deltaSqr1 = deltaSqr;
		if (deltaSqr2 <= deltaSqr1 && deltaSqr1 <= deltaSqr || deltaSqr2 == deltaSqr)
			--stopEquals;
		else
			stopEquals = 2 * TracingPeriod;

		if (stopEquals <= 0)
		{
			log << "equals break" << endl;
			break;
		}
	}

	cout << "stop k = " << k << endl;
	cout << "time = " << omp_get_wtime() - start << endl;
	log << "stop k = " << k << endl;
	log << "time = " << omp_get_wtime() - start << endl;
	log.close();
	{
		ofstream fout("f/final.txt");
		SaveResults(w, sizeX, sizeY, fout);
		fout.close();
	}

	#ifdef SHOWDELTAGRAPHIC
	{
		deltaLog.close();
	}
	#endif // SHOWDELTAGRAPHIC

	//Освобождение памяти
	for (int i = 0; i < sizeX; ++i)
	{
		delete[] w[i];
		delete[] r[i];
		delete[] a[i];
		delete[] b[i];
		delete[] F[i];
	}
	delete[] w;
	delete[] r;
	delete[] a;
	delete[] b;
	delete[] F;

	return 0;
}

void ReadParameters(int argc, char **argv)
{
	for (int i = 1; i < argc; i++)
	{
		char c;
		std::istringstream iss(argv[i]);

		iss >> c;
		if (c == 'N')
		{
			iss >> c >> N;
		}
		else if (c == 'M')
		{
			iss >> c >> M;
		}
		else if (c == 'P')
		{
			iss >> c >> TracingPeriod;
		}
	}

}

double CalculateA(double x, double y, double h1, double h2)
{
	if ((x <= P0.X) || (y <= P0.Y) || (x > P1.X) || (y > P1.Y))
		return 0;

	x = x - h1 / 2;
	if (x <= C.X)
		return 1.0;

	double yCB = 9 - 3 * x;
	double aEps = y + h2 / 2 - yCB;
	if (aEps < 0)//под CB
		return 1.0;

	double aOne = yCB - (y - h2 / 2);
	if (aOne < 0)//над CB
		return 1.0 / epsilon;

	return (aEps / epsilon + aOne) / h2;
}

double CalculateB(double x, double y, double h1, double h2)
{
	if ((x <= P0.X) || (y <= P0.Y) || (x > P1.X) || (y > P1.Y))
		return 0;

	if (x + h1 / 2 <= C.X)
		return 1.0;

	double xCB = 3 - (y - h2 / 2) / 3;
	double bEps = x + h1 / 2 - xCB;
	if (bEps < 0)//левее CB
		return 1.0;

	double bOne = xCB - (x - h1 / 2);
	if (bOne < 0)//правее CB
		return 1.0 / epsilon;

	return (bOne + bEps / epsilon) / h1;
}

double CalculateF(double x, double y, double h1, double h2)
{
	if ((x <= P0.X) || (y <= P0.Y) || (x >= P1.X) || (y >= P1.Y))
		return 0;

	double xLeft = x - h1 / 2;
	double xRight = x + h1 / 2;
	if (xRight <= C.X)//правее CB
		return 1.0;

	double yTop = y + h2 / 2;
	if (yTop > C.Y)//не должно произойти, так как вызываем от правильных точек
		yTop = C.Y;

	double yDown = y - h2 / 2;
	if (yDown < B.Y)//не должно произойти, так как вызываем от правильных точек
		yDown = B.Y;

	double yCBTop = yTop;
	double yCBDown = yDown;

	double xCDLeft = 3 - (yTop) / 3;
	if (xCDLeft >= xRight)//левее CB
		return 1.0;

	double xCDRight = 3 - (yDown) / 3;
	if (xCDRight <= xLeft)//правее CB
		return 0;

	if (xCDLeft < xLeft)//Пересечение с левой гранью
	{
		xCDLeft = xLeft;
		yCBTop = 9 - 3 * xCDLeft;
	}
	if (xCDRight > xRight)//Пересечение с правой гранью
	{
		xCDRight = xRight;
		yCBDown = 9 - 3 * xCDRight;
	}

	double dx1 = xCDLeft - xLeft;
	if (dx1 < 0)//По идее, никогда не должно такого произойти
		dx1 = 0;
	double dx2 = xCDRight - xCDLeft;

	double dy1 = yCBTop - yCBDown;
	double dy2 = yCBDown - yDown;

	double s1 = dy1 * (dx1 + dx2 / 2) + dy2 * (dx1 + dx2);

	return (s1) / (h1 * h2);
}

double MainFunction(double **w, int i, int j, int M, int N, double **a, double **b, double h1, double h2)
{
	double center, left, right, top, down;
	if (i <= 0 || i >= M || j <= 0 || j >= N)
		return 0;

	center = w[i][j];
	left = i - 1 == 0 ? 0 : w[i - 1][j];
	right = i + 1 == M ? 0 : w[i + 1][j];
	top = j + 1 == N ? 0 : w[i][j + 1];
	down = j - 1 == 0 ? 0 : w[i][j - 1];

	double dx = (a[i + 1][j] * (right - center) - a[i][j] * (center - left)) / (h1 * h1);
	double dy = (b[i][j + 1] * (top - center) - b[i][j] * (center - down)) / (h2 * h2);
	return -(dx + dy);
}

void SaveResults(double **w, int N, int M, ofstream &fileoutput)
{
	for (int j = 0; j < N; ++j)
	{
		for (int i = 0; i < M; ++i)
		{
			fileoutput << w[i][j] << " ";
		}
		fileoutput << ";" << endl;
	}
}
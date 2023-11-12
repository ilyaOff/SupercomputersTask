#include <iostream>
#include "Point.h"
#include "Task2.h"


using namespace std;

//Параметры области
const Point A = Point(0.0, 0.0);
const Point B = Point(3.0, 0.0);
const Point C = Point(2.0, 3.0);
const Point D = Point(0.0, 3.0);
const Point P0 = Point(0.0, 0.0);
const Point P1 = Point(3.0, 3.0);

const double DELTA = 0.000001;
const double epsilon = 0.001;
const int KMAX = 100000;

int N = 100;
int M = 100;

int main(int argc, char **argv)
{
	ReadParameters(argc, argv);
	//return 0;
	double h1 = (P1.X - P0.X) / (M);
	double h2 = (P1.Y - P0.Y) / (N);
	++N;
	++M;


	double **w = new double *[M];
	double **wNew = new double *[M];
	double **r = new double *[M];
	double **a = new double *[M];
	double **b = new double *[M];
	double **F = new double *[M];
	for (int i = 0; i < M; ++i)
	{
		w[i] = new double[N];
		wNew[i] = new double[N];
		r[i] = new double[N];
		a[i] = new double[N];
		b[i] = new double[N];
		F[i] = new double[N];
	}

	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; j++)
		{
			w[i][j] = 0;
			wNew[i][j] = 0;
			r[i][j] = 0;
			a[i][j] = CalculateA(i, j, h1, h2);
			b[i][j] = CalculateB(i, j, h1, h2);
			F[i][j] = CalculateF(i, j, h1, h2);
		}
	}

	double tau = 0.0;
	double rA = 0.0;
	double tauNumerator = 0.0, tauDenominator = 0.0;
	double deltaResult = 0.0;
	int k = 0;
	for (; k < KMAX; ++k)
	{
		//посчитать невязку r
		for (int i = 1; i < M; ++i)
		{
			for (int j = 1; j < N; ++j)
			{
				double tmp = MainFunction(w, i, j, M, N, a, b, h1, h2);
				tmp -= F[i][j];
				r[i][j] = tmp;
			}
		}
		//посчитать невязку
		for (int i = 1; i < M; ++i)
		{
			for (int j = 1; j < N; ++j)
			{
				rA = MainFunction(r, i, j, M, N, a, b, h1, h2);
				tauNumerator += rA * r[i][j];
				tauDenominator += rA * rA;
			}
		}
		tau = tauNumerator / tauDenominator;
		deltaResult = 0.0;
		//посчитать w(k+1)
		//посчитать точность
		for (int i = 1; i < M; ++i)
		{
			for (int j = 1; j < N; ++j)
			{
				double tmp = tau *r[i][j];
				wNew[i][j] = w[i][j] - tmp;

				if (tmp < 0)
					tmp = -tmp;

				if(deltaResult < tmp)
					deltaResult = tmp;
			}
		}

		if (deltaResult < 0)
			deltaResult = -deltaResult;
	
		if (k % 100 == 0)
		{
			cout << "delta = " << deltaResult;
			cout << " tau = " << tau;
			cout << " tauNumerator = " << tauNumerator;
			cout << " tauDenominator = " << tauDenominator << endl;
		}

		if (deltaResult < DELTA * DELTA)
			break;

		double **swap = w;
		w = wNew;
		wNew = swap;
	}
	cout << "sto[ k = " << k << endl;
	for (int i = 0; i < M; ++i)
	{
		delete[] w[i];
		delete[] wNew[i];
		delete[] r[i];
		delete[] a[i];
		delete[] b[i];
		delete[] F[i];
	}
	delete[] w;
	delete[] wNew;
	delete[] r;
	delete[] a;
	delete[] b;
	delete[] F;

	return 0;
}

void ReadParameters(int argc, char **argv)
{
	if (argc == 1)
	{ // если в аргументах только имя программы
		cout << "no arguments!" << endl; // выводим, что нет аргументов
	}
	else
	{
		// иначе выводим все аргументы, которые переданы
		for (int i = 1; i < argc; i++)
		{
			cout << "argv[" << i << "] - " << argv[i] << endl;
		}
	}
}

double CalculateA(double x, double y, double h1, double h2)
{
	if (x < P0.X)
		return 0;
	if (y < P0.Y)
		return 0;
	if (x > P1.X)
		return 0;
	if (y > P1.Y)
		return 0;

	x = x - h1 / 2;

	if (x <= C.X)
	{
		return 1 / h2;
	}

	double yCD = 9 - 3 * (x - h1 / 2);
	double aEps = y + h2 / 2 - yCD;
	aEps = aEps > 0 ? aEps / epsilon : 0;
	double aOne = yCD - (y - h2 / 2);
	aOne = aOne > 0 ? aOne : 0;
	return (aEps + aOne) / h2;
}

double CalculateB(double x, double y, double h1, double h2)
{
	if (x < P0.X)
		return 0;
	if (y < P0.Y)
		return 0;
	if (x > P1.X)
		return 0;
	if (y > P1.Y)
		return 0;

	if (x + h1 / 2 <= C.X)
	{
		return 1 / h1;
	}

	double xCD = 3 - (y - h2 / 2) / 3;
	double bEps = x + h1 / 2 - xCD;
	bEps = bEps > 0 ? bEps / epsilon : 0;
	double bOne = xCD - (x - h1 / 2);
	bOne = bOne > 0 ? bOne : 0;
	return (bOne + bEps) / h1;
}

double CalculateF(double x, double y, double h1, double h2)
{
	if (x <= P0.X)
		return 0;
	if (y <= P0.Y)
		return 0;
	if (x >= P1.X)
		return 0;
	if (y >= P1.Y)
		return 0;

	if (x + h1 / 2 <= C.X)
		return 1 / (h1 * h2);

	double yTop = y + h2 / 2;
	if (yTop > C.Y)//не должно произойти, так как вызываем от правильных точек
		yTop = C.Y;

	double yDown = y - h2 / 2;
	if (yDown < B.Y)//не должно произойти, так как вызываем от правильных точек
		yDown = B.Y;

	double xCDLeft = 3 - (yTop) / 3;//Всегда >= C.X
	double xCDRight = 3 - (yDown) / 3;//Всегда <= B.X и >= C.X
	double xLeft = x - h1 / 2;
	if (xLeft > C.X)
		xLeft = C.X;

	double dx1 = xCDLeft - xLeft;
	if (dx1 < 0)//По идее, никогда не должно такого произойти
		dx1 = 0;

	double height = yTop - yDown;
	double dx2 = xCDRight - xCDLeft;
	return height * (dx1 + dx2 / 2) / (h1 * h2);
}

double MainFunction(double **w, int i, int j, int M, int N, double **a, double **b, double h1, double h2)
{
	double center, left, right, top, down;
	if (i <= 0 || i >= M - 1 || j <= 0 || j >= N)
		return 0;

	center = w[i][j];
	left = i - 1 == 0 ? 0 : w[i - 1][j];
	right = i + 1 == M ? 0 : w[i + 1][j];
	top = j + 1 == N ? 0 : w[i][j + 1];
	down = j - 1 == N ? 0 : w[i][j - 1];

	double dx = (a[i + 1][j] * (left - center) - a[i][j] * (center - right)) / (h1 * h1);
	double dy = (b[i][j + 1] * (top - center) - a[i][j] * (center - down)) / (h2 * h2);
	return -(dx + dy);
}
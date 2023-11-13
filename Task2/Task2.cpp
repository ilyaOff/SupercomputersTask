﻿#include <iostream>
#include <string> 
#include <sstream>
#include "Point.h"
#include "Task2.h"
//#define WRITEFILE
#define SHOWINFO

using namespace std;

//Параметры области
const Point A = Point(0.0, 0.0);
const Point B = Point(3.0, 0.0);
const Point C = Point(2.0, 3.0);
const Point D = Point(0.0, 3.0);
const Point P0 = Point(0.0, 0.0);
const Point P1 = Point(3.0, 3.0);

const int PERIOD = 10000;
const double DELTA = 0.000001;
const double epsilon = 0.001;
const int KMAX = 100000000;

int N = 40;
int M = 40;

#include<fstream>
using namespace std;

int main(int argc, char **argv)
{
	//ReadParameters(argc, argv);
	//return 0;

	double h1 = (P1.X - P0.X) / (M);
	double h2 = (P1.Y - P0.Y) / (N);

	double sizeX = M + 1;
	double sizeY = N + 1;
	double **w = new double *[sizeX];
	double **wNew = new double *[sizeX];
	double **r = new double *[sizeX];
	double **a = new double *[sizeX];
	double **b = new double *[sizeX];
	double **F = new double *[sizeX];

	for (int i = 0; i < sizeX; ++i)
	{
		w[i] = new double[sizeY];
		wNew[i] = new double[sizeY];
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
			wNew[i][j] = 0;
			r[i][j] = 0;
			double x = P0.X + i * h1;
			double y = P0.Y + j * h2;
			a[i][j] = CalculateA(x, y, h1, h2);
			b[i][j] = CalculateB(x, y, h1, h2);
			F[i][j] = 0.0;
		}
	}

	double tau = 0.0;
	double rA = 0.0;
	double tauNumerator = 0.0, tauDenominator = 0.0;
	double deltaResult2 = 0.0, deltaResult1 = 0.0, deltaResult = 0.0;
	int k = 1, stopEquals = 2 * PERIOD;
	int M1 = M;
	int N1 = N;
	ofstream fout("f/F.txt");
	for (int i = 0; i < M1 + 1; ++i)
	{
		for (int j = 0; j < N1 + 1; ++j)
		{
			double x = P0.X + i * h1;
			double y = P0.Y + j * h2;
			F[i][j] = CalculateF(x, y, h1, h2);
			fout << F[i][j] << " ";
		}
		fout << endl;
	}
	fout.close();
	for (; k < KMAX; ++k)
	{
		//посчитать невязку r
		for (int i = 1; i < M1; ++i)
		{
			for (int j = 1; j < N1; ++j)
			{
				double tmp = MainFunction(w, i, j, M, N, a, b, h1, h2);
				tmp -= F[i][j];
				r[i][j] = tmp;
			}
		}
		//посчитать невязку
		for (int i = 1; i < M1; ++i)
		{
			for (int j = 1; j < N1; ++j)
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
		for (int i = 1; i < M1; ++i)
		{
			for (int j = 1; j < N1; ++j)
			{
				double tmp = tau * r[i][j];
				wNew[i][j] = w[i][j] - tmp;

				if (tmp < 0)
					tmp = -tmp;

				if (deltaResult < tmp)
					deltaResult = tmp;
			}
		}

		/*if (deltaResult < 0)
			deltaResult = -deltaResult;*/

		if (k % PERIOD == 0)
		{
		#ifdef SHOWINFO

			cout << k << ")";
			cout << " delta = " << deltaResult;
			cout << " delta1 = " << deltaResult1;
			cout << " delta2 = " << deltaResult2;
			cout << " tau = " << tau;
			cout << " tauNumerator = " << tauNumerator;
			cout << " tauDenominator = " << tauDenominator << endl;
		#endif // SHOWINFO

		#ifdef WRITEFILE
			std::ostringstream oss;
			oss << "f/result" << k << ".txt";
			string fileName = oss.str();
			SaveResults(w, N, M, fileName.c_str());
		#endif
		}

		if (deltaResult < DELTA)
			break;
		deltaResult2 = deltaResult1;
		deltaResult1 = deltaResult;
		if (deltaResult2 <= deltaResult1 && deltaResult1 <= deltaResult || deltaResult2 == deltaResult)
			--stopEquals;
		else
			stopEquals = 2 * PERIOD;

		if (stopEquals <= 0)
			break;

		double **swap = w;
		w = wNew;
		wNew = swap;
	}
	cout << "stop k = " << k << endl;
	std::ostringstream oss;
	oss << "f/final.txt";
	string fileName = oss.str();
	SaveResults(w, sizeX, sizeY, fileName.c_str());

	for (int i = 0; i < sizeX; ++i)
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
	if (x <= P0.X)
		return 0;
	if (y <= P0.Y)
		return 0;
	/*if (x >= P1.X)
		return 0;
	if (y >= P1.Y)
		return 0;*/

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
	if (x <= P0.X)
		return 0;
	if (y <= P0.Y)
		return 0;
	/*if (x >= P1.X)
		return 0;
	if (y >= P1.Y)
		return 0;*/

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

	double xLeft = x - h1 / 2;
	double xRight = x + h1 / 2;
	if (xRight <= C.X)
		return 1 / (h1 * h2);

	double yTop = y + h2 / 2;
	if (yTop > C.Y)//не должно произойти, так как вызываем от правильных точек
		yTop = C.Y;
	double yCDTop = yTop;

	double yDown = y - h2 / 2;
	if (yDown < B.Y)//не должно произойти, так как вызываем от правильных точек
		yDown = B.Y;
	double yCDDown = yDown;

	double xCDLeft = 3 - (yTop) / 3;
	if (xCDLeft >= xRight)
		return 1 / (h1 * h2);	

	double xCDRight = 3 - (yDown) / 3;//Всегда <= B.X и >= C.X
	if (xCDRight <= xLeft)
		return 0;

	if (xCDLeft < xLeft)
	{
		xCDLeft = xLeft;
		yCDTop = 9 - 3 * xCDLeft;
	}
	if (xCDRight > xRight)
	{
		xCDRight = xRight;
		yCDDown = 9 - 3 * xCDRight;
	}

	double dx1 = xCDLeft - xLeft;
	if (dx1 < 0)//По идее, никогда не должно такого произойти
		dx1 = 0;
	double dx2 = xCDRight - xCDLeft;

	double dy1 = yCDTop - yCDDown;
	double dy2 = yCDDown - yDown;

	double s1 = dy1 * (dx1 + dx2 / 2) + dy2 * (dx1 + dx2);
	double s2 = h1 * h2 - s1;

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

void SaveResults(double **w, int N, int M, const char *fileName)
{
	ofstream fout(fileName);
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			fout << w[i][j] << " ";
		}
		fout << ";" << endl;
	}

	fout.close();
}
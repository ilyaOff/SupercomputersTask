﻿#include <iostream>
#include <string> 
#include <sstream>
#include <fstream>
#include <omp.h>
#include <mpi.h>
#include "Point.h"
#include "Task2.h"
#include "MyMacroses.h"
//#define WRITEFILE
#define WRITEFILER
//#define SHOWINFO
//#define SHOWDELTAGRAPHIC
#define SHOWCOUNT
#define SHOWERRORGRAPHIC

using namespace std;

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
	int numtasks, rank;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	cout << "HELLO MPI. id process = " << rank << " from " << numtasks << " processes" << endl;	

	MPI_Finalize();
	return 0;
	if (argc == 1)
	{
		cout << "no arguments!" << endl;
		return -1;
	}
	else
	{
		ReadParameters(argc, argv);
	}

	cout << "M = " << M << " N = " << N << " Period = " << TracingPeriod << endl;

	if (M <= 0 || N <= 0)
	{
		cout << "invalid parametres (M, N)";
		return -1;
	}

	#ifdef  SHOWDELTAGRAPHIC
	ofstream deltaLog("f/DeltaLog.txt");
	#endif //  SHOWDELTAGRAPHIC

	int sizeX = M + 1;
	int sizeY = N + 1;

	//Выделение памяти под массивы
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

	#ifdef SHOWERRORGRAPHIC
	double **err = new double *[sizeX];
	double **rAGrid = new double *[sizeX];

	for (int i = 0; i < sizeX; ++i)
	{
		err[i] = new double[sizeY];
		rAGrid[i] = new double[sizeY];
	}
	#endif // SHOWERRORGRAPHIC

	//Установка параметров OpenMP
	//omp_set_nested(1);
	//omp_set_dynamic(1);

	cout << "start" << endl;
	double start = omp_get_wtime();

	double h1 = (P1.X - P0.X) / (M);
	double h2 = (P1.Y - P0.Y) / (N);
	epsilon = h1;
	if (h2 > h1)
		epsilon = h2;
	epsilon = epsilon * epsilon;
	
	#pragma omp parallel for
	for (int i = 0; i < sizeX; ++i)
	{
		#pragma omp parallel for
		for (int j = 0; j < sizeY; ++j)
		{
			w[i][j] = 0;
			r[i][j] = 0;
			double x = P0.X + i * h1;
			double y = P0.Y + j * h2;
			a[i][j] = CalculateA(x, y, h1, h2);
			b[i][j] = CalculateB(x, y, h1, h2);
			F[i][j] = CalculateF(x, y, h1, h2);
			#ifdef SHOWERRORGRAPHIC
			err[i][j] = 0;
			#endif
		}
	}
	w[sizeX / 2][sizeY / 2] = 1;

	double norma2R = 0.0;
	double tau = 0.0;
	double rA = 0.0;
	double tauNumerator = 0.0, tauDenominator = 0.0;
	double deltaSqr = 0.0;
	int k = 1;
	int i, j;

	#ifdef WRITEFILER
	ofstream foutR("f/normaR.txt");
	#endif
	//Вывод коэффициентов рассчёта
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
	#endif

	//Основной цикл
	#pragma omp parallel private(i, j, rA, tau)
	for (; k < KMAX; )
	{
		//посчитать невязку r
		#pragma omp  for  collapse(2) schedule(static) 
		for (i = 1; i < M; ++i)
		{
			for (j = 1; j < N; ++j)
			{
				MainFunctionParallel2(r[i][j], -F[i][j], w, i, j, M, N, a, b, h1, h2);
			}
		}

		#pragma omp single nowait
		{
			tauNumerator = 0.0, tauDenominator = 0.0;
			deltaSqr = 0.0;
			++k;
			#ifdef WRITEFILER
			norma2R = 0.0;
			#endif // WRITEFILER
		}
		#pragma omp barrier

		//посчитать итерационный параметр
		#pragma omp for collapse(2) schedule(static) reduction(+:tauNumerator, tauDenominator)
		for (i = 1; i < M; ++i)
		{
			for (j = 1; j < N; ++j)
			{
				MainFunctionParallel2(rA, 0, r, i, j, M, N, a, b, h1, h2);

				#ifdef SHOWERRORGRAPHIC
				rAGrid[i][j] = rA;
				#endif // SHOWERRORGRAPHIC

				tauNumerator += rA * r[i][j];
				tauDenominator += rA * rA;
			}
		}

		//#pragma omp single
		{
			/*if (tauDenominator == 0 || isnan(tauDenominator))
				tau = 0;
			else*/
			tau = tauNumerator / tauDenominator;
		}

		//посчитать w(k+1)
		//посчитать точность
		#pragma omp for schedule(static) collapse(2) nowait reduction(+:deltaSqr)
		for (i = 1; i < M; ++i)
		{
			for (j = 1; j < N; ++j)
			{
				double step = tau * r[i][j];
				w[i][j] = w[i][j] - step;

				deltaSqr += step * step;
			}
		}

		#ifdef SHOWERRORGRAPHIC
		if (k % TracingPeriod == 0)
		{
			#pragma omp for schedule(static) collapse(2)
			for (i = 1; i < M; ++i)
			{
				for (j = 1; j < N; ++j)
				{
					err[i][j] = -tau * r[i][j];
				}
			}

			#pragma omp single nowait
			{
				std::ostringstream oss;
				oss << "f/RAA/errorGrid" << k << ".txt";
				ofstream fout(oss.str());
				SaveResults(rAGrid, sizeX, sizeY, fout);
				fout.close();
			}

			#pragma omp single nowait
			{
				std::ostringstream oss;
				oss << "f/err/errorGrid" << k << ".txt";
				ofstream fout(oss.str());
				SaveResults(err, sizeX, sizeY, fout);
				fout.close();
			}

			#pragma omp single nowait
			{
				std::ostringstream oss;
				oss << "f/Rerr/errorGrid" << k << ".txt";
				ofstream fout(oss.str());
				SaveResults(r, sizeX, sizeY, fout);
				fout.close();
			}
		}
		#endif // SHOWERRORGRAPHIC

		#if defined SHOWCOUNT || defined SHOWDELTAGRAPHIC
		#pragma omp barrier
		#if defined SHOWCOUNT
		#pragma omp single nowait
		{
			if (k % TracingPeriod == 0)
			{
				cout << k << ")";
				cout << " delta^2 = " << deltaSqr;
				cout << " tau = " << tau;
				/*cout << " tauNumerator = " << tauNumerator;
				cout << " tauDenominator = " << tauDenominator;*/
				cout << endl;


				#ifdef WRITEFILE
				std::ostringstream oss;
				oss << "f/result" << k << ".txt";
				ofstream fout(oss.str());
				SaveResults(w, sizeX, sizeY, fout);
				fout.close();

				#endif

				cout << "timeStep = " << (omp_get_wtime() - start) << endl;
			}
		}
		#endif // SHOWINFO

		#ifdef WRITEFILER
		#pragma omp barrier
		if (k % TracingPeriod == 0)
		{
			#pragma omp for schedule(static) collapse(2) reduction(+:norma2R)
			for (i = 1; i < M; ++i)
			{
				for (j = 1; j < N; ++j)
				{
					norma2R += r[i][j] * r[i][j];
				}
			}

			#pragma omp single
			{
				foutR << norma2R << ";" << endl;
				norma2R = 0;
			}
		}
		#endif// WRITEFILER

		#ifdef SHOWDELTAGRAPHIC
		#pragma omp single nowait
		{
			deltaLog << deltaSqr << endl;
		}
		#endif // SHOWDELTAGRAPHIC
		#endif // OR DEFINED

		#pragma omp barrier

		if (deltaSqr < DELTA * DELTA)
			break;
	}

	//#ifdef SHOWCOUNT
	cout << "stop k = " << k << endl;
	cout << "time = " << (omp_get_wtime() - start) << endl;
	//#endif // SHOWCOUT

	#ifdef SHOWDELTAGRAPHIC
	#pragma omp single nowait
	{
		deltaLog << deltaSqr << endl;
		deltaLog.close();
	}
	#endif // SHOWDELTAGRAPHIC

	#ifdef WRITEFILER
	{
		norma2R = 0.0;
		for (i = 1; i < M; ++i)
		{
			for (j = 1; j < N; ++j)
			{
				norma2R += r[i][j] * r[i][j];
			}
		}
		foutR << norma2R << endl;
		foutR.close();
	}
	#endif

	#ifdef SHOWERRORGRAPHIC
	{
		std::ostringstream oss;
		oss << "f/RAA/errorGrid" << k << ".txt";
		ofstream fout(oss.str());
		SaveResults(rAGrid, sizeX, sizeY, fout);
		fout.close();
	}
	{
		std::ostringstream oss;
		oss << "f/Rerr/errorGrid" << k << ".txt";
		ofstream fout(oss.str());
		SaveResults(r, sizeX, sizeY, fout);
		fout.close();
	}
	{
		std::ostringstream oss;
		oss << "f/err/errorGrid" << k << ".txt";
		ofstream fout(oss.str());
		SaveResults(err, sizeX, sizeY, fout);
		fout.close();
	}
	#endif // SHOWERRORGRAPHIC

	//Вывод результата в файл
	{
		#ifdef SHOWINFO
		ofstream fout("f/final.txt");
		SaveResults(w, sizeX, sizeY, fout);
		fout.close();
		#else
		cout << endl << "result:" << endl;
		//SaveResults(w, sizeX, sizeY);
		cout << endl;
		#endif // SHOWINFO

	}

	//Освобождение памяти
	for (int i = 0; i < sizeX; ++i)
	{
		delete[] w[i];
		delete[] r[i];
		delete[] a[i];
		delete[] b[i];
		delete[] F[i];
		#ifdef SHOWERRORGRAPHIC
		delete[] err[i];
		delete[] rAGrid[i];
		#endif // SHOWERRORGRAPHIC

	}
	delete[] w;
	delete[] r;
	delete[] a;
	delete[] b;
	delete[] F;

	#ifdef SHOWERRORGRAPHIC
	delete[] err;
	delete[] rAGrid;
	#endif // SHOWERRORGRAPHIC

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

void SaveResults(double **w, int N, int M)
{
	for (int j = 0; j < N; ++j)
	{
		for (int i = 0; i < M; ++i)
		{
			cout << w[i][j] << " ";
		}
		cout << ";" << endl;
	}
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

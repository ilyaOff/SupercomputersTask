#include <iostream>
#include <string> 
#include <sstream>
#include <fstream>
#include <omp.h>
#include <mpi.h>
#include "Point.h"
#include "Task2.h"
#include "MyMacroses.h"

//#define WRITEFILE
//#define WRITEFILER
//#define SHOWINFO
#define RESULTINFILE
//#define SHOWDELTAGRAPHIC
#define SHOWCOUNT
//#define SHOWERRORGRAPHIC

#define SENDW 10

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
	MPI_Init(&argc, &argv);

	int numtasks, rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	if (rank == 0)
	{
		if (argc == 1)
		{
			cout << "no arguments!" << endl;
			M = N = 0;
		}
		else
		{
			ReadParameters(argc, argv);
		}

		cout << "M = " << M << " N = " << N << " Period = " << TracingPeriod << endl;
	}

	int size[3];
	size[0] = M;
	size[1] = N;
	size[2] = TracingPeriod;
	MPI_Bcast(size, 3, MPI_INT, 0, MPI_COMM_WORLD);
	M = size[0];
	N = size[1];
	TracingPeriod = size[2];

	if (M <= 0 || N <= 0)
	{
		cout << "invalid parametres (M, N)";
		MPI_Finalize();
		return -1;
	}

	double h1 = (P1.X - P0.X) / (M);
	double h2 = (P1.Y - P0.Y) / (N);
	epsilon = h1;
	if (h2 > h1)
		epsilon = h2;
	epsilon = epsilon * epsilon;

	MPI_Comm vu;
	int dims[2];
	CreateGridCommunicator(numtasks, vu, dims);

	int topNode, downNode, leftNode, rightNode;
	MPI_Cart_shift(vu, 1, -1, &topNode, &downNode);
	MPI_Cart_shift(vu, 0, 1, &leftNode, &rightNode);

	int coord[2];
	MPI_Cart_coords(vu, rank, 2, coord);

	int sizeX = CalculateSize(M, dims[0], coord[0]);
	int sizeY = CalculateSize(N, dims[1], coord[1]);

	cout << rank << " UP " << topNode << " down " << downNode << " left " << leftNode << " right " << rightNode << endl;
	cout << rank << " coord = " << coord[0] << ":" << coord[1] << endl;
	cout << rank << " elements = (" << sizeX << "; " << sizeY << ")" << endl;

	#ifdef  SHOWDELTAGRAPHIC
	ofstream deltaLog("f/DeltaLog.txt");
	#endif //  SHOWDELTAGRAPHIC

	//Выделение памяти под массивы
	double **w = new double *[sizeX];
	double **r = new double *[sizeX];
	double **a = new double *[sizeX];
	double **b = new double *[sizeX];
	double **F = new double *[sizeX];
	double *sharedWbyX = new double[sizeX];
	double *sharedWbyX2 = new double[sizeX];

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
	cout << "MaxThereads = " <<omp_get_max_threads() << endl;
	double start = MPI_Wtime();

	for (int i = 0; i < sizeX; ++i)
	{
		sharedWbyX[i] = 0;
	}

	int shiftX = coord[0] * GetCountElementInRow(M, dims[0]);
	if (coord[0] != 0)
		shiftX -= 1;
	int shiftY = coord[1] * GetCountElementInRow(N, dims[1]);
	if (coord[1] != 0)
		shiftY -= 1;
	cout << rank << " shift: " << shiftX << " " << shiftY << endl;
	for (int i = 0; i < sizeX; ++i)
	{
		for (int j = 0; j < sizeY; ++j)
		{
			w[i][j] = 0;
			r[i][j] = 0;
			double x = P0.X + (i + shiftX) * h1;
			double y = P0.Y + (j + shiftY) * h2;
			a[i][j] = CalculateA(x, y, h1, h2);
			b[i][j] = CalculateB(x, y, h1, h2);
			F[i][j] = CalculateF(x, y, h1, h2);
			#ifdef SHOWERRORGRAPHIC
			err[i][j] = 0;
			#endif
		}
	}

	/*int middleX = (M + 1) / 2 - shiftX;
	int middleY = (N + 1) / 2 - shiftY;
	if (middleX >= 0 && middleX < sizeX && middleY >= 0 && middleY < sizeY)
		w[middleX][middleY] = 1;*/

	double norma2R = 0.0;
	double tau = 0.0;
	double rA = 0.0;
	double tauNumerator = 0.0, tauDenominator = 0.0;
	double tau4Send[2], tau4Recive[2];
	double deltaSqr = 0.0, deltaSqr4Recive = 0.0;
	int k = 1;
	int i, j;

	#ifdef WRITEFILER
	ofstream foutR("f/normaR.txt");
	#endif
	//Вывод коэффициентов рассчёта
	#ifdef SHOWINFO

	Write2File("f/F", rank, F, sizeX, sizeY);
	Write2File("f/A", rank, a, sizeX, sizeY);
	Write2File("f/B", rank, b, sizeX, sizeY);

	#endif

	int Mfor = sizeX - 1;
	int Nfor = sizeY - 1;
	int startI = 1;
	int startJ = 1;

	cout << Mfor << Nfor << endl;
	//Основной цикл
	#pragma omp parallel private(i, j, rA, tau)
	for (; k < KMAX; )
	{
		//записать значения массива от соседей
		#pragma omp single
		{
			BorderPointExchange(w, sharedWbyX, sharedWbyX2, sizeX, sizeY, coord, rightNode, leftNode, downNode, topNode, vu);
		}

		//посчитать невязку r
		#pragma omp for collapse(2) schedule(static)
		for (i = startI; i < Mfor; ++i)
		{
			for (j = startJ; j < Nfor; ++j)
			{
				MainFunctionParallel2(r[i][j], -F[i][j], w, i, j, Mfor, Nfor, a, b, h1, h2);
			}
		}

		#pragma omp single
		{
			tauNumerator = 0.0, tauDenominator = 0.0;
			deltaSqr = 0.0;
			++k;
			#ifdef WRITEFILER
			norma2R = 0.0;
			#endif // WRITEFILER
		}


		#pragma omp single
		{
			BorderPointExchange(r, sharedWbyX, sharedWbyX2, sizeX, sizeY, coord, rightNode, leftNode, downNode, topNode, vu);
		}

		//посчитать итерационный параметр
		#pragma omp for collapse(2) schedule(static) reduction(+:tauNumerator, tauDenominator)
		for (i = startI; i < Mfor; ++i)
		{
			for (j = startJ; j < Nfor; ++j)
			{
				MainFunctionParallel2(rA, 0, r, i, j, Mfor, Nfor, a, b, h1, h2);

				#ifdef SHOWERRORGRAPHIC
				rAGrid[i][j] = rA;
				#endif // SHOWERRORGRAPHIC

				tauNumerator += rA * r[i][j];
				tauDenominator += rA * rA;
			}
		}

		#pragma omp single
		{
			//if (rank == 0)
			//	cout << "before " << tauNumerator / tauDenominator << "split " << tauNumerator << " " << tauDenominator << endl;
			tau4Send[0] = tauNumerator;
			tau4Send[1] = tauDenominator;
			MPI_Allreduce(tau4Send, tau4Recive, 2, MPI_DOUBLE, MPI_SUM, vu);
			tauNumerator = tau4Recive[0];
			tauDenominator = tau4Recive[1];

			
			//if (rank == 0)
			//	cout << "after " << tau << endl;
		}
		tau = tauNumerator /( tauDenominator);

		//посчитать w(k+1)
		//посчитать точность
		#pragma omp for schedule(static) collapse(2) reduction(+:deltaSqr)
		for (i = startI; i < Mfor; ++i)
		{
			for (j = startJ; j < Nfor; ++j)
			{
				double step = tau * r[i][j];
				w[i][j] = w[i][j] - step;

				deltaSqr += step * step;
			}
		}

		#pragma omp single
		{		
			MPI_Allreduce(&deltaSqr, &deltaSqr4Recive, 1, MPI_DOUBLE, MPI_SUM, vu);
			deltaSqr = deltaSqr4Recive;
		}

		#ifdef SHOWERRORGRAPHIC
		if (k % TracingPeriod == 0)
		{
			#pragma omp for schedule(static) collapse(2)
			for (i = 1; i < Mfor; ++i)
			{
				for (j = 1; j < Nfor; ++j)
				{
					err[i][j] = -tau * r[i][j];
				}
			}

			#pragma omp single nowait
			{
				Write2FileWithStep("f/RAA/errorGrid", rank, k, rAGrid, sizeX, sizeY);
			}
			#pragma omp single nowait
			{
				Write2FileWithStep("f/err/errorGrid", rank, k, err, sizeX, sizeY);
			}
			#pragma omp single nowait
			{
				Write2FileWithStep("f/Rerr/errorGrid", rank, k, r, sizeX, sizeY);
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
				cout << k << " node " << rank << ")";
				cout << " delta^2 = " << deltaSqr;
				cout << " tau = " << tau;
				cout << endl;

				#ifdef WRITEFILE
				Write2FileWithStep("f/result", rank, k, w, sizeX, sizeY);
				#endif

				cout << "timeStep = " << (MPI_Wtime() - start) << endl;
			}
		}
		#endif // SHOWINFO

		#ifdef WRITEFILER
		#pragma omp barrier
		if (k % TracingPeriod == 0)
		{
			#pragma omp for schedule(static) collapse(2) reduction(+:norma2R)
			for (i = 1; i < Mfor; ++i)
			{
				for (j = 1; j < Nfor; ++j)
				{
					norma2R += r[i][j] * r[i][j];
				}
			}

			MPI_Reduce(&norma2R, &deltaSqr4Recive, 1, MPI_DOUBLE, MPI_SUM, 0, vu);
			norma2R = deltaSqr4Recive;
			#pragma omp single nowait
			{
				if (rank == 0)
				{
					cout << "start write norma2R" << endl;
					foutR << norma2R << ";" << endl;
					cout << "end write norma2R" << endl;
				}
				norma2R = 0;
			}
		}
		#endif// WRITEFILER

		#ifdef SHOWDELTAGRAPHIC
		#pragma omp single nowait
		{
			if (rank == 0)
			{
				deltaLog << deltaSqr << endl;
			}
		}
		#endif // SHOWDELTAGRAPHIC
		#endif // OR DEFINED

		#pragma omp barrier

		if (deltaSqr < DELTA * DELTA)
		{
			break;
		}
	}

	//#ifdef SHOWCOUNT
	cout << "stop k = " << k << endl;
	cout << "time = " << (MPI_Wtime() - start) << endl;
	//#endif // SHOWCOUT
	#ifdef SHOWDELTAGRAPHIC
	//#pragma omp single nowait
	{
		deltaLog << deltaSqr << endl;
		deltaLog.close();
	}
	#endif // SHOWDELTAGRAPHIC

	#ifdef WRITEFILER
	{
		norma2R = 0.0;
		for (i = 1; i < Mfor; ++i)
		{
			for (j = 1; j < Nfor; ++j)
			{
				norma2R += r[i][j] * r[i][j];
			}
		}
		MPI_Reduce(&norma2R, &deltaSqr4Recive, 1, MPI_DOUBLE, MPI_SUM, 0, vu);
		norma2R = deltaSqr4Recive;
		if (rank == 0)
			foutR << norma2R << endl;
		foutR.close();
	}
	#endif

	#ifdef SHOWERRORGRAPHIC

	Write2FileWithStep("f/RAA/errorGrid", rank, k, rAGrid, sizeX, sizeY);
	Write2FileWithStep("f/Rerr/errorGrid", rank, k, r, sizeX, sizeY);
	Write2FileWithStep("f/err/errorGrid", rank, k, err, sizeX, sizeY);

	#endif // SHOWERRORGRAPHIC

	//Вывод результата в файл
	{
		cout << rank << " start write file" << endl;
		#ifdef RESULTINFILE

		Write2File("f/final", rank, w, sizeX, sizeY);

		#else
		cout << endl << "result:" << endl;
		//SaveResults(w, sizeX, sizeY);
		cout << endl;
		#endif // RESULTINFILE
		cout << rank << " end write file" << endl;
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

	delete[] sharedWbyX;
	delete[] sharedWbyX2;

	#ifdef SHOWERRORGRAPHIC
	delete[] err;
	delete[] rAGrid;
	#endif // SHOWERRORGRAPHIC

	MPI_Finalize();
	return 0;
}

void BorderPointExchange(double **w, double *sharedWbyX, double *sharedWbyX2, int sizeX, int sizeY, int(&coord)[2], int rightNode, int leftNode, int downNode, int topNode, const MPI_Comm &vu)
{
	MPI_Status status;

	//if (coord[1] % 2 == 0)
	{
		if (downNode >= 0)
		{
			for (int i = 0; i < sizeX; ++i)
			{
				sharedWbyX[i] = w[i][1];
			}
			MPI_Sendrecv(sharedWbyX, sizeX, MPI_DOUBLE, downNode, SENDW,
						 sharedWbyX2, sizeX, MPI_DOUBLE, downNode, SENDW, vu, &status);
			for (int i = 0; i < sizeX; ++i)
			{
				w[i][0] = sharedWbyX2[i];
			}
		}

		if (topNode >= 0)
		{
			for (int i = 0; i < sizeX; ++i)
			{
				sharedWbyX[i] = w[i][sizeY - 2];
			}
			MPI_Sendrecv(sharedWbyX, sizeX, MPI_DOUBLE, topNode, SENDW,
						 sharedWbyX2, sizeX, MPI_DOUBLE, topNode, SENDW, vu, &status);
			for (int i = 0; i < sizeX; ++i)
			{
				w[i][sizeY - 1] = sharedWbyX2[i];
			}
		}
	}

	//if (coord[0] % 2 == 0)
	{
		if (rightNode >= 0)
		{
			MPI_Sendrecv(w[sizeX - 2], sizeY, MPI_DOUBLE, rightNode, SENDW,
						 w[sizeX - 1], sizeY, MPI_DOUBLE, rightNode, SENDW, vu, &status);
		}

		if (leftNode >= 0)
		{
			MPI_Sendrecv(w[1], sizeY, MPI_DOUBLE, leftNode, SENDW,
						 w[0], sizeY, MPI_DOUBLE, leftNode, SENDW, vu, &status);
		}
	}
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

void CreateGridCommunicator(int numtasks, MPI_Comm &vu, int *dims)
{
	int period[2], reorder;


	dims[0] = 0; dims[1] = 0;
	MPI_Dims_create(numtasks, 2, dims);

	period[0] = false; period[1] = false;
	reorder = false;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, &vu);
}

int GetCountElementInRow(int lengthBigGrid, int maxElemets)
{
	return (lengthBigGrid) / maxElemets;
}

int CalculateSize(int lengthBigGrid, int maxElemets, int gridCoordinate)
{
	int size = GetCountElementInRow(lengthBigGrid, maxElemets);
	if (gridCoordinate == maxElemets - 1)
		size = (lengthBigGrid + 1) - size * gridCoordinate;

	size += 2;//Для обмена с соседними узлами на сетке
	if (gridCoordinate == 0)//у первого узла на соседа меньше
		size -= 1;
	//не ifelse, так как элемент может быть один по данному измерению
	if (gridCoordinate == maxElemets - 1)// у последнего узла на соседа меньше
		size -= 1;
	return size;
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

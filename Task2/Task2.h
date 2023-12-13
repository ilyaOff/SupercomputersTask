#pragma once

void ReadParameters(int argc, char **argv);

void CreateGridCommunicator(int numtasks, MPI_Comm &vu, int *dims);
int GetCountElementInRow(int lengthBigGrid, int maxElemets);
int CalculateSize(int lengthBigGrid, int maxElemets, int gridCoordinate);

double CalculateA(double x, double y, double h1, double h2);
double CalculateB(double x, double y, double h1, double h2);
double CalculateF(double x, double y, double h1, double h2);

double MainFunction(double **w, int i, int j, int M, int N, double **a, double **b, double h1, double h2);
void BorderPointExchange(double **w, double *sharedWbyX, int sizeX, int sizeY, int (&coord)[2], int rightNode, int leftNode, int downNode, int topNode, const MPI_Comm &vu);

void SaveResults(double **w, int N, int M, std::ofstream &fileoutput);
void SaveResults(double **w, int N, int M);
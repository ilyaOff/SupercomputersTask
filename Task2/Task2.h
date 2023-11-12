#pragma once

void ReadParameters(int argc, char **argv);

double CalculateA(double x, double y, double h1, double h2);
double CalculateB(double x, double y, double h1, double h2);
double CalculateF(double x, double y, double h1, double h2);

double MainFunction(double **w, int i, int j, int M, int N, double **a, double **b, double h1, double h2);

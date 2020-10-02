#ifndef INC_4ISLAKI_COMMON_H
#define INC_4ISLAKI_COMMON_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

double **CreatingMatrix(int n, int m);

double **CreatingMatrix(int n);

double *GaussMethod(double **matrix, double *vec, int size);

void DeletingMatrix(double **X, int x);

void MatrixMultiply(double **matrix, double *roots, double *res, int size);

#endif //INC_4ISLAKI_COMMON_H

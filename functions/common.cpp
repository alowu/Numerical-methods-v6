#include "common.h"

double **CreatingMatrix(int n, int m) {
    auto **A = new double *[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[m];
    return A;
}

double **CreatingMatrix(int n) {
    auto **A = new double *[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];
    return A;
}

double *GaussMethod(double **matrix, double *vec, int size) {
    double *x, max;
    int k, index;
    const double eps = 1e-6;
    x = new double[size];
    k = 0;
    while (k < size) {
        max = abs(matrix[k][k]);
        index = k;
        for (int i = k + 1; i < size; i++) {
            if (abs(matrix[i][k]) > max) {
                max = abs(matrix[i][k]);
                index = i;
            }
        }
        double temp;
        for (int j = 0; j < size; j++) {
            temp = matrix[k][j];
            matrix[k][j] = matrix[index][j];
            matrix[index][j] = temp;
        }
        temp = vec[k];
        vec[k] = vec[index];
        vec[index] = temp;

        for (int i = k; i < size; i++) {
            try {
                temp = matrix[i][k];
                //if (abs(temp) < eps) continue;
                if (abs(temp) < eps) throw "IER = 1";
            } catch (const char *msg) {
                cout << msg << '\n';
                exit(1);
            }
            for (int j = 0; j < size; j++) {
                matrix[i][j] /= temp;
            }
            vec[i] /= temp;
            if (i == k) continue;
            for (int j = 0; j < size; j++) {
                matrix[i][j] = matrix[i][j] - matrix[k][j];
            }
            vec[i] = vec[i] - vec[k];
        }
        k++;
    }

    for (k = size - 1; k >= 0; k--) {
        x[k] = vec[k];
        for (int i = 0; i < k; i++)
            vec[i] -= matrix[i][k] * x[k];
    }
    return x;
}

void DeletingMatrix(double **X, int x) {
    for (int i = 0; i < x; i++) {
        delete X[i];
    }
    delete[] X;
}

void MatrixMultiply(double **matrix, double *roots, double *res, int size) {
    for (int i = 0; i < size; i++) {
        double temp = 0;
        for (int j = 0; j < size; j++) {
            temp += matrix[i][j] * roots[j];
        }
        res[i] = temp;
    }
}
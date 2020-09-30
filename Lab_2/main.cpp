#include <iostream>
#include <iomanip>

using namespace std;

const double M = 0.5;
const int NIT = 100;

double **CreatingMatrix(int n, int m) {
    auto **A = new double *[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[m];
    return A;
}

void DeletingMatrix(double **X, int x) {
    for (int i = 0; i < x; i++) {
        delete X[i];
    }
    delete[] X;
}

double f1(double x1, double x2) {
    return (2 * x1 * x1 * x1 - x2 * x2 - 1);
    //return x1 * x1 * x1 + x2 * x2 * x2 - 6 * x1 + 3;//v5
}

double f2(double x1, double x2) {
    return (x1 * x2 * x2 * x2 - x2 - 4);
    //return x1 * x1 * x1 - x2 * x2 * x2 - 6 * x2 + 2;//v5
}

double func11(double x1, double x2) {
    return (f1(x1 + M * x1, x2) - f1(x1, x2)) / M * x1;
}

double func12(double x1, double x2) {
    return (f1(x1, x2 + M * x2) - f1(x1, x2)) / M * x2;
}

double func21(double x1, double x2) {
    return (f2(x1 + M * x1, x2) - f2(x1, x2)) / M * x1;
}

double func22(double x1, double x2) {
    return (f2(x1, x2 + M * x2) - f2(x1, x2)) / M * x2;
}

void nev(double *F, double x1, double x2) {
    F[0] = -f1(x1, x2);
    F[1] = -f2(x1, x2);
}

void JIter(double **matrix, double x1, double x2) {
    matrix[0][0] = func11(x1, x2);
    matrix[0][1] = func12(x1, x2);
    matrix[1][0] = func21(x1, x2);
    matrix[1][1] = func22(x1, x2);
}

void JDiff(double **matrix, double x1, double x2) {
    matrix[0][0] = 6 * x1 * x1;
    matrix[0][1] = -2 * x2;
    matrix[1][0] = x2 * x2 * x2;
    matrix[1][1] = 3 * x1 * x2 * x2 - 1;
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

double *NewtonMethod(const int SIZE, double x1, double x2) {
    const double EPS = 1e-9;

    double **Jak = CreatingMatrix(SIZE, SIZE);

    auto *F = new double[SIZE];
    auto *delta = new double[SIZE];
    auto *roots = new double[SIZE];
    auto *b = new double[SIZE];

    roots[0] = x1;
    roots[1] = x2;

    double delta1, delta2;
    int k = 1;
    printf("-----------------------------------------------------------\n");
    printf("|      x1   |      x2   |    delta1   |      delta2    | k \n");
    printf("-----------------------------------------------------------\n");

    do {
        nev(F, x1, x2);
        JIter(Jak, x1, x2);

        for (int i = 0; i < SIZE; ++i) b[i] = F[i];

        delta = GaussMethod(Jak, b, SIZE);

        for (int i = 0; i < SIZE; ++i) {
            roots[i] += delta[i];
        }

        double max1 = 0;
        double max2 = 0;
        for (int i = 0; i < SIZE; ++i) {
            if (abs(F[i]) > max1) max1 = abs(F[i]);
            if (abs(roots[i]) < 1) {
                if (abs(delta[i]) > max2) max2 = abs(delta[i]);
            } else {
                if (abs(delta[i] / roots[i]) > max2) max2 = abs(delta[i] / roots[i]);
            }
        }

        delta1 = max1;
        delta2 = max2;

        x1 = roots[0];
        x2 = roots[1];
        for (int i = 0; i < SIZE; i++) {
            cout << '|' << setw(10) << roots[i] << ' ';
        }
        cout << '|' << setw(13) << delta1 << '|' << setw(13) << delta2 << "   " << '|' << setw(2) << k;
        cout << '\n';

        ++k;
        try {
            if (k >= NIT) throw "IER = 2";
        } catch (const char *msg) {
            cout << msg << '\n';
            exit(2);
        }

    } while (delta1 > EPS && delta2 > EPS);

    DeletingMatrix(Jak, SIZE);
    delete[] F;
    delete[] delta;
    delete[] b;

    return roots;
}

int main() {
    const int SIZE = 2;
    double x1 = 1, x2 = 1;

    auto *roots = new double[SIZE];
    roots = NewtonMethod(SIZE, x1, x2);
    printf("------------------------ Roots ----------------------------\n");
    for (int i = 0; i < SIZE; ++i) {
        cout << "X" << i + 1 << " = " << roots[i] << '\n';
    }
    delete[] roots;
}

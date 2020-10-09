#include "../functions/common.h"

const double M = 0.01;
const int NIT = 100;

double f1(double x1, double x2) {
    return (2 * x1 * x1 * x1 - x2 * x2 - 1);// v6
}

double f2(double x1, double x2) {
    return (x1 * x2 * x2 * x2 - x2 - 4);//v6
}

double func11(double x1, double x2) {
    return (f1(x1 + M, x2) - f1(x1, x2)) / M;
}

double func12(double x1, double x2) {
    return (f1(x1, x2 + M) - f1(x1, x2)) / M;
}

double func21(double x1, double x2) {
    return (f2(x1 + M, x2) - f2(x1, x2)) / M;
}

double func22(double x1, double x2) {
    return (f2(x1, x2 + M) - f2(x1, x2)) / M;
}

void nev(double *F, double x1, double x2) {
    F[0] = -1. * f1(x1, x2);
    F[1] = -1. * f2(x1, x2);
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

    while (true) {
        nev(F, x1, x2);
        JDiff(Jak, x1, x2);

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
        if (delta1 <= EPS && delta2 <= EPS) break;
    }

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
#include "../functions/common.h"

/*
 * --!tests!--
21.547 -95.520 -96.121 -49.930
10.223 -91.065 -7.343 -12.465
51.218 12.264 86.457 60.812

2.30 5.70 -0.80 -6.49
3.50 -2.70 5.30 19.20
1.70 2.30 -1.80 -5.09

 21.547 -95.510 -96.121 -49.930
 10.223 -91.065 -7.343 -12.465
 51.218 12.264 86.457 60.812

 2.21 3.65 1.69 6.99 -8.35
 8.30 2.62 4.10 1.90 -10.65
 3.92 8.45 7.78 2.46 12.21
 3.77 7.21 8.04 2.28 15.45

 0.1582 1.1675 0.1768 0.1871 1.6471
 0.1968 0.2071 1.2168 0.2271 1.7471
 0.2368 0.2471 0.2568 1.2671 1.8471
 1.1161 0.1254 0.1397 0.1490 1.5471

 2 3 11 5 2
 1 1 5 2 1
 2 1 3 2 -3
 1 1 3 4 -3
 */

void InputingOfTheSystem(double **matrix, double **copyM, double *vec, double *copyV, int size) {
    for (int i = 0; i < size; ++i) {
        cout << "input " << i + 1 << " row" << '\n';
        for (int j = 0; j < size; ++j) {
            cin >> matrix[i][j];
            copyM[i][j] = matrix[i][j];
        }
        cin >> vec[i];
        copyV[i] = vec[i];
    }
}

void OutputingOfTheSystem(double **matrix, double *vec, int size) {
    cout << "----------------------------------------------\n";
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cout << setw(10) << matrix[i][j];
        }
        cout << setw(2) << '|' << setw(10) << vec[i] << '\n';
    }
    cout << "----------------------------------------------\n";
}

void OutputingRoots(double *roots, int size) {
    for (int i = 0; i < size; ++i) {
        cout << i + 1 << " root is " << roots[i] << '\n';
    }
}

int main() {
    int size;
    cout << "Enter the number of X: ";
    cin >> size;

    double **matrix = CreatingMatrix(size);
    double **copyOfTheOrigMatrix = CreatingMatrix(size);

    auto *factor = new double[size];
    auto *copyOfTheOrigFactor = new double[size];
    auto *roots = new double[size];

    InputingOfTheSystem(matrix, copyOfTheOrigMatrix, factor, copyOfTheOrigFactor, size);
    OutputingOfTheSystem(matrix, factor, size);

    roots = GaussMethod(matrix, factor, size);
    OutputingRoots(roots, size);

    auto *inconspicuous = new double[size];

    MatrixMultiply(copyOfTheOrigMatrix, roots, inconspicuous, size);

    cout << "----------------------------------------------\n";
    for (int i = 0; i < size; ++i) {
        inconspicuous[i] -= copyOfTheOrigFactor[i];
        cout << "Norma vectora " << i + 1 << " : " << inconspicuous[i] << '\n';
    }
    cout << "----------------------------------------------\n";

    auto *newB = new double [size];
    MatrixMultiply(copyOfTheOrigMatrix, roots, newB, size);
    auto *newRoots = new double[size];

    newRoots = GaussMethod(copyOfTheOrigMatrix, newB, size);

    cout << "----------------------------------------------\n";
    for (int i = 0; i < size; ++i) {
        double pogr = abs(newRoots[i] - roots[i]) / abs(roots[i]);
        cout << "Otnositelnaya pogreshnost' " << i + 1 << " : " << pogr << '\n';
    }
    cout << "----------------------------------------------\n";

    DeletingMatrix(matrix, size);
    DeletingMatrix(copyOfTheOrigMatrix, size);
    delete[] factor;
    delete[] copyOfTheOrigFactor;
    delete[] roots;
    delete[] inconspicuous;
    delete[] newB;
    delete[] newRoots;
}
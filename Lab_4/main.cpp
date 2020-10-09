#include <iostream>
//include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

//степень аппроксимирующего полинома; число измерений
const int m = 1, N = 7;

double **createMatrix(int x) {
    auto **A = new double *[x];
    for (int i = 0; i < x; i++)
        A[i] = new double[x + 1];
    return A;
}

void deleteMatrix(double **X, int x) {
    for (int i = 0; i < x; i++) {
        delete X[i];
    }
    delete[] X;
}

void enterMatrix(double *X, int x, ifstream &in) {
    for (int i = 0; i < x; i++) {
        double element;
        in >> element;
        X[i] = element;
    }
}

void copyMatrix(double **X, double **copyX, int x) {
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < x + 1; j++) {
            copyX[i][j] = X[i][j];
        }
    }
}

void viewAnswer(double *X, int x) {
    for (int i = 0; i < x; i++) {
        cout << "a[" << i << "] = " << X[i] << endl;
    }
}

bool Gauss(double *An, double **X, int x) {

    for (int k = 0; k < x; k++) {
        double max = fabs(X[k][k]);
        int remeber = k;
        for (int i = k + 1; i < x; i++) {
            if (max < fabs(X[i][k])) {
                max = fabs(X[i][k]);
                remeber = i;
            }
        }

        if (fabs(max - 0) < 1e-6) {
            exit(0);
        }

        if (k != remeber) {
            double *temp = X[k];
            X[k] = X[remeber];
            X[remeber] = temp;
        }

        double lead = X[k][k];
        for (int r = k; r < x + 1; r++) {
            X[k][r] /= lead;
        }

        for (int i = k + 1; i < x; i++) {
            double temp = X[i][k];
            for (int j = k; j < x + 1; j++) {
                X[i][j] -= X[k][j] * temp;
            }
        }

    }

    An[x - 1] = X[x - 1][x + 1 - 1];
    for (int i = x - 2; i >= 0; i--) {
        An[i] = X[i][x + 1 - 1];
        for (int j = i + 1; j < x + 1 - 1; j++) {
            An[i] -= X[i][j] * An[j];
        }
    }
    return true;
}

double valueSigma(double *Datax, const double *Datay, const double *An, double x) {
    double sigma = 0, temp;
    for (int i = 0; i < N; i++) {
        temp = Datay[i];
        for (int j = 0; j < x; j++) {
            temp -= An[j] * pow(Datax[i], j);
        }
        sigma += temp * temp;
    }
    return sigma /= (N - m - 1);
}

int main() {
    auto *Datax = new double[N];
    auto *Datay = new double[N];
    ifstream inx(R"(H:\CLionLabs\Numerical-methods-v6\Lab_4\data_x.txt)");
    ifstream iny(R"(H:\CLionLabs\Numerical-methods-v6\Lab_4\data_y.txt)");
    enterMatrix(Datax, N, inx);
    enterMatrix(Datay, N, iny);

    ofstream foun("H:/CLionLabs/Numerical-methods-v6/Lab_4/n.dat", ios_base::out | ios_base::trunc | ios_base::binary);
    int n = N;
    foun.write((char *)&n, sizeof n);
    foun.close();

    double **A = createMatrix(m + 1);

    for (int i = 0; i < m + 1; i++) {
        for (int j = 0; j < m + 1; j++) {
            A[i][j] = 0;
            for (int k = 0; k < N; k++) {
                A[i][j] += pow(Datax[k], i + j);
            }
        }
    }
    A[0][0] = N;

    for (int i = 0; i < m + 1; i++) {
        for (int j = 0; j < N; j++) {
            A[i][m + 1] += Datay[j] * pow(Datax[j], i);
        }
    }

    auto *An = new double[m + 1];
    double **copyA = createMatrix(m + 1);
    copyMatrix(A, copyA, m + 1);

    Gauss(An, copyA, m + 1);
    viewAnswer(An, m + 1);

    cout << "Sigma = " << sqrt(valueSigma(Datax, Datay, An, m + 1)) << endl;

    double w;
    ofstream foutx("H:/CLionLabs/Numerical-methods-v6/Lab_4/x.dat", ios_base::out | ios_base::trunc | ios_base::binary);
    for (int i = 0; i < N; ++i) {
        w = Datax[i];
        foutx.write((char*)&w, sizeof w);
    }
    for (int i = 0; i < N; ++i) {
        w = Datay[i];
        foutx.write((char*)&w, sizeof w);
    }
    auto* newDatay = new double[N];
    for (int i = 0; i < N; ++i) {
        newDatay[i] = An[1] * Datax[i] + An[0];
    }
    for (int i = 0; i < N; ++i) {
        w = newDatay[i];
        foutx.write((char*)&w, sizeof w);
    }
    foutx.close();

    delete[] An;
    delete[] Datax;
    delete[] Datay;
    delete[] newDatay;
    deleteMatrix(A, m + 1);
    deleteMatrix(copyA, m + 1);
}
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x) {
    return exp(x / 2) / sqrt(x + 1);
}

double f(double x, double y) {
    return x * x / (1 + y * y);
}

double TrapeziumMethod(double a, double b, int &iter) {
    double sum = 0, prevSum = 0;
    double h = b - a;
    double eps = 1e-4;

    do {
        prevSum = sum;
        sum = 0;
        for (int i = 1; i <= (b - a) / h - 1; ++i) {
            sum += 2 * f(a + h * i);
        }
        sum += f(a) + f(b);
        sum *= h / 2;
        h /= 2;
        ++iter;
    } while (abs(sum - prevSum) > 3 * eps);

    return sum;
}

double SimpsonMethod(double a, double b, int &iter) {
    double sum = 0, prevSum = 0;
    double h = (b - a) / 2;
    double eps = 1e-5;

    do {
        prevSum = sum;
        sum = 0;
        for (int i = 1; i <= (b - a) / h; i += 2) {
            sum += 4 * f(a + h * i);
        }
        for (int i = 2; i < (b - a) / h - 1; i += 2) {
            sum += 2 * f(a + h * i);
        }
        sum += f(a) + f(b);
        sum *= h / 3;
        h /= 2;
        ++iter;
    } while (abs(sum - prevSum) > 15 * eps);

    return sum;
}

double CubatureSimpson(double a, double b, double c, double d) {
    int m = 10;
    int n = m * 2;

    double hx = (b - a) / (2 * n);
    double hy = (d - c) / (2 * m);
    double sum = 0;

    double xi = a;
    double yi = c;

    auto *Xi = new double[2 * n + 1];
    Xi[0] = xi;

    for (int i = 1; i <= 2 * n; ++i) {
        Xi[i] = Xi[i - 1] + hx;
    }

    auto *Yi = new double[2 * m + 1];
    Yi[0] = yi;

    for (int i = 1; i <= 2 * m; ++i) {
        Yi[i] = Yi[i - 1] + hy;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            sum += f(Xi[2 * i], Yi[2 * j]);
            sum += 4 * f(Xi[2 * i + 1], Yi[2 * j]);
            sum += f(Xi[2 * i + 2], Yi[2 * j]);
            sum += 4 * f(Xi[2 * i], Yi[2 * j + 1]);
            sum += 16 * f(Xi[2 * i + 1], Yi[2 * j + 1]);
            sum += 4 * f(Xi[2 * i + 2], Yi[2 * j + 1]);
            sum += f(Xi[2 * i], Yi[2 * j + 2]);
            sum += 4 * f(Xi[2 * i + 1], Yi[2 * j + 2]);
            sum += f(Xi[2 * i + 2], Yi[2 * j + 2]);
        }
    }
    sum *= hx * hy / 9;
    return sum;

}

int main() {
    double a = 0; //cin >> a;
    double b = 1.047; //cin >> b;

    int countOfIteration = 0;

    double resultTrapeziumMethod = TrapeziumMethod(a, b, countOfIteration);
    cout << "================== Trapezium Method ==================\n";
    cout << "Result: " << setw(10) << resultTrapeziumMethod <<
         "\nNumber of method iterations: " << setw(3) << countOfIteration << '\n';

    countOfIteration = 0;
    double resultSimpsonMethod = SimpsonMethod(a, b, countOfIteration);
    cout << "=================== Simpson Method ===================\n";
    cout << "Result: " << setw(10) << resultSimpsonMethod <<
         "\nNumber of method iterations: " << setw(3) << countOfIteration << '\n';

    a = 0;
    b = 4;
    double c = 1;
    double d = 2;

    countOfIteration = 0;
    double resultCubatureSimpsonMethod = CubatureSimpson(a, b, c, d);
    cout << "============= Cubature Simpson Method ================\n";
    cout << "Result: " << setw(10) << resultCubatureSimpsonMethod;

}
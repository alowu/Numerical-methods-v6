#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
const double de = 1e-9;
const int itr_max = 100;

//выделение памяти с возращением указателя
double **createMatrix(int x) {
    auto **A = new double *[x];
    for (int i = 0; i < x; i++)
        A[i] = new double[x + 1];
    return A;
}

//освобождение памяти
void deleteMatrix(double **X, int x) {
    for (int i = 0; i < x; i++)
        delete X[i];
    delete[] X;
}

//копирования для сохранения оригинала
void copyMatrix(double **X, double **copyX, int x) {
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < x + 1; j++)
            copyX[i][j] = X[i][j];
    }
}

//отображение в консоли
void viewMatrix(double **X, int x) {
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < x + 1; j++)
            cout << setw(12) << X[i][j];
        cout << endl;
    }
    cout << endl << endl;
}

//отображение в консоли ответов
void viewAnswer(double *X, int x) {
    for (int i = 0; i < x; i++)
        cout << setw(12) << X[i];
    cout << endl;
}

//весь метод Гаусса
bool Gauss(double *An, double **X, int x) {
    for (int k = 0; k < x; k++) {
        double max = fabs(X[k][k]);
        int remeber = k;        //запоминаем строку, чтобы не поменять саму себя
        for (int i = k + 1; i < x; i++) {
            if (max < fabs(X[i][k]))        //находим максимальный по модулю элемент в столбце
            {
                max = fabs(X[i][k]);
                remeber = i;
            }
        }

        if (fabs(max - 0) < 1e-6) {
            return 0;
        }

        if (k != remeber)                //меняем строки местами
        {
            double *temp = X[k];
            X[k] = X[remeber];
            X[remeber] = temp;
        }

        //viewMatrix(X, x);

        double lead = X[k][k];            //запоминаем ведущий элемент
        for (int r = k; r < x + 1; r++) {
            X[k][r] /= lead;
        }
        //начиная со следующей строки приводим исходную матрицу к диагональному виду
        for (int i = k + 1; i < x; i++) {
            double temp = X[i][k];
            for (int j = k; j < x + 1; j++) {
                X[i][j] -= X[k][j] * temp;
            }
        }
        //viewMatrix(X, x);
    }

    An[x - 1] = X[x - 1][x + 1 - 1];                //обратный ход
    for (int i = x - 2; i >= 0; i--) {
        An[i] = X[i][x + 1 - 1];
        for (int j = i + 1; j < x + 1 - 1; j++) {
            An[i] -= X[i][j] * An[j];
        }
    }
    return 1;
}

//первое уравнение системы
double f1(double x1, double x2) {
    //return cos(0.4*x2 + x1*x1) + x2*x2 + x1*x1 - 1.6;
    //return 1.5*x1*x1*x1 - x2*x2 - 1;
    return sin(x1 - x2) - 1.32;
    //return x1*x1*x2*x2 - 3 * x1*x1 - 6 * x2*x2*x2 + 8; //проверка Вар3
    //return x1 - x2 - 6 * log10(x1) - 1;
}

//второе уравнение системы
double f2(double x1, double x2) {
    //return 1.5*x1*x1 - x2*x2 / 0.36 - 1;
    //return x1*x2*x2*x2 - x2 - 3;
    //return x1*x1*x1*x1 - 9 * x2 + 2;   //проверка Вар3
    //return x1 - 3 * x2 - 6 * log10(x2) - 2;
    return cos(x2 - x1) + 0.85;
}

typedef double(*pf)(double, double); //pf указатель на функцию возвращающую double
//подсчет производной через дифференциал (как lim)
//Аргументы - (уравнение системы, переменные, по какой переменной ищем)
double Differential(pf f, double x1, double x2, int x) {
    if (x == 1)                                            //по первой неизвестной
        return (f(x1 + de, x2) - f(x1, x2)) / de;
    else                                                //по второй неизвестной
        return (f(x1, x2 + de) - f(x1, x2)) / de;
}

//метод Ньютона
double Newton(double *function, double *approximate, int n) {
    double **F = createMatrix(n);
    auto *An = new double[n];
    double e = 1e-9, b1 = 0, b2 = 0;

    cout << endl << setw(15) << "k" << setw(15) << "b1      " << setw(15) << "b2      " << endl << endl;
    int itr = 1;

    do {
        //заполнение матрицы Якоби в ячейки F[0][0], F[0][1], F[1][0], F[1][1]
        F[0][0] = Differential(f1, approximate[0], approximate[1], 1);
        F[0][1] = Differential(f1, approximate[0], approximate[1], 2);
        F[0][2] = -f1(approximate[0], approximate[1]);                    //вектор столбец с -
        F[1][0] = Differential(f2, approximate[0], approximate[1], 1);
        F[1][1] = Differential(f2, approximate[0], approximate[1], 2);
        F[1][2] = -f2(approximate[0], approximate[1]);                    //вектор столбец с -

        double **copyF = createMatrix(n);                    // копия для сохранения оригинала
        copyMatrix(F, copyF, n);

        if (!Gauss(An, copyF, n))                //Проверка метода Гаусса
        {
            deleteMatrix(copyF, n);        //очистка памяти
            deleteMatrix(F, n);
            delete[] An;
            //	cout << itr;
            return 0;
        }

        approximate[0] += An[0];                //уточнение решения
        approximate[1] += An[1];

        if (fabs(f1(approximate[0], approximate[1])) >
            fabs(f2(approximate[0], approximate[1]))) //подсчет первой погрешности
            b1 = fabs(f1(approximate[0], approximate[1]));
        else
            b1 = fabs(f2(approximate[0], approximate[1]));

        for (int i = 0; i < n; i++)                                        //подсчет второй погрешности
        {
            if (fabs(An[i]) < 1)
                b2 = fabs(An[i]);
            else if (fabs(An[i]) >= 1)
                b2 = fabs(An[i] / approximate[i]);
        }

        cout << setw(15) << itr << setw(15) << b1 << setw(15) << b2 << endl;
        itr++;
    } while ((b1 > e || b2 > e) && (itr < itr_max));

    function[n - 2] = approximate[0];
    function[n - 1] = approximate[1];

    //deleteMatrix(copyF, n);
    deleteMatrix(F, n);                //очистка памяти
    delete[] An;

    return *function;    //возвращаем решение СНЛАУ
}

int main() {
    const int n = 2;
    auto *Approximate = new double[n];
    auto *Function = new double[n];
    Approximate[0] = 1;
    Approximate[1] = 0;
    *Function = Newton(Function, Approximate, n);
    if (*Function != 0) {
        cout << "\n\n Answer: ";
        viewAnswer(Function, n);
    } else
        cout << "\n\n There is no solution!\n\n";

    /*Approximate[0] = -1; Approximate[1] = 1;
    *Function = Newton(Function, Approximate, n);
    if (*Function != 0)
    {
        cout << "\n\n Answer: ";
        viewAnswer(Function, n);
    }
    else
        cout << "\n\n There is no solution!\n\n";*/
}
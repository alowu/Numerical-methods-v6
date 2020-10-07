#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;
//������� ����������������� ��������; ����� ���������

//const int m = 3;  const int N = 21;  //var 7
//const int m = 2; const int N = 10; //var 8
//const int m = 1; const int N = 5; //var 10
const int m = 1 , N = 7;

//��������� ������ � ����������� ���������
double** createMatrix(int x)
{
    double **A = new double*[x];
    for (int i = 0; i < x; i++)
        A[i] = new double[x + 1];
    return A;
}
//������������ ������
void deleteMatrix(double **X, int x)
{
    for (int i = 0; i < x; i++)
        delete X[i];
    delete[] X;
}
//���������� ������� ��������
void enterMatrix(double *X, int x, ifstream &in)
{
    for (int i = 0; i < x; i++)
    {
        double element; in >> element;
        X[i] = element;
    }
}
//����������� � �������
void viewMatrix(double **X, int x)
{
    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < x + 1; j++)
            cout << setw(15) << X[i][j];
        cout << endl;
    }
    cout << endl << endl;
}
//����������� ��� ���������� ���������
void copyMatrix(double **X, double **copyX, int x)
{
    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < x + 1; j++)
            copyX[i][j] = X[i][j];
    }
}
//�������� �������
void zeroMatrix(double **X, int x)
{
    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < x + 1; j++)
        {
            X[i][j] = 0.0;
        }
    }
}
//����������� � ������� �������
void viewAnswer(double *X, int x)
{
    for (int i = 0; i < x; i++)
        cout << "a[" << i << "] = " << X[i] << endl;
}
//����� ������ � �������� ���
bool Gauss(double *An, double **X, int x)
{

    for (int k = 0; k < x; k++)
    {
        double max = fabs(X[k][k]);
        int remeber = k;		//���������� ������, ����� �� �������� ���� ����
        for (int i = k + 1; i < x; i++)
        {
            if (max < fabs(X[i][k]))		//������� ������������ �� ������ ������� � �������
            {
                max = fabs(X[i][k]);
                remeber = i;
            }
        }

        if (fabs(max - 0) < 1e-6)
        {
            return 0;
        }

        if (k != remeber)				//������ ������ �������
        {
            double *temp = X[k];
            X[k] = X[remeber];
            X[remeber] = temp;
        }

        //viewMatrix(X, row, column);

        double lead = X[k][k];			//���������� ������� �������
        for (int r = k; r < x + 1; r++)
        {
            X[k][r] /= lead;
        }
        //������� �� ��������� ������ ��������� �������������� ������
        for (int i = k + 1; i < x; i++)
        {
            double temp = X[i][k];
            for (int j = k; j < x + 1; j++)
            {
                X[i][j] -= X[k][j] * temp;
            }
        }
        //viewMatrix(A, n, m);
    }

    An[x - 1] = X[x - 1][x + 1 - 1];				//�������� ���
    for (int i = x - 2; i >= 0; i--)
    {
        An[i] = X[i][x + 1 - 1];
        for (int j = i + 1; j < x + 1 - 1; j++)
        {
            An[i] -= X[i][j] * An[j];
        }
    }
    return 1;
}

//�������������������� ����������
double valueSigma(double *Datax, double *Datay, double *An, double x)
{
    double sigma = 0, temp;
    for (int i = 0; i < N; i++)	{
        temp = Datay[i];
        for (int j = 0; j < x; j++)
            temp -= An[j] * pow(Datax[i], j);  // ��������� 4.5
        sigma += temp*temp;// � � �������
    }
    return sigma /= (N - m - 1);
}

int main()
{
    double *Datax = new double[N];
    double *Datay = new double[N];
    ifstream inx(R"(H:\CLionLabs\Numerical-methods-v6\Lab_4\data_x.txt)");
    ifstream iny(R"(H:\CLionLabs\Numerical-methods-v6\Lab_4\data_y.txt)");
    enterMatrix(Datax, N, inx);
    enterMatrix(Datay, N, iny);

//	for (int j = 0; j < N; j++)
//			cout << Datax[j] << "   " << Datay[j] << endl;

    double **A = createMatrix(m + 1);
    zeroMatrix(A, m + 1);

    for (int i = 0; i < m + 1; i++)   ///  ���������� �������
    {
        for (int j = 0; j < m + 1; j++)
        {
            for (int k = 0; k < N; k++)
                A[i][j] += pow(Datax[k], i + j);;
        }
    }
    A[0][0] = N;

    for (int i = 0; i < m + 1; i++) // ������ �����
    {
        for (int j = 0; j < N; j++)
            A[i][m + 1] += Datay[j] * pow(Datax[j], i);
    }

    double *An = new double[m + 1];
    double **copyA = createMatrix(m + 1);
    copyMatrix(A, copyA, m + 1);

    Gauss(An, copyA, m + 1);
    viewAnswer(An, m + 1);

    cout << "Sigma = " << sqrt(valueSigma(Datax, Datay, An, m + 1)) << endl;
    // �� ������� ������ ����  ���������� � ��������



    delete An;
    delete Datax;
    delete Datay;
    deleteMatrix(A, m + 1);
    deleteMatrix(copyA, m + 1);
}
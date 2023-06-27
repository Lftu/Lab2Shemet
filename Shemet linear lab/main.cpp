#include <bits/stdc++.h>
#include <locale>;
#include "windows.h";

using namespace std;

const int n = 7;
const int N = 7;

void inversion(double **A)
{
    double temp;

    double **E = new double *[N];

    for (int i = 0; i < N; i++)
        E[i] = new double [N];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    for (int k = 0; k < N; k++)
    {
        temp = A[k][k];

        for (int j = 0; j < N; j++)
        {
            A[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < N; i++)
        {
            temp = A[i][k];

            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = A[i][k];

            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = E[i][j];

    for (int i = 0; i < N; i++)
        delete [] E[i];

    delete [] E;
}

void PrintMatrix(double **A, double *F)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << A[i][j] << " ";
        cout << F[i];
        cout << endl;
    }
}
void PrintMatrix(double **A)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
}
void PrintVec(double *A)
{
    for (int i = 0; i < n; i++)
    {
        cout << A[i] << " ";
    }
    cout << endl;
}
void Print(double *A)
{
    for (int i = 0; i < n; i++)
    {
        cout <<"x"<<i+1<<"="<< A[i] << " ";
    }
    cout << endl;
}
void Jacobi(double**A, double* F, double* X)
{
    double c = 0;
    double e;
    cout << "Введіть точність розрахунку: ";
    cin >> e;
    double norm, *TempX = new double[n];
    for (int k = 0; k < n; k++)
        TempX[k] = X[k];
    int cnt = 0;
    do
    {
        for (int i = 0; i < n; i++)
        {
            TempX[i] = F[i];
            for (int g = 0; g < n; g++)
                if (i != g)
                    TempX[i] -= A[i][g] * X[g];
            TempX[i] /= A[i][i];
        }
        norm = abs(X[0] - TempX[0]);
        for (int h = 0; h < n; h++)
        {
            if (abs(X[h] - TempX[h]) > norm)
                norm = abs(X[h] - TempX[h]);
            X[h] = TempX[h];
        }
        cnt++;
        c++;
        cout << "Норма №" << c << " = " << norm << endl;
    }
    while (norm > e);
    double *nev = new double[n];
    cout << "Вектор нев'язки: (";
    for (int i = 0; i < n; i++)
    {
        nev[i] = -F[i];
        for (int k = 0; k < n; k++)
            nev[i] += A[i][k] * TempX[k];
        cout << nev[i];
        if (i != n - 1)
            cout << ", ";
        else
            cout << ")" << endl;
    }
    cout << "Кількість ітерацій методу Якобі = " << cnt << endl;
    delete[] TempX;
}
void input(double **&A, double *&F) /// A - Matrix
{
    F = new double[n];
    A = new double *[n];
    for (int i = 0; i < n; i++)
        A[i] = new double[n];
    for (int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            A[i][j] = 0;
            if(i == j)
            {
                A[i][j] = 4.0 + 5.0 / double(i + 1);
            }
            if(abs(i - j) == 1)
            {
                A[i][j] = 1;
            }
        }
        if(i == 0)
        {
            F[i] = n;
        }
        else
        {
            F[i] = n + 1 - 2 * i;
        }
        A[i][i] *= 10;
    }
    A[0][n - 1] = 5;
    A[n - 1][0] = 5;
    F[n - 1] = 1;
}

bool checkMatrix(double **A, double *F) /// Достатня умова збіжності
{
    for(int i = 0; i < n; i++)
    {
        int sum = 0;
        for(int j = 0; j < n; j++)
        {
            if(i != j)
            {
                sum += abs(A[i][j]);
            }
        }
        if(abs(A[i][i]) < sum)
        {
            return false;
        }
    }
    return true;
}

void print(int n, int t, double **z)// вывести матрицу
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < t; j++)
        {
            cout << setw(15) << z[i][j];
        }
        cout << endl;
    }
}

double** transp(double** a, int n, int m)// транспонировать матрицу
{
    int i, j;
    double **b;
    b = new double *[n];
    for (i = 0; i< n; i++)
    {
        b[i] = new double[m];
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
            b[j][i] = a[i][j];
    }
    return b;
}

int sign(double z)
{
    if (z > 0)
        return 1;
    else
        return -1;
}

double Norm(double **a)
{
    double temp = 0, norm_m = 0;
    for(int i = 0; i < n; i++)
    {
        temp = 0;
        for(int j = 0; j < n; j++)
            temp += fabs(a[i][j]);
        if(temp > norm_m)
            norm_m = temp;
    }
    return norm_m;
}

void SquareRoot(double **a, double *b)
{
    double **at = new double*[n], **rab = new double*[n], **s = new double*[n], **d = new double*[n];
    double opred = 1, *nev = new double[n], *y = new double[n], *x = new double[n], *b1 = new double[n];
    int i, j, k;
    bool flag = true;
    double kst;
    for (i = 0; i < n; i++)
    {
        d[i] = new double[n];
        s[i] = new double[n];
        rab[i] = new double[n];
        at[i] = new double[n];
        for (j = 0; j < n; j++)
        {
            if (a[i][j] != a[j][i])
                flag = false;
            d[i][j] = 0;
            s[i][j] = 0;
        }
    }
    if (!flag)
    {
        cout << "Матриця не симетрична, необхідно домножити на транспоновану до неї" << endl;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                at[j][i] = a[i][j];
        }
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                rab[i][j] = 0;
                b1[i] = 0;
                for (k = 0; k < n; k++)
                {
                    rab[i][j] += at[i][k] * a[k][j];
                    b1[i] += at[i][k] * b[k];
                }
            }
        }
    }
    else
    {
        cout << "Матриця симетрична" << endl;
        for (i = 0; i < n; i++)
        {
            b1[i] = b[i];
            for (j = 0; j < n; j++)
                rab[i][j] = a[i][j];
        }
    }
    cout << "Робоча матриця:" << endl;
    print(n, n, rab);
    for (int i = 0; i< n; i++)
    {
        for (int k = 0; k < (i + 1); k++)
        {
            double sum = 0;
            for (int j = 0; j < k; j++)
                sum += s[i][j] * s[k][j]*d[j][j];
            if (i == k)
            {
                s[i][k] = sqrt(abs(rab[i][i] - sum));
                d[i][k] = sign(rab[i][i] - sum);
            }
            else
                s[i][k]= (1.0 / s[k][k] * (rab[i][k] - sum));
        }
    }
    s = transp(s, n, n);
    cout << "Mатриця S:" << endl;
    print(n, n, s);
    cout << "Вектор Y: (";
    for (i = 0; i < n; i++)
    {
        double sum = 0;
        for (int k = 0; k <= i - 1; k++)
            sum += y[k] * s[k][i];
        y[i] = (b1[i] - sum) / s[i][i];
        cout << y[i];
        if (i != n - 1)
            cout << ", ";
        else
            cout << ")";
    }
    cout << endl;
    cout << "Вектор X: (";
    for (i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int k = i + 1; k <= n - 1; k++)
            sum += s[i][k] * x[k];
        x[i] = (y[i] - sum) / s[i][i];
    }
    for (i = 0; i < n; i++)
    {
        opred *= s[i][i] * s[i][i] * d[i][i];
        cout << x[i];
        if (i != n - 1)
            cout << ", ";
        else
            cout << ")";
    }
    cout << endl;
    cout << "Визначник матриці: "<<opred<<endl;
    cout << "Вектор нев'язки: (";
    for (i = 0; i < n; i++)
    {
        nev[i] = -b[i];
        for (k = 0; k < n; k++)
            nev[i] += a[i][k] * x[k];
        cout << nev[i];
        if (i != n - 1)
            cout << ", ";
        else
            cout << ")" << endl;
    }
    double n1 = Norm(a);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            rab[i][j] = a[i][j];
    }
    inversion(rab);
    double n2 = Norm(rab);
    cout << "Матриця A^-1:" << endl;
    print(n, n, rab);
    cout << "Норма матриці A:" << endl;
    cout << n1 << endl;
    cout << "Норма матриці A^-1:" << endl;
    cout << n2 << endl;
    cout << "Число обумовленості матриці A:" << endl;
    cout << n1 * n2 << endl;
}

int main()
{
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    setlocale (LC_CTYPE, "ukr");
    double **Matrix, *b, *y, *x;
    input(Matrix, b);

    cout << "1) Метод квадратного кореня: " << endl;

    SquareRoot(Matrix, b);

    cout << "2) Метод Якобі: " << endl;
    if(checkMatrix(Matrix, b))
    {
        x = new double[n]; /// Початкові оцінки
        for (int i = 0; i < n; i++)
            x[i] = 1.0;
        y = new double[n];
        Jacobi(Matrix, b, x);
        cout << "Результат: ";
        Print(x);
    }
    else
    {
        cout << "Достатня умова збіжності методу Якобі не виконується" << endl;
    }
    system("pause");
    return 0;
}

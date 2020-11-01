#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Precision 1e-12
#define MAXK 1100
#define N 501

double C[6][N+1] = { 0 };

int s = 2;

//幂法所需函数
double Func1(double u[]);
double Func2(double a[], double b[]);
double PowerMethod();
//Doolittle分解所需函数
double sum1(int k, int j);
double sum2(int k, int i);
double sum3(int i, double t[]);
double sum4(int i, double x[]);
double *Doolittle(double y[]);
//反幂法
double InversePower();

int main()
{
	//将A中的元素存入C中
	int temp = 0;
	for (int i=1,j = 1; i <= N-2; i++)
	{
		j = i;
		temp = i - j + s + 1;
		C[temp][j] = ((1.64 - 0.024 * i) * sin(0.2 * i)) - 0.64 * exp(0.1 / i);
		C[temp-1][j + 1] = 0.16;
		C[temp-2][j + 2] = -0.064;
		C[temp + 1][j] = 0.16;
		C[temp + 2][j] = -0.064;
	}
	int i = N-1, j = N-1;
	temp = i - j + s + 1;
	C[temp][j] = ((1.64 - 0.024 * i) * sin(0.2 * i)) - 0.64 * exp(0.1 / i);
	C[temp-1][j + 1] = 0.16;
	C[temp + 1][j] = 0.16;
	i++;
	j++;
	C[temp][j] = ((1.64 - 0.024 * i) * sin(0.2 * i)) - 0.64 * exp(0.1 / i);


	//求按模最大的特征值-执行幂法
	double lambda1 = PowerMethod();
	

	//求另一个最大（小）的特征值
	//原点平移
	for (int i = 1; i <= N; i++)
	{
		C[temp][i] -= lambda1;
	}

	//再一次执行幂法
	double lambda2 = PowerMethod();
	double lambdan = lambda1 + lambda2;
	printf("最小的特征值lambda 1 = %.12e\n", lambda1);
	printf("最大的特征值lambda 501 = %.12e\n", lambdan);


	//求模最小的特征值
	double lambdas = InversePower();
	printf("模最小的特征值lambda s = %.12e\n", 1/lambdas);

	system("PAUSE");
	return 0;
}


double Func1(double u[])
{
	double result = 0;
	for (int i = 1; i <= N; i++)
	{
		result += u[i] * u[i];
	}
	return sqrt(result);
}

double Func2(double a[], double b[])
{
	double result = 0;
	for (int i = 1; i <= N; i++)
	{
		result += a[i] * b[i];
	}
	return result;
}

double PowerMethod()
{
	double u[N + 1] = { 1 ,1 };
	u[N] = 1;
	double yita, y[N + 1];
	int k = 0;
	double b[2] = { 0 };
	do {
		b[0] = b[1];
		k++;
		yita = Func1(u);
		for (int i = 1; i <= N; i++)
		{
			y[i] = u[i] / yita;
		}

		double temp;
		int border;
		int ctemp;
		for (int i = 1; i <= N; i++)
		{
			temp = 0;//记录u每个元素的结果
			border = i + 2;
			ctemp = 0;
			for (int j = i - 2; j <= border; j++)
			{
				if (j < 1)
				{
					j = 1;
				}
				if (border > N)
				{
					border = N;
				}
				ctemp = i - j + s + 1;
				temp += C[ctemp][j] * y[j];
			}
			u[i] = temp;
		}

		b[1] = Func2(y, u);
		

	} while (fabs(b[1] - b[0]) / fabs(b[1]) > Precision || k == 1);
	printf("%d %.12e\n", k, b[1]);

	return b[1];
}

double sum1(int k, int j)
{
	double result = 0;
	int temp1,temp2;
	int tmax = fmax(1, k - s);
	for (int t = 1; t <= k - 1; t++)
	{
		temp1 = k - t + s + 1;
		temp2 = t - j + s + 1;
		result += C[temp1][t] * C[temp2][j];
	}
	return result;
}

double sum2(int k, int i)
{
	double result = 0;
	int temp1, temp2;
	for (int t = 1; t <= k - 1; t++)
	{
		temp1 = i - t + s + 1;
		temp2 = t - k + s + 1;
		result += C[temp1][t] * C[temp2][k];
	}
	return result;
}

double sum3(int i,double t[])
{
	double result = 0;
	int temp;
	for (int z = 1; z <= i - 1; z++)
	{
		temp = i - z + s + 1;
		result += C[temp][z] * t[z];
	}
	return result;
}

double sum4(int i,double x[])
{
	double result = 0;
	int temp;
	for (int z = i + 1; z <= N; z++)
	{
		temp = i - z + s + 1;
		result += C[temp][z] * x[z];
	}
	return result;
}

double *Doolittle(double y[])
{
	double x[N + 1];
	double t[N + 1];
	int temp;
	//Doolittle分解
	for (int k = 1; k <= N; k++)
	{
		for (int j = k; j <= fmin(k+s,N); j++)
		{
			temp = k - j + s + 1;
			C[temp][j] = C[temp][j] - sum1(k, j);
		}
		if (k < N)
		{
			for (int i = k + 1; i <= fmin(k+s,N); i++)
			{
				temp = i - k + s + 1;
				C[temp][k] = (C[temp][k] - sum2(k, i)) / C[s+1][k];
			}
		}

	}

	//回代求解
	t[1] = y[1];
	for (int i = 2; i <= N; i++)
	{
		t[i] = y[i] - sum3(i,t);
	}

	x[N] = t[N] / C[s+1][N];
	for (int i = N - 1; i >= 1; i--)
	{
		x[i] = (t[i] - sum4(i,x)) / C[s+1][i];
	}
	return x;
}

double InversePower()
{
	double *u = new double[N+1];
	u[1] = u[2] = 1;
	u[N] = 1;
	double yita, y[N + 1];
	int k = 0;
	double b[2] = { 0 };
	do {
		b[0] = b[1];
		k++;
		yita = Func1(u);
		for (int i = 1; i <= N; i++)
		{
			y[i] = u[i] / yita;
		}

		u = Doolittle(y);
		b[1] = Func2(y, u);
	} while (fabs(b[1] - b[0]) / fabs(b[1]) > Precision || k == 1);

	delete []u;
	return b[1];
}

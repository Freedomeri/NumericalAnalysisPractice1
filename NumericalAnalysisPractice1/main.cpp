#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Precision 1e-12
#define MAXK 1100
#define N 501

double C[6][N+1] = { 0 };
double A[6][N + 1] = { 0 };//�洢C�ĸ���

int s = 2;

//����C��Ԫ�ص�A
void Copy();
//�ݷ����躯��
double Func1(double u[]);
double Func2(double a[], double b[]);
double PowerMethod();
//Doolittle�ֽ����躯��
double sum1(int k, int j);
double sum2(int k, int i);
double sum3(int i, double t[]);
double sum4(int i, double x[]);
void Doolittle();
//���ݷ�
double InversePower();

int main()
{
	//��A�е�Ԫ�ش���C��
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

	Copy();
	Doolittle();
	//��ģ��������ֵ-ִ���ݷ�
	double lambda1 = PowerMethod();
	

	//����һ�����С��������ֵ
	//ԭ��ƽ��
	for (int i = 1; i <= N; i++)
	{
		C[temp][i] -= lambda1;
	}

	//��һ��ִ���ݷ�
	double lambda2 = PowerMethod();
	double lambdan = lambda1 + lambda2;
	printf("��С������ֵlambda 1 = %.12e\n", lambda1);
	printf("��������ֵlambda 501 = %.12e\n", lambdan);

	//ƽ�ƻ�ԭ����
	for (int i = 1; i <= N; i++)
	{
		C[temp][i] += lambda1;
	}

	//��ģ��С������ֵ
	double lambdas = InversePower();
	printf("ģ��С������ֵlambda s = %.12e\n", 1/lambdas);

	//����uk��ӽ�������ֵ
	double u[40],lambda[40];
	for (int k = 1; k <= 39; k++)
	{
		u[k] = lambda1 + k / 40.0 * (lambdan - lambda1);
		for (int t = 1; t <= N; t++)
		{
			C[temp][t] -= u[k];
		}
		lambda[k] = InversePower();
		printf("u%d = %.12e����u%d��ӽ�������ֵΪ��%.12e\n",k, u[k],k, 1/lambda[k]+u[k]);
		//��ԭ
		for (int t = 1; t <= N; t++)
		{
			C[temp][t] += u[k];
		}
	}

	//��A���׷���������ʽֵ
	printf("A���׷���cond(A)2Ϊ��%.12e\n", fabs(lambda1) * fabs(lambdas));
	
	//printf("A������ʽֵdetAΪ��%.12e\n", );

	system("PAUSE");
	return 0;
}

void Copy()
{
	for (int i = 1; i <= s+s+1; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			A[i][j] = C[i][j];
		}
	}
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
			temp = 0;//��¼uÿ��Ԫ�صĽ��
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
	//printf("%d %.12e\n", k, b[1]);

	return b[1];
}

double sum1(int k, int j)
{
	double result = 0;
	int temp1,temp2;
	int tmax = 1>(k-s)?1:(k-s);
	if (tmax < j - s)
	{
		tmax = j - s;
	}
	for (int t = tmax; t <= k - 1; t++)
	{
		temp1 = k - t + s + 1;
		temp2 = t - j + s + 1;
		result += A[temp1][t] * A[temp2][j];
	}
	return result;
}

double sum2(int k, int i)
{
	double result = 0;
	int temp1, temp2;
	int tmax = 1 > (i - s) ? 1 : (i - s);
	if (tmax < k - s)
	{
		tmax = k - s;
	}
	for (int t = tmax; t <= k - 1; t++)
	{
		temp1 = i - t + s + 1;
		temp2 = t - k + s + 1;
		result += A[temp1][t] * A[temp2][k];
	}
	return result;
}

double sum3(int i,double t[])
{
	double result = 0;
	int temp;
	int tmax = 1 > (i - s) ? 1 : (i - s);
	for (int z = tmax; z <= i - 1; z++)
	{
		temp = i - z + s + 1;
		result += A[temp][z] * t[z];
	}
	return result;
}

double sum4(int i,double x[])
{
	double result = 0;
	int temp;
	int tmin = (i + s) < N ? (i + s) : N;
	for (int z = i + 1; z <= tmin; z++)
	{
		temp = i - z + s + 1;
		result += A[temp][z] * x[z];
	}
	return result;
}

void Doolittle()
{
	//double x[N + 1];
	//double t[N + 1];
	int temp;
	//Doolittle�ֽ�
	int jmin;
	for (int k = 1; k <= N; k++)
	{
		jmin = ((k + s)) < N ? (k + s) : N;
		for (int j = k; j <= jmin; j++)
		{
			temp = k - j + s + 1;
			A[temp][j] = A[temp][j] - sum1(k, j);
		}
		if (k < N)
		{
			for (int i = k + 1; i <= jmin; i++)
			{
				temp = i - k + s + 1;
				A[temp][k] = (A[temp][k] - sum2(k, i)) / A[s+1][k];
			}
		}

	}
}

double InversePower()
{
	
	double *u = new double[N+1];
	u[1]  = 1;
	double t[N + 1];
	/*for (int i = 1; i <= N; i++)
	{
		u[i] = 1;
	}*/
	double yita, y[N + 1];
	int k = 0;
	double b[2] = { 0 };
	
	do {
		//Copy();
		b[0] = b[1];
		k++;
		yita = Func1(u);
		for (int i = 1; i <= N; i++)
		{
			y[i] = u[i] / yita;
		}

		//�ش����
		t[1] = y[1];
		for (int i = 2; i <= N; i++)
		{
			t[i] = y[i] - sum3(i, t);
		}

		u[N] = t[N] / A[s + 1][N];
		for (int i = N - 1; i >= 1; i--)
		{
			u[i] = (t[i] - sum4(i, u)) / A[s + 1][i];
		}

		b[1] = Func2(y, u);
	} while (fabs(b[1] - b[0]) / fabs(b[1]) > Precision || k == 1);

	//printf("%d\n", k);
	return b[1];
}

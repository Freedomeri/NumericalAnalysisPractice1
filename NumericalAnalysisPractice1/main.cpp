#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double C[6][502];
int s = 2;

int main()
{
	//将A中的元素存入C中
	int temp = 0;
	for (int i=1,j = 1; i <= 499; i++)
	{
		j = i;
		temp = i - j + s + 1;
		C[temp][j] = ((1.64 - 0.024 * i) * sin(0.2 * i)) - 0.64 * exp(0.1 / i);
		C[temp][j + 1] = 0.16;
		C[temp][j + 2] = -0.064;
		C[temp + 1][j] = 0.16;
		C[temp + 2][j] = -0.064;
	}
	int i = 500, j = 500;
	temp = i - j + s + 1;
	C[temp][j] = ((1.64 - 0.024 * i) * sin(0.2 * i)) - 0.64 * exp(0.1 / i);
	C[temp][j + 1] = 0.16;
	C[temp + 1][j] = 0.16;
	i++;
	j++;
	C[temp][j] = ((1.64 - 0.024 * i) * sin(0.2 * i)) - 0.64 * exp(0.1 / i);


}

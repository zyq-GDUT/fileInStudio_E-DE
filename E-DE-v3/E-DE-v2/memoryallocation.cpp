#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
*Ϊ�ַ���ָ������ڴ�
*/
char* malloc1Char(int D) {
	char* str = (char*)malloc(D);
	if (str == NULL) {
		printf("malloc error!");
		exit(1);
	}
	return str;

}
void malloc1D(double*& a, int D) {
	a = (double*)malloc(D * sizeof(double));
	if (a == NULL)
		perror("malloc");
}

void malloc1E(int*& a, int D) {
	a = (int*)malloc(D* sizeof(int));
	if (a == NULL)
		perror("malloc");
}

void malloc2D(int**& a, int xDim, int yDim)
{
	a = (int**)malloc(xDim * sizeof(int*));//a�Ƕ�ά���飬����ά����a����xDim���ܴ�����ָ����ڴ�
	a[0] = (int*)malloc(xDim * yDim * sizeof(int));//�ȸ���άָ��ĵ�һ��Ԫ�ط����ʼ��ַ���õ�ַ��ָ��������ά������ڴ�
	for (int i = 1; i < xDim; i++)
	{
		a[i] = a[i - 1] + yDim;//��ʣ�µ�Ԫ�ظ�ֵ
	}
	if (a == NULL)
		perror("malloc");
}

void malloc2E(double**& a, int xDim, int yDim)
{
	a = (double**)malloc(xDim * sizeof(double*));
	a[0] = (double*)malloc(xDim * yDim * sizeof(double));//a[0]��ʾָ���ά����a���׵�ַ
	for (int i = 1; i < xDim; i++)
	{
		a[i] = a[i - 1] + yDim;//a�Ƕ�ά���飬a[i]��ʾ��ά������ÿ��Ԫ�ص��׵�ַ����ά�����ÿ��Ԫ�ض���һά���顣
	}
	if (a == NULL)
		perror("malloc");
}
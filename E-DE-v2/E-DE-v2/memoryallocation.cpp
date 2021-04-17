#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
*为字符串指针分配内存
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
	a = (int**)malloc(xDim * sizeof(int*));//a是二维数组，给二维数组a分配xDim个能存整型指针的内存
	a[0] = (int*)malloc(xDim * yDim * sizeof(int));//先给二维指针的第一个元素分配初始地址，该地址是指向整个二维数组的内存
	for (int i = 1; i < xDim; i++)
	{
		a[i] = a[i - 1] + yDim;//给剩下的元素赋值
	}
	if (a == NULL)
		perror("malloc");
}

void malloc2E(double**& a, int xDim, int yDim)
{
	a = (double**)malloc(xDim * sizeof(double*));
	a[0] = (double*)malloc(xDim * yDim * sizeof(double));//a[0]表示指向二维数组a的首地址
	for (int i = 1; i < xDim; i++)
	{
		a[i] = a[i - 1] + yDim;//a是二维数组，a[i]表示二维数组中每个元素的首地址，二维数组的每个元素都是一维数组。
	}
	if (a == NULL)
		perror("malloc");
}
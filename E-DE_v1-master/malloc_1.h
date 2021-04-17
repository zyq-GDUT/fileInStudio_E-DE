#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double **loadData(int *d, int *n);
void malloc1D(double *&a, int D);
void malloc1E(int *&a, int D);
void malloc2D(int **&a, int xDim, int yDim);
void malloc2E(double **&a, int xDim, int yDim);

/*****************************************************************
* func_name: load the testing data as you need 
* input: dimension and data number
* descript: load data from the txt file
*d是指数据的维数，n是指数据样本的总数
*****************************************************************/
double **loadData(int *d, int *n)
{
	int i, j;
	double **arraydata; 
	FILE *fp;

	if ((fp = fopen("testing_data.txt", "r")) == NULL)    
			fprintf(stderr, "cannot open data.txt!\n");
			
	//fscanf(fp, "D=%d,N=%d\n", d, n)的返回值是从文件读取的数据个数。这里的D=%d，N=%d只是为了方便可读性,实际上意思是读取文件中两个整型数据
	if (fscanf(fp, "D=%d,N=%d\n", d, n) != 2)   fprintf(stderr, "load error!\n");

	malloc2E(arraydata, *n, *d); //给二维数组arraydata分配一个n行d列的数组，该数组刚好能把整个数据样本保存起来

	for (i = 0; i<*n; i++)
		for (j = 0; j<*d; j++)
			fscanf(fp, "%lf", &arraydata[i][j]); //从保存数据样本的文件中读取数据到二维数组arraydata

	return arraydata;
}
// double **loadData2(int &d, int &n)
// {
// 	int i, j;
// 	double **arraydata; 
// 	FILE *fp;

// 	if ((fp = fopen("testing_data.txt", "r")) == NULL)    
// 			fprintf(stderr, "cannot open data.txt!\n");
			
// 	fscanf(fp, "%d", &d);
// 	fscanf(fp, "%d", &n);
// 	malloc2E(arraydata, n, d); //给二维数组arraydata分配一个n行d列的数组，该数组刚好能把整个数据样本保存起来

// 	for (i = 0; i<n; i++)
// 		for (j = 0; j<d; j++)
// 			fscanf(fp, "%lf", &arraydata[i][j]); //从保存数据样本的文件中读取数据到二维数组arraydata

// 	return arraydata;
// }

/*****************************************************************
* func_name: array function
* input: array name and dimension
* descript: allocate the int or double array
*****************************************************************/
void malloc1D(double *&a, int D) {
	a = (double *)malloc(D * sizeof(double));
	if (a == NULL)
		perror("malloc");
}

void malloc1E(int *&a, int D) {
	a = (int *)malloc(D 
	* sizeof(int));
	if (a == NULL)
		perror("malloc");
}

void malloc2D(int **&a, int xDim, int yDim)
{
	a = (int **)malloc(xDim * sizeof(int *));//a是二维数组，给二维数组a分配xDim个能存整型指针的内存
	a[0] = (int *)malloc(xDim * yDim * sizeof(int));//先给二维指针的第一个元素分配初始地址，该地址是指向整个二维数组的内存
	for (int i = 1; i<xDim; i++)
	{
		a[i] = a[i - 1] + yDim;//给剩下的元素赋值
	}
	if (a == NULL)
		perror("malloc");
}

void malloc2E(double **&a, int xDim, int yDim)
{
	a = (double **)malloc(xDim * sizeof(double *));
	a[0] = (double *)malloc(xDim * yDim * sizeof(double));//a[0]表示指向二维数组a的首地址
	for (int i = 1; i<xDim; i++)
	{
		a[i] = a[i - 1] + yDim;//a是二维数组，a[i]表示二维数组中每个元素的首地址，二维数组的每个元素都是一维数组。
	}
	if (a == NULL)
		perror("malloc");
}
// int main(){
// double **t;
// malloc2E(t,3,4);
// t[1][2]=1;
// printf("%lf\n",t[1][2]);

// return 0;

// }
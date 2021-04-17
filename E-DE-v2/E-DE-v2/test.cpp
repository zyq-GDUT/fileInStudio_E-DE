
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include "ede.h"
#include <time.h>
#include<math.h>
//#include <iostream>
//using namespace std;

//保存当前个体每个簇的均值向量，每一行都是对应一个簇的均值向量，若当前个体k个簇则一共是k行。
double** mean_cluster;
//保存当前个体每个簇对应的样本数
int* dataNum_Incluster;

//需要特殊处理的文件


/**
 * 计算当前个体每个簇的均值向量和每个簇对应的样本数
 * k_num是指总簇数
 * **/
double* count_meanV1(double* data, int* in_clu, int k_num, int N, int D) {
	double* meanV;//用来保存每个簇的中心坐标（均值向量）
	int* numInclu;//用来保存每个簇对应的样本数
	malloc1D(meanV, k_num * D);
	malloc1E(numInclu, k_num);
	memset(meanV, 0, sizeof(*meanV) * k_num * D);
	memset(numInclu, 0, sizeof(*numInclu) * k_num);
	for (int i = 0; i < N; i++) {
		int cn = in_clu[i];
		numInclu[cn]++;
		for (int j=0; j < D; j++) {
			*(meanV + cn * D + j) += *(data + i * D + j);
		}
	}
	//求均值向量
	for (int k = 0; k < k_num; k++) {
		for (int i = 0; i < D; i++) {
			*(meanV + k * D + i) = *(meanV + k * D + i) / numInclu[k];
		}
	}

	printf("\n输出均值向量\n");
	for (int i = 0; i < k_num * D; i++) {
		for (int j = 0; j < k_num; j++) {
			printf("%f", *(meanV + i * D + j));
		}
		printf("\t||\t");
		if (i % D == 0) {
			printf("\n");
		}

	}
	return meanV;
}

/**i指标函数
 * k_num是指当前个体的簇数
 * data保存了所有样本。
 * p_index保存了当前个体的所有簇中心坐标
 * in_clu记录了每个样本对应的簇下标
 * N是指样本数
 * D是指样本的维数
 * **/
double I_index1(double k_num, double* data, double* p_index, int* in_clu, int N, int D) {
	int i, j, k;
	double E = 0;
	//计算每个簇与簇中心的距离之和，即计算E值
	for (i = 0; i < N; i++) {
		for (j = 0; j < D; j++) {
			int cnum = in_clu[i];
			E += fabs(*(data + i * D + j) - *(p_index + cnum * D + j));
		}
	}
	//找出所有簇中距离最远的两个簇的距离。即计算Dmax值
	double Dmax = 0.0;
	double distance;

	for (i = 0; i  < k_num; i++) {
		for (j = i + 1; j < k_num; j++) {
			distance = 0;
			for (k = 0; k < D; k++) {
				distance += fabs(*(p_index + i * D + k) - *(p_index + j * D + k));
			}
			if (distance > Dmax) {
				Dmax = distance;
			}
		}
	}

	//计算I值
	double I_index = pow((1 / k_num) * (N / E) * Dmax, N);

	return I_index;
}




//DB指标函数
double DB_index1(double k_num, double* p_pars, double* p_index, int* in_clu, int N, int D) {

	return 0.1;
}


/*
统计文件中的样本数和维度数
filename是指要统计的文件名
n用来返回文件中的行数，即样本数
d用来返回样本的维度数
*/
void count_N_D1(const char* filename, int& n, int& d) {
	FILE* fp;
	char buffer[1000];
	char c;
	int bufferlen;
	int i;
	n = 0;
	d = 0;

	if ((fp = fopen(filename, "rb")) == NULL) {
		printf("文件不能打开\n");
		exit(0);
	}
	while (fgets(buffer, 1000, fp) != NULL) {
		bufferlen = strlen(buffer);
		//跳过空白行
		if (bufferlen == 2 && buffer[0] == 13 && buffer[1] == 10) {
			continue;
		}

		n++;
		//在文件的第二行统计当前数据集的维数
		if (n == 2) {
			for (i = 0; i < bufferlen; i++) {
				c = buffer[i];
				if (c == ',') {
					d++;
				}
			}
		}

	}
	fclose(fp);
	printf("行数：%d\n", n);
	printf("维数：%d\n", d);
}

double** loaddata2(const char* filename, int& n, int& d)
{
	int i, j;
	double** arraydata;
	char buffer[100];
	FILE* fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("file open error\n");
		exit(1);
	}
	count_N_D1(filename, n, d);

	malloc2E(arraydata, n, d);
	char* bu[100];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < d; j++) {
			const char* s = "%lf";
			fscanf(fp, "%lf", &arraydata[i][j]);
			printf("%lf\t", arraydata[i][j]);
			if (j == d - 1) {
				printf("\n");
				fgets(buffer, 100, fp);
			}
		}
		// fscanf(fp, "%lf,%lf,%lf,%lf", &a[i][0], &a[i][1], &a[i][2],&a[i][3]);
		// printf("%lf,%lf,%lf,%lf\n", a[i][0],a[i][1],a[i][2],a[i][3]);
	}
	fclose(fp);
	return arraydata;
}


/*
若当前行一个单词不是数字则跳过
返回字符串buffer忽略第一个单词后的第一个下标
注意这里只能识别用逗号，空格，制表符进行分开的字符串
*/

int skip_fw(char* buffer) {
	
	int t = 0;
	int flag = 0;
	while (buffer[t] != 44 && buffer[t] != 32 && buffer[t] != 9) {
		//有一个英文字母则认为是单词
		if (buffer[t] >= 65 && buffer[t] <= 90 || buffer[t] >= 97 && buffer[t] <= 122) {
			flag = 1;
		}
		t++;
	}
	//不需要跳过首个单词
	if (flag == 0) {
		t = 0;
	}
	return t;

}
/*
若当前行最后一个单词不是数字则跳过
返回字符串buffer忽略最后一个单词后的下标
注意这里只能识别用逗号，空格，制表符进行分开的字符串
*/
int skip_lw(char* buffer, int bufferlen) {
	int flag = 0;
	int w = bufferlen - 1;
	while (buffer[w] != 44 && buffer[w] != 32 && buffer[w] != 9) {
		if (buffer[w] >= 65 && buffer[w] <= 90 || buffer[w] >= 97 && buffer[w] <= 122) {
			flag = 1;
		}
		w--;
	}
	//不需要跳过最后一个单词
	if (flag == 0) {
		w = bufferlen;
	}
	return w;

}

/*
*跳过字符串最后一个单词，
*返回字符串buffer忽略最后一个单词后的下标
*注意这里只能识别用逗号，空格，制表符进行分开的字符串
*/
int skipLastWord(char* buffer, int bufferlen) {
	int w = bufferlen - 1;
	while (buffer[w] != 44 && buffer[w] != 32 && buffer[w] != 9) {
		w--;
	}
	return w;
}



/*
从字符串buffer的[t,w)范围区间中，读取数字型数据.
将读取到的数据保存到数组middle中，num是middle数组的长度
*/
void  readInNumberData(char* buffer,int bufferlen, int t, int w,double* middle,int &num) {
	char next = '\0';
	char c;
	int k = 0;
	char sk[20] = { 0 };
	for (int i = t; i < w; i++) {
		c = buffer[i];
		if (i + 1 < bufferlen) {
			next = buffer[i + 1];
		}
		if (c <= 57 && c >= 48 || c == 46) {
			sk[k] = c;
			k++;
			if (i + 1 >= bufferlen || next > 57 || (next < 48 && next != 46)) {
				// if(i+1>=bufferlen||next==32||next==44||next==13||next==10){
					// printf("%s\t",sk);
					 //printf("next=%c=%d\t",next,next);
				char str1[20] = { 0 };
				strncat(str1, sk, k);
				// printf("%s\t",str1);
				//printf("%lf\n", atof(str1));
				middle[num] = atof(str1);
				num++;
				k = 0;//sk字符串重新开始存数字
			}
		}

	}
}

/*
* 将一维数组转成n行num/n列的二维数组，并将其返回
* k是一维数组的个数
*/
double** arrConver(double* middle, int num, int n) {
	if (num % n != 0) {
		printf("行数不能整除个数，无法转换");
		exit(-1);
	}
	int d = num / n;
	double** arraydata;
	malloc2E(arraydata, n, d);
	printf("\n");
	for (int i = 0; i < n; i++) {
		printf("第%d行:\t", i + 1);
		for (int j = 0; j < d; j++) {
			arraydata[i][j] = middle[i * d + j];
			printf("%lf\t ", arraydata[i][j]);
		}
		printf("\n");
	}

	printf("个数：%d\n", num);

	
	printf("行数：%d\n", n);
	// printf("维数：%d\n",d);
	return arraydata;
	

}

//需要进行特殊处理的文件
char specialfile[10][20] = { "data.txt", "seeds.txt", "iris.txt" ,"test.txt" };

/*
若当前文件名在数组specialfile中，则返回1
*/
int checkfile(const char* filename) {
	int k = 0;
	int special = 0;
	while (*(specialfile[k]) != '\0') {
		if (!strcmp(specialfile[k], filename)) {
			//printf("specialfile[%d]=%s,filename=%s\n", k, specialfile[k], filename);
			special = 1;
			break;
		}

		k++;
	}
	return special;
}



double** loaddata3(const char* filename, int& n, int& d)
{
	FILE* fp;
	char buffer[1000];
	int bufferlen;
	//int i;
	int num = 0;
	n = 0,d=0;
	int t;
	int w;
	//int flag;
	//double middle[10000] = { 0 };//用于保存数据集中的数字
	double* middle;//用于保存数据集中的数字
	malloc1D(middle, 40000);//这种方式分配内存，与new方法一致，分配的都是动态内存，即堆内存。若使用数组初始化方式，则分配的是静态内存即栈内存。一般栈内存不能分配过大
	if ((fp = fopen(filename, "rb")) == NULL) {
		printf("文件不能打开\n");
		exit(0);
	}
	int special;//特殊处理标志位，由于有些数据集在文件中的簇名是数字，因此不能把簇名当做数据处理。
	/*int k = 0;
	while(*(s[k])!='\0'){
		if (!strcmp(s[k], filename)) {
			printf("s[%d]=%s,filename=%s\n",k,s[k],filename);
			special = 1;
			break;
		}

		k++;
	}*/

	special=checkfile(filename);

	while (fgets(buffer, 1000, fp) != NULL) {
		bufferlen = strlen(buffer);
		//跳过空白行
		if (bufferlen == 2 && buffer[0] == 13 && buffer[1] == 10) {
			continue;
		}
		n++;//计算非空行数，即可用数据的个数
	
		//t = 0;
		//flag = 0;
		//while(buffer[t]!=44&&buffer[t]!=32&&buffer[t]!=9){
		//	if(buffer[t]>=65&&buffer[t]<=90||buffer[t]>=97&&buffer[t]<=122){
		// 		flag=1;
		//	}
		//	t++;
		//}
		// //不需要跳过首个单词
		//if(flag==0){
		//	t=0;
		//}
		t = skip_fw(buffer);

		//flag=0;
		//w = bufferlen - 1;
		//while(buffer[w]!=44&&buffer[w]!=32&&buffer[w]!=9){
		//	if(buffer[w]>=65&&buffer[w]<=90||buffer[w]>=97&&buffer[w]<=122){
		// 		flag=1;
		//	}
		//	w--;
		//}
		////不需要跳过最后一个单词
		//if(flag==0){
		//	w=bufferlen;
		//}
		if (special == 1) {
			w=skipLastWord(buffer, bufferlen);
			//printf("bufferlen=%d,,w=%d\n",bufferlen, w);
		}else {
			w = skip_lw(buffer, bufferlen);
		}
		

		//char next='\0';
		//char c;
		//int k = 0;
		//char sk[20] = { 0 };
		//for (i = t; i < w ; i++) {
		//	c = buffer[i];
		//	if (i + 1 < bufferlen) {
		//		next = buffer[i + 1];
		//	}
		//	if (c <= 57 && c >= 48 || c == 46) {
		//		sk[k] = c;
		//		k++;
		//		if (i + 1 >= bufferlen || next > 57 || (next < 48 && next != 46)) {
		//			// if(i+1>=bufferlen||next==32||next==44||next==13||next==10){
		//				// printf("%s\t",sk);
		//				 //printf("next=%c=%d\t",next,next);
		//			char str1[20] = { 0 };
		//			strncat(str1, sk, k);
		//			// printf("%s\t",str1);
		//			//printf("%lf\n", atof(str1));
		//			middle[num] = atof(str1);
		//			num++;
		//			k = 0;//sk字符串重新开始存数字
		//		}
		//	}

		//}
		readInNumberData(buffer, bufferlen, t, w, middle, num);


	}
	
	if (num % n != 0) {
		printf("有些数据集的维度不一致\n");
		exit(-1);
	}
	d = num / n;//计算数据集的维度
	//printf("\n维度是%d", d);
	//double** arraydata;
	//malloc2E(arraydata, n, d);
	//printf("\n");
	//for (int i = 0; i < n; i++) {
	//	printf("第%d行:\t",i+1);
	//	for (int j = 0; j < d; j++) {
	//		arraydata[i][j]=middle[i * d + j];
	//		printf("%lf\t ", arraydata[i][j]);
	//	}
	//	printf("\n");
	//}

	//printf("个数：%d\n", num);

	//fclose(fp);
	//printf("行数：%d\n", n);
	//// printf("维数：%d\n",d);
	double** arraydata;
	arraydata=arrConver(middle, num, n);
	free(middle);//释放动态分配的内存
	fclose(fp);
	return arraydata;
}



/*
生成n个不同的随机数
随机数范围0-rang-1
*/
int* rand_generator1(int n, int rang) {

	if (n >= rang) {
		exit(1);
	}
	int* a;
	int* b;
	int rand_d;
	malloc1E(a, n);
	malloc1E(b, rang);
	for (int i = 0; i < rang; i++) {
		b[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		rand_d = (int)rand() % rang;
		while (b[rand_d] != 0) {
			rand_d = (int)rand() % rang;
		}
		b[rand_d] = 1;
		a[i] = rand_d;
	}
	return a;


}
int  main()
{
	// 	srand((unsigned int)(time(NULL)));
	// double** a;
	// int d;
	// int n;
	// 	int q,e;
	// 	const char *filename="1.txt";
	//     loaddata2(filename,q,e);
	// 	// count_N_D("iris.txt",q,e);
	// 	printf("\n n=%d,d=%d",q,e);

	// int n=10,rang=11;
	// int *a=rand_generator(n,rang);
	// for(int i=0;i<n;i++){
	// 	printf("%d\n",a[i]);
	// }
	// double f=fabs(-0.001);
	// printf("f=%lf",f);
		// int k_num=10;
		// int D=4;
		// double *meanV;//用来保存每个簇的中心坐标（均值向量）
		// int *numInclu;//用来保存每个簇对应的样本数
		// malloc1D(meanV,k_num*D);
		// malloc1E(numInclu,k_num);
		// memset(meanV,0,sizeof(*meanV)*k_num*D);
		// memset(numInclu,0,sizeof(*numInclu)*k_num);
		// for(int i=0;i<k_num;i++){
		// 	printf("meanV[%d]=%d\n",i,numInclu[i]);

		// }

	// 	 char str1[101] = { 0 };
	//     char str2[50] = { 0 };
	//    str2[0]='3';
	//     str2[1]='.';
	// 	str2[2]='1';
	// 	char d='d';
	// 	char *p;
	//     strncat(str1, str2,3);
	//     printf("%s",str1);
	//    printf("\n%lf",atof(str1));

	int q, e;
	const char* filename = "seeds.txt";
	loaddata3(filename, q, e);
	// const char *a="test";
	// int t=strcmp(a,"dddd");



	//char s[10][20] = { "data.txt", "seeds.data", "iris.txt" };
	//int k = 0;
	//int special=0;
	//const char* file = "data.txt";
	//while(*(s[k])!='\0'){
	//	//printf("%s\n", s[k]);
	//	if (!strcmp(s[k], file)) {
	//		printf("s[%d]=%s,filename=%s\n",k,s[k],file);
	//		special = 1;
	//		break;
	//	}
	//
	//	k++;
	//}
	//printf("special=%d\n", special);


	// printf("%d",t);
	 //count_N_D("iris.txt",q,e);
	 printf("\n n=%d,d=%d",q,e);

	return 0;
}
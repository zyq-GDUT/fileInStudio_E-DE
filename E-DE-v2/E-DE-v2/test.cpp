
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include "ede.h"
#include <time.h>
#include<math.h>
//#include <iostream>
//using namespace std;

//���浱ǰ����ÿ���صľ�ֵ������ÿһ�ж��Ƕ�Ӧһ���صľ�ֵ����������ǰ����k������һ����k�С�
double** mean_cluster;
//���浱ǰ����ÿ���ض�Ӧ��������
int* dataNum_Incluster;

//��Ҫ���⴦����ļ�


/**
 * ���㵱ǰ����ÿ���صľ�ֵ������ÿ���ض�Ӧ��������
 * k_num��ָ�ܴ���
 * **/
double* count_meanV1(double* data, int* in_clu, int k_num, int N, int D) {
	double* meanV;//��������ÿ���ص��������꣨��ֵ������
	int* numInclu;//��������ÿ���ض�Ӧ��������
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
	//���ֵ����
	for (int k = 0; k < k_num; k++) {
		for (int i = 0; i < D; i++) {
			*(meanV + k * D + i) = *(meanV + k * D + i) / numInclu[k];
		}
	}

	printf("\n�����ֵ����\n");
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

/**iָ�꺯��
 * k_num��ָ��ǰ����Ĵ���
 * data����������������
 * p_index�����˵�ǰ��������д���������
 * in_clu��¼��ÿ��������Ӧ�Ĵ��±�
 * N��ָ������
 * D��ָ������ά��
 * **/
double I_index1(double k_num, double* data, double* p_index, int* in_clu, int N, int D) {
	int i, j, k;
	double E = 0;
	//����ÿ����������ĵľ���֮�ͣ�������Eֵ
	for (i = 0; i < N; i++) {
		for (j = 0; j < D; j++) {
			int cnum = in_clu[i];
			E += fabs(*(data + i * D + j) - *(p_index + cnum * D + j));
		}
	}
	//�ҳ����д��о�����Զ�������صľ��롣������Dmaxֵ
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

	//����Iֵ
	double I_index = pow((1 / k_num) * (N / E) * Dmax, N);

	return I_index;
}




//DBָ�꺯��
double DB_index1(double k_num, double* p_pars, double* p_index, int* in_clu, int N, int D) {

	return 0.1;
}


/*
ͳ���ļ��е���������ά����
filename��ָҪͳ�Ƶ��ļ���
n���������ļ��е���������������
d��������������ά����
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
		printf("�ļ����ܴ�\n");
		exit(0);
	}
	while (fgets(buffer, 1000, fp) != NULL) {
		bufferlen = strlen(buffer);
		//�����հ���
		if (bufferlen == 2 && buffer[0] == 13 && buffer[1] == 10) {
			continue;
		}

		n++;
		//���ļ��ĵڶ���ͳ�Ƶ�ǰ���ݼ���ά��
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
	printf("������%d\n", n);
	printf("ά����%d\n", d);
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
����ǰ��һ�����ʲ�������������
�����ַ���buffer���Ե�һ�����ʺ�ĵ�һ���±�
ע������ֻ��ʶ���ö��ţ��ո��Ʊ�����зֿ����ַ���
*/

int skip_fw(char* buffer) {
	
	int t = 0;
	int flag = 0;
	while (buffer[t] != 44 && buffer[t] != 32 && buffer[t] != 9) {
		//��һ��Ӣ����ĸ����Ϊ�ǵ���
		if (buffer[t] >= 65 && buffer[t] <= 90 || buffer[t] >= 97 && buffer[t] <= 122) {
			flag = 1;
		}
		t++;
	}
	//����Ҫ�����׸�����
	if (flag == 0) {
		t = 0;
	}
	return t;

}
/*
����ǰ�����һ�����ʲ�������������
�����ַ���buffer�������һ�����ʺ���±�
ע������ֻ��ʶ���ö��ţ��ո��Ʊ�����зֿ����ַ���
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
	//����Ҫ�������һ������
	if (flag == 0) {
		w = bufferlen;
	}
	return w;

}

/*
*�����ַ������һ�����ʣ�
*�����ַ���buffer�������һ�����ʺ���±�
*ע������ֻ��ʶ���ö��ţ��ո��Ʊ�����зֿ����ַ���
*/
int skipLastWord(char* buffer, int bufferlen) {
	int w = bufferlen - 1;
	while (buffer[w] != 44 && buffer[w] != 32 && buffer[w] != 9) {
		w--;
	}
	return w;
}



/*
���ַ���buffer��[t,w)��Χ�����У���ȡ����������.
����ȡ�������ݱ��浽����middle�У�num��middle����ĳ���
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
				k = 0;//sk�ַ������¿�ʼ������
			}
		}

	}
}

/*
* ��һά����ת��n��num/n�еĶ�ά���飬�����䷵��
* k��һά����ĸ���
*/
double** arrConver(double* middle, int num, int n) {
	if (num % n != 0) {
		printf("�������������������޷�ת��");
		exit(-1);
	}
	int d = num / n;
	double** arraydata;
	malloc2E(arraydata, n, d);
	printf("\n");
	for (int i = 0; i < n; i++) {
		printf("��%d��:\t", i + 1);
		for (int j = 0; j < d; j++) {
			arraydata[i][j] = middle[i * d + j];
			printf("%lf\t ", arraydata[i][j]);
		}
		printf("\n");
	}

	printf("������%d\n", num);

	
	printf("������%d\n", n);
	// printf("ά����%d\n",d);
	return arraydata;
	

}

//��Ҫ�������⴦����ļ�
char specialfile[10][20] = { "data.txt", "seeds.txt", "iris.txt" ,"test.txt" };

/*
����ǰ�ļ���������specialfile�У��򷵻�1
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
	//double middle[10000] = { 0 };//���ڱ������ݼ��е�����
	double* middle;//���ڱ������ݼ��е�����
	malloc1D(middle, 40000);//���ַ�ʽ�����ڴ棬��new����һ�£�����Ķ��Ƕ�̬�ڴ棬�����ڴ档��ʹ�������ʼ����ʽ���������Ǿ�̬�ڴ漴ջ�ڴ档һ��ջ�ڴ治�ܷ������
	if ((fp = fopen(filename, "rb")) == NULL) {
		printf("�ļ����ܴ�\n");
		exit(0);
	}
	int special;//���⴦���־λ��������Щ���ݼ����ļ��еĴ��������֣���˲��ܰѴ����������ݴ���
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
		//�����հ���
		if (bufferlen == 2 && buffer[0] == 13 && buffer[1] == 10) {
			continue;
		}
		n++;//����ǿ����������������ݵĸ���
	
		//t = 0;
		//flag = 0;
		//while(buffer[t]!=44&&buffer[t]!=32&&buffer[t]!=9){
		//	if(buffer[t]>=65&&buffer[t]<=90||buffer[t]>=97&&buffer[t]<=122){
		// 		flag=1;
		//	}
		//	t++;
		//}
		// //����Ҫ�����׸�����
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
		////����Ҫ�������һ������
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
		//			k = 0;//sk�ַ������¿�ʼ������
		//		}
		//	}

		//}
		readInNumberData(buffer, bufferlen, t, w, middle, num);


	}
	
	if (num % n != 0) {
		printf("��Щ���ݼ���ά�Ȳ�һ��\n");
		exit(-1);
	}
	d = num / n;//�������ݼ���ά��
	//printf("\nά����%d", d);
	//double** arraydata;
	//malloc2E(arraydata, n, d);
	//printf("\n");
	//for (int i = 0; i < n; i++) {
	//	printf("��%d��:\t",i+1);
	//	for (int j = 0; j < d; j++) {
	//		arraydata[i][j]=middle[i * d + j];
	//		printf("%lf\t ", arraydata[i][j]);
	//	}
	//	printf("\n");
	//}

	//printf("������%d\n", num);

	//fclose(fp);
	//printf("������%d\n", n);
	//// printf("ά����%d\n",d);
	double** arraydata;
	arraydata=arrConver(middle, num, n);
	free(middle);//�ͷŶ�̬������ڴ�
	fclose(fp);
	return arraydata;
}



/*
����n����ͬ�������
�������Χ0-rang-1
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
		// double *meanV;//��������ÿ���ص��������꣨��ֵ������
		// int *numInclu;//��������ÿ���ض�Ӧ��������
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
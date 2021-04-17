
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include "malloc_1.h"
#include <time.h>
//#include <iostream>
//using namespace std;

//保存当前个体每个簇的均值向量，每一行都是对应一个簇的均值向量，若当前个体k个簇则一共是k行。
double **mean_cluster;
//保存当前个体每个簇对应的样本数
int* dataNum_Incluster;

/**
 * 计算当前个体每个簇的均值向量和每个簇对应的样本数
 * k_num是指总簇数
 * **/
double* count_meanV(double* data,int *in_clu,int k_num,int N,int D){
	double *meanV;//用来保存每个簇的中心坐标（均值向量）
	int *numInclu;//用来保存每个簇对应的样本数
	malloc1D(meanV,k_num*D);
	malloc1E(numInclu,k_num);
	memset(meanV,0,sizeof(*meanV)*k_num*D);
	memset(numInclu,0,sizeof(*numInclu)*k_num);
	for(int i=0;i<N;i++){
		int cn=in_clu[i];
		numInclu[cn]++;
		for(int j;j<D;j++){
			*(meanV+cn*D+j)+=*(data+i*D+j);
		}
	}
	//求均值向量
	for(int k=0;k<k_num;k++){
		for(int i=0;i<D;i++){
			*(meanV+k*D+i)=*(meanV+k*D+i)/numInclu[k];
		}
	}

	printf("\n输出均值向量\n");
	for(int i=0;i<k_num*D;i++){
		for(int j=0;j<k_num;j++){
			printf("%f",*(meanV+i*D+j));
		}
		printf("\t||\t");
		if(i%D==0){
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
double I_index(double k_num, double *data, double *p_index, int *in_clu, int N, int D){
	int i,j,k;
	double E=0;
	//计算每个簇与簇中心的距离之和，即计算E值
	for(i=0;i<N;i++){
		for( j=0;j<D;j++){
			int cnum=in_clu[i];
			E+=fabs(*(data+i*D+j)-*(p_index+cnum*D+j));
		}
	}
	//找出所有簇中距离最远的两个簇的距离。即计算Dmax值
	double Dmax=0.0;
	double distance;
	
	for(i=0;i<k_num;i++){
		for(j=i+1;j<k_num;j++){
				distance=0;
				for(k=0;k<D;k++){
						distance+=fabs(*(p_index+i*D+k)-*(p_index+j*D+k));
				}
				if(distance>Dmax){
					Dmax=distance;
				}
		}
	}

	//计算I值
	double I_index=pow((1/k_num)*(N/E)*Dmax,N);

	return I_index;
}


double I_index_test( double *data){
	for(int i=0;i<4;i++){
		for(int j=0;j<3;j++){
			printf("输出:%lf\t",*(data+i*3+j));
			
		}
		printf("\n");
		}

}

//DB指标函数
double DB_index(double k_num, double *p_pars, double *p_index, int *in_clu, int N, int D){


}


/*
统计文件中的样本数和维度数
filename是指要统计的文件名
n用来返回文件中的行数，即样本数
d用来返回样本的维度数
*/
void count_N_D(const char *filename,int &n,int &d){
	FILE *fp;
	char buffer[1000];
	char c;
	int bufferlen;
	int i;
	n=0;
	d=0;
 
	if((fp=fopen(filename,"rb"))==NULL){
		printf("文件不能打开\n");
		exit(0);
	}
	while(fgets(buffer,1000,fp)!=NULL){
		bufferlen=strlen(buffer);
		//跳过空白行
		if(bufferlen==2&&buffer[0]==13&&buffer[1]==10){
			continue;
		}
		
		n++;
		//在文件的第二行统计当前数据集的维数
		if(n==2){
			for(i=0;i<bufferlen;i++){
				c=buffer[i];
				if(c==','){
					d++;
				}
			}
		}
		
	}
	fclose(fp);
	printf("行数：%d\n",n);
	printf("维数：%d\n",d);
}

double **loaddata2(const char* filename,int &n,int &d)
{
	int i, j;
	double **arraydata; 
	char buffer[100];
 	FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("file open error\n");
        exit(1);
    }
	 count_N_D(filename,n,d);

    malloc2E(arraydata,n,d);
	char *bu[100];
    for (int i = 0; i < n; i++)
    {
		for(int j=0;j<d;j++){
			const char *s="%lf";
			fscanf(fp,"%lf",&arraydata[i][j]);
			printf("%lf\t",arraydata[i][j]);
			if(j==d-1){
				printf("\n");
				fgets(buffer,100,fp);
			}
		}
        // fscanf(fp, "%lf,%lf,%lf,%lf", &a[i][0], &a[i][1], &a[i][2],&a[i][3]);
        // printf("%lf,%lf,%lf,%lf\n", a[i][0],a[i][1],a[i][2],a[i][3]);
    }
    fclose(fp);
	return arraydata;
}

double **loaddata3(const char* filename,int &n,int &d)
{
	FILE *fp;
	char buffer[1000];
	int bufferlen;
	int i;
	n=0;
	d=0;
 	int num=0;
	int t;
	int w;
	int flag;
	if((fp=fopen(filename,"rb"))==NULL){
		printf("文件不能打开\n");
		exit(0);
	}
	while(fgets(buffer,1000,fp)!=NULL){
		bufferlen=strlen(buffer);
		//跳过空白行
		if(bufferlen==2&&buffer[0]==13&&buffer[1]==10){
			continue;
		}
		n++;
		// //如果当前行第一个字符不是数字，则忽略该行的第一个单词
		// if(buffer[0]>57||buffer[0]<48){
		// 		  t=0;
		// 		for(int j=0;j<bufferlen;j++){
		// 			//是逗号，制表，空格则退出循环
		// 			if(buffer[j]==44||buffer[j]==32||buffer[j]==9){
		// 				break;
		// 			}
		// 			t++;
					
		// 		}
		// }
		   t=0;
		   flag=0;
		// while(buffer[t]!=44||buffer[t]!=32||buffer[t]!=9){
		// 	if(buffer[t]>=65&&buffer[t]<=90||buffer[t]>=97&&buffer[t]<=122){
		// 		flag=1;
		// 	}
		// 	t++;
		// }
		// //不需要跳过首个单词
		// if(flag==0){
		// 	t=0;
		// }

		// flag=0;
		 w=bufferlen-1;
		// while(buffer[bufferlen-1]!=44||buffer[bufferlen-1]!=32||buffer[bufferlen-1]!=9){
		// 	if(buffer[w]>=65&&buffer[w]<=90||buffer[w]>=97&&buffer[w]<=122){
		// 		flag=1;
		// 	}
		// 	w--;
		// }
		// //不需要跳过最后一个单词
		// if(flag==0){
		// 	w=bufferlen;
		// }

		char next;
		char c;
		int k=0;
		char sk[20] = { 0 };
		for(i=t;i<w+1;i++){
			c=buffer[i];
			if(i+1<bufferlen){
				next=buffer[i+1];
			}
			if(c<=57&&c>=48||c==46){		
				sk[k]=c;
				k++;	
				if(i+1>=bufferlen||next>57||(next<48&&next!=46)){
				// if(i+1>=bufferlen||next==32||next==44||next==13||next==10){
					// printf("%s\t",sk);
					 //printf("next=%c=%d\t",next,next);
 					char str1[20] = { 0 };
 					strncat(str1, sk,k);
					// printf("%s\t",str1);
					printf("%lf\n",atof(str1));
					num++;
					k=0;
				}
			}
		
	   }


	}
			printf("个数数：%d\n",num);

	fclose(fp);
	printf("行数：%d\n",n);
	// printf("维数：%d\n",d);
}



/*
生成n个不同的随机数
随机数范围0-rang-1
*/
int* rand_generator(int n,int rang){

	if(n>=rang){
		exit(1);
	}
	int *a;
	int *b;
	int rand_d;
	malloc1E(a,n);
	malloc1E(b,rang);
	for(int i=0;i<rang;i++){
		b[i]=0;
	}
	for(int i=0;i<n;i++){
		rand_d=(int)rand()%rang;
		while(b[rand_d]!=0){
			rand_d=(int)rand()%rang;
		}
		b[rand_d]=1;
		a[i]=rand_d;
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

	int q,e;
	const char *filename="zoo.data";
    loaddata3(filename,q,e);
	// const char *a="test";
	// int t=strcmp(a,"dddd");

	// printf("%d",t);
	// count_N_D("iris.txt",q,e);
	// printf("\n n=%d,d=%d",q,e);
   return 0;
}
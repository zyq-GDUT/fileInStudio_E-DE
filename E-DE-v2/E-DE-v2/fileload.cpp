#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include"ede.h"
/*****************************************************************
* func_name: load the testing data as you need
* input: dimension and data number
* descript: load data from the txt file
*d是指数据的维数，n是指数据样本的总数
*****************************************************************/

//需要进行特殊处理的文件，即如果文件中的数据，每行内容最后的一个数字是表示簇编号的的话即作为特殊文件进行处理
//处理特殊文件时，读取样本数据时不会把簇编号当做样本的坐标。如果是正常文件，即会把每行内容的所有数字都作为样本的坐标来处理
char specialfile[10][20] = {  "seeds.txt"  };


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
void  readInNumberData(char* buffer, int bufferlen, int t, int w, double* middle, int& num) {
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
			
				char str1[20] = { 0 };
				strncat(str1, sk, k);
				 //printf("%s\t",str1);
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
	//printf("\n");
	for (int i = 0; i < n; i++) {
		//printf("第%d行:\t", i + 1);
		for (int j = 0; j < d; j++) {
			arraydata[i][j] = middle[i * d + j];
			//printf("%lf\t ", arraydata[i][j]);
		}
		//printf("\n");
	}



	return arraydata;


}


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


double** loadData2(const char* filename, int& n, int& d)
{
	FILE* fp;
	char buffer[1000];
	int bufferlen;
	int num = 0;
	n = 0, d = 0;
	int t;
	int w;
	double* middle;//用于保存数据集中的数字
	int special;//特殊处理标志位，由于有些数据集在文件中的簇名是数字，因此不能把簇名当做数据处理。
	malloc1D(middle, 40000);//这种方式分配内存，与new方法一致，分配的都是动态内存，即堆内存。若使用数组初始化方式，则分配的是静态内存即栈内存。一般栈内存不能分配过大
	if ((fp = fopen(filename, "rb")) == NULL) {
		printf("文件不能打开\n");
		exit(0);
	}
	special = checkfile(filename);

	while (fgets(buffer, 1000, fp) != NULL) {
		bufferlen = strlen(buffer);
		//跳过空白行
		if (bufferlen == 2 && buffer[0] == 13 && buffer[1] == 10) {
			continue;
		}
		n++;//计算非空行数，即可用数据的个数

		t = skip_fw(buffer);

		if (special == 1) {
			w = skipLastWord(buffer, bufferlen);
			//printf("bufferlen=%d,,w=%d\n",bufferlen, w);
		}
		else {
			w = skip_lw(buffer, bufferlen);
		}
		readInNumberData(buffer, bufferlen, t, w, middle, num);

	}

	if (num % n != 0) {
		printf("有些数据集的维度不一致\n");
		exit(-1);
	}
	d = num / n;//计算数据集的维度
	double** arraydata;
	arraydata = arrConver(middle, num, n);
	free(middle);//释放动态分配的内存
	fclose(fp);
	printf("\n行数：%d\n", n);
	 printf("维数：%d\n",d);
	 printf("个数：%d\n", num);
	return arraydata;
}





/*
* 该方法返回记录了数据集本身的聚类情况的数组
* filename表示数据集文件的文件名，N表示该文件的行数，即样本个数
* 所加载的数据集，每一行第一个字符都是数字表示该样本所在的簇编号（数据集本身已经分好类）
*/
int* getDataOriginClusterInfo(const char* filename,int N) {
	FILE* fp;
	char buffer[1000];
	int bufferlen;
	int* originClusterInfo;
	int lineNum = 0;//用于记录文件的行数
	int c;
	int clusterNum=0;//用来记录簇数，样本中的簇编号都是从1开始依次递增的，因此最大的簇编号就是簇数
	malloc1E(originClusterInfo, N+1);//最后一位用来保存簇编号
	if ((fp = fopen(filename, "rb")) == NULL) {
		printf("文件不能打开\n");
		exit(0);
	}
	//读取每行的内容到缓冲区buffer
	while (fgets(buffer, 1000, fp) != NULL) {
		bufferlen = strlen(buffer);
		//跳过空白行
		if (bufferlen == 2 && buffer[0] == 13 && buffer[1] == 10) {
			continue;
		}
		lineNum++;
		//不是如果当前行的第一个字符不是数字则无法处理
		if (buffer[0] > 57 || buffer[0] < 48) {
			printf("文件%s第%d行没有在行头标出簇编号！\n", filename, lineNum);
			exit(2);
		}
		 c=atoi(&buffer[0]);
		 //找最大的簇编号
		 if (c > clusterNum) {
			 clusterNum = c;
		 }
		 originClusterInfo[lineNum - 1] = c-1;//簇编号从0开始
		 //printf("第%d个样本的簇编号为%d\n", lineNum, originClusterInfo[lineNum-1]);
	}
	originClusterInfo[N] = clusterNum;//在数组的最后保存总簇数
	if (N != lineNum) {
		printf("文件%s数据有误！", filename);
		exit(1);
	}
	fclose(fp);
	return originClusterInfo;
}


/*
统计文件中的样本数和维度数
filename是指要统计的文件名
n用来返回文件中的行数，即样本数
d用来返回样本的维度数
*/
void count_N_D(const char* filename, int& n, int& d) {
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

double** loadData(const char* filename, int& n, int& d)
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
	count_N_D(filename, n, d);

	malloc2E(arraydata, n, d);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < d; j++) {

			fscanf(fp, "%lf,", &arraydata[i][j]);
			// printf("%lf\t",arraydata[i][j]);

			if (j == d - 1) {
				// printf("\n");
				fgets(buffer, 100, fp);
			}
		}
		// fscanf(fp, "%lf,%lf,%lf,%lf", &a[i][0], &a[i][1], &a[i][2],&a[i][3]);
		// printf("%lf,%lf,%lf,%lf\n", a[i][0],a[i][1],a[i][2],a[i][3]);
	}
	fclose(fp);
	return arraydata;
}





//int main() {
//	int q, e;
//	const char* filename = "seeds.txt";
//	loadData2(filename, q, e);
//	getDataOriginClusterInfo(filename, q);
//
//
//
//}





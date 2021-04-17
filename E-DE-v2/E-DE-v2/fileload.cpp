#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include"ede.h"
/*****************************************************************
* func_name: load the testing data as you need
* input: dimension and data number
* descript: load data from the txt file
*d��ָ���ݵ�ά����n��ָ��������������
*****************************************************************/

//��Ҫ�������⴦����ļ���������ļ��е����ݣ�ÿ����������һ�������Ǳ�ʾ�ر�ŵĵĻ�����Ϊ�����ļ����д���
//���������ļ�ʱ����ȡ��������ʱ����Ѵر�ŵ������������ꡣ����������ļ��������ÿ�����ݵ��������ֶ���Ϊ����������������
char specialfile[10][20] = {  "seeds.txt"  };


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
	//printf("\n");
	for (int i = 0; i < n; i++) {
		//printf("��%d��:\t", i + 1);
		for (int j = 0; j < d; j++) {
			arraydata[i][j] = middle[i * d + j];
			//printf("%lf\t ", arraydata[i][j]);
		}
		//printf("\n");
	}



	return arraydata;


}


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


double** loadData2(const char* filename, int& n, int& d)
{
	FILE* fp;
	char buffer[1000];
	int bufferlen;
	int num = 0;
	n = 0, d = 0;
	int t;
	int w;
	double* middle;//���ڱ������ݼ��е�����
	int special;//���⴦���־λ��������Щ���ݼ����ļ��еĴ��������֣���˲��ܰѴ����������ݴ���
	malloc1D(middle, 40000);//���ַ�ʽ�����ڴ棬��new����һ�£�����Ķ��Ƕ�̬�ڴ棬�����ڴ档��ʹ�������ʼ����ʽ���������Ǿ�̬�ڴ漴ջ�ڴ档һ��ջ�ڴ治�ܷ������
	if ((fp = fopen(filename, "rb")) == NULL) {
		printf("�ļ����ܴ�\n");
		exit(0);
	}
	special = checkfile(filename);

	while (fgets(buffer, 1000, fp) != NULL) {
		bufferlen = strlen(buffer);
		//�����հ���
		if (bufferlen == 2 && buffer[0] == 13 && buffer[1] == 10) {
			continue;
		}
		n++;//����ǿ����������������ݵĸ���

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
		printf("��Щ���ݼ���ά�Ȳ�һ��\n");
		exit(-1);
	}
	d = num / n;//�������ݼ���ά��
	double** arraydata;
	arraydata = arrConver(middle, num, n);
	free(middle);//�ͷŶ�̬������ڴ�
	fclose(fp);
	printf("\n������%d\n", n);
	 printf("ά����%d\n",d);
	 printf("������%d\n", num);
	return arraydata;
}





/*
* �÷������ؼ�¼�����ݼ�����ľ������������
* filename��ʾ���ݼ��ļ����ļ�����N��ʾ���ļ�������������������
* �����ص����ݼ���ÿһ�е�һ���ַ��������ֱ�ʾ���������ڵĴر�ţ����ݼ������Ѿ��ֺ��ࣩ
*/
int* getDataOriginClusterInfo(const char* filename,int N) {
	FILE* fp;
	char buffer[1000];
	int bufferlen;
	int* originClusterInfo;
	int lineNum = 0;//���ڼ�¼�ļ�������
	int c;
	int clusterNum=0;//������¼�����������еĴر�Ŷ��Ǵ�1��ʼ���ε����ģ�������Ĵر�ž��Ǵ���
	malloc1E(originClusterInfo, N+1);//���һλ��������ر��
	if ((fp = fopen(filename, "rb")) == NULL) {
		printf("�ļ����ܴ�\n");
		exit(0);
	}
	//��ȡÿ�е����ݵ�������buffer
	while (fgets(buffer, 1000, fp) != NULL) {
		bufferlen = strlen(buffer);
		//�����հ���
		if (bufferlen == 2 && buffer[0] == 13 && buffer[1] == 10) {
			continue;
		}
		lineNum++;
		//���������ǰ�еĵ�һ���ַ������������޷�����
		if (buffer[0] > 57 || buffer[0] < 48) {
			printf("�ļ�%s��%d��û������ͷ����ر�ţ�\n", filename, lineNum);
			exit(2);
		}
		 c=atoi(&buffer[0]);
		 //�����Ĵر��
		 if (c > clusterNum) {
			 clusterNum = c;
		 }
		 originClusterInfo[lineNum - 1] = c-1;//�ر�Ŵ�0��ʼ
		 //printf("��%d�������Ĵر��Ϊ%d\n", lineNum, originClusterInfo[lineNum-1]);
	}
	originClusterInfo[N] = clusterNum;//���������󱣴��ܴ���
	if (N != lineNum) {
		printf("�ļ�%s��������", filename);
		exit(1);
	}
	fclose(fp);
	return originClusterInfo;
}


/*
ͳ���ļ��е���������ά����
filename��ָҪͳ�Ƶ��ļ���
n���������ļ��е���������������
d��������������ά����
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





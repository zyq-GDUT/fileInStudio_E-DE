#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include <float.h>
#include "ede.h"


/*
* ��ȫ�ֱ���traildata������ֵ
*/
void set_traildataInfo(trailData& traildata, int correctClusterNum, const char* filename, int usedindex) {
	if (sizeof(traildata) <= 0 || correctClusterNum <= 0 || filename==NULL || usedindex <0) {
		printf("����traildata����ʧ��\n");
		exit(-1);
	}
	traildata.correctClusterNum = correctClusterNum;
	traildata.filename = filename;
	traildata.usedindex = usedindex;
}

/*
* ���ظ�����
* D��ʾ���ݵ�ά��
* timesά���ı�����ά�����Ը�ֵ��Ϊ������
*/
int getNP(int D, int times) {
	int NP = 0;
	if (D < 2 || times < 1) {
		exit(-1);
	}
	NP = D * times;
	return NP;
}

/*
* ���ص�������
* NP��ʾ������
* divisor/NP��Ϊ��������
*/
int getGmax(int NP,int divisor) {
	if (NP > divisor) {
		exit(-1);
	}
	int Gmax = 0;
	Gmax = 1000000 / NP;
	return Gmax;
}
/*
* ��ʼ�����Ͷ�ά����
* x������y����
* value��ʼ����ֵ
*/
void initial_two_Dim_intarr(int** arr,int x,int y,int value) {
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			arr[i][j] = value;
		}
	}
}

/*
* �ҵ����ݼ���ÿά�ȵ����ֵ����Сֵ
* ����uk����ÿά�ȵ����ֵ
* ����lk����ÿά�ȵ���Сֵ
* data������ÿ�����������ݣ������ݼ�
*/
void getDataDim_max_min(double** data, double* uk, double* lk, int D,int N) {
	double data_min, data_max;
	int i, j;
	for (i = 0; i < D; i++) {
		data_min = DBL_MAX, data_max = DBL_MIN;
		for (j = 0; j < N; j++) {
			if (data[j][i] > data_max)
				data_max = data[j][i];
			if (data[j][i] < data_min)
				data_min = data[j][i];
		}

		uk[i] = data_max;
		lk[i] = data_min;
	}
}




void initial_cluster_center(double** data, double* popul_index, double* popul_param, int N, int D, int max_cluster_num, int min_cluster_num) {
	int k_rand = rand() % (max_cluster_num - min_cluster_num + 1) + min_cluster_num;   //the cluster num is between[minclusternum,maxclusternum]
	popul_param[0] = (double)k_rand; //popul_paramÿһ�ж������У���һ���Ǳ�����ж�Ӧ�ĴصĴ������ڶ��ж�Ӧ���Ǹ��ж�Ӧ�Ĵص���Ӧ��
	int* rand_arr = rand_generator(k_rand, N);
	for (int j = 0; j < k_rand; j++) {
		for (int k = 0; k < D; k++) {
			// popul_index[i][j*D + k] = lk[k] + URAND *(uk[k] - lk[k]);//��ʼ�������������k_rand���أ�������ÿ���ص�ÿһά�ȶ������ʼ��һ��ȡֵ��Χ�ڵ�ֵ
			popul_index[j * D + k] = data[rand_arr[j]][k];
		}
	}
	free(rand_arr);
}

void initial_incluster(double** data, double* popul_index, double* popul_param, int* in_cluster, int N, int D) {
	for (int j = 0; j < N; j++) {

		double min = DBL_MAX;

		for (int l = 0; l < ((int)popul_param[0]); l++) {
			double distance = getDistance(data[j], &popul_index[l * D], D);
			// printf("\n����%d���%d�ľ���Ϊ%lf",j,l,distance);
			if (distance < min) {
				min = distance;
				in_cluster[j] = l;    //record the data index��������¼��ǰ�����ڵ�ǰ������±�Ϊl�Ĵ��С�
			}
		}

	}

}

/*
* ��ʼ�����壬����ʼ������Ĵ�����
* ͬʱ��Ҫ��ʼ����������Ӧ�Ĵ����;������
* popul_index���浱ǰ����Ĵ�����
* popul_param�������浱ǰ����Ĵ�������Ӧ��
* in_cluster�������浱ǰ����ľ����������ÿ�����������Ĵر��
*/
void initial_individual(double** data,double* popul_index,double* popul_param,int* in_cluster,int N,int D,int max_cluster_num,int min_cluster_num) {
	//int k_rand = rand() % (max_cluster_num - min_cluster_num + 1) + min_cluster_num;   //the cluster num is between[minclusternum,maxclusternum]
	//popul_param[0] = (double)k_rand; //popul_paramÿһ�ж������У���һ���Ǳ�����ж�Ӧ�ĴصĴ������ڶ��ж�Ӧ���Ǹ��ж�Ӧ�Ĵص���Ӧ��
	//int* rand_arr = rand_generator(k_rand, N);
	//for (int j = 0; j < k_rand; j++) {
	//	for (int k = 0; k < D; k++) {
	//		// popul_index[i][j*D + k] = lk[k] + URAND *(uk[k] - lk[k]);//��ʼ�������������k_rand���أ�������ÿ���ص�ÿһά�ȶ������ʼ��һ��ȡֵ��Χ�ڵ�ֵ
	//		popul_index[j * D + k] = data[rand_arr[j]][k];
	//	}
	//}
	initial_cluster_center( data,popul_index,  popul_param,  N,  D, max_cluster_num,  min_cluster_num);

	//evaluate the fitness of the initial population
	//for (int j = 0; j < N; j++) {

	//	double min = DBL_MAX;

	//	for (int l = 0; l < ((int)popul_param[0]); l++) {
	//		double distance = getDistance(data[j], &popul_index[l * D], D);
	//		// printf("\n����%d���%d�ľ���Ϊ%lf",j,l,distance);
	//		if (distance < min) {
	//			min = distance;
	//			in_cluster[j] = l;    //record the data index��������¼��ǰ�����ڵ�ǰ������±�Ϊl�Ĵ��С�
	//		}
	//	}

	//}
	initial_incluster(data, popul_index, popul_param, in_cluster, N, D);
	dealwith_emptyCluster(in_cluster, N, popul_param[0], data, D, popul_index);
}

#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include <float.h>
#include "ede.h"


/*
* 给全局变量traildata赋部分值
*/
void set_traildataInfo(trailData& traildata, int correctClusterNum, const char* filename, int usedindex) {
	if (sizeof(traildata) <= 0 || correctClusterNum <= 0 || filename==NULL || usedindex <0) {
		printf("保存traildata数据失败\n");
		exit(-1);
	}
	traildata.correctClusterNum = correctClusterNum;
	traildata.filename = filename;
	traildata.usedindex = usedindex;
}

/*
* 返回个体数
* D表示数据的维数
* times维数的倍数，维数乘以该值即为个体数
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
* 返回迭代代数
* NP表示个体数
* divisor/NP即为迭代代数
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
* 初始化整型二维数组
* x行数，y列数
* value初始化的值
*/
void initial_two_Dim_intarr(int** arr,int x,int y,int value) {
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			arr[i][j] = value;
		}
	}
}

/*
* 找到数据集中每维度的最大值和最小值
* 数组uk保存每维度的最大值
* 数组lk保存每维度的最小值
* data保存了每个样本的数据，即数据集
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
	popul_param[0] = (double)k_rand; //popul_param每一行都有两列，第一列是保存该行对应的簇的簇数，第二列对应的是该行对应的簇的适应度
	int* rand_arr = rand_generator(k_rand, N);
	for (int j = 0; j < k_rand; j++) {
		for (int k = 0; k < D; k++) {
			// popul_index[i][j*D + k] = lk[k] + URAND *(uk[k] - lk[k]);//初始化操作，随机了k_rand个簇，这里是每个簇的每一维度都随机初始化一个取值范围内的值
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
			// printf("\n样本%d与簇%d的距离为%lf",j,l,distance);
			if (distance < min) {
				min = distance;
				in_cluster[j] = l;    //record the data index，用来记录当前样本在当前个体的下标为l的簇中。
			}
		}

	}

}

/*
* 初始化个体，即初始化个体的簇中心
* 同时还要初始化个体所对应的簇数和聚类情况
* popul_index保存当前个体的簇中心
* popul_param用来保存当前个体的簇数和适应度
* in_cluster用来保存当前个体的聚类情况，即每个样本所属的簇编号
*/
void initial_individual(double** data,double* popul_index,double* popul_param,int* in_cluster,int N,int D,int max_cluster_num,int min_cluster_num) {
	//int k_rand = rand() % (max_cluster_num - min_cluster_num + 1) + min_cluster_num;   //the cluster num is between[minclusternum,maxclusternum]
	//popul_param[0] = (double)k_rand; //popul_param每一行都有两列，第一列是保存该行对应的簇的簇数，第二列对应的是该行对应的簇的适应度
	//int* rand_arr = rand_generator(k_rand, N);
	//for (int j = 0; j < k_rand; j++) {
	//	for (int k = 0; k < D; k++) {
	//		// popul_index[i][j*D + k] = lk[k] + URAND *(uk[k] - lk[k]);//初始化操作，随机了k_rand个簇，这里是每个簇的每一维度都随机初始化一个取值范围内的值
	//		popul_index[j * D + k] = data[rand_arr[j]][k];
	//	}
	//}
	initial_cluster_center( data,popul_index,  popul_param,  N,  D, max_cluster_num,  min_cluster_num);

	//evaluate the fitness of the initial population
	//for (int j = 0; j < N; j++) {

	//	double min = DBL_MAX;

	//	for (int l = 0; l < ((int)popul_param[0]); l++) {
	//		double distance = getDistance(data[j], &popul_index[l * D], D);
	//		// printf("\n样本%d与簇%d的距离为%lf",j,l,distance);
	//		if (distance < min) {
	//			min = distance;
	//			in_cluster[j] = l;    //record the data index，用来记录当前样本在当前个体的下标为l的簇中。
	//		}
	//	}

	//}
	initial_incluster(data, popul_index, popul_param, in_cluster, N, D);
	dealwith_emptyCluster(in_cluster, N, popul_param[0], data, D, popul_index);
}

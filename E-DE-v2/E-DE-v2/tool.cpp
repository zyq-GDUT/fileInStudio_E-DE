#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include <float.h>
#include"ede.h"

/*****************************************************************
* func_name: distance calculation
* input: the two testing vertors and the dimension length
* descript: for calculating the Euclidean distance 欧几里得距离
*****************************************************************/
double getDistance(double* avector, double* bvector, int n)
{
	int i;
	double sum = 0;
	for (i = 0; i < n; i++)
		sum = sum + pow(*(avector + i) - *(bvector + i), 2); //两个向量的对应维度的平方和再开方即是欧几里得距离

	return sqrt(sum);
}

/*
*数组incluster保存了每个样本对应的簇编号，N是数组长度，也是样本总数
* k是指总的簇数，簇编号是0~K-1
* 本方法返回每个簇对应的样本数
*/
int* getDataNumInCluster(int* incluster, int N, int K) {
	if (!incluster) {
		printf("cluster指针为空!");
		exit(1);
	}
	int* clusterArr;
	malloc1E(clusterArr, K);
	memset(clusterArr, 0, sizeof(int) * K);
	for (int i = 0; i < N; i++) {
		clusterArr[incluster[i]]++;
	}
	return clusterArr;
}

/**
*用来判断是否有空簇
*样本数小于2的簇为空。
*N为数组incluster的长度，即数据集的个数
*数组每个下标都表示对应的样本编号，对应的值为该样本所在簇编号
*k为簇的个数。
*返回值1表示有空簇。
*/
int isExistEmptyCluster(int* incluster,int N,int K) {
	if (!incluster) {
		printf("cluster指针为空!");
		exit(1);
	}
	int flag = 0;
	//int c[100] = { 0 };//用来记录每个簇的样本数
	//for (int i = 0; i < N; i++) {
	//	c[incluster[i]]++;
	//}
	int* c = getDataNumInCluster(incluster, N, K);
	for (int j = 0; j < K; j++) {
		if (c[j] < 2) {
			flag = 1;
			break;
		}
	}
	free(c);
	return flag;


}

/*
返回一个D维的簇中心.
N是样本个数
K是簇数
D是样本维度
*/
double* getDdimCluCenter(double** data, int N,int K ,int D) {
	int num = N / K;//用来生成簇中心的样本个数
	double* center;
	malloc1D(center, D);
	memset(center, 0, sizeof(double) * D);
	for (int i = 0; i < num; i++) {
		int s=rand() % N;//这里随机选择的num个样本可能会出现重复，但是不影响总体结果
		for (int j = 0; j < D; j++) {
			center[j] += data[s][j];
		}
	}
	//求平均
	for (int i = 0; i < D; i++) {
		center[i] = center[i] / num;
	}
	return center;
}

/*
//根据新的簇中心重新聚类,聚类结果保存到数组in_cluster中
*/
void dataInCluster(double** data, int N, int K, double* cluster,int D,int* in_cluster) {
	//int* Incluster;
	//malloc1E(Incluster, N);
	double Min_dist;
	int index=0;
	for (int i = 0; i < N; i++) {
		Min_dist= DBL_MAX;
		for (int j = 0; j < K; j++) {
			double d=getDistance(data[i], cluster+(j * D), D);
			if (d < Min_dist) {
				Min_dist = d;
				index = j;
			}
		}
		//*(Incluster+i) = index;
		in_cluster[i] = index;
	}
	
}

/*
由于每个簇的样本数至少为2，因此数据集的样本数至少为2*最大簇数
若样本数小于2*最大簇数，因此最大簇数就只能是样本数/2
返回最大簇数
*/
int getMaxClusterNum(int N, int maxClusterNum) {

	if (N < 2 * maxClusterNum) {
		maxClusterNum = N / 2;
	}
	return maxClusterNum;
}

/*
找到初始种群中最小的适应度或对应的簇数
flag是1则返回最小适应度，0则返回对应的簇数
*/
double getMinFitness(double* param,int np,int flag) {
	double minFitness = DBL_MAX;
	int clusterNum = 0;
	for (int i = 0; i < np; i++) {
		if (*(param+i*2+1) < minFitness) {
			minFitness = *(param + i * 2 + 1);
			clusterNum = *(param + i * 2 + 0);
		}

	}
	if (flag == 1) {
		return minFitness;
	}
	else {
		return clusterNum;
	}
	
}

/*
找到初始种群中最大的适应度或对应的簇数
flag是1则返回最小适应度，0则返回对应的簇数
*/
double getMaxFitness(double* param, int np, int flag) {
	double maxFitness = DBL_MIN;
	int clusterNum = 0;
	for (int i = 0; i < np; i++) {
		if (*(param + i * 2 + 1) > maxFitness) {
			maxFitness = *(param + i * 2 + 1);
			clusterNum = *(param + i * 2 + 0);
		}

	}
	if (flag == 1) {
		return maxFitness;
	}
	else {
		return clusterNum;
	}

}
/*
该方法是用于杂交后试验变量的簇数小于2或者大于20的情况，生成新的试验向量
该方法会将新生成的试验向量直接保存到next_index中，即没有返回值。
该种方式获得试验向量可能会存在相同的簇中心,但是该方法之后的流程会对该种情况作出处理，因此这里就不会对这种情况进行处理
p_index保存了当前处理个体的簇中心
p_param保存了当前个体的簇数
next_index保存了当前变异个体的簇中心
next_param保存了当前变异个体的簇数
D是指维度
*/
void dealWithClusterNumNotCorrect(double* p_index, double* p_param, double* next_index, double* next_param, int D) {
	double* trailVector;
	double min;
	double distance=0;
	int min_index=0;
	int i, j, k;
	//试验向量的簇数与目标向量的簇数保存一致
	malloc1D(trailVector,p_param[0]*D);
	for ( i = 0; i < p_param[0]; i++) {
		min = DBL_MAX;
		for ( j = 0; j < next_param[0]; j++) {
			distance=getDistance(p_index + i * D, next_index + j * D, D);
			if (distance < min) {
				min = distance;
				min_index = j;
			}
		}
		for (int k = 0; k < D; k++) {
			*(trailVector + i * D + k) = *(next_index + min_index * D + k);
		}
	}
	//将试验向量保存到数组next_index中
	for (int l= 0; l < p_param[0]; l++) {
		for (int u = 0; u < D; u++) {
			*(next_index + l * D + u) = *(trailVector + l * D + u);
		}
	}
	next_param[0] = p_param[0];

}





/*
* 处理空簇,处理后当前个体的每个簇中心对应的样本数都大于等于2
*/

void dealwith_emptyCluster(int* in_cluster, int N, double& k_num, double** data, int D,double* p_index) {
	int count = 0;
	int* c;
	while (isExistEmptyCluster(in_cluster, N, (int)k_num)) {
		
		//printf("当前是第%d代，处理空簇的次数共为%d\n",k,++fixEmptyClusterNum);
		 c = getDataNumInCluster(in_cluster, N, (int)k_num);
		//对空簇进行处理
		for (int l = 0; l < (int)k_num; l++) {
			if (c[l] < 2) {
				double* center = getDdimCluCenter(data, N, (int)k_num, D);
				//把新生成的簇中心替换到空簇的簇中心
				for (int j = 0; j < D; j++) {
					p_index[l * D + j] = center[j];
				}
				free(center);
			}
		}
		free(c);
		//free(in_cluster2[i]);
		//根据新的簇中心重新聚类
		//in_cluster2[i] = dataInCluster(data, N, next_param[i][0],next_index[i],D);
		dataInCluster(data, N, (int)k_num, p_index, D, in_cluster);
		count++;
		//如果循环了5次都不能消除空簇的情况则该个体重新初始化
		if (count > 5) {
			k_num--;
			if (k_num < 2) k_num = 2;
			int* rand_arr = rand_generator((int)k_num, N);
			for (int j = 0; j < (int)k_num; j++) {
				for (int k = 0; k < D; k++) {
					p_index[j * D + k] = data[rand_arr[j]][k];
				}
			}
			dataInCluster(data, N, (int)k_num, p_index, D, in_cluster);
			count = 0;
		}
	}
}



/*
* 释放double类型的二维数组的内存空间
*/
void freeTwoDimArr_double(double** &tdd,int x) {
	/*for (int i = 0; i < x; i++) {
		free(tdd[i]);
	}
	free(tdd);*/
	free(tdd[0]);
	free(tdd);
}

/*
* 释放double类型的二维数组的内存空间
*/
void freeTwoDimArr_int(int** &tdi, int x) {
	/*for (int i = 0; i < x; i++) {
		free(tdi[i]);
	}
	free(tdi);*/
	free(tdi[0]);
	free(tdi);
}









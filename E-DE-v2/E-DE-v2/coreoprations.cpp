#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctime>
#include <malloc.h>
#include <float.h>
#include<string.h>
#include<string>
#include"ede.h"

/*
* 通过评价指标获取最优簇数和适应度
*/
void getBestValByIndex(double& bestFitness,double &bestClusterNum,double* popul_param,int NP,int usedIndex){
	//使用DB指标
	if (usedIndex == DB_INDEX) {
		bestFitness = getMinFitness(popul_param, NP, GET_FITNESS);
		bestClusterNum = getMinFitness(popul_param, NP, GET_CLUSTER_NUM);
	}
	//使用I指标
	else if (usedIndex == I_INDEX|| usedIndex == SIHOUETTES_INDEX) {
		bestFitness = getMaxFitness(popul_param, NP, GET_FITNESS);
		bestClusterNum = getMaxFitness(popul_param, NP, GET_CLUSTER_NUM);
	}

}

/*
* 初始化用于变异，杂交操作的几个核心变量
* next_param用于保存变异向量的簇数，杂交操作后用于保存试验向量的簇数和适应度
* next_index用于保存变异向量的簇中心，杂交操作后用于保存试验向量的簇中心
* in_cluster2用于保存试验向量的聚类情况，即每个下标表示样本编号值表示对应的簇编号
* in_area用于保存杂交操作时从变异向量选出的簇中心
* out_area用于保存从目标个体选出的向量和in_area中的向量，然后该数组会把值赋给保存试验向量的数组
*/
void initialCoreVriable(double** next_param,int** in_cluster2,double** out_area,double** in_area,double** next_index,int NP,int max_cluster_num,int D,int N) {
	for (int i = 0; i < NP; i++) {
		for (int j = 0; j < 2; j++) {
			next_param[i][j] = 0;
		}
		for (int j = 0; j < N; j++) {
			in_cluster2[i][j] = 0;
		}
		for (int j = 0; j < max_cluster_num * D; j++) {
			out_area[i][j] = 0;
			in_area[i][j] = 0;
			next_index[i][j] = 0; //用来保存变异个体
		}
	}
}

/*
* 初始化数组dist和center
* dist用于保存杂交子空间每维度的距离
* center用于保存杂交子空间的中心坐标
* D表示数据的维度
*/
void initial_dist_center(double* dist, double* center,int D) {
	for (int i = 0; i < D; i++) {
		dist[i] = 0; //距离数组，用来记录杂交操作时，杂交子空间中，每一维的取值范围
		center[i] = 0; //用来记录子空间中的簇中心
	}
}

/*
* 生成三个随机数，并且四个数互不相同
* NP是随机数的范围，即0-(NP-1)
*/
void getRandomNum(int& r1,int& r2,int& r3,int i ,int NP) {
	//Mutation operator
			//step 1 of the mutation: random select 3 individuals
	do {
		r1 = (int)(NP * URAND); //在变异操作中，变异个体是基于第r1个个体产生的，r2,r3个体提供差分的信息，然后在r1上做出相应改变，即得到了当前目标个体的变异个体
	} while (r1 == i);
	do {
		r2 = (int)(NP * URAND);
	} while (r2 == i || r2 == r1);
	do {
		r3 = (int)(NP * URAND);
	} while (r3 == i || r3 == r1 || r3 == r2);
}

/*
* 越界处理
* 数值越界则取边界值
*/
int cross_border_process(int num,int ub,int lb) {
	int new_num = num;
	if (new_num > ub) {
		new_num = ub;
	}
	if (new_num < lb) {
		new_num = lb;
	}
	return new_num;
}

/*
* 执行变异操作
* mutate_u 变异个体的簇数
* mutate_d1,mutate_d2,mutate_d3从目标向量中随机选取的三个个体的簇数
* r1,r2,r3三个随机选取个体的个体编号
* D 数据的维度
* N 数据总数
* next_index 用于保存变异个体的簇中心
* popul_index 用于保存目标个体的簇中心
*/
//void mutate_process(double* popul_index,double* next_index,int mutate_u, int mutate_d1, int mutate_d2, int mutate_d3, int r1, int r2, int r3, int D, int N) {
//	if () {
//
//	}
//	else if () {
//
//	}
//	else {
//
//	}
//}


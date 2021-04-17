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
* ͨ������ָ���ȡ���Ŵ�������Ӧ��
*/
void getBestValByIndex(double& bestFitness,double &bestClusterNum,double* popul_param,int NP,int usedIndex){
	//ʹ��DBָ��
	if (usedIndex == DB_INDEX) {
		bestFitness = getMinFitness(popul_param, NP, GET_FITNESS);
		bestClusterNum = getMinFitness(popul_param, NP, GET_CLUSTER_NUM);
	}
	//ʹ��Iָ��
	else if (usedIndex == I_INDEX|| usedIndex == SIHOUETTES_INDEX) {
		bestFitness = getMaxFitness(popul_param, NP, GET_FITNESS);
		bestClusterNum = getMaxFitness(popul_param, NP, GET_CLUSTER_NUM);
	}

}

/*
* ��ʼ�����ڱ��죬�ӽ������ļ������ı���
* next_param���ڱ�����������Ĵ������ӽ����������ڱ������������Ĵ�������Ӧ��
* next_index���ڱ�����������Ĵ����ģ��ӽ����������ڱ������������Ĵ�����
* in_cluster2���ڱ������������ľ����������ÿ���±��ʾ�������ֵ��ʾ��Ӧ�Ĵر��
* in_area���ڱ����ӽ�����ʱ�ӱ�������ѡ���Ĵ�����
* out_area���ڱ����Ŀ�����ѡ����������in_area�е�������Ȼ���������ֵ����������������������
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
			next_index[i][j] = 0; //��������������
		}
	}
}

/*
* ��ʼ������dist��center
* dist���ڱ����ӽ��ӿռ�ÿά�ȵľ���
* center���ڱ����ӽ��ӿռ����������
* D��ʾ���ݵ�ά��
*/
void initial_dist_center(double* dist, double* center,int D) {
	for (int i = 0; i < D; i++) {
		dist[i] = 0; //�������飬������¼�ӽ�����ʱ���ӽ��ӿռ��У�ÿһά��ȡֵ��Χ
		center[i] = 0; //������¼�ӿռ��еĴ�����
	}
}

/*
* ��������������������ĸ���������ͬ
* NP��������ķ�Χ����0-(NP-1)
*/
void getRandomNum(int& r1,int& r2,int& r3,int i ,int NP) {
	//Mutation operator
			//step 1 of the mutation: random select 3 individuals
	do {
		r1 = (int)(NP * URAND); //�ڱ�������У���������ǻ��ڵ�r1����������ģ�r2,r3�����ṩ��ֵ���Ϣ��Ȼ����r1��������Ӧ�ı䣬���õ��˵�ǰĿ�����ı������
	} while (r1 == i);
	do {
		r2 = (int)(NP * URAND);
	} while (r2 == i || r2 == r1);
	do {
		r3 = (int)(NP * URAND);
	} while (r3 == i || r3 == r1 || r3 == r2);
}

/*
* Խ�紦��
* ��ֵԽ����ȡ�߽�ֵ
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
* ִ�б������
* mutate_u �������Ĵ���
* mutate_d1,mutate_d2,mutate_d3��Ŀ�����������ѡȡ����������Ĵ���
* r1,r2,r3�������ѡȡ����ĸ�����
* D ���ݵ�ά��
* N ��������
* next_index ���ڱ���������Ĵ�����
* popul_index ���ڱ���Ŀ�����Ĵ�����
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


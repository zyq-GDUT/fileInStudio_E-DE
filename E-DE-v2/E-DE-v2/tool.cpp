#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include <float.h>
#include"ede.h"

/*****************************************************************
* func_name: distance calculation
* input: the two testing vertors and the dimension length
* descript: for calculating the Euclidean distance ŷ����þ���
*****************************************************************/
double getDistance(double* avector, double* bvector, int n)
{
	int i;
	double sum = 0;
	for (i = 0; i < n; i++)
		sum = sum + pow(*(avector + i) - *(bvector + i), 2); //���������Ķ�Ӧά�ȵ�ƽ�����ٿ�������ŷ����þ���

	return sqrt(sum);
}

/*
*����incluster������ÿ��������Ӧ�Ĵر�ţ�N�����鳤�ȣ�Ҳ����������
* k��ָ�ܵĴ������ر����0~K-1
* ����������ÿ���ض�Ӧ��������
*/
int* getDataNumInCluster(int* incluster, int N, int K) {
	if (!incluster) {
		printf("clusterָ��Ϊ��!");
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
*�����ж��Ƿ��пմ�
*������С��2�Ĵ�Ϊ�ա�
*NΪ����incluster�ĳ��ȣ������ݼ��ĸ���
*����ÿ���±궼��ʾ��Ӧ��������ţ���Ӧ��ֵΪ���������ڴر��
*kΪ�صĸ�����
*����ֵ1��ʾ�пմء�
*/
int isExistEmptyCluster(int* incluster,int N,int K) {
	if (!incluster) {
		printf("clusterָ��Ϊ��!");
		exit(1);
	}
	int flag = 0;
	//int c[100] = { 0 };//������¼ÿ���ص�������
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
����һ��Dά�Ĵ�����.
N����������
K�Ǵ���
D������ά��
*/
double* getDdimCluCenter(double** data, int N,int K ,int D) {
	int num = N / K;//�������ɴ����ĵ���������
	double* center;
	malloc1D(center, D);
	memset(center, 0, sizeof(double) * D);
	for (int i = 0; i < num; i++) {
		int s=rand() % N;//�������ѡ���num���������ܻ�����ظ������ǲ�Ӱ��������
		for (int j = 0; j < D; j++) {
			center[j] += data[s][j];
		}
	}
	//��ƽ��
	for (int i = 0; i < D; i++) {
		center[i] = center[i] / num;
	}
	return center;
}

/*
//�����µĴ��������¾���,���������浽����in_cluster��
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
����ÿ���ص�����������Ϊ2��������ݼ�������������Ϊ2*������
��������С��2*�������������������ֻ����������/2
����������
*/
int getMaxClusterNum(int N, int maxClusterNum) {

	if (N < 2 * maxClusterNum) {
		maxClusterNum = N / 2;
	}
	return maxClusterNum;
}

/*
�ҵ���ʼ��Ⱥ����С����Ӧ�Ȼ��Ӧ�Ĵ���
flag��1�򷵻���С��Ӧ�ȣ�0�򷵻ض�Ӧ�Ĵ���
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
�ҵ���ʼ��Ⱥ��������Ӧ�Ȼ��Ӧ�Ĵ���
flag��1�򷵻���С��Ӧ�ȣ�0�򷵻ض�Ӧ�Ĵ���
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
�÷����������ӽ�����������Ĵ���С��2���ߴ���20������������µ���������
�÷����Ὣ�����ɵ���������ֱ�ӱ��浽next_index�У���û�з���ֵ��
���ַ�ʽ��������������ܻ������ͬ�Ĵ�����,���Ǹ÷���֮������̻�Ը���������������������Ͳ��������������д���
p_index�����˵�ǰ�������Ĵ�����
p_param�����˵�ǰ����Ĵ���
next_index�����˵�ǰ�������Ĵ�����
next_param�����˵�ǰ�������Ĵ���
D��ָά��
*/
void dealWithClusterNumNotCorrect(double* p_index, double* p_param, double* next_index, double* next_param, int D) {
	double* trailVector;
	double min;
	double distance=0;
	int min_index=0;
	int i, j, k;
	//���������Ĵ�����Ŀ�������Ĵ�������һ��
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
	//�������������浽����next_index��
	for (int l= 0; l < p_param[0]; l++) {
		for (int u = 0; u < D; u++) {
			*(next_index + l * D + u) = *(trailVector + l * D + u);
		}
	}
	next_param[0] = p_param[0];

}





/*
* ����մ�,�����ǰ�����ÿ�������Ķ�Ӧ�������������ڵ���2
*/

void dealwith_emptyCluster(int* in_cluster, int N, double& k_num, double** data, int D,double* p_index) {
	int count = 0;
	int* c;
	while (isExistEmptyCluster(in_cluster, N, (int)k_num)) {
		
		//printf("��ǰ�ǵ�%d��������մصĴ�����Ϊ%d\n",k,++fixEmptyClusterNum);
		 c = getDataNumInCluster(in_cluster, N, (int)k_num);
		//�Կմؽ��д���
		for (int l = 0; l < (int)k_num; l++) {
			if (c[l] < 2) {
				double* center = getDdimCluCenter(data, N, (int)k_num, D);
				//�������ɵĴ������滻���մصĴ�����
				for (int j = 0; j < D; j++) {
					p_index[l * D + j] = center[j];
				}
				free(center);
			}
		}
		free(c);
		//free(in_cluster2[i]);
		//�����µĴ��������¾���
		//in_cluster2[i] = dataInCluster(data, N, next_param[i][0],next_index[i],D);
		dataInCluster(data, N, (int)k_num, p_index, D, in_cluster);
		count++;
		//���ѭ����5�ζ����������մص������ø������³�ʼ��
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
* �ͷ�double���͵Ķ�ά������ڴ�ռ�
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
* �ͷ�double���͵Ķ�ά������ڴ�ռ�
*/
void freeTwoDimArr_int(int** &tdi, int x) {
	/*for (int i = 0; i < x; i++) {
		free(tdi[i]);
	}
	free(tdi);*/
	free(tdi[0]);
	free(tdi);
}









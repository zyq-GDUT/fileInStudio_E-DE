#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctime>
#include <malloc.h>
#include <float.h>
#include<string.h>
#include"ede.h"
/*
输出种群信息
p存放了每个个体的信息，即每种聚类方案的簇中心坐标
p2保存了每个个体的簇数
D是数据的维度
np是种群规模，即个体数
maxClusterNum是最大的簇数
*/
void outputPopulationInfo(double* p,double* p2,int D,int np,int maxClusterNum,int** incluster,int N) {
	
	for (int i = 0; i < np; i++) {
		
		int k_num = p2[i * 2];
		int* c = getDataNumInCluster(incluster[i],N,k_num );
		printf("\n个体%d,%d个簇，适应度为%g：\n",i,k_num,p2[i*2+1]);
		for (int j = 0; j < k_num; j++) {
			printf("  %d号簇,%d个样本：",j,c[j]);
			printf("[");
			for (int l = 0; l < D; l++) {
				printf("%g", *(p + i * maxClusterNum * D + j * D + l));
				if (l != D - 1) {
					printf(",");
				}
				else{
					printf("],");
				}
				
			}
			if ((j+1) % 3 == 0) {
				printf("\n");
			}
	

		}
		
		free(c);
	}
	
}

double* getclusterCenter(double cluster[][7],int num,int d) {
	double* meanV;//用来保存每个簇的中心坐标（均值向量）
	//int* numInclu;//用来保存每个簇对应的样本数
	malloc1D(meanV, d);
	//malloc1E(numInclu, k_num);
	memset(meanV, 0, sizeof(*meanV) *d);
	//memset(numInclu, 0, sizeof(*numInclu) * k_num);
	for (int i = 0; i <num ; i++) {
		
		for (int j = 0; j < d; j++) {
			meanV[j] += cluster[i][j];
		}
	}

	//求均值向量
	
		for (int i = 0; i < d; i++) {
			meanV[i] = meanV[i] / num;
		}
	
	
	return meanV;


}


void testfitness() {
	const char* filename = "ecoli.txt";
	double** data;
	int n, d;
	double clu_cp[142][7] = { 0 };
	double clu_im[77][7] = { 0 };
	double clu_imS[2][7] = { 0 };
	double clu_imL[2][7] = { 0 };
	double clu_imu[35][7] = { 0 };
	double clu_om[20][7] = { 0 };
	double clu_oml[5][7] = { 0 };
	double clu_pp[52][7] = { 0 };

	double* clu_cp_center,*clu_im_center, *clu_imS_center, *clu_imL_center, *clu_imu_center, *clu_om_center, *clu_oml_center, *clu_pp_center;
	int* in_cluster;
	data=loadData2(filename,n,d);
	malloc1E(in_cluster, n);
	for (int i = 0; i < n; i++) {
		printf("第%d行:\t", i + 1);
		if (i <142) {
			for (int j = 0; j < d; j++) {
				clu_cp[i][j] = data[i][j];
				in_cluster[i] = 0;
				printf("%g\t", clu_cp[i][j]);
			}
			
		}
		else if (i >= 142 && i < 219) {
			for (int j = 0; j < d; j++) {
				clu_im[i - 142][j] = data[i][j];
				in_cluster[i] = 1;
				printf("%g\t", clu_im[i - 142][j]);

			}
		}
		else if (i >= 219 && i < 221) {
			for (int j = 0; j < d; j++) {
				clu_imS[i - 219][j] = data[i][j];
				in_cluster[i] = 2;
				printf("%g\t", clu_imS[i - 219][j]);
			}
		}
		else if (i >= 221 && i < 223) {
			for (int j = 0; j < d; j++) {
				clu_imL[i - 221][j] = data[i][j];
				in_cluster[i] = 3;
				printf("%g\t", clu_imL[i - 221][j]);
			}
		}
		else if (i >= 223 && i < 258) {
			for (int j = 0; j < d; j++) {
				clu_imu[i - 223][j] = data[i][j];
				in_cluster[i] = 4;
				printf("%g\t", clu_imu[i - 223][j]);
			}
		}
		else if (i >= 258 && i < 278) {
			for (int j = 0; j < d; j++) {
				clu_om[i - 258][j] = data[i][j];
				in_cluster[i] = 5;
				printf("%g\t", clu_om[i - 258][j]);
			}
		}
		else if (i >= 278 && i < 283) {
			for (int j = 0; j < d; j++) {
				clu_oml[i - 278][j] = data[i][j];
				in_cluster[i] = 6;
				printf("%g\t", clu_oml[i - 278][j]);
			}
		}
		else if (i >= 283 ) {
			for (int j = 0; j < d; j++) {
				clu_pp[i - 283][j] = data[i][j];
				in_cluster[i] = 7;
				printf("%g\t", clu_pp[i - 283][j]);
			}
		}
		printf("\n");
		
		
	}
	clu_cp_center = getclusterCenter(clu_cp, 142, 7);
	clu_im_center = getclusterCenter(clu_im, 77, 7);
	clu_imS_center = getclusterCenter(clu_imS, 2, 7);
	clu_imL_center = getclusterCenter(clu_imL, 2, 7);
	clu_imu_center = getclusterCenter(clu_imu, 35, 7);
	clu_om_center = getclusterCenter(clu_om, 20, 7);
	clu_oml_center = getclusterCenter(clu_oml, 5, 7);
	clu_pp_center = getclusterCenter(clu_pp, 52, 7);
	double p_index[7 * 8] = { 0 };
	for (int i = 0; i < 7; i++) {
		p_index[i] = clu_cp_center[i];
		p_index[7+i] = clu_im_center[i];
		p_index[14 + i] = clu_imS_center[i];
		p_index[21 + i] = clu_imL_center[i];
		p_index[28 + i] = clu_imu_center[i];
		p_index[35 + i] = clu_om_center[i];
		p_index[42 + i] = clu_oml_center[i];
		p_index[49 + i] = clu_pp_center[i];
	
	}
	double f=DB_index(8, *data, p_index, in_cluster, n, d);
	printf("DB值=%g\n", f);
	free(clu_cp_center);
	free(clu_im_center);
	free(clu_imS_center);
	free(clu_imL_center);
	free(clu_imu_center);
	free(clu_om_center);
	free(clu_oml_center);
	free(clu_pp_center);
	free(in_cluster);
}

/*
* gap是指每隔gap代输出一次
*/
void outputPopulationInfo2(double* p, double* p2, int D, int np, int maxClusterNum, int** incluster, int N,int generation,double bestfitness,double bestclusternum,int gap) {
	printf("\n---------第%d代-最优适应度为%g-对应的簇数为%d------------------\n", generation + 1, bestfitness, (int)bestclusternum);
	if (generation % gap == 0) {
		outputPopulationInfo(p, p2, D, np, maxClusterNum, incluster, N);

	}
}



//int main() {
//	testfitness();	
//	return 0;
//}
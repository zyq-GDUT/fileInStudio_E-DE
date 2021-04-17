#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include <float.h>
#include "ede.h"

/**
 * 计算当前个体每个簇的均值向量和每个簇对应的样本数
 * k_num是指总簇数
 * **/
double* count_meanV(double* data, int* in_clu, int k_num, int N, int D) {
	double* meanV;//用来保存每个簇的中心坐标（均值向量）
	int* numInclu;//用来保存每个簇对应的样本数
	malloc1D(meanV, k_num * D);
	malloc1E(numInclu, k_num);
	memset(meanV, 0, sizeof(*meanV) * k_num * D);
	memset(numInclu, 0, sizeof(*numInclu) * k_num);
	for (int i = 0; i < N; i++) {
		int cn = in_clu[i];
		numInclu[cn]++;
		for (int j = 0; j < D; j++) {
			*(meanV + cn * D + j) += *(data + i * D + j);
		}
	}

	//求均值向量
	for (int k = 0; k < k_num; k++) {
		for (int i = 0; i < D; i++) {
			*(meanV + k * D + i) = *(meanV + k * D + i) / numInclu[k];
		}
	}
	free(numInclu);
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
double I_index(int k_num, double* data, double* p_index, int* in_clu, int N, int D) {

	//p_index = count_meanV(data, in_clu, k_num, N, D);
	int i, j, k;
	double E = 0;
	//计算每个样本与簇中心的距离之和，即计算E值
	for (i = 0; i < N; i++) {
		double d = 0;
		int cnum = in_clu[i];
		for (j = 0; j < D; j++) {
			//计算样本与其簇中心的欧氏距离
			d += pow(*(data + i * D + j) - *(p_index + cnum * D + j),2);
		}
		E += sqrt(d);
	}
	// printf("E值=%lf\n",E);
	//找出所有簇中距离最远的两个簇的距离。即计算Dmax值
	double Dmax = 0.0;
	double distance;
	for (i = 0; i < k_num; i++) {
		for (j = i + 1; j < k_num; j++) {
			distance = 0.0;
			for (k = 0; k < D; k++) {
				distance+=pow(*(p_index + i * D + k) - *(p_index + j * D + k),2);
			}
			distance = sqrt(distance);
			if (distance > Dmax) {
				Dmax = distance;
			}
		}
	}

	// printf("k_num=%d\tN=%d\tE=%lf\tD=%lf\n",k_num,N,E,Dmax);
	//计算I值,要注意若除数是整数则除法表达式就相当于是整除。
	double i_index = pow((1.0 / (double)k_num) * ((double)N / E) * Dmax, 2);
	//printf("i值=%g\n", i_index);
	//printf("\n------\n");
	//free(p_index);
	return i_index;

}

//DB指标函数
double DB_index(int k_num, double* data, double* p_index, int* in_clu, int N, int D) {
	int c;
	double dist;
	double* e;//用来保存每个簇与簇中心的平均误差
	int* numInclu;//用来记录每个簇的样本数
	//p_index = count_meanV(data, in_clu, k_num, N, D);//求每个簇的中心坐标
	malloc1D(e, k_num);
	malloc1E(numInclu, k_num);
	memset(e, 0, sizeof(double) * k_num);
	memset(numInclu, 0, sizeof(int) * k_num);

	for (int i = 0; i < N; i++) {
		c = in_clu[i];
		numInclu[c] = numInclu[c] + 1;
		dist = 0;
		//计算样本与簇中心的欧氏距离的平方
		dist=getDistance(data + i * D, p_index + c * D, D);
		dist = dist * dist;
		e[c] += dist;

	}
	//求均值
	for (int i = 0; i < k_num; i++) {
		//printf("\n当前的簇编号是%d，样本个数为%d  e[%d]=%lf\n",i, numInclu[i],i, e[i]);
		//如果某个簇对应的样本数是0，则下面求e值时会出现0/0的情况，求出来的e值是一个nan(ind)值，即不是数值。因此在前面
		e[i] = e[i] / numInclu[i];
		
	}
	

	double** ed;
	double dist_i_j;
	malloc2E(ed, k_num, k_num);
	for (int i = 0; i < k_num; i++) {
		for (int j = i+1; j < k_num; j++) {
			dist_i_j = 0.0;
			//计算两个簇中心的欧氏距离的平方。即论文公式（1）中的Dij
			for (int k = 0; k < D; k++) {
				dist_i_j += pow(*(p_index + i * D + k) - *(p_index + j * D + k),2);
			}
			ed[i][j] = (e[i] + e[j]) / dist_i_j;
			ed[j][i] = ed[i][j];
		}
	}
	double maxe;
	//求ed矩阵中每行的最大值
	for (int i = 0; i < k_num; i++) {
		maxe = DBL_MIN;
		for (int j = 0; j < k_num; j++) {
			
			if (i == j) {
				continue;
			}
			if (maxe < ed[i][j]) {
				maxe = ed[i][j];
				ed[i][i] =(double) j;
			}
			//printf("ed[%d][%d]=%g\t", i, j, ed[i][j]);

		}
		//printf("\n");
	}

	//求DB值
	double DB = 0;
	for (int i = 0; i < k_num; i++) {
		int num = (int)ed[i][i];
		DB += ed[i][num];
	}


	DB = DB / k_num;


	//printf("DB值=%g\n", DB);
	//free(p_index);
	//free(ed);
	freeTwoDimArr_double(ed,k_num);
	free(e);
	free(numInclu);
	return DB;

}


/*
* 求簇内距离，即每个簇的样本与中心的平均距离之和
* data保存了所有样本的信息
* p_index保存了本次测试中最优个体
* p_param保存了最优个体的簇数
* in_cluster保存了样本对应的簇编号
* N 样本总数
* D 数据的维度
*/
double getIntraDist(double** data, double* p_index, double* p_param, int* in_cluster, int N, int D) {
	double Intra_dist = 0;
	int** clusteringInfo;//每行表示一个簇，每行的第一个元素保存的是该簇的样本个数，从下标1开始保存的是样本的编号
	double* IntraDistPreCluster;//用来保存每个簇的簇内平均距离
	malloc2D(clusteringInfo, p_param[0], N);
	malloc1D(IntraDistPreCluster, p_param[0]);
	for (int i = 0; i < p_param[0]; i++) {
		memset(clusteringInfo[i],0, sizeof(int) * N);

	}
	memset(IntraDistPreCluster, 0, sizeof(double) * p_param[0]);
	for (int i = 0; i < N; i++) {
		int c = in_cluster[i];
		clusteringInfo[c][0]++;
		clusteringInfo[c][clusteringInfo[c][0]]=i;
	}
	for (int i = 0; i < p_param[0]; i++) {
		int num1 = clusteringInfo[i][0];
		for (int j = 1; j <=num1; j++) {
			for (int k = j + 1; k <= num1; k++) {
				int n1 = clusteringInfo[i][j];
				int n2 = clusteringInfo[i][k];
				//计算两个样本的距离
				double dist=getDistance(data[n1], data[n2], D);
				dist = dist * dist;//计算距离的平方
				IntraDistPreCluster[i] += dist;
			}
		}
		double num2 = (double)num1;
		IntraDistPreCluster[i] = 2.0 * IntraDistPreCluster[i] / (num2 * (num2 - 1));
	}
	for (int i = 0; i < p_param[0]; i++) {
		Intra_dist += IntraDistPreCluster[i];
	}
	freeTwoDimArr_int(clusteringInfo, p_param[0]);
	free(IntraDistPreCluster);
	return Intra_dist;

}

/*
* 求簇间距离,即任意两个簇中心距离的平均值
*/
double getInterDist(double* p_index, double* p_param,int D) {
	double Inter_dist = 0;
	for (int i = 0; i < p_param[0]; i++) {
		for (int j = i + 1; j < p_param[0]; j++) {
			double distance = getDistance(p_index + i * D, p_index + j * D, D);
			distance = distance * distance;
			Inter_dist += distance;
		}
	}
	double clusterNum = (double)p_param[0];
	Inter_dist = 2.0 * Inter_dist / (clusterNum * (clusterNum - 1));
	return Inter_dist;
}

/*
* 返回每个簇对应了所有样本编号的二维数组
* incluster保存了每个样本对应的簇编号，数组的下标为样本编号，值为簇编号
* N是样本的总数
* k_num是簇数
*/
int** getClusterSambles(int* incluster, int N, int k_num) {
	int** clusterSambles;//数组每行的第一个单元保存的是当前簇的样本数
	int c;
	malloc2D(clusterSambles, k_num, N);
	for (int i = 0; i < k_num; i++) {
		clusterSambles[i][0] = 0;
	}
	for (int i = 0; i < N; i++) {
		 c = incluster[i];
		 clusterSambles[c][++clusterSambles[c][0]] = i;
		 
	}
	return clusterSambles;
}
/*
* 求组合Cnm，n是下标，m是上标
* 若n<m则返回0
*/
double Combination_n_m(int n, int m) {
	if (n < m) {
		return 0.0;
	}
	double s_up = 1.0;
	double s_down = 1.0;
	for (int i = 0; i < m; i++) {
		s_up = s_up * (n - i);
		s_down = s_down * (m - i);
	}
	return (s_up / s_down);
}
/*
* 求ARI值
* originClusterInfo记录了数据集原来的聚类情况
* in_cluster记录了E-DE算法的聚类情况
* k_num是算法聚类后的簇数
* N是样本的总数
*/
double getARI(int* originClusterInfo, int* in_cluster, int k_num, int N) {
	int origin_k_num = originClusterInfo[N];//获取本身的聚类数
	int** contingency_table;
	int** originClusterSambles;
	double ARI = 0;
	malloc2D(contingency_table, origin_k_num+1, k_num+1);
	for (int i = 0; i < origin_k_num + 1; i++) {
		memset(contingency_table[i], 0, sizeof(int) * (k_num + 1));
	}
	originClusterSambles = getClusterSambles(originClusterInfo, N, origin_k_num);
	for (int i = 0; i < origin_k_num; i++) {
		for (int j = 1; j <= originClusterSambles[i][0]; j++) {
			int samble = originClusterSambles[i][j];
			int cn = in_cluster[samble];
			contingency_table[i][cn]++;
			contingency_table[i][k_num]++;
			contingency_table[origin_k_num][cn]++;
			contingency_table[origin_k_num][k_num]++;
		}
	}

	int s_a = 0;//协和矩阵中个数大于2的组合之和
	for (int i = 0; i < origin_k_num; i++) {
		for (int j = 0; j < k_num; j++) {
			s_a += Combination_n_m(contingency_table[i][j], 2);
		}
	}
	double s_lastColumn = 0;//协和矩阵最右边那列的组合之和
	double s_lastRank = 0;//协和矩阵最下边一行的组合之和
	for (int i = 0; i < origin_k_num; i++) {
		s_lastColumn += Combination_n_m(contingency_table[i][k_num],2);
	}
	for (int j = 0; j < k_num; j++) {
		s_lastRank += Combination_n_m(contingency_table[origin_k_num][j], 2);
	}
	double s_n;
	s_n = Combination_n_m(contingency_table[origin_k_num][k_num], 2);

	double ARI_up = s_a - (s_lastRank * s_lastColumn) / s_n;
	double ARI_down = (s_lastRank + s_lastColumn) / 2 -  (s_lastRank * s_lastColumn) / s_n;
	freeTwoDimArr_int(contingency_table,origin_k_num);
	freeTwoDimArr_int(originClusterSambles, origin_k_num);
	ARI = ARI_up / ARI_down;
	return ARI;
}







//int  main() {
//	int c;
//	double dist;
//	int k_num = 5;
//	double** e;//用来保存每个簇与簇中心的平均误差与对应簇的样本个数
//	malloc2E(e, k_num, 2);
//	for (int i = 0; i < k_num; i++) {
//		memset(e[i], 0, sizeof(double) * 2);
//	}
//
//	for (int i = 0; i < k_num; i++) {
//		for (int j = 0; j < 2; j++) {
//			printf("%lf\t", e[i][j]);
//		}
//		printf("\n");
//	}
//
//	system("pause");
//	return 0;
//
//}
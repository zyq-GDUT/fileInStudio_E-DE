#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include <float.h>
#include "ede.h"

/**
 * ���㵱ǰ����ÿ���صľ�ֵ������ÿ���ض�Ӧ��������
 * k_num��ָ�ܴ���
 * **/
double* count_meanV(double* data, int* in_clu, int k_num, int N, int D) {
	double* meanV;//��������ÿ���ص��������꣨��ֵ������
	int* numInclu;//��������ÿ���ض�Ӧ��������
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

	//���ֵ����
	for (int k = 0; k < k_num; k++) {
		for (int i = 0; i < D; i++) {
			*(meanV + k * D + i) = *(meanV + k * D + i) / numInclu[k];
		}
	}
	free(numInclu);
	return meanV;
}

/**iָ�꺯��
 * k_num��ָ��ǰ����Ĵ���
 * data����������������
 * p_index�����˵�ǰ��������д���������
 * in_clu��¼��ÿ��������Ӧ�Ĵ��±�
 * N��ָ������
 * D��ָ������ά��
 * **/
double I_index(int k_num, double* data, double* p_index, int* in_clu, int N, int D) {

	//p_index = count_meanV(data, in_clu, k_num, N, D);
	int i, j, k;
	double E = 0;
	//����ÿ������������ĵľ���֮�ͣ�������Eֵ
	for (i = 0; i < N; i++) {
		double d = 0;
		int cnum = in_clu[i];
		for (j = 0; j < D; j++) {
			//����������������ĵ�ŷ�Ͼ���
			d += pow(*(data + i * D + j) - *(p_index + cnum * D + j),2);
		}
		E += sqrt(d);
	}
	// printf("Eֵ=%lf\n",E);
	//�ҳ����д��о�����Զ�������صľ��롣������Dmaxֵ
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
	//����Iֵ,Ҫע����������������������ʽ���൱����������
	double i_index = pow((1.0 / (double)k_num) * ((double)N / E) * Dmax, 2);
	//printf("iֵ=%g\n", i_index);
	//printf("\n------\n");
	//free(p_index);
	return i_index;

}

//DBָ�꺯��
double DB_index(int k_num, double* data, double* p_index, int* in_clu, int N, int D) {
	int c;
	double dist;
	double* e;//��������ÿ����������ĵ�ƽ�����
	int* numInclu;//������¼ÿ���ص�������
	//p_index = count_meanV(data, in_clu, k_num, N, D);//��ÿ���ص���������
	malloc1D(e, k_num);
	malloc1E(numInclu, k_num);
	memset(e, 0, sizeof(double) * k_num);
	memset(numInclu, 0, sizeof(int) * k_num);

	for (int i = 0; i < N; i++) {
		c = in_clu[i];
		numInclu[c] = numInclu[c] + 1;
		dist = 0;
		//��������������ĵ�ŷ�Ͼ����ƽ��
		dist=getDistance(data + i * D, p_index + c * D, D);
		dist = dist * dist;
		e[c] += dist;

	}
	//���ֵ
	for (int i = 0; i < k_num; i++) {
		//printf("\n��ǰ�Ĵر����%d����������Ϊ%d  e[%d]=%lf\n",i, numInclu[i],i, e[i]);
		//���ĳ���ض�Ӧ����������0����������eֵʱ�����0/0��������������eֵ��һ��nan(ind)ֵ����������ֵ�������ǰ��
		e[i] = e[i] / numInclu[i];
		
	}
	

	double** ed;
	double dist_i_j;
	malloc2E(ed, k_num, k_num);
	for (int i = 0; i < k_num; i++) {
		for (int j = i+1; j < k_num; j++) {
			dist_i_j = 0.0;
			//�������������ĵ�ŷ�Ͼ����ƽ���������Ĺ�ʽ��1���е�Dij
			for (int k = 0; k < D; k++) {
				dist_i_j += pow(*(p_index + i * D + k) - *(p_index + j * D + k),2);
			}
			ed[i][j] = (e[i] + e[j]) / dist_i_j;
			ed[j][i] = ed[i][j];
		}
	}
	double maxe;
	//��ed������ÿ�е����ֵ
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

	//��DBֵ
	double DB = 0;
	for (int i = 0; i < k_num; i++) {
		int num = (int)ed[i][i];
		DB += ed[i][num];
	}


	DB = DB / k_num;


	//printf("DBֵ=%g\n", DB);
	//free(p_index);
	//free(ed);
	freeTwoDimArr_double(ed,k_num);
	free(e);
	free(numInclu);
	return DB;

}


/*
* ����ھ��룬��ÿ���ص����������ĵ�ƽ������֮��
* data������������������Ϣ
* p_index�����˱��β��������Ÿ���
* p_param���������Ÿ���Ĵ���
* in_cluster������������Ӧ�Ĵر��
* N ��������
* D ���ݵ�ά��
*/
double getIntraDist(double** data, double* p_index, double* p_param, int* in_cluster, int N, int D) {
	double Intra_dist = 0;
	int** clusteringInfo;//ÿ�б�ʾһ���أ�ÿ�еĵ�һ��Ԫ�ر�����Ǹôص��������������±�1��ʼ������������ı��
	double* IntraDistPreCluster;//��������ÿ���صĴ���ƽ������
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
				//�������������ľ���
				double dist=getDistance(data[n1], data[n2], D);
				dist = dist * dist;//��������ƽ��
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
* ��ؼ����,���������������ľ����ƽ��ֵ
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
* ����ÿ���ض�Ӧ������������ŵĶ�ά����
* incluster������ÿ��������Ӧ�Ĵر�ţ�������±�Ϊ������ţ�ֵΪ�ر��
* N������������
* k_num�Ǵ���
*/
int** getClusterSambles(int* incluster, int N, int k_num) {
	int** clusterSambles;//����ÿ�еĵ�һ����Ԫ������ǵ�ǰ�ص�������
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
* �����Cnm��n���±꣬m���ϱ�
* ��n<m�򷵻�0
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
* ��ARIֵ
* originClusterInfo��¼�����ݼ�ԭ���ľ������
* in_cluster��¼��E-DE�㷨�ľ������
* k_num���㷨�����Ĵ���
* N������������
*/
double getARI(int* originClusterInfo, int* in_cluster, int k_num, int N) {
	int origin_k_num = originClusterInfo[N];//��ȡ����ľ�����
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

	int s_a = 0;//Э�;����и�������2�����֮��
	for (int i = 0; i < origin_k_num; i++) {
		for (int j = 0; j < k_num; j++) {
			s_a += Combination_n_m(contingency_table[i][j], 2);
		}
	}
	double s_lastColumn = 0;//Э�;������ұ����е����֮��
	double s_lastRank = 0;//Э�;������±�һ�е����֮��
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
//	double** e;//��������ÿ����������ĵ�ƽ��������Ӧ�ص���������
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
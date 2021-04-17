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
	else {
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
*�����������е���Сֵ
*/
int getMin(int a,int b) {
	if (a < b) {
		return a;
	}
	else {
		return b;
	}
}

/*
*��������a���Ƿ����v
* len������a�ĳ���
*�����򷵻�1���������򷵻�0
*/
int isContain_v(int* a, int len, int v) {
	int flag=0;
	for (int i = 0; i < len; i++) {
		if (a[i] == v) {
			flag = 1;
			break;
		}
	}
	return flag;
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
void mutate_process(double* p,double* next_index,int mutate_u, int mutate_d1, int mutate_d2, int mutate_d3, int r1, int r2, int r3, int D, int N) {

	if (mutate_u==mutate_d1) {
		//�������Ĵ���������Ŀ�����r1�Ĵ���������ȵ����
		for (int j = 0; j < mutate_d1 * D; j++) //r1Ŀ���������д�����ֱ�Ӹ��Ƶ����������
			next_index[j] = *(p + r1 * 20 * D + j);
			
	}
	else if (mutate_u>mutate_d1) {
		 for (int j = 0; j < mutate_d1 * D; j++) //�Ȱ�r1����Ĵ����ĸ��Ƶ����������
			next_index[j] = *(p + r1 * 20 * D + j);
		 int rand_num = (mutate_u - mutate_d1);
		 int rand_num3=getMin(mutate_d3, rand_num/2);
		 int* rand3=rand_generator(rand_num3, mutate_d3);
		 int* rand2 = rand_generator(rand_num - rand_num3, mutate_d2);
		 int index = mutate_d1*D;
		 for (int i = 0; i < rand_num3; i++) {
			 int c = rand3[i];
			 for(int j=0;j<D;j++){
				 next_index[index + j] = *(p + r3 * 20 * D + c * D + j);
			 }
			 index = index + D;
		 }
		 for (int i = 0; i < (rand_num - rand_num3); i++) {
			 int c = rand2[i];
			 for (int j = 0; j < D; j++) {
				 next_index[index + j] = *(p + r2 * 20 * D + c * D + j);
			 }
			 index = index + D;
		 }
		 free(rand3);
		 free(rand2);

	}
	else {//mutate_u<mutate_d1�����
		int delnum = mutate_d1 - mutate_u;
		int* rand_delnum = rand_generator(delnum, mutate_d1);
		int index = 0;
		for (int i = 0; i < mutate_d1; i++) {
			if (isContain_v(rand_delnum, delnum, i) ){
				continue;
			}
			for (int j = 0; j < D; j++) {
				next_index[index + j] = *(p + r1 * 20 * D + i * D + j);
			}
			index = index + D;

		}
	}
}
/*
* ������尴��һ������������������̫�ֲ����Ŷ�
* next_index�ǵ�ǰ�������
* k_num�Ǳ������Ĵ���
* ub,lb��ָ�Ͻ���½�
* D��ָ���ݵ�ά��
*/
void add_normal_distribution_disturb(double* next_index,int k_num,int D,double* Xu,double* Xl) {
	for (int j = 0; j <k_num; j++) { //��ǰ��ʼ��������ÿһ�������ģ���ÿһά������Ӹ�˹�Ŷ�
		for (int l = 0; l < D; l++) {
			double gauss_index = next_index[j * D + l] + sampleNormal(0, 0.1) * (Xu[l] - Xl[l]);//sampleNormal�����ɸ�˹�������������0�Ǿ�ֵ��0.1�Ǳ�׼���˹�ֲ��������ʾ���ʣ��������ʾ��ֵ�����ɸú���ͼ����Եó���������������0��0.1���������������
			if (gauss_index > Xu[l]) {
				next_index[j * D + l] = Xu[l]; //Խ���ȡ�߽�ֵ
			}
			else if (gauss_index < Xl[l]) {
				next_index[j * D + l] = Xl[l];
			}
			else {
				next_index[j * D + l] = gauss_index;
			}
		}
	}
}

/*
* �ӽ�����,�ڱ���������ѡ���ӽ����壬�����������ӽ��Ĵ���
* next_index�ǵ�ǰ�ı������
* k_num�Ǳ������Ĵ���
* in_area�Ǳ���ӵ�ǰ��������л�ȡ�Ĵ�����
* D������ά��
* CR�Ǳ������
* �÷������شӱ��������ѡ����ӽ�����
*/
int selectClustersFromMutateInd(double* next_index,int k_num,double* in_area,int D,double CR) {
	int cr_lenth = 0, swap = 0, swap1 = 0; 
	int cr_n;
	do {
		cr_lenth = cr_lenth + 1;
	} while (URAND < CR && cr_lenth < k_num);//���ӽ�ǰ���ȼ�������������������Ĵ�������
	cr_n = (int)(k_num * URAND); //���һ���ӽ��㣬ȡֵ��Χ��0~rand_basic-1  
	//step 2 of the crossover: determine the subspace of crossover
	int rand1;
	   //the number of points that outside and inside the swap area
	if ((cr_n + cr_lenth - 1) < k_num) {//rand_basic�ǵ�ǰ�������Ĵ���������cr_n+cr_length�ǵ�ǰ����������һ���ص�λ��
		for (int j = cr_n; j < cr_n + cr_lenth; j++) {
			for (int l = 0; l < D; l++)
				in_area[(j - cr_n) * D + l] = next_index[j * D + l]; //in_area��ά���鱣�����ÿ�����������ӽ��Ĵ����ģ������Ǵ�0��ʼ��������
		}
	}
	else {//������ӽ���ʼ��cr_n��ʼ�ӽ�һֱ�ӽ�cr_length���ӽ����������ӽ��������Ѿ������˵�ǰ�������Ĵ���  
		for (int j = cr_n; j < k_num; j++) {
			for (int l = 0; l < D; l++)
				in_area[(j - cr_n) * D + l] = next_index[j * D + l];
		}
		for (int j = 0; j < (cr_n + cr_lenth - k_num); j++) {
			for (int l = 0; l < D; l++)
				in_area[(j + k_num - cr_n) * D + l] = next_index[j * D + l];
		}
	}
	return cr_lenth;
}



/*
* ��out_area�������Ϣ���Ƶ�next_index�����еõ���������
* k_num��ָout_area�еĴ���
*/
void get_trailVector(double* out_area, int k_num,double* next_index,double*next_param,int D) {
	for (int j = 0; j < k_num * D; j++)
		next_index[j] = out_area[j]; //��out_area�����еĴظ��Ƶ��������б����������飬���õ��������������顣
	next_param[0] = k_num;
}

/*
* �����ӿռ��ӽ�
* cr_lenth �Ǳ�������ӿռ�Ĵ���
* in_area �Ǵӱ��������ȡ���������ӽ��Ĵ�
* out_area �ӽ���Ĵ��ȱ����ڸ����飬����ٱ��浽next_index����
* p_index �����˵�ǰĿ�����
* p_param �����˵�ǰ����Ĵ�������Ӧ��
* next_index ���ڱ����ӽ���Ĵأ�����Ϊ��������
* next_param ��������ʵ�������Ĵ���
* D ���ݵ�ά��
* 
*/
void subspace_cross(int cr_lenth, double* in_area, double* out_area, double* p_index, double* p_param, double* next_index, double* next_param,int D) {
	double* min_b;
	double* max_b;
	double minb, maxb;
	malloc1D(min_b, D);
	malloc1D(max_b, D);
	if (cr_lenth == 1) {//��Ҫ�ӽ��Ĵ�ֻ��һ���������
		for (int j = 0; j < 20 * D; j++)
			out_area[j] =p_index[j]; //�Ȱѵ�ǰĿ���������дر��浽����out_area�С�
		//�����λ��Ŀ������������������20���أ���ô�ӿռ佻��ʱ��ֻ��Ҫ����Ҫ�ӽ��Ĵ��滻��λĿ����������һ���ء�  
		//if ((int)*(p2 + i * 2) == 20) {
		int rand1 = (int)((int)p_param[0] * URAND);
		for (int j = 0; j < D; j++)
			out_area[rand1 * D + j] = in_area[j]; //���cr_lenth==1,������Ĵ�ֻ��һ������ô����in_area[i]ֻ������һ���أ�Ȼ���������滻��out_area��������һ��λ��
		//for (int j = 0; j < 20 * D; j++)
		//	next_index[j] = out_area[j]; //��out_area�����еĴظ��Ƶ��������б����������飬���õ��������������顣
		//next_param[0] = p_param[0]; //next_param���鱣����ǵ�i���������Ĵ�����������Ϊ�ӽ��Ĵ�ֻ��һ����˵�i��������壨�����������Ĵ������λ��Ŀ�����Ĵ�������һ��
		get_trailVector(out_area, (int)p_param[0], next_index, next_param, D);
	}
	else {
		for (int j = 0; j <D; j++) {
			minb = DBL_MAX;
			maxb = DBL_MIN;
			for (int i = 0; i < cr_lenth; i++) {//����Ҫ�ӽ������дأ�����������ÿά��ƽ��ֵ֮�󣬾Ϳ��Եõ���ֵ����
				if (in_area[i * D + j] < minb) {
					minb = in_area[i * D + j];
				}
				if (in_area[i * D + j]> maxb) {
					maxb = in_area[i * D + j];
				}
			}
			min_b[j] = minb;
			max_b[j] = maxb;
		}

		int swap_if;
		int swap1 = 0, swap;
		for (int j = 0; j < (int)p_param[0]; j++) {
			swap_if = 0;
			for (int l = 0; l < D; l++) {
				//������������˵����ǰ�Ĵز����ӽ��ӿռ䷶Χ��
				if (p_index[j * D + l] < min_b[l]|| p_index[j*D+l]>max_b[l]) {   //the center can be changed
					swap_if = 1;
					break;
				}
			}
			if (swap_if == 1) {//�Ѳ���Ҫ�����ǵĴر������������ⲿ�ִ���Ŀ������в���Ҫ����������滻�Ĵء�����ⲿ�ִ����������б�ѡ�е��ӽ��ؽ�ϳ���������
				for (int l = 0; l < D; l++)
					out_area[swap1 * D + l] = p_index[j*D+l];
				swap1 = swap1 + 1;
			}
		}
		swap = cr_lenth + swap1;//swap1��Ŀ������в���Ҫ�������Ĵ�����cr_lenth��ָ��������ӿռ䣨����in_area���Ĵ��������������ǵĴ���

		//step 3 of the crossover : subarea swap
		if (swap <= MAX_CLUSTER && swap >= MIN_CLUSTER) {
			for (int j = 0; j < cr_lenth * D; j++)
				next_index[j] = in_area[j]; //�Ȱѱ�������ӿռ��еĴأ����뵽����������
			for (int j = cr_lenth * D; j < swap * D; j++)
				next_index[j] = out_area[j - cr_lenth * D];//�ٰ�Ŀ������в���Ҫ�����ǵĴز��뵽���������У��ò�����ɼ�ԭ�������������������ͱ�ɱ������������ˡ�
			next_param[0] = swap;
		}else {
			/*���ӽ���Ĵ�������20��С��2��ʱ�򣬽������´���*/
			dealWithClusterNumNotCorrect(p_index, p_param, next_index, next_param, D);
		}
	}

	free(min_b);
	free(max_b);
}
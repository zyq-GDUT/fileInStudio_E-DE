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
	else {
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
*返回两个数中的最小值
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
*返回数组a中是否包含v
* len是数组a的长度
*包含则返回1，不包含则返回0
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
* 执行变异操作
* mutate_u 变异个体的簇数
* mutate_d1,mutate_d2,mutate_d3从目标向量中随机选取的三个个体的簇数
* r1,r2,r3三个随机选取个体的个体编号
* D 数据的维度
* N 数据总数
* next_index 用于保存变异个体的簇中心
* popul_index 用于保存目标个体的簇中心
*/
void mutate_process(double* p,double* next_index,int mutate_u, int mutate_d1, int mutate_d2, int mutate_d3, int r1, int r2, int r3, int D, int N) {

	if (mutate_u==mutate_d1) {
		//变异个体的簇中心数与目标个体r1的簇中心数相等的情况
		for (int j = 0; j < mutate_d1 * D; j++) //r1目标个体的所有簇中心直接复制到变异个体中
			next_index[j] = *(p + r1 * 20 * D + j);
			
	}
	else if (mutate_u>mutate_d1) {
		 for (int j = 0; j < mutate_d1 * D; j++) //先把r1个体的簇中心复制到变异个体中
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
	else {//mutate_u<mutate_d1的情况
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
* 变异个体按照一定概率来产生满足正太分布的扰动
* next_index是当前变异个体
* k_num是变异个体的簇数
* ub,lb是指上界和下界
* D是指数据的维度
*/
void add_normal_distribution_disturb(double* next_index,int k_num,int D,double* Xu,double* Xl) {
	for (int j = 0; j <k_num; j++) { //当前初始变异个体的每一个簇中心，即每一维都会添加高斯扰动
		for (int l = 0; l < D; l++) {
			double gauss_index = next_index[j * D + l] + sampleNormal(0, 0.1) * (Xu[l] - Xl[l]);//sampleNormal是生成高斯随机数生成器，0是均值，0.1是标准差。高斯分布纵坐标表示概率，横坐标表示数值，而由该函数图像可以得出，该生成器生成0±0.1附近的数概率最大。
			if (gauss_index > Xu[l]) {
				next_index[j * D + l] = Xu[l]; //越界就取边界值
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
* 杂交操作,在变异向量中选择杂交个体，并返回用于杂交的簇数
* next_index是当前的变异个体
* k_num是变异个体的簇数
* in_area是保存从当前变异个体中获取的簇中心
* D是数据维度
* CR是变异概率
* 该方法返回从变异个体中选择的杂交簇数
*/
int selectClustersFromMutateInd(double* next_index,int k_num,double* in_area,int D,double CR) {
	int cr_lenth = 0, swap = 0, swap1 = 0; 
	int cr_n;
	do {
		cr_lenth = cr_lenth + 1;
	} while (URAND < CR && cr_lenth < k_num);//在杂交前，先计算变异个体中用来交叉的簇中心数
	cr_n = (int)(k_num * URAND); //随机一个杂交点，取值范围是0~rand_basic-1  
	//step 2 of the crossover: determine the subspace of crossover
	int rand1;
	   //the number of points that outside and inside the swap area
	if ((cr_n + cr_lenth - 1) < k_num) {//rand_basic是当前变异个体的簇中心数，cr_n+cr_length是当前变异个体最后一个簇的位置
		for (int j = cr_n; j < cr_n + cr_lenth; j++) {
			for (int l = 0; l < D; l++)
				in_area[(j - cr_n) * D + l] = next_index[j * D + l]; //in_area二维数组保存的是每个变异个体待杂交的簇中心，而且是从0开始连续保存
		}
	}
	else {//如果从杂交开始点cr_n开始杂交一直杂交cr_length个杂交向量，而杂交结束点已经超出了当前变异个体的簇数  
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
* 将out_area数组的信息复制到next_index数组中得到试验向量
* k_num是指out_area中的簇数
*/
void get_trailVector(double* out_area, int k_num,double* next_index,double*next_param,int D) {
	for (int j = 0; j < k_num * D; j++)
		next_index[j] = out_area[j]; //把out_area数组中的簇复制到保存所有变异个体的数组，即得到试验向量的数组。
	next_param[0] = k_num;
}

/*
* 进行子空间杂交
* cr_lenth 是变异个体子空间的簇数
* in_area 是从变异个体中取到的用于杂交的簇
* out_area 杂交后的簇先保存于该数组，最后再保存到next_index数组
* p_index 保存了当前目标个体
* p_param 保存了当前个体的簇数和适应度
* next_index 用于保存杂交后的簇，即作为试验向量
* next_param 用来保存实验向量的簇数
* D 数据的维度
* 
*/
void subspace_cross(int cr_lenth, double* in_area, double* out_area, double* p_index, double* p_param, double* next_index, double* next_param,int D) {
	double* min_b;
	double* max_b;
	double minb, maxb;
	malloc1D(min_b, D);
	malloc1D(max_b, D);
	if (cr_lenth == 1) {//需要杂交的簇只有一个的情况，
		for (int j = 0; j < 20 * D; j++)
			out_area[j] =p_index[j]; //先把当前目标个体的所有簇保存到数组out_area中。
		//如果对位的目标个体是满簇情况，即20个簇，那么子空间交叉时，只需要把需要杂交的簇替换对位目标个体随机的一个簇。  
		//if ((int)*(p2 + i * 2) == 20) {
		int rand1 = (int)((int)p_param[0] * URAND);
		for (int j = 0; j < D; j++)
			out_area[rand1 * D + j] = in_area[j]; //如果cr_lenth==1,即交叉的簇只有一个，那么数组in_area[i]只保存了一个簇，然后把这个簇替换到out_area数组的随机一个位置
		//for (int j = 0; j < 20 * D; j++)
		//	next_index[j] = out_area[j]; //把out_area数组中的簇复制到保存所有变异个体的数组，即得到试验向量的数组。
		//next_param[0] = p_param[0]; //next_param数组保存的是第i个变异个体的簇数，这里因为杂交的簇只有一个因此第i个变异个体（试验向量）的簇数与对位的目标个体的簇数保持一致
		get_trailVector(out_area, (int)p_param[0], next_index, next_param, D);
	}
	else {
		for (int j = 0; j <D; j++) {
			minb = DBL_MAX;
			maxb = DBL_MIN;
			for (int i = 0; i < cr_lenth; i++) {//将需要杂交的所有簇，叠加起来，每维求平均值之后，就可以得到均值向量
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
				//满足条件，则说明当前的簇不在杂交子空间范围内
				if (p_index[j * D + l] < min_b[l]|| p_index[j*D+l]>max_b[l]) {   //the center can be changed
					swap_if = 1;
					break;
				}
			}
			if (swap_if == 1) {//把不需要被覆盖的簇保存起来，即这部分簇是目标个体中不需要被变异个体替换的簇。随后这部分簇与变异个体中被选中的杂交簇结合成试验向量
				for (int l = 0; l < D; l++)
					out_area[swap1 * D + l] = p_index[j*D+l];
				swap1 = swap1 + 1;
			}
		}
		swap = cr_lenth + swap1;//swap1是目标个体中不需要被交换的簇数，cr_lenth是指变异个体子空间（数组in_area）的簇数，即用来覆盖的簇数

		//step 3 of the crossover : subarea swap
		if (swap <= MAX_CLUSTER && swap >= MIN_CLUSTER) {
			for (int j = 0; j < cr_lenth * D; j++)
				next_index[j] = in_area[j]; //先把变异个体子空间中的簇，插入到变异向量中
			for (int j = cr_lenth * D; j < swap * D; j++)
				next_index[j] = out_area[j - cr_lenth * D];//再把目标个体中不需要被覆盖的簇插入到变异向量中，该操作完成即原本保存变异向量的数组就变成保存试验向量了。
			next_param[0] = swap;
		}else {
			/*当杂交后的簇数大于20或小于2的时候，进行如下处理*/
			dealWithClusterNumNotCorrect(p_index, p_param, next_index, next_param, D);
		}
	}

	free(min_b);
	free(max_b);
}
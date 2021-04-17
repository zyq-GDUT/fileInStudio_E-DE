/*******************************************************************************/
/* This is a simple implementation of Elastic Differential Evolution (E-DE).   */
/* The codes are written in C.                                                 */
/* For any questions, please contact J. Chen (junxianchen001@gmail.com).       */
/* E-DE_1.0, Edited on August, 2019.                                           */
/*******************************************************************************/

// add your header files here 
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
#include"testtool.h"




int min_cluster_num = MIN_CLUSTER;
int max_cluster_num = MAX_CLUSTER;


int fixEmptyClusterNum = 0;

//用来保存每次测试的数据,trailData定义在了头文件ede.h
trailData traildata = { "", {0},{0},{0}, {0},{0},0,0 };


//// change any of these parameters to match your needs 
//#define URAND  ((double)rand()/((double)RAND_MAX+1.0)) //产生的随机数是[0,1),rand()产生0-RAND_MAX范围的数

// declaration of functions used by this algorithm 

// 1. load and arrange the data 
double** data;


// 2. functions to implement the E-DE algorithm 
double func(double k_num, double* p_pars, double* p_index, int* in_clu, int N, int D,int useIndex);
//double getDistance(double* avector, double* bvector, int n);
int e_de(double* p, double* p2, int N, int D, int Gmax, double* Xl, double* Xu,int usedindex, int** in_cluster);

/********************************************************************
* func_name: fitness calculation
* input: the cluster number, population info, data info
* descript: calculating the fitness of an individual
* this takes a user defined function
********************************************************************/

double func(double k_num, double* p_pars, double* p_index, int* in_clu, int N, int D,int usedIndex) {
	/* define your testing function here */
	double fes;
	if (usedIndex==I_INDEX) {
		fes = I_index((int)k_num, p_pars, p_index, in_clu, N, D);
	}
	else {
		fes = DB_index((int)k_num, p_pars, p_index, in_clu, N, D);
	}
	

	return fes;
}

///*****************************************************************
//* func_name: the gaussian generator
//* input: mean and standard deviation
//* descript: for the crossover operation
//*****************************************************************/
//double sampleNormal(double mean, double sigma)
//{
//	 // apply the unix time to set the seed of rand
//	 static boost::mt19937 rng(static_cast< unsigned >(std::time(0)));
//
//	// select the gaussian random distribution
//	 boost::normal_distribution< double > norm_dist(mean, sigma);
//
//	 // generate the random generator
//	 boost::variate_generator< boost::mt19937&, boost::normal_distribution< double > > normal_sampler(rng, norm_dist);
//	
//	 return normal_sampler();
//	
//	//double gassrandN = gaussrand(0.0, 0.1);
//	//return gassrandN;
//}




/*****************************************************************
* func_name: E-DE algorithm
* input: the evolution population, the data number N, dimension D
*        the lower and upper bound Xl, Xu
* descript: do the mutation and crossover opertors, and select the
*           best results for the next generation
*****************************************************************/
int e_de(double* p, double* p2, int N, int D, int Gmax, double* Xl, double* Xu,int usedIndex,int **in_cluster) {

	
	int i, k, r1, r2, r3, r4, r1_full, r3_full;
	int* r1_rand, * r2_rand, * r3_rand, ** in_cluster2;
	int numofE = 0, index = 0,NP, swap, swap1, swap_if, cr_lenth, cr_n;
	double F = 0.5, CR = 0.4, MU, gauss_index;
	double** next_index, ** next_param, distance;
	double** out_area, ** in_area, * dist, * center,  min;
	double  best_val ;
	double best_val2 ;
	NP = getNP(D, 10);
	/*outputPopulationInfo(p, p2, D, NP, max_cluster_num);*/

	////使用DB指标
	//if (usedIndex == DB_INDEX) {
	//     best_val = getMinFitness(p2, NP, GET_FITNESS);
	//	 best_val2 = getMinFitness(p2, NP, GET_CLUSTER_NUM);
	//}
	////使用I指标
	//else {
	//	 best_val = getMaxFitness(p2, NP, GET_FITNESS);
	//	 best_val2 = getMaxFitness(p2, NP, GET_CLUSTER_NUM);
	//}
	//根据当前的评价指标获取最优适应度和最优簇数，分别保存到best_val和best_val2
	getBestValByIndex(best_val, best_val2, p2, NP, usedIndex);
	
	malloc1D(dist, D);
	malloc1D(center, D);
	malloc1E(r1_rand, max_cluster_num);  //最大20个簇，这三个数组与簇数相关。因此这三个随机数组大小设置为20. 
	malloc1E(r2_rand, max_cluster_num);
	malloc1E(r3_rand, max_cluster_num);

	malloc2E(out_area, NP, max_cluster_num * D);
	malloc2E(in_area, NP, max_cluster_num * D);
	malloc2E(next_index, NP, max_cluster_num * D);//NP行20*D列，用来保存每个变异个体的所有簇中心
	malloc2E(next_param, NP, 2); //NP行2列，用来保存每个变异个体的簇中心数
	malloc2D(in_cluster2, NP, N);//用来标记每个个体中，样本所属的簇，每个个体都是一种聚类方案，而个体只是保存了一种聚类方案的簇中心。因此需要用一个二维数组来表明每种聚类方案的样本归属情况

	//for (i = 0; i < NP; i++) {
	//	for (int j = 0; j < 2; j++) {
	//		next_param[i][j] = 0;
	//	}
	//	for (int j = 0; j < N; j++) {
	//		in_cluster2[i][j] = 0;
	//	}
	//	for (int j = 0; j <max_cluster_num * D; j++) {
	//		out_area[i][j] = 0;
	//		in_area[i][j] = 0;
	//		next_index[i][j] = 0; //用来保存变异个体
	//	}
	//}
	//
	initialCoreVriable(next_param, in_cluster2, out_area, in_area, next_index, NP, max_cluster_num, D, N);
	//for (i = 0; i < D; i++) {
	//	dist[i] = 0; //距离数组，用来记录杂交操作时，杂交子空间中，每一维的取值范围
	//	center[i] = 0; //用来记录子空间中的簇中心
	//}
	initial_dist_center(dist, center, D);
	outputPopulationInfo2(p, p2, D, NP, max_cluster_num, in_cluster, N, 1, best_val, best_val2, 20000);
	//一共迭代Gmax次
	Gmax = 1000;
	for (k = 0; k < Gmax; k++)      //Gmax denotes the maximum iteration times
	{
		
		for (i = 0; i < NP; i++)
		{
			getRandomNum(r1, r2, r3, i, NP);
			//step 2 of the mutation: determine the cluster number of mutant vector
			//这里应该是获取随机三个个体的簇中心数，这里p2应该也是二维数组，一共有np行和2列，用来保存种群中每个个体的簇数。
			int mutate_d1 = (int)*(p2 + r1 * 2);  //获取目标个体r1的簇中心数
			int mutate_d2 = (int)*(p2 + r2 * 2);	//获取目标个体r2的簇中心数
			int mutate_d3 = (int)*(p2 + r3 * 2);	//获取目标个体r3的簇中心数
			int mutate_u = mutate_d1 + (int)(F * ((double)mutate_d2 - (double)mutate_d3));
			mutate_u = cross_border_process(mutate_u, max_cluster_num, min_cluster_num);
			next_param[i][0] = mutate_u; //变异向量中有NP个个体，该数组是保存变异向量中每个个体的簇中心数
			//生成mutate_u个簇中心保存到next_index中
			mutate_process(p, next_index[i], mutate_u, mutate_d1, mutate_d2, mutate_d3, r1, r2, r3, D, N);
			//step 4 of the mutation: fine tune the cluster centroids
			MU = 0.1 - 0.06 * k / ((double)Gmax - 1.0); //上面的操作是得到初始的变异个体，该步骤是对初始变异个体按一定概率进行处理，每轮迭代添加扰动处理的变异个体数大概是NP*MU
			if (URAND < MU) {//每次都有概率MU对初始变异个体添加高斯扰动
				add_normal_distribution_disturb( next_index[i], next_param[i][0], D, Xu, Xl);
			}
			int rand1;
			swap = 0, swap1 = 0;    //the number of points that outside and inside the swap area
			/*杂交操作*/
			//Crossover operator
			//step 1 of the crossover: determine the length of crossover
			int rand_basic = (int)next_param[i][0];//获取第i个变异个体的簇中心数
			cr_lenth = 0;
			do {
				cr_lenth = cr_lenth + 1;
			} while (URAND < CR && cr_lenth < rand_basic);//在杂交前，先计算变异个体中用来交叉的簇中心数
			cr_n = (int)(rand_basic * URAND); //随机一个杂交点，取值范围是0~rand_basic-1  
			//step 2 of the crossover: determine the subspace of crossover
			
			if ((cr_n + cr_lenth - 1) < rand_basic) {//rand_basic是当前变异个体的簇中心数，cr_n+cr_length是当前变异个体最后一个簇的位置
				for (int j = cr_n; j < cr_n + cr_lenth; j++) {
					for (int l = 0; l < D; l++)
						in_area[i][(j - cr_n) * D + l] = next_index[i][j * D + l]; //in_area二维数组保存的是每个变异个体待杂交的簇中心，而且是从0开始连续保存
				}
			}
			else {//如果从杂交开始点cr_n开始杂交一直杂交cr_length个杂交向量，而杂交结束点已经超出了当前变异个体的簇数  
				for (int j = cr_n; j < rand_basic; j++) {
					for (int l = 0; l < D; l++)
						in_area[i][(j - cr_n) * D + l] = next_index[i][j * D + l];
				}
				for (int j = 0; j < (cr_n + cr_lenth - rand_basic); j++) {
					for (int l = 0; l < D; l++)
						in_area[i][(j + rand_basic - cr_n) * D + l] = next_index[i][j * D + l];
				}
			}
			//int cr_lenth=selectClustersFromMutateInd(next_index[i], (int)next_param[i][0],  in_area[i],  D, CR);

			/*进行子空间交叉*/
			if (cr_lenth == 1) {//需要杂交的簇只有一个的情况，
				for (int j = 0; j < 20 * D; j++)
					out_area[i][j] = *(p + i * 20 * D + j); //先把当前目标个体的所有簇保存到数组out_area中。
				//如果对位的目标个体是满簇情况，即20个簇，那么子空间交叉时，只需要把需要杂交的簇替换对位目标个体随机的一个簇。  
				//if ((int)*(p2 + i * 2) == 20) {
					rand1 = (int)(((int)*(p2 + i * 2)) * URAND);
					for (int j = 0; j < D; j++)
						out_area[i][rand1 * D + j] = in_area[i][j]; //如果cr_lenth==1,即交叉的簇只有一个，那么数组in_area[i]只保存了一个簇，然后把这个簇替换到out_area数组的随机一个位置
					for (int j = 0; j < 20 * D; j++)
						next_index[i][j] = out_area[i][j]; //把out_area数组中的簇复制到保存所有变异个体的数组，即得到试验向量的数组。
					next_param[i][0] = *(p2 + i * 2); //next_param数组保存的是第i个变异个体的簇数，这里因为杂交的簇只有一个因此第i个变异个体（试验向量）的簇数与对位的目标个体的簇数保持一致
				//}
				////对位的目标个体的簇数，小于20个，即没满簇，那么子空间杂交，只需要把待杂交的簇插入到目标个体的末尾。
				//else {
				//	for (int j = 0; j < D; j++)
				//		out_area[i][((int)*(p2 + i * 2)) * D + j] = in_area[i][j];
				//	for (int j = 0; j < 20 * D; j++)
				//		next_index[i][j] = out_area[i][j];
				//	next_param[i][0] = *(p2 + i * 2) + 1;
				//}
			}
			//需要杂交的簇数大于1，即需要使用子空间杂交法，首先需要得到子空间的均值向量（中心向量）
			else {
				for (int j = 0; j < D; j++)
					center[j] = 0; //子空间的中心向量初始值为全0

				//average all node to get the swap center
				for (int j = 0; j < cr_lenth; j++) {
					for (int l = 0; l < D; l++) {//将需要杂交的所有簇，叠加起来，每维求平均值之后，就可以得到均值向量
						center[l] = center[l] + in_area[i][j * D + l];
					}
				}

				for (int j = 0; j < D; j++) {//保存子空间中的第一个簇与最后一个簇的对位距离
					dist[j] = fabs(in_area[i][0 * D + j] - in_area[i][(cr_lenth - 1) * D + j]) / 2; //first and last node
					center[j] = center[j] / cr_lenth;  //average of all node
				}


				for (int j = 0; j < (int)*(p2 + i * 2); j++) {
					swap_if = 0;
					for (int l = 0; l < D; l++) {
						//满足条件，则说明当前的簇不在杂交子空间范围内
						if (fabs(*(p + i * 20 * D + j * D + l) - center[l]) > dist[l]) {   //the center can be changed
							swap_if = 1;
							break;
						}
					}
					if (swap_if == 1) {//把不需要被覆盖的簇保存起来，即这部分簇是目标个体中不需要被变异个体替换的簇。随后这部分簇与变异个体中被选中的杂交簇结合成试验向量
						for (int l = 0; l < D; l++)
							out_area[i][swap1 * D + l] = *(p + i * 20 * D + j * D + l);
						swap1 = swap1 + 1;
					}
				}
				swap = cr_lenth + swap1;//swap1是目标个体中不需要被交换的簇数，cr_lenth是指变异个体子空间（数组in_area）的簇数，即用来覆盖的簇数

				//step 3 of the crossover : subarea swap
				//因为每个簇至少要有2个样本，因此如果出现平均每个簇的样本数小于2，则该簇数就认为过大。
				if (swap <=max_cluster_num && swap >= min_cluster_num) {
					for (int j = 0; j < cr_lenth * D; j++)
						next_index[i][j] = in_area[i][j]; //先把变异个体子空间中的簇，插入到变异向量中
					for (int j = cr_lenth * D; j < swap * D; j++)
						next_index[i][j] = out_area[i][j - cr_lenth * D];//再把目标个体中不需要被覆盖的簇插入到变异向量中，该操作完成即原本保存变异向量的数组就变成保存试验向量了。
					next_param[i][0] = swap;
				}
				else {
					/*当杂交后的簇数大于20或小于2的时候，进行如下处理*/
					dealWithClusterNumNotCorrect(p + i * max_cluster_num * D, p2 + i * 2, next_index[i], next_param[i], D);
				}
				
			}
			//subspace_cross(cr_lenth, in_area[i], out_area[i], (p + i * MAX_CLUSTER * D), (p2 + i * 2), next_index[i], next_param[i], D);

			//Selection operator,将样本分配到具体的簇中
			//step 1 of the selection: assign the data object
			for (int j = 0; j < N; j++) {
				min = DBL_MAX;
				for (int l = 0; l < ((int)next_param[i][0]); l++) {
					distance = getDistance(data[j], &next_index[i][l * D], D);
					if (distance < min) {
						min = distance;
						in_cluster2[i][j] = l;    //record the data index
					}
				}
			}

			dealwith_emptyCluster(in_cluster2[i], N, next_param[i][0], data, D, next_index[i]);

			//step 2 of the selection: calculate the fitness，next_param[i][0]记录了当前试验向量的簇数，data是样本数据集，next_index[i]保存了当前这个试验向量的各个簇中心，in_cluster2[i]数据集中每个数据属于哪个簇。
			next_param[i][1] = func(next_param[i][0], *data, next_index[i], in_cluster2[i], N, D,usedIndex);
			numofE = numofE + 1;

	
			//I指标适应度越大越好
			if (usedIndex == I_INDEX) {
				if (next_param[i][1] > *(p2 + i * 2 + 1)) {
					for (int j = 0; j < 20 * D; j++) {
						*(p + i * 20 * D + j) = next_index[i][j];
					}
					*(p2 + i * 2 + 0) = next_param[i][0];
					*(p2 + i * 2 + 1) = next_param[i][1];
					for (int k = 0; k < N; k++) {
						in_cluster[i][k] = in_cluster2[i][k];
					}
				}
				
				if (*(p2 + i * 2 + 1) > best_val) {
					best_val = *(p2 + i * 2 + 1);
					best_val2 = *(p2 + i * 2 + 0);
					index = i;//记录最优个体的位置
				}
			}
			//DB指标适应度越小越好
			else {
				if (next_param[i][1] < *(p2 + i * 2 + 1)) {
					for (int j = 0; j < 20 * D; j++) {
						*(p + i * 20 * D + j) = next_index[i][j];
					}
					*(p2 + i * 2 + 0) = next_param[i][0];
					*(p2 + i * 2 + 1) = next_param[i][1];
					for (int k = 0; k < N; k++) {
						in_cluster[i][k] = in_cluster2[i][k];
					}
				}
				
				if (*(p2 + i * 2 + 1) < best_val) {
					best_val = *(p2 + i * 2 + 1);//最优的适应度
					best_val2 = *(p2 + i * 2 + 0);//最优的簇数
					index = i;//记录最优个体的位置
				}

			}
			
		}
	/*	printf("\n---------第%d代-最优适应度为%g-对应的簇数为%d------------------\n", k+1, best_val, (int)best_val2);
		outputPopulationInfo(p, p2, D, NP, max_cluster_num, in_cluster, N);*/
		outputPopulationInfo2(p, p2, D, NP, max_cluster_num, in_cluster, N,k+1,best_val,best_val2,20000);
		
	}
	


	free(dist);
	free(center);
	free(r1_rand);
	free(r2_rand);
	free(r3_rand);
	freeTwoDimArr_double(out_area, NP);
	freeTwoDimArr_double(in_area, NP);
	freeTwoDimArr_double(next_index, NP);
	freeTwoDimArr_double(next_param, NP);
	freeTwoDimArr_int(in_cluster2, NP);
	return index;
}



void ede_initialization(const char* filename,int usedIndex) {
	//srand((unsigned int)(time(NULL)));
	int i,  D, N, Gmax, NP, best = 0, * popul_rand, ** in_cluster;
	double   * uk, * lk, ** popul_index, ** popul_param;
	int* originDataClusterInfo;//用来保存数据集本身的聚类情况
	data = loadData2(filename, N, D);   //load data from the text file,将数据集存到了二维数组data中
	originDataClusterInfo=getDataOriginClusterInfo(filename, N);
	set_traildataInfo(traildata, originDataClusterInfo[N], filename, usedIndex);
	printf("文件%s共%d个数据，维度是%d\n", filename, N, D);
	max_cluster_num = getMaxClusterNum(N, max_cluster_num);//每种簇数的聚类方案都要满足，每个簇至少有两个样本。因此最大的簇数20满足了该要求，2-19这些簇数都能都满足该要求
	NP = getNP(D, 10);//维数越大，问题的解空间就越大，因此种群的个数就应该更大
	Gmax = getGmax(NP, 1000000);//每一代都需要进行NP次评估，因此该式子能得到最大的迭代代数
	printf("The times of iteration(Gmax):%d\n", Gmax);
	malloc1D(uk, D);
	malloc1D(lk, D);
	malloc1E(popul_rand, N);
	malloc2D(in_cluster, NP, N);        //population index，二维数组in_cluster用来记录当前样本在哪个簇的情况
	malloc2E(popul_index, NP,max_cluster_num * D);  //population cluster centroids
	malloc2E(popul_param, NP, 2);       //population info(cluster number, fitness value)
	initial_two_Dim_intarr(in_cluster, NP, N, 0);
	getDataDim_max_min(data, uk, lk, D, N);//获取数据每个维度的最大值最小值并保存到数组uk,lk
	//population initialization
	for (i = 0; i < NP; i++) {
		initial_individual(data, popul_index[i], popul_param[i], in_cluster[i], N, D, max_cluster_num, min_cluster_num);
		popul_param[i][1] = func(popul_param[i][0], *data, popul_index[i], in_cluster[i], N, D,usedIndex);//计算当前个体的适应度	
	}
	best = e_de(*popul_index, *popul_param, N, D, Gmax, lk, uk,usedIndex,in_cluster); //run the E-DE algorithm

	//计算本次聚类的簇内距离，簇间距离，最优适应度值，标准差，R指标
	popul_index[best]= count_meanV(*data, in_cluster[best], popul_param[best][0], N, D);//重新计算获得真正的簇中心
	double ARI = getARI(originDataClusterInfo, in_cluster[best], popul_param[best][0], N);
	double IntraDistance = getIntraDist(data, popul_index[best], popul_param[best],in_cluster[best],N,D);
	double InterDistance = getInterDist(popul_index[best], popul_param[best], D);

	recordTraildata(traildata, ARI, InterDistance, IntraDistance, popul_param[best][1], popul_param[best][0]);
	free(uk);
	free(lk);
	free(popul_rand);
	free(originDataClusterInfo);
	freeTwoDimArr_int(in_cluster,NP);
	freeTwoDimArr_double(popul_param, NP);
	freeTwoDimArr_double(popul_index, NP);
	freeTwoDimArr_double(data, N);

}






void test_ede(int usedIndex,const char* filename){
	char* outputfile;
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		 outputfile = getOutputFilename(filename, usedIndex,i+1);	
		 freopen((const char* )outputfile, "w", stdout);
		 printf("当前操作的文件为:%s\n", filename);
		 ede_initialization(filename,usedIndex);
		
	}
	free(outputfile);


}
void ede_two_index(const char* filename) {
	for (int i = 0; i < 2; i++) {
		char* indexname = getIndexName(i);
		freopen("CON", "w", stdout);
		printf("正在使用指标%s进行处理...\n", indexname);
		test_ede(i, filename);
		free(indexname);
	}
	

}

char* showmenu() {
	 char* filename = malloc1Char(100);
	int i;
	
	printf("**************************\n");
	printf("***** 1. ecoli.txt *****\n");
	printf("***** 2. iris.txt *****\n");
	printf("***** 3. seeds.txt *****\n");
	printf("***** 4. zoo.txt *****\n");
	printf("**************************\n");
	printf("请输入你要处理的文件编号：\n");
	scanf("%d", &i);
	while (i < 1 || i>4) {
		printf("输入有误请重新输入(1~4)\n");
		scanf("%d", &i);
	}
	switch (i) {
		case 1:
			strcpy(filename, "ecoli.txt");
			break;
		case 2:
			strcpy(filename, "iris.txt");
			break;
		case 3:
			strcpy(filename, "seeds.txt");
			break;
		case 4:
			strcpy(filename, "zoo.txt");
			break;
		default:
			break;

	}
	printf("当前处理的文件是%s\n程序执行中,请勿关闭当前窗口....\n", filename);
	return filename;
}


int  main()
{

	////依次用两个指标进行处理
	//srand((unsigned int)(time(NULL)));
	char* filename = showmenu();
	ede_two_index(filename);
	free(filename);
	freopen("CON", "w", stdout);
	printf("程序执行完毕，");
	fclose(stdout);
	system("pause");
	return 0;
}



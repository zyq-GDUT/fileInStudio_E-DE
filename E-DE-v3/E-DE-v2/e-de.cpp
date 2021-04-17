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

//��������ÿ�β��Ե�����,trailData��������ͷ�ļ�ede.h
trailData traildata = { "", {0},{0},{0}, {0},{0},0,0 };


//// change any of these parameters to match your needs 
//#define URAND  ((double)rand()/((double)RAND_MAX+1.0)) //�������������[0,1),rand()����0-RAND_MAX��Χ����

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

	////ʹ��DBָ��
	//if (usedIndex == DB_INDEX) {
	//     best_val = getMinFitness(p2, NP, GET_FITNESS);
	//	 best_val2 = getMinFitness(p2, NP, GET_CLUSTER_NUM);
	//}
	////ʹ��Iָ��
	//else {
	//	 best_val = getMaxFitness(p2, NP, GET_FITNESS);
	//	 best_val2 = getMaxFitness(p2, NP, GET_CLUSTER_NUM);
	//}
	//���ݵ�ǰ������ָ���ȡ������Ӧ�Ⱥ����Ŵ������ֱ𱣴浽best_val��best_val2
	getBestValByIndex(best_val, best_val2, p2, NP, usedIndex);
	
	malloc1D(dist, D);
	malloc1D(center, D);
	malloc1E(r1_rand, max_cluster_num);  //���20���أ������������������ء������������������С����Ϊ20. 
	malloc1E(r2_rand, max_cluster_num);
	malloc1E(r3_rand, max_cluster_num);

	malloc2E(out_area, NP, max_cluster_num * D);
	malloc2E(in_area, NP, max_cluster_num * D);
	malloc2E(next_index, NP, max_cluster_num * D);//NP��20*D�У���������ÿ�������������д�����
	malloc2E(next_param, NP, 2); //NP��2�У���������ÿ���������Ĵ�������
	malloc2D(in_cluster2, NP, N);//�������ÿ�������У����������Ĵأ�ÿ�����嶼��һ�־��෽����������ֻ�Ǳ�����һ�־��෽���Ĵ����ġ������Ҫ��һ����ά����������ÿ�־��෽���������������

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
	//		next_index[i][j] = 0; //��������������
	//	}
	//}
	//
	initialCoreVriable(next_param, in_cluster2, out_area, in_area, next_index, NP, max_cluster_num, D, N);
	//for (i = 0; i < D; i++) {
	//	dist[i] = 0; //�������飬������¼�ӽ�����ʱ���ӽ��ӿռ��У�ÿһά��ȡֵ��Χ
	//	center[i] = 0; //������¼�ӿռ��еĴ�����
	//}
	initial_dist_center(dist, center, D);
	outputPopulationInfo2(p, p2, D, NP, max_cluster_num, in_cluster, N, 1, best_val, best_val2, 20000);
	//һ������Gmax��
	Gmax = 1000;
	for (k = 0; k < Gmax; k++)      //Gmax denotes the maximum iteration times
	{
		
		for (i = 0; i < NP; i++)
		{
			getRandomNum(r1, r2, r3, i, NP);
			//step 2 of the mutation: determine the cluster number of mutant vector
			//����Ӧ���ǻ�ȡ�����������Ĵ�������������p2Ӧ��Ҳ�Ƕ�ά���飬һ����np�к�2�У�����������Ⱥ��ÿ������Ĵ�����
			int mutate_d1 = (int)*(p2 + r1 * 2);  //��ȡĿ�����r1�Ĵ�������
			int mutate_d2 = (int)*(p2 + r2 * 2);	//��ȡĿ�����r2�Ĵ�������
			int mutate_d3 = (int)*(p2 + r3 * 2);	//��ȡĿ�����r3�Ĵ�������
			int mutate_u = mutate_d1 + (int)(F * ((double)mutate_d2 - (double)mutate_d3));
			mutate_u = cross_border_process(mutate_u, max_cluster_num, min_cluster_num);
			next_param[i][0] = mutate_u; //������������NP�����壬�������Ǳ������������ÿ������Ĵ�������
			//����mutate_u�������ı��浽next_index��
			mutate_process(p, next_index[i], mutate_u, mutate_d1, mutate_d2, mutate_d3, r1, r2, r3, D, N);
			//step 4 of the mutation: fine tune the cluster centroids
			MU = 0.1 - 0.06 * k / ((double)Gmax - 1.0); //����Ĳ����ǵõ���ʼ�ı�����壬�ò����ǶԳ�ʼ������尴һ�����ʽ��д���ÿ�ֵ�������Ŷ�����ı�������������NP*MU
			if (URAND < MU) {//ÿ�ζ��и���MU�Գ�ʼ���������Ӹ�˹�Ŷ�
				add_normal_distribution_disturb( next_index[i], next_param[i][0], D, Xu, Xl);
			}
			int rand1;
			swap = 0, swap1 = 0;    //the number of points that outside and inside the swap area
			/*�ӽ�����*/
			//Crossover operator
			//step 1 of the crossover: determine the length of crossover
			int rand_basic = (int)next_param[i][0];//��ȡ��i���������Ĵ�������
			cr_lenth = 0;
			do {
				cr_lenth = cr_lenth + 1;
			} while (URAND < CR && cr_lenth < rand_basic);//���ӽ�ǰ���ȼ�������������������Ĵ�������
			cr_n = (int)(rand_basic * URAND); //���һ���ӽ��㣬ȡֵ��Χ��0~rand_basic-1  
			//step 2 of the crossover: determine the subspace of crossover
			
			if ((cr_n + cr_lenth - 1) < rand_basic) {//rand_basic�ǵ�ǰ�������Ĵ���������cr_n+cr_length�ǵ�ǰ����������һ���ص�λ��
				for (int j = cr_n; j < cr_n + cr_lenth; j++) {
					for (int l = 0; l < D; l++)
						in_area[i][(j - cr_n) * D + l] = next_index[i][j * D + l]; //in_area��ά���鱣�����ÿ�����������ӽ��Ĵ����ģ������Ǵ�0��ʼ��������
				}
			}
			else {//������ӽ���ʼ��cr_n��ʼ�ӽ�һֱ�ӽ�cr_length���ӽ����������ӽ��������Ѿ������˵�ǰ�������Ĵ���  
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

			/*�����ӿռ佻��*/
			if (cr_lenth == 1) {//��Ҫ�ӽ��Ĵ�ֻ��һ���������
				for (int j = 0; j < 20 * D; j++)
					out_area[i][j] = *(p + i * 20 * D + j); //�Ȱѵ�ǰĿ���������дر��浽����out_area�С�
				//�����λ��Ŀ������������������20���أ���ô�ӿռ佻��ʱ��ֻ��Ҫ����Ҫ�ӽ��Ĵ��滻��λĿ����������һ���ء�  
				//if ((int)*(p2 + i * 2) == 20) {
					rand1 = (int)(((int)*(p2 + i * 2)) * URAND);
					for (int j = 0; j < D; j++)
						out_area[i][rand1 * D + j] = in_area[i][j]; //���cr_lenth==1,������Ĵ�ֻ��һ������ô����in_area[i]ֻ������һ���أ�Ȼ���������滻��out_area��������һ��λ��
					for (int j = 0; j < 20 * D; j++)
						next_index[i][j] = out_area[i][j]; //��out_area�����еĴظ��Ƶ��������б����������飬���õ��������������顣
					next_param[i][0] = *(p2 + i * 2); //next_param���鱣����ǵ�i���������Ĵ�����������Ϊ�ӽ��Ĵ�ֻ��һ����˵�i��������壨�����������Ĵ������λ��Ŀ�����Ĵ�������һ��
				//}
				////��λ��Ŀ�����Ĵ�����С��20������û���أ���ô�ӿռ��ӽ���ֻ��Ҫ�Ѵ��ӽ��Ĵز��뵽Ŀ������ĩβ��
				//else {
				//	for (int j = 0; j < D; j++)
				//		out_area[i][((int)*(p2 + i * 2)) * D + j] = in_area[i][j];
				//	for (int j = 0; j < 20 * D; j++)
				//		next_index[i][j] = out_area[i][j];
				//	next_param[i][0] = *(p2 + i * 2) + 1;
				//}
			}
			//��Ҫ�ӽ��Ĵ�������1������Ҫʹ���ӿռ��ӽ�����������Ҫ�õ��ӿռ�ľ�ֵ����������������
			else {
				for (int j = 0; j < D; j++)
					center[j] = 0; //�ӿռ������������ʼֵΪȫ0

				//average all node to get the swap center
				for (int j = 0; j < cr_lenth; j++) {
					for (int l = 0; l < D; l++) {//����Ҫ�ӽ������дأ�����������ÿά��ƽ��ֵ֮�󣬾Ϳ��Եõ���ֵ����
						center[l] = center[l] + in_area[i][j * D + l];
					}
				}

				for (int j = 0; j < D; j++) {//�����ӿռ��еĵ�һ���������һ���صĶ�λ����
					dist[j] = fabs(in_area[i][0 * D + j] - in_area[i][(cr_lenth - 1) * D + j]) / 2; //first and last node
					center[j] = center[j] / cr_lenth;  //average of all node
				}


				for (int j = 0; j < (int)*(p2 + i * 2); j++) {
					swap_if = 0;
					for (int l = 0; l < D; l++) {
						//������������˵����ǰ�Ĵز����ӽ��ӿռ䷶Χ��
						if (fabs(*(p + i * 20 * D + j * D + l) - center[l]) > dist[l]) {   //the center can be changed
							swap_if = 1;
							break;
						}
					}
					if (swap_if == 1) {//�Ѳ���Ҫ�����ǵĴر������������ⲿ�ִ���Ŀ������в���Ҫ����������滻�Ĵء�����ⲿ�ִ����������б�ѡ�е��ӽ��ؽ�ϳ���������
						for (int l = 0; l < D; l++)
							out_area[i][swap1 * D + l] = *(p + i * 20 * D + j * D + l);
						swap1 = swap1 + 1;
					}
				}
				swap = cr_lenth + swap1;//swap1��Ŀ������в���Ҫ�������Ĵ�����cr_lenth��ָ��������ӿռ䣨����in_area���Ĵ��������������ǵĴ���

				//step 3 of the crossover : subarea swap
				//��Ϊÿ��������Ҫ��2������������������ƽ��ÿ���ص�������С��2����ô�������Ϊ����
				if (swap <=max_cluster_num && swap >= min_cluster_num) {
					for (int j = 0; j < cr_lenth * D; j++)
						next_index[i][j] = in_area[i][j]; //�Ȱѱ�������ӿռ��еĴأ����뵽����������
					for (int j = cr_lenth * D; j < swap * D; j++)
						next_index[i][j] = out_area[i][j - cr_lenth * D];//�ٰ�Ŀ������в���Ҫ�����ǵĴز��뵽���������У��ò�����ɼ�ԭ�������������������ͱ�ɱ������������ˡ�
					next_param[i][0] = swap;
				}
				else {
					/*���ӽ���Ĵ�������20��С��2��ʱ�򣬽������´���*/
					dealWithClusterNumNotCorrect(p + i * max_cluster_num * D, p2 + i * 2, next_index[i], next_param[i], D);
				}
				
			}
			//subspace_cross(cr_lenth, in_area[i], out_area[i], (p + i * MAX_CLUSTER * D), (p2 + i * 2), next_index[i], next_param[i], D);

			//Selection operator,���������䵽����Ĵ���
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

			//step 2 of the selection: calculate the fitness��next_param[i][0]��¼�˵�ǰ���������Ĵ�����data���������ݼ���next_index[i]�����˵�ǰ������������ĸ��������ģ�in_cluster2[i]���ݼ���ÿ�����������ĸ��ء�
			next_param[i][1] = func(next_param[i][0], *data, next_index[i], in_cluster2[i], N, D,usedIndex);
			numofE = numofE + 1;

	
			//Iָ����Ӧ��Խ��Խ��
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
					index = i;//��¼���Ÿ����λ��
				}
			}
			//DBָ����Ӧ��ԽСԽ��
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
					best_val = *(p2 + i * 2 + 1);//���ŵ���Ӧ��
					best_val2 = *(p2 + i * 2 + 0);//���ŵĴ���
					index = i;//��¼���Ÿ����λ��
				}

			}
			
		}
	/*	printf("\n---------��%d��-������Ӧ��Ϊ%g-��Ӧ�Ĵ���Ϊ%d------------------\n", k+1, best_val, (int)best_val2);
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
	int* originDataClusterInfo;//�����������ݼ�����ľ������
	data = loadData2(filename, N, D);   //load data from the text file,�����ݼ��浽�˶�ά����data��
	originDataClusterInfo=getDataOriginClusterInfo(filename, N);
	set_traildataInfo(traildata, originDataClusterInfo[N], filename, usedIndex);
	printf("�ļ�%s��%d�����ݣ�ά����%d\n", filename, N, D);
	max_cluster_num = getMaxClusterNum(N, max_cluster_num);//ÿ�ִ����ľ��෽����Ҫ���㣬ÿ��������������������������Ĵ���20�����˸�Ҫ��2-19��Щ�������ܶ������Ҫ��
	NP = getNP(D, 10);//ά��Խ������Ľ�ռ��Խ�������Ⱥ�ĸ�����Ӧ�ø���
	Gmax = getGmax(NP, 1000000);//ÿһ������Ҫ����NP����������˸�ʽ���ܵõ����ĵ�������
	printf("The times of iteration(Gmax):%d\n", Gmax);
	malloc1D(uk, D);
	malloc1D(lk, D);
	malloc1E(popul_rand, N);
	malloc2D(in_cluster, NP, N);        //population index����ά����in_cluster������¼��ǰ�������ĸ��ص����
	malloc2E(popul_index, NP,max_cluster_num * D);  //population cluster centroids
	malloc2E(popul_param, NP, 2);       //population info(cluster number, fitness value)
	initial_two_Dim_intarr(in_cluster, NP, N, 0);
	getDataDim_max_min(data, uk, lk, D, N);//��ȡ����ÿ��ά�ȵ����ֵ��Сֵ�����浽����uk,lk
	//population initialization
	for (i = 0; i < NP; i++) {
		initial_individual(data, popul_index[i], popul_param[i], in_cluster[i], N, D, max_cluster_num, min_cluster_num);
		popul_param[i][1] = func(popul_param[i][0], *data, popul_index[i], in_cluster[i], N, D,usedIndex);//���㵱ǰ�������Ӧ��	
	}
	best = e_de(*popul_index, *popul_param, N, D, Gmax, lk, uk,usedIndex,in_cluster); //run the E-DE algorithm

	//���㱾�ξ���Ĵ��ھ��룬�ؼ���룬������Ӧ��ֵ����׼�Rָ��
	popul_index[best]= count_meanV(*data, in_cluster[best], popul_param[best][0], N, D);//���¼����������Ĵ�����
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
		 printf("��ǰ�������ļ�Ϊ:%s\n", filename);
		 ede_initialization(filename,usedIndex);
		
	}
	free(outputfile);


}
void ede_two_index(const char* filename) {
	for (int i = 0; i < 2; i++) {
		char* indexname = getIndexName(i);
		freopen("CON", "w", stdout);
		printf("����ʹ��ָ��%s���д���...\n", indexname);
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
	printf("��������Ҫ������ļ���ţ�\n");
	scanf("%d", &i);
	while (i < 1 || i>4) {
		printf("������������������(1~4)\n");
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
	printf("��ǰ������ļ���%s\n����ִ����,����رյ�ǰ����....\n", filename);
	return filename;
}


int  main()
{

	////����������ָ����д���
	//srand((unsigned int)(time(NULL)));
	char* filename = showmenu();
	ede_two_index(filename);
	free(filename);
	freopen("CON", "w", stdout);
	printf("����ִ����ϣ�");
	fclose(stdout);
	system("pause");
	return 0;
}



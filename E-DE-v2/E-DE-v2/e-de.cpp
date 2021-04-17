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
#include "boost/random.hpp"
#include"ede.h"
#include"testtool.h"




int min_cluster_num = MIN_CLUSTER;
int max_cluster_num = MAX_CLUSTER;

FILE* fp;//�����������
double** distarr;//���ڼ�������ϵ��
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
	}else if(usedIndex==DB_INDEX) {
		fes = DB_index((int)k_num, p_pars, p_index, in_clu, N, D);
	}else if (usedIndex == SIHOUETTES_INDEX) {
		fes = getSilhouettes(in_clu, N, k_num, D, distarr);
	}else {
		printf("�����ָ������\n");
		exit(1);
	}
	

	return fes;
}

/*****************************************************************
* func_name: the gaussian generator
* input: mean and standard deviation
* descript: for the crossover operation
*****************************************************************/
double sampleNormal(double mean, double sigma)
{
	 // apply the unix time to set the seed of rand
	 static boost::mt19937 rng(static_cast< unsigned >(std::time(0)));

	// select the gaussian random distribution
	 boost::normal_distribution< double > norm_dist(mean, sigma);

	 // generate the random generator
	 boost::variate_generator< boost::mt19937&, boost::normal_distribution< double > > normal_sampler(rng, norm_dist);
	
	 return normal_sampler();
	
	//double gassrandN = gaussrand(0.0, 0.1);
	//return gassrandN;
}




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
	double F = 0.5, CR = 0.5, MU, gauss_index;
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
	outputPopulationInfo2(p, p2, D, NP, max_cluster_num, in_cluster, N, 1, best_val, best_val2, 20000,fp);
	//һ������Gmax��
	Gmax = 2000;
	for (k = 0; k < Gmax; k++)      //Gmax denotes the maximum iteration times
	{
		printf("���ڴ����%d����Ⱥ...\n", k + 1);
		
		for (i = 0; i < NP; i++)
		{
			////Mutation operator
			////step 1 of the mutation: random select 3 individuals
			//do {
			//	r1 = (int)(NP * URAND); //�ڱ�������У���������ǻ��ڵ�r1����������ģ�r2,r3�����ṩ��ֵ���Ϣ��Ȼ����r1��������Ӧ�ı䣬���õ��˵�ǰĿ�����ı������
			//} while (r1 == i);
			//do {
			//	r2 = (int)(NP * URAND);
			//} while (r2 == i || r2 == r1);
			//do {
			//	r3 = (int)(NP * URAND);
			//} while (r3 == i || r3 == r1 || r3 == r2);
			getRandomNum(r1, r2, r3, i, NP);
			//step 2 of the mutation: determine the cluster number of mutant vector
			//����Ӧ���ǻ�ȡ�����������Ĵ�������������p2Ӧ��Ҳ�Ƕ�ά���飬һ����np�к�2�У�����������Ⱥ��ÿ������Ĵ�����
			int mutate_d1 = (int)*(p2 + r1 * 2);  //��ȡĿ�����r1�Ĵ�������
			int mutate_d2 = (int)*(p2 + r2 * 2);	//��ȡĿ�����r2�Ĵ�������
			int mutate_d3 = (int)*(p2 + r3 * 2);	//��ȡĿ�����r3�Ĵ�������

			int mutate_u = mutate_d1 + (int)(F * ((double)mutate_d2 - (double)mutate_d3));

			//if (mutate_u > max_cluster_num) //����������һ������Ĵ���������ȡֵ��Χ��2~20
			//	mutate_u = max_cluster_num;
			//if (mutate_u < min_cluster_num)
			//	mutate_u = min_cluster_num;
			mutate_u = cross_border_process(mutate_u, max_cluster_num, min_cluster_num);
			next_param[i][0] = mutate_u; //������������NP�����壬�������Ǳ������������ÿ������Ĵ�������

			for (int j = 0; j < max_cluster_num; j++) { //������������������ı��촦��
				r1_rand[j] = 0;
				r2_rand[j] = 0;
				r3_rand[j] = 0;
			}

			//step 3 of the mutation: produce the initial mutant vector
			//pӦ���Ǳ�������Ŀ���������飬����pҲ�Ƕ�ά����p���׵�ַ��r1��ָ��r1�����壬
			//ÿ�����嶼Ԥ����20���صĴ�С�����20*D��ָÿ������ĳ��ȣ�r1*20*D���ǵ�r1��������׵�ַ
			if (mutate_u == mutate_d1) {      //�������Ĵ���������Ŀ�����r1�Ĵ���������ȵ����
				for (int j = 0; j < mutate_d1 * D; j++) //r1Ŀ���������д�����ֱ�Ӹ��Ƶ����������
					next_index[i][j] = *(p + r1 * 20 * D + j);
			}
			else if (mutate_u > mutate_d1) {//�������Ĵ�������Ŀ�����������������ʱ�򣬸���r2�Ĵ��������Ǳȸ���r3�Ĵ�������Ҫ��������ߵĲ�ֵ���ܱ�r3�Ĵ�������Ҫ�󣬵�Ȼ��ֵ�϶��Ǳ�r2�Ĵ�������ҪС��
				for (int j = 0; j < mutate_d1 * D; j++)//�Ȱ�Ŀ�����r1�����д����ĸ��Ƶ��������
					next_index[i][j] = *(p + r1 * 20 * D + j);

				//select from the r2 or r3 individual
				//ʣ�µĴ����ĴӸ���r2��r3�����ѡȡ
				r3_full = 0;
				for (int j = mutate_d1; j < mutate_u; j++) { //��mutate_d1��mutate_u���Ǳ�������������Ŀ�����r1�����Ĳ�ֵ
					if (URAND < 0.5) {   //select from the r2
						r4 = (int)(mutate_d2 * URAND);
						while (r2_rand[r4] == 1) {
							r4 = (int)(mutate_d2 * URAND);
						}
						for (int l = 0; l < D; l++) {
							next_index[i][j * D + l] = *(p + r2 * 20 * D + r4 * D + l);
						}
						r2_rand[r4] = 1;
					}
					else {   //select from the r3��r3�Ĵ����������ܻ�Ȳ�ֵС�������������r3_full����Ǵ�r3���ѡ�еĴ��������������ֵ������r3�Ĵ�������mutate_d3���ʾr3�е����и����Ѿ���ѡ�ϡ�
						if (r3_full < mutate_d3) {
							r4 = (int)(mutate_d3 * URAND);
							while (r3_rand[r4] == 1) {
								r4 = (int)(mutate_d3 * URAND);
							}
							for (int l = 0; l < D; l++) {
								next_index[i][j * D + l] = *(p + r3 * 20 * D + r4 * D + l);
							}
							r3_rand[r4] = 1;
							r3_full = r3_full + 1; //��r3ѡ���˸��壬r3_full�ͼ�1.
						}
						else { //���r3�Ĵ��������Ѿ���ѡ���˻�����������Ĳ�ֵ����ʣ�µĴ����ľʹӸ���r2��ѡ��
							r4 = (int)(mutate_d2 * URAND);
							while (r2_rand[r4] == 1) {
								r4 = (int)(mutate_d2 * URAND);
							}
							for (int l = 0; l < D; l++) {
								next_index[i][j * D + l] = *(p + r2 * 20 * D + r4 * D + l);
							}
							r2_rand[r4] = 1;
						}
					}
				}
			}
			else {  //delete from the r1 individual.�������Ĵ���������Ŀ������������ҪС�����
				int mutate_num = mutate_d1 - mutate_u;//�������У�mutate_d1�϶��Ǵ���mutate_u�ģ���mutate_u����2~20��Χ�ڣ���mutate_u�Ǵ���0�ģ����mutate_num�϶���С��mutate_d1�ģ�Ҳ�����ڸ���r1��ɾ��mutate_num�������ĺ󣬸���r1��Ȼ�ᱣ������2�������ġ�
				while (mutate_num != 0) {
					r4 = (int)(mutate_d1 * URAND);
					while (r1_rand[r4] == 1) { //����r1_rand��ֵΪ1������Ӧ�Ĵؾ�����Ҫ��Ŀ������r1ɾ���Ĵ�
						r4 = (int)(mutate_d1 * URAND);
					}
					r1_rand[r4] = 1;
					mutate_num = mutate_num - 1;
				}

				r1_full = 0;
				for (int j = 0; j < mutate_d1; j++) { //j��ʾ����r1�ĵ�j����
					if (r1_rand[j] == 1)
						continue;
					//next_index�Ǳ����ʼ�����������飬���������Ĵ�����Ŀ���������ٵ�ʱ�򣬱����������д�������Ŀ��������ɾ��mutate_num�������ĺ����µĴ�����
					for (int l = 0; l < D; l++) {//ÿ�θ���һ����j,r1_full�������Ʊ�������еĴ�����λ�ã�Ŀ�����øø��������еĴ����Ķ������ı�����һ��
						next_index[i][r1_full * D + l] = *(p + r1 * 20 * D + j * D + l);//j���Ǳ�ʾĿ���������Ҫ�����Ƶ���������д�������ʼ�±꣬����mutate_num�������ģ�����Ĵ����Ķ���Ҫ������
					}
					r1_full = r1_full + 1;
				}
			}

			//step 4 of the mutation: fine tune the cluster centroids
			MU = 0.10 - 0.06 * k / ((double)Gmax - 1.0); //����Ĳ����ǵõ���ʼ�ı�����壬�ò����ǶԳ�ʼ������尴һ�����ʽ��д���ÿ�ֵ�����Ӹ�˹�Ŷ�����ı�������������NP*MU
			 //MU = 0.1 + 0.3 * k / ((double)Gmax - 1.0);
			if (URAND < MU) {//ÿ�ζ��и���MU�Գ�ʼ���������Ӹ�˹�Ŷ�
				for (int j = 0; j < mutate_u; j++) { //��ǰ��ʼ��������ÿһ�������ģ���ÿһά������Ӹ�˹�Ŷ�
					for (int l = 0; l < D; l++) {
						gauss_index = next_index[i][j * D + l] + sampleNormal(0, 0.1) * (Xu[l] - Xl[l]);//sampleNormal�����ɸ�˹�������������0�Ǿ�ֵ��0.1�Ǳ�׼���˹�ֲ��������ʾ���ʣ��������ʾ��ֵ�����ɸú���ͼ����Եó���������������0��0.1���������������
						if (gauss_index > Xu[l]) {
							next_index[i][j * D + l] = Xu[l]; //Խ���ȡ�߽�ֵ
						}
						else if (gauss_index < Xl[l]) {
							next_index[i][j * D + l] = Xl[l];
						}
						else {
							next_index[i][j * D + l] = gauss_index;
						}
					}
				}
			}
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
			int rand1;
			swap = 0, swap1 = 0;    //the number of points that outside and inside the swap area
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
			/*�����ӿռ佻��*/
			if (cr_lenth == 1) {//��Ҫ�ӽ��Ĵ�ֻ��һ���������
				for (int j = 0; j < 20 * D; j++)
					out_area[i][j] = *(p + i * 20 * D + j); //�Ȱѵ�ǰĿ���������дر��浽����out_area�С�
				//�����λ��Ŀ������������������20���أ���ô�ӿռ佻��ʱ��ֻ��Ҫ����Ҫ�ӽ��Ĵ��滻��λĿ����������һ���ء�  
				if ((int)*(p2 + i * 2) == 20) {
					rand1 = (int)(((int)*(p2 + i * 2)) * URAND);
					for (int j = 0; j < D; j++)
						out_area[i][rand1 * D + j] = in_area[i][j]; //���cr_lenth==1,������Ĵ�ֻ��һ������ô����in_area[i]ֻ������һ���أ�Ȼ���������滻��out_area��������һ��λ��
					for (int j = 0; j < 20 * D; j++)
						next_index[i][j] = out_area[i][j]; //��out_area�����еĴظ��Ƶ��������б����������飬���õ��������������顣
					next_param[i][0] = *(p2 + i * 2); //next_param���鱣����ǵ�i���������Ĵ�����������Ϊ�ӽ��Ĵ�ֻ��һ����˵�i��������壨�����������Ĵ������λ��Ŀ�����Ĵ�������һ��
				}
				//��λ��Ŀ�����Ĵ�����С��20������û���أ���ô�ӿռ��ӽ���ֻ��Ҫ�Ѵ��ӽ��Ĵز��뵽Ŀ������ĩβ��
				else {
					for (int j = 0; j < D; j++)
						out_area[i][((int)*(p2 + i * 2)) * D + j] = in_area[i][j];
					for (int j = 0; j < 20 * D; j++)
						next_index[i][j] = out_area[i][j];
					next_param[i][0] = *(p2 + i * 2) + 1;
				}
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
				int swap_flag, swap_num=0;

				////�ӱ��������л�ȡ�ӿռ��ڵĴ�
				//for (int j = 0; j < (int)next_param[i][0]; j++) {
				//	swap_flag = 1;
				//	for (int l = 0; l < D; l++) {
				//		//������������˵����ǰ�Ĵز����ӽ��ӿռ䷶Χ��
				//		if (fabs(next_index[i][j*D+l] - center[l]) > (dist[l]+0.00001)) {   //the center can be changed
				//			swap_flag = 0;
				//			break;
				//		}
				//	}
				//	if (swap_flag == 1) {//�Ѳ���Ҫ�����ǵĴر������������ⲿ�ִ���Ŀ������в���Ҫ����������滻�Ĵء�����ⲿ�ִ����������б�ѡ�е��ӽ��ؽ�ϳ���������
				//		for (int l = 0; l < D; l++)
				//			in_area[i][swap_num * D + l] = next_index[i][j * D + l];
				//		swap_num = swap_num + 1;
				//	}
				//}
				//cr_lenth = swap_num;

				for (int j = 0; j < (int)*(p2 + i * 2); j++) {
					swap_if = 0;
					for (int l = 0; l < D; l++) {
						//������������˵����ǰ�Ĵز����ӽ��ӿռ䷶Χ��
						if (fabs(*(p + i * 20 * D + j * D + l) - center[l]) > (dist[l] + 0.00001)) {   //the center can be changed
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
				}else {
					/*���ӽ���Ĵ�������20��С��2��ʱ�򣬽������´���*/
					dealWithClusterNumNotCorrect(p + i * max_cluster_num * D, p2 + i * 2, next_index[i], next_param[i], D);
				}

				

			}

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

	
			//Iָ�������ϵ��Խ��Խ��
			if (usedIndex == I_INDEX||usedIndex==SIHOUETTES_INDEX) {
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
		outputPopulationInfo2(p, p2, D, NP, max_cluster_num, in_cluster, N,k+1,best_val,best_val2,20000,fp);
		
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
	srand((unsigned int)(time(NULL)));
	int i, j, D, N, Gmax, NP, best = 0, * popul_rand, ** in_cluster;
	double data_min, data_max, min, distance, * uk, * lk, ** popul_index, ** popul_param;
	int* originDataClusterInfo;//�����������ݼ�����ľ������
	data = loadData2(filename, N, D);   //load data from the text file,�����ݼ��浽�˶�ά����data��
	distarr = getDataDistance( data,  N,D);
	originDataClusterInfo=getDataOriginClusterInfo(filename, N);
	set_traildataInfo(traildata, originDataClusterInfo[N], filename, usedIndex);
	fprintf(fp,"�ļ�%s��%d�����ݣ�ά����%d\n", filename, N, D);
	max_cluster_num = getMaxClusterNum(N, max_cluster_num);//ÿ�ִ����ľ��෽����Ҫ���㣬ÿ��������������������������Ĵ���20�����˸�Ҫ��2-19��Щ�������ܶ������Ҫ��
	NP = getNP(D, 10);//ά��Խ������Ľ�ռ��Խ�������Ⱥ�ĸ�����Ӧ�ø���
	Gmax = getGmax(NP, 1000000);//ÿһ������Ҫ����NP����������˸�ʽ���ܵõ����ĵ�������
	fprintf(fp,"The times of iteration(Gmax):%d\n", Gmax);
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
	
	//getSilhouettes(in_cluster[0], N, (int)popul_param[0][0], D, distarr);


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
	freeTwoDimArr_double(distarr, N);

}






void test_ede(int usedIndex,const char* filename){
	char* outputfile;
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		printf("#####���ڽ��е�%d�β���..\n",i);
		outputfile = getOutputFilename(filename, usedIndex,i+1);	
		 if ((fp = fopen(outputfile, "w")) == NULL) {
			 printf("�ļ����ܴ�\n");
			 exit(0);
		 }
		 //freopen((const char* )outputfile, "w", stdout);
		 fprintf(fp,"��ǰ�������ļ�Ϊ:%s\n", filename);
		 ede_initialization(filename,usedIndex);
		 fclose(fp);
		 free(outputfile);
	}


}
void ede_two_index(const char* filename) {
	for (int i = 2; i < 3; i++) {
		char* indexname = getIndexName(i);
		//freopen("CON", "w", stdout);
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
	//freopen("CON", "w", stdout);
	printf("����ִ����ϣ�");
	//fclose(stdout);
	free(filename);
	system("pause");
	return 0;
}



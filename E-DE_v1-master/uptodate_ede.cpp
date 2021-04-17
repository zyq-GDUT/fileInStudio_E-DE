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
 //#include "boost/random.hpp"
#include"gaussrand.h"
// change any of these parameters to match your needs 

#define URAND  ((double)rand()/((double)RAND_MAX+1.0)) //产生的随机数是[0,1)


// declaration of functions used by this algorithm 

// 1. load and arrange the data 
double **data;
double **loadData(int *d, int *n);
void malloc1D(double *&a, int D);
void malloc1E(int *&a, int D);
void malloc2D(int **&a, int xDim, int yDim);
void malloc2E(double **&a, int xDim, int yDim);

// 2. functions to implement the E-DE algorithm 
double func(double k_num, double *p_pars, double *p_index, int *in_clu, int N, int D);
double getDistance(double *avector, double *bvector, int n);
int e_de(double *p, double *p2, int N, int D, int Gmax, double *Xl, double *Xu);

/********************************************************************
* func_name: fitness calculation
* input: the cluster number, population info, data info
* descript: calculating the fitness of an individual    
* this takes a user defined function                            
********************************************************************/
/**
 * 计算当前个体每个簇的均值向量和每个簇对应的样本数
 * k_num是指总簇数
 * **/
double* count_meanV(double* data,int *in_clu,int k_num,int N,int D){
	double *meanV;//用来保存每个簇的中心坐标（均值向量）
	int *numInclu;//用来保存每个簇对应的样本数
	malloc1D(meanV,k_num*D);
	malloc1E(numInclu,k_num);
	memset(meanV,0,sizeof(*meanV)*k_num*D);
	memset(numInclu,0,sizeof(*numInclu)*k_num);
	for(int i=0;i<N;i++){
		int cn=in_clu[i];
		numInclu[cn]++;
		for(int j=0;j<D;j++){
			*(meanV+cn*D+j)+=*(data+i*D+j);
		}
	}

	//求均值向量
	for(int k=0;k<k_num;k++){
		for(int i=0;i<D;i++){
			*(meanV+k*D+i)=*(meanV+k*D+i)/numInclu[k];
		}
	}


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
double I_index(int k_num, double *data, double *p_index, int *in_clu, int N, int D){
	
	p_index=count_meanV(data,in_clu,k_num,N,D);
	int i,j,k;
	double E=0;
		//计算每个簇与簇中心的距离之和，即计算E值
	for(i=0;i<N;i++){
		for( j=0;j<D;j++){
			int cnum=in_clu[i];
			E+=fabs(*(data+i*D+j)-*(p_index+cnum*D+j));
		}
	}
	// printf("E值=%lf\n",E);
	//找出所有簇中距离最远的两个簇的距离。即计算Dmax值
	double Dmax=0.0;
	double distance;
	for(i=0;i<k_num;i++){
		for(j=i+1;j<k_num;j++){
				distance=0.0;
				for(k=0;k<D;k++){
					double d=*(p_index+i*D+k)-*(p_index+j*D+k);
					distance+=fabs(d);
				}
				if(distance>Dmax){
					Dmax=distance;
				}
		}
	}

	// printf("k_num=%d\tN=%d\tE=%lf\tD=%lf\n",k_num,N,E,Dmax);
	//计算I值
	double i_index=pow((1.0/(double)k_num)*((double)N/E)*Dmax,N);
	printf("i值=%g\n",i_index);
	printf("\n------\n");
	return i_index;

}

//DB指标函数
double DB_index(double k_num, double *p_pars, double *p_index, int *in_clu, int N, int D){


}
double func(double k_num, double *p_pars, double *p_index, int *in_clu, int N, int D) {
	double fes;  
	/* define your testing function here */
	
	 fes=I_index((int)k_num,p_pars,p_index,in_clu,N,D);
	return fes;
}

/*****************************************************************
* func_name: the gaussian generator
* input: mean and standard deviation
* descript: for the crossover operation
*****************************************************************/
double sampleNormal(double mean, double sigma)
{
	// // apply the unix time to set the seed of rand
	// static boost::mt19937 rng(static_cast< unsigned >(std::time(0)));

	// // select the gaussian random distribution
	// boost::normal_distribution< double > norm_dist(mean, sigma);

	// // generate the random generator
	// boost::variate_generator< boost::mt19937&, boost::normal_distribution< double > > normal_sampler(rng, norm_dist);

	// return normal_sampler();
	double gassrandN=gaussrand(0,0.1);
	return gassrandN;
}


/*****************************************************************
* func_name: load the testing data as you need 
* input: dimension and data number
* descript: load data from the txt file
*d是指数据的维数，n是指数据样本的总数
*****************************************************************/
/*
统计文件中的样本数和维度数
filename是指要统计的文件名
n用来返回文件中的行数，即样本数
d用来返回样本的维度数
*/
void count_N_D(const char *filename,int &n,int &d){
	FILE *fp;
	char buffer[1000];
	char c;
	int bufferlen;
	int i;
	n=0;
	d=0;
 
	if((fp=fopen(filename,"rb"))==NULL){
		printf("文件不能打开\n");
		exit(0);
	}
	while(fgets(buffer,1000,fp)!=NULL){
		bufferlen=strlen(buffer);
		//跳过空白行
		if(bufferlen==2&&buffer[0]==13&&buffer[1]==10){
			continue;
		}
		
		n++;
		//在文件的第二行统计当前数据集的维数
		if(n==2){
			for(i=0;i<bufferlen;i++){
				c=buffer[i];
				if(c==','){
					d++;
				}
			}
		}
		
	}
	fclose(fp);
	printf("行数：%d\n",n);
	printf("维数：%d\n",d);
 
	
}

double **loadData(const char* filename,int &n,int &d)
{
	int i, j;
	double **arraydata; 
	char buffer[100];
 	FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("file open error\n");
        exit(1);

    }
	 count_N_D(filename,n,d);

    malloc2E(arraydata,n,d);

    for (int i = 0; i < n; i++)
    {
		for(int j=0;j<d;j++){
			
			fscanf(fp,"%lf,",&arraydata[i][j]);
			// printf("%lf\t",arraydata[i][j]);
			
			if(j==3){
			// printf("\n");
			fgets(buffer,100,fp);
			}
		}
        // fscanf(fp, "%lf,%lf,%lf,%lf", &a[i][0], &a[i][1], &a[i][2],&a[i][3]);
        // printf("%lf,%lf,%lf,%lf\n", a[i][0],a[i][1],a[i][2],a[i][3]);
    }
    fclose(fp);
	return arraydata;
}

/*****************************************************************
* func_name: array function
* input: array name and dimension
* descript: allocate the int or double array
*****************************************************************/
void malloc1D(double *&a, int D) {
	a = (double *)malloc(D * sizeof(double));
	if (a == NULL)
		perror("malloc");
}

void malloc1E(int *&a, int D) {
	a = (int *)malloc(D 
	* sizeof(int));
	if (a == NULL)
		perror("malloc");
}

void malloc2D(int **&a, int xDim, int yDim)
{
	a = (int **)malloc(xDim * sizeof(int *));//a是二维数组，给二维数组a分配xDim个能存整型指针的内存
	a[0] = (int *)malloc(xDim * yDim * sizeof(int));//先给二维指针的第一个元素分配初始地址，该地址是指向整个二维数组的内存
	for (int i = 1; i<xDim; i++)
	{
		a[i] = a[i - 1] + yDim;//给剩下的元素赋值
	}
	if (a == NULL)
		perror("malloc");
}

void malloc2E(double **&a, int xDim, int yDim)
{
	a = (double **)malloc(xDim * sizeof(double *));
	a[0] = (double *)malloc(xDim * yDim * sizeof(double));//a[0]表示指向二维数组a的首地址
	for (int i = 1; i<xDim; i++)
	{
		a[i] = a[i - 1] + yDim;//a是二维数组，a[i]表示二维数组中每个元素的首地址，二维数组的每个元素都是一维数组。
	}
	if (a == NULL)
		perror("malloc");
}

/*****************************************************************
* func_name: distance calculation
* input: the two testing vertors and the dimension length
* descript: for calculating the Euclidean distance 欧几里得距离
*****************************************************************/
double getDistance(double *avector, double *bvector, int n)
{
	int i;
	double sum = 0;
	for (i = 0; i<n; i++)
		sum = sum + pow(*(avector + i) - *(bvector + i), 2); //两个向量的对应维度的平方和再开方即是欧几里得距离

	return sqrt(sum);
}

/*****************************************************************
* func_name: E-DE algorithm
* input: the evolution population, the data number N, dimension D
*        the lower and upper bound Xl, Xu
* descript: do the mutation and crossover opertors, and select the
*           best results for the next generation
*****************************************************************/
int e_de(double *p, double *p2, int N, int D, int Gmax, double *Xl, double *Xu) {
	int i, k, r1, r2, r3, r4, r1_full, r3_full;
	int *r1_rand, *r2_rand, *r3_rand, **in_cluster2;
	int NP = 10 * D, numofE = 0, index = 0, swap, swap1, swap_if, cr_lenth, cr_n;

	double F = 0.5, CR = 0.4, MU, gauss_index;
	double **next_index, **next_param, distance;
	double **out_area, **in_area, *dist, *center, best_val = DBL_MIN, best_val2 = 0, min;

	malloc1D(dist, D);
	malloc1D(center, D);
	malloc1E(r1_rand, 20);  //最大20个簇，这三个数组与簇数相关。因此这三个随机数组大小设置为20. 
	malloc1E(r2_rand, 20);
	malloc1E(r3_rand, 20);

	malloc2E(out_area, NP, 20 * D);   
	malloc2E(in_area, NP, 20 * D);     
	malloc2E(next_index, NP, 20 * D);//NP行20*D列，用来保存每个变异个体的所有簇中心
	malloc2E(next_param, NP, 2); //NP行2列，用来保存每个变异个体的簇中心数
	malloc2D(in_cluster2, NP, N);//用来标记每个个体中，样本所属的簇，每个个体都是一种聚类方案，而个体只是保存了一种聚类方案的簇中心。因此需要用一个二维数组来表明每种聚类方案的样本归属情况

	for (i = 0; i<NP; i++) {
		for (int j = 0; j<2; j++) {
			next_param[i][j] = 0;
		}
		for (int j = 0; j<N; j++) {
			in_cluster2[i][j] = 0;
		}
		for (int j = 0; j<20 * D; j++) {
			out_area[i][j] = 0;
			in_area[i][j] = 0;
			next_index[i][j] = 0; //用来保存变异个体
		}
	}

	for (i = 0; i<D; i++) {
		dist[i] = 0; //距离数组，用来记录杂交操作时，杂交子空间中，每一维的取值范围
		center[i] = 0; //用来记录子空间中的簇中心
	}
	//一共迭代Gmax次
	for (k = 0; k<Gmax; k++)      //Gmax denotes the maximum iteration times
	{
		for (i = 0; i<NP; i++)   
		{
			//Mutation operator
			//step 1 of the mutation: random select 3 individuals
			do{
				r1 = (int)(NP*URAND); //在变异操作中，变异个体是基于第r1个个体产生的，r2,r3个体提供差分的信息，然后在r1上做出相应改变，即得到了当前目标个体的变异个体
			} while (r1 == i);
			do{
				r2 = (int)(NP*URAND);
			} while (r2 == i || r2 == r1);
			do{
				r3 = (int)(NP*URAND);
			} while (r3 == i || r3 == r1 || r3 == r2);

			//step 2 of the mutation: determine the cluster number of mutant vector
			//这里应该是获取随机三个个体的簇中心数，这里p2应该也是二维数组，一共有np行和2列，用来保存种群中每个个体的簇数。
			int mutate_d1 = (int) *(p2 + r1 * 2);  //获取目标个体r1的簇中心数
			int mutate_d2= (int) *(p2 + r2 * 2);	//获取目标个体r2的簇中心数
			int mutate_d3 =  (int) *(p2 + r3 * 2);	//获取目标个体r3的簇中心数

			int mutate_u = mutate_d1 + (int)( F *(mutate_d2 - mutate_d3));

			if (mutate_u > 20) //变异向量中一个个体的簇中心数的取值范围是2~20
				mutate_u = 20;
			if (mutate_u < 2)
				mutate_u = 2;

			next_param[i][0] = mutate_u; //变异向量中有NP个个体，该数组是保存变异向量中每个个体的簇中心数

			for (int j = 0; j<20; j++) { //这三个数组用于下面的变异处理
				r1_rand[j] = 0; 
				r2_rand[j] = 0;
				r3_rand[j] = 0;
			}
			
			//step 3 of the mutation: produce the initial mutant vector
			//p应该是保存所有目标个体的数组，这里p也是二维数组p的首地址，r1是指第r1个个体，
			//每个个体都预留了20个簇的大小，因此20*D是指每个个体的长度，r1*20*D就是第r1个个体的首地址
			if (mutate_u == mutate_d1) {      //变异个体的簇中心数与目标个体r1的簇中心数相等的情况
				for (int j = 0; j<mutate_d1*D; j++) //r1目标个体的所有簇中心直接复制到变异个体中
					next_index[i][j] = *(p + r1 * 20 * D + j); 
			}
			else if (mutate_u > mutate_d1) {//变异个体的簇数大于目标个体簇数的情况，这时候，个体r2的簇中心数是比个体r3的簇中心数要大。因此两者的差值可能比r3的簇中心数要大，当然差值肯定是比r2的簇中心数要小的
				for (int j = 0; j<mutate_d1*D; j++)//先把目标个体r1的所有簇中心复制到变异个体
					next_index[i][j] = *(p + r1 * 20 * D + j);

				//select from the r2 or r3 individual
				//剩下的簇中心从个体r2或r3中随机选取
				r3_full = 0;
				for (int j = mutate_d1; j<mutate_u; j++) { //从mutate_d1到mutate_u就是变异向量簇数与目标个体r1簇数的差值
					if (URAND < 0.5) {   //select from the r2
						r4 = (int)(mutate_d2*URAND);   
						while (r2_rand[r4] == 1) {   
							r4 = (int)(mutate_d2*URAND);   
						}
						for (int l = 0; l<D; l++) {
							next_index[i][j*D + l] = *(p + r2 * 20 * D + r4 * D + l);
						}                 
						r2_rand[r4] = 1;
					}
					else {   //select from the r3，r3的簇中心数可能会比差值小，因此这里用了r3_full来标记从r3随机选中的簇中心数，如果该值等于了r3的簇中心数mutate_d3则表示r3中的所有个体已经被选上。
						if (r3_full < mutate_d3) {
							r4 = (int)(mutate_d3*URAND);   
							while (r3_rand[r4] == 1) {   
								r4 = (int)(mutate_d3*URAND);   
							}
							for (int l = 0; l<D; l++) {
								next_index[i][j*D + l] = *(p + r3 * 20 * D + r4 * D + l);
							}
							r3_rand[r4] = 1;
							r3_full = r3_full + 1; //从r3选择了个体，r3_full就加1.
						}
						else { //如果r3的簇中心数已经被选完了还不满足所需的差值，则剩下的簇中心就从个体r2中选择。
							r4 = (int)(mutate_d2*URAND);   
							while (r2_rand[r4] == 1) {   
								r4 = (int)(mutate_d2*URAND);  
							}
							for (int l = 0; l<D; l++) {
								next_index[i][j*D + l] = *(p + r2 * 20 * D + r4 * D + l);
							}
							r2_rand[r4] = 1;
						}
					}
				}
			}
			else {  //delete from the r1 individual.变异个体的簇中心数比目标个体簇中心数要小的情况
				int mutate_num = mutate_d1 - mutate_u;//该条件中，mutate_d1肯定是大于mutate_u的，而mutate_u是在2~20范围内，即mutate_u是大于0的，因此mutate_num肯定是小于mutate_d1的，也就是在个体r1中删除mutate_num个簇中心后，个体r1仍然会保留至少2个簇中心。
				while (mutate_num != 0) {
					r4 = (int)(mutate_d1*URAND);
					while (r1_rand[r4] == 1) { //数组r1_rand中值为1的所对应的簇就是需要在目标向量r1删除的簇
						r4 = (int)(mutate_d1*URAND);
					}
					r1_rand[r4] = 1;
					mutate_num = mutate_num - 1;
				}

				r1_full = 0;
				for (int j = 0; j<mutate_d1; j++) { //j表示个体r1的第j个簇
					if (r1_rand[j] == 1)
						continue;
					//next_index是保存初始变异个体的数组，当变异个体的簇数比目标个体簇数少的时候，变异个体的所有簇中心是目标个体随机删除mutate_num个簇中心后留下的簇中心
					for (int l = 0; l<D; l++) {//每次复制一个簇j,r1_full用来控制变异个体中的簇中心位置，目的是让该个体中所有的簇中心都连续的保存在一起
						next_index[i][r1_full*D + l] = *(p + r1 * 20 * D + j * D + l);//j则是表示目标个体中需要被复制到变异个体中簇中心起始下标，除了mutate_num个簇中心，其余的簇中心都需要被复制
					}
					r1_full = r1_full + 1;
				}
			}

			//step 4 of the mutation: fine tune the cluster centroids
			MU = 0.1 - 0.06 * k / (Gmax - 1); //上面的操作是得到初始的变异个体，该步骤是对初始变异个体按一定概率进行处理，每轮迭代添加高斯扰动处理的变异个体数大概是NP*MU
			if (URAND < MU) {//每次都有概率MU对初始变异个体添加高斯扰动
				for (int j = 0; j<mutate_u; j++) { //当前初始变异个体的每一个簇中心，即每一维都会添加高斯扰动
					for (int l = 0; l<D; l++) { 
						gauss_index = next_index[i][j*D + l] + sampleNormal(0,0.1)*(Xu[l] - Xl[l]);//sampleNormal是生成高斯随机数生成器，0是均值，0.1是标准差。高斯分布纵坐标表示概率，横坐标表示数值，而由该函数图像可以得出，该生成器生成0±0.1附近的数概率最大。
						if (gauss_index > Xu[l]) {
							next_index[i][j*D + l] = Xu[l]; //越界就取边界值
						}
						else if (gauss_index < Xl[l]) {
							next_index[i][j*D + l] = Xl[l];
						}
						else {
							next_index[i][j*D + l] = gauss_index;
						}
					}
				}
			}
			/*杂交操作*/
			//Crossover operator
			//step 1 of the crossover: determine the length of crossover
			int rand_basic = (int)next_param[i][0];//获取第i个变异个体的簇中心数
			cr_lenth = 0;
			do {
				cr_lenth = cr_lenth + 1;
			} while (URAND < CR && cr_lenth < rand_basic);//在杂交前，先计算变异个体中用来交叉的簇中心数
			cr_n = (int)(rand_basic*URAND); //随机一个杂交点，取值范围是0~rand_basic-1  
			//step 2 of the crossover: determine the subspace of crossover
			int rand1;
			swap = 0, swap1 = 0;    //the number of points that outside and inside the swap area
			if ((cr_n + cr_lenth - 1) < rand_basic) {//rand_basic是当前变异个体的簇中心数，cr_n+cr_length是当前变异个体最后一个簇的位置
				for (int j = cr_n; j<cr_n + cr_lenth; j++) {
					for (int l = 0; l<D; l++)
						in_area[i][(j - cr_n)*D + l] = next_index[i][j*D + l]; //in_area二维数组保存的是每个变异个体待杂交的簇中心，而且是从0开始连续保存
				}
			}
			else {//如果从杂交开始点cr_n开始杂交一直杂交cr_length个杂交向量，而杂交结束点已经超出了当前变异个体的簇数  
				for (int j = cr_n; j<rand_basic; j++) {
					for (int l = 0; l<D; l++)
						in_area[i][(j - cr_n)*D + l] = next_index[i][j*D + l];
				}
				for (int j = 0; j<(cr_n + cr_lenth - rand_basic); j++) { 
					for (int l = 0; l<D; l++)
						in_area[i][(j + rand_basic - cr_n)*D + l] = next_index[i][j*D + l];
				}
			}
			/*进行子空间交叉*/
			if (cr_lenth == 1) {//需要杂交的簇只有一个的情况，
				for (int j = 0; j<20 * D; j++)
					out_area[i][j] = *(p + i * 20 * D + j); //先把当前目标个体的所有簇保存到数组out_area中。
				//如果对位的目标个体是满簇情况，即20个簇，那么子空间交叉时，只需要把需要杂交的簇替换对位目标个体随机的一个簇。  
				if ((int)*(p2 + i * 2) == 20) {    
					rand1 = (int)(((int) *(p2 + i * 2)) * URAND);
					for (int j = 0; j<D; j++)
						out_area[i][rand1*D + j] = in_area[i][j]; //如果cr_lenth==1,即交叉的簇只有一个，那么数组in_area[i]只保存了一个簇，然后把这个簇替换到out_area数组的随机一个位置
					for (int j = 0; j<20 * D; j++)
						next_index[i][j] = out_area[i][j]; //把out_area数组中的簇复制到保存所有变异个体的数组，即得到试验向量的数组。
					next_param[i][0] = *(p2 + i * 2); //next_param数组保存的是第i个变异个体的簇数，这里因为杂交的簇只有一个因此第i个变异个体（试验向量）的簇数与对位的目标个体的簇数保持一致
				}
				//对位的目标个体的簇数，小于20个，即没满簇，那么子空间杂交，只需要把待杂交的簇插入到目标个体的末尾。
				else { 
					for (int j = 0; j<D; j++)
						out_area[i][((int)*(p2 + i * 2)) *D + j] = in_area[i][j];
					for (int j = 0; j<20 * D; j++)
						next_index[i][j] = out_area[i][j];
					next_param[i][0] = *(p2 + i * 2) + 1;
				}
			}
			//需要杂交的簇数大于1，即需要使用子空间杂交法，首先需要得到子空间的均值向量（中心向量）
			else {    
				for (int j = 0; j<D; j++)
					center[j] = 0; //子空间的中心向量初始值为全0

				//average all node to get the swap center
				for (int j = 0; j<cr_lenth; j++) {
					for (int l = 0; l<D; l++) {//将需要杂交的所有簇，叠加起来，每维求平均值之后，就可以得到均值向量
						center[l] = center[l] + in_area[i][j*D + l];
					}
				}

				for (int j = 0; j<D; j++) {//保存子空间中的第一个簇与最后一个簇的对位距离
					dist[j] = fabs(in_area[i][0 * D + j] - in_area[i][(cr_lenth - 1)*D + j]) / 2; //first and last node
					center[j] = center[j] / cr_lenth;  //average of all node
				}


				for (int j = 0; j<(int) *(p2 + i * 2); j++) {
					swap_if = 0;
					for (int l = 0; l<D; l++) {
						//满足条件，则说明当前的簇不在杂交子空间范围内
						if (fabs(*(p + i * 20 * D + j * D + l) - center[l]) > dist[l]) {   //the center can be changed
							swap_if = 1;
							break;
						}
					}
					if (swap_if == 1) {//把不需要被覆盖的簇保存起来，即这部分簇是目标个体中不需要被变异个体替换的簇。随后这部分簇与变异个体中被选中的杂交簇结合成试验向量
						for (int l = 0; l<D; l++)
							out_area[i][swap1*D + l] = *(p + i * 20 * D + j * D + l);
						swap1 = swap1 + 1;
					}
				}
				swap = cr_lenth + swap1;//swap1是目标个体中不需要被交换的簇数，cr_lenth是指变异个体子空间（数组in_area）的簇数，即用来覆盖的簇数

				//step 3 of the crossover : subarea swap
				if (swap <= 20 && swap >= 2) {
					for (int j = 0; j<cr_lenth*D; j++)
						next_index[i][j] = in_area[i][j]; //先把变异个体子空间中的簇，插入到变异向量中
					for (int j = cr_lenth * D; j<swap*D; j++)
						next_index[i][j] = out_area[i][j - cr_lenth * D];//再把目标个体中不需要被覆盖的簇插入到变异向量中，该操作完成即原本保存变异向量的数组就变成保存试验向量了。
					next_param[i][0] = swap;
				}
				/*当杂交后的簇数大于20或小于2的时候，这里是直接用变异向量作为试验向量*/
			}

			//Selection operator,将样本分配到具体的簇中
			//step 1 of the selection: assign the data object
			for (int j = 0; j<N; j++) {
				min = DBL_MAX;    
				for (int l = 0; l<((int)next_param[i][0]); l++) {   
					distance = getDistance(data[j], &next_index[i][l*D], D);
					if (distance < min) {
						min = distance;
						in_cluster2[i][j] = l;    //record the data index
					}
				}
			}
			/*这里应该还要对样本数小于2的簇进行处理。*/

			//step 2 of the selection: calculate the fitness，next_param[i][0]记录了当前试验向量的簇数，data是样本数据集，next_index[i]保存了当前这个试验向量的各个簇中心，in_cluster2[i]数据集中每个数据属于哪个簇。
			next_param[i][1] = func(next_param[i][0], *data, next_index[i], in_cluster2[i], N, D);	
			numofE = numofE + 1;

			//step 3 of the selection: preserve the better one (depends on the testing index)
			if (next_param[i][1] > *(p2 + i * 2 + 1)) {  
				for (int j = 0; j<20 * D; j++) {
					*(p + i * 20 * D + j) = next_index[i][j];
				}
				*(p2 + i * 2 + 0) = next_param[i][0];
				*(p2 + i * 2 + 1) = next_param[i][1];
			}
			//适应度越大越好，符合i指标，若用DB指标则是越小越好
			if (*(p2 + i * 2 + 1) > best_val) {     
				best_val = *(p2 + i * 2 + 1);
				best_val2 = *(p2 + i * 2 + 0);
				index = i;
			}
		}
	}
	// return index;
	return 1;
}


/*
生成n个不同的随机数
随机数范围0-rang-1
*/
int* rand_generator(int n,int rang){
	if(n>=rang){
		exit(1);
	}
	int *a;
	int *b;
	int rand_d;
	malloc1E(a,n);
	malloc1E(b,rang);
	for(int i=0;i<rang;i++){
		b[i]=0;
	}
	for(int i=0;i<n;i++){
		rand_d=(int)rand()%rang;
		while(b[rand_d]!=0){
			rand_d=(int)rand()%rang;
		}
		b[rand_d]=1;
		a[i]=rand_d;
	}  
	return a;
	
}

void printf_data(double **data,int D,int N){
   
    printf("\n输出所有样本的坐标\n");
   
    for(int i=0;i<N;i++){
            printf("第%d个样本坐标：",i);  
         for(int j=0;j<D;j++){
            printf("%lf\t",data[i][j]);
         }
            printf("\n");
    }

}

void printf_pop(double **pop,double **p,int D,int np){
   
    printf("\n输出所有簇中心的坐标\n");
   
    for(int i=0;i<np;i++){
            printf("第%d个个体坐标：",i);  
         for(int j=0;j<p[i][0]*D;j++){
            if(j%D==0&&j!=0){
            printf(" || ");
            }
            
            printf("%lf\t",pop[i][j]);

         }
            printf("\n");
    }

}

void printf_inclu(int**inclu,double **p,int np,int N){
     printf("\n输出聚类情况\n");
    for(int i=0;i<np;i++){
          printf("第%d个个体的聚类情况,一共有%g个簇\t",i,p[i][0]);
        for (int j = 0; j<N; j++)
        {
          printf("样本%d->编号:%d\t",j,inclu[i][j]);    
        }
          printf("\n");
    }
}


int  main()
{
	srand((unsigned int)(time(NULL)));

	int i, j, D, N, Gmax, NP, best = 0, *popul_rand, **in_cluster;
	double data_min, data_max, min, distance, *uk, *lk, **popul_index, **popul_param;
		const char *filename="iris.txt";
	data = loadData(filename,N, D);   //load data from the text file,将数据集存到了二维数组data中
		
	NP = 10 * D; //维数越大，问题的解空间就越大，因此种群的个数就应该更大
	Gmax = 1000000/ NP;//每一代都需要进行NP次评估，因此该式子能得到最大的迭代代数
	printf("The times of iteration(Gmax):%d\n", Gmax);

	malloc1D(uk, D);
	malloc1D(lk, D);
	malloc1E(popul_rand, N);

	malloc2D(in_cluster, NP, N);        //population index，二维数组in_cluster用来记录当前样本在哪个簇的情况
	malloc2E(popul_index, NP, 20 * D);  //population cluster centroids
	malloc2E(popul_param, NP, 2);       //population info(cluster number, fitness value)

	for (i = 0; i < NP; i++) {
		for (j = 0; j < N; j++) {
			in_cluster[i][j] = 0;
		}
	}
	//data是保存所有样本的二维数组。该循环是求每一维度的最大值和最小值，然后将该维度的最大值最小值分别设为该维度的上界和下界
	for (i = 0; i<D; i++) {
		data_min = DBL_MAX, data_max = DBL_MIN;
		for (j = 0; j<N; j++) {
			if (data[j][i] > data_max)
				data_max = data[j][i];
			if (data[j][i] < data_min)
				data_min = data[j][i];
		}
		
		uk[i] = data_max;
		lk[i] = data_min;
	}

	//population initialization
	for (i = 0; i<NP; i++) {     
		int k_rand = rand() % 19 + 2;   //the cluster num is between[2,20]
		popul_param[i][0] = (double)k_rand; //popul_param每一行都有两列，第一列是保存该行对应的簇的簇数，第二列对应的是该行对应的簇的适应度
		 int *rand_arr=rand_generator(k_rand,N);
		for (j = 0; j<k_rand; j++) {  
			for (int k = 0; k < D; k++) {
				// popul_index[i][j*D + k] = lk[k] + URAND *(uk[k] - lk[k]);//初始化操作，随机了k_rand个簇，这里是每个簇的每一维度都随机初始化一个取值范围内的值
				popul_index[i][j*D + k]=data[rand_arr[j]][k];
			}		
		}

		//evaluate the fitness of the initial population
        for (int j = 0; j<N; j++) {
             
			min = DBL_MAX;   
			
            for (int l = 0; l<((int)popul_param[i][0]); l++) { 
				distance = getDistance(data[j], &popul_index[i][l*D], D);    
				// printf("\n样本%d与簇%d的距离为%lf",j,l,distance);
				if (distance < min) {
					min = distance;
					in_cluster[i][j] = l;    //record the data index，用来记录当前样本在当前个体的下标为l的簇中。
                }
			}
            	
		}
       
		popul_param[i][1] = func(popul_param[i][0], *data, popul_index[i], in_cluster[i], N, D);//计算当前个体的适应度	
	}
	    
	best = e_de(*popul_index, *popul_param, N, D, Gmax, lk, uk); //run the E-DE algorithm

	// //output the results as your need
	printf("E-DE run successful!\n");

	return 0;
}



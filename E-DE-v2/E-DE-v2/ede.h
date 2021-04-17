//#ifndef _CLUSTER_NUM
//#define MIN_CLUSTER 2  
//#define	MAX_CLUSTER 20
//#endif

#ifndef MIN_CLUSTER
#define MIN_CLUSTER 2  
#endif // !MIN_CLUSTER

#ifndef MAX_CLUSTER
#define	MAX_CLUSTER 20
#endif

#ifndef MAX_TEST_NUM
#define MAX_TEST_NUM 5
#endif // !MAX_TEST_NUM

#ifndef DB_INDEX
#define DB_INDEX 0
#endif 

#ifndef I_INDEX
#define I_INDEX 1
#endif // !I_INDEX

#ifndef SIHOUETTES_INDEX
#define SIHOUETTES_INDEX 2
#endif  

#ifndef GET_CLUSTER_NUM
#define GET_CLUSTER_NUM 0
#endif // !GET_CLUSTER_NUM

#ifndef GET_FITNESS 
#define GET_FITNESS 1
#endif

#ifndef TRAIL_DATA
#define  TRAIL_DATA

#ifndef URAND
// change any of these parameters to match your needs 
#define URAND  ((double)rand()/((double)RAND_MAX+1.0)) //产生的随机数是[0,1),rand()产生0-RAND_MAX范围的数
#endif
typedef struct{
	const char* filename;
	double bestfitness[MAX_TEST_NUM];
	int bestClusterNum[MAX_TEST_NUM];
	double ARI_arr[MAX_TEST_NUM];
	double Interdist_arr[MAX_TEST_NUM];
	double Intradist_arr[MAX_TEST_NUM];
	double correctClusterNum;//保存当前数据集正确的簇数
	int length;
	int usedindex;//记录当前使用的评价指标
} trailData;
#endif // !1



#ifndef _GAUSS_RAND
#define  _GAUSS_RAND

/*产生满足均值为0，标准差为1的正态分布随机数*/
double gaussrand_NORMAL();
/*产生满足均值为mean，标准差为stdc的正态分布随机数*/
double gaussrand(double mean, double stdc);
#endif

#ifndef _RAND_GENERATOR
#define _RAND_GENERATOR
/*生成n个不同的随机数,随机数范围0-rang-1*/
int* rand_generator(int n, int rang);
#endif


#ifndef _MEMORY_ALLOCATION
#define	_MEMORY_ALLOCATION
/*给一维double型数组分配长度为D个单元的内存*/
void malloc1D(double*& a, int D);
/*给一维int型数组分配长度为D个单元的内存*/
void malloc1E(int*& a, int D);
/*给二维double型数组分配长度为xDim行yDim列个单元的内存*/
void malloc2D(int**& a, int xDim, int yDim);
/*给二维int型数组分配长度为xDim行yDim列个单元的内存*/
void malloc2E(double**& a, int xDim, int yDim);
char* malloc1Char( int D);
#endif


#ifndef _FITNESS_FUNCS
#define _FITNESS_FUNCS
/*计算每个簇的均值向量*/
double* count_meanV(double* data, int* in_clu, int k_num, int N, int D);
/*I指标函数*/
double I_index(int k_num, double* data, double* p_index, int* in_clu, int N, int D);
/*DB指标函数*/
double DB_index(int k_num, double* p_pars, double* p_index, int* in_clu, int N, int D);
/*求ARI值*/
double getARI(int* originClusterInfo, int* in_cluster, int k_num, int N);
double getIntraDist(double** data, double* p_index, double* p_param, int* in_cluster, int N, int D);
double getInterDist(double* p_index, double* p_param, int D);
double** getDataDistance(double** data, int N, int D);
void readyforSilhouettes(double* sihouettes, int* adjacent_cluster, int* incluster,  int N, int k_num, int D, double** distarr);
double getSilhouettes(int* incluster, int N, int k_num, int D, double** distarr);
#endif


#ifndef _FILE_LOAD
#define _FILE_LOAD
/*统计文件中的样本数和维度数*/
void count_N_D(const char* filename, int& n, int& d);

double** loadData(const char* filename, int& n, int& d);
double** loadData2(const char* filename, int& n, int& d);
int* getDataOriginClusterInfo(const char* filename, int N);

#endif

#ifndef _TOOL
#define _TOOL
int isExistEmptyCluster(int* incluster, int N, int K);
int* getDataNumInCluster(int* incluster, int N, int K);
double* getDdimCluCenter(double** data, int N, int K, int D);
double getDistance(double* avector, double* bvector, int n);
void dataInCluster(double** data, int N, int K, double* cluster, int D,int* in_cluster);
int getMaxClusterNum(int N, int maxClusterNum);
double getMinFitness(double* param, int np, int flag);
double getMaxFitness(double* param, int np, int flag);
void dealWithClusterNumNotCorrect(double* p_index, double* p_param, double* next_index, double* next_param, int D);
void freeTwoDimArr_double(double** &tdd, int x);
void freeTwoDimArr_int(int** &tdi, int x);
void dealwith_emptyCluster(int* in_cluster, int N, double& k_num, double** data, int D, double* p_index);

#endif

#ifndef RECORD_TRAIL_DATA
#define RECORD_TRAIL_DATA
// 记录多次测试的实验数据
void recordTraildata(trailData& traildata, double ARI, double interDist, double intraDist, double bestFitness, double bestClusterNum);
char* cutFileNameSuffix(const char* filename, const char* suffix);
char* getOutputFilename(const char* filename, int index, int i);
char* getIndexName(int indexNum);
#endif 


#ifndef INITIALIZATION
#define INITIALIZATION
void set_traildataInfo(trailData& traildata, int correctClusterNum, const char* filename, int usedindex);
int getNP(int D, int times);
int getGmax(int NP, int divisor);
void initial_two_Dim_intarr(int** arr, int x, int y, int value);
void getDataDim_max_min(double** data, double* uk, double* lk, int D, int N);
void initial_individual(double** data, double* popul_index, double* popul_param, int* in_cluster, int N, int D, int max_cluster_num, int min_cluster_num);
#endif 
#ifndef CORE_OPERATION
#define CORE_OPERATION
void getBestValByIndex(double& bestFitness, double& bestClusterNum, double* popul_param, int NP, int usedIndex);
void initialCoreVriable(double** next_param, int** in_cluster2, double** out_area, double** in_area, double** next_index, int NP, int max_cluster_num, int D, int N);
void initial_dist_center(double* dist, double* center, int D);
void getRandomNum(int& r1, int& r2, int& r3, int i, int NP);
int cross_border_process(int num, int ub, int lb);
#endif 
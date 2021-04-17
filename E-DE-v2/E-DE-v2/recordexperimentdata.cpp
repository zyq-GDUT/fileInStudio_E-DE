#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include<math.h>
#include <float.h>
#include <direct.h>
#include <conio.h>
#include "ede.h"
#include <corecrt_io.h>


char* creatFolderWithPrefix_experimentdata(char* folder);
char* cutFileNameSuffix(const char* filename, const char* suffix);
char* getIndexName(int indexNum);
char* cutFileNameSuffix(const char* filename, const char* suffix);
void processOutputFilename(char* outputfile, int i, const char* suffix);
char* creatTwoFolder(char* firstFolder, int index);
/*
* 将实验数据输出到文件中
* 该方法输出单次测试的数据
*/
void outputExperimentdata1(trailData& traildata,FILE* fp) {
	
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		fprintf(fp,"\t\t第%d次  ",i+1);
	}
	fprintf(fp, "\n");
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		if (i == 0) {
			fprintf(fp, "bestClusterNum:  ");
		}
		fprintf(fp, "%12d ", traildata.bestClusterNum[i]);

	}
	fprintf(fp, "\n");
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		if (i == 0) {
			fprintf(fp, "bestFitness：\t");
		}
		fprintf(fp, "%10.4f ", traildata.bestfitness[i]);

	}
	fprintf(fp, "\n");
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		if (i == 0) {
			fprintf(fp, "IntraDistance：\t");
		}
		fprintf(fp, "%10.6g ", traildata.Intradist_arr[i]);

	}
	fprintf(fp, "\n");
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		if (i == 0) {
			fprintf(fp, "InterDistance：\t");
		}
		fprintf(fp, "%10.6g ", traildata.Interdist_arr[i]);

	}
	fprintf(fp, "\n");
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		if (i == 0) {
			fprintf(fp, "ARI：\t\t");
		}
		fprintf(fp, "%10.6f ", traildata.ARI_arr[i]);

	}
	fprintf(fp, "\n");
	
}
/*
* 该方法输出多次测试的数据到文件，如均值和标准差
*/
void outputExperimentdata2(double avgARI, double avgInterDist, double avgIntraDist, double avgfitness, int* clusterCount, double stdARI, double stdInterDist, double stdIntraDist, double stdfitness,FILE* fp) {
	fprintf(fp, "\n%d次测试的聚类情况[簇数(出现次数)]：", MAX_TEST_NUM);
	for (int i = 0; i < MAX_CLUSTER+1; i++) {
		if (clusterCount[i] == 0) continue;
		fprintf(fp, "\t%d(%d) ", i, clusterCount[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "平均适应度(AVG.Fitness): %g\t标准差: ±%g\n",avgfitness,stdfitness);
	fprintf(fp, "平均簇内距离(Intra Distance): %g\t标准差: ±%g\n", avgIntraDist,stdIntraDist);
	fprintf(fp, "平均簇间距离(Inter Distance): %g\t标准差: ±%g\n", avgInterDist,stdInterDist);
	fprintf(fp, "平均ARI值(AVG.ARI): \t%10g\t\t标准差: ±%g\n", avgARI,stdARI);
	
}

/*
* 计算数据的平均值标准差等
*/
void manageExperimentData(trailData& traildata) {
	double avgARI = 0.0;
	double stdARI = 0.0;
	double avgInterDist = 0.0;
	double stdInterDist = 0.0;
	double avgIntraDist = 0.0;
	double stdIntraDist = 0.0;
	double avgfitness = 0.0;
	double stdfitness = 0.0;
	double CSR = 0.0;//三十次测试的平均聚类成功率，该成功率是指找到正确簇数的成功率
	int clusterCount[MAX_CLUSTER+1] = { 0 };//下标表示簇数，值表示该簇数的出现的次数
	const char* recordFilename = "epData.txt";
	char* str=cutFileNameSuffix(traildata.filename, ".txt");
	char* folderName=creatFolderWithPrefix_experimentdata(str);
	char* indexname=getIndexName(traildata.usedindex);
	strcat(indexname, recordFilename);
	strcat(str, indexname);
	//给文件加上文件夹
	strcat(folderName, str);
	//输出到文件
	//freopen(folderName, "w", stdout);
	//求平均数
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		avgARI += traildata.ARI_arr[i];
		avgInterDist += traildata.Interdist_arr[i];
		avgIntraDist += traildata.Intradist_arr[i];
		avgfitness += traildata.bestfitness[i];
		clusterCount[traildata.bestClusterNum[i]]++;
	}
	avgARI = avgARI / MAX_TEST_NUM;
	avgInterDist = avgInterDist / MAX_TEST_NUM;
	avgIntraDist = avgIntraDist / MAX_TEST_NUM;
	avgfitness = avgfitness / MAX_TEST_NUM;
	//求标准差
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		stdARI += pow(traildata.ARI_arr[i] - avgARI, 2);
		stdInterDist += pow(traildata.Interdist_arr[i] - avgInterDist, 2);
		stdIntraDist += pow(traildata.Intradist_arr[i] - avgIntraDist, 2);
		stdfitness += pow(traildata.bestfitness[i] - avgfitness, 2);
	}
	stdARI = sqrt(stdARI / MAX_TEST_NUM);
	stdInterDist = sqrt(stdInterDist / MAX_TEST_NUM);
	stdIntraDist = sqrt(stdIntraDist / MAX_TEST_NUM);
	stdfitness = sqrt(stdfitness / MAX_TEST_NUM);
	FILE* fp;
	if ((fp = fopen(folderName, "w")) == NULL) {
		printf("文件%s不能打开\n",folderName);
		exit(0);
	}
	fprintf(fp,"处理的文件：%s\n一共测试了%d次\n\n", traildata.filename, MAX_TEST_NUM);
	outputExperimentdata1(traildata,fp);
	outputExperimentdata2(avgARI, avgInterDist, avgIntraDist, avgfitness, clusterCount, stdARI, stdInterDist, stdIntraDist, stdfitness,fp);
	fclose(fp);
	free(str);
	free(folderName);
	free(indexname);
}

/*
* 记录多次测试的实验数据
* ARI,interDist,intraDist,bestFitness,bestCluster都是一次测试下的的数据
*/
void recordTraildata(trailData& traildata, double ARI, double interDist, double intraDist, double bestFitness, double bestClusterNum) {
	

	
	traildata.ARI_arr[traildata.length] = ARI;
	traildata.Interdist_arr[traildata.length] = interDist;
	traildata.Intradist_arr[traildata.length] = intraDist;
	traildata.bestfitness[traildata.length] = bestFitness;
	traildata.bestClusterNum[traildata.length] = bestClusterNum;
	traildata.length++;

	
	//如果运行到最后一次则处理实验数据并输出到文件
	if (traildata.length == MAX_TEST_NUM) {
		manageExperimentData(traildata);
		traildata.length = 0;
	}
	
}

/*
*通过I指标或DB指标的编号返回字符串
*/
 char* getIndexName(int indexNum) {
	 char* indexname = malloc1Char(100);
	if (indexNum==I_INDEX) {
		strcpy(indexname,"_I_index_");
	}
	else if (indexNum == DB_INDEX) {
		strcpy(indexname, "_DB_index_");
	
	}
	else if (indexNum == SIHOUETTES_INDEX) {
		strcpy(indexname, "_SIHOUETTES_index_");

	}
	else {
		free(indexname);
		printf("返回字符串出错！");
		exit(2);
	}
	return indexname;
}



/*
* 去掉文件名的.txt后缀,返回去掉后缀后的字符串
*/
char*  cutFileNameSuffix( const char* filename,const char* suffix) {
    char* str =malloc1Char(100);
	int filenamelen=strlen(filename);
	int suffixlen = strlen(suffix);
	if (suffix[0]!='.'||filenamelen <= suffixlen) {
		printf("文件名没有后缀或输入的后缀没有加\".\"，无法处理");
		exit(1);
	}
	strncpy(str, filename,  filenamelen - suffixlen);
	str[filenamelen -suffixlen] = '\0';
	return str;

}

/*
* 给输出文件名添加后缀
*/
void processOutputFilename(char* outputfile, int i,const char* suffix) {
	if (outputfile == NULL || strlen(suffix) <= 0 || suffix[0] != '.') {
		exit(3);
	}
	char str[10] = { 0 };
	_itoa(i, str, 10);//10是十进制
	strcat(str, suffix);
	strcat(outputfile, str);
}



/*
* 创建文件夹并返回文件夹名称
* 
* prefix是文件夹的前缀
* prefix+folder就是最终的文件夹名称
*/
char* creatFolderName(char* folder, const char* prefix) {
	//char* folderName = (char*)malloc(50);
	char* folderName = malloc1Char(150);
	strcpy(folderName, prefix);
	strcat(folderName, folder);
	// 文件夹不存在则创建文件夹,0表示文件是否存在若存在则返回0，否则返回-1
	if (_access(folderName,0) == -1)
	{
		_mkdir(folderName);
	}
	//int f = _mkdir(folderName);//创建失败返回-1，否则返回0
	strcat(folderName, "/");
	return folderName;
}

/*
* folder是待创建的文件夹名
* 通过该方法创建的文件夹会带有“experiment_data”这样的前缀
* 返回创建后文件夹的名称
*/
char* creatFolderWithPrefix_experimentdata(char* folder) {
	return creatFolderName(folder, "experiment_data_");
}



/*
* 根据评价指标创建二级文件夹并返回二级文件夹的路径
*若文件夹存在则只返回该二级文件夹的路径。
*/

char* creatTwoFolder(char* firstFolder, int index) {
	char* secondFolderName = malloc1Char(100);
	if (index == I_INDEX) {
		strcpy(secondFolderName, "I_index_eve_generation_info");
	}
	else if(index == DB_INDEX) {
		strcpy(secondFolderName, "DB_index_eve_generation_info");
	}
	else if(index == SIHOUETTES_INDEX) {
		strcpy(secondFolderName, "SIHOUETTES_index_eve_generation_info");
	}
	//返回一级文件夹的名称，若该文件夹不存在还会创建
	char* folderName = creatFolderWithPrefix_experimentdata(firstFolder);
	strcat(folderName, secondFolderName);
	free(secondFolderName);
	return creatFolderName(folderName, "");

}

/*
* 获取输出文件名
* filename当前处理的数据集文件名
* index当前使用的评价指标
* i是文件名要添加的编号
*/
char* getOutputFilename(const char* filename,int index,int i) {
	char* outputfilename = cutFileNameSuffix(filename, ".txt");
	//新建文件夹,若存在则不创建
	//char* folderName = creatFolderWithPrefix_experimentdata(outputfilename);
	char* folderName = creatTwoFolder(outputfilename, index);
	//构造输出文件的名称
	char* indexname=getIndexName(index);
	strcat(outputfilename, indexname);
	processOutputFilename(outputfilename, i, ".txt");
	//在文件名前添加文件夹名
	strcat(folderName, outputfilename);
	free(indexname);
	free(outputfilename);
	//返回要收集数据的文件的完整路径
	return folderName;
}


//int main() {
//
//	//char folderName[50] = "testFolder11";
//	//creatFolderName(folderName, "");//若目录已经存在则返回-1，否则返回0
//	//
//	//const char* filename = "test.txt";
//	//char* outputfilename=getOutputFilename(filename, 1, 1);
//	//freopen((const char*)outputfilename,"a",stdout);
//	//printf("hello world!\n");
//	//freopen("CON", "w", stdout);
//	const char* filename = "testsss.txt";
//	char* t=getOutputFilename(filename, 1, 1);
//	printf("%s", t);
//	system("pause");
//	return 0;
//}

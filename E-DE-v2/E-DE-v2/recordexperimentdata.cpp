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
* ��ʵ������������ļ���
* �÷���������β��Ե�����
*/
void outputExperimentdata1(trailData& traildata,FILE* fp) {
	
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		fprintf(fp,"\t\t��%d��  ",i+1);
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
			fprintf(fp, "bestFitness��\t");
		}
		fprintf(fp, "%10.4f ", traildata.bestfitness[i]);

	}
	fprintf(fp, "\n");
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		if (i == 0) {
			fprintf(fp, "IntraDistance��\t");
		}
		fprintf(fp, "%10.6g ", traildata.Intradist_arr[i]);

	}
	fprintf(fp, "\n");
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		if (i == 0) {
			fprintf(fp, "InterDistance��\t");
		}
		fprintf(fp, "%10.6g ", traildata.Interdist_arr[i]);

	}
	fprintf(fp, "\n");
	for (int i = 0; i < MAX_TEST_NUM; i++) {
		if (i == 0) {
			fprintf(fp, "ARI��\t\t");
		}
		fprintf(fp, "%10.6f ", traildata.ARI_arr[i]);

	}
	fprintf(fp, "\n");
	
}
/*
* �÷��������β��Ե����ݵ��ļ������ֵ�ͱ�׼��
*/
void outputExperimentdata2(double avgARI, double avgInterDist, double avgIntraDist, double avgfitness, int* clusterCount, double stdARI, double stdInterDist, double stdIntraDist, double stdfitness,FILE* fp) {
	fprintf(fp, "\n%d�β��Եľ������[����(���ִ���)]��", MAX_TEST_NUM);
	for (int i = 0; i < MAX_CLUSTER+1; i++) {
		if (clusterCount[i] == 0) continue;
		fprintf(fp, "\t%d(%d) ", i, clusterCount[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "ƽ����Ӧ��(AVG.Fitness): %g\t��׼��: ��%g\n",avgfitness,stdfitness);
	fprintf(fp, "ƽ�����ھ���(Intra Distance): %g\t��׼��: ��%g\n", avgIntraDist,stdIntraDist);
	fprintf(fp, "ƽ���ؼ����(Inter Distance): %g\t��׼��: ��%g\n", avgInterDist,stdInterDist);
	fprintf(fp, "ƽ��ARIֵ(AVG.ARI): \t%10g\t\t��׼��: ��%g\n", avgARI,stdARI);
	
}

/*
* �������ݵ�ƽ��ֵ��׼���
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
	double CSR = 0.0;//��ʮ�β��Ե�ƽ������ɹ��ʣ��óɹ�����ָ�ҵ���ȷ�����ĳɹ���
	int clusterCount[MAX_CLUSTER+1] = { 0 };//�±��ʾ������ֵ��ʾ�ô����ĳ��ֵĴ���
	const char* recordFilename = "epData.txt";
	char* str=cutFileNameSuffix(traildata.filename, ".txt");
	char* folderName=creatFolderWithPrefix_experimentdata(str);
	char* indexname=getIndexName(traildata.usedindex);
	strcat(indexname, recordFilename);
	strcat(str, indexname);
	//���ļ������ļ���
	strcat(folderName, str);
	//������ļ�
	//freopen(folderName, "w", stdout);
	//��ƽ����
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
	//���׼��
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
		printf("�ļ�%s���ܴ�\n",folderName);
		exit(0);
	}
	fprintf(fp,"������ļ���%s\nһ��������%d��\n\n", traildata.filename, MAX_TEST_NUM);
	outputExperimentdata1(traildata,fp);
	outputExperimentdata2(avgARI, avgInterDist, avgIntraDist, avgfitness, clusterCount, stdARI, stdInterDist, stdIntraDist, stdfitness,fp);
	fclose(fp);
	free(str);
	free(folderName);
	free(indexname);
}

/*
* ��¼��β��Ե�ʵ������
* ARI,interDist,intraDist,bestFitness,bestCluster����һ�β����µĵ�����
*/
void recordTraildata(trailData& traildata, double ARI, double interDist, double intraDist, double bestFitness, double bestClusterNum) {
	

	
	traildata.ARI_arr[traildata.length] = ARI;
	traildata.Interdist_arr[traildata.length] = interDist;
	traildata.Intradist_arr[traildata.length] = intraDist;
	traildata.bestfitness[traildata.length] = bestFitness;
	traildata.bestClusterNum[traildata.length] = bestClusterNum;
	traildata.length++;

	
	//������е����һ������ʵ�����ݲ�������ļ�
	if (traildata.length == MAX_TEST_NUM) {
		manageExperimentData(traildata);
		traildata.length = 0;
	}
	
}

/*
*ͨ��Iָ���DBָ��ı�ŷ����ַ���
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
		printf("�����ַ�������");
		exit(2);
	}
	return indexname;
}



/*
* ȥ���ļ�����.txt��׺,����ȥ����׺����ַ���
*/
char*  cutFileNameSuffix( const char* filename,const char* suffix) {
    char* str =malloc1Char(100);
	int filenamelen=strlen(filename);
	int suffixlen = strlen(suffix);
	if (suffix[0]!='.'||filenamelen <= suffixlen) {
		printf("�ļ���û�к�׺������ĺ�׺û�м�\".\"���޷�����");
		exit(1);
	}
	strncpy(str, filename,  filenamelen - suffixlen);
	str[filenamelen -suffixlen] = '\0';
	return str;

}

/*
* ������ļ�����Ӻ�׺
*/
void processOutputFilename(char* outputfile, int i,const char* suffix) {
	if (outputfile == NULL || strlen(suffix) <= 0 || suffix[0] != '.') {
		exit(3);
	}
	char str[10] = { 0 };
	_itoa(i, str, 10);//10��ʮ����
	strcat(str, suffix);
	strcat(outputfile, str);
}



/*
* �����ļ��в������ļ�������
* 
* prefix���ļ��е�ǰ׺
* prefix+folder�������յ��ļ�������
*/
char* creatFolderName(char* folder, const char* prefix) {
	//char* folderName = (char*)malloc(50);
	char* folderName = malloc1Char(150);
	strcpy(folderName, prefix);
	strcat(folderName, folder);
	// �ļ��в������򴴽��ļ���,0��ʾ�ļ��Ƿ�����������򷵻�0�����򷵻�-1
	if (_access(folderName,0) == -1)
	{
		_mkdir(folderName);
	}
	//int f = _mkdir(folderName);//����ʧ�ܷ���-1�����򷵻�0
	strcat(folderName, "/");
	return folderName;
}

/*
* folder�Ǵ��������ļ�����
* ͨ���÷����������ļ��л���С�experiment_data��������ǰ׺
* ���ش������ļ��е�����
*/
char* creatFolderWithPrefix_experimentdata(char* folder) {
	return creatFolderName(folder, "experiment_data_");
}



/*
* ��������ָ�괴�������ļ��в����ض����ļ��е�·��
*���ļ��д�����ֻ���ظö����ļ��е�·����
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
	//����һ���ļ��е����ƣ������ļ��в����ڻ��ᴴ��
	char* folderName = creatFolderWithPrefix_experimentdata(firstFolder);
	strcat(folderName, secondFolderName);
	free(secondFolderName);
	return creatFolderName(folderName, "");

}

/*
* ��ȡ����ļ���
* filename��ǰ��������ݼ��ļ���
* index��ǰʹ�õ�����ָ��
* i���ļ���Ҫ��ӵı��
*/
char* getOutputFilename(const char* filename,int index,int i) {
	char* outputfilename = cutFileNameSuffix(filename, ".txt");
	//�½��ļ���,�������򲻴���
	//char* folderName = creatFolderWithPrefix_experimentdata(outputfilename);
	char* folderName = creatTwoFolder(outputfilename, index);
	//��������ļ�������
	char* indexname=getIndexName(index);
	strcat(outputfilename, indexname);
	processOutputFilename(outputfilename, i, ".txt");
	//���ļ���ǰ����ļ�����
	strcat(folderName, outputfilename);
	free(indexname);
	free(outputfilename);
	//����Ҫ�ռ����ݵ��ļ�������·��
	return folderName;
}


//int main() {
//
//	//char folderName[50] = "testFolder11";
//	//creatFolderName(folderName, "");//��Ŀ¼�Ѿ������򷵻�-1�����򷵻�0
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

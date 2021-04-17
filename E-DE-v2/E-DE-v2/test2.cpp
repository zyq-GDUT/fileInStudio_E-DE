#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctime>
#include <malloc.h>
#include <float.h>
#include<string.h>

struct stu
{
	int bestClusterNum[30];
	int length;
}t = { {1},2 };

int main() {
	printf("%d,%d\n", t.bestClusterNum[0], t.length);
	system("pause");
	return 0;
}
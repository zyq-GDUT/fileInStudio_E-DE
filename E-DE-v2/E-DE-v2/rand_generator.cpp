#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"ede.h"
#define T  1;

/*产生满足均值为0，标准差为1的正态分布随机数*/
double gaussrand_NORMAL() {
	static double V1, V2, S;
	static int phase = 0;
	double X;
	if (phase == 0) {
		do {
			double U1 = (float)rand() / RAND_MAX;
			double U2 = (float)rand() / RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);
		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);
	phase = 1 - phase;
	return X;
}

/*产生满足均值为mean，标准差为stdc的正态分布随机数*/
double gaussrand(double mean, double stdc) {
	return mean + gaussrand_NORMAL() * stdc;
}




/*
生成n个不同的随机数
随机数范围0-rang-1
*/
int* rand_generator(int n, int rang) {
	if (n >= rang) {
		exit(1);
	}
	int* a;
	int* b;
	int rand_d;
	malloc1E(a, n);
	malloc1E(b, rang);
	for (int i = 0; i < rang; i++) {
		b[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		rand_d = (int)rand() % rang;
		while (b[rand_d] != 0) {
			rand_d = (int)rand() % rang;
		}
		b[rand_d] = 1;
		a[i] = rand_d;
	}
	free(b);
	return a;

}
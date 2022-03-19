#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#pragma warning(disable:4996)
#define ITERATIONS 3;

// funtion
long double f(long double x) {
	return x * x * x - exp(x) + 1;
}

long double g(long double x) {
	long double res = x * x * x - exp(x) + 1;
	return fabs(x * x * x - exp(x) + 1);
}

// File and print function

// fill file with values
void FillFile(FILE* fpdata, long double* mass, int size) {
	for (int i = 0; i < size; i++) {
		fprintf(fpdata, "%.20lf ", (mass[i]));
	}
	fprintf(fpdata, "\n");
}

// print array
void PrintMass(long double* mass, int size) {
	for (int i = 0; i < size; i++) {
		printf("%lf ", mass[i]);
	}
	printf("\n");
}

long double* FillH(long double a, long double b, long double from, long double to, int amount) {
	long double* hmass = (long double*)malloc(sizeof(long double) * amount);
	int k = 0;
	for (int i = from; i <= to; i = i * 2, k++) {
		hmass[k] = (b - a) / i;
	}
	return hmass;
}

//interval Runge functions
long double* SplitInterval(long double a, long double b, long long n) {
	long double* xmass = (long double*)malloc(sizeof(long double) * (n + 1));
	long double h = (b - a) / n;
	xmass[0] = a;
	xmass[n] = b;
	for (int i = 1; i < n; i++) {
		xmass[i] = a + i * h;
	}
	return xmass;
}

// simson functions
long double CountH(long double a, long double b, long long n) {
	return (b - a) / n;
}

//long double CountUnevenSummPrev(long double* xmass, double(*f)(double), long long n) {
//	long double summ = 0;
//	for (int k = 1; k <= (n / 2); k++) {
//		summ += f(xmass[2 * k - 1]);
//	}
//	return summ;
//}
//
//long double CountEvenSummPrev(long double* xmass, double(*f)(double), long long n) {
//	long double summ = 0;
//	for (int k = 1; k <= ((n / 2) - 1); k++) {
//		summ += f(xmass[2 * k]);
//	}
//	return summ;
//}
//
//long double Simpson(long double* xmass, long long n, long double (*f)(long double)) {
//	long double h = CountH(xmass[0], xmass[n], n);
//	long double evenSumm = CountEvenSummPrev(xmass,f, n);
//	long double unevenSumm = CountUnevenSummPrev(xmass,f, n);
//	long double res = (h / 3) * (f(xmass[0]) + f(xmass[n]) + 4 * unevenSumm + 2 * evenSumm);
//	return res;
//}

long double CountUnevenSumm(long double a, long double h, long long n, long double(*f)(long double)) {
	long double summ = 0;
	long double x = a;
	for (int k = 1; k <= (n / 2); k++) {
		x = a + (2 * k - 1) * h;
		summ += f(x);
	}
	return summ;
}

long double CountEvenSumm(long double a, long double h, long long n, long double(*f)(long double)) {
	long double summ = 0;
	long double x = a;
	for (int k = 1; k <= ((n / 2) - 1); k++) {
		x = a + 2 * k * h;
		summ += f(x);
	}
	return summ;
}

long double Simpson(long double a, long double b, long double (*f)(long double), long double eps, long double* res) {
	long long n = 1;
	int m = 4;
	long double I1 = 0, I2 = 0;

	long double h = CountH(a, b, n);
	long double evenSumm = CountEvenSumm(a, h, n, f);
	long double unevenSumm = CountUnevenSumm(a, h, n, f);
	I2 = (h / 3) * (f(a) + f(b) + 4 * unevenSumm + 2 * evenSumm);
	//printf("%d: %lf", n, I2);
	do {
		I1 = I2;
		n *= 2;
		h = CountH(a, b, n);
		evenSumm = CountEvenSumm(a, h, n, f);
		unevenSumm = CountUnevenSumm(a, h, n, f);
		I2 = (h / 3) * (f(a) + f(b) + 4 * unevenSumm + 2 * evenSumm);
	} while (fabs(I1 - I2) / (powl(2, m) - 1) > eps);

	printf("\n%d", n);
	*res = I2;
	return h;
}


int main() {
#ifdef ITERATIONS
	char* FileNameSmooth = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab3_Data.csv";
	char* FileNameBroken = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab3_DataModule.csv";
	FILE* dataFile = fopen(FileNameSmooth, "a");
	long double a = 1.2;
	long double b = 2;
	long double h = 0;
	long double eps = 1e-1;
	long long n = 0;  // amount intervals on fixed gap
	int amount = 15;
	long double res;
	long double* hmass = (long double*)malloc(sizeof(long double) * amount);
	long double* resmass = (long double*)malloc(sizeof(long double) * amount);
	for (int k = 0; k < amount; k++) {
		h = Simpson(a, b, f, eps, &res);
		printf("\n%lf\n", res);
		resmass[k] = res;
		hmass[k] = h;
		eps = eps / 10;

	}
	printf("\n\neps = %.15lf", eps);
	FillFile(dataFile, hmass, amount);
	FillFile(dataFile, resmass, amount);
	fclose(dataFile);
#endif // ITERATIONS

#ifndef ITERATIONS
	char* FileNameSmooth = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab3_Data_ChangeN.csv";
	char* FileNameBroken = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab3_DataModule_ChangeN.csv";
	FILE* dataFile = fopen(FileNameBroken, "a");
	long double a = 1.2;
	long double b = 2;
	int from = 2, to = 65536;
	long long n = 0;  // amount intervals on fixed gap
	int amount = log(to) / log(from);
	double* hmass = FillH(a, b, from, to, amount);
	double* resmass = (double*)malloc(sizeof(double) * amount);
	FillFile(dataFile, hmass, amount);
	for (int k = 0, n = from; n <= to; n = n * 2, k++) {
		double* xmass = SplitInterval(a, b, n);
		double res = Simpson(xmass, n, g);
		resmass[k] = res;

	}
	FillFile(dataFile, resmass, amount);
	fclose(dataFile);

#endif // ITERATIONS
}

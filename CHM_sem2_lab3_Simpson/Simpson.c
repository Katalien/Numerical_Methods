#include<stdio.h>
#include <math.h>
#include<stdlib.h>
#pragma warning(disable:4996)
//#define ITERATIONS 3;

// funtion
double f(double x) {
	return x * x * x - exp(x) + 1;
}

double g(double x) {
	return fabs(x * x * x - exp(x) + 1);
}

// File and print function

// fill file with values
void FillFile(FILE* fpdata, double* mass, int size) {
	for (int i = 0; i < size; i++) {
		fprintf(fpdata, "%.15lf ", (mass[i]));
	}
	fprintf(fpdata, "\n");
}

// print array
void PrintMass(double* mass, int size) {
	for (int i = 0; i < size; i++) {
		printf("%lf ", mass[i]);
	}
	printf("\n");
}

double* FillH(double a, double b, double from, double to, int amount) {
	double* hmass = (double*)malloc(sizeof(double) * amount);
	int k = 0;
	for (int i = from; i <= to; i = i * 2, k++) {
		hmass[k] = (b - a) / i;
		printf("\n%d: %lf", i, hmass[k]);
	}
	return hmass;
}

//interval Runge functions
double* SplitInterval(double a, double b, int n) {
	double* xmass = (double*)malloc(sizeof(double) * (n + 1));
	double h = (b - a) / n;
	xmass[0] = a;
	xmass[n] = b;
	for (int i = 1; i < n; i++) {
		xmass[i] = a + i * h;
	}
	return xmass;
}

// simson functions
double CountH(double a, double b, int n) {
	return (b - a) / n;
}

double CountUnevenSumm(double* xmass, int n) {
	double summ = 0;
	for (int k = 1; k <= (n / 2); k++) {
		summ += f(xmass[2 * k - 1]);
	}
	return summ;
}

double CountEvenSumm(double* xmass, int n) {
	double summ = 0;
	for (int k = 1; k <= ((n / 2) - 1); k++) {
		summ += f(xmass[2 * k]);
	}
	return summ;
}

double Simpson(double* xmass, int n, double (*f)(double)) {
	double h = CountH(xmass[0], xmass[n], n);
	double evenSumm = CountEvenSumm(xmass, n);
	double unevenSumm = CountUnevenSumm(xmass, n);
	double res = (h / 3) * (f(xmass[0]) + f(xmass[n]) + 4 * unevenSumm + 2 * evenSumm);
	return res;
}


double* FillY(double* xmass, double(*f)(double), int n) {
	double* ymass = (double*)malloc(sizeof(double) * (n + 1));
	for (int i = 0; i < (n + 1); i++) {
		ymass[i] = f(xmass[i]);
	}
	return ymass;
}

int main() {
#ifdef ITERATIONS

#endif // ITERATIONS

#ifndef ITERATIONS
	char* hFileName = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab3_Data.csv";
	char* FileNameBroken = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab3_DataModule.csv";
	FILE* dataFile = fopen(FileNameBroken, "a");
	double a = 1.2;
	double b = 2;
	int from = 2, to = 50000;
	int n = 0;  // amount intervals on fixed gap
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

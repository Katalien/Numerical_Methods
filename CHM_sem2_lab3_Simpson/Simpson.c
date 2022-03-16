#include<stdio.h>
#include <math.h>
#pragma warning(disable:4996)
//#define ITERATIONS 3;

// funtion
double f(double x) {
	return sqrt(1 + 2 * x * x - x * x * x);
}

// File and print function

// read and fill table function x & y
void FillData(FILE* fpdata, double* xmass, double* ymass, int amount) {
	for (int i = 0; i < amount; i++) {
		fscanf(fpdata, "%lf", &(xmass[i]));
	}
	for (int i = 0; i < amount; i++) {
		fscanf(fpdata, "%lf", &(ymass[i]));
	}
}


// read x values
void ReadFile(FILE* fpdata, double* mass, int amount) {
	for (int i = 0; i < amount; i++) {
		fscanf(fpdata, "%lf", &(mass[i]));
	}
}

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

double CountH(double a, double b, int n) {
	return (b - a) / n;
}

double CountUnevenSumm(double* xmass, int n) {
	double summ = 0;
	printf("\nuneven:");
	for (int k = 1; k <= (n / 2); k++) {
		printf("%lf ", f(xmass[2 * k - 1]));
		summ += f(xmass[2 * k - 1]);
	}
	return summ;
}

double CountEvenSumm(double* xmass, int n) {
	double summ = 0;
	printf("\neven:");
	for (int k = 1; k <= ((n / 2) - 1); k++) {
		printf("%lf ", f(xmass[2 * k]));
		summ += f(xmass[2 * k]);
	}
	return summ;
}

double Simpson(double* xmass, int n) {
	double h = CountH(xmass[0], xmass[n], n);
	printf("%lf\n", h);
	double evenSumm = CountEvenSumm(xmass, n);
	double unevenSumm = CountUnevenSumm(xmass, n);
	double res = (h / 3) * (f(xmass[0]) + f(xmass[n]) + 4 * unevenSumm + 2 * evenSumm);
	return res;
}

int main() {
#ifdef ITERATIONS

#endif // ITERATIONS

#ifndef ITERATIONS
	double a = -2;
	double b = 2;
	int N = 4;  // 2N intervals on fixed gap
	double xmass[5] = { 1.2, 1.4, 1.6, 1.8, 2 };
	printf("\n%lf", Simpson(xmass, N));
#endif // ITERATIONS
}

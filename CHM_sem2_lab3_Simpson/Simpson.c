#include<stdio.h>
#include <math.h>
#pragma warning(disable:4996)
//#define ITERATIONS 3;

// funtion
double f(double x) {
	return fabs(x*x*x - exp(x) + 1);
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

double CountH(double a, double b) {
	return (b - a) / 2;
}


double Simpson(double* xmass) {
	double h = CountH(xmass[0], xmass[2]);
	double res = (h / 3) * (f(xmass[0]) + 4 * f(xmass[1]) + f(xmass[2]));
	return res;
}

int main() {
#ifdef ITERATIONS

#endif // ITERATIONS

#ifndef ITERATIONS
	double xmass[3] = { -0.2, 0, 0.2 };
	printf("%lf, %lf, %lf", f(xmass[0]), f(xmass[1]), f(xmass[2]));
	printf("\n%lf", Simpson(xmass));
#endif // ITERATIONS



}

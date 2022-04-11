#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#pragma warning(disable:4996)

long double f(long double x) {
	return x * x * x - exp(x) + 1;
}

long double g(long double x) {
	return fabs(x * x * x - exp(x) + 1);
}

void FillFile(FILE* fpdata, long double* mass, int size) {
	for (int i = 0; i < size; i++) {
		fprintf(fpdata, "%.20lf ", (mass[i]));
	}
	fprintf(fpdata, "\n");
}

void PrintMass(long double* mass, int size) {
	for (int i = 0; i < size; i++) {
		printf("%lf ", mass[i]);
	}
	printf("\n");
}

void ProjectX(long double* xMass, long double a, long double b, long double* newMass) {
	long double x = xMass[1];
	newMass[0] = a;
	newMass[2] = b;
	newMass[1] = (b + a) / 2 + (b - a) * x / 2;
}

void ProjectKoeff(long double* koeffMass, long double a, long double b, long double* newMass) {
	long double tmp = (b - a) / 2;
	for (int i = 0; i < 3; i++) {
		newMass[i] = koeffMass[i] * tmp;
	}
}

long double CountH(long double a, long double b, int n) {
	return (b - a) / n;
}

long double LobattoFormula(long double* xMass, long double* koeffMass, long double (*f)(long double)) {
	return koeffMass[0] * f(xMass[0]) + koeffMass[1] * f(xMass[1]) + koeffMass[2] * f(xMass[2]);
}

long double Lobatto(long double* xMass, long double* koeffMass, long double (*f)(long double), long double eps , int* N, long double* H) {
	long double m = 2.5;
	int n = 2;
	long double h = CountH(xMass[0], xMass[2], n);
	long double* xMassNext = (long double*)malloc(sizeof(long double) * 3);
	long double* koeffMassNext = (long double*)malloc(sizeof(long double) * 3);
	if (!xMassNext || !koeffMassNext) {
		printf("error in memory");
	}
	long double I1 = 0, I2 = 0, Itmp = 0;
	I2 = LobattoFormula(xMass, koeffMass, f);	
	do {
		I1 = I2;
		I2 = 0;
		n *= 2;
		h = CountH(xMass[0], xMass[2], n);
		for (int i = 0, k = 0; i < n / 2; i++, k += 2) {
			ProjectX(xMass, xMass[0] + k * h, xMass[0] + (k + 2) * h, xMassNext );
			ProjectKoeff(koeffMass, xMass[0] + k * h, xMass[0] + (k + 2) * h, koeffMassNext);
			I2 += LobattoFormula(xMassNext, koeffMassNext, f);
		}
	} while (fabs(I1 - I2) / (powl(2, m) - 1) > eps);
	*N = n / 2;
	*H = h;
 	return I2;
}

int main() {
	char* FileNameSmooth = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab4_DataSmooth.csv";
	char* FileNameBroken = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab4_DataBroken.csv";
	FILE* dataFile = fopen(FileNameSmooth, "a");

	long double koeffMass[3] = { 1.0/3.0, 4.0/3.0, 1.0/3.0 };
	long double xMass[3] = { -1.0, 0.0, 1.0 };
	long double eps = 1e-1;
	int amount = 10;
	long double res;
	int N = 0;
	long double H = 0;
 	long double* Nmass = (long double*)malloc(sizeof( long double) * amount);
	long double* Hmass = (long double*)malloc(sizeof(long double) * amount);
	long double* resmass = (long double*)malloc(sizeof(long double) * amount);
	for (int k = 0; k < amount; k++) {
		res = Lobatto(xMass, koeffMass, f, eps, &N, &H);
		//printf("\n%lf\n", res);
		resmass[k] = res;
		Nmass[k] = N;
		Hmass[k] = H;
		eps = eps / 10;
	}
	FillFile(dataFile, resmass, amount);
	FillFile(dataFile, Nmass, amount);
	FillFile(dataFile, Hmass, amount);
	fclose(dataFile);
}
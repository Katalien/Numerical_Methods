#include <stdio.h>
#include <math.h>
#pragma warning(disable:4996)
//H(x) = SUMM (yi* fi(x) + y'i*gi(x))
# define ITERATIONS

void FillData(FILE* fpdata, double* xmass, double* ymass, double* ydermass, int amount) {
	for (int i = 0; i < amount; i++) {
		fscanf(fpdata, "%lf", &(xmass[i]));
	}
	for (int i = 0; i < amount; i++) {
		fscanf(fpdata, "%lf", &(ymass[i]));
	}
	for (int i = 0; i < amount; i++) {
		fscanf(fpdata, "%lf", &(ydermass[i]));
	}
}

void ReadFile(FILE* fpdata, double* mass, int amount) {
	for (int i = 0; i < amount; i++) {
		fscanf(fpdata, "%lf", &(mass[i]));
	}
}

void FillFile(FILE* fpdata, double* mass, int size) {
	for (int i = 0; i < size; i++) {
		fprintf(fpdata, "%.15lf ", (mass[i]));
	}
	fprintf(fpdata, "\n");
}

double H(double x, double* xmass, double* ymass, double* ydermass, int size) {
	double res = 0, sum = 0, mult = 1;
	double* A = (double*)malloc(sizeof(double) * size);
	double* B = (double*)malloc(sizeof(double) * size);

	for (int i = 0; i < size; i++) {
		// mult
		for (int k = 0; k < size; k++) {
			if (xmass[i] != xmass[k]) {
				double temp = (x - xmass[k]) / (xmass[i] - xmass[k]) * (x - xmass[k]) / (xmass[i] - xmass[k]);
				mult = mult * temp;
			}
		}

		//summ
		for (int k = 0; k < size; k++) {
			if (xmass[i] != xmass[k]) {
				double temp = (x - xmass[i]) / (xmass[i] - xmass[k]);
				sum = sum + temp;
			}
		}

		double add = (1 - 2 * sum) * mult * ymass[i] + ((x - xmass[i]) * mult * ydermass[i]);
		//res += ((x - points[j]) * yder[j] + (1 - 2 * sum) * ymass[j]) * prod;
		res = res + add;
		mult = 1;
		sum = 0;
	}
	return res;
}

void CountValues(double* xmass, double* ymass, double* ydermass, double* points, double* y, int amount, int size) {
	for (int i = 0; i < amount; i++) {
		y[i] = H(points[i], xmass, ymass, ydermass, size);
	}
}

void PrintMass(double* mass, int size) {
	for (int i = 0; i < size; i++) {
		printf("%lf ", mass[i]);
	}
	printf("\n");
}

int main() {
#ifdef ITERATIONS
	int iter = 13;
	int n = 3;
	int amount = 100;
	char* filePoints = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab1_points.csv";
	//char* fileTable = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab1_compareTableFunc.csv";	
	//char* fileValues = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab1_compareTableValues.csv";
	char* fileTable = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab1_compareModTableFunc.csv";
	char* fileValues = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab1_compareModuleValues.csv";

	FILE* tableFile = fopen(fileTable, "r");
	FILE* pointsFile = fopen(filePoints, "r");
	FILE* valuesFile = fopen(fileValues, "a");
	if (tableFile == NULL || pointsFile == NULL || valuesFile == NULL) {
		printf("\nerror in file\n");
		return -1;
	}
	double* points = (double*)malloc(sizeof(double) * amount);
	ReadFile(pointsFile, points, amount);

	for (int i = 0; i < iter; i++) {
		int size = n + 1;
		//int size = n;
		double* xmass = (double*)malloc(sizeof(double) * size);
		double* ymass = (double*)malloc(sizeof(double) * size);
		double* ydermass = (double*)malloc(sizeof(double) * size);
		double* values = (double*)malloc(sizeof(double) * amount);
		//FillData(tableFile, xmass, ymass, ydermass, size);

		/*printf("n = %d\n", n);
		PrintMass(xmass, n);
		PrintMass(ymass, n);
		PrintMass(ydermass, n);
		printf("\n\n");*/

		CountValues(xmass, ymass, ydermass, points, values, amount, size);

		printf("n = %d\n", n);
		PrintMass(values, amount);

		FillFile(valuesFile, values, amount);
		free(xmass);
		free(ymass);
		free(ydermass);
		free(values);
		n++;
	}
	free(points);
	fclose(tableFile);
	fclose(pointsFile);
	fclose(valuesFile);
#endif
#ifndef ITERATIONS
		int n = 6;
	int amount = 100;
	char* fileTable = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab1_compareTableFunc.csv";
	char* filePoints = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab1_points.csv";
	char* fileValues = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab1_compareTableFunc.csv";
	FILE* tableFile = fopen(fileTable, "r");
	FILE* pointsFile = fopen(filePoints, "r");
	FILE* valuesFile = fopen(fileValues, "a");
	if (tableFile == NULL || pointsFile == NULL || valuesFile == NULL) {
		printf("\nerror in file\n");
		return -1;
	}
	double* points = (double*)malloc(sizeof(double) * amount);
	ReadFile(pointsFile, points, amount);

	int size = n;
	double* xmass = (double*)malloc(sizeof(double) * size);
	double* ymass = (double*)malloc(sizeof(double) * size);
	double* ydermass = (double*)malloc(sizeof(double) * size);
	double* values = (double*)malloc(sizeof(double) * amount);
	FillData(tableFile, xmass, ymass, ydermass, size);
	//
	/*printf("n = %d\n", n);
	PrintMass(xmass, n);
	PrintMass(ymass, n);
	PrintMass(ydermass, n);*/
	printf("\n\n");
	PrintMass(ydermass, n);
	printf("\n\n");
	CountValues(xmass, ymass, ydermass, points, values, amount, size);
	//
	printf("n = %d\n", n);
	PrintMass(values, amount);
	//
	FillFile(valuesFile, values, amount);

	free(xmass);
	free(ymass);
	free(ydermass);
	free(values);
	free(points);
	fclose(tableFile);
	fclose(pointsFile);
#endif
}

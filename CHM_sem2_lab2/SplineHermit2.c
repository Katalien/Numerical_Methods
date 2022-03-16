#include <stdio.h>
#include <math.h>
#pragma warning(disable:4996)
#define ITERATIONS 3;

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

void FillDataWithDer(FILE* fpdata, double* xmass, double* ymass, double* ydermass, int amount) {
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

// Count additional values
double Abs(double x) {
	if (x < 0) {
		return x * (-1);
	}
	else {
		return x;
	}
}

void CountDividedDiff(double* xmass, double* ymass, double* divDif, int amount) {
	for (int i = 0; i < amount; i++) {
		divDif[i] = (ymass[i + 1] - ymass[i]) / (xmass[i + 1] - xmass[i]);
	}
}

void CountEndDerivative( double* xmass, double* ymass, double* der1, double* der2, int amount) {
	double h1, h2, delta1, delta2;
	// first der:
	int i = 0;
	h1 = xmass[i + 1] - xmass[i];
	h2 = xmass[i + 2] - xmass[i + 1];
	delta1 = (ymass[i + 1] - ymass[i]) / (xmass[i + 1] - xmass[i]);
	delta2 = (ymass[i + 2] - ymass[i+1]) / (xmass[i + 2] - xmass[i+1]);
	*der1 = ((2 * h1 + h2) * delta1 - h1 * delta2) / (h1 + h2);
	if ((*der1) * delta1 < 0) {
		*der1 = 0;
	}
	else if (delta1 * delta2 < 0 && Abs(*der1) > Abs(3 * delta1)) {
		*der1 = 3 * delta1;
	}
	//last der:
	i = amount - 2;
	h1 = xmass[i + 1] - xmass[i];
	h2 = xmass[i ] - xmass[i - 1];
	delta1 = (ymass[i + 1] - ymass[i]) / (xmass[i + 1] - xmass[i]);
	delta2 = (ymass[i ] - ymass[i -1]) / (xmass[i ] - xmass[i -1]);
	*der2 = ((2 * h1 + h2) * delta1 - h1 * delta2) / (h1 + h2);
	if ((*der2) * delta1 < 0) {
		*der2 = 0;
	}
	else if (delta1 * delta2 < 0 && Abs(*der2) > Abs(3 * delta1)) {
		*der2 = 3 * delta1;
	}
}

void CountInteriorDerivative(double* xmass, double* ymass, double* ydermass, int amount) {
	double tmp = 0;
	double* divDiff = (double*)malloc((amount - 1) * sizeof(double));   // array n-2 points
	CountDividedDiff(xmass, ymass, divDiff, amount - 1);
	// count interior points
	for (int i = 1; i < amount - 1; i++) {   // fill from 1 to n-1
		if (divDiff[i - 1] * divDiff[i] <= 0) {
			ydermass[i] = 0;
		}
		
		else if(divDiff[i - 1] * divDiff[i] > 0 ) {   // ???
			tmp = 0.5*(1 / divDiff[i - 1] + 1 / divDiff[i]);
			ydermass[i] = 1 / tmp;
		}
	}
}

// Hermit function
double Multiply(double x, double* xmass, int i, int size) {
	double res = 1, temp = 0;
	for (int j = 0; j < size; j++) {
		if (xmass[i] != xmass[j]) {
			temp = (x - xmass[j]) / (xmass[i] - xmass[j]);
			res = res * temp * temp;
		}
	}
	return res;
}


double Summ(double x, double* xmass, int i, int size) {
	double res = 0;
	for (int k = 0; k < size; k++) {
		if (xmass[i] != xmass[k]) {
			double tmp = (x - xmass[i]) / (xmass[i] - xmass[k]);
			res = res + tmp;
		}
	}
	return res;
}

double HermitSpline(double x, double* xmass, double* ymass, double* ydermass, int size) {
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
		res = res + add;
		mult = 1;
		sum = 0;
	}
	return res;
}

void FillAdditionalArray(double* xmass, double* ymass, double* ydermass, double* addxmass, double* addymass, double* addydermass, int from) {
		addxmass[0] = xmass[from];
		addymass[0] = ymass[from];
		addydermass[0] = ydermass[from];
		addxmass[1] = xmass[from + 1];
		addymass[1] = ymass[from + 1];
		addydermass[1] = ydermass[from + 1];
}

void CountValues(double* xmass, double* ymass, double* ydermass, double* points, double* y, int amount, int n) {
	int i = 0;
	int j = 0;
	double curxmass[2], curymass[2], curydermass[2];
	for (j; j < n-1; j++) {    // цикл по табличным точкам
		FillAdditionalArray(xmass, ymass, ydermass, curxmass, curymass, curydermass, j);
		for (i; i < amount; i++) {    //цикл по 100 точкам
			if (points[i] <= xmass[j+1]) {
				y[i] = HermitSpline(points[i], curxmass, curymass, curydermass, 2);
			}
			else {
				break;
			}

		}
	}
}

int main() {
#ifdef ITERATIONS
		int iter = 13;
		int n = 4;
		int amount = 100;
		char* filePoints = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab2_points.csv";
		//char* fileTable = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab2_graphTableFunc.csv";
		//char* fileValues = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab2_graphValues.csv";
		char* fileTable = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab2_WithDerTableFunc.csv";
		char* fileValues = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab2_WithDerValues.csv";
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
			double* xmass = (double*)malloc(sizeof(double) * n);
			double* ymass = (double*)malloc(sizeof(double) * n);
			double* ydermass = (double*)malloc(sizeof(double) * n);
			double* values = (double*)malloc(sizeof(double) * amount);
			//FillData(tableFile, xmass, ymass, n);
			/*CountInteriorDerivative(xmass, ymass, ydermass, n);
			CountEndDerivative(xmass, ymass, &(ydermass[0]), &(ydermass[n - 1]), n);*/
			FillDataWithDer(tableFile, xmass, ymass, ydermass, n);
			PrintMass(xmass, n);
			printf("\n");
			PrintMass(ymass, n);
			printf("\n");
			PrintMass(ydermass, n);
			printf("\n");
			CountValues(xmass, ymass, ydermass, points, values, amount, n);
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
	//int amount = 100;
	//int n = 16;
	//double* xmass = (double*)malloc(sizeof(double) * n);
	//double* ymass = (double*)malloc(sizeof(double) * n);
	//double* ydermass = (double*)malloc(sizeof(double) * n);
	//double* points = (double*)malloc(sizeof(double) * amount);
	//double* values = (double*)malloc(sizeof(double) * amount);
	//char* fileTable = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab2_tableFunc.csv";
	//char* filePoints = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab2_points.csv";
	//char* fileValues = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab2_values.csv";
	//FILE* tableFile = fopen(fileTable, "r");
	//FILE* pointsFile = fopen(filePoints, "r");
	//FILE* valuesFile = fopen(fileValues, "a");
	//if (tableFile == NULL || pointsFile == NULL || valuesFile == NULL) {
	//	printf("\nerror in file\n");
	//	return -1;
	//}
	//FillData(tableFile, xmass, ymass, n);
	//ReadFile(pointsFile, points, amount);
	//CountInteriorDerivative(xmass, ymass, ydermass, n);
	//CountEndDerivative(xmass, ymass, &(ydermass[0]), &(ydermass[n - 1]), n);
	//CountValues(xmass, ymass, ydermass, points, values, amount, n);
	////FillFile(valuesFile, values, amount);
	//fclose(tableFile);
	//fclose(pointsFile);
	//fclose(valuesFile);
	
#endif
}

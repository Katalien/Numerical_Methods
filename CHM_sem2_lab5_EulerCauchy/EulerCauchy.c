#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#pragma warning(disable:4996)
//#define ITERATIONS 1
#//define TEST 2
//#define N_EPS 3
#define TEST_OTCHET 2

long double f(long double x, long double y) {
	return powl(y, 2) * exp(x) - 2 * y;
}

long double fExact(long double x) {
	return exp(-x);
}

void FillFile(FILE* fpdata, long double* mass, int size) {
	for (int i = 0; i < size; i++) {
		fprintf(fpdata, "%.20lf ", (mass[i]));
	}
	fprintf(fpdata, "\n");
}

void PrintMass(long double* mass, int size) {
	for (int i = 0; i < size; i++) {
		printf("%.15lf ", mass[i]);
	}
	printf("\n");
}

void CopyMass(long double* from, long double* to, int n) {
	for (int i = 0; i < n; i++) {
		to[i] = from[i];
	}
}

long double CountH(long double a, long double b, int n) {
	return (long double)((b - a) / n);
}

long double NextX_(long double a, long double h, int k) {
	return a + k * h;
}

long double NextX(long double x, long double h) {
	return x + h;
}

// f: counter + 2
long double CountNextY(long double x, long double y, long double h, long double (*f)(long double, long double)) {
	long double xNext = NextX(x, h);
	long double f_val_i = f(x, y);
	long double yWaveNext = y + h * f_val_i;
	long double y_next = y + h * (f_val_i + f(xNext, yWaveNext)) / 2;
	return y_next;
}

long double* FindLocalErr(long double a, long double y_a, long double b, long double eps, long double (*f)(long double, long double), int amount) {
	int n = 2;
	long double h_max = CountH(a, b, n), h = h_max;;
	long double cur_x = a, x_h, x_h_div;
	long double cur_y = y_a, y_h = 0, y_h_div = 0;
	long double* errMass = (long double*)malloc(sizeof(long double) * amount);
	// считаем локальную погрешность
	for (int i = 0; i < amount; i++) {
		x_h = NextX(a, h);
		y_h = CountNextY(a, cur_y, h, f);
		h = h / 2;
		if (h < 1e-18) {
			printf("ooops, h is tooo smalll\n");
			break;
		}
		printf("\n%.55lf", x_h);
		errMass[i] = fabs(fExact(x_h) - y_h);
	}
	printf("\n\n\n");
	return errMass;
}

long double* FindGlobalErr(long double a, long double y_a, long double b, long double eps, long double (*f)(long double, long double), int amount) {
	int n = 2;
	long double h_max = CountH(a, b, n), h = h_max;;
	long double cur_x = a, x_h, x_h_div, x, next_x, tmp_x;
	long double cur_y = y_a, y_h = 0, y_h_div = 0, tmp_y;
	long double* errMass = (long double*)malloc(sizeof(long double) * amount);


	// считаем глабольную погрешность
	for (int i = 0; i < amount; i++) {
		cur_x = a;
		cur_y = y_a;
		do {
			tmp_x = NextX(cur_x, h);
			tmp_y = CountNextY(cur_x, cur_y, h, f);
			cur_x = tmp_x;
			cur_y = tmp_y;
		} while (cur_x != b - h);
		x_h = NextX(cur_x, h);
		y_h = CountNextY(cur_x, cur_y, h, f);
		h = h / 2;
		printf("\n%lf, %lf", x_h, y_h);
		errMass[i] = fabs(fExact(b) - y_h);
	}
	return errMass;
}

// p = 2!!!!!
int FindSolution(long double a, long double y_a, long double b, long double eps, long double (*f)(long double, long double), long double** xMass, long double** yMass) {
	int p = 2;
	int counter = 0; // ampont of points
	int iterations = 0;  // amount of f call
	long double x_next, x_cur = a, x_middle, x_div_next;
	long double y_next, y_cur = y_a, y_middle, y_div_next;
	long double H = CountH(a, b, 1);
	long double h = H;
	long double* bufferX = (long double*)malloc(sizeof(long double) * 140000000);
	long double* bufferY = (long double*)malloc(sizeof(long double) * 140000000);
	bufferY[counter] = y_a;
	bufferX[counter] = a;
	do {
		do {
			h /= 2;
			x_next = NextX(x_cur, h);
			y_next = CountNextY(x_cur, y_cur, h, f);

			x_middle = NextX(x_cur, h / 2);
			y_middle = CountNextY(x_cur, y_cur, h / 2, f);

			x_div_next = NextX(x_middle, h / 2);
			y_div_next = CountNextY(x_middle, y_middle, h / 2, f);

			iterations += 2;

		} while (fabs(y_next - y_div_next) / (powl(2, p) - 1) >= eps);
		printf("%.8lf %.8lf %.8lf\n", x_div_next, y_div_next, fExact(x_div_next));
		x_cur = x_next;
		y_cur = y_div_next;   ///
		h *= 2;
		counter++;
		bufferY[counter] = y_cur;
		bufferX[counter] = x_cur;

	} while (x_cur != b);
	long double* xRes = (long double*)malloc(sizeof(long double) * counter);
	long double* yRes = (long double*)malloc(sizeof(long double) * counter);

	CopyMass(bufferX, xRes, counter);
	CopyMass(bufferY, yRes, counter);
	*xMass = xRes;
	*yMass = yRes;
	free(bufferX);
	free(bufferY);
	return iterations;
}


int main() {
#ifdef ITERATIONS
	char* FileNameError = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_GlobalLocalErrEps.csv";
	char* FileNameData = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_Data.csv";
	FILE* dataFile = fopen(FileNameError, "a");

	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	long double eps = 1e-1;
	int amount = 15;
	int size = 0;
	for (int i = 0; i < amount; i++) {
		long double* bufferY = (long double*)malloc(sizeof(long double) * 140000);
		long double* bufferX = (long double*)malloc(sizeof(long double) * 140000);
		long double* resGlobal = (long double*)malloc(sizeof(long double) * amount);
		size = FindSolution(a, y_a, b, eps, f, bufferX, bufferY);
		long double* yMass = (long double*)malloc(sizeof(long double) * size);
		long double* xMass = (long double*)malloc(sizeof(long double) * size);
		CopyMass(bufferX, xMass, size);
		CopyMass(bufferY, yMass, size);
		fprintf(dataFile, "%d ", size);
		FillFile(dataFile, xMass, size);
		FillFile(dataFile, yMass, size);
		free(bufferX);
		free(bufferY);
		free(xMass);
		free(yMass);
	}
#endif ITERATIONS

#ifdef N_EPS
	char* FileName = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_err_eps.csv";
	FILE* dataFile = fopen(FileName, "a");
	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	long double eps = 1e-1;
	long double yGlobal = fExact(1);
	long double yLocal = fExact(-1);
	long double epsMass[3] = { 1e-3, 1e-7, 1e-15 };
	int amount = 15;
	int size;
	long double* globalY = (long double*)malloc(sizeof(long double) * amount);
	long double* localY = (long double*)malloc(sizeof(long double) * amount);
	long double* globalErr = (long double*)malloc(sizeof(long double) * amount);
	long double* localErr = (long double*)malloc(sizeof(long double) * amount);
	int* sizeMass = (int*)malloc(sizeof(int) * amount);
	for (int i = 0; i < amount; i++) {
		long double* bufferY = (long double*)malloc(sizeof(long double) * 140000);
		long double* bufferX = (long double*)malloc(sizeof(long double) * 140000);
		size = FindSolution(a, y_a, b, eps, f, bufferX, bufferY);
		long double* yMass = (long double*)malloc(sizeof(long double) * size);
		long double* xMass = (long double*)malloc(sizeof(long double) * size);
		CopyMass(bufferX, xMass, size);
		CopyMass(bufferY, yMass, size);
		localY[i] = yMass[0];
		globalY[i] = yMass[size - 1];
		localErr[i] = fabs(yMass[1] - fExact(xMass[1]));
		globalErr[i] = fabs(yMass[size - 1] - fExact(1));
		printf("\n");
		printf("%.7lf %.7lf\n", localErr[i], globalErr[i]);
		free(bufferX);
		free(bufferY);
		free(xMass);
		free(yMass);
		eps /= 10;
	}
	//PrintMass(globalErr, amount);
	FillFile(dataFile, localErr, amount);
	FillFile(dataFile, globalErr, amount);
	//fclose(dataFile);
#endif N_EPS

#ifdef TEST
	char* FileNameError = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_eps1.csv";
	char* FileNameData = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_eps3.csv";
	FILE* dataFile = fopen(FileNameError, "a");

	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	int n = 2;
	long double eps = 1e-15;
	long double* bufferY = (long double*)malloc(sizeof(long double) * 140000);
	long double* bufferX = (long double*)malloc(sizeof(long double) * 140000);
	long double* localErrMass = NULL;
	long double* globalErrMass = NULL;
	int amount = 20;
	int size = 0;
	size = FindSolution(a, y_a, b, 1e-1, f, bufferX, bufferY);
	printf("%d\n", size);

	long double* yMass = (long double*)malloc(sizeof(long double) * size);
	long double* xMass = (long double*)malloc(sizeof(long double) * size);
	CopyMass(bufferX, xMass, size);
	CopyMass(bufferY, yMass, size);
	PrintMass(xMass, size);
	FillFile(dataFile, xMass, size);
	FillFile(dataFile, yMass, size);



	//size = FindSolution(a, y_a, b, eps, f, bufferX, bufferY);
	//long double* yMass = (long double*)malloc(sizeof(long double) * size);
	//long double* xMass = (long double*)malloc(sizeof(long double) * size);
	//CopyMass(bufferX, xMass, size);
	//CopyMass(bufferY, yMass, size);
	//printf("\n%d ", size);

	/*FillFile(dataFile, xMass, size);
	FillFile(dataFile, yMass, size);
	free(bufferX);
	free(bufferY);
	free(xMass);
	free(yMass);*/
	fclose(dataFile);
#endif // TEST
#ifdef TEST_OTCHET

	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	int n = 2;
	long double eps = 1e-1;
	long double* bufferY = (long double*)malloc(sizeof(long double) * 140000);
	long double* bufferX = (long double*)malloc(sizeof(long double) * 140000);
	long double* localErrMass = NULL;
	long double* globalErrMass = NULL;
	long double* xMass = NULL;
	long double* yMass = NULL;
	int amount = 20;
	int size = 0;
	size = FindSolution(a, y_a, b, eps, f, &xMass, &yMass);
	printf("\n%d\n", size);
#endif
}


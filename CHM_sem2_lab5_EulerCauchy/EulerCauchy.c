#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#pragma warning(disable:4996)
//#define ITERATIONS 1
#define TEST 2
//#define N_EPS 3

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
		printf("%lf ", mass[i]);
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
	return a + k*h;
}

long double NextX(long double x, long double h) {
	return x + h;
}

long double CountWaveY(long double x, long double y, long double h, long double (*f)(long double, long double)) {
	long double tmp = y + h * f(x, y);
	return tmp;
}

long double CountNextY(long double x, long double y, long double h,long double (*f)(long double, long double)) {
	long double xNext = NextX(x, h);
	long double yWaveNext = CountWaveY(x, y, h, f);
	long double tmp = y + h * (f(x, y) + f(xNext, yWaveNext)) / 2;
	return y + h * (f(x, y) + f(xNext, yWaveNext)) / 2;
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
		x_h_div = NextX(cur_x, h / 2);
		y_h = CountNextY(a, cur_y, h, f);
		y_h = CountNextY(cur_x, cur_y, h, f);
		y_h_div = CountNextY(cur_x, cur_y, h / 2, f);
		h = h / 2;
		if (h < 1e-18) {
			printf("ooops, h is tooo smalll\n");
			break;
		}
		printf("%d: %.15lf, %.20lf\n", i, h, fabs(fExact(a) - y_h));
		errMass[i] = fabs(fExact(a) - y_h);
	}
	return errMass;
}

long double* FindGlobalErr(long double a, long double y_a, long double b, long double eps, long double (*f)(long double, long double), int amount) {
	int n = 2;
	long double h_max = CountH(a, b, n), h = h_max;;
	long double cur_x = a, x_h, x_h_div, x, next_x;
	long double cur_y = y_a, y_h = 0, y_h_div = 0;
	long double* errMass = (long double*)malloc(sizeof(long double) * amount);
	
	
	// считаем глабольную погрешность
	for (int i = 0; i < amount; i++) {
		do {
			cur_x = NextX(cur_x, h);
		} while (cur_x != b-h);		
		x_h_div = NextX(cur_x, h / 2);
		y_h = CountNextY(cur_x, cur_y, h, f);
		y_h = CountNextY(cur_x, cur_y, h, f);
		y_h_div = CountNextY(cur_x, cur_y, h / 2, f);
		h = h / 2;

		if (h < 1e-18) {
			printf("ooops, h is tooo smalll\n");
			break;
		}
		printf("%d: %.15lf, %.20lf\n", i, h, fabs(fExact(a) - y_h));
		errMass[i] = fabs(fExact(a) - y_h);
	}
	return errMass;
}

//int FindSolution(long double a, long double y_a, long double b, long double eps, long double (*f)(long double, long double), long double* xMass, long double* yMass) {
//	int p = 1;
//	int counter = 0;
//	long double x_next, x_cur = a, x_middle, x_div_next;
//	long double y_next, y_cur = y_a, y_middle, y_div_next;
//	long double H = CountH(a, b, 1);
//	long double h = H ;
//	long double* bufferY = (long double*)malloc(sizeof(long double) * 140000);
//	long double* bufferX = (long double*)malloc(sizeof(long double) * 140000);
//	bufferY[counter] = y_a;
//	bufferX[counter] = a;
//	do {
//		do {
//			h /= 2;
//			x_next = NextX(x_cur, h);
//			y_next = CountNextY(x_cur, y_cur, h, f);
//			x_middle = NextX(x_cur, h / 2);
//			y_middle = CountNextY(x_cur, y_cur, h / 2, f);
//			x_div_next = NextX(x_middle, h/2);
//			y_div_next = CountNextY(x_middle, y_middle, h / 2, f);
//		} while (fabs(y_next - y_div_next) / (powl(2, p) - 1) > eps);
//		x_cur = x_next;
//		y_cur = y_next;   ///
//		h *= 2;
//		counter++;
//		bufferY[counter] = y_cur;
//		bufferX[counter] = x_cur;
//		//printf("%.10lf\n", fabs(y_cur - fExact(x_cur)));
//	} while (x_cur != b);
//	long double* resMassY = (long double*)malloc(sizeof(long double) * counter);
//	long double* resMassX = (long double*)malloc(sizeof(long double) * counter);
//	CopyMass(bufferY, xMass, counter);
//	CopyMass(bufferX, yMass, counter);
//	printf("%d", counter);
//	free(bufferY);
//	free(bufferX);
//	return counter;
//}

int FindSolution(long double a, long double y_a, long double b, long double eps, long double (*f)(long double, long double), long double* xMass, long double* yMass) {
	int p = 1;
	int counter = 0;
	long double x_next, x_cur = a, x_middle, x_div_next;
	long double y_next, y_cur = y_a, y_middle, y_div_next;
	long double H = CountH(a, b, 1);
	long double h = H;
	//long double* bufferY = (long double*)malloc(sizeof(long double) * 140000);
	//long double* bufferX = (long double*)malloc(sizeof(long double) * 140000);
	yMass[counter] = y_a;
	xMass[counter] = a;
	do {
		do {
			h /= 2;
			x_next = NextX(x_cur, h);
			y_next = CountNextY(x_cur, y_cur, h, f);
			x_middle = NextX(x_cur, h / 2);
			y_middle = CountNextY(x_cur, y_cur, h / 2, f);
			x_div_next = NextX(x_middle, h / 2);
			y_div_next = CountNextY(x_middle, y_middle, h / 2, f);
			printf("\n %.7lf",fabs(1 - y_div_next));
} while (fabs(y_next - y_div_next) / (powl(2, p) - 1) > eps);
x_cur = x_next;
y_cur = y_next;   ///
h *= 2;
counter++;
yMass[counter] = y_cur;
xMass[counter] = x_cur;
//printf("%.10lf\n", fabs(y_cur - fExact(x_cur)));
	} while (x_cur != b);
	return counter;
}


int main() {
#ifdef ITERATIONS
	char* FileNameError = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_LocalGlobalErr.csv";
	char* FileNameData = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_Data.csv";
	FILE* dataFile = fopen(FileNameData, "a");

	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	long double eps = 1e-1;
	int amount = 15;
	int size = 0;
	for (int i = 0; i < amount; i++) {
		long double* bufferY = (long double*)malloc(sizeof(long double) * 140000);
		long double* bufferX = (long double*)malloc(sizeof(long double) * 140000);
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
	char* FileName = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_ChislennoeReshenie_3eps.csv";
	//FILE* dataFile = fopen(FileName, "a");
	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	long double eps;
	long double epsMass[3] = { 1e-3, 1e-7, 1e-15 };
	int amount = 15;
	int size;
	long double ywave =0 + f(0, 0);
	long double y = 0 + 0.5*(f(0.0, 0.0) + f(1, ywave));
	printf("%lf, %lf", ywave, y);
	//int* sizeMass = (int*)malloc(sizeof(int) * amount);
	
		/*eps = epsMass[0];
		long double* bufferY = (long double*)malloc(sizeof(long double) * 140000);
		long double* bufferX = (long double*)malloc(sizeof(long double) * 140000);
		size = FindSolution(a, y_a, b, eps, f, bufferX, bufferY);
		long double* yMass = (long double*)malloc(sizeof(long double) * size);
		long double* xMass = (long double*)malloc(sizeof(long double) * size);
		CopyMass(bufferX, xMass, size);
		CopyMass(bufferY, yMass, size);
		printf("%d\n", size);
		FillFile(dataFile, xMass, size);
		FillFile(dataFile, yMass, size);
		free(bufferX);
		free(bufferY);
		free(xMass);
		free(yMass);		
	
	fclose(dataFile);*/
#endif N_EPS

#ifdef TEST
	char* FileNameError = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_LocalGlobalErr.csv";
	char* FileNameData = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab5_Eps15.csv";
	FILE* dataFile = fopen(FileNameData, "a");

	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	int n = 2;
	long double eps = 1e-5;
	long double* bufferY = (long double*)malloc(sizeof(long double) * 140000);
	long double* bufferX = (long double*)malloc(sizeof(long double) * 140000);
	long double* localErrMass = NULL;
	long double* globalErrMass = NULL;
	int amount = 55;
	int size = 0;
	printf("%.7lf", fExact(0));
	size = FindSolution(a, y_a, b, eps, f, bufferX, bufferY);
	long double* yMass = (long double*)malloc(sizeof(long double) * size);
	long double* xMass = (long double*)malloc(sizeof(long double) * size);
	CopyMass(bufferX, xMass, size);
	CopyMass(bufferY, yMass, size);
	PrintMass(xMass, size);
	printf("\n\n");
	PrintMass(yMass, size);
	printf("\n\n%d", size);
	FillFile(dataFile, xMass, size);
	FillFile(dataFile, yMass, size);
	free(bufferX);
	free(bufferY);
	free(xMass);
	free(yMass);
#endif // TEST
}


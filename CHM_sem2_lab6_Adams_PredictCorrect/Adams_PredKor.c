#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#pragma warning(disable:4996)

//#define ITERATIONS 1
//#define TEST 2
#define TEST_OTCHET

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
		printf("%.20lf ", mass[i]);
	}
	printf("\n");
}

void CopyMass(long double* from, long double* to, int n) {
	for (int i = 0; i < n; i++) {
		to[i] = from[i];
	}
}

long double NextX(long double x, long double h) {
	return x + h;
}

long double CountH(long double a, long double b, int n) {
	return (long double)((b - a) / n);
}

// считаем разгонные точки
long double CountWaveY(long double x, long double y, long double h, long double (*f)(long double, long double)) {
	long double tmp = y + h * f(x, y);
	return tmp;
}

long double CountAddY(long double x, long double y, long double h, long double (*f)(long double, long double)) {
	long double xNext = NextX(x, h);
	long double yWaveNext = CountWaveY(x, y, h, f);
	//printf("ywave = %lf\n\n", yWaveNext);
	long double tmp = y + h * (f(x, y) + f(xNext, yWaveNext)) / 2;
	return y + h * (f(x, y) + f(xNext, yWaveNext)) / 2;
}

long double* FindAddPoints(long double a, long double y_a, long double b, long double h, long double (*f)(long double, long double)) {
	int p = 2;
	long double* points = (long double*)malloc(sizeof(long double) * 2);
	long double x_next, x_cur = a, x_middle, x_div_next;
	long double y_next, y_cur = y_a, y_middle, y_div_next;

	x_next = NextX(x_cur, h);
	y_next = CountAddY(x_cur, y_cur, h, f);

	x_middle = NextX(x_cur, h / 2);
	y_middle = CountAddY(x_cur, y_cur, h / 2, f);

	x_div_next = NextX(x_middle, h / 2);
	y_div_next = CountAddY(x_middle, y_middle, h / 2, f);

	points[0] = x_div_next;
	points[1] = y_div_next;
	return points;
}
// Закончили считать разгонные точки


// f: +3 counter
long double PredictCorrectY(long double x, long double y, long double x_prev, long double y_prev, long double h, long double (*f)(long double, long double)) {
	long double x_next = NextX(x, h);
	long double f_val_i = f(x, y);
	long double f_val_prev = f(x_prev, y_prev);
	long double y_pred = y + h * 0.5 * (3 * f_val_i - f(x_prev, y_prev));
	long double y_corr = y + h * 0.5 * (f(x_next, y_pred) + f_val_i);
	return y_corr;
}

// Закончили уравнения системы

//long double* FindLocalErr(long double a, long double y_a, long double b, long double eps, long double (*f)(long double, long double), int amount, long double* hMass) {
//	long double p = 2;
//	int counter = 1;
//	long double h = CountH(a, b, 2);
//
//	long double x_prev = a, x_cur, x_next, x_middle, x_div_next;
//	long double y_prev = y_a, y_cur, y_next, y_middle, y_div_next;
//	long double* errMass = (long double*)malloc(sizeof(long double) * amount);
//
//	//long double* points = FindAddPoints(a, y_a, b, eps, f);
//
//
//	// считаем локальную погрешность
//	for (int i = 0; i < amount; i++) {
//		x_cur = NextX(x_prev, h);
//		y_cur = CountAddY(x_prev, y_prev, h, f);
//		hMass[i] = h;
//		x_next = NextX(x_cur, h);
//		y_next = PredictCorrectY(x_cur, y_cur, x_prev, y_prev, h, f);
//		//printf("\n%d %.10lf ", i, h);
//		errMass[i] = fabs(fExact(x_next) - y_next);
//		h = h / 2;
//	}
//	printf("\n\n\n");
//	return errMass;
//}
//
//
//long double* FindGlobalErr(long double a, long double y_a, long double b, long double eps, long double (*f)(long double, long double), int amount, long double* hMass) {
//	long double p = 2;
//	int counter = 1;
//	long double h = CountH(a, b, 2);
//	long double* errMass = (long double*)malloc(sizeof(long double) * amount);
//
//
//	long double x_prev = a, x_cur, x_next, x_next_tmp, x_h;
//	long double y_prev = y_a, y_cur, y_next, y_next_tmp, y_h;
//	int n = 2;
//	//long double* points = FindAddPoints(a, y_a, b, eps, f);
//	x_cur = NextX(x_prev, h);
//	y_cur = CountAddY(x_prev, y_prev, h, f);
//	//x_cur = points[0];
//	//y_cur = points[1];
//
//
//	// считаем глбольную погрешность
//	for (int i = 0; i < amount; i++) {
//		h /= 2;
//		x_prev = a;
//		y_prev = y_a;
//		x_cur = NextX(x_prev, h);
//		y_cur = CountAddY(x_prev, y_prev, h, f);
//
//		printf("%d %.10lf  \n", i, h);
//		hMass[i] = h;
//		do {
//			x_next_tmp = NextX(x_cur, h);
//			y_next_tmp = PredictCorrectY(x_cur, y_cur, x_prev, y_prev, h, f);
//			x_prev = x_cur;
//			y_prev = y_cur;
//			x_cur = x_next_tmp;
//			y_cur = y_next_tmp;
//		} while (x_cur != b - h);
//		x_h = NextX(x_cur, h);
//		y_h = PredictCorrectY(x_cur, y_cur, x_prev, y_prev, h, f);
//		//printf("\n%d %.10lf ",i, h);
//		errMass[i] = fabs(fExact(b) - y_h);
//	}
//	return errMass;
//}


long double* FindSolution(long double a, long double b, long double y_a, long double eps, long double (*f)(long double, long double), long double** xMass, long double** yMass) {
	long double p = 2;
	int counter = 1;
	int iterations = 0;
	int iterations_test = 0;
	long double h = CountH(a, b, 1);

	long double* bufferX = (long double*)malloc(sizeof(long double) * 140000000);
	long double* bufferY = (long double*)malloc(sizeof(long double) * 140000000);

	long double x_prev = a, x_cur, x_next, x_middle, x_div_next;
	long double y_prev = y_a, y_cur, y_next, y_middle, y_div_next;

	do {
		h /= 2;
		long double* points = FindAddPoints(a, y_a, b, h, f);
		x_cur = points[0];
		y_cur = points[1];
		free(points);

		x_next = NextX(x_cur, h);
		y_next = PredictCorrectY(x_cur, y_cur, x_prev, y_prev, h, f);

		x_middle = NextX(x_cur, h / 2);
		y_middle = PredictCorrectY(x_cur, y_cur, x_prev, y_prev, h / 2, f);


		x_div_next = NextX(x_middle, h / 2);
		y_div_next = PredictCorrectY(x_middle, y_middle, x_cur, y_cur, h / 2, f);

		iterations += 3;

	} while (fabs(y_next - y_div_next) / (powl(2, p) - 1) >= eps);
	bufferX[0] = a;
	bufferY[0] = y_a;
	bufferX[1] = x_cur;
	bufferY[1] = y_cur;
	counter++;
	iterations_test += 3;
	bufferX[counter] = x_div_next;
	bufferY[counter] = y_div_next;

	x_prev = x_cur;
	y_prev = y_cur;
	x_cur = x_div_next;
	y_cur = y_div_next;

	do {
		do {

			x_next = NextX(x_cur, h);
			y_next = PredictCorrectY(x_cur, y_cur, x_prev, y_prev, h, f);

			x_middle = NextX(x_cur, h / 2);
			y_middle = PredictCorrectY(x_cur, y_cur, x_prev, y_prev, h / 2, f);


			x_div_next = NextX(x_middle, h / 2);
			y_div_next = PredictCorrectY(x_middle, y_middle, x_cur, y_cur, h / 2, f);

			iterations += 3;
			
			h /= 2;
			
		} while (fabs(y_next - y_div_next) / (powl(2, p) - 1) >= eps);
		x_prev = x_cur;
		y_prev = y_cur;
		x_cur = x_div_next;
		y_cur = y_div_next;
		bufferX[counter] = x_cur;
		bufferY[counter] = y_cur;
		counter++;
		printf("%d %d\n", iterations_test, iterations);
		h *= 2;
	} while (x_cur != b);

	long double* xRes = (long double*)malloc(sizeof(long double) * counter);
	long double* yRes = (long double*)malloc(sizeof(long double) * counter);

	CopyMass(bufferX, xRes, counter);
	CopyMass(bufferY, yRes, counter);
	*xMass = xRes;
	*yMass = yRes;
	free(bufferX);
	free(bufferY);
	printf("\n\n%d %d\n", iterations_test, iterations);
	return iterations;
}


int main() {
#ifdef TEST

	char* FileNameData = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab6_localH_again_.csv";
	FILE* dataFile = fopen(FileNameData, "a");
	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	long double eps = 1e-3;
	long double* xMass = NULL;
	long double* yMass = NULL;
	int amount = 20;
	long double* localH = (long double*)malloc(sizeof(long double) * amount);
	long double* globalH = (long double*)malloc(sizeof(long double) * amount);
	long double* localErr = FindLocalErr(a, y_a, b, eps, f, amount, localH);
	// long double* globalErr = FindGlobalErr(a, y_a, b, eps, f, amount, globalH);




	FillFile(dataFile, localH, amount);
	FillFile(dataFile, localErr, amount);



	free(xMass);
	free(yMass);
	fclose(dataFile);
	return 0;

#endif // 
#ifdef ITERATIONS

	char* FileNameData = "C:/Users/z.kate/Desktop/4 семестр/CHM/data/lab6_eps_err_global.csv";
	FILE* dataFile = fopen(FileNameData, "a");
	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	long double eps = 1e-1;
	printf("%d \n", 3);
	int amount = 25;
	int size = 0;
	long double* xMass = NULL;
	long double* yMass = NULL;
	long double* err = (long double*)malloc(sizeof(long double) * 15);
	for (int i = 0; i < 15; i++) {

		size = FindSolution(a, b, y_a, eps, f, &xMass, &yMass);
		printf("%d \n", i);
		err[i] = fabs(fExact(1) - yMass[size - 1]);
		eps /= 10;
	}
	FillFile(dataFile, err, 15);
	fclose(dataFile);
	return 0;
#endif // ITERATIONS
#ifdef TEST_OTCHET
	long double a = -1;
	long double b = 1;
	long double y_a = exp(1);
	long double eps = 1e-10;
	long double h = 0.5;
	long double x_first = a;
	long double y_first = y_a;
	long double x_second = -0.5;
	long double* xMass = NULL;
	long double* yMass = NULL;
	int iterations = FindSolution(a, b, y_a, eps, f, &xMass, &yMass);
	printf("\n\n%d", iterations);
#endif

}
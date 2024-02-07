//
// Created by cest on 1/31/24.
//

#include "numpy.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void linspace(double *arr, double start, double end, int n) {
    double step = (end - start) / (n - 1);
    for(int i = 0; i < n; i++) {
        arr[i] = start + i * step;
    }
}

void meshgrid(double *x, double *y, int n, double **X, double **Y) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            X[i][j] = x[j];
            Y[i][j] = y[i];
        }
    }
}

double** zeros(int rows, int cols) {
    double** matrix = (double**) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*) calloc(cols, sizeof(double));  // Use calloc to initialize to 0
    }
    return matrix;
}

double** full(int rows, int cols, int value) {
    double** matrix = (double**) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = value;
        }
    }
    return matrix;
}

double** ones(int rows, int cols) {
    return full(rows, cols, 1);
}

void freeMatrix(double** matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

Array arange(int start, int stop, int step) {
    Array array;
    array.length = (int)((stop - start) / step) + 1;
    array.array = (double *)malloc(array.length * sizeof(double));
    for (int i = 0; i < array.length; i++) {
        array.array[i] = step * (i + 1);
    }
    return array;
}

typedef int (*MathFunc)(float, int);

int do_math(float arg1, int arg2) {
    return arg2;
}

int call_a_func(MathFunc call_this) {
    int output = call_this(5.5, 7);

    return output;
}

#define q	3		/* for 2^3 points */
#define N	(1<<q)		/* N-point FFT, iFFT */

typedef float real;
typedef struct{real Re; real Im;} complex;

#ifndef PI
# define PI	3.14159265358979323846264338327950288
#endif
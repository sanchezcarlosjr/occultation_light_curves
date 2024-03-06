//
// Created by cest on 1/31/24.
//

#include "numpy.h"



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

Matrix zeros(int rows, int cols) {
    double** matrix = (double**) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*) calloc(cols, sizeof(double));  // Use calloc to initialize to 0
    }
    return matrix;
}

Matrix full(int rows, int cols, int value) {
    double** matrix = (double**) malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = value;
        }
    }
    return matrix;
}

Matrix ones(int rows, int cols) {
    return full(rows, cols, 1);
}

void freeMatrix(Matrix matrix, int rows) {
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


// Function to initialize a 2D array of complex numbers with zeros
void zerosComplex(ComplexMatrix array, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            array[i][j] = 0.0 + 0.0 * I; // Initialize to complex zero
        }
    }
}

ComplexMatrix allocateComplex2DArray(int rows, int cols) {
    ComplexMatrix array = (ComplexMatrix) malloc(rows * sizeof(double complex *));
    if (!array) return NULL;

    for (int i = 0; i < rows; i++) {
        array[i] = (double complex *) malloc(cols * sizeof(double complex));
        if (!array[i]) {
            // Free previously allocated memory in case of failure
            for (int j = 0; j < i; j++) free(array[j]);
            free(array);
            return NULL;
        }
    }
    return array;
}

// Function to free a dynamically allocated 2D array of complex numbers
void freeComplex2DArray(ComplexMatrix array, int rows) {
    for (int i = 0; i < rows; i++) free(array[i]);
    free(array);
}

typedef int (*MathFunc)(float, int);

int do_math(float arg1, int arg2) {
    return arg2;
}

int call_a_func(MathFunc call_this) {
    int output = call_this(5.5, 7);

    return output;
}

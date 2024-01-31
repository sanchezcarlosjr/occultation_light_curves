//
// Created by cest on 1/31/24.
//

#include "numpy.h"
#include <stdlib.h>


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
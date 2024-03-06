//
// Created by cest on 1/31/24.
//

#ifndef OCCULTATION_LIGHT_CURVES_NUMPY_H
#define OCCULTATION_LIGHT_CURVES_NUMPY_H
typedef struct {
    int length;
    double* array;
} Array;

typedef double** Matrix;

void linspace(double *arr, double start, double end, int n);
void meshgrid(double *x, double *y, int n, double **X, double **Y);
Matrix zeros(int rows, int cols);
Matrix full(int rows, int cols, int value);
Matrix ones(int rows, int cols);
void freeMatrix(double** matrix, int rows);
Array arange(int start, int stop, int step);

#endif //OCCULTATION_LIGHT_CURVES_NUMPY_H

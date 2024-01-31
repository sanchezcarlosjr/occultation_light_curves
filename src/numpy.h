//
// Created by cest on 1/31/24.
//

#ifndef OCCULTATION_LIGHT_CURVES_NUMPY_H
#define OCCULTATION_LIGHT_CURVES_NUMPY_H

void linspace(double *arr, double start, double end, int n);
void meshgrid(double *x, double *y, int n, double **X, double **Y);
double** zeros(int rows, int cols);
double** full(int rows, int cols, int value);
double** ones(int rows, int cols);
void freeMatrix(double** matrix, int rows);


#endif //OCCULTATION_LIGHT_CURVES_NUMPY_H

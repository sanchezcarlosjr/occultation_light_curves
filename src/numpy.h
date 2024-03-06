//
// Created by cest on 1/31/24.
//
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>

#ifndef OCCULTATION_LIGHT_CURVES_NUMPY_H
#define OCCULTATION_LIGHT_CURVES_NUMPY_H
typedef struct {
    int length;
    double* array;
} Array;

typedef double** Matrix;
typedef double complex** ComplexMatrix;

void linspace(double *arr, double start, double end, int n);
void meshgrid(double *x, double *y, int n, double **X, double **Y);
Matrix zeros(int rows, int cols);
Matrix full(int rows, int cols, int value);
Matrix ones(int rows, int cols);
void freeMatrix(double** matrix, int rows);
Array arange(int start, int stop, int step);
void zerosComplex(ComplexMatrix array, int rows, int cols);
ComplexMatrix allocateComplex2DArray(int rows, int cols);
void freeComplex2DArray(ComplexMatrix array, int rows);

#endif //OCCULTATION_LIGHT_CURVES_NUMPY_H

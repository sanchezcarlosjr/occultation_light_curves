//
// Created by cest on 1/31/24.
//

#ifndef OCCULTATION_LIGHT_CURVES_DIFFRACTION_H
#define OCCULTATION_LIGHT_CURVES_DIFFRACTION_H
#include <gsl/gsl_matrix.h>

double SNR_TAOS2(double mV);
double calcPlano(double d, double lmda, double ua);
gsl_matrix* pupilCO(int M, double D, double d);

#endif //OCCULTATION_LIGHT_CURVES_DIFFRACTION_H

//
// Created by cest on 1/31/24.
//

#ifndef OCCULTATION_LIGHT_CURVES_SPECTRA_H
#define OCCULTATION_LIGHT_CURVES_SPECTRA_H
#include <stddef.h>
#include <gsl/gsl_matrix.h>
#include "hdf5.h"

gsl_matrix* readHDF5File(const char* filename, size_t rows, size_t cols);

#endif //OCCULTATION_LIGHT_CURVES_SPECTRA_H

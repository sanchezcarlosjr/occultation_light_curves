//
// Created by cest on 3/6/24.
//

#ifndef OCCULTATION_LIGHT_CURVES_DATA_WRITERS_H
#define OCCULTATION_LIGHT_CURVES_DATA_WRITERS_H
#include <gsl/gsl_matrix.h>
#include "hdf5.h"

void writeHDF5File(gsl_matrix *gslMatrix, const char *output, const char *dataset);

#endif //OCCULTATION_LIGHT_CURVES_DATA_WRITERS_H
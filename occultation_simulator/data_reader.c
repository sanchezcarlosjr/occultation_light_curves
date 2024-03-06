//
// Created by cest on 1/31/24.
//

#include "data_reader.h"
#include "stdio.h"
#include <gsl/gsl_matrix.h>


/*
 * Python:
import numpy as np

# Example 2D NumPy array
arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float64)

# Save to a binary file
arr.tofile("matrix.bin")

C:

gsl_matrix* matrix = readNumpyBinaryFile("matrix.bin", 3, 3);

 */
gsl_matrix* readNumpyBinaryFile(const char* filename, size_t rows, size_t cols) {
    FILE* file = fopen(filename, "rb");
    if (!file) {
        perror("Failed to open file");
        return NULL;
    }
    // Allocate a GSL matrix
    gsl_matrix *m = gsl_matrix_alloc(rows, cols);

    // Read the binary data directly into the GSL matrix's data block
    if (fread(m->data, sizeof(double), rows * cols, file) != rows * cols) {
        fprintf(stderr, "Error reading matrix data from the binary file.\n");
        gsl_matrix_free(m);
        fclose(file);
        return NULL;
    }

    return m;
}
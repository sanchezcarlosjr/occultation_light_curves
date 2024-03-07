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
gsl_matrix* readHDF5File(const char* filename, size_t rows, size_t cols) {
    hid_t file_id, dataset_id;
    gsl_matrix* matrix = gsl_matrix_alloc(rows, cols);
    unsigned mode = H5F_ACC_RDONLY;

    if ((file_id = H5Fopen(filename, mode, H5P_DEFAULT)) < 0) {
        fprintf(stderr, "unable to open the file: %s\n", filename);
        gsl_matrix_free(matrix);
        return NULL;
    }

    if ((dataset_id = H5Dopen2(file_id, "dataset", H5P_DEFAULT)) < 0) {
        fprintf(stderr, "unable to open the dataset in file: %s\n", filename);
        H5Fclose(file_id);
        gsl_matrix_free(matrix);
        return NULL;
    }

    if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix->data) < 0) {
        fprintf(stderr, "unable to read the dataset from file: %s\n", filename);
        H5Dclose(dataset_id);
        H5Fclose(file_id);
        gsl_matrix_free(matrix);
        return NULL;
    }

    // Successfully read the data, now close the dataset and file handles
    H5Dclose(dataset_id);
    H5Fclose(file_id);

    return matrix;
}
//
// Created by cest on 3/6/24.
//

#include "data_writers.h"

void writeHDF5File(gsl_matrix *gslMatrix, const char *output, const char *dataset) {
    // Identifiers for the file, dataset, and dataspace. Initialize to invalid values.
    hid_t file_id = -1, dataset_id = -1, dataspace_id = -1;

    // Create a new HDF5 file (H5F_ACC_TRUNC ensures the file is overwritten if it exists)
    file_id = H5Fcreate(output, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error creating file: %s\n", output);
        goto cleanup;
    }

    // Define the dimensions of the dataset
    hsize_t dims[2] = {gslMatrix->size1, gslMatrix->size2};

    // Create a simple dataspace with the defined dimensions
    dataspace_id = H5Screate_simple(2, dims, NULL);
    if (dataspace_id < 0) {
        fprintf(stderr, "Error creating dataspace\n");
        goto cleanup;
    }

    // Create the dataset to hold the matrix data
    dataset_id = H5Dcreate2(file_id, dataset, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error creating dataset\n");
        goto cleanup;
    }

    // Write the matrix data to the dataset
    if (H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, gslMatrix->data) < 0) {
        fprintf(stderr, "Error writing to dataset\n");
    }

    cleanup:
    // Attempt to close the dataspace if it was successfully created
    if (dataspace_id >= 0) {
        H5Sclose(dataspace_id);
    }
    // Attempt to close the dataset if it was successfully created
    if (dataset_id >= 0) {
        H5Dclose(dataset_id);
    }
    // Attempt to close the file if it was successfully created
    if (file_id >= 0) {
        H5Fclose(file_id);
    }
}

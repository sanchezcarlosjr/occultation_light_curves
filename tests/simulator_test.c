#include <unity.h>
#include "diffraction.h"
#include "data_reader.h"
#include "numpy.h"

double M=2048;
double lambda=600e-9;
double nLamb=10;
double d=3000;
double ua=45;

void setUp(void) {
    // Initialize your test setup here
}

void tearDown(void) {
    // Clean up your test setup here
}

void test_it_calculates_the_optimal_plane_size(void)
{
    double D = calcPlano(d, lambda, ua);
    TEST_ASSERT_FLOAT_WITHIN(0.000000000001, D, 142112.63138792414);
}

void test_it_generates_a_circular_obstruction(void) {
    gsl_matrix* expected_matrix = readHDF5File("circle_object_provided_by_pupil_co.h5", M, M);
    TEST_ASSERT_EQUAL(expected_matrix->size1, M);
    TEST_ASSERT_EQUAL(expected_matrix->size2, M);
    TEST_ASSERT_EQUAL(gsl_matrix_get(expected_matrix, 0, 0), 1);
    double D = calcPlano(d, lambda, ua);
    gsl_matrix* actual_matrix = pupilCO(M, D, d);
    int result = gsl_matrix_equal_with_tolerance(expected_matrix, actual_matrix, 0.000000000001);
    TEST_ASSERT_EQUAL(result, 1);
    gsl_matrix_free(expected_matrix);
    gsl_matrix_free(actual_matrix);
}

void test_it_generates_a_contact_binary(void) {
    gsl_matrix* expected_matrix = readHDF5File("contact_binary_provided_by_pupil_doble.h5", M, M);
    TEST_ASSERT_EQUAL(expected_matrix->size1, M);
    TEST_ASSERT_EQUAL(expected_matrix->size2, M);
    TEST_ASSERT_EQUAL(gsl_matrix_get(expected_matrix, 0, 0), 1);
    gsl_matrix* actual_matrix = pupilDoble(M, calcPlano(d, lambda, ua), d);
    TEST_ASSERT_EQUAL(gsl_matrix_get(actual_matrix, 0, 0), 1);
    TEST_ASSERT_EQUAL(gsl_matrix_get(actual_matrix, 1024, 1024), 0);
    int result = gsl_matrix_equal_with_tolerance(expected_matrix, actual_matrix, 0);
    TEST_ASSERT_EQUAL(result, 1);
    gsl_matrix_free(expected_matrix);
    gsl_matrix_free(actual_matrix);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_it_calculates_the_optimal_plane_size);
    RUN_TEST(test_it_generates_a_circular_obstruction);
    RUN_TEST(test_it_generates_a_contact_binary);

    return UNITY_END();
}
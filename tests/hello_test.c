#include <unity.h>
#include "diffraction.h"

void setUp(void) {
    // Initialize your test setup here
}

void tearDown(void) {
    // Clean up your test setup here
}

void test_it_calculates_the_optimal_plane_size(void)
{
    double d = 3000;
    double lambda = 600e-9;
    double ua = 45;
    double D = calcPlano(d, lambda, ua);
    TEST_ASSERT_FLOAT_WITHIN(0.000000000001, D, 142112.63138792414);
}

int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_it_calculates_the_optimal_plane_size);

    return UNITY_END();
}
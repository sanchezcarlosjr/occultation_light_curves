#include <gtest/gtest.h>
extern "C" {
   #include "numpy.h"
}

TEST(HelloTest, BasicAssertions) {
    EXPECT_STRNE("hello", "world");
    EXPECT_EQ(7 * 6, 42);
}
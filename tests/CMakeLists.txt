cmake_minimum_required(VERSION 3.6)

# Fetch Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/release-1.10.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Define your test executable
add_executable(your_test_executable tests.c)

# Link your test executable against Google Test and your library
target_link_libraries(your_test_executable gtest_main your_library_name)

# Add tests to CMake's test runner
include(GoogleTest)
gtest_discover_tests(your_test_executable)
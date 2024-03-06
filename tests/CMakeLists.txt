cmake_minimum_required(VERSION 3.27)
project(OccultationSimulatorTesting VERSION 1.0.0 LANGUAGES C CXX)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/31993dfa6b47e11c7a6ef67cfa8af90892b9bd1c.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(hello_test hello_test.cc)
target_link_libraries(hello_test PUBLIC GTest::gtest_main ${LIBRARY_NAME})

include(GoogleTest)
gtest_discover_tests(hello_test)
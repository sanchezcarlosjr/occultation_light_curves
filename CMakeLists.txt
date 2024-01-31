cmake_minimum_required(VERSION 3.27)
project(occultation_light_curves C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file(GLOB SRC_FILES "src/*.c")

add_executable(occultation_light_curves ${SRC_FILES})

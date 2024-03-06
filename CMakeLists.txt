cmake_minimum_required(VERSION 3.27)
project(occultation_light_curves C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(occultation_light_curves src/main.c src/numpy.c src/numpy.h src/diffraction.c src/diffraction.h
        src/spectra.c
        src/spectra.h)


find_package(FFTW3 REQUIRED COMPONENTS DOUBLE_LIB)

# Include FFTW headers
include_directories(${FFTW3_INCLUDE_DIRS})

target_link_libraries(occultation_light_curves m ${FFTW3_LIBRARIES})
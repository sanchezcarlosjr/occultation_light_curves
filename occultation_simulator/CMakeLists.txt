add_library(${LIBRARY_NAME} STATIC numpy.c diffraction.c data_reader.c)
find_package(GSL REQUIRED)
find_package(FFTW3 REQUIRED COMPONENTS DOUBLE_LIB)

target_include_directories(${LIBRARY_NAME} PUBLIC "./")
target_link_libraries(${LIBRARY_NAME} PUBLIC gsl gslcblas m fftw3)

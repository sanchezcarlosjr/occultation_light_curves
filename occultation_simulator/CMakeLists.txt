add_library(${LIBRARY_NAME} STATIC numpy.c diffraction.c data_readers.c data_writers.c)
find_package(GSL REQUIRED)
find_package(FFTW3 REQUIRED COMPONENTS DOUBLE_LIB)
find_package(HDF5 REQUIRED COMPONENTS C HL)

target_include_directories(${LIBRARY_NAME} PUBLIC "./")
target_link_libraries(${LIBRARY_NAME} PUBLIC gsl gslcblas m fftw3 ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

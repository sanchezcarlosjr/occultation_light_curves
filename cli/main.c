#include "numpy.h"
#include "cmdline.h"
#include "diffraction.h"
#include "data_writers.h"
#include <stdlib.h>

int main(int argc, char *argv[]) {
    struct gengetopt_args_info args_info;

    if (cmdline_parser(argc, argv, &args_info) != 0) {
        exit(1);
    }

    double D = calcPlano(args_info.d_arg, args_info.lamb_arg, args_info.ua_arg);
    gsl_matrix* O1 = pupilCO(args_info.M_arg, D, args_info.d_arg);
    writeHDF5File(O1, args_info.output_arg, "O1");

    gsl_matrix_free(O1);


    cmdline_parser_free(&args_info);
    return 0;
}
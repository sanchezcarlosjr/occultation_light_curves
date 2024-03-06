#include "cmdline.h"
#include <stdlib.h>

int main(int argc, char *argv[]) {
    struct gengetopt_args_info args_info;

    if (cmdline_parser(argc, argv, &args_info) != 0) {
        exit(1);
    }

    if (args_info.verbose_given) {
        // Verbose output
    }

    if (args_info.output_given) {
        // Handle output file
    }

    cmdline_parser_free(&args_info);
    return 0;
}
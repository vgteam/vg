/** \file version_main.cpp
 *
 * Defines the "vg version" subcommand, which evaluates graphs and alignments.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../version.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_version(char** argv) {
    cerr << "usage: " << argv[0] << " version" << endl
         << "options: " << endl
         << "  -s, --slug           print only the one-line, whitespace-free version string" << endl
         << "  -h, --help           print this help message to stderr and exit" << endl
         << endl;
}

int main_version(int argc, char** argv){

    bool slug_only = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"slug", no_argument, 0, 's'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}};
        int option_index = 0;
        c = getopt_long(argc, argv, "sh?",
                        long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 's':
            slug_only = true;
            break;
        case 'h':
        case '?':
        default:
            help_version(argv);
            exit(1);
        }
    }

    if (argc > optind) {
        help_version(argv);
        return 1;
    }

    cout << (slug_only ? Version::get_version() : Version::get_long()) << endl;
    return 0;
}
// Register subcommand
static Subcommand vg_version("version", "version information", DEVELOPMENT, main_version);


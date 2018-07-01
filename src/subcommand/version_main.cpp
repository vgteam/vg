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

void help_version(char** argv){
    cerr << "usage: " << argv[0] << " version" << endl
         << "options: " << endl
         << endl;
}

int main_version(int argc, char** argv){

    if (argc != 2) {
        help_version(argv);
        return 1;
    }

    cout << Version::get_long() << endl;
    return 0;
}
// Register subcommand
static Subcommand vg_version("version", "version information", DEVELOPMENT, main_version);


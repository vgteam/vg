#include <getopt.h>

#include <string>
#include <vector>

#include <subcommand.hpp>

#include "../primer_filter.hpp"
#include "../snarl_distance_index.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_primers(char** argv) {
    cerr << "usage: " << argv[0] << " primers [options] input.primer3 > filtered_primers.out" << endl
         << endl
         << "options:" << endl
         << "    -z, --zero-variance        Allow no variance in the product" << endl
         << "    -l, --tolerance INT        Allow this much difference between minimum and maximum sizes compared to the linear product size" << endl
         << "    -n, --minimum-size INT     Minimum product size allowed (has precedence over --tolerance)" << endl
         << "    -m, --maximum-size INT     Maximum product size allowed (has precedence over --tolerance)" << endl;
}

int main_primers(int argc, char** argv) {
    
    if (argc == 2) {
        help_priemrs(argv)
    }
}
/*
 * Define the "vg remove_duplicate" subcommand, which remove the duplicate PCRs
 *
 *
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"
#include "../vg.hpp"
#include "../algorithms/gfa_to_handle.hpp"
#include "../algorithms/id_sort.hpp"
#include <vg/io/vpkg.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_rmvdup(char **argv){
// TODO: add whatever option is needed to this list. Change long_option and getopt_long if want to add an option.
// TODO: see what input and output file formats is possible
    cerr << "usage: " << argv[0] << " rmvdup [options] inputfile.gam > output.gam " << endl
            << "Remove duplicate PCRs from the input file" << endl
            << "  -p, --progress               Show progress." << endl
            << "    -t, --threads N            number of threads to use" << endl;

}
int main_rmvdup(int argc, char *argv[]){
    string filename;
    bool show_progress = false;
    int threads = 1;

    int c;
    optind = 2;  // force optind past command positional argument
    while(true){
        static struct option long_options[] ={
                {"help", no_argument, 0, 'h'},
                {"progress", no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hpt",
                         long_options, &option_index);

        if (c == -1) break;

        switch (c) {

            case 'p':
                show_progress = true;
                break;

            case 't':
                threads = parse<int>(optarg);
                omp_set_num_threads(threads);
                break;

            case 'h':
            case '?':
                help_rmvdup(argv);
                exit(1);
                break;

            default:
                help_rmvdup(argv);
                abort();
        }
    }

    if (optind >= argc) {
        cerr << "[vg view] error: no filename given" << endl;
        exit(1);
    }
    string file_name = get_input_file_name(optind, argc, argv);
//    cout << file_name << endl;


    return 0;

}

// Register subcommand
static Subcommand vg_removeduplicate("rmvdup", "Remove duplicate PCRs from the input file", WIDGET, main_rmvdup);
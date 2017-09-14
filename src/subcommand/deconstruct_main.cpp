/** \file deconstruct_main.cpp
 *
 * Defines the "vg deconstruct" subcommand, which turns graphs back into VCFs.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../deconstructor.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_deconstruct(char** argv){
    cerr << "usage: " << argv[0] << " deconstruct [options] -p <PATH> <my_graph>.vg" << endl
         << "Outputs VCF records for Snarls present in a graph (relative to a chosen reference path)." << endl
         << "options: " << endl
         << "--path / -p     REQUIRED: A reference path to deconstruct against." << endl
         << endl;
}

int main_deconstruct(int argc, char** argv){
    //cerr << "WARNING: EXPERIMENTAL" << endl;
    if (argc <= 2) {
        help_deconstruct(argv);
        return 1;
    }

    vector<string> refpaths;
    string graphname;
    string outfile = "";
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"path", required_argument, 0, 'p'},
                {0, 0, 0, 0}

            };

            int option_index = 0;
            c = getopt_long (argc, argv, "hp:",
                    long_options, &option_index);

            // Detect the end of the options.
            if (c == -1)
                break;

            switch (c)
            {
                case 'p':
                    refpaths = split(optarg, ",");
                    break;
                case '?':
                case 'h':
                    help_deconstruct(argv);
                    return 1;
                default:
                    help_deconstruct(argv);
                    abort();
            }

        }
        graphname = argv[optind];
        vg::VG* graph;
        if (!graphname.empty()){
            ifstream gstream(graphname);
            graph = new vg::VG(gstream);
        }

        // load graph

        // Deconstruct
        Deconstructor dd;
        dd.deconstruct(refpaths, graph);
    return 0;
}

// Register subcommand
static Subcommand vg_deconstruct("deconstruct", "convert a graph into VCF relative to a reference", main_deconstruct);


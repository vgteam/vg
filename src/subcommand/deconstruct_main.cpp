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
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_deconstruct(char** argv){
    cerr << "usage: " << argv[0] << " deconstruct [options] -p <PATH> <my_graph>.vg" << endl
         << "Outputs VCF records for Snarls present in a graph (relative to a chosen reference path)." << endl
         << "options: " << endl
         << "    -p, --path NAME        REQUIRED: A reference path to deconstruct against (comma-separated list accepted)." << endl
         << "    -r, --snarls FILE      Snarls file (from vg snarls) to avoid recomputing." << endl
         << "    -e, --path-traversals  Only consider traversals that correspond to paths in the grpah." << endl
         << "    -t, --threads N        Use N threads" << endl
         << "    -v, --verbose          Print some status messages" << endl
         << endl;
}

int main_deconstruct(int argc, char** argv){
    if (argc <= 2) {
        help_deconstruct(argv);
        return 1;
    }

    vector<string> refpaths;
    string graphname;
    string snarl_file_name;
    bool path_restricted_traversals = false;
    bool show_progress = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"path", required_argument, 0, 'p'},
                {"snarls", required_argument, 0, 'r'},
                {"path-traversals", no_argument, 0, 'e'},
                {"threads", required_argument, 0, 't'},
                {"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}

            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hp:r:et:v",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            refpaths = split(optarg, ",");
            break;
        case 'r':
            snarl_file_name = optarg;
            break;
        case 'e':
            path_restricted_traversals = true;
            break;
        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;
        case 'v':
            show_progress = true;
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
    
    string graph_file_name = get_input_file_name(optind, argc, argv);

    vg::VG* graph;
    get_input_file(graph_file_name, [&](istream& in) {
            if (show_progress) {
                cerr << "Loading graph" << endl;
            }
            graph = new VG(in);
    });

    // Load or compute the snarls
    unique_ptr<SnarlManager> snarl_manager;    
    if (!snarl_file_name.empty()) {
        ifstream snarl_file(snarl_file_name.c_str());
        if (!snarl_file) {
            cerr << "[vg deconstruct]: Unable to load snarls file: " << snarl_file_name << endl;
            return 1;
        }
        if (show_progress) {
            cerr << "Loading snarls" << endl;
        }
        snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
    } else {
        CactusSnarlFinder finder(*graph);
        if (show_progress) {
            cerr << "Finding snarls" << endl;
        }
        snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls())));
    }

    // Deconstruct
    Deconstructor dd;
    if (show_progress) {
        cerr << "Decsontructing top-level snarls" << endl;
    }
    dd.deconstruct(refpaths, graph, snarl_manager.get(), path_restricted_traversals);
    return 0;
}

// Register subcommand
static Subcommand vg_deconstruct("deconstruct", "convert a graph into VCF relative to a reference", main_deconstruct);


/** \file sort_main.cpp
 *
 * Defines the "vg sort" subcommand, which sorts graph nodes.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../gfa.hpp"
#include "../flow_sort.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_sort(char** argv){
    cerr << "usage: " << argv[0] << " sort [options] -i <input_file> -r <reference_name> > sorted.vg " << endl
         << "options: " << endl
         << "           -g, --gfa              input in GFA format" << endl
         << "           -i, --in               input file" << endl
         << "           -r, --ref              reference name" << endl
         << "           -w, --without-grooming no grooming mode" << endl
         << "           -f, --fast             sort using Eades algorithm, otherwise max-flow sorting is used" << endl   
         << endl;
}

int main_sort(int argc, char *argv[]) {

    //default input format is vg
    bool gfa_input = false;
    string file_name = "";
    string reference_name = "";
    bool without_grooming = false;
    bool use_fast_algorithm = false;
    int c;
    while (true) {
        static struct option long_options[] =
            {
                {"gfa", no_argument, 0, 'g'},
                {"in", required_argument, 0, 'i'},
                {"ref", required_argument, 0, 'r'},
                {"without-grooming", no_argument, 0, 'w'},
                {"fast", no_argument, 0, 'f'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "i:r:gwf",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'g':
            gfa_input = true;
            break;
        case 'r':
            reference_name = optarg;
            break;
        case 'i':
            file_name = optarg;
            break;
        case 'w':
            without_grooming = true;
            break;
        case 'f':
            use_fast_algorithm = true;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_sort(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }
  
    if (reference_name.empty() || file_name.empty()) {
        help_sort(argv);
        exit(1);
    }
    
    ifstream in;
    std::unique_ptr<VG> graph;
    {
        in.open(file_name.c_str());        
        if (gfa_input) {
            graph.reset(new VG());
            if (!gfa_to_graph(in, graph.get())) {
                // GFA loading has failed because the file is invalid
                exit(1);
            }
        } else {
            graph.reset(new VG(in));
        }
    }
    FlowSort flow_sort(*graph.get());
    if (use_fast_algorithm) {
        flow_sort.fast_linear_sort(reference_name, !without_grooming);
    } else {
        flow_sort.max_flow_sort(reference_name);
    }
    
    graph->serialize_to_ostream(std::cout);
    in.close();
    return 0;
}

// Register subcommand
static Subcommand vg_sort("sort", "sort variant graph using max flow algorithm or Eades fast heuristic algorithm", main_sort);


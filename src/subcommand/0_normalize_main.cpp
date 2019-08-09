// mod.cpp: define the "normalize" subcommand, which realigns snarls to produce more
// efficient representations of snarls.

#include <getopt.h>
#include <omp.h>
#include <unistd.h>

#include "subcommand.hpp"

#include "../../include/sglib/hash_graph.hpp"
#include "../../include/vg/io/vpkg.hpp"
#include "../algorithms/0_draft_haplotype_realignment.hpp"
#include "../algorithms/0_draft_snarl_normalization_evaluation.cpp"
#include "../gbwt_helper.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_normalize(char **argv) {
    cerr << "usage: " << argv[0] << " normalize [options] <graph.vg> >[mod.vg]" << endl
         << "Modifies snarls, outputs modified on stdout." << endl
         << endl
         << "options:" << endl
         << "    -g, --gbwt       gbwt corresponding to hashgraph." << endl
         << "    -s, --snarls       snarls file corresponding to hashgraph." << endl;
}

int main_normalize(int argc, char **argv) {

    if (argc == 2) {
        help_normalize(argv);
        return 1;
    }

    bool evaluate = false;
    bool normalize = false;
    string gbwt;
    string snarls;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {{"help", no_argument, 0, 'h'},
             {"gbwt", required_argument, 0, 'g'},
             {"snarls", required_argument, 0, 's'},
             {"evaluate", no_argument, 0, 'e'},
             {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "g:s:eh", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {

        case 'g':
            gbwt = optarg;
            normalize = true;
            break;

        case 'e':
            evaluate = true;
            break;

        case 's':
            snarls = optarg;
            break;

        default:
            abort();
        }
    }

    sglib::HashGraph *graph;
    get_input_file(optind, argc, argv,
                   [&](istream &in) { graph = new sglib::HashGraph(in); });

    if (normalize) {
        /// Build the gbwt:
        ifstream gbwt_stream;
        gbwt_stream.open(gbwt);

        // Load the GBWT from its container
        unique_ptr<gbwt::GBWT> gbwt;
        gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
        GBWTGraph haploGraph = vg::GBWTGraph(*gbwt, *graph);

        std::ifstream snarl_stream;
        string snarl_file = snarls;
        snarl_stream.open(snarl_file);

        if (!snarl_stream) {
            cerr << "error:[vg mod] Cannot open Snarls file " << snarl_file << endl;
            exit(1);
        }

        // run test code on all snarls in graph.
        disambiguate_top_level_snarls(*graph, haploGraph, snarl_stream);
    }

    if (evaluate) {
        std::ifstream snarl_stream;
        string snarl_file = snarls;
        snarl_stream.open(snarl_file);
        cerr << "about to evaluate normalized snarls" << endl;
        vg::evaluate_normalized_snarls(snarl_stream);
    }

    // TODO: NOTE: this may be cumbersome code if we decide to add more argument types.
    // Consider changing.
    if (normalize) {
        graph->serialize(std::cout);
    }
    delete graph;

    return 0;
}

// Register subcommand
static Subcommand vg_normalize("normalize",
                               "edit snarls to reduce information duplication", TOOLKIT,
                               main_normalize);
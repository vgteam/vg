// mod.cpp: define the "normalize" subcommand, which realigns snarls to produce more
// efficient representations of snarls.

#include <getopt.h>
#include <omp.h>
#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>

#include "subcommand.hpp"

#include "../../include/bdsg/hash_graph.hpp"
#include "../../include/vg/io/vpkg.hpp"
// #include "../algorithms/0_draft_haplotype_realignment.hpp"
#include "../algorithms/0_draft_snarl_normalization_evaluation.cpp"
#include "../algorithms/0_oo_normalize_snarls.hpp"
#include "../gbwt_helper.hpp"

#include <chrono> // for high_resolution_clock

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_normalize(char **argv) {
    cerr
        << "usage: " << argv[0] << " normalize [options] <graph.vg> >[mod.vg]" << endl
        << "Modifies snarls, outputs modified on stdout." << endl
        << endl
        << "options:" << endl
        << "    -g, --gbwt       gbwt corresponding to hashgraph." << endl
        << "    -s, --snarls       snarls file corresponding to hashgraph." << endl
        << "    -m, --max_alignment_size       limits the number of threads that will "
           "be aligned in any snarl. If exceeded, program skips snarl. Default is 200 "
           "threads. If you don't want to skip any snarls based on thread count, enter 0."
        << endl
        << "    -s, --snarls       snarls file corresponding to hashgraph." << endl;
}

int main_normalize(int argc, char **argv) {

    if (argc == 2) {
        help_normalize(argv);
        return 1;
    }

    bool evaluate = false;
    bool normalize = false;
    int max_alignment_size = 200; // default cutoff is 200 threads in a snarl.
    string gbwt;
    string snarls;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {{"help", no_argument, 0, 'h'},
             {"gbwt", required_argument, 0, 'g'},
             {"snarls", required_argument, 0, 's'},
             {"max_alignment_size", optional_argument, 0, 'm'},
             {"evaluate", no_argument, 0, 'e'},
             {0, 0, 0, 0}};

        int option_index = 0;
        c = getopt_long(argc, argv, "hg:s:m:e", long_options, &option_index);

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

        case 'm':
            max_alignment_size = parse<int>(optarg);
            // if max_alignment_size is 0, then that signifies that it should actually be
            // infinite, i.e. that we should not exclude any snarls.
            if (max_alignment_size == 0) {
                max_alignment_size = INT_MAX;
            }

        case 's':
            snarls = optarg;
            break;

        default:
            abort();
        }
    }

    // //get a VG graph:
    // cerr << "getting a vg graph" << endl;
    // unique_ptr<VG> graph;
    // get_input_file(optind, argc, argv, [&](istream &in) {
    //     graph = vg::io::VPKG::load_one<VG>(in);
    // });
    // cerr << "got a vg graph" << endl;

    // //using vpkg for just HashGraph:
    // unique_ptr<MutablePathDeletableHandleGraph> graph;
    // get_input_file(optind, argc, argv, [&](istream &in) {
    //     graph = vg::io::VPKG::load_one<bdsg::HashGraph>(in);
    // });

    //getting graph of any type, requires vpkg wrapping:
    unique_ptr<MutablePathDeletableHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream &in) {
        graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
    });

    // cerr << "getting hashgraph" << endl;
    // bdsg::HashGraph *graph;
    // get_input_file(optind, argc, argv,
    //                [&](istream &in) { graph = new bdsg::HashGraph(in); });
    // cerr << "got hashgraph" << endl;

    if (normalize) {
        cerr << "running normalize!" << endl;

        /// Build the gbwt:
        ifstream gbwt_stream;
        gbwt_stream.open(gbwt);

        // Load the GBWT from its container
        unique_ptr<gbwt::GBWT> gbwt;
        gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
        gbwtgraph::GBWTGraph haploGraph = gbwtgraph::GBWTGraph(*gbwt, *graph);

        std::ifstream snarl_stream;
        string snarl_file = snarls;
        snarl_stream.open(snarl_file);

        if (!snarl_stream) {
            cerr << "error:[vg mod] Cannot open Snarls file " << snarl_file << endl;
            exit(1);
        }
        // Record start time
        auto start = chrono::high_resolution_clock::now();

        SnarlNormalizer normalizer =
            SnarlNormalizer(*graph, haploGraph, max_alignment_size);

        // run test code on all snarls in graph.
        normalizer.normalize_top_level_snarls(snarl_stream);

        // // run test code on all snarls in graph. (non obj-oriented code)
        // disambiguate_top_level_snarls(*graph, haploGraph, snarl_stream,
        // max_alignment_size);

        // Record end time
        auto finish = std::chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = finish - start;
        cerr << "Elapsed time: " << elapsed.count() << " s\n";
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
        vg::io::VPKG::save(*dynamic_cast<bdsg::HashGraph *>(graph.get()), cout);

        // graph->serialize(std::cout);
    }
    // delete graph;

    return 0;
}

// Register subcommand
static Subcommand vg_normalize("normalize",
                               "edit snarls to reduce information duplication", TOOLKIT,
                               main_normalize);
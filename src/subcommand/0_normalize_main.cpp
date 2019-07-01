// mod.cpp: define the "normalize" subcommand, which realigns snarls to produce more efficient representations of snarls.

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include "../hash_graph.hpp"
#include "../algorithms/0_draft_haplotype_realignment.hpp"
#include "../gbwt_helper.hpp"
#include "../../include/vg/io/vpkg.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_normalize(char** argv) {
    cerr << "usage: " << argv[0] << " normalize [options] <graph.vg> >[mod.vg]" << endl
         << "Modifies snarls, outputs modified on stdout." << endl
         << endl
         << "options:" << endl
         << "    -n, --normalize       normalizes a currently-hardcoded snarl from a graph." << endl;
}

int main_normalize(int argc, char** argv) {

    if (argc == 2) {
        help_normalize(argv);
        return 1;
    }

    bool normalize = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

        {
            {"help", no_argument, 0, 'h'},
            {"normalize", no_argument, 0, 'n'}, //TODO: change no_argument to required_argument, assuming we want one.
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "n", //TODO: change to "n:" later, when we have something to specify.
                long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'n':
            normalize = true;
        }
    }

    HashGraph* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new HashGraph(in);
    });

    if ( normalize ) 
    {
        /// Build the gbwt:
        ifstream gbwt_stream;
        string gbwt_name = "test/robin_haplotypes/threads_in_middle_example/chr10_subgraph_0_new.gbwt"; //Nodes 23493 to 23505
        gbwt_stream.open(gbwt_name);

        // Load the GBWT from its container
        unique_ptr<gbwt::GBWT> gbwt;
        gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
        GBWTGraph haploGraph = vg::GBWTGraph(*gbwt, *graph);

        std::ifstream snarl_stream;
        string snarl_file = "test/robin_haplotypes/threads_in_middle_example/chr10_subgraph_0_new.snarls";
        snarl_stream.open(snarl_file);
        
        if (!snarl_stream) {
            cerr << "error:[vg mod] Cannot open Snarls file " << snarl_file << endl;
            exit(1);
        }

        // run test code on all snarls in graph.
        disambiguate_top_level_snarls(*graph, haploGraph, snarl_stream);


        // /// Run test code on a single snarl:
        // vg::id_t source = 23493; vg::id_t sink = 23505;
        // disambiguate_snarl(*graph, haploGraph, source, sink);

    }

    graph->serialize(std::cout);
    delete graph;

    return 0;
}

// Register subcommand
static Subcommand vg_normalize("normalize", "edit snarls to reduce information duplication", TOOLKIT, main_normalize);

















//TODO: Remove JUNK:

        // vg::id_t source = 23251;//for robin_haplotypes/simple
        // vg::id_t sink = 23257;//for robin_haplotypes/simple
        // /Testing gbwt_helper.hpp's for_each_kmer function. This issue is that I don't know how to construct a gbwt::GBWT haplotypes object. Nor do I know how to determine what size k I should use.
        // vg::id_t source = 23251;//for robin_haplotypes/simple
        // vg::id_t sink = 23257;//for robin_haplotypes/simple
        // clean_snarl_from_haplotypes(*graph, source, sink);
        // cerr << "done!" << endl;
        // vg::handle_t source_handle = graph->get_handle(source);
        // vg::handle_t sink_handle = graph->get_handle(sink);

        // vector<string> haplotypes = depth_first_haplotypes_to_strings(*graph, source, sink);
        // cerr << "finished depth_first, now on to reference." << endl;
        // vector<string> reference = get_paths(*graph, source_handle, sink_handle);

        // haplotypes.insert(end(haplotypes), begin(reference), end(reference));

        // cerr << "here goes!" << endl;
        // for(string haplotype : haplotypes) {
            
        //     cerr << haplotype << endl;
        // }
        // cerr << "done" << endl;
















        //     std::ifstream snarl_stream;
        //     snarl_stream.open(demo_0);
            
        //     if (!snarl_stream) {
        //         cerr << "error:[vg mod] Cannot open Snarls file " << demo_0 << endl;
        //         exit(1);
        //     }

        //     clean_all_snarls(*graph, snarl_stream);

        // string gbwt_name = "test/robin_haplotypes/simple/chr10_subgraph_2dels-shift-729006.gbwt";


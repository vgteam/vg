/** \file prune_main.cpp
 *
 * Defines the "vg prune" subcommand, which prunes the complex regions of the
 * graph for GCSA2 indexing.
 *
 * By default, pruning removes the nodes touched by paths of length
 * --kmer-length crossing more than --edge-max non-trivial edges. Graph
 * regions shorter than --kmer-length are also removed.
 * TODO: Optionally, nodes on the paths stored in the XG index are not removed.
 * TODO: Optionally, replace each pruned region with a set of disjoint paths
 * corresponding to all distinct XG paths and GBWT threads in that region.
 * TODO: We also need to output a mapping from duplicated node ids to original
 * node ids for GCSA2 construction.
 */

#include "subcommand.hpp"
#include "../vg.hpp"
#include "../xg.hpp"

#include <gbwt/gbwt.h>

#include <cstdlib>
#include <iostream>
#include <string>

#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;


struct PruningParameters
{
    const static int LENGTH = 16;
    const static int EDGE_MAX = 4;
};

void help_prune(char** argv) {
    std::cerr << "usage: " << argv[0] << " prune [options] <graph.vg> >[output.vg]" << std::endl;
    std::cerr << "Prunes the complex regions of the graph for GCSA2 indexing. Note that pruning" << std::endl;
    std::cerr << "removes paths from the graph." << std::endl;
    std::cerr << "pruning options:" << std::endl;
    std::cerr << "    -k, --kmer-length N    kmer length used for pruning (default: " << PruningParameters::LENGTH << ")" << std::endl;
    std::cerr << "    -e, --edge-max N       only consider paths making > N edge choices (default: " << PruningParameters::EDGE_MAX << ")" << std::endl;
    std::cerr << "path options:" << std::endl;
    std::cerr << "    -x, --xg-name FILE     preserve the paths in this XG index" << std::endl;
    std::cerr << "    -g, --gbwt-name FILE   unfold the threads in this GBWT index in the complex" << std::endl;
    std::cerr << "                           regions instead (requires -x)" << std::endl;
    std::cerr << "    -m, --mapping FILE     store the node mapping from -g in this file" << std::endl;
}

int main_prune(int argc, char** argv) {

    if (argc == 2) {
        help_prune(argv);
        return 1;
    }

    int kmer_length = PruningParameters::LENGTH;
    int edge_max = PruningParameters::EDGE_MAX;
    std::string xg_name, gbwt_name, mapping_name;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "kmer-length", required_argument, 0, 'k' },
            { "edge-max", required_argument, 0, 'e' },
            { "xg-name", required_argument, 0, 'x' },
            { "gbwt-name", required_argument, 0, 'g' },
            { "mapping", required_argument, 0, 'm' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "k:e:x:g:m:", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'k':
            kmer_length = stoi(optarg);
            break;
        case 'e':
            edge_max = stoi(optarg);
            break;
        case 'x':
            xg_name = optarg;
            break;
        case 'g':
            gbwt_name = optarg;
            break;
        case 'm':
            mapping_name = optarg;
            break;

        case 'h':
        case '?':
            help_prune(argv);
            return 1;
        default:
            std::abort();
        }
    }
    if (!(kmer_length > 0 && edge_max > 0)) {
        cerr << "[vg prune]: --kmer-length and --edge-max must be positive" << endl;
        return 1;
    }

    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });
    graph->paths.clear(); // Updating the paths can be expensive, so we just get rid of them.

    // Prune complex regions.
    if (xg_name.empty()) {
        graph->prune_complex_with_head_tail(kmer_length, edge_max);
    } else {
        std::cerr << "Not implemented yet!" << std::endl;
        // find complex regions
        // cluster complex regions into connected components
        // remove complex regions, except if they are on XG paths
        // alternatively unfold the complex regions
    }
    graph->prune_short_subgraphs(kmer_length);

    graph->serialize_to_ostream(std::cout);
    if (!mapping_name.empty()) {
        std::cerr << "Not implemented yet!" << std::endl;
        // (write mapping)
    }

    delete graph; graph = 0;

    return 0;
}

// Register subcommand
static Subcommand vg_prune("prune", "prune the graph for GCSA2 indexing", TOOLKIT, main_prune);


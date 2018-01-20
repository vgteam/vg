/** \file prune_main.cpp
 *
 * Defines the "vg prune" subcommand, which prunes the complex regions of the
 * graph for GCSA2 indexing.
 *
 * By default, pruning removes the nodes touched by paths of length
 * --kmer-length crossing more than --edge-max non-trivial edges. Graph
 * regions shorter than --subgraph_min are also removed. Pruning either removes
 * all embedded paths or preserves both the paths and the edges on non-alt
 * paths.
 * TODO: Optionally, replace each pruned region with a set of disjoint paths
 * corresponding to all distinct XG paths and GBWT threads in that region.
 * TODO: We also need to output a mapping from duplicated node ids to original
 * node ids for GCSA2 construction.
 */

#include "subcommand.hpp"
#include "../vg.hpp"
#include "../xg.hpp"

#include <gbwt/gbwt.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

#include <getopt.h>
#include <omp.h>

using namespace vg;
using namespace vg::subcommand;


struct PruningParameters
{
    const static int    KMER_LENGTH = 16;
    const static int    EDGE_MAX = 4;
    const static size_t SUBGRAPH_MIN = 33;
};

void help_prune(char** argv) {
    std::cerr << "usage: " << argv[0] << " prune [options] <graph.vg> >[output.vg]" << std::endl;
    std::cerr << "Prunes the complex regions of the graph for GCSA2 indexing. By default," << std::endl;
    std::cerr << "pruning removes the embedded paths." << std::endl;
    std::cerr << "pruning options:" << std::endl;
    std::cerr << "    -k, --kmer-length N    kmer length used for pruning (default: " << PruningParameters::KMER_LENGTH << ")" << std::endl;
    std::cerr << "    -e, --edge-max N       prune paths making > N edge choices (default: " << PruningParameters::EDGE_MAX << ")" << std::endl;
    std::cerr << "    -s, --subgraph-min N   prune subgraphs of < N bases (default: " << PruningParameters::SUBGRAPH_MIN << ")" << std::endl;
    std::cerr << "path options:" << std::endl;
    std::cerr << "    -p, --preserve-paths   preserve the embedded non-alt paths in the graph" << std::endl;
    std::cerr << "    -x, --xg-name FILE     use this XG index" << std::endl;
    std::cerr << "    -g, --gbwt-name FILE   unfold the complex regions by using the threads in" << std::endl;
    std::cerr << "                           this GBWT index (requires -x)" << std::endl;
    std::cerr << "    -m, --mapping FILE     store the node mapping from -g in this file" << std::endl;
    std::cerr << "other options:" << std::endl;
    std::cerr << "    -t, --threads N        use N threads (default: " << omp_get_max_threads() << ")" << std::endl;
}

int main_prune(int argc, char** argv) {

    if (argc == 2) {
        help_prune(argv);
        return 1;
    }

    int kmer_length = PruningParameters::KMER_LENGTH;
    int edge_max = PruningParameters::EDGE_MAX;
    size_t subgraph_min = PruningParameters::SUBGRAPH_MIN;
    int threads = omp_get_max_threads();
    bool preserve_paths = false;
    std::string xg_name, gbwt_name, mapping_name;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "kmer-length", required_argument, 0, 'k' },
            { "edge-max", required_argument, 0, 'e' },
            { "subgraph-min", required_argument, 0, 's' },
            { "preserve-paths", no_argument, 0, 'p' },
            { "xg-name", required_argument, 0, 'x' },
            { "gbwt-name", required_argument, 0, 'g' },
            { "mapping", required_argument, 0, 'm' },
            { "threads", required_argument, 0, 't' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "k:e:s:px:g:m:t:", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'k':
            kmer_length = stoi(optarg);
            break;
        case 'e':
            edge_max = stoi(optarg);
            break;
        case 's':
            subgraph_min = stoul(optarg);
            break;
        case 'p':
            preserve_paths = true;
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
        case 't':
            threads = stoi(optarg);
            threads = std::min(threads, omp_get_max_threads());
            threads = std::max(threads, 1);
            omp_set_num_threads(threads);
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

    // Prune complex regions.
    if (!preserve_paths) {
        graph->paths.clear();
    }
    graph->prune_complex_with_head_tail(kmer_length, edge_max, preserve_paths);
    // TODO: Another approach
    // load XG, GBWT
    // find complex regions
    // cluster complex regions into connected components
    // remove complex regions, except if they are on VG paths
    // unfold the complex regions
    graph->prune_short_subgraphs(subgraph_min);

    graph->serialize_to_ostream(std::cout);
    if (!mapping_name.empty()) {
        std::cerr << "Not implemented yet!" << std::endl;
        // (write mapping)
    }

    delete graph; graph = nullptr;

    return 0;
}

// Register subcommand
static Subcommand vg_prune("prune", "prune the graph for GCSA2 indexing", TOOLKIT, main_prune);


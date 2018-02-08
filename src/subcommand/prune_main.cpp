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
 *
 * With --gbwt-name and --xg-name, pruning unfolds the paths and threads in
 * the complex regions instead of pruning them. If --mapping is specified,
 * the mapping from duplicate node identifiers to the original identifiers is
 * stored in a file. This file can be used for building a GCSA2 index that
 * maps to the original graph.
 */

#include "../phase_unfolder.hpp"
#include "subcommand.hpp"

#include <gbwt/gbwt.h>

#include <cstdlib>
#include <iostream>
#include <regex>
#include <set>
#include <string>

#include <getopt.h>
#include <omp.h>

using namespace vg;
using namespace vg::subcommand;


struct PruningParameters
{
    const static int    KMER_LENGTH = 16;
    const static int    KMER_LENGTH_UNFOLDING = 16;
    const static int    EDGE_MAX = 4;
    const static int    EDGE_MAX_UNFOLDING = 4;
    const static size_t SUBGRAPH_MIN = 33;
    const static size_t SUBGRAPH_MIN_UNFOLDING = 33;
};

void help_prune(char** argv) {
    std::cerr << "usage: " << argv[0] << " prune [options] <graph.vg> >[output.vg]" << std::endl;
    std::cerr << "Prunes the complex regions of the graph for GCSA2 indexing. By default," << std::endl;
    std::cerr << "pruning removes the embedded paths." << std::endl;
    std::cerr << "pruning options:" << std::endl;
    std::cerr << "    -k, --kmer-length N    kmer length used for pruning" << std::endl;
    std::cerr << "                           (default: " << PruningParameters::KMER_LENGTH << "; " << PruningParameters::KMER_LENGTH_UNFOLDING << " with unfolding)" << std::endl;
    std::cerr << "    -e, --edge-max N       remove the edges on kmers making > N edge choices" << std::endl;
    std::cerr << "                           (default: " << PruningParameters::EDGE_MAX << "; " << PruningParameters::EDGE_MAX_UNFOLDING << " with unfolding)" << std::endl;
    std::cerr << "    -s, --subgraph-min N   remove subgraphs of < N bases" << std::endl;
    std::cerr << "                           (default: " << PruningParameters::SUBGRAPH_MIN << "; " << PruningParameters::SUBGRAPH_MIN_UNFOLDING << " with unfolding)" << std::endl;
    std::cerr << "    -P, --preserve-paths   preserve the non-alt paths and the edges on them" << std::endl;
    std::cerr << "unfolding options:" << std::endl;
    std::cerr << "    -x, --xg-name FILE     unfold the paths in this XG index (ignores -P)" << std::endl;
    std::cerr << "    -g, --gbwt-name FILE   unfold the threads in this GBWT index (requires -x)" << std::endl;
    std::cerr << "    -m, --mapping FILE     store the node mapping for duplicates in this file" << std::endl;
    std::cerr << "    -a, --append-mapping   append to the existing node mapping (requires -m)" << std::endl;
    std::cerr << "other options:" << std::endl;
    std::cerr << "    -p, --progress         show progress" << std::endl;
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
    bool preserve_paths = false, append_mapping = false, show_progress = false;
    std::string gbwt_name, xg_name, mapping_name;

    // Derived variables.
    bool kmer_length_set = false, edge_max_set = false, subgraph_min_set = false, unfold_paths = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "kmer-length", required_argument, 0, 'k' },
            { "edge-max", required_argument, 0, 'e' },
            { "subgraph-min", required_argument, 0, 's' },
            { "preserve-paths", no_argument, 0, 'P' },
            { "xg-name", required_argument, 0, 'x' },
            { "gbwt-name", required_argument, 0, 'g' },
            { "mapping", required_argument, 0, 'm' },
            { "append-mapping", no_argument, 0, 'a' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "k:e:s:Px:g:m:apt:", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'k':
            kmer_length = stoi(optarg);
            kmer_length_set = true;
            break;
        case 'e':
            edge_max = stoi(optarg);
            edge_max_set = true;
            break;
        case 's':
            subgraph_min = stoul(optarg);
            subgraph_min_set = true;
            break;
        case 'P':
            preserve_paths = true;
            break;
        case 'x':
            xg_name = optarg;
            unfold_paths = true;
            break;
        case 'g':
            gbwt_name = optarg;
            break;
        case 'm':
            mapping_name = optarg;
            break;
        case 'a':
            append_mapping = true;
            break;
        case 'p':
            show_progress = true;
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
        std::cerr << "[vg prune]: --kmer-length and --edge-max must be positive" << std::endl;
        return 1;
    }
    if (append_mapping && mapping_name.empty()) {
        std::cerr << "[vg prune]: parameter --append-mapping requires --mapping" << std::endl;
        return 1;
    }
    if (unfold_paths) {
        preserve_paths = false;
        if (!kmer_length_set) {
            kmer_length = PruningParameters::KMER_LENGTH_UNFOLDING;
        }
        if (!edge_max_set) {
            edge_max = PruningParameters::EDGE_MAX_UNFOLDING;
        }
        if (!subgraph_min_set) {
            subgraph_min = PruningParameters::SUBGRAPH_MIN_UNFOLDING;
        }
    }

    // Handle input.
    VG* graph;
    get_input_file(optind, argc, argv, [&](std::istream& in) {
        graph = new VG(in);
    });
    vg::id_t max_node_id = graph->max_node_id();
    if (show_progress) {
        std::cerr << "Loaded VG graph " << argv[optind - 1] << ": "
                  << graph->node_count() << " nodes, " << graph->edge_count() << " edges" << std::endl;
    }

    // Remove unnecessary paths.
    if (preserve_paths) {
        regex is_alt("_alt_.+_[0-9]+");
        std::set<std::string> to_remove;
        graph->paths.for_each([&](const Path& path) {
            if (regex_match(path.name(), is_alt)) {
                to_remove.insert(path.name());
            }
        });
        graph->paths.remove_paths(to_remove);
        if (show_progress) {
            std::cerr << "Removed non-alt paths" << std::endl;
        }
    } else {
        graph->paths.clear();
        if (show_progress) {
            std::cerr << "Removed all paths" << std::endl;
        }
    }

    // Prune the graph.
    graph->prune_complex_with_head_tail(kmer_length, edge_max, preserve_paths);
    if (show_progress) {
        std::cerr << "Pruned complex regions: "
                  << graph->node_count() << " nodes, " << graph->edge_count() << " edges" << std::endl;
    }
    graph->prune_short_subgraphs(subgraph_min);
    if (show_progress) {
        std::cerr << "Removed small subgraphs: "
                  << graph->node_count() << " nodes, " << graph->edge_count() << " edges" << std::endl;
    }

    // Unfold phase threads.
    if (unfold_paths) {
        xg::XG xg_index;
        get_input_file(xg_name, [&](std::istream& in) {
           xg_index.load(in); 
        });
        gbwt::GBWT gbwt_index;
        if (!gbwt_name.empty()) {
            get_input_file(gbwt_name, [&](std::istream& in) {
               gbwt_index.load(in);
            });
        }
        PhaseUnfolder unfolder(xg_index, gbwt_index, max_node_id + 1);
        if (append_mapping) {
            unfolder.read_mapping(mapping_name);
        }
        unfolder.unfold(*graph, show_progress);
        if (!mapping_name.empty()) {
            unfolder.write_mapping(mapping_name);
        }
    }

    // Serialize.
    graph->serialize_to_ostream(std::cout);
    if (show_progress) {
        std::cerr << "Serialized the graph: "
                  << graph->node_count() << " nodes, " << graph->edge_count() << " edges" << std::endl;
    }

    delete graph; graph = nullptr;

    return 0;
}

// Register subcommand
static Subcommand vg_prune("prune", "prune the graph for GCSA2 indexing", TOOLKIT, main_prune);


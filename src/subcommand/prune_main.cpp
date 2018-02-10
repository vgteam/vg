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

enum PruningMode { mode_prune, mode_preserve, mode_restore, mode_unfold };

std::string mode_name(PruningMode mode) {
    std::string result = "unknown";
    switch (mode) {
    case mode_prune:
        result = "(plain)";
        break;
    case mode_preserve:
        result = "--preserve-paths";
        break;
    case mode_restore:
        result = "--restore-paths";
        break;
    case mode_unfold:
        result = "--unfold-paths";
        break;
    }
    return result;
}

void help_prune(char** argv) {
    std::cerr << "usage: " << argv[0] << " prune [options] <graph.vg> >[output.vg]" << std::endl;
    std::cerr << "Prunes the complex regions of the graph for GCSA2 indexing. By default," << std::endl;
    std::cerr << "pruning removes the embedded paths." << std::endl;
    std::cerr << "pruning parameters:" << std::endl;
    std::cerr << "    -k, --kmer-length N    kmer length used for pruning" << std::endl;
    std::cerr << "                           (default: " << PruningParameters::KMER_LENGTH << "; " << PruningParameters::KMER_LENGTH_UNFOLDING << " with unfolding)" << std::endl;
    std::cerr << "    -e, --edge-max N       remove the edges on kmers making > N edge choices" << std::endl;
    std::cerr << "                           (default: " << PruningParameters::EDGE_MAX << "; " << PruningParameters::EDGE_MAX_UNFOLDING << " with unfolding)" << std::endl;
    std::cerr << "    -s, --subgraph-min N   remove subgraphs of < N bases" << std::endl;
    std::cerr << "                           (default: " << PruningParameters::SUBGRAPH_MIN << "; " << PruningParameters::SUBGRAPH_MIN_UNFOLDING << " with unfolding)" << std::endl;
    std::cerr << "pruning modes (-P, -r, and -u are mutually exclusive):" << std::endl;
    std::cerr << "    -P, --preserve-paths   preserve the non-alt paths and the edges on them" << std::endl;
    std::cerr << "    -r, --restore-paths    restore the edges on XG paths (requires -x)" << std::endl;
    std::cerr << "    -u, --unfold-paths     unfold XG paths (requires -x)" << std::endl;
    std::cerr << "    -x, --xg-name FILE     use this XG index" << std::endl;
    std::cerr << "    -v, --verify-paths     verify that the path exist after pruning" << std::endl;
    std::cerr << "unfolding options:" << std::endl;
    std::cerr << "    -g, --gbwt-name FILE   also unfold GBWT threads" << std::endl;
    std::cerr << "    -m, --mapping FILE     store the node mapping for duplicates in this file" << std::endl;
    std::cerr << "    -a, --append-mapping   append to the existing node mapping (requires -m)" << std::endl;
    std::cerr << "other options:" << std::endl;
    std::cerr << "    -p, --progress         show progress" << std::endl;
    std::cerr << "    -t, --threads N        use N threads (default: " << omp_get_max_threads() << ")" << std::endl;
    std::cerr << "    -d, --dry-run          determine the validity of the parameter combination" << std::endl;
}

int main_prune(int argc, char** argv) {

    if (argc == 2) {
        help_prune(argv);
        return 1;
    }

    // Command-line options.
    int kmer_length = PruningParameters::KMER_LENGTH;
    int edge_max = PruningParameters::EDGE_MAX;
    size_t subgraph_min = PruningParameters::SUBGRAPH_MIN;
    PruningMode mode = mode_prune;
    int threads = omp_get_max_threads();
    bool verify_paths = false, append_mapping = false, show_progress = false, dry_run = false;
    std::string vg_name, gbwt_name, xg_name, mapping_name;

    // Derived variables.
    bool kmer_length_set = false, edge_max_set = false, subgraph_min_set = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "kmer-length", required_argument, 0, 'k' },
            { "edge-max", required_argument, 0, 'e' },
            { "subgraph-min", required_argument, 0, 's' },
            { "preserve-paths", no_argument, 0, 'P' },
            { "restore-paths", no_argument, 0, 'r' },
            { "unfold-paths", no_argument, 0, 'u' },
            { "xg-name", required_argument, 0, 'x' },
            { "verify-paths", no_argument, 0, 'v' },
            { "gbwt-name", required_argument, 0, 'g' },
            { "mapping", required_argument, 0, 'm' },
            { "append-mapping", no_argument, 0, 'a' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { "dry-run", no_argument, 0, 'd' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "k:e:s:Prux:vg:m:apt:d", long_options, &option_index);
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
            mode = mode_preserve;
            break;
        case 'r':
            mode = mode_restore;
            break;
        case 'u':
            mode = mode_unfold;
            break;
        case 'x':
            xg_name = optarg;
            break;
        case 'v':
            verify_paths = true;
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
        case 'd':
            dry_run = true;
            break;

        case 'h':
        case '?':
            help_prune(argv);
            return 1;
        default:
            std::abort();
        }
    }
    if (optind >= argc) {
        std::cerr << "[vg prune]: graph name was not specified" << std::endl;
        return 1;
    }
    vg_name = argv[optind];
    if (!(kmer_length > 0 && edge_max > 0)) {
        std::cerr << "[vg prune]: --kmer-length and --edge-max must be positive" << std::endl;
        return 1;
    }
    if (append_mapping && mapping_name.empty()) {
        std::cerr << "[vg prune]: parameter --append-mapping requires --mapping" << std::endl;
        return 1;
    }

    // Mode-specific checks.
    if (mode == mode_prune) {
        if (verify_paths) {
            std::cerr << "[vg prune]: mode " << mode_name(mode) << " does not have paths to verify" << std::endl;
            return 1;
        }
        if (!(xg_name.empty() && gbwt_name.empty() && mapping_name.empty())) {
            std::cerr << "[vg prune]: mode " << mode_name(mode) << " does not use additional files" << std::endl;
            return 1;
        }
    }
    if (mode == mode_preserve) {
        if (!(xg_name.empty() && gbwt_name.empty() && mapping_name.empty())) {
            std::cerr << "[vg prune]: mode " << mode_name(mode) << " does not use additional files" << std::endl;
            return 1;
        }
    }
    if (mode == mode_restore) {
        if (xg_name.empty()) {
            std::cerr << "[vg prune]: mode " << mode_name(mode) << " requires --xg-name" << std::endl;
            return 1;
        }
        if (!(gbwt_name.empty() && mapping_name.empty())) {
            std::cerr << "[vg prune]: mode " << mode_name(mode) << " does not use --gbwt-name or --mapping-name" << std::endl;
            return 1;
        }
    }
    if (mode == mode_unfold) {
        if (xg_name.empty()) {
            std::cerr << "[vg prune]: mode " << mode_name(mode) << " requires --xg-name" << std::endl;
        }
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

    // Dry run.
    if (dry_run) {
        std::cerr << "Pruning mode:   " << mode_name(mode) << std::endl;
        std::cerr << "Parameters:     --kmer-length " << kmer_length << " --edge-max " << edge_max << " --subgraph-min " << subgraph_min << std::endl;
        std::cerr << "Options:        --threads " << omp_get_max_threads();
        if (verify_paths) {
            std::cerr << " --verify-paths";
        }
        if (append_mapping) {
            std::cerr << " --append_mapping";
        }
        if (show_progress) {
            std::cerr << " --progress";
        }
        if (dry_run) {
            std::cerr << " --dry-run";
        }
        std::cerr << std::endl;
        if (!vg_name.empty()) {
            std::cerr << "VG:             " << vg_name << std::endl;
        }        
        if (!xg_name.empty()) {
            std::cerr << "XG:             " << xg_name << std::endl;
        }
        if (!gbwt_name.empty()) {
            std::cerr << "GBWT:           " << gbwt_name << std::endl;
        }
        if (!mapping_name.empty()) {
            std::cerr << "Mapping:        " << mapping_name << std::endl;
        }
        return 0;
    }

    // Handle the input.
    VG* graph;
    get_input_file(vg_name, [&](std::istream& in) {
        graph = new VG(in);
    });
    vg::id_t max_node_id = graph->max_node_id();
    if (show_progress) {
        std::cerr << "Original graph " << vg_name << ": " << graph->node_count() << " nodes, " << graph->edge_count() << " edges" << std::endl;
    }

    // Remove unnecessary paths.
    if (mode == mode_preserve) {
        regex is_alt("_alt_.+_[0-9]+");
        std::set<std::string> to_remove;
        graph->paths.for_each([&](const Path& path) {
            if (regex_match(path.name(), is_alt)) {
                to_remove.insert(path.name());
            }
        });
        graph->paths.remove_paths(to_remove);
        if (show_progress) {
            std::cerr << "Removed alt paths" << std::endl;
        }
    } else {
        graph->paths.clear();
        if (show_progress) {
            std::cerr << "Removed all paths" << std::endl;
        }
    }

    // Prune the graph.
    graph->prune_complex_with_head_tail(kmer_length, edge_max, (mode == mode_preserve));
    if (show_progress) {
        std::cerr << "Pruned complex regions: "
                  << graph->node_count() << " nodes, " << graph->edge_count() << " edges" << std::endl;
    }
    graph->prune_short_subgraphs(subgraph_min);
    if (show_progress) {
        std::cerr << "Removed small subgraphs: "
                  << graph->node_count() << " nodes, " << graph->edge_count() << " edges" << std::endl;
    }

    // Preserve the VG paths.
    if (mode == mode_preserve) {
        if (verify_paths) {
            // TODO: implement
        }
    }

    // Restore the XG paths.
    if (mode == mode_restore) {
        xg::XG xg_index;
        get_input_file(xg_name, [&](std::istream& in) {
           xg_index.load(in); 
        });
        gbwt::GBWT gbwt_index;
        PhaseUnfolder unfolder(xg_index, gbwt_index, max_node_id + 1);
        unfolder.restore_paths(*graph, show_progress);
        if (verify_paths) {
            // TODO: implement
        }
    }

    // Unfold the XG paths and the GBWT threads.
    if (mode == mode_unfold) {
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
        if (verify_paths) {
            // TODO: implement
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


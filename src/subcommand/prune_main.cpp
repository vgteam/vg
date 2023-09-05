/** \file prune_main.cpp
 *
 * Defines the "vg prune" subcommand, which prunes the complex regions of the
 * graph for GCSA2 indexing.
 *
 * By default, pruning removes the nodes touched by paths of length
 * --kmer-length crossing more than --edge-max non-trivial edges. Graph
 * regions shorter than --subgraph_min are also removed. Pruning also removes
 * all embedded paths.
 *
 * For very complex graphs, there is an option to remove high-degree nodes
 * before pruning. Otherwise enumerating the k bp paths would take too long.
 *
 * With --restore-paths, the nodes and edges on non-alt paths are added back
 * after pruning.
 *
 * With --unfold-paths, pruning unfolds the non-alt paths and GBWT threads in
 * the complex regions. If --mapping is specified, the mapping from duplicate
 * node identifiers to the original identifiers is stored in a file. This file
 * can be used for building a GCSA2 index that maps to the original graph.
 */

#include "../phase_unfolder.hpp"
#include <vg/io/vpkg.hpp>
#include "subcommand.hpp"
#include "xg.hpp"
#include "../algorithms/prune.hpp"
#include "../io/save_handle_graph.hpp"

#include <gbwt/gbwt.h>

#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <string>

#include <getopt.h>
#include <omp.h>

using namespace vg;
using namespace vg::subcommand;


enum PruningMode { mode_prune, mode_restore, mode_unfold };

struct PruningParameters
{
    static std::map<PruningMode, int>    kmer_length;
    static std::map<PruningMode, int>    edge_max;
    static std::map<PruningMode, size_t> subgraph_min;
    static std::map<PruningMode, int>    max_degree;
};

std::map<PruningMode, int> PruningParameters::kmer_length { { mode_prune, 24 }, { mode_restore, 24 }, { mode_unfold, 24 } };
std::map<PruningMode, int> PruningParameters::edge_max { { mode_prune, 3 }, { mode_restore, 3 }, { mode_unfold, 3 } };
std::map<PruningMode, size_t> PruningParameters::subgraph_min { { mode_prune, 33 }, { mode_restore, 33 }, { mode_unfold, 33 } };
std::map<PruningMode, int> PruningParameters::max_degree { { mode_prune, 0 }, { mode_restore, 0 }, { mode_unfold, 0 } };

std::string mode_name(PruningMode mode) {
    std::string result = "(unknown)";
    switch (mode) {
    case mode_prune:
        result = "--prune";
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

std::string short_mode_name(PruningMode mode) {
    std::string result = "???";
    switch (mode) {
    case mode_prune:
        result = "-P";
        break;
    case mode_restore:
        result = "-r";
        break;
    case mode_unfold:
        result = "-u";
        break;
    }
    return result;
}

template<class ValueType>
void print_defaults(const std::map<PruningMode, ValueType>& defaults) {
    std::cerr << "defaults: ";
    bool first = true;
    for (auto setting : defaults) {
        if (!first) {
            std::cerr << "; ";
        }
        std::cerr << setting.second << " with " << short_mode_name(setting.first);
        first = false;
    }
    std::cerr << std::endl;
}

void help_prune(char** argv) {
    std::cerr << "usage: " << argv[0] << " prune [options] <graph.vg> >[output.vg]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Prunes the complex regions of the graph for GCSA2 indexing. Pruning the graph" << std::endl;
    std::cerr << "removes embedded paths." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Pruning parameters:" << std::endl;
    std::cerr << "    -k, --kmer-length N    kmer length used for pruning" << std::endl;
    std::cerr << "                           "; print_defaults(PruningParameters::kmer_length);
    std::cerr << "    -e, --edge-max N       remove the edges on kmers making > N edge choices" << std::endl;
    std::cerr << "                           "; print_defaults(PruningParameters::edge_max);
    std::cerr << "    -s, --subgraph-min N   remove subgraphs of < N bases" << std::endl;
    std::cerr << "                           "; print_defaults(PruningParameters::subgraph_min);
    std::cerr << "    -M, --max-degree N     if N > 0, remove nodes with degree > N before pruning" << std::endl;
    std::cerr << "                           "; print_defaults(PruningParameters::max_degree);
    std::cerr << std::endl;
    std::cerr << "Pruning modes (-P, -r, and -u are mutually exclusive):" << std::endl;
    std::cerr << "    -P, --prune            simply prune the graph (default)" << std::endl;
    std::cerr << "    -r, --restore-paths    restore the edges on non-alt paths" << std::endl;
    std::cerr << "    -u, --unfold-paths     unfold non-alt paths and GBWT threads" << std::endl;
    std::cerr << "    -v, --verify-paths     verify that the paths exist after pruning" << std::endl;
    std::cerr << "                           (potentially very slow)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Unfolding options:" << std::endl;
    std::cerr << "    -g, --gbwt-name FILE   unfold the threads from this GBWT index" << std::endl;
    std::cerr << "    -m, --mapping FILE     store the node mapping for duplicates in this file (required with -u)" << std::endl;
    std::cerr << "    -a, --append-mapping   append to the existing node mapping" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "    -p, --progress         show progress" << std::endl;
    std::cerr << "    -t, --threads N        use N threads (default: " << omp_get_max_threads() << ")" << std::endl;
    std::cerr << "    -d, --dry-run          determine the validity of the combination of options" << std::endl;
    std::cerr << std::endl;
}

int main_prune(int argc, char** argv) {

    if (argc == 2) {
        help_prune(argv);
        return 1;
    }

    // Command-line options.
    int kmer_length = 0;
    int edge_max = 0;
    size_t subgraph_min = 0;
    int max_degree = 0;
    PruningMode mode = mode_prune;
    int threads = omp_get_max_threads();
    bool verify_paths = false, append_mapping = false, show_progress = false, dry_run = false;
    std::string vg_name, gbwt_name, mapping_name;

    // Derived variables.
    bool kmer_length_set = false, edge_max_set = false, subgraph_min_set = false, max_degree_set = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "kmer-length", required_argument, 0, 'k' },
            { "edge-max", required_argument, 0, 'e' },
            { "subgraph-min", required_argument, 0, 's' },
            { "max-degree", required_argument, 0, 'M' },
            { "prune", no_argument, 0, 'P' },
            { "restore-paths", no_argument, 0, 'r' },
            { "unfold-paths", no_argument, 0, 'u' },
            { "verify-paths", no_argument, 0, 'v' },
            { "xg-name", no_argument, 0, 'x' }, // no longer needed
            { "gbwt-name", required_argument, 0, 'g' },
            { "mapping", required_argument, 0, 'm' },
            { "append-mapping", no_argument, 0, 'a' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { "dry-run", no_argument, 0, 'd' },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "k:e:s:M:Pruvx:g:m:apt:dh", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'k':
            kmer_length = parse<int>(optarg);
            kmer_length_set = true;
            break;
        case 'e':
            edge_max = parse<int>(optarg);
            edge_max_set = true;
            break;
        case 's':
            subgraph_min = parse<size_t>(optarg);
            subgraph_min_set = true;
            break;
        case 'M':
            max_degree = parse<int>(optarg);
            max_degree_set = true;
            break;
        case 'P':
            mode = mode_prune;
            break;
        case 'r':
            mode = mode_restore;
            break;
        case 'u':
            mode = mode_unfold;
            break;
        case 'v':
            verify_paths = true;
            break;
        case 'x': // no longer needed
            std::cerr << "warning: [vg prune] option --xg-name is no longer needed" << std::endl;
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
            threads = parse<int>(optarg);
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
    
    if (optind < argc) {
        // There's an input file specified.
        vg_name = get_input_file_name(optind, argc, argv);
    } else {
        // Assume they want stdin
        vg_name = "-";
    }
    
    if (!kmer_length_set) {
        kmer_length = PruningParameters::kmer_length[mode];
    }
    if (!edge_max_set) {
        edge_max = PruningParameters::edge_max[mode];
    }
    if (!subgraph_min_set) {
        subgraph_min = PruningParameters::subgraph_min[mode];
    }
    if (!max_degree_set) {
        max_degree = PruningParameters::max_degree[mode];
    }
    if (!(kmer_length > 0 && edge_max > 0)) {
        std::cerr << "error: [vg prune] --kmer-length and --edge-max must be positive" << std::endl;
        return 1;
    }

    // Mode-specific checks.
    if (mode == mode_prune) {
        if (verify_paths) {
            std::cerr << "error: [vg prune] mode " << mode_name(mode) << " does not have paths to verify" << std::endl;
            return 1;
        }
        if (!(gbwt_name.empty() && mapping_name.empty())) {
            std::cerr << "error: [vg prune] mode " << mode_name(mode) << " does not use additional files" << std::endl;
            return 1;
        }
    }
    if (mode == mode_restore) {
        if (!(gbwt_name.empty() && mapping_name.empty())) {
            std::cerr << "error: [vg prune] mode " << mode_name(mode) << " does not use additional files" << std::endl;
            return 1;
        }
    }
    if (mode == mode_unfold) {
        if (mapping_name.empty()) {
            std::cerr << "error: [vg prune] mode --unfold requires a node mapping file specified with --mapping" << std::endl;
            return 1;
        }
    }

    // Dry run.
    if (dry_run) {
        std::cerr << "Pruning mode:   " << mode_name(mode) << std::endl;
        std::cerr << "Parameters:     --kmer-length " << kmer_length << " --edge-max " << edge_max << " --subgraph-min " << subgraph_min << " --max-degree " << max_degree << std::endl;
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
            std::cerr << "VG:             " << (vg_name == "-" ? "(stdin)" : vg_name) << std::endl;
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
    std::unique_ptr<MutablePathDeletableHandleGraph> graph;
    graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(vg_name);
    xg::XG xg_index;
    std::unique_ptr<gbwt::GBWT> gbwt_index;
    
    
    vg::id_t max_node_id = graph->max_node_id();
    if (show_progress) {
        std::cerr << "Original graph " << vg_name << ": " << graph->get_node_count() << " nodes, " << graph->get_edge_count() << " edges" << std::endl;
    }

    // Remove the paths and build an XG index if needed.
    if (mode == mode_restore || mode == mode_unfold) {
        vector<path_handle_t> alt_path_handles;
        graph->for_each_path_handle([&](path_handle_t path_handle) {
            if (Paths::is_alt(graph->get_path_name(path_handle))) {
                alt_path_handles.push_back(path_handle);
            }
        });
        for (auto& alt_path_handle : alt_path_handles) {
            graph->destroy_path(alt_path_handle);
        }
        xg_index.from_path_handle_graph(*graph);
        if (show_progress) {
            std::cerr << "Built a temporary XG index" << std::endl;
        }
    }
    
    // Destroy all remaining paths
    vector<path_handle_t> path_handles;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        path_handles.push_back(path_handle);
    });
    for (auto path_handle : path_handles) {
        graph->destroy_path(path_handle);
    }
    
    if (show_progress) {
        std::cerr << "Removed all paths" << std::endl;
    }

    // Remove high-degree nodes.
    if (max_degree > 0) {
        algorithms::remove_high_degree_nodes(*graph, max_degree);
        if (show_progress) {
            std::cerr << "Removed high-degree nodes: "
                      << graph->get_node_count() << " nodes, " << graph->get_edge_count() << " edges" << std::endl;
        }
    }

    // Prune the graph.
    algorithms::prune_complex_with_head_tail(*graph, kmer_length, edge_max);
    if (show_progress) {
        std::cerr << "Pruned complex regions: "
                  << graph->get_node_count() << " nodes, " << graph->get_edge_count() << " edges" << std::endl;
    }
    algorithms::prune_short_subgraphs(*graph, subgraph_min);
    if (show_progress) {
        std::cerr << "Removed small subgraphs: "
                  << graph->get_node_count() << " nodes, " << graph->get_edge_count() << " edges" << std::endl;
    }

    // Restore the non-alt paths.
    if (mode == mode_restore) {
        // Make an empty GBWT index to pass along
        gbwt::GBWT empty_gbwt;
        PhaseUnfolder unfolder(xg_index, empty_gbwt, max_node_id + 1);
        unfolder.restore_paths(*graph, show_progress);
        if (verify_paths) {
            size_t failures = unfolder.verify_paths(*graph, show_progress);
            if (failures > 0) {
                std::cerr << "warning: [vg prune] verification failed for " << failures << " paths" << std::endl;
            }
        }
    }

    // Unfold the XG paths and the GBWT threads.
    if (mode == mode_unfold) {
        if (!gbwt_name.empty()) {
            get_input_file(gbwt_name, [&](std::istream& in) {
                gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(in);
                if (gbwt_index.get() == nullptr) {
                    std::cerr << "[vg prune]: could not load GBWT" << std::endl;
                    exit(1);
                }
            });
        } else {
            // The PhaseUnfolder can't deal with having no GBWT at all; we need to give it an empty one.
            gbwt_index = unique_ptr<gbwt::GBWT>(new gbwt::GBWT());
            // TODO: Let us pass in null pointers instead.
        }
        PhaseUnfolder unfolder(xg_index, *gbwt_index, max_node_id + 1);
        if (append_mapping) {
            unfolder.read_mapping(mapping_name);
        }
        unfolder.unfold(*graph, show_progress);
        if (!mapping_name.empty()) {
            unfolder.write_mapping(mapping_name);
        }
        if (verify_paths) {
            size_t failures = unfolder.verify_paths(*graph, show_progress);
            if (failures > 0) {
                std::cerr << "warning: [vg prune] verification failed for " << failures << " paths" << std::endl;
            }
        }
    }

    // Serialize.
    
    vg::io::save_handle_graph(graph.get(), std::cout);
    if (show_progress) {
        std::cerr << "Serialized the graph: "
                  << graph->get_node_count() << " nodes, " << graph->get_edge_count() << " edges" << std::endl;
    }

    return 0;
}

// Register subcommand
static Subcommand vg_prune("prune", "prune the graph for GCSA2 indexing", TOOLKIT, main_prune);


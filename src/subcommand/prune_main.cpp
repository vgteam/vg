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
 * Optionally, replace each pruned region with a set of disjoint paths
 * corresponding to all distinct XG paths and GBWT threads in that region.
 * TODO: We also need to output a mapping from duplicated node ids to original
 * node ids for GCSA2 construction.
 */

#include "subcommand.hpp"
#include "../vg.hpp"

#include <gbwt/gbwt.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <list>
#include <regex>
#include <set>
#include <stack>
#include <string>
#include <vector>

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
    std::cerr << "    -P, --preserve-paths   preserve the embedded non-alt paths in the graph" << std::endl;
    std::cerr << "    -g, --gbwt-name FILE   unfold the complex regions by using the threads in" << std::endl;
    std::cerr << "                           this GBWT index (requires -g, avoid -P)" << std::endl;
    std::cerr << "    -x, --xg_name FILE     use this XG index with -g" << std::endl;
    std::cerr << "    -m, --mapping FILE     store the node mapping from -g in this file" << std::endl;
    std::cerr << "other options:" << std::endl;
    std::cerr << "    -p, --progress         show progress" << std::endl;
    std::cerr << "    -t, --threads N        use N threads (default: " << omp_get_max_threads() << ")" << std::endl;
}

void create_state(const gbwt::GBWT& gbwt_index, vg::id_t from, bool is_reverse,
                  std::stack<std::pair<gbwt::SearchState, std::vector<gbwt::node_type>>>& intermediate) {
    gbwt::node_type node = gbwt::Node::encode(from, is_reverse);
    gbwt::SearchState state = gbwt_index.find(node);
    if (state.empty()) {
        return;
    }
    std::vector<gbwt::node_type> path { node };
    intermediate.push(std::make_pair(state, path));
}

bool extend_state(const gbwt::GBWT& gbwt_index, gbwt::SearchState state, std::vector<gbwt::node_type> path,
                  vg::id_t to, bool is_reverse,
                  std::stack<std::pair<gbwt::SearchState, std::vector<gbwt::node_type>>>& intermediate) {
    gbwt::node_type node = gbwt::Node::encode(to, is_reverse);
    state = gbwt_index.extend(state, node);
    if (state.empty()) {
        return false;
    }
    path.push_back(node);
    intermediate.push(std::make_pair(state, path));
    return true;
}

void insert_path(std::set<std::vector<gbwt::node_type>>& paths, std::vector<gbwt::node_type>& path) {
    if (path.size() < 2) {
        return;
    }
    std::vector<gbwt::node_type> reverse_complement(path.size(), 0);
    for (size_t i = 0; i < path.size(); i++) {
        reverse_complement[path.size() - 1 - i] = gbwt::Node::reverse(path[i]);
    }
    paths.insert(std::min(path, reverse_complement));
}

void list_paths(VG& component, const std::set<vg::id_t> border, const gbwt::GBWT& gbwt_index, vg::id_t from, std::set<std::vector<gbwt::node_type>>& paths) {
    std::stack<std::pair<gbwt::SearchState, std::vector<gbwt::node_type>>> intermediate;
    create_state(gbwt_index, from, false, intermediate);
    create_state(gbwt_index, from, true, intermediate);

    while (!intermediate.empty()) {
        gbwt::SearchState state;
        std::vector<gbwt::node_type> path;
        std::tie(state, path) = intermediate.top(); intermediate.pop();
        vg::id_t node = gbwt::Node::id(path.back());
        bool is_reverse = gbwt::Node::is_reverse(path.back());

        std::vector<Edge*> edges = component.edges_of(component.get_node(node));
        bool has_successors = false;
        for (Edge* edge : edges) {
            if (edge->from() == node && edge->from_start() == is_reverse) {
                has_successors |= extend_state(gbwt_index, state, path, edge->to(), edge->to_end(), intermediate);
            }
            else if (edge->to() == node && edge->to_end() != is_reverse) {
                has_successors |= extend_state(gbwt_index, state, path, edge->from(), !edge->from_start(), intermediate);
            }
        }

        if (!has_successors || border.find(node) != border.end()) {
            insert_path(paths, path);
        }
    }
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
    bool preserve_paths = false, show_progress = false;
    std::string gbwt_name, xg_name, mapping_name;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "kmer-length", required_argument, 0, 'k' },
            { "edge-max", required_argument, 0, 'e' },
            { "subgraph-min", required_argument, 0, 's' },
            { "preserve-paths", no_argument, 0, 'P' },
            { "gbwt-name", required_argument, 0, 'g' },
            { "xg-name", required_argument, 0, 'x' },
            { "mapping", required_argument, 0, 'm' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "k:e:s:Pg:x:m:pt:", long_options, &option_index);
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
        case 'P':
            preserve_paths = true;
            break;
        case 'g':
            gbwt_name = optarg;
            break;
        case 'x':
            xg_name = optarg;
            break;
        case 'm':
            mapping_name = optarg;
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
        cerr << "[vg prune]: --kmer-length and --edge-max must be positive" << endl;
        return 1;
    }
    if (!gbwt_name.empty() && xg_name.empty()) {
        cerr << "[vg prune]: parameter --gbwt-name requires --xg-name" << endl;
        return 1;
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
    if (!gbwt_name.empty()) {
        gbwt::GBWT gbwt_index;
        get_input_file(gbwt_name, [&](std::istream& in) {
           gbwt_index.load(in);
        });
        xg::XG xg_index;
        get_input_file(xg_name, [&](std::istream& in) {
           xg_index.load(in); 
        });

        // TODO: Move this to PhaseDuplicator.
        VG complement;
        for (gbwt::comp_type comp = 1; comp < gbwt_index.effective(); comp++) {
            gbwt::node_type gbwt_node = gbwt_index.toNode(comp);
            std::vector<gbwt::edge_type> outgoing = gbwt_index.edges(gbwt_node);
            for (gbwt::edge_type outedge : outgoing) {
                if (outedge.first == gbwt::ENDMARKER) {
                    continue;
                }
                Edge candidate;
                candidate.set_from(gbwt::Node::id(gbwt_node));
                candidate.set_to(gbwt::Node::id(outedge.first));
                candidate.set_from_start(gbwt::Node::is_reverse(gbwt_node));
                candidate.set_to_end(gbwt::Node::is_reverse(outedge.first));
                if (!graph->has_edge(candidate)) {
                    complement.add_node(xg_index.node(candidate.from()));
                    complement.add_node(xg_index.node(candidate.to()));
                    complement.add_edge(candidate);
                }
            }
        }
        std::list<VG> components;
        complement.disjoint_subgraphs(components);
        if (show_progress) {
            std::cerr << "Complement graph: "
                      << complement.node_count() << " nodes, " << complement.edge_count() << " edges in "
                      << components.size() << " components" << std::endl;
        }

        size_t haplotype_paths = 0;
        VG unfolded;
        for (VG& component : components) {
            std::set<vg::id_t> border;   // Nodes that exist in the component and the pruned graph.
            component.for_each_node([&](Node* node) {
                if (graph->has_node(node->id())) {
                    border.insert(node->id());
                }
            });
            std::set<std::vector<gbwt::node_type>> paths;
            for (vg::id_t start_node : border) {
                list_paths(component, border, gbwt_index, start_node, paths);
            }
            haplotype_paths += paths.size();
            for (const std::vector<gbwt::node_type>& path : paths) {
                Node prev = xg_index.node(gbwt::Node::id(path.front()));
                unfolded.add_node(prev);
                for (size_t i = 1; i < path.size(); i++) {
                    Node node = xg_index.node(gbwt::Node::id(path[i]));
                    if (i + 1 < path.size()) {
                        node.set_id(max_node_id + 1); max_node_id++; // TODO: We need to store the mapping.
                    }
                    unfolded.add_node(node);
                    Edge edge;
                    edge.set_from(prev.id());
                    edge.set_to(node.id());
                    edge.set_from_start(gbwt::Node::is_reverse(path[i - 1]));
                    edge.set_to_end(gbwt::Node::is_reverse(path[i]));
                    unfolded.add_edge(edge);
                    prev = node;
                }
            }
        }
        if (show_progress) {
            std::cerr << "Unfolded graph: "
                      << unfolded.node_count() << " nodes, " << unfolded.edge_count() << " edges on "
                      << haplotype_paths << " paths" << std::endl;
        }

        graph->extend(unfolded);
    }

    // Serialize.
    graph->serialize_to_ostream(std::cout);
    if (show_progress) {
        std::cerr << "Serialized the graph: "
                  << graph->node_count() << " nodes, " << graph->edge_count() << " edges" << std::endl;
    }
    if (!mapping_name.empty()) {
        std::cerr << "Not implemented yet!" << std::endl;
        // (write mapping)
    }

    std::cerr << "Max node: " << graph->max_node_id() << ", next node: " << graph->current_id << std::endl;

    delete graph; graph = nullptr;

    return 0;
}

// Register subcommand
static Subcommand vg_prune("prune", "prune the graph for GCSA2 indexing", TOOLKIT, main_prune);


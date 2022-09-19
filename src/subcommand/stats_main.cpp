/** \file stats_main.cpp
 *
 * Defines the "vg stats" subcommand, which evaluates graphs and alignments.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>

#include "subcommand.hpp"
#include "../algorithms/distance_to_head.hpp"
#include "../algorithms/distance_to_tail.hpp"
#include "../handle.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../annotation.hpp"
#include "../snarl_distance_index.hpp"

#include "../path.hpp"
#include "../statistics.hpp"
#include "../genotypekit.hpp"

#include "xg.hpp"
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"
#include "../io/converted_hash_graph.hpp"
#include "../io/save_handle_graph.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;
using namespace vg::algorithms;

void help_stats(char** argv) {
    cerr << "usage: " << argv[0] << " stats [options] [<graph file>]" << endl
         << "options:" << endl
         << "    -z, --size             size of graph" << endl
         << "    -N, --node-count       number of nodes in graph" << endl
         << "    -E, --edge-count       number of edges in graph" << endl
         << "    -l, --length           length of sequences in graph" << endl
         << "    -L, --self-loops       number of self-loops" << endl
         << "    -s, --subgraphs        describe subgraphs of graph" << endl
         << "    -H, --heads            list the head nodes of the graph" << endl
         << "    -T, --tails            list the tail nodes of the graph" << endl
         << "    -e, --nondeterm        list the nondeterministic edge sets" << endl
         << "    -c, --components       print the strongly connected components of the graph" << endl
         << "    -A, --is-acyclic       print if the graph is acyclic or not" << endl
         << "    -n, --node ID          consider node with the given id" << endl
         << "    -d, --to-head          show distance to head for each provided node" << endl
         << "    -t, --to-tail          show distance to head for each provided node" << endl
         << "    -a, --alignments FILE  compute stats for reads aligned to the graph" << endl
         << "    -r, --node-id-range    X:Y where X and Y are the smallest and largest "
        "node id in the graph, respectively" << endl
         << "    -o, --overlap PATH    for each overlapping path mapping in the graph write a table:" << endl
         << "                              PATH, other_path, rank1, rank2" << endl
         << "                          multiple allowed; limit comparison to those provided" << endl
         << "    -O, --overlap-all     print overlap table for the cartesian product of paths" << endl
         << "    -R, --snarls          print statistics for each snarl" << endl
         << "    -F, --format          graph format from {VG-Protobuf, PackedGraph, HashGraph, ODGI, XG}. " <<
        "Can't detect Protobuf if graph read from stdin" << endl
         << "    -D, --degree-dist     print degree distribution of the graph." << endl
         << "    -b, --dist-snarls FILE print the sizes and depths of the snarls in a given distance index." << endl
         << "    -p, --threads N       number of threads to use [all available]" << endl
         << "    -v, --verbose         output longer reports" << endl;
}

int main_stats(int argc, char** argv) {

    if (argc == 2) {
        help_stats(argv);
        return 1;
    }

    bool stats_size = false;
    bool stats_length = false;
    bool stats_self_loops = false;
    bool stats_subgraphs = false;
    bool stats_heads = false;
    bool stats_tails = false;
    bool stats_nondeterm = false;
    bool show_sibs = false;
    bool show_components = false;
    bool head_distance = false;
    bool tail_distance = false;
    bool node_count = false;
    bool edge_count = false;
    bool verbose = false;
    bool is_acyclic = false;
    bool stats_range = false;
    set<vg::id_t> ids;
    // What alignments GAM file should we read and compute stats on with the
    // graph?
    string alignments_filename;
    vector<string> paths_to_overlap;
    bool overlap_all_paths = false;
    bool snarl_stats = false;
    bool format = false;
    bool degree_dist = false;
    string distance_index_filename;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"size", no_argument, 0, 'z'},
            {"node-count", no_argument, 0, 'N'},
            {"edge-count", no_argument, 0, 'E'},
            {"length", no_argument, 0, 'l'},
            {"self-loops", no_argument, 0, 'L'},
            {"subgraphs", no_argument, 0, 's'},
            {"heads", no_argument, 0, 'H'},
            {"tails", no_argument, 0, 'T'},
            {"nondeterm", no_argument, 0, 'e'},
            {"help", no_argument, 0, 'h'},
            {"components", no_argument, 0, 'c'},
            {"to-head", no_argument, 0, 'd'},
            {"to-tail", no_argument, 0, 't'},
            {"node", required_argument, 0, 'n'},
            {"alignments", required_argument, 0, 'a'},
            {"is-acyclic", no_argument, 0, 'A'},
            {"node-id-range", no_argument, 0, 'r'},
            {"verbose", no_argument, 0, 'v'},
            {"overlap", no_argument, 0, 'o'},
            {"overlap-all", no_argument, 0, 'O'},
            {"snarls", no_argument, 0, 'R'},
            {"format", no_argument, 0, 'F'},
            {"degree-dist", no_argument, 0, 'D'},
            {"dist-snarls", required_argument, 0, 'b'},
            {"threads", required_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hzlLsHTecdtn:NEa:vAro:ORFDb:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'z':
            stats_size = true;
            break;

        case 'N':
            node_count = true;
            break;

        case 'E':
            edge_count = true;
            break;

        case 'l':
            stats_length = true;
            break;

        case 'L':
            stats_self_loops = true;
            break;

        case 's':
            stats_subgraphs = true;
            break;

        case 'H':
            stats_heads = true;
            break;

        case 'T':
            stats_tails = true;
            break;

        case 'e':
            stats_nondeterm = true;
            break;

        case 'S':
            show_sibs = true;
            break;

        case 'c':
            show_components = true;
            break;

        case 'd':
            head_distance = true;
            break;

        case 't':
            tail_distance = true;
            break;

        case 'n':
            ids.insert(parse<vg::id_t>(optarg));
            break;

        case 'A':
            is_acyclic = true;
            break;

        case 'a':
            alignments_filename = optarg;
            break;

        case 'r':
            stats_range = true;
            break;

        case 'o':
            paths_to_overlap.push_back(optarg);
            break;

        case 'O':
            overlap_all_paths = true;
            break;
            
        case 'R':
            snarl_stats = true;
            break;

        case 'v':
            verbose = true;
            break;

        case 'F':
            format = true;
            break;

        case 'D':
            degree_dist = true;
            break;
        case 'b':
            distance_index_filename = optarg;
            break;
        case 'p':
        {
            int num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg stats] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }

        case 'h':
        case '?':
            help_stats(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    unique_ptr<PathHandleGraph> graph;
    string graph_file_name;
    if (have_input_file(optind, argc, argv)) {
        // We have an (optional, because we can just process alignments) graph input file.
        // TODO: we can load any PathHandleGraph, but some operations still require a VG
        // In those cases, we convert back to vg::VG
        graph_file_name = get_input_file_name(optind, argc, argv);
        graph = vg::io::VPKG::load_one<PathHandleGraph>(graph_file_name);
    }
    
    // We have function to make sure the graph was passed and complain if not
    auto require_graph = [&graph]() {
        if (graph.get() == nullptr) {
            cerr << "error[vg stats]: The selected operation requires passing a graph file to work on" << endl;
            exit(1);
        }
    };
    

    if (stats_size) {
        require_graph();
        cout << "nodes" << "\t" << graph->get_node_count() << endl
            << "edges" << "\t" << graph->get_edge_count() << endl;
    }

    if (node_count) {
        require_graph();
        cout << graph->get_node_count() << endl;
    }

    if (edge_count) {
        require_graph();
        cout << graph->get_edge_count() << endl;
    }

    if (stats_length) {
        require_graph();
        cout << "length" << "\t" << graph->get_total_length() << endl;
    }

    if (stats_self_loops) {
        require_graph();
        size_t total = 0;
        graph->for_each_edge([&](const edge_t& edge) {
            if (graph->get_id(edge.first) == graph->get_id(edge.second)) {
                total++;
            }
        });
        cout << "self-loops" << "\t" << total << endl;
    }

    if (stats_heads) {
        require_graph();
        vector<handle_t> heads = handlealgs::head_nodes(graph.get());
        cout << "heads" << "\t";
        for (auto& h : heads) {
            cout << graph->get_id(h) << " ";
        }
        cout << endl;
    }

    if (stats_tails) {
        require_graph();
        vector<handle_t> tails = handlealgs::tail_nodes(graph.get());
        cout << "tails" << "\t";
        for (auto& t : tails) {
            cout << graph->get_id(t) << " ";
        }
        cout << endl;
    }

    if (stats_nondeterm) {
        require_graph();
        graph->for_each_handle([&](const handle_t& handle) {
            nid_t id = graph->get_id(handle);
            for (bool is_reverse : { false, true }) {
                std::map<char, std::vector<handle_t>> edges;
                graph->follow_edges(graph->get_handle(id, is_reverse), false, [&](const handle_t& to) {
                    edges[graph->get_base(to, 0)].push_back(to);
                });
                for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
                    if (iter->second.size() > 1) {
                        std::cout << "nondeterministic\t" << id << (is_reverse ? "-" : "+");
                        for (const handle_t& to : iter->second) {
                            std::cout << "\t" << graph->get_id(to) << (graph->get_is_reverse(to) ? "-" : "+");
                        }
                        std::cout << std::endl;
                    }
                }
            }
        });
    }

    if (stats_subgraphs) {
        require_graph();
        
        // TODO: Pretty sure "subgraphs" means "weakly connected components",
        // but this isn't really explained.
        
        vector<pair<unordered_set<nid_t>, vector<handle_t>>> subgraphs_with_tips =
            handlealgs::weakly_connected_components_with_tips(graph.get());
        
        for (auto& subgraph_and_tips : subgraphs_with_tips) {
            // For each subgraph set and its inward tip handles
            auto& subgraph = subgraph_and_tips.first;
            auto& tips = subgraph_and_tips.second;
            
            // Decide if we need a comma before us or not
            bool first = true;
            for (handle_t& tip : tips) {
                // Print all the IDs of heads
                if (graph->get_is_reverse(tip)) {
                    // Heads are locally forward, so this isn't one.
                    // TODO: subgraphs with only tails get no identification.
                    continue;
                }
                if (!first) {
                    cout << ",";
                } else {
                    first = false;
                }
                cout << graph->get_id(tip);
            }
            cout << "\t";
            
            // Now we need the total length. TODO: can we do a batch lookup?
            size_t total_length = 0;
            for (auto& id : subgraph) {
                total_length += graph->get_length(graph->get_handle(id));
            }
            
            cout << total_length << endl;
        }
    }

    if (stats_range) {
        require_graph();
        cout << "node-id-range\t" << graph->min_node_id() << ":" << graph->max_node_id() << endl;
    }

    if (show_components) {
        require_graph();
        for (auto& c : handlealgs::strongly_connected_components(graph.get())) {
            for (auto& id : c) {
                cout << id << ", ";
            }
            cout << endl;
        }
    }

    if (is_acyclic) {
        require_graph();
        if (handlealgs::is_acyclic(graph.get())) {
            cout << "acyclic" << endl;
        } else {
            cout << "cyclic" << endl;
        }
    }

    if (head_distance) {
        require_graph();
        for (auto id : ids) {
            auto n = graph->get_handle(id, false);
            cout << id << " to head:\t"
                 << distance_to_head(n, 1000, graph.get()) << endl;
        }
    }

    if (tail_distance) {
        require_graph();
        for (auto id : ids) {
            auto n = graph->get_handle(id, false);
            cout << id << " to tail:\t"
                << distance_to_tail(n, 1000, graph.get()) << endl;
        }
    }

    if (format) {
        require_graph();
        string format_string;
        if (dynamic_cast<xg::XG*>(graph.get()) != nullptr) {
            format_string = "XG";
        } else if (dynamic_cast<GFAHandleGraph*>(graph.get()) != nullptr) {
            // important this check comes before PackedGraph
            format_string = "GFA";
        } else if (dynamic_cast<bdsg::PackedGraph*>(graph.get()) != nullptr) {
            format_string = "PackedGraph";
        } else if (dynamic_cast<vg::io::ConvertedHashGraph*>(graph.get()) != nullptr) {
            // Was Protobuf but we're using a HashGraph internally
            format_string = "VG-Protobuf";
        } else if (dynamic_cast<bdsg::HashGraph*>(graph.get()) != nullptr) {
            format_string = "HashGraph";
        } else if (dynamic_cast<bdsg::ODGI*>(graph.get()) != nullptr) {
            format_string = "ODGI";
        } else {
            format_string = "Unknown";
        }
        cout << "format: " << format_string << endl;
    }

    if (degree_dist) {
        require_graph();
        // compute degrees
        map<size_t, tuple<size_t, size_t, size_t, size_t>> degree_to_count;
        graph->for_each_handle([&degree_to_count, &graph](handle_t handle) {
                size_t left_degree = graph->get_degree(handle, true);
                size_t right_degree = graph->get_degree(handle, false);
                // update sides count
                ++get<0>(degree_to_count[left_degree]);
                ++get<0>(degree_to_count[right_degree]);
                // update min count
                ++get<1>(degree_to_count[std::min(left_degree, right_degree)]);
                // update max count
                ++get<2>(degree_to_count[std::max(left_degree, right_degree)]);
                // update total count
                ++get<3>(degree_to_count[left_degree + right_degree]);
            });
        // print degrees
        cout << "Degree\tSides\tNodes(min)\tNodes(max)\tNodes(total)" << endl;
        for (const auto& dg : degree_to_count) {
            cout << dg.first << "\t" << get<0>(dg.second) << "\t" <<get<1>(dg.second) << "\t" << get<2>(dg.second) << "\t" << get<3>(dg.second) << endl;
        }
    }

    if (!paths_to_overlap.empty() || overlap_all_paths) {
        require_graph();
        
        VG* vg_graph = dynamic_cast<VG*>(graph.get());
        if (vg_graph == nullptr) {
            // TODO: This path overlap code can be handle-ified, and should be.
            vg_graph = new vg::VG();
            handlealgs::copy_path_handle_graph(graph.get(), vg_graph);
            // Give the unique_ptr ownership and delete the graph we loaded.
            graph.reset(vg_graph);
            // Make sure the paths are all synced up
            vg_graph->paths.to_graph(vg_graph->graph);
        }
    
        auto cb = [&](const Path& p1, const Path& p2) {
            // sparse storage of the correspondence matrix
            // map from ranks in first to ranks in second
            map<vg::id_t, vector<int64_t> > p1_ranks;
            map<vg::id_t, vector<int64_t> > p2_ranks;
            //map<int64_t, int64_t> relative_rank;
            for (int64_t i = 0; i < p1.mapping_size(); ++i) {
                auto& mapping = p1.mapping(i);
                p1_ranks[mapping.position().node_id()].push_back(mapping.rank());
            }
            for (int64_t i = 0; i < p2.mapping_size(); ++i) {
                auto& mapping = p2.mapping(i);
                p2_ranks[mapping.position().node_id()].push_back(mapping.rank());
            }
            // intersect
            set<pair<int64_t, int64_t> > seen;
            for (auto& p : p1_ranks) {
                auto f = p2_ranks.find(p.first);
                if (f != p2_ranks.end()) {
                    for (auto& id1 : p.second) {
                        for (auto& id2 : f->second) {
                            if (seen.count(make_pair(id1, id2))) continue;
                            cout << p1.name() << "____" << p2.name() << "\t" << id1 << "\t" << id2 << endl;
                            seen.insert(make_pair(id1, id2));
                        }
                    }
                } else {
                    // non-overlap
                    for (auto& id1 : p.second) {
                        if (seen.count(make_pair(id1, 0))) continue;
                        cout << p1.name() << "____" << p2.name() << "\t" << id1 << "\t" << 0 << endl;
                        seen.insert(make_pair(id1, 0));
                    }
                }
            }
            // get the non-overlapping bits of the second sequence
            for (auto& p : p2_ranks) {
                auto f = p1_ranks.find(p.first);
                if (f != p1_ranks.end()) {
                    for (auto& id2 : p.second) {
                        if (seen.count(make_pair(0, id2))) continue;
                        cout << p1.name() << "____" << p2.name() << "\t" << 0 << "\t" << id2 << endl;
                        seen.insert(make_pair(0, id2));
                    }
                }
            }

        };
        cout << "comparison" << "\t" << "x" << "\t" << "y" << endl;
        vector<string> path_names;
        if (overlap_all_paths) {
            path_names = vg_graph->paths.all_path_names();
        } else {
            path_names = paths_to_overlap;
        }
        for (auto& p1_name : path_names) {
            Path p1 = vg_graph->paths.path(p1_name);
            for (auto& p2_name : path_names) {
                if (p1_name == p2_name) {
                    continue;
                }
                Path p2 = vg_graph->paths.path(p2_name);
                cb(p1, p2);
            }
        }
    }

    if (!alignments_filename.empty()) {
        // Read in the given GAM
        ifstream alignment_stream(alignments_filename);

        // We need some allele parsing functions

        // This one gets the site name from an allele path name
        auto path_name_to_site = [](const string& path_name) -> string {
            auto last_underscore = path_name.rfind('_');
            assert(last_underscore != string::npos);
            return path_name.substr(0, last_underscore);
        };

        // This one gets the allele name from an allele path name
        auto path_name_to_allele = [](const string& path_name) -> string {
            auto last_underscore = path_name.rfind('_');
            assert(last_underscore != string::npos);
            return path_name.substr(last_underscore + 1);
        };
        
        // In order to do stats across multiple threads, we define these add-able bundles of stats.
        struct ReadStats {
            // These are the general stats we will compute.
            size_t total_alignments = 0;
            size_t total_aligned = 0;
            size_t total_primary = 0;
            size_t total_secondary = 0;
            size_t total_perfect = 0; // Number of reads with no indels or substitutions relative to their paths
            size_t total_gapless = 0; // Number of reads with no indels relative to their paths

            // These are for tracking which nodes are covered and which are not
            map<vg::id_t, size_t> node_visit_counts;

            // And for counting indels
            // Inserted bases also counts softclips
            size_t total_insertions = 0;
            size_t total_inserted_bases = 0;
            size_t total_deletions = 0;
            size_t total_deleted_bases = 0;
            // And substitutions
            size_t total_substitutions = 0;
            size_t total_substituted_bases = 0;
            // And softclips
            size_t total_softclips = 0;
            size_t total_softclipped_bases = 0;
            // And pairing
            size_t total_paired = 0;
            size_t total_proper_paired = 0;

            // In verbose mode we want to report details of insertions, deletions,
            // and substitutions, and soft clips.
            vector<pair<vg::id_t, Edit>> insertions;
            vector<pair<vg::id_t, Edit>> deletions;
            vector<pair<vg::id_t, Edit>> substitutions;
            vector<pair<vg::id_t, Edit>> softclips;
            
            // This is going to be indexed by site
            // ("_alt_f6d951572f9c664d5d388375aa8b018492224533") and then by allele
            // ("0"). A read only counts if it visits a node that's on one allele
            // and not any others in that site.
            map<string, map<string, size_t>> reads_on_allele;
        
            inline ReadStats& operator+=(const ReadStats& other) {
                total_alignments += other.total_alignments;
                total_aligned += other.total_aligned;
                total_primary += other.total_primary;
                total_secondary += other.total_secondary;
                total_perfect += other.total_perfect;
                total_gapless += other.total_gapless;
                
                for (auto& kv : other.node_visit_counts) {
                    node_visit_counts[kv.first] += kv.second;
                }
                
                total_insertions += other.total_insertions;
                total_inserted_bases += other.total_inserted_bases;
                total_deletions += other.total_deletions;
                total_deleted_bases += other.total_deleted_bases;
                total_substitutions += other.total_substitutions;
                total_substituted_bases += other.total_substituted_bases;
                total_softclips += other.total_softclips;
                total_softclipped_bases += other.total_softclipped_bases;
                total_paired += other.total_paired;
                total_proper_paired += other.total_proper_paired;
                
                std::copy(other.insertions.begin(), other.insertions.end(), std::back_inserter(insertions));
                std::copy(other.deletions.begin(), other.deletions.end(), std::back_inserter(deletions));
                std::copy(other.substitutions.begin(), other.substitutions.end(), std::back_inserter(substitutions));
                std::copy(other.softclips.begin(), other.softclips.end(), std::back_inserter(softclips));
                
                for (auto& kv : other.reads_on_allele) {
                    auto& dest = reads_on_allele[kv.first];
                    for (auto& kv2 : kv.second) {
                        dest[kv2.first] += kv2.second;
                    }
                }
                
                return *this;
            }
        };

        // Before we go over the reads, we need to make a map that tells us what
        // nodes are unique to what allele paths. Stores site and allele parts
        // separately.
        map<vg::id_t, pair<string, string>> allele_path_for_node;

        // Create a combined ReadStats accumulator. We need to pre-populate its
        // reads_on_allele with 0s when we look at the alleles so we know which
        // sites actually have 2 alleles and which only have 1 in the graph.
        ReadStats combined;

        if (graph.get() != nullptr) {
            // We have a graph to work on

            // For each pair of allele paths in the graph, we need to find out
            // whether the coverage imbalance between them among primary alignments
            // is statistically significant. For this, we need to track how many
            // reads overlap the distinct parts of allele paths.

            graph->for_each_handle([&](handle_t node) {
                // For every node in parallel

                // We want a unique allele path on it
                string allele_path;
                
                graph->for_each_step_on_handle(node, [&](const step_handle_t& step) -> bool {
                    // Get the name of every patht hat goes here (some may repeat)
                    auto path_name = graph->get_path_name(graph->get_path_handle_of_step(step));
                    
                    if(Paths::is_alt(path_name) && path_name != allele_path) {
                        // If it's a new/distinct allele path
                        if(allele_path.empty()) {
                            // It's the first. Take it.
                            allele_path = path_name;
                            // Check for more overlappin alt paths
                            return true;
                        } else {
                            // It's a subsequent one. This node is not uniquely part
                            // of any allele path. So we want to skip the node.
                            allele_path.clear();
                            return false;
                        }
                    }
                    
                    // If not an alt, keep going
                    return true;
                });
                
                if (allele_path.empty()) {
                    // We did not find a unique overlapping allele path.
                    // Skip the node.
                    return;
                }
                
                // We found an allele path for this node

                // Get its site and allele so we can count it as a biallelic
                // site. Note that sites where an allele has no unique nodes
                // (pure indels, for example) can't be handled and will be
                // ignored.
                auto site = path_name_to_site(allele_path);
                auto allele = path_name_to_allele(allele_path);


                #pragma omp critical (allele_path_for_node)
                allele_path_for_node[graph->get_id(node)] = make_pair(site, allele);

                #pragma omp critical (reads_on_allele)
                combined.reads_on_allele[site][allele] = 0;
            }, true);
        }

        // Allocate per-thread storage for stats
        size_t thread_count = vg::get_thread_count();
        vector<ReadStats> read_stats;
        read_stats.resize(thread_count); 

        // when we get each read, process it into the current thread's stats
        function<void(Alignment&)> lambda = [&](Alignment& aln) {
            int tid = omp_get_thread_num();
            auto& stats = read_stats.at(tid);
            // We ought to be able to do many stats on the alignments.

            // Now do all the non-mapping stats
            stats.total_alignments++;
            if(aln.is_secondary()) {
                stats.total_secondary++;
            } else {
                stats.total_primary++;
                bool has_alignment = aln.score() > 0;
                if (has_alignment) {
                    // We only count aligned primary reads in "total aligned";
                    // the primary can't be unaligned if the secondary is
                    // aligned.
                    stats.total_aligned++;
                }
                
                if (aln.has_fragment_next() || aln.has_fragment_prev() || has_annotation(aln, "proper_pair")) {
                    stats.total_paired++;
                    if (has_annotation(aln, "proper_pair") && get_annotation<bool>(aln, "proper_pair")) {
                        stats.total_proper_paired++;
                    }
                }

                // Which sites and alleles does this read support. TODO: if we hit
                // unique nodes from multiple alleles of the same site, we should...
                // do something. Discard the read? Not just count it on both sides
                // like we do now.
                set<pair<string, string>> alleles_supported;
                
                // We check if the read has non-softclip indels, or any edits at all.
                bool has_non_match_edits = false;
                bool has_non_softclip_indel_edits = false;

                for(size_t i = 0; i < aln.path().mapping_size(); i++) {
                    // For every mapping...
                    auto& mapping = aln.path().mapping(i);
                    vg::id_t node_id = mapping.position().node_id();

                    if(allele_path_for_node.count(node_id)) {
                        // We hit a unique node for this allele. Add it to the set,
                        // in case we hit another unique node for it later in the
                        // read.
                        alleles_supported.insert(allele_path_for_node.at(node_id));
                    }

                    // Record that there was a visit to this node.
                    stats.node_visit_counts[node_id]++;

                    for(size_t j = 0; j < mapping.edit_size(); j++) {
                        // Go through edits and look for each type.
                        auto& edit = mapping.edit(j);

                        if(edit.to_length() > edit.from_length()) {
                            // This is an insert or softclip and not a match
                            has_non_match_edits = true;
                            if((j == 0 && i == 0) || (j == mapping.edit_size() - 1 && i == aln.path().mapping_size() - 1)) {
                                // We're at the very end of the path, so this is a soft clip.
                                stats.total_softclipped_bases += edit.to_length() - edit.from_length();
                                stats.total_softclips++;
                                if(verbose) {
                                    // Record the actual insertion
                                    stats.softclips.push_back(make_pair(node_id, edit));
                                }
                            } else {
                                // This is not a softclip
                                has_non_softclip_indel_edits = true;
                                
                                // Record this insertion
                                stats.total_inserted_bases += edit.to_length() - edit.from_length();
                                stats.total_insertions++;
                                if(verbose) {
                                    // Record the actual insertion
                                    stats.insertions.push_back(make_pair(node_id, edit));
                                }
                            }

                        } else if(edit.from_length() > edit.to_length()) {
                            // This is a deletion and not a match
                            has_non_match_edits = true;
                            
                            // This is not a softclip either
                            has_non_softclip_indel_edits = true;
                            
                            // Record this deletion
                            stats.total_deleted_bases += edit.from_length() - edit.to_length();
                            stats.total_deletions++;
                            if(verbose) {
                                // Record the actual deletion
                                stats.deletions.push_back(make_pair(node_id, edit));
                            }
                        } else if(!edit.sequence().empty()) {
                            // This is a substitution and not a match
                            has_non_match_edits = true;
                        
                            // Record this substitution
                            // TODO: a substitution might also occur as part of a deletion/insertion above!
                            stats.total_substituted_bases += edit.from_length();
                            stats.total_substitutions++;
                            if(verbose) {
                                // Record the actual substitution
                                stats.substitutions.push_back(make_pair(node_id, edit));
                            }
                        }

                    }
                }

                for(auto& site_and_allele : alleles_supported) {
                    // This read is informative for an allele of a site.
                    // Up the reads on that allele of that site.
                    stats.reads_on_allele[site_and_allele.first][site_and_allele.second]++;
                }
            
                // If there's no non-match edits, call it a perfect alignment
                stats.total_perfect += !has_non_match_edits && has_alignment;
                
                // If there's no non-softclip indel edits, the alignment is gapless
                stats.total_gapless += !has_non_softclip_indel_edits && has_alignment;
            
            }

        };

        // Actually go through all the reads and count stuff up.
        vg::io::for_each_parallel(alignment_stream, lambda);
        
        // Now combine into a single ReadStats object (for which we pre-populated reads_on_allele with 0s).
        for (auto& per_thread : read_stats) {
            combined += per_thread;
        }
        read_stats.clear();

        // Go through all the nodes again and sum up unvisited nodes
        size_t unvisited_nodes = 0;
        // And unvisited base count
        size_t unvisited_node_bases = 0;
        // And nodes that are visited by only one thing (which is useful if
        // we're checking diploid assembly pairs).
        size_t single_visited_nodes = 0;
        size_t single_visited_node_bases = 0;
        // If we're in verbose mode, collect IDs too.
        set<vg::id_t> unvisited_ids;
        set<vg::id_t> single_visited_ids;
        // Note that you need to subtract out substituted-away and deleted bases
        // from the sum of 2 * double- and single-visited bases to get the bases
        // actually present in reads, because deleted bases are still "visited"
        // as many times as their nodes are touched. Also note that we ignore
        // edge effects and a read that stops before the end of a node will
        // visit the whole node.
        
        // These are for counting significantly allele-biased hets
        size_t total_hets = 0;
        size_t significantly_biased_hets = 0;

        if (graph.get() != nullptr) {

            // Calculate stats about the reads per allele data
            for(auto& site_and_alleles : combined.reads_on_allele) {
                // For every site
                if(site_and_alleles.second.size() == 2) {
                    // If it actually has 2 alleles with unique nodes in the
                    // graph (so we can use the binomial)

                    // We'll fill this with the counts for the two present alleles.
                    vector<size_t> counts;

                    for(auto& allele_and_count : site_and_alleles.second) {
                        // Collect all the counts
                        counts.push_back(allele_and_count.second);
                    }

                    if(counts[0] > counts[1]) {
                        // We have a 50% underlying probability so we can just put
                        // the rarer allele first.
                        swap(counts[0], counts[1]);
                    }

                    // What's the log prob for the smaller tail?
                    auto tail_logprob = binomial_cmf_ln(prob_to_logprob(0.5),  counts[1] + counts[0], counts[0]);

                    // Double it to get the two-tailed test
                    tail_logprob += prob_to_logprob(2);

#ifdef debug
                    cerr << "Site " << site_and_alleles.first << " has " << counts[0]
                        << " and " << counts[1] << " p=" << logprob_to_prob(tail_logprob) << endl;
#endif

                    if(tail_logprob < prob_to_logprob(0.05)) {
                        significantly_biased_hets++;
                    }
                    total_hets++;

                }
            }

            graph->for_each_handle([&](handle_t node) {
                // For every node
                
                // Look up its stats
                nid_t id = graph->get_id(node);
                size_t length = graph->get_length(node);
                
                if(!combined.node_visit_counts.count(id) || combined.node_visit_counts.at(id) == 0) {
                    // If we never visited it with a read, count it.
                    #pragma omp critical (unvisited_nodes)
                    unvisited_nodes++;
                    #pragma omp critical (unvisited_node_bases)
                    unvisited_node_bases += length;
                    if(verbose) {
                        #pragma omp critical (unvisited_ids)
                        unvisited_ids.insert(id);
                    }
                } else if(combined.node_visit_counts.at(id) == 1) {
                    // If we visited it with only one read, count it.
                    #pragma omp critical (single_visited_nodes)
                    single_visited_nodes++;
                    #pragma omp critical (single_visited_node_bases)
                    single_visited_node_bases += length;
                    if(verbose) {
                        #pragma omp critical (single_visited_ids)
                        single_visited_ids.insert(id);
                    }
                }
            });
            
        }

        cout << "Total alignments: " << combined.total_alignments << endl;
        cout << "Total primary: " << combined.total_primary << endl;
        cout << "Total secondary: " << combined.total_secondary << endl;
        cout << "Total aligned: " << combined.total_aligned << endl;
        cout << "Total perfect: " << combined.total_perfect << endl;
        cout << "Total gapless (softclips allowed): " << combined.total_gapless << endl;
        cout << "Total paired: " << combined.total_paired << endl;
        cout << "Total properly paired: " << combined.total_proper_paired << endl;

        cout << "Insertions: " << combined.total_inserted_bases << " bp in " << combined.total_insertions << " read events" << endl;
        if(verbose) {
            for(auto& id_and_edit : combined.insertions) {
                cout << "\t" << id_and_edit.second.from_length() << " -> " << id_and_edit.second.sequence()
                    << " on " << id_and_edit.first << endl;
            }
        }
        cout << "Deletions: " << combined.total_deleted_bases << " bp in " << combined.total_deletions << " read events" << endl;
        if(verbose) {
            for(auto& id_and_edit : combined.deletions) {
                cout << "\t" << id_and_edit.second.from_length() << " -> " << id_and_edit.second.to_length()
                    << " on " << id_and_edit.first << endl;
            }
        }
        cout << "Substitutions: " << combined.total_substituted_bases << " bp in " << combined.total_substitutions << " read events" << endl;
        if(verbose) {
            for(auto& id_and_edit : combined.substitutions) {
                cout << "\t" << id_and_edit.second.from_length() << " -> " << id_and_edit.second.sequence()
                    << " on " << id_and_edit.first << endl;
            }
        }
        cout << "Softclips: " << combined.total_softclipped_bases << " bp in " << combined.total_softclips << " read events" << endl;
        if(verbose) {
            for(auto& id_and_edit : combined.softclips) {
                cout << "\t" << id_and_edit.second.from_length() << " -> " << id_and_edit.second.sequence()
                    << " on " << id_and_edit.first << endl;
            }
        }
        
        if (graph.get() != nullptr) {
            cout << "Unvisited nodes: " << unvisited_nodes << "/" << graph->get_node_count()
                << " (" << unvisited_node_bases << " bp)" << endl;
            if(verbose) {
                for(auto& id : unvisited_ids) {
                    cout << "\t" << id << endl;
                }
            }

            cout << "Single-visited nodes: " << single_visited_nodes << "/" << graph->get_node_count()
                << " (" << single_visited_node_bases << " bp)" << endl;
            if(verbose) {
                for(auto& id : single_visited_ids) {
                    cout << "\t" << id << endl;
                }
            }

            cout << "Significantly biased heterozygous sites: " << significantly_biased_hets << "/" << total_hets;
            if(total_hets > 0) {
                cout << " (" << (double)significantly_biased_hets / total_hets * 100 << "%)";
            }
            cout << endl;
        }


    }
    
    if (snarl_stats) {
        // We will go through all the snarls and compute stats.
        
        require_graph();
        
        // First compute the snarls
        auto manager = IntegratedSnarlFinder(*graph).find_snarls_parallel();
        
        // We will track depth for each snarl
        unordered_map<const Snarl*, size_t> depth;

        // TSV header
        cout << "Start\tStart-Reversed\tEnd\tEnd-Reversed\tUltrabubble\tUnary\tShallow-Nodes\tShallow-Edges\tShallow-bases\tDeep-Nodes\tDeep-Edges\tDeep-Bases\tDepth\tChildren\tChains\tChains-Children\tNet-Graph-Size\n";
        
        manager.for_each_snarl_preorder([&](const Snarl* snarl) {
            // Loop over all the snarls and print stats.

            // snarl
            cout << snarl->start().node_id() << "\t" << snarl->start().backward() << "\t";
            cout << snarl->end().node_id() << "\t" << snarl->end().backward() << "\t";
            
            // Snarl metadata
            cout << (snarl->type() == ULTRABUBBLE) << "\t";
            cout << (snarl->type() == UNARY) << "\t";

            // Snarl size not including boundary nodes
            pair<unordered_set<vg::id_t>, unordered_set<vg::edge_t> > contents = manager.shallow_contents(snarl, *graph, false);
            size_t num_bases = 0;
            for (vg::id_t node_id : contents.first) {
                num_bases += graph->get_length(graph->get_handle(node_id));
            }
            cout << contents.first.size() << "\t";
            cout << contents.second.size() << "\t";
            cout << num_bases << "\t";
            contents = manager.deep_contents(snarl, *graph, false);
            num_bases = 0;
            for (vg::id_t node_id : contents.first) {
                num_bases += graph->get_length(graph->get_handle(node_id));
            }
            cout << contents.first.size() << "\t";
            cout << contents.second.size() << "\t";
            cout << num_bases << "\t";
            
            // Compute depth
            auto parent = manager.parent_of(snarl);
            
            if (parent == nullptr) {
                depth[snarl] = 0;
            } else {
                depth[snarl] = depth[parent] + 1;
            }
            cout << depth[snarl] << "\t";
            
            // Number of children (looking inside chains)
            cout << manager.children_of(snarl).size() << "\t";
            
            // Number of chains (including unary child snarls)
            // Will be 0 for leaves
            auto chains = manager.chains_of(snarl);
            cout << chains.size() << "\t";

            for (size_t i = 0; i < chains.size(); ++i) {
                // Number of children in each chain
                cout << chains[i].size();
                if (i < chains.size() - 1) {
                    cout << ",";
                }
            }
            if (chains.empty()) {
                cout << "0";
            }
            cout << "\t";
            
            // Net graph info
            // Internal connectivity not important, we just want the size.
            auto netGraph = manager.net_graph_of(snarl, graph.get(), false);
            cout << netGraph.get_node_count() << endl;
        });
        
    }

    if (!distance_index_filename.empty()) {
        //Print snarl stats from a distance index
        ifstream infile;
        infile.open(distance_index_filename);
        auto distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(infile);
        distance_index->print_snarl_stats();
        infile.close();
    }

    return 0;

}

// Register subcommand
static Subcommand vg_stats("stats", "metrics describing graph and alignment properties", TOOLKIT, main_stats);


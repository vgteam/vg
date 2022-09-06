/** \file depth_main.cpp
 *
 * Estimate sequencing depth from a (packed) alignment.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <algorithm>
#include <iostream>

#include "subcommand.hpp"

#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include "../handle.hpp"
#include <bdsg/overlays/overlay_helper.hpp>
#include "../utility.hpp"
#include "../packer.hpp"
#include "algorithms/coverage_depth.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_depth(char** argv) {
    cerr << "usage: " << argv[0] << " depth [options] <graph>" << endl
         << "options:" << endl
         << "  packed coverage depth (print 1-based positional depths along path):" << endl
         << "    -k, --pack FILE        supports created from vg pack for given input graph" << endl
         << "    -d, --count-dels       count deletion edges within the bin as covering reference positions" << endl
         << "  GAM/GAF coverage depth (print <mean> <stddev> for depth):" << endl
         << "    -g, --gam FILE         read alignments from this GAM file (could be '-' for stdin)" << endl
         << "    -a, --gaf FILE         read alignments from this GAF file (could be '-' for stdin)" << endl
         << "    -n, --max-nodes N      maximum nodes to consider [1000000]" << endl
         << "    -s, --random-seed N    random seed for sampling nodes to consider" << endl
         << "    -Q, --min-mapq N       ignore alignments with mapping quality < N [0]" << endl
         << "  path coverage depth (print 1-based positional depths along path):" << endl
         << "     activate by specifiying -p without -k" << endl
         << "    -c, --count-cycles     count each time a path steps on a position (by default paths are only counted once)" << endl
         << "  common options:" << endl
         << "    -p, --ref-path NAME    reference path to call on (multipile allowed.  defaults to all paths)" << endl
         << "    -P, --paths-by STR     select the paths with the given name prefix" << endl        
         << "    -b, --bin-size N       bin size (in bases) [1] (2 extra columns printed when N>1: bin-end-pos and stddev)" << endl
         << "    -m, --min-coverage N   ignore nodes with less than N coverage depth [1]" << endl
         << "    -t, --threads N        number of threads to use [all available]" << endl;
}

int main_depth(int argc, char** argv) {

    if (argc == 2) {
        help_depth(argv);
        return 1;
    }

    string pack_filename;
    unordered_set<string> ref_paths_input_set;
    vector<string> path_prefixes;
    size_t bin_size = 1;
    bool count_dels = false;
    
    string gam_filename;
    string gaf_filename;
    size_t max_nodes = 1000000;
    int random_seed = time(NULL);
    size_t min_mapq = 0;
    bool count_cycles = false;

    size_t min_coverage = 1;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {

        static const struct option long_options[] = {
            {"pack", required_argument, 0, 'k'},            
            {"ref-path", required_argument, 0, 'p'},
            {"paths-by", required_argument, 0, 'P'}, 
            {"bin-size", required_argument, 0, 'b'},
            {"count-dels", no_argument, 0, 'd'},
            {"gam", required_argument, 0, 'g'},
            {"gaf", no_argument, 0, 'a'},
            {"max-nodes", required_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"min-mapq", required_argument, 0, 'Q'},
            {"min-coverage", required_argument, 0, 'm'},
            {"count-cycles", no_argument, 0, 'c'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:p:P:b:dg:a:n:s:m:ct:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'k':
            pack_filename = optarg;
            break;
        case 'p':
            ref_paths_input_set.insert(optarg);
            break;
        case 'P':
            path_prefixes.push_back(optarg);
            break;
        case 'b':
            bin_size = parse<size_t>(optarg);
            break;
        case 'd':
            count_dels = true;
            break;            
        case 'g':
            gam_filename = optarg;
            break;
        case 'a':
            gaf_filename = optarg;
            break;
        case 'n':
            max_nodes = parse<size_t>(optarg);
            break;
        case 's':
            random_seed = parse<size_t>(optarg);
            break;
        case 'Q':
            min_mapq = parse<size_t>(optarg);
            break;
        case 'm':
            min_coverage = parse<size_t>(optarg);
            break;
        case 'c':
            count_cycles = true;
            break;
        case 't':
        {
            int num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg depth] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_depth(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_depth(argv);
        return 1;
    }

    size_t input_count = pack_filename.empty() ? 0 : 1;
    if (!gam_filename.empty()) ++input_count;
    if (!gaf_filename.empty()) ++input_count;
    if (input_count > 1) {                                          
        cerr << "error:[vg depth] At most one of a pack file (-k), a GAM file (-g), or a GAF file (-a) must be given" << endl;
        exit(1);
    }

    // Read the graph
    unique_ptr<PathHandleGraph> path_handle_graph;
    string path_handle_graph_filename = get_input_file_name(optind, argc, argv);
    path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(path_handle_graph_filename);
    PathHandleGraph* graph = path_handle_graph.get();
    
    // Apply the overlay if necessary
    bdsg::PathVectorizableOverlayHelper overlay_helper;
    if (!pack_filename.empty()) {
        graph = dynamic_cast<PathHandleGraph*>(overlay_helper.apply(path_handle_graph.get()));
        assert(graph != nullptr);
    }

    // Process the pack (or paths)
    unique_ptr<Packer> packer;
    if (!pack_filename.empty() || input_count == 0) {
        if (!pack_filename.empty()) {
            // Load our packed supports (they must have come from vg pack on graph)
            packer = unique_ptr<Packer>(new Packer(graph));
            packer->load_from_file(pack_filename);
        }

        // we want our paths sorted by the subpath parse so the output is sorted
        map<pair<string, int64_t>, string> ref_paths;
        unordered_set<string> base_path_set;
        
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                string path_name = graph->get_path_name(path_handle);
                tuple<bool, string, size_t, size_t> subpath_parse = Paths::parse_subpath_name(path_name);
                const string& base_name = get<0>(subpath_parse) ? get<1>(subpath_parse) : path_name;
                base_path_set.insert(base_name);
                // just take anything if no selection
                bool use_it = !Paths::is_alt(path_name) && path_prefixes.empty() && ref_paths_input_set.empty();

                // then look in the input paths -p
                if (!use_it && ref_paths_input_set.count(base_name)) {
                    use_it = true;
                }
                
                // then look in the prefixes
                for (size_t i = 0; i < path_prefixes.size() && !use_it; ++i) {
                    if (path_name.substr(0, path_prefixes[i].length()) == path_prefixes[i]) {
                        use_it = true;
                    }
                }
                if (use_it) {
                    auto coord = make_pair(base_name, get<2>(subpath_parse));
                    assert(!ref_paths.count(coord));
                    ref_paths[coord] = path_name;
                }
            });
        
        for (const auto& ref_name : ref_paths_input_set) {
            if (!base_path_set.count(ref_name)) {
                cerr << "error:[vg depth] Path \"" << ref_name << "\" not found in graph" << endl;
            }
        }

        for (const auto& ref_coord_path : ref_paths) {
            const string& ref_path = ref_coord_path.second;
            const string& base_path = ref_coord_path.first.first;
            const size_t subpath_offset = ref_coord_path.first.second;
            
            if (bin_size > 1) {
                vector<tuple<size_t, size_t, double, double>> binned_depth;
                if (!pack_filename.empty()) {
                    binned_depth = algorithms::binned_packed_depth(*packer, ref_path, bin_size, min_coverage, count_dels);
                } else {
                    binned_depth = algorithms::binned_path_depth(*graph, ref_path, bin_size, min_coverage, count_cycles);
                }
                for (auto& bin_cov : binned_depth) {
                    // bins can ben nan if min_coverage filters everything out.  just skip
                    if (!isnan(get<3>(bin_cov))) {
                        cout << base_path << "\t" << (get<0>(bin_cov) + 1 + subpath_offset)<< "\t" << (get<1>(bin_cov) + 1 + subpath_offset) << "\t" << get<2>(bin_cov)
                             << "\t" << sqrt(get<3>(bin_cov)) << endl;
                    }
                }
            } else {
                if (!pack_filename.empty()) {
                    algorithms::packed_depths(*packer, ref_path, min_coverage, cout);
                } else {
                    algorithms::path_depths(*graph, ref_path, min_coverage, count_cycles, cout);
                }
            }
        }
    }

    // Process the gam
    if (!gam_filename.empty() || !gaf_filename.empty()) {
        const string& mapping_filename = !gam_filename.empty() ? gam_filename : gaf_filename;
        pair<double, double> mapping_cov;
        mapping_cov = algorithms::sample_mapping_depth(*graph, mapping_filename, max_nodes, random_seed,
                                                       min_coverage, min_mapq, !gam_filename.empty() ? "GAM" : "GAF");
        cout << mapping_cov.first << "\t" << sqrt(mapping_cov.second) << endl;
    }
        
    return 0;

}

// Register subcommand
static Subcommand vg_depth("depth", "estimate sequencing depth", main_depth);


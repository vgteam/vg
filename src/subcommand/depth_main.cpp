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
         << "  packed coverage depth (print positional depths along path):" << endl
         << "    -k, --pack FILE        supports created from vg pack for given input graph" << endl
         << "    -d, --count-dels       count deletion edges within the bin as covering reference positions" << endl
         << "  GAM/GAF coverage depth (print <mean> <stddev> for depth):" << endl
         << "    -g, --gam FILE         read alignments from this GAM file (could be '-' for stdin)" << endl
         << "    -a, --gaf FILE         read alignments from this GAF file (could be '-' for stdin)" << endl
         << "    -n, --max-nodes N      maximum nodes to consider [1000000]" << endl
         << "    -s, --random-seed N    random seed for sampling nodes to consider" << endl
         << "    -Q, --min-mapq N       ignore alignments with mapping quality < N [0]" << endl
         << "  path coverage depth (print positional depths along path):" << endl
         << "     activate by specifiying -p without -k" << endl
         << "  common options:" << endl
         << "    -p, --ref-path NAME    reference path to call on (multipile allowed.  defaults to all paths)" << endl
         << "    -b, --bin-size N       bin size (in bases) [1] (2 extra columns printed when N>1: bin-end-pos and stddev)" << endl
         << "    -m, --min-coverage N   ignore nodes with less than N coverage [1]" << endl
         << "    -t, --threads N        number of threads to use [all available]" << endl;
}

int main_depth(int argc, char** argv) {

    if (argc == 2) {
        help_depth(argv);
        return 1;
    }

    string pack_filename;
    vector<string> ref_paths;
    size_t bin_size = 1;
    bool count_dels = false;
    
    string gam_filename;
    string gaf_filename;
    size_t max_nodes = 1000000;
    int random_seed = time(NULL);
    size_t min_mapq = 0;

    size_t min_coverage = 1;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {

        static const struct option long_options[] = {
            {"pack", required_argument, 0, 'k'},            
            {"ref-path", required_argument, 0, 'p'},
            {"bin-size", required_argument, 0, 'b'},
            {"count-dels", no_argument, 0, 'd'},
            {"gam", required_argument, 0, 'g'},
            {"gaf", no_argument, 0, 'a'},
            {"max-nodes", required_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"min-mapq", required_argument, 0, 'Q'},
            {"min-coverage", required_argument, 0, 'm'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:p:c:b:dg:a:n:s:m:t:",
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
            ref_paths.push_back(optarg);
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
    get_input_file(optind, argc, argv, [&](istream& in) {
            path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(in);
        });
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

        // All paths if none given
        if (ref_paths.empty()) {
            graph->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = graph->get_path_name(path_handle);
                    if (!Paths::is_alt(path_name)) {
                        ref_paths.push_back(path_name);
                    }
                });
        } else {
            for (const string& ref_name : ref_paths) {
                if (!graph->has_path(ref_name)) {
                    cerr << "error:[vg depth] Path \"" << ref_name << "\" not found in graph" << endl;
                }
            }
        }
        
        for (const string& ref_path : ref_paths) {
            if (bin_size > 1) {
                vector<tuple<size_t, size_t, double, double>> binned_depth;
                if (!pack_filename.empty()) {
                    binned_depth = algorithms::binned_packed_depth(*packer, ref_path, bin_size, min_coverage, count_dels);
                } else {
                    binned_depth = algorithms::binned_path_depth(*graph, ref_path, bin_size, min_coverage);
                }
                for (auto& bin_cov : binned_depth) {
                    // bins can ben nan if min_coverage filters everything out.  just skip
                    if (!isnan(get<3>(bin_cov))) {
                        cout << ref_path << "\t" << (get<0>(bin_cov) + 1)<< "\t" << (get<1>(bin_cov) + 1) << "\t" << get<2>(bin_cov)
                             << "\t" << sqrt(get<3>(bin_cov)) << endl;
                    }
                }
            } else {
                if (!pack_filename.empty()) {
                    algorithms::packed_depths(*packer, ref_path, min_coverage, cout);
                } else {
                    algorithms::path_depths(*graph, ref_path, min_coverage, cout);
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


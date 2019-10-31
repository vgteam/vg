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
#include <bdsg/overlay_helper.hpp>
#include "../utility.hpp"
#include "../packer.hpp"
#include "algorithms/coverage_depth.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_depth(char** argv) {
    cerr << "usage: " << argv[0] << " depth [options] <graph>" << endl
         << "options:" << endl
         << "  packed coverage depth:" << endl
         << "    -k, --pack FILE       Supports created from vg pack for given input graph" << endl
         << "    -p, --ref-path NAME   Reference path to call on (multipile allowed.  defaults to all paths)" << endl
         << "    -c, --context-size N  Context size (steps) for expanding bin subgraphs [50]" << endl
         << "    -b, --bin-size N      Bin size (in bases) [10000000]" << endl
         << "  GAM coverage depth:" << endl
         << "    -g, --gam FILE        read alignments from this file (could be '-' for stdin)" << endl
         << "    -n, --max-nodes N     maximum nodes to consider [1000000]" << endl
         << "    -s, --random-seed N   random seed for sampling nodes to consider" << endl
         << "  common options:" << endl
         << "    -m, --min-coverage N  ignore nodes with less than N coverage [1]" << endl
         << "    -t, --threads N       Number of threads to use [all available]" << endl;
}

int main_depth(int argc, char** argv) {

    if (argc == 2) {
        help_depth(argv);
        return 1;
    }

    string pack_filename;
    vector<string> ref_paths;
    size_t context_steps = 50;
    size_t bin_size = 10000000;
    
    string gam_filename;
    size_t max_nodes = 1000000;
    int random_seed = time(NULL);

    size_t min_coverage = 1;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {

        static const struct option long_options[] = {
            {"pack", required_argument, 0, 'k'},            
            {"ref-path", required_argument, 0, 'p'},
            {"context-size", required_argument, 0, 'c'},
            {"bin-size", required_argument, 0, 'b'},
            {"gam", required_argument, 0, 'g'},
            {"max-nodes", required_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"min-coverage", required_argument, 0, 'm'},
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:p:c:b:g:n:s:m:t:",
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
        case 'c':
            context_steps = parse<size_t>(optarg);
            break;
        case 'b':
            bin_size = parse<size_t>(optarg);
            break;
        case 'g':
            gam_filename = optarg;
            break;
        case 'n':
            max_nodes = parse<size_t>(optarg);
            break;
        case 's':
            random_seed = parse<size_t>(optarg);
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

    if (pack_filename.empty() == gam_filename.empty() ) {
        cerr << "error:[vg depth] Either a pack file (-k) or a gam file (-g) must be given" << endl;
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

    // Process the pack
    unique_ptr<Packer> packer;
    if (!pack_filename.empty()) {        
        // Load our packed supports (they must have come from vg pack on graph)
        packer = unique_ptr<Packer>(new Packer(graph));
        packer->load_from_file(pack_filename);

        // All paths if none given
        if (ref_paths.empty()) {
            graph->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = graph->get_path_name(path_handle);
                    if (!Paths::is_alt(path_name)) {
                        ref_paths.push_back(path_name);
                    }
                });
        }

        for (const string& ref_path : ref_paths) {
            map<size_t, double> binned_depth = algorithms::binned_packed_depth(*graph, *packer, ref_path, bin_size, get_thread_count());
            for (auto& bin_cov : binned_depth) {
                cerr << ref_path << "\t" << bin_cov.first << "\t" << bin_cov.second << endl;
            }
        }
    }

    // Process the gam
    if (!gam_filename.empty()) {
        double gam_cov;
        get_input_file(gam_filename, [&] (istream& gam_stream) {
                gam_cov = algorithms::sample_gam_depth(*graph, gam_stream, max_nodes, random_seed, min_coverage);
            });
        cerr << "gam-coverage\t" << gam_cov << endl;
    }

    return 0;

}

// Register subcommand
static Subcommand vg_depth("depth", "estimate sequencing depth", main_depth);


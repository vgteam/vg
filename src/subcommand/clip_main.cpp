#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../io/save_handle_graph.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <vg/io/alignment_emitter.hpp>
#include "../clip.hpp"
#include <bdsg/overlays/overlay_helper.hpp>

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;
using namespace vg::io;

void help_clip(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph>" << endl
       << "Chop out path intervals from a vg graph" << endl
       << endl
       << "bed clipping options: " << endl
       << "    -b, --bed FILE            Clip out alt-alleles in snarls that are contained in a region from given BED file" << endl
       << "depth clipping options: " << endl
       << "    -d, --depth N             Clip out alt-alleles with average depth below N" << endl
       << "    -P, --path-prefix STRING  Do not clip out alleles on paths beginning with given prefix." << endl
       << "general options: " << endl
       << "    -r, --snarls FILE         Snarls (from vg snarls) to avoid recomputing " << endl
       << "    -t, --threads N           number of threads to use (only used to computing snarls) [default: all available]" << endl
       << "    -v, --verbose             Print some logging messages" << endl
       << endl;
}    

int main_clip(int argc, char** argv) {

    string bed_path;
    string snarls_path;
    string ref_prefix;
    int64_t min_depth = 0;
    bool verbose = false;
    int input_count = 0;

    if (argc == 2) {
        help_clip(argv);
        return 1;
    }
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bed", required_argument, 0, 'b'},
            {"depth", required_argument, 0, 'd'},
            {"path-prefix", required_argument, 0, 'P'},
            {"snarls", required_argument, 0, 'r'},
            {"threads", required_argument, 0, 't'},
            {"verbose", required_argument, 0, 'v'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hb:d:P:r:t:v",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case '?':
        case 'h':
            help_clip(argv);
            return 0;
        case 'b':
            bed_path = optarg;
            ++input_count;
            break;
        case 'd':
            min_depth = stol(optarg);
            ++input_count;
            break;
        case 'P':
            ref_prefix = optarg;
            break;            
        case 'r':
            snarls_path = optarg;
            break;
        case 'v':
            verbose = true;
            break;
        case 't':
        {
            int num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg clip] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }            
        default:
            abort();
        }
    }

    if (input_count != 1) {
        cerr << "error:[vg clip] Exactly one of -b or -d must be specified to select what to clip" << endl;
        return 1;
    }

    // load the graph
    string graph_path = get_input_file_name(optind, argc, argv);
    unique_ptr<MutablePathMutableHandleGraph> graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(graph_path);

    // Load or compute the snarls
    unique_ptr<SnarlManager> snarl_manager;
    if (!snarls_path.empty()) {
        ifstream snarl_file(snarls_path.c_str());
        if (!snarl_file) {
            cerr << "Error [vg clip]: Unable to load snarls file: " << snarls_path << endl;
            return 1;
        }
        snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
        if (verbose) {
            cerr << "[vg clip]: Loaded " << snarl_manager->num_snarls() << " snarls" << endl;
        }
    } else {
        IntegratedSnarlFinder finder(*graph);
        snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls_parallel())));
        if (verbose) {
            cerr << "[vg clip]: Computed " << snarl_manager->num_snarls() << " snarls" << endl;
        }
    }

    if (!bed_path.empty()) {
        // load the bed file
        vector<Region> bed_regions;
        parse_bed_regions(bed_path, bed_regions);
        if (verbose) {
            cerr << "[vg clip]: Loaded " << bed_regions.size() << " BED regions" << endl;
        }

        // need path positions
        bdsg::PathPositionOverlayHelper overlay_helper;
        PathPositionHandleGraph* pp_graph = overlay_helper.apply(graph.get());

        // run the clipping
        clip_contained_snarls(graph.get(), pp_graph, bed_regions, *snarl_manager, false, verbose);
    }

    if (min_depth > 0) {
        // run the clipping
        //clip_low_depth_traversals(graph.get(), min_depth, *snarl_manager, hardmask);
    }
        
    // write the graph
    vg::io::save_handle_graph(graph.get(), std::cout);
    
    return 0;
}


// Register subcommand
static Subcommand vg_clip("clip", "remove BED regions (other other nodes from their snarls) from a graph", main_clip);

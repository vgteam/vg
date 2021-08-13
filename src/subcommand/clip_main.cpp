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
       << "Chop out variation within path intervals of a vg graph" << endl
       << endl
       << "input options: " << endl
       << "    -b, --bed FILE            BED regions corresponding to path intervals of the graph to target" << endl
       << "    -r, --snarls FILE         Snarls from vg snarls (recomputed if not given).  Snarls used to identify subgraphs within target intervals" << endl
       << "depth clipping options: " << endl
       << "    -d, --depth N             Clip out nodes with average depth below N" << endl
       << "    -P, --path-prefix STRING  Do not clip out alleles on paths beginning with given prefix (references must be specified either with -P or -b). " << endl
       << "general options: " << endl
       << "    -m, --min-fragment-len N  Don't write novel path fragment if it less than N bp long" << endl
       << "    -t, --threads N           number of threads to use (only used to computing snarls) [default: all available]" << endl
       << "    -v, --verbose             Print some logging messages" << endl
       << endl;
}    

int main_clip(int argc, char** argv) {

    string bed_path;
    string snarls_path;
    string ref_prefix;
    int64_t min_depth = -1;
    int64_t min_fragment_len = 0;
    bool verbose = false;
    bool depth_clipping = false;

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
            {"min-fragment-len", required_argument, 0, 'm'},
            {"threads", required_argument, 0, 't'},
            {"verbose", required_argument, 0, 'v'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hb:d:P:r:m:t:v",
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
            break;
        case 'd':
            min_depth = parse<size_t>(optarg);
            break;
        case 'P':
            ref_prefix = optarg;
            break;            
        case 'r':
            snarls_path = optarg;
            break;
        case 'm':
            min_fragment_len = parse<int>(optarg);
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

    if (min_depth >= 0) {
        if (bed_path.empty() == ref_prefix.empty()) {
            cerr << "error:[vg-clip] Depth clipping (-d) requires reference intervals specified with one of -b or -P" << endl;
            return 1;
        }
    } else if (bed_path.empty()) {
        cerr << "error:[vg-clip] BED intervals must be specified with -b" << endl;
        return 1;
    }

    // load the graph
    string graph_path = get_input_file_name(optind, argc, argv);
    unique_ptr<MutablePathMutableHandleGraph> graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(graph_path);

    // optional overlay only needed with bed regions input
    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* pp_graph = nullptr;

    unique_ptr<SnarlManager> snarl_manager;
    vector<Region> bed_regions;
    
    if (!bed_path.empty()) {
        // need path positions to intersect with regions
        pp_graph = overlay_helper.apply(graph.get());
        
        // Load or compute the snarls which are required for targetting bed regions
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
        
        // load the bed file
        parse_bed_regions(bed_path, bed_regions);
        if (verbose) {
            cerr << "[vg clip]: Loaded " << bed_regions.size() << " BED regions" << endl;
        }
    }

    if (min_depth >= 0) {
        // run the depth clipping       
        if (bed_regions.empty()) {            
            // do the whole graph
            clip_low_depth_nodes(graph.get(), min_depth, ref_prefix, min_fragment_len, verbose);
        } else {
            // do the contained snarls
            clip_contained_low_depth_nodes(graph.get(), pp_graph, bed_regions, *snarl_manager, false, min_depth, min_fragment_len, verbose);
        }
        
    } else {
        // run the alt-allele clipping
        clip_contained_snarls(graph.get(), pp_graph, bed_regions, *snarl_manager, false, min_fragment_len, verbose);
    }
        
    // write the graph
    vg::io::save_handle_graph(graph.get(), std::cout);
    
    return 0;
}


// Register subcommand
static Subcommand vg_clip("clip", "remove BED regions (other other nodes from their snarls) from a graph", main_clip);

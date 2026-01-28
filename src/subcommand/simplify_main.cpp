// simplify_main.cpp: define the "vg simplify" subcommand, which removes small variation

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include <vg/io/vpkg.hpp>

#include "../vg.hpp"
#include "../small_snarl_simplifier.hpp"
#include "../rare_variant_simplifier.hpp"
#include "../io/save_handle_graph.hpp"
#include "../traversal_clusters.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_simplify(char** argv) {
    cerr << "usage: " << argv[0] << " simplify [options] old.vg >new.vg" << endl
         << "general options:" << endl
         << "  -h, --help              print this help message to stderr and exit" << endl
         << "  -a, --algorithm NAME    simplification algorithm (small, rare) [small]" << endl
         << "  -t, --threads N         use N threads to construct graph [numCPUs]" << endl
         << "  -p, --progress          show progress" << endl
         << "  -b, --bed-in FILE       BED file with the cordinates of the original paths" << endl
         << "  -B, --bed-out FILE      output transformed features with new path coordinates" << endl
         << "path snarl simplifier options:" << endl
         << "  -P, --path-prefix STR   [NECESSARY TO SCALE PAST TINY GRAPHS]" << endl
         << "                          all paths with this prefix selected as reference paths" << endl
         << "                          (default: all reference-sense paths)" << endl
         << "small snarl simplifier options:" << endl       
         << "  -m, --min-size N        remove leaf sites with fewer than N bases" << endl
         << "                          (with -P, uses max allele length) involved [10]" << endl
         << "  -i, --max-iterations N  perform up to N iterations of simplification [10]" << endl
         << "  -L, --cluster F         cluster traversals whose (handle) Jaccard coefficient" << endl
         << "                          is >= F together [1.0]" << endl
         << "  -k, --keep-paths        non-reference (-P) paths are removed by default." << endl
         << "                          use this flag to keep them (the resulting graph will" << endl
         << "                          be more complex and possibly more difficult to load)" << endl
         << "rare variant simplifier options:" << endl
         << "  -v, --vcf FILE          use this VCF to determine variant frequency (required)" << endl
         << "  -f, --min-freq FLOAT    remove variants with total alt frequency <FLOAT [0]" << endl
         << "  -c, --min-count N       remove variants with total alt occurrence count <N [0]" << endl;
}

int main_simplify(int argc, char** argv) {
    Logger logger("vg simplify");

    if (argc == 2) {
        help_simplify(argv);
        return 1;
    }

    // What algorithm should we use for simplification ("path", "small" or "rare").
    string algorithm = "small";

    // General options
    string bed_in_filename;
    string bed_out_filename;
    bool show_progress = false;

    // for simplifying based on path traversals
    double cluster_threshold = 1.0;
    string ref_path_prefix;
    bool keep_nonref_paths = false;
    
    // For simplifying small variants
    size_t min_size = 10;
    size_t max_iterations = 10;

    // For simplifying rare variants
    string vcf_filename;
    double min_frequency = 0;
    size_t min_count = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"algorithm", required_argument, 0, 'a'},
                {"progress",  no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {"bed-in", required_argument, 0, 'b'},
                {"bed-out", required_argument, 0, 'B'},
                {"min-size", required_argument, 0, 'm'},
                {"max-iterations", required_argument, 0, 'i'},
                {"vcf", required_argument, 0, 'v'},
                {"min-freq", required_argument, 0, 'f'},
                {"min-count", required_argument, 0, 'c'},
                {"cluster", required_argument, 0, 'L'},
                {"keep-paths", no_argument, 0, 'k'},
                {"path-prefix", required_argument, 0, 'P'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "a:pt:b:B:m:i:v:f:c:L:kP:h?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'a':
            algorithm = optarg;
            break;

        case 'p':
            show_progress = true;
            break;

        case 't':
            set_thread_count(logger, optarg);
            break;
            
        case 'b':
            bed_in_filename = require_exists(logger, optarg);
            break;
        
        case 'B':
            bed_out_filename = ensure_writable(logger, optarg);
            break;

        case 'm':
            min_size = parse<int>(optarg);
            break;

        case 'L':
            cluster_threshold = parse<double>(optarg);
            break;

        case 'k':
            keep_nonref_paths = true;
            break;            

        case 'P':
            ref_path_prefix = optarg;
            break;
            
        case 'i':
            max_iterations = parse<int>(optarg);
            break;

        case 'v':
            vcf_filename = require_exists(logger, optarg);
            break;

        case 'f':
            min_frequency = parse<double>(optarg);
            break;

        case 'c':
            min_count = parse<size_t>(optarg);
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_simplify(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }

    // Do preliminary options checks
    if (!bed_out_filename.empty() && bed_in_filename.empty()) {
        // Don't allow writing out a BED without reading one
        logger.error() << "Cannot output a BED (-B) unless a BED is read in first (-b)" << endl;
    }
    
    if (algorithm != "small" && !ref_path_prefix.empty()) {
        logger.error() << "Path simplification (-P) can only be used with -a small" << endl;
    }

    if (algorithm == "small") {
        if (ref_path_prefix.empty()) {
        logger.warn() << "By not specifying a reference path (-P) you are using old logic which requires "
                      << "protobuf input, and scales very poorly" << endl;
        }
        if (!vcf_filename.empty()) {
            logger.error() << "A VCF file (-v) cannot be used with small snarl simplification" << endl;
        }
    }

    if (algorithm == "rare") {
        if (vcf_filename.empty()) {
            logger.error() << "The \"rare\" simplification algorithm requires a VCF (-v)" << endl;
        }
    }
    
    // Load the graph
    unique_ptr<handlegraph::MutablePathDeletableHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        if (!ref_path_prefix.empty()) {
            graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
        } else {
            try {
                graph = unique_ptr<MutablePathDeletableHandleGraph>(new VG(in, show_progress));
            } catch(...) {
                logger.error() << "Error loading input Protobuf graph." << endl;
            }
        }
    });

    if (graph == nullptr) {
        logger.error() << "Could not load graph." << endl;
    }

    // This will hold BED features if we are tracking those
    unique_ptr<FeatureSet> features;
    if (!bed_in_filename.empty()) {
        // Go and load up the BED features
        get_input_file(bed_in_filename, [&](istream& bed_stream) {
            features = unique_ptr<FeatureSet>(new FeatureSet());
            features->load_bed(bed_stream);
        });
    }

    if (!ref_path_prefix.empty()) {
        if (!vcf_filename.empty() || !bed_in_filename.empty() || !bed_out_filename.empty()) {
            logger.error()<< "-v/-b/-B options cannot be used with path-based simplification (-P)" << endl;
        }

        nid_t min_id = graph->min_node_id();
        
        simplify_graph_using_traversals(dynamic_cast<MutablePathMutableHandleGraph*>(graph.get()),
                                        ref_path_prefix, min_size, cluster_threshold, max_iterations, 100000);

        if (!keep_nonref_paths) {
            vector<path_handle_t> to_destroy;
            graph->for_each_path_of_sense({PathSense::REFERENCE, PathSense::GENERIC}, [&](const path_handle_t path_handle) {
                if (graph->get_path_name(path_handle).compare(0, ref_path_prefix.length(), ref_path_prefix) != 0) {
                    to_destroy.push_back(path_handle);
                }
            });
            graph->destroy_paths(to_destroy);            
        }
        
        handlealgs::unchop(*graph);

        if (graph->min_node_id() < min_id) {
            // we want chromosome graphs to keep disjoint node ids, so we preserve the min_id            
            graph->increment_node_ids(min_id - graph->min_node_id());
        }
        
    } else if (algorithm == "small") {
        // Make a SmallSnarlSimplifier for the graph and copy over settings.
        SmallSnarlSimplifier simplifier(*dynamic_cast<VG*>(graph.get()));
        simplifier.show_progress = show_progress;
        simplifier.max_iterations = max_iterations;
        simplifier.min_size = min_size;
        simplifier.features = features.get();
         
        // Do the simplification
        simplifier.simplify();
    } else if (algorithm == "rare") {
        // We are going to remove rare variants as noted in a VCF

        // Load the VCF
        vcflib::VariantCallFile variant_file;
        variant_file.parseSamples = false; // Major speedup if there are many samples.
        variant_file.open(vcf_filename);
        if (!variant_file.is_open()) {
            logger.error() << "could not open " << vcf_filename << endl;
        }

        // Buffer it
        VcfBuffer buffer(&variant_file);

        // Make a RareVariantSimplifier for the graph
        RareVariantSimplifier simplifier(*graph, buffer);

        // Set its settings
        simplifier.min_frequency_to_keep = min_frequency;
        simplifier.min_count_to_keep = min_count;

        // Run it
        simplifier.simplify();
    } else {
        logger.error() << "Unknown algorithm \"" 
                       << algorithm << "\"; use \"small\" or \"rare\"." << endl;
    }

    // Serialize the graph
    if (!ref_path_prefix.empty()) {
        vg::io::save_handle_graph(graph.get(), std::cout);
    } else {
        dynamic_cast<VG*>(graph.get())->serialize_to_ostream(std::cout);
    }
        
    if (!bed_out_filename.empty()) {
        // Save BED features
        assert(features.get() != nullptr);
        ofstream bed_stream(bed_out_filename.c_str());
        features->save_bed(bed_stream);
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_simplify("simplify", "graph simplification", main_simplify);


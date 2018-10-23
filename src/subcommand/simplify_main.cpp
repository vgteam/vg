// simplify_main.cpp: define the "vg simplify" subcommand, which removes small variation

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../small_snarl_simplifier.hpp"
#include "../rare_variant_simplifier.hpp"



using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_simplify(char** argv) {
    cerr << "usage: " << argv[0] << " simplify [options] old.vg >new.vg" << endl
         << "general options:" << endl
         << "    -a, --algorithm NAME   simplify using the given algorithm (small, rare; default: small)" << endl
         << "    -t, --threads N        use N threads to construct graph (defaults to numCPUs)" << endl
         << "    -p, --progress         show progress" << endl
         << "    -b, --bed-in           read in the given BED file in the cordinates of the original paths" << endl
         << "    -B, --bed-out          output transformed features in the coordinates of the new paths" << endl
         << "small snarl simplifier options:" << endl
         << "    -m, --min-size N       remove leaf sites with fewer than N bases involved (default: 10)" << endl
         << "    -i, --max-iterations N perform up to N iterations of simplification (default: 10)" << endl
         << "rare variant simplifier options:" << endl
         << "    -v, --vcf FILE         use the given VCF file to determine variant frequency (required)" << endl
         << "    -f, --min-freq FLOAT   remove variants with total alt frequency under FLOAT (default: 0)" << endl
         << "    -c, --min-count N      remove variants with total alt occurrence count under N (default: 0)" << endl;
}

int main_simplify(int argc, char** argv) {

    if (argc == 2) {
        help_simplify(argv);
        return 1;
    }

    // What algorithm should we use for simplification ("small" or "rare").
    string algorithm = "small";

    // General options
    string bed_in_filename;
    string bed_out_filename;
    bool show_progress = false;

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
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "a:pt:b:B:m:i:v:f:c:h?",
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
            omp_set_num_threads(parse<int>(optarg));
            break;
            
        case 'b':
            bed_in_filename = optarg;
            break;
        
        case 'B':
            bed_out_filename = optarg;
            break;

        case 'm':
            min_size = parse<int>(optarg);
            break;
            
        case 'i':
            max_iterations = parse<int>(optarg);
            break;

        case 'v':
            vcf_filename = optarg;
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
        cerr << "error[vg simplify]: Cannot output a bed (-B) unless a BED is read in first (-b)" << endl;
        exit(1);
    }
    
    // Load the graph
    unique_ptr<VG> graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = unique_ptr<VG>(new VG(in, show_progress));
    });
    
    if (graph == nullptr) {
        cerr << "error:[vg simplify]: Could not load graph" << endl;
        exit(1);
    }

    // This will hold BED features if we are tracking those
    unique_ptr<FeatureSet> features;
    if (!bed_in_filename.empty()) {
        // Go and ,oad up the BED features
        get_input_file(bed_in_filename, [&](istream& bed_stream) {
            features = unique_ptr<FeatureSet>(new FeatureSet());
            features->load_bed(bed_stream);
        });
    }

    if (algorithm == "small") {
        if (!vcf_filename.empty()) {
            cerr << "error[vg simplify]: A VCF file (-v) cannot be used with small snarl simplification" << endl;
            exit(1);
        }

        // Make a SmallSnarlSimplifier for the graph and copy over settings.
        SmallSnarlSimplifier simplifier(*graph);
        simplifier.show_progress = show_progress;
        simplifier.max_iterations = max_iterations;
        simplifier.min_size = min_size;
        simplifier.features = features.get();
         
        // Do the simplification
        simplifier.simplify();
    } else if (algorithm == "rare") {
        // We are going to remove rare variants as noted in a VCF
        if (vcf_filename.empty()) {
            cerr << "error[vg simplify]: \"rare\" simplification algorithm requires a VCF (-v)" << endl;
            exit(1);
        }

        // Load the VCF
        vcflib::VariantCallFile variant_file;
        variant_file.parseSamples = false; // Major speedup if there are many samples.
        variant_file.open(vcf_filename);
        if (!variant_file.is_open()) {
            cerr << "error:[vg simplify] could not open" << vcf_filename << endl;
            exit(1);
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
        cerr << "error[vg simplify]: Unknown algorithm \"" << algorithm << "\"; use \"small\" or \"rare\"." << endl;
        exit(1);
    }

    // Serialize the graph
    graph->serialize_to_ostream(std::cout);
        
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


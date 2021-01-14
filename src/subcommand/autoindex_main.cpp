/** \file autoindex_main.cpp
 *
 * Defines the "vg autoindex" subcommand, which produces indexes needed for other subcommands
 */
#include <getopt.h>
#include <iostream>

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include "subcommand.hpp"
#include "index_registry.hpp"
#include "utility.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

bool vcf_is_phased(const string& filepath) {
    
    if (IndexingParameters::verbose) {
        cerr << "[IndexRegistry]: Checking for phasing in VCF." << endl;
    }
    
    // check about 30k variants before concluding that the VCF isn't phased
    // TODO: will there be contig ordering biases that make this a bad assumption?
    constexpr int vars_to_check = 1 << 15;
    
    htsFile* file = hts_open(filepath.c_str(),"rb");
    bcf_hdr_t* header = bcf_hdr_read(file);
    int phase_set_id = bcf_hdr_id2int(header, BCF_DT_ID, "PS");
    if (phase_set_id < 0) {
        // no PS tag means no phasing
        bcf_hdr_destroy(header);
        hts_close(file);
        return false;
    }
    
    // iterate over records
    bcf1_t* bcf_rec = bcf_init();
    int iter = 0;
    bool found_phased = false;
    while (bcf_read(file, header, bcf_rec) >= 0 && iter < vars_to_check && !found_phased)
    {
        bcf_unpack(bcf_rec, BCF_UN_ALL);
        if (phase_set_id == BCF_HT_INT) {
            // phase sets are integers
            int num_phase_set_arr = 0;
            int32_t* phase_sets = NULL;
            int num_phase_sets = bcf_get_format_int32(header, bcf_rec, "PS", &phase_sets, &num_phase_set_arr);
            for (int i = 0; i < num_phase_sets && !found_phased; ++i) {
                found_phased = phase_sets[i] != 0;
            }
            free(phase_sets);
        }
        else if (phase_set_id == BCF_HT_STR) {
            // phase sets are strings
            int num_phase_set_arr = 0;
            char** phase_sets = NULL;
            int num_phase_sets = bcf_get_format_string(header, bcf_rec, "PS", &phase_sets, &num_phase_set_arr);
            for (int i = 0; i < num_phase_sets && !found_phased; ++i) {
                found_phased = strcmp(phase_sets[i], ".") != 0;
            }
            if (phase_sets) {
                // all phase sets are concatenated in one malloc's char*, pointed to by the first pointer
                free(phase_sets[0]);
            }
            // free the array of pointers
            free(phase_sets);
        }
        else {
            cerr << "warning: unrecognized PS type in VCF file " << filepath << ", interpreting as unphased." << endl;
            break;
        }
        ++iter;
    }
    // clean up
    bcf_destroy(bcf_rec);
    bcf_hdr_destroy(header);
    hts_close(file);
    
    return found_phased;
}

void help_autoindex(char** argv) {
    cerr
    << "usage: " << argv[0] << " autoindex [options]" << endl
    << "options:" << endl
    << "  output:" << endl
    << "    -p, --prefix PREFIX   prefix to use for all output (default: index)" << endl
    << "    -w, --workflow NAME   workflow to produce indexes for: map, mpmap, giraffe (default: map)" << endl
    << "  input data:" << endl
    << "    -r, --ref-fasta FILE  FASTA file containing the reference sequence" << endl
    << "    -v, --vcf FILE        VCF file with sequence names matching -r" << endl
    << "    -i, --ins-fasta FILE  FASTA file with sequences of INS variants from -v" << endl
    << "    -g, --gfa FILE        GFA file to make a graph from" << endl
    << "  logging and computation:" << endl
    << "    -T, --tmp-dir DIR     temporary directory to use for intermediate files" << endl
    << "    -V, --verbose         log progress to stderr" << endl
    << "    -d, --dot             print the dot-formatted graph of index recipes and exit" << endl
    << "    -h, --help            print this help message to stderr and exit" << endl;
}

int main_autoindex(int argc, char** argv) {

    if (argc == 2) {
        help_autoindex(argv);
        return 1;
    }
    
    // load the registry
    IndexRegistry registry = VGIndexes::get_vg_index_registry();
    bool print_dot = false;
    vector<string> targets;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"prefix", required_argument, 0, 'p'},
            {"workflow", required_argument, 0, 'w'},
            {"ref-fasta", required_argument, 0, 'r'},
            {"vcf", required_argument, 0, 'v'},
            {"ins-fasta", required_argument, 0, 'i'},
            {"gfa", required_argument, 0, 'g'},
            {"tmp-dir", required_argument, 0, 'T'},
            {"verbose", no_argument, 0, 'V'},
            {"dot", no_argument, 0, 'd'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "p:w:r:v:i:g:T:dVh",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'p':
                registry.set_prefix(optarg);
                break;
            case 'w':
                if (optarg == string("map")) {
                    targets = VGIndexes::get_default_map_indexes();
                }
                else if (optarg == string("mpmap")) {
                    targets = VGIndexes::get_default_mpmap_indexes();
                    cerr << "mpmap indexing not yet implemented" << endl;
                    return 1;
                }
                else if (optarg == string("giraffe")) {
                    targets = VGIndexes::get_default_giraffe_indexes();
                    cerr << "giraffe indexing not yet implemented" << endl;
                    return 1;
                }
                else {
                    cerr << "error: Unrecognized workflow (-w): " << optarg << endl;
                    return 1;
                }
                break;
            case 'r':
                registry.provide("Reference FASTA", optarg);
                break;
            case 'v':
                if (vcf_is_phased(optarg)) {
                    registry.provide("Phased VCF", optarg);
                }
                else {
                    registry.provide("VCF", optarg);
                }
                break;
            case 'i':
                registry.provide("Insertion Sequence FASTA", optarg);
                break;
            case 'g':
                registry.provide("Reference GFA", optarg);
                break;
            case 'T':
                temp_file::set_dir(optarg);
                break;
            case 'V':
                IndexingParameters::verbose = true;
                break;
            case 'd':
                print_dot = true;
                break;
            case 'h':
                help_autoindex(argv);
                return 0;
            default:
                abort ();
        }
    }

    if (print_dot) {
        // don't index, just visualize the plan
        cout << registry.to_dot(targets);
        return 0;
    }
    
    if (targets.empty()) {
        // default to vg map indexes
        targets = VGIndexes::get_default_map_indexes();
    }
    
    registry.make_indexes(targets);
    
    return 0;

}

// Register subcommand
static Subcommand vg_autoindex("autoindex", "produce indexes for other subcommands", TOOLKIT, main_autoindex);


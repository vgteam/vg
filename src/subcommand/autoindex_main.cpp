/** \file autoindex_main.cpp
 *
 * Defines the "vg autoindex" subcommand, which produces indexes needed for other subcommands
 */
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <omp.h>

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
    
    
    htsFile* file = hts_open(filepath.c_str(), "rb");
    bcf_hdr_t* hdr = bcf_hdr_read(file);
    int phase_set_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "PS");
    // note: it seems that this is not necessary for expressing phasing after all
//    if (phase_set_id < 0) {
//        // no PS tag means no phasing
//        bcf_hdr_destroy(hdr);
//        hts_close(file);
//        return false;
//    }
    
    // iterate over records
    bcf1_t* line = bcf_init();
    int iter = 0;
    bool found_phased = false;
    while (bcf_read(file, hdr, line) >= 0 && iter < vars_to_check && !found_phased)
    {
        if (phase_set_id >= 0) {
            if (phase_set_id == BCF_HT_INT) {
                // phase sets are integers
                int num_phase_set_arr = 0;
                int32_t* phase_sets = NULL;
                int num_phase_sets = bcf_get_format_int32(hdr, line, "PS", &phase_sets, &num_phase_set_arr);
                for (int i = 0; i < num_phase_sets && !found_phased; ++i) {
                    found_phased = phase_sets[i] != 0;
                }
                free(phase_sets);
            }
            else if (phase_set_id == BCF_HT_STR) {
                // phase sets are strings
                int num_phase_set_arr = 0;
                char** phase_sets = NULL;
                int num_phase_sets = bcf_get_format_string(hdr, line, "PS", &phase_sets, &num_phase_set_arr);
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
        }
        
        // init a genotype array
        int32_t* genotypes = nullptr;
        int arr_size = 0;
        // and query it
        int num_genotypes = bcf_get_genotypes(hdr, line, &genotypes, &arr_size);
        if (num_genotypes >= 0) {
            // we got genotypes, check to see if they're phased
            int num_samples = bcf_hdr_nsamples(hdr);
            int ploidy = num_genotypes / num_samples;
            for (int i = 0; i < num_genotypes && !found_phased; i += ploidy) {
                for (int j = 0; j < ploidy && !found_phased; ++j) {
                    if (genotypes[i + j] == bcf_int32_vector_end) {
                        // sample has lower ploidy
                        break;
                    }
                    if (bcf_gt_is_missing(genotypes[i + j])) {
                        continue;
                    }
                    if (bcf_gt_is_phased(genotypes[i + j])) {
                        // the VCF expresses phasing, we can
                        found_phased = true;;
                    }
                }
            }
        }
        
        free(genotypes);
        ++iter;
    }
    // clean up
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    hts_close(file);
    return found_phased;
}

void help_autoindex(char** argv) {
    cerr
    << "usage: " << argv[0] << " autoindex [options]" << endl
    << "options:" << endl
    << "  output:" << endl
    << "    -p, --prefix PREFIX    prefix to use for all output (default: index)" << endl
    << "    -w, --workflow NAME    workflow to produce indexes for, can be provided multiple" << endl
    << "                           times. options: map, mpmap, giraffe (default: map)" << endl
    << "  input data:" << endl
    << "    -r, --ref-fasta FILE   FASTA file containing the reference sequence (may repeat)" << endl
    << "    -v, --vcf FILE         VCF file with sequence names matching -r (may repeat)" << endl
    << "    -i, --ins-fasta FILE   FASTA file with sequences of INS variants from -v" << endl
    << "    -g, --gfa FILE         GFA file to make a graph from" << endl
    << "    -x, --tx-gff FILE      GTF/GFF file with transcript annotations (may repeat)" << endl
    << "  configuration:" << endl
    << "    -f, --gff-feature STR  GTF/GFF feature type (col. 3) to add to graph (default: " << IndexingParameters::gff_feature_name << ")" << endl
    << "    -a, --gff-tx-tag STR   GTF/GFF tag (in col. 9) for transcript ID (default: " << IndexingParameters::gff_transcript_tag << ")" << endl
    << "  logging and computation:" << endl
    << "    -T, --tmp-dir DIR      temporary directory to use for intermediate files" << endl
    << "    -t, --threads NUM      number of threads (default: all available)" << endl
    << "    -V, --verbose          log progress to stderr" << endl
    << "    -d, --dot              print the dot-formatted graph of index recipes and exit" << endl
    << "    -h, --help             print this help message to stderr and exit" << endl;
}

int main_autoindex(int argc, char** argv) {

    if (argc == 2) {
        help_autoindex(argv);
        return 1;
    }
    
#define OPT_KEEP_INTERMEDIATE 1000
#define OPT_FORCE_UNPHASED 1001
#define OPT_FORCE_PHASED 1002
    
    // load the registry
    IndexRegistry registry = VGIndexes::get_vg_index_registry();
    bool print_dot = false;
    vector<IndexName> targets;
    vector<string> vcf_names;
    bool force_unphased = false;
    bool force_phased = false;
    
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
            {"tx-gff", required_argument, 0, 'x'},
            {"gff-feature", required_argument, 0, 'f'},
            {"gff-tx-tag", required_argument, 0, 'a'},
            {"tmp-dir", required_argument, 0, 'T'},
            {"threads", required_argument, 0, 't'},
            {"verbose", no_argument, 0, 'V'},
            {"dot", no_argument, 0, 'd'},
            {"help", no_argument, 0, 'h'},
            {"keep-intermediate", no_argument, 0, OPT_KEEP_INTERMEDIATE},
            {"force-unphased", no_argument, 0, OPT_FORCE_UNPHASED},
            {"force-phased", no_argument, 0, OPT_FORCE_PHASED},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "p:w:r:v:i:g:x:a:f:T:t:dVh",
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
                    for (auto& target : VGIndexes::get_default_map_indexes()) {
                        targets.emplace_back(move(target));
                    }
                }
                else if (optarg == string("mpmap")) {
                    for (auto& target : VGIndexes::get_default_mpmap_indexes()) {
                        targets.emplace_back(move(target));
                    }
                }
                else if (optarg == string("giraffe")) {
                    for (auto& target : VGIndexes::get_default_giraffe_indexes()) {
                        targets.emplace_back(move(target));
                    }                    cerr << "giraffe indexing not yet implemented" << endl;
                    return 1;
                }
                else {
                    cerr << "error: Unrecognized workflow (-w): " << optarg << endl;
                    return 1;
                }
                break;
            case 'r':
                registry.provide({"Reference FASTA"}, optarg);
                break;
            case 'v':
                vcf_names.push_back(optarg);
                break;
            case 'i':
                registry.provide({"Insertion Sequence FASTA"}, optarg);
                break;
            case 'g':
                registry.provide({"Reference GFA"}, optarg);
                break;
            case 'x':
                registry.provide({"GTF/GFF"}, optarg);
                break;
            case 'f':
                IndexingParameters::gff_feature_name = optarg;
                break;
            case 'a':
                IndexingParameters::gff_transcript_tag = optarg;
                break;
            case 'T':
                temp_file::set_dir(optarg);
                break;
            case 't':
                omp_set_num_threads(parse<int>(optarg));
                break;
            case 'V':
                IndexingParameters::verbose = true;
                break;
            case 'd':
                print_dot = true;
                break;
            case OPT_KEEP_INTERMEDIATE:
                registry.set_intermediate_file_keeping(true);
                break;
            case OPT_FORCE_UNPHASED:
                force_unphased = true;
                break;
            case OPT_FORCE_PHASED:
                force_phased = true;
                break;
            case 'h':
                help_autoindex(argv);
                return 0;
            default:
                abort ();
        }
    }
    
    assert(!(force_phased && force_unphased));
    
    // we have special logic for VCFs to make it friendly to both phased
    // and unphased VCF files
    if (!vcf_names.empty()) {
        // we interpret it as a phased VCF if any of the VCFs have phasing
        bool phased = force_phased;
        if (!force_unphased) {
            for (size_t i = 0; i < vcf_names.size() && !phased; ++i) {
                phased = vcf_is_phased(vcf_names[i]);
            }
        }
        
        for (auto& vcf_name : vcf_names) {
            if (phased) {
                registry.provide({"VCF w/ Phasing"}, vcf_name);
            }
            else {
                registry.provide({"VCF"}, vcf_name);
            }
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
    // deduplicate
    sort(targets.begin(), targets.end());
    targets.resize(unique(targets.begin(), targets.end()) - targets.begin());
    
    try {
        registry.make_indexes(targets);
    }
    catch (InsufficientInputException ex) {
        cerr << "error:[vg autoindex] Input is not sufficient to create indexes" << endl;
        cerr << ex.what();
        return 1;
    }
    
    return 0;

}

// Register subcommand
static Subcommand vg_autoindex("autoindex", "produce indexes for other subcommands", TOOLKIT, main_autoindex);


/** \file autoindex_main.cpp
 *
 * Defines the "vg autoindex" subcommand, which produces indexes needed for other subcommands
 */
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <unistd.h>
#include <getopt.h>

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <gbwt/utils.h>
#include <omp.h>

#include "subcommand.hpp"
#include "index_registry.hpp"
#include "utility.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

int64_t parse_memory_usage(const string& mem_arg) {
    if (mem_arg.empty()) {
        cerr << "error:[vg autoindex] target memory usage arg is empty" << endl;
        exit(1);
    }
    string mem = mem_arg;
    if (mem.back() == 'B') {
        mem.pop_back();
    }
    int64_t base;
    if (mem.back() == 'k') {
        base = 1024;
        mem.pop_back();
    }
    else if (mem.back() == 'M') {
        base = 1024 * 1024;
        mem.pop_back();
    }
    else if (mem.back() == 'G') {
        base = 1024 * 1024 * 1024;
        mem.pop_back();
    }
    else if (isdigit(mem.back())) {
        base = 1;
    }
    else {
        cerr << "error:[vg autoindex] unrecognized unit for target memory usage: " << mem.back() << endl;
        exit(1);
    }
    return parse<int64_t>(mem) * base;
}

string mem_usage_string(int64_t mem) {
    assert(mem > 0);
    stringstream strm;
    strm.precision(1);
    if (mem >= 1024 * 1024 * 1024) {
        strm << double(mem) / (1024 * 1024 * 1024) << "GB";
    }
    else if (mem >= 1024 * 1024) {
        strm << double(mem) / (1024 * 1024) << "MB";
    }
    else if (mem >= 1024) {
        strm << double(mem) / (1024) << "kB";
    }
    else {
        strm << double(mem) << "B";
    }
    return strm.str();
};

// expects a string of form "Index Registry Name:filepath1,filepath2,filepath3"
pair<string, vector<string>> parse_provide_string(const string& str) {
    
    pair<string, vector<string>> return_val;
    
    size_t i = str.find(':');
    if (i >= str.size()) {
        cerr << "error: Couldn't parse index provide string: " << str << endl;
        exit(1);
    }
    return_val.first = str.substr(0, i);
    while (i < str.size()) {
        size_t end = str.find(',', i + 1);
        return_val.second.emplace_back(str.substr(i + 1, end - i - 1));
        i = end;
    }
    if (return_val.second.empty()) {
        cerr << "error: Couldn't parse index provide string: " << str << endl;
        exit(1);
    }
    return return_val;
}

void help_autoindex(char** argv) {
    cerr << "usage: " << argv[0] << " autoindex [options]" << endl
         << "output:" << endl
         << "  -p, --prefix PREFIX    prefix to use for all output [index]" << endl
         << "  -w, --workflow NAME    workflow to produce indexes for (may repeat) [map]" << endl
         << "                         {map, mpmap, rpvg, giraffe, sr-giraffe, lr-giraffe}" << endl
         << "input data:" << endl
         << "  -r, --ref-fasta FILE   FASTA file with the reference sequence (may repeat)" << endl
         << "  -v, --vcf FILE         VCF file with sequence names matching -r (may repeat)" << endl
         << "  -i, --ins-fasta FILE   FASTA file with sequences of INS variants from -v" << endl
         << "  -g, --gfa FILE         GFA file to make a graph from" << endl
         << "  -G, --gbz FILE         GBZ file to make indexes from" << endl
         << "  -x, --tx-gff FILE      GTF/GFF file with transcript annotations (may repeat)" << endl
         << "  -H, --hap-tx-gff FILE  GTF/GFF file with transcript annotations " << endl
         << "                         of a named haplotype (may repeat)" << endl
         << "  -n, --no-guessing      do not guess that pre-existing files are indexes" << endl
         << "                         i.e. force-regenerate any index not explicitly provided" << endl
         << "configuration:" << endl
         << "  -f, --gff-feature STR  GTF/GFF feature type (col. 3) to add to graph "
                                      << "[" << IndexingParameters::gff_feature_name << "]" << endl
         << "  -a, --gff-tx-tag STR   GTF/GFF tag (in col. 9) for ID "
                                      << "[" << IndexingParameters::gff_transcript_tag << "]" << endl
         << "logging and computation:" << endl
         << "  -T, --tmp-dir DIR      temporary directory to use for intermediate files" << endl
         << "  -M, --target-mem MEM   target max memory usage (not exact, formatted INT[kMG])" << endl
         << "                         [1/2 of available]" << endl
// TODO: hiding this now that we have rewinding options, since detailed args aren't really in the spirit of this subcommand
//    << "  --gbwt-buffer-size NUM GBWT construction buffer size in millions of nodes; may need to be" << endl
//    << "                              increased for graphs with long haplotypes "
//                                   << "[" << IndexingParameters::gbwt_insert_batch_size / gbwt::MILLION << "]" << endl
//    << "  --gcsa-size-limit NUM limit on size of GCSA2 temporary files on disk in bytes" << endl
         << "  -t, --threads NUM      number of threads [all available]" << endl
         << "  -V, --verbosity NUM    log to stderr {0 = none, 1 = basic, 2 = debug}"
         <<                           "[" << (int) IndexingParameters::verbosity << "]" << endl
       //<< "  -d, --dot              print the dot-formatted graph of index recipes and exit" << endl
         << "  -h, --help             print this help message to stderr and exit" << endl;
}

int main_autoindex(int argc, char** argv) {
    
    if (argc == 2) {
        help_autoindex(argv);
        return 1;
    }
    
    constexpr int OPT_KEEP_INTERMEDIATE =  1000;
    constexpr int OPT_FORCE_UNPHASED = 1001;
    constexpr int OPT_FORCE_PHASED = 1002;
    constexpr int OPT_GBWT_BUFFER_SIZE = 1003;
    constexpr int OPT_GCSA_SIZE_LIMIT = 1004;
    
    // load the registry
    IndexRegistry registry = VGIndexes::get_vg_index_registry();
    bool print_dot = false;
    vector<IndexName> targets;
    vector<string> vcf_names;
    bool force_unphased = false;
    bool force_phased = false;
    int64_t target_mem_usage = IndexRegistry::get_system_memory() / 2;
    bool allow_guessing = true;
    
    string gfa_name;
    string gbz_name;
    
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
            {"gbz", required_argument, 0, 'G'},
            {"tx-gff", required_argument, 0, 'x'},
            {"hap-tx-gff", required_argument, 0, 'H'},
            {"no-guessing", no_argument, 0, 'n'},
            {"gff-feature", required_argument, 0, 'f'},
            {"gff-tx-tag", required_argument, 0, 'a'},
            {"provide", required_argument, 0, 'P'},
            {"request", required_argument, 0, 'R'},
            {"target-mem", required_argument, 0, 'M'},
            {"gbwt-buffer-size", required_argument, 0, OPT_GBWT_BUFFER_SIZE},
            {"gcsa-size-limit", required_argument, 0, OPT_GCSA_SIZE_LIMIT},
            {"tmp-dir", required_argument, 0, 'T'},
            {"threads", required_argument, 0, 't'},
            {"verbosity", required_argument, 0, 'V'},
            {"dot", no_argument, 0, 'd'},
            {"help", no_argument, 0, 'h'},
            {"keep-intermediate", no_argument, 0, OPT_KEEP_INTERMEDIATE},
            {"force-unphased", no_argument, 0, OPT_FORCE_UNPHASED},
            {"force-phased", no_argument, 0, OPT_FORCE_PHASED},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "p:w:r:v:i:g:G:x:H:na:P:R:f:M:T:t:dV:h?",
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
                        targets.emplace_back(std::move(target));
                    }
                }
                else if (optarg == string("mpmap")) {
                    for (auto& target : VGIndexes::get_default_mpmap_indexes()) {
                        targets.emplace_back(std::move(target));
                    }
                }
                else if (optarg == string("giraffe") || optarg == string("sr-giraffe")) {
                    for (auto& target : VGIndexes::get_default_short_giraffe_indexes()) {
                        targets.emplace_back(std::move(target));
                    }
                }
                else if (optarg == string("lr-giraffe")) {
                    for (auto& target : VGIndexes::get_default_long_giraffe_indexes()) {
                        targets.emplace_back(std::move(target));
                    }
                }
                else if (optarg == string("rpvg")) {
                    for (auto& target : VGIndexes::get_default_rpvg_indexes()) {
                        targets.emplace_back(std::move(target));
                    }
                }
                else {
                    cerr << "error: Unrecognized workflow (-w): " << optarg << endl;
                    return 1;
                }
                break;
            case 'r':
                if (ends_with(optarg, vg::GZ_SUFFIX)) {
                        cerr << "[autoindex] ERROR: FASTA file " << optarg
                             << " appears to be gzipped. Decompress it before use." << endl;    
                    return 1;
                }
                registry.provide("Reference FASTA", optarg);
                break;
            case 'v':
                vcf_names.push_back(optarg);
                break;
            case 'i':
                if (ends_with(optarg, vg::GZ_SUFFIX)) {
                        cerr << "[autoindex] ERROR: FASTA file " << optarg
                             << " appears to be gzipped. Decompress it before use." << endl;    
                    return 1;
                }
                registry.provide("Insertion Sequence FASTA", optarg);
                break;
            case 'g':
                gfa_name = optarg;
                break;
            case 'G':
                gbz_name = optarg;
                break;
            case 'x':
                registry.provide("GTF/GFF", optarg);
                break;
            case 'H':
                registry.provide("Haplotype GTF/GFF", optarg);
                break;
            case 'n':
                allow_guessing = false;
                break;
            case 'f':
                IndexingParameters::gff_feature_name = optarg;
                break;
            case 'a':
                IndexingParameters::gff_transcript_tag = optarg;
                break;
            case 'P':
            {
                auto parsed = parse_provide_string(optarg);
                registry.provide(parsed.first, parsed.second);
                break;
            }
            case 'R':
                targets.emplace_back(optarg);
                break;
            case 'M':
                target_mem_usage = parse_memory_usage(optarg);
                break;
            case OPT_GBWT_BUFFER_SIZE:
                IndexingParameters::gbwt_insert_batch_size = std::max(parse<size_t>(optarg), 1ul) * gbwt::MILLION;
                break;
            case 'T':
                temp_file::set_dir(optarg);
                break;
            case 't':
                omp_set_num_threads(parse<int>(optarg));
                break;
            case 'V':
            {
                int verbosity = parse<int>(optarg);
                if (verbosity < IndexingParameters::None || verbosity > IndexingParameters::Debug) {
                    cerr << "error: Verbosity (-V) must be integer in {0, 1, 2}: " << optarg << endl;
                    return 1;
                }
                IndexingParameters::verbosity = (IndexingParameters::Verbosity) verbosity;
                break;
            }
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
            case OPT_GCSA_SIZE_LIMIT:
                IndexingParameters::gcsa_size_limit = parse<int64_t>(optarg);
                break;
            case '?':
            case 'h':
            default:
                help_autoindex(argv);
                return 1;
        }
    }
    
    if (IndexingParameters::verbosity >= IndexingParameters::Basic) {
        cerr << "[vg autoindex] Executing command:";
        for (int i = 0; i < argc; ++i) {
            cerr << " " << argv[i];
        }
        cerr << endl;
    }
    
    assert(!(force_phased && force_unphased));
    
    // we have special logic for VCFs to make it friendly to both phased
    // and unphased VCF files
    if (!vcf_names.empty()) {
        // we interpret it as a phased VCF if any of the VCFs have phasing
        bool phased = force_phased;
        if (!force_unphased) {
            for (size_t i = 0; i < vcf_names.size() && !phased; ++i) {
                phased = IndexRegistry::vcf_is_phased(vcf_names[i]);
            }
        }
        
        for (auto& vcf_name : vcf_names) {
            if (phased) {
                registry.provide("VCF w/ Phasing", vcf_name);
            }
            else {
                registry.provide("VCF", vcf_name);
            }
        }
    }
    
    if (!gfa_name.empty()) {
        if (IndexRegistry::gfa_has_haplotypes(gfa_name)) {
            registry.provide("Reference GFA w/ Haplotypes", gfa_name);
        }
        else {
            registry.provide("Reference GFA", gfa_name);
        }
    }
    if (!gbz_name.empty()) {
        registry.provide("GBZ", gbz_name);
        // Also mark this as the Giraffe GBZ so if we index for Giraffe we won't duplicate it.
        // TODO: Why are these separate?
        registry.provide("Giraffe GBZ", gbz_name);
    }


    if (print_dot) {
        // don't index, just visualize the plan
        cout << registry.to_dot(targets);
        return 0;
    }
    
    registry.set_target_memory_usage(target_mem_usage);
    
    if (targets.empty()) {
        // default to vg map indexes
        targets = VGIndexes::get_default_map_indexes();
    }
    // deduplicate
    sort(targets.begin(), targets.end());
    targets.resize(unique(targets.begin(), targets.end()) - targets.begin());

    //Check if we can automatically load other indexes in the plan based on the names
    if (allow_guessing) {
        for (const IndexName& target : targets) {
            if (!registry.available(target)) {
                vector<string> inferred_file_names = registry.get_possible_filenames(target);
                for (const string& filename : inferred_file_names) {
                    if (ifstream(filename).is_open()) {
                        cerr << "[vg autoindex] Guessing that " << filename << " is " << target << endl;
                        registry.provide(target, filename);
                        break;
                    }
                }
            }
        }
    }
    
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
static Subcommand vg_autoindex("autoindex", "mapping tool-oriented index construction from interchange formats", 
                               PIPELINE, 1, main_autoindex);


/**
 * \file mpmap_main.cpp: multipath mapping of reads to a graph
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include "../multipath_mapper.hpp"
#include "../path.hpp"

//#define record_read_run_times

#ifdef record_read_run_times
#define READ_TIME_FILE "_read_times.tsv"
#include <ctime>
#include <iostream>
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_mpmap(char** argv) {
    cerr
    << "usage: " << argv[0] << " mpmap [options] -x index.xg -g index.gcsa [-f reads1.fq [-f reads2.fq] | -G reads.gam] > aln.gamp" << endl
    << "Multipath align reads to a graph." << endl
    << endl
    << "basic options:" << endl
    << "graph/index:" << endl
    << "  -x, --xg-name FILE        use this xg index (required)" << endl
    << "  -g, --gcsa-name FILE      use this GCSA2/LCP index pair (required; both FILE and FILE.lcp)" << endl
    << "  -H, --gbwt-name FILE      use this GBWT haplotype index for population-based MAPQs" << endl
    << "input:" << endl
    << "  -f, --fastq FILE          input FASTQ (possibly compressed), can be given twice for paired ends (for stdin use -)" << endl
    << "  -G, --gam-input FILE      input GAM (for stdin, use -)" << endl
    << "  -i, --interleaved         FASTQ or GAM contains interleaved paired ends" << endl
    << "  -N, --sample NAME         add this sample name to output GAM" << endl
    << "  -R, --read-group NAME     add this read group to output GAM" << endl
    << "  -e, --same-strand         read pairs are from the same strand of the DNA molecule" << endl
    << "algorithm:" << endl
    << "  -S, --single-path-mode    produce single-path alignments (GAM) instead of multipath alignments (GAMP) (ignores -sua)" << endl
    << "  -s, --snarls FILE         align to alternate paths in these snarls" << endl
    << "scoring:" << endl
    << "  -A, --no-qual-adjust      do not perform base quality adjusted alignments (required if input does not have base qualities)" << endl
    << endl
    << "advanced options:" << endl
    << "algorithm:" << endl
    << "  -X, --snarl-max-cut INT   do not align to alternate paths in a snarl if an exact match is at least this long (0 for no limit) [5]" << endl
    << "  -a, --alt-paths INT       align to (up to) this many alternate paths in between MEMs or in snarls [4]" << endl
    << "  -n, --unstranded          use lazy strand consistency when clustering MEMs" << endl
    << "  -b, --frag-sample INT     look for this many unambiguous mappings to estimate the fragment length distribution [1000]" << endl
    << "  -I, --frag-mean           mean for fixed fragment length distribution" << endl
    << "  -D, --frag-stddev         standard deviation for fixed fragment length distribution" << endl
    << "  -B, --no-calibrate        do not auto-calibrate mismapping dectection" << endl
    << "  -P, --max-p-val FLOAT     background model p value must be less than this to avoid mismapping detection [0.00001]" << endl
    << "  -v, --mq-method OPT       mapping quality method: 0 - none, 1 - fast approximation, 2 - adaptive, 3 - exact [2]" << endl
    << "  -Q, --mq-max INT          cap mapping quality estimates at this much [60]" << endl
    << "  -p, --band-padding INT    pad dynamic programming bands in inter-MEM alignment by this much [2]" << endl
    << "  -u, --map-attempts INT    perform (up to) this many mappings per read (0 for no limit) [48]" << endl
    << "  -O, --max-paths INT       consider (up to) this many paths per alignment when scoring by population consistency [1]" << endl
    << "  -M, --max-multimaps INT   report (up to) this many mappings per read [1]" << endl
    << "  -r, --reseed-length INT   reseed SMEMs for internal MEMs if they are at least this long (0 for no reseeding) [28]" << endl
    << "  -W, --reseed-diff FLOAT   require internal MEMs to have length within this much of the SMEM's length [0.45]" << endl
    << "  -k, --min-mem-length INT  minimum MEM length to anchor multipath alignments [1]" << endl
    << "  -K, --clust-length INT    minimum MEM length form clusters [automatic]" << endl
    << "  -c, --hit-max INT         ignore MEMs that occur greater than this many times in the graph (0 for no limit) [256]" << endl
    << "  -d, --max-dist-error INT  maximum typical deviation between distance on a reference path and distance in graph [8]" << endl
    << "  -w, --approx-exp FLOAT    let the approximate likelihood miscalculate likelihood ratios by this power [6.5]" << endl
    << "  -C, --drop-subgraph FLOAT drop alignment subgraphs whose MEMs cover this fraction less of the read than the best subgraph [0.2]" << endl
    << "  -U, --prune-exp FLOAT     prune MEM anchors if their approximate likelihood is this root less than the optimal anchors [1.25]" << endl
    << "scoring:" << endl
    << "  -q, --match INT           use this match score [1]" << endl
    << "  -z, --mismatch INT        use this mismatch penalty [4]" << endl
    << "  -o, --gap-open INT        use this gap open penalty [6]" << endl
    << "  -y, --gap-extend INT      use this gap extension penalty [1]" << endl
    << "  -L, --full-l-bonus INT    add this score to alignments that use the full length of the read [5]" << endl
    << "  -m, --remove-bonuses      remove full length alignment bonuses in reported scores" << endl
    << "computational parameters:" << endl
    << "  -t, --threads INT         number of compute threads to use" << endl
    << "  -Z, --buffer-size INT     buffer this many alignments together (per compute thread) before outputting to stdout [100]" << endl;
    
}

int main_mpmap(int argc, char** argv) {

    if (argc == 2) {
        help_mpmap(argv);
        return 1;
    }

    // initialize parameters with their default options
    string xg_name;
    string gcsa_name;
    string gbwt_name;
    string snarls_name;
    string fastq_name_1;
    string fastq_name_2;
    string gam_file_name;
    int match_score = default_match;
    int mismatch_score = default_mismatch;
    int gap_open_score = default_gap_open;
    int gap_extension_score = default_gap_extension;
    int full_length_bonus = 5;
    bool interleaved_input = false;
    int snarl_cut_size = 5;
    int max_map_attempts = 48;
    int population_max_paths = 1;
    int max_rescue_attempts = 32;
    int max_num_mappings = 1;
    int buffer_size = 100;
    int hit_max = 256;
    int min_mem_length = 1;
    int min_clustering_mem_length = 0;
    int reseed_length = 28;
    double reseed_diff = 0.45;
    double reseed_exp = 0.065;
    bool use_adaptive_reseed = true;
    double cluster_ratio = 0.2;
    bool qual_adjusted = true;
    bool strip_full_length_bonus = false;
    MappingQualityMethod mapq_method = Adaptive;
    int band_padding = 2;
    int max_dist_error = 8;
    int num_alt_alns = 4;
    double suboptimal_path_exponent = 1.25;
    double likelihood_approx_exp = 6.5;
    bool single_path_alignment_mode = false;
    int max_mapq = 60;
    size_t frag_length_sample_size = 1000;
    double frag_length_robustness_fraction = 0.95;
    double frag_length_mean = NAN;
    double frag_length_stddev = NAN;
    bool same_strand = false;
    bool auto_calibrate_mismapping_detection = true;
    double max_mapping_p_value = 0.00001;
    size_t num_calibration_simulations = 250;
    size_t calibration_read_length = 150;
    bool unstranded_clustering = false;
    size_t order_length_repeat_hit_max = 3000;
    size_t sub_mem_count_thinning = 4;
    size_t sub_mem_thinning_burn_in = 16;
    double secondary_rescue_score_diff = 0.8;
    size_t rescue_only_min = numeric_limits<size_t>::max(); // disabling this for now
    size_t rescue_only_anchor_max = 16;
    string sample_name = "";
    string read_group = "";
    bool prefilter_redundant_hits = true;
    bool precollapse_order_length_hits = true;
    int max_sub_mem_recursion_depth = 1;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gcsa-name", required_argument, 0, 'g'},
            {"gbwt-name", required_argument, 0, 'H'},
            {"fastq", required_argument, 0, 'f'},
            {"gam-input", required_argument, 0, 'G'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"interleaved", no_argument, 0, 'i'},
            {"same-strand", no_argument, 0, 'e'},
            {"single-path-mode", no_argument, 0, 'S'},
            {"snarls", required_argument, 0, 's'},
            {"snarl-max-cut", required_argument, 0, 'X'},
            {"alt-paths", required_argument, 0, 'a'},
            {"unstranded", no_argument, 0, 'n'},
            {"frag-sample", required_argument, 0, 'b'},
            {"frag-mean", required_argument, 0, 'I'},
            {"frag-stddev", required_argument, 0, 'D'},
            {"no-calibrate", no_argument, 0, 'B'},
            {"max-p-val", required_argument, 0, 'P'},
            {"mq-method", required_argument, 0, 'v'},
            {"mq-max", required_argument, 0, 'Q'},
            {"band-padding", required_argument, 0, 'p'},
            {"map-attempts", required_argument, 0, 'u'},
            {"max-paths", required_argument, 0, 'O'},
            {"max-multimaps", required_argument, 0, 'M'},
            {"reseed-length", required_argument, 0, 'r'},
            {"reseed-diff", required_argument, 0, 'W'},
            {"min-mem-length", required_argument, 0, 'k'},
            {"clustlength", required_argument, 0, 'K'},
            {"hit-max", required_argument, 0, 'c'},
            {"max-dist-error", required_argument, 0, 'd'},
            {"approx-exp", required_argument, 0, 'w'},
            {"drop-subgraph", required_argument, 0, 'C'},
            {"prune-exp", required_argument, 0, 'U'},
            {"match", required_argument, 0, 'q'},
            {"mismatch", required_argument, 0, 'z'},
            {"gap-open", required_argument, 0, 'o'},
            {"gap-extend", required_argument, 0, 'y'},
            {"full-l-bonus", required_argument, 0, 'L'},
            {"remove-bonuses", no_argument, 0, 'm'},
            {"no-qual-adjust", no_argument, 0, 'A'},
            {"threads", required_argument, 0, 't'},
            {"buffer-size", required_argument, 0, 'Z'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:H:f:G:N:R:ieSs:u:O:a:nb:I:D:BP:v:Q:p:M:r:W:k:K:c:d:w:C:R:q:z:o:y:L:mAt:Z:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                if (xg_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide XG file with -x." << endl;
                    exit(1);
                }
                break;
                
            case 'g':
                gcsa_name = optarg;
                if (gcsa_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide GCSA file with -g." << endl;
                    exit(1);
                }
                break;
                
            case 'H':
                gbwt_name = optarg;
                break;
                
            case 'f':
                if (fastq_name_1.empty()) {
                    fastq_name_1 = optarg;
                    if (fastq_name_1.empty()) {
                        cerr << "error:[vg mpmap] Must provide FASTQ file with -f" << endl;
                        exit(1);
                    }
                }
                else if (fastq_name_2.empty()) {
                    fastq_name_2 = optarg;
                    if (fastq_name_2.empty()) {
                        cerr << "error:[vg mpmap] Must provide FASTQ file with -f" << endl;
                        exit(1);
                    }
                }
                else {
                    cerr << "error:[vg mpmap] Cannot specify more than two FASTQ files" << endl;
                    exit(1);
                }
                break;
                
            case 'G':
                gam_file_name = optarg;
                if (gam_file_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide GAM file with -G." << endl;
                    exit(1);
                }
                break;
                
            case 'N':
                sample_name = optarg;
                if (sample_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide sample name file with -N." << endl;
                    exit(1);
                }
                break;
                
            case 'R':
                read_group = optarg;
                if (read_group.empty()) {
                    cerr << "error:[vg mpmap] Must provide read group with -R." << endl;
                    exit(1);
                }
                break;
                
            case 'i':
                interleaved_input = true;
                break;
                
            case 'e':
                same_strand = true;
                break;
                
            case 'S':
                single_path_alignment_mode = true;
                break;
                
            case 's':
                snarls_name = optarg;
                if (snarls_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide snarl file with -s." << endl;
                    exit(1);
                }
                break;
                
            case 'X':
                snarl_cut_size = atoi(optarg);
                break;
                
            case 'a':
                num_alt_alns = atoi(optarg);
                break;
                
            case 'n':
                unstranded_clustering = true;
                break;
                
            case 'b':
                frag_length_sample_size = atoi(optarg);
                break;
                
            case 'I':
                frag_length_mean = atof(optarg);
                break;
                
            case 'D':
                frag_length_stddev = atof(optarg);
                break;
                
            case 'B':
                auto_calibrate_mismapping_detection = false;
                break;
                
            case 'P':
                max_mapping_p_value = atof(optarg);
                break;
                
            case 'v':
            {
                int mapq_arg = atoi(optarg);
                if (mapq_arg == 0) {
                    mapq_method = None;
                }
                else if (mapq_arg == 1) {
                    mapq_method = Approx;
                }
                else if (mapq_arg == 2) {
                    mapq_method = Adaptive;
                }
                else if (mapq_arg == 3) {
                    mapq_method = Exact;
                }
                else {
                    cerr << "error:[vg mpmap] Unrecognized mapping quality (-v) option: " << mapq_arg << ". Choose from {0, 1, 2}." << endl;
                    exit(1);
                }
            }
                break;
                
            case 'Q':
                max_mapq = atoi(optarg);
                break;
                
            case 'p':
                band_padding = atoi(optarg);
                break;
                
            case 'u':
                max_map_attempts = atoi(optarg);
                break;
                
            case 'O':
                population_max_paths = atoi(optarg);
                break;
                
            case 'M':
                max_num_mappings = atoi(optarg);
                break;
                
            case 'r':
                reseed_length = atoi(optarg);
                break;
                
            case 'W':
                reseed_diff = atof(optarg);
                break;
                
            case 'k':
                min_mem_length = atoi(optarg);
                break;
                
            case 'K':
                min_clustering_mem_length = atoi(optarg);
                break;
                
            case 'c':
                hit_max = atoi(optarg);
                break;
                
            case 'd':
                max_dist_error = atoi(optarg);
                break;
                
            case 'w':
                likelihood_approx_exp = atof(optarg);
                break;
                
            case 'C':
                cluster_ratio = atof(optarg);
                break;
                
            case 'U':
                suboptimal_path_exponent = atof(optarg);
                break;
                
            case 'q':
                match_score = atoi(optarg);
                break;
                
            case 'z':
                mismatch_score = atoi(optarg);
                break;
                
            case 'o':
                gap_open_score = atoi(optarg);
                break;
                
            case 'y':
                gap_extension_score = atoi(optarg);
                break;
                
            case 'm':
                strip_full_length_bonus = true;
                break;
                
            case 'L':
                full_length_bonus = atoi(optarg);
                break;
                
            case 'A':
                qual_adjusted = false;
                break;
                
            case 't':
            {
                int num_threads = atoi(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg mpmap] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'Z':
                buffer_size = atoi(optarg);
                break;
                
            case 'h':
            case '?':
            default:
                help_mpmap(argv);
                exit(1);
                break;
        }
    }
    
    // check for valid parameters
    
    if (std::isnan(frag_length_mean) != std::isnan(frag_length_stddev)) {
        cerr << "error:[vg mpmap] Cannot specify only one of fragment length mean (-I) and standard deviation (-D)." << endl;
        exit(1);
    }
    
    if (!std::isnan(frag_length_mean) && frag_length_mean < 0) {
        cerr << "error:[vg mpmap] Fragment length mean (-I) must be nonnegative." << endl;
        exit(1);
    }
    
    if (!std::isnan(frag_length_stddev) && frag_length_stddev < 0) {
        cerr << "error:[vg mpmap] Fragment length standard deviation (-D) must be nonnegative." << endl;
        exit(1);
    }
    
    if (interleaved_input && !fastq_name_2.empty()) {
        cerr << "error:[vg mpmap] Cannot designate both interleaved paired ends (-i) and separate paired end file (-f)." << endl;
        exit(1);
    }
    
    if (!fastq_name_1.empty() && !gam_file_name.empty()) {
        cerr << "error:[vg mpmap] Cannot designate both FASTQ input (-f) and GAM input (-G) in same run." << endl;
        exit(1);
    }
    
    if (fastq_name_1.empty() && gam_file_name.empty()) {
        cerr << "error:[vg mpmap] Must designate reads to map from either FASTQ (-f) or GAM (-G) file." << endl;
        exit(1);
    }
    
    if (!interleaved_input && fastq_name_2.empty() && same_strand) {
        cerr << "warning:[vg mpmap] Ignoring same strand parameter (-d) because no paired end input provided." << endl;
    }
    
    if (num_alt_alns <= 0) {
        cerr << "error:[vg mpmap] Number of alternate snarl paths (-a) set to " << num_alt_alns << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (frag_length_sample_size <= 0) {
        cerr << "error:[vg mpmap] Fragment length distribution sample size (-b) set to " << frag_length_sample_size << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (snarl_cut_size < 0) {
        cerr << "error:[vg mpmap] Max snarl cut size (-U) set to " << snarl_cut_size << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (max_mapq <= 0 && mapq_method != None) {
        cerr << "error:[vg mpmap] Maximum mapping quality (-Q) set to " << max_mapq << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (band_padding < 0) {
        cerr << "error:[vg mpmap] Band padding (-p) set to " << band_padding << ", must set to a nonnegative integer." << endl;
        exit(1);
    }
    
    if (max_map_attempts < 0) {
        cerr << "error:[vg mpmap] Maximum number of mapping attempts (-u) set to " << max_map_attempts << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (population_max_paths < 1) {
        cerr << "error:[vg mpmap] Maximum number of paths per alignment for population scoring (-O) set to " << population_max_paths << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (population_max_paths != 1 && gbwt_name.empty()) {
        cerr << "error:[vg mpmap] Maximum number of paths per alignment for population scoring (-O) set when population database (-H) not provided." << endl;
        exit(1);
    }
    
    if (max_num_mappings > max_map_attempts && max_map_attempts != 0) {
        cerr << "warning:[vg mpmap] Reporting up to " << max_num_mappings << " mappings, but only computing up to " << max_map_attempts << " mappings." << endl;
    }
    
    if (max_num_mappings <= 0) {
        cerr << "error:[vg mpmap] Maximum number of mappings per read (-M) set to " << max_num_mappings << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (reseed_length < 0) {
        cerr << "error:[vg mpmap] Reseeding length (-r) set to " << reseed_length << ", must set to a positive integer or 0 for no reseeding." << endl;
        exit(1);
    }
    
    if ((reseed_diff <= 0 || reseed_diff >= 1.0) && reseed_length != 0) {
        cerr << "error:[vg mpmap] Reseeding length difference (-W) set to " << reseed_diff << ", must set to a number between 0.0 and 1.0." << endl;
        exit(1);
    }
    
    if (reseed_exp < 0 && reseed_length != 0 && use_adaptive_reseed) {
        cerr << "error:[vg mpmap] Reseeding exponent set to " << reseed_exp << ",  must set to a nonnegative number." << endl;
        exit(1);
    }
    
    if (hit_max < 0) {
        cerr << "error:[vg mpmap] MEM hit max (-c) set to " << hit_max << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (max_dist_error < 0) {
        cerr << "error:[vg mpmap] Maximum distance approximation error (-d) set to " << max_dist_error << ", must set to a nonnegative integer." << endl;
        exit(1);
    }
    
    if (likelihood_approx_exp < 1.0) {
        cerr << "error:[vg mpmap] Likelihood approximation exponent (-w) set to " << likelihood_approx_exp << ", must set to at least 1.0." << endl;
        exit(1);
    }
    
    if (cluster_ratio < 0.0 || cluster_ratio >= 1.0) {
        cerr << "error:[vg mpmap] Cluster drop ratio (-C) set to " << cluster_ratio << ", must set to a number between 0.0 and 1.0." << endl;
        exit(1);
    }
    
    if (suboptimal_path_exponent < 1.0) {
        cerr << "error:[vg mpmap] Suboptimal path likelihood root (-R) set to " << suboptimal_path_exponent << ", must set to at least 1.0." << endl;
        exit(1);
    }
    
    if (min_mem_length <= 0) {
        cerr << "error:[vg mpmap] Minimum MEM length (-k) set to " << min_mem_length << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (match_score > std::numeric_limits<int8_t>::max() || mismatch_score > std::numeric_limits<int8_t>::max()
        || gap_open_score > std::numeric_limits<int8_t>::max() || gap_extension_score > std::numeric_limits<int8_t>::max()
        || full_length_bonus > std::numeric_limits<int8_t>::max() || match_score < 0 || mismatch_score < 0
        || gap_open_score < 0 || gap_extension_score < 0 || full_length_bonus < 0) {
        cerr << "error:[vg mpmap] All alignment scoring parameters (-qzoyL) must be between 0 and " << (int) std::numeric_limits<int8_t>::max() << endl;
        exit(1);
    }
    
    if (buffer_size <= 0) {
        cerr << "error:[vg mpmap] Buffer size (-Z) set to " << buffer_size << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    // adjust parameters that produce irrelevant extra work in single path mode
    
    if (single_path_alignment_mode && population_max_paths == 1) {
        // TODO: I don't like having these constants floating around in two different places, but it's not very risky, just a warning
        if (!snarls_name.empty()) {
            cerr << "warning:[vg mpmap] Snarl file (-s) is ignored in single path mode (-S) without multipath population scoring (-O)." << endl;
        }
        
        if (snarl_cut_size != 5) {
            cerr << "warning:[vg mpmap] Snarl cut limit (-u) is ignored in single path mode (-S) without multipath population scoring (-O)." << endl;
        }
        
        if (num_alt_alns != 4) {
            cerr << "warning:[vg mpmap] Number of alternate alignments (-a) is ignored in single path mode (-S) without multipath population scoring (-O)." << endl;
        }
        num_alt_alns = 1;
        
    }
    
    // ensure required parameters are provided
    
    if (xg_name.empty()) {
        cerr << "error:[vg mpmap] Multipath mapping requires an XG index, must provide XG file" << endl;
        exit(1);
    }
    
    if (gcsa_name.empty()) {
        cerr << "error:[vg mpmap] Multipath mapping requires a GCSA2 index, must provide GCSA2 file" << endl;
        exit(1);
    }
    
    // create in-memory objects
    
    ifstream xg_stream(xg_name);
    if (!xg_stream) {
        cerr << "error:[vg mpmap] Cannot open XG file " << xg_name << endl;
        exit(1);
    }
    
    ifstream gcsa_stream(gcsa_name);
    if (!xg_stream) {
        cerr << "error:[vg mpmap] Cannot open GCSA2 file " << gcsa_name << endl;
        exit(1);
    }
    
    string lcp_name = gcsa_name + ".lcp";
    ifstream lcp_stream(lcp_name);
    if (!xg_stream) {
        cerr << "error:[vg mpmap] Cannot open LCP file " << lcp_name << endl;
        exit(1);
    }
    
    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Configure its temp directory to the system temp directory
    gcsa::TempFile::setDirectory(find_temp_dir());
    
    xg::XG xg_index(xg_stream);
    gcsa::GCSA gcsa_index;
    gcsa_index.load(gcsa_stream);
    gcsa::LCPArray lcp_array;
    lcp_array.load(lcp_stream);
    
    gbwt::GBWT* gbwt = nullptr;
    haplo::ScoreProvider* haplo_score_provider = nullptr;
    if (!gbwt_name.empty()) {
        ifstream gbwt_stream(gbwt_name);
        if (!gbwt_stream) {
            cerr << "error:[vg mpmap] Cannot open GBWT file " << gbwt_name << endl;
            exit(1);
        }
        gbwt = new gbwt::GBWT();
        gbwt->load(gbwt_stream);
        
        // We have the GBWT available for scoring haplotypes
        haplo_score_provider = new haplo::GBWTScoreProvider<gbwt::GBWT>(*gbwt);
    }
    // TODO: Allow using haplo::XGScoreProvider?
    
    SnarlManager* snarl_manager = nullptr;
    if (!snarls_name.empty()) {
        ifstream snarl_stream(snarls_name);
        if (!snarl_stream) {
            cerr << "error:[vg mpmap] Cannot open Snarls file " << snarls_name << endl;
            exit(1);
        }
        snarl_manager = new SnarlManager(snarl_stream);
    }
        
    MultipathMapper multipath_mapper(&xg_index, &gcsa_index, &lcp_array, haplo_score_provider, snarl_manager);
    
    // set alignment parameters
    multipath_mapper.set_alignment_scores(match_score, mismatch_score, gap_open_score, gap_extension_score, full_length_bonus);
    multipath_mapper.adjust_alignments_for_base_quality = qual_adjusted;
    multipath_mapper.strip_bonuses = strip_full_length_bonus;
    multipath_mapper.band_padding = band_padding;
    
    // set mem finding parameters
    multipath_mapper.hit_max = hit_max;
    multipath_mapper.mem_reseed_length = reseed_length;
    multipath_mapper.fast_reseed = true;
    multipath_mapper.fast_reseed_length_diff = reseed_diff;
    multipath_mapper.sub_mem_count_thinning = sub_mem_count_thinning;
    multipath_mapper.sub_mem_thinning_burn_in = sub_mem_thinning_burn_in;
    multipath_mapper.order_length_repeat_hit_max = order_length_repeat_hit_max;
    multipath_mapper.min_mem_length = min_mem_length;
    multipath_mapper.adaptive_reseed_diff = use_adaptive_reseed;
    multipath_mapper.adaptive_diff_exponent = reseed_exp;
    multipath_mapper.use_approx_sub_mem_count = false;
    multipath_mapper.prefilter_redundant_hits = prefilter_redundant_hits;
    multipath_mapper.precollapse_order_length_hits = precollapse_order_length_hits;
    multipath_mapper.max_sub_mem_recursion_depth = max_sub_mem_recursion_depth;
    multipath_mapper.max_mapping_p_value = max_mapping_p_value;
    if (min_clustering_mem_length) {
        multipath_mapper.min_clustering_mem_length = min_clustering_mem_length;
    }
    else {
        multipath_mapper.set_automatic_min_clustering_length();
    }
    
    // set mapping quality parameters
    multipath_mapper.mapping_quality_method = mapq_method;
    multipath_mapper.max_mapping_quality = max_mapq;
    multipath_mapper.use_population_mapqs = (haplo_score_provider != nullptr);
    multipath_mapper.population_max_paths = population_max_paths;
    
    // set pruning and clustering parameters
    multipath_mapper.max_expected_dist_approx_error = max_dist_error;
    multipath_mapper.mem_coverage_min_ratio = cluster_ratio;
    multipath_mapper.log_likelihood_approx_factor = likelihood_approx_exp;
    multipath_mapper.num_mapping_attempts = max_map_attempts ? max_map_attempts : numeric_limits<int>::max();
    multipath_mapper.unstranded_clustering = unstranded_clustering;
    
    // set pair rescue parameters
    multipath_mapper.secondary_rescue_score_diff = secondary_rescue_score_diff;
    multipath_mapper.max_rescue_attempts = max_rescue_attempts;
    multipath_mapper.rescue_only_min = rescue_only_min;
    multipath_mapper.rescue_only_anchor_max = rescue_only_anchor_max;
    
    // set multipath alignment topology parameters
    multipath_mapper.max_snarl_cut_size = snarl_cut_size;
    multipath_mapper.num_alt_alns = num_alt_alns;
    multipath_mapper.max_suboptimal_path_score_ratio = suboptimal_path_exponent;
    
    // if directed to, auto calibrate the mismapping detection to the graph
    if (auto_calibrate_mismapping_detection) {
        multipath_mapper.calibrate_mismapping_detection(num_calibration_simulations, calibration_read_length);
    }
    
    // set computational paramters
    int thread_count = get_thread_count();
    multipath_mapper.set_alignment_threads(thread_count);
    
    // are we doing paired ends?
    if (interleaved_input || !fastq_name_2.empty()) {
        // make sure buffer size is even (ensures that output will be interleaved)
        if (buffer_size % 2 == 1) {
            buffer_size++;
        }

        if (!std::isnan(frag_length_mean) && !std::isnan(frag_length_stddev)) {
            // Force a fragment length distribution
            multipath_mapper.force_fragment_length_distr(frag_length_mean, frag_length_stddev);
        }
        else {
            // choose the sample size and tail-fraction for estimating the fragment length distribution
            multipath_mapper.set_fragment_length_distr_params(frag_length_sample_size, frag_length_sample_size,
                                                              frag_length_robustness_fraction);
        }

    }
    
#ifdef record_read_run_times
    ofstream read_time_file(READ_TIME_FILE);
#endif
    
    // a buffer to hold read pairs that can't be unambiguously mapped before the fragment length distribution
    // is estimated
    // note: sufficient to have only one buffer because multithreading code enforces single threaded mode
    // during distribution estimation
    vector<pair<Alignment, Alignment>> ambiguous_pair_buffer;
    
    vector<vector<Alignment> > single_path_output_buffer(thread_count);
    vector<vector<MultipathAlignment> > multipath_output_buffer(thread_count);
    
    // write unpaired multipath alignments to stdout buffer
    auto output_multipath_alignments = [&](vector<MultipathAlignment>& mp_alns) {
        auto& output_buf = multipath_output_buffer[omp_get_thread_num()];
        
        // move all the alignments over to the output buffer
        for (MultipathAlignment& mp_aln : mp_alns) {
            output_buf.emplace_back(move(mp_aln));
            
            // label with read group and sample name
            if (!read_group.empty()) {
                output_buf.back().set_read_group(read_group);
            }
            if (!sample_name.empty()) {
                output_buf.back().set_sample_name(sample_name);
            }
        }
        
        stream::write_buffered(cout, output_buf, buffer_size);
    };
    
    // convert to unpaired single path alignments and write stdout buffer
    auto output_single_path_alignments = [&](vector<MultipathAlignment>& mp_alns) {
        auto& output_buf = single_path_output_buffer[omp_get_thread_num()];
        // add optimal alignments to the output buffer
        for (MultipathAlignment& mp_aln : mp_alns) {
            output_buf.emplace_back();
            optimal_alignment(mp_aln, output_buf.back());
            // compute the Alignment identity to make vg call happy
            output_buf.back().set_identity(identity(output_buf.back().path()));
            
            // label with read group and sample name
            if (!read_group.empty()) {
                output_buf.back().set_read_group(read_group);
            }
            if (!sample_name.empty()) {
                output_buf.back().set_sample_name(sample_name);
            }
        }
        
        stream::write_buffered(cout, output_buf, buffer_size);
    };
    
    // write paired multipath alignments to stdout buffer
    auto output_multipath_paired_alignments = [&](vector<pair<MultipathAlignment, MultipathAlignment>>& mp_aln_pairs) {
        auto& output_buf = multipath_output_buffer[omp_get_thread_num()];
        
        // move all the alignments over to the output buffer
        for (pair<MultipathAlignment, MultipathAlignment>& mp_aln_pair : mp_aln_pairs) {
            output_buf.emplace_back(move(mp_aln_pair.first));
            
            // label with read group and sample name
            if (!read_group.empty()) {
                output_buf.back().set_read_group(read_group);
            }
            if (!sample_name.empty()) {
                output_buf.back().set_sample_name(sample_name);
            }
            
            // switch second read back to the opposite strand if necessary
            if (same_strand) {
                output_buf.emplace_back(move(mp_aln_pair.second));
            }
            else {
                output_buf.emplace_back();
                rev_comp_multipath_alignment(mp_aln_pair.second,
                                             [&](vg::id_t node_id) { return xg_index.node_length(node_id); },
                                             output_buf.back());
            }
            
            // label with read group and sample name
            if (!read_group.empty()) {
                output_buf.back().set_read_group(read_group);
            }
            if (!sample_name.empty()) {
                output_buf.back().set_sample_name(sample_name);
            }
        }
        
        stream::write_buffered(cout, output_buf, buffer_size);
    };
    
    // convert to paired single path alignments and write stdout buffer
    auto output_single_path_paired_alignments = [&](vector<pair<MultipathAlignment, MultipathAlignment>>& mp_aln_pairs) {
        auto& output_buf = single_path_output_buffer[omp_get_thread_num()];
        
        // add optimal alignments to the output buffer
        for (pair<MultipathAlignment, MultipathAlignment>& mp_aln_pair : mp_aln_pairs) {
            
            output_buf.emplace_back();
            optimal_alignment(mp_aln_pair.first, output_buf.back());
            // compute the Alignment identity to make vg call happy
            output_buf.back().set_identity(identity(output_buf.back().path()));
            
            // label with read group and sample name
            if (!read_group.empty()) {
                output_buf.back().set_read_group(read_group);
            }
            if (!sample_name.empty()) {
                output_buf.back().set_sample_name(sample_name);
            }
            // arbitrarily decide that this is the "previous" fragment
            output_buf.back().mutable_fragment_next()->set_name(mp_aln_pair.second.name());
            
            output_buf.emplace_back();
            optimal_alignment(mp_aln_pair.second, output_buf.back());
            // compute identity again
            output_buf.back().set_identity(identity(output_buf.back().path()));
            
            // switch second read back to the opposite strand if necessary
            if (!same_strand) {
                reverse_complement_alignment_in_place(&output_buf.back(),
                                                      [&](vg::id_t node_id) { return xg_index.node_length(node_id); });
            }
            
            // label with read group and sample name
            if (!read_group.empty()) {
                output_buf.back().set_read_group(read_group);
            }
            if (!sample_name.empty()) {
                output_buf.back().set_sample_name(sample_name);
            }
            // arbitrarily decide that this is the "next" fragment
            output_buf.back().mutable_fragment_prev()->set_name(mp_aln_pair.first.name());
        }
        stream::write_buffered(cout, output_buf, buffer_size);
    };
    
    // do unpaired multipath alignment and write to buffer
    function<void(Alignment&)> do_unpaired_alignments = [&](Alignment& alignment) {
#ifdef record_read_run_times
        clock_t start = clock();
#endif
        vector<MultipathAlignment> mp_alns;
        multipath_mapper.multipath_map(alignment, mp_alns, max_num_mappings);
        if (single_path_alignment_mode) {
            output_single_path_alignments(mp_alns);
        }
        else {
            output_multipath_alignments(mp_alns);
        }
#ifdef record_read_run_times
        clock_t finish = clock();
#pragma omp critical
        read_time_file << alignment.name() << "\t" << double(finish - start) / CLOCKS_PER_SEC << endl;
#endif
    };
    
    // do paired multipath alignment and write to buffer
    function<void(Alignment&, Alignment&)> do_paired_alignments = [&](Alignment& alignment_1, Alignment& alignment_2) {
        // get reads on the same strand so that oriented distance estimation works correctly
        // but if we're clearing the ambiguous buffer we already RC'd these on the first pass
#ifdef record_read_run_times
        clock_t start = clock();
#endif
        if (!same_strand) {
            // remove the path so we won't try to RC it (the path may not refer to this graph)
            alignment_2.clear_path();
            reverse_complement_alignment_in_place(&alignment_2, [&](vg::id_t node_id) { return xg_index.node_length(node_id); });
        }
                
        vector<pair<MultipathAlignment, MultipathAlignment>> mp_aln_pairs;
        multipath_mapper.multipath_map_paired(alignment_1, alignment_2, mp_aln_pairs, ambiguous_pair_buffer, max_num_mappings);
        if (single_path_alignment_mode) {
            output_single_path_paired_alignments(mp_aln_pairs);
        }
        else {
            output_multipath_paired_alignments(mp_aln_pairs);
        }
#ifdef record_read_run_times
        clock_t finish = clock();
#pragma omp critical
        read_time_file << alignment_1.name() << "\t" << alignment_2.name() << "\t" << double(finish - start) / CLOCKS_PER_SEC << endl;
#endif
    };
    
    // for streaming paired input, don't spawn parallel tasks unless this evalutes to true
    function<bool(void)> multi_threaded_condition = [&](void) {
        return multipath_mapper.has_fixed_fragment_length_distr();
    };
    
    
    // FASTQ input
    if (!fastq_name_1.empty()) {
        if (interleaved_input) {
            fastq_paired_interleaved_for_each_parallel_after_wait(fastq_name_1, do_paired_alignments,
                                                                  multi_threaded_condition);
        }
        else if (fastq_name_2.empty()) {
            fastq_unpaired_for_each_parallel(fastq_name_1, do_unpaired_alignments);
        }
        else {
            fastq_paired_two_files_for_each_parallel_after_wait(fastq_name_1, fastq_name_2, do_paired_alignments,
                                                                multi_threaded_condition);
        }
    }
    
    // GAM input
    if (!gam_file_name.empty()) {
        function<void(istream&)> execute = [&](istream& gam_in) {
            if (!gam_in) {
                cerr << "error:[vg mpmap] Cannot open GAM file " << gam_file_name << endl;
                exit(1);
            }
            if (interleaved_input) {
                stream::for_each_interleaved_pair_parallel_after_wait(gam_in, do_paired_alignments,
                                                                      multi_threaded_condition);
            }
            else {
                stream::for_each_parallel(gam_in, do_unpaired_alignments);
            }
        };
        get_input_file(gam_file_name, execute);
    }

    // take care of any read pairs that we couldn't map unambiguously before the fragment length distribution
    // had been estimated
    if (!ambiguous_pair_buffer.empty()) {
        if (multipath_mapper.has_fixed_fragment_length_distr()) {
#pragma omp parallel for
            for (size_t i = 0; i < ambiguous_pair_buffer.size(); i++) {
                pair<Alignment, Alignment>& aln_pair = ambiguous_pair_buffer[i];
                // we reverse complemented the alignment on the first pass, so switch back so we don't break
                // the alignment functions expectations
                // TODO: slightly wasteful, inelegant
                if (!same_strand) {
                    reverse_complement_alignment_in_place(&aln_pair.second,
                                                          [&](vg::id_t node_id) { return xg_index.node_length(node_id); });
                }
                do_paired_alignments(aln_pair.first, aln_pair.second);
            }
        }
        else {
            cerr << "warning:[vg mpmap] Could not find " << frag_length_sample_size << " unambiguous read pair mappings to estimate fragment length ditribution. Mapping read pairs as independent single-ended reads. Consider decreasing sample size (-b)." << endl;
            
#pragma omp parallel for
            for (size_t i = 0; i < ambiguous_pair_buffer.size(); i++) {
                pair<Alignment, Alignment>& aln_pair = ambiguous_pair_buffer[i];
                // we reverse complemented the alignment on the first pass, so switch back so we don't break
                // the alignment function's expectations
                // TODO: slightly wasteful, inelegant
                if (!same_strand) {
                    reverse_complement_alignment_in_place(&aln_pair.second,
                                                          [&](vg::id_t node_id) { return xg_index.node_length(node_id); });
                }
                do_unpaired_alignments(aln_pair.first);
                do_unpaired_alignments(aln_pair.second);
            }
        }
    }
    
    // flush output buffers
    for (int i = 0; i < thread_count; i++) {
        vector<Alignment>& single_path_buffer = single_path_output_buffer[i];
        stream::write_buffered(cout, single_path_buffer, 0);
        
        vector<MultipathAlignment>& multipath_buffer = multipath_output_buffer[i];
        stream::write_buffered(cout, multipath_buffer, 0);
    }
    cout.flush();
    
#ifdef record_read_run_times
    read_time_file.close();
#endif
    
    //cerr << "MEM length filtering efficiency: " << ((double) OrientedDistanceClusterer::MEM_FILTER_COUNTER) / OrientedDistanceClusterer::MEM_TOTAL << " (" << OrientedDistanceClusterer::MEM_FILTER_COUNTER << "/" << OrientedDistanceClusterer::MEM_TOTAL << ")" << endl;
    //cerr << "MEM cluster filtering efficiency: " << ((double) OrientedDistanceClusterer::PRUNE_COUNTER) / OrientedDistanceClusterer::CLUSTER_TOTAL << " (" << OrientedDistanceClusterer::PRUNE_COUNTER << "/" << OrientedDistanceClusterer::CLUSTER_TOTAL << ")" << endl;
    //cerr << "subgraph filtering efficiency: " << ((double) MultipathMapper::PRUNE_COUNTER) / MultipathMapper::SUBGRAPH_TOTAL << " (" << MultipathMapper::PRUNE_COUNTER << "/" << MultipathMapper::SUBGRAPH_TOTAL << ")" << endl;
    
    if (snarl_manager != nullptr) {
        delete snarl_manager;
    }
   
    if (haplo_score_provider != nullptr) {
        delete haplo_score_provider;
    }
   
    if (gbwt != nullptr) {
        delete gbwt;
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_mpmap("mpmap", "multipath alignments of reads to a graph", main_mpmap);



    /**
 * \file mpmap_main.cpp: multipath mapping of reads to a graph
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include "../multipath_mapper.hpp"
#include "../path.hpp"
#include "../watchdog.hpp"

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
    << "  -x, --xg-name FILE            use this xg index (required)" << endl
    << "  -g, --gcsa-name FILE          use this GCSA2/LCP index pair (required; both FILE and FILE.lcp)" << endl
    << "  -H, --gbwt-name FILE          use this GBWT haplotype index for population-based MAPQs" << endl
    << "  -d, --dist-name FILE          use this snarl distance index for clustering (if given, requires matched snarls from -s)" << endl
    << "      --linear-index FILE       use this sublinear Li and Stephens index file for population-based MAPQs" << endl
    << "      --linear-path PATH        use the given path name as the path that the linear index is against" << endl
    << "input:" << endl
    << "  -f, --fastq FILE              input FASTQ (possibly compressed), can be given twice for paired ends (for stdin use -)" << endl
    << "  -G, --gam-input FILE          input GAM (for stdin, use -)" << endl
    << "  -i, --interleaved             FASTQ or GAM contains interleaved paired ends" << endl
    << "  -N, --sample NAME             add this sample name to output GAMP" << endl
    << "  -R, --read-group NAME         add this read group to output GAMP" << endl
    << "  -e, --same-strand             read pairs are from the same strand of the DNA molecule" << endl
    << "algorithm:" << endl
    << "  -S, --single-path-mode        produce single-path alignments (GAM) instead of multipath alignments (GAMP) (ignores -sua)" << endl
    << "  -s, --snarls FILE             align to alternate paths in these snarls" << endl
    << "scoring:" << endl
    << "  -A, --no-qual-adjust          do not perform base quality adjusted alignments (required if input does not have base qualities)" << endl
    << "  -E, --long-read-scoring       set alignment scores to long-read defaults: -q1 -z1 -o1 -y1 -L0 (can be overridden)" << endl
    << endl
    << "advanced options:" << endl
    << "algorithm:" << endl
    << "  -v, --tvs-clusterer           use the target value search based clusterer (requies a distance index from -d)" << endl
    << "  -X, --snarl-max-cut INT       do not align to alternate paths in a snarl if an exact match is at least this long (0 for no limit) [5]" << endl
    << "  -a, --alt-paths INT           align to (up to) this many alternate paths in between MEMs or in snarls [4]" << endl
    << "  -n, --unstranded              use lazy strand consistency when clustering MEMs" << endl
    << "  -b, --frag-sample INT         look for this many unambiguous mappings to estimate the fragment length distribution [1000]" << endl
    << "  -I, --frag-mean               mean for fixed fragment length distribution" << endl
    << "  -D, --frag-stddev             standard deviation for fixed fragment length distribution" << endl
    << "  -B, --no-calibrate            do not auto-calibrate mismapping dectection" << endl
    << "  -P, --max-p-val FLOAT         background model p value must be less than this to avoid mismapping detection [0.00001]" << endl
    << "  -Q, --mq-max INT              cap mapping quality estimates at this much [60]" << endl
    << "  -p, --padding-mult FLOAT      pad dynamic programming bands in inter-MEM alignment FLOAT * sqrt(read length) [1.0]" << endl
    << "  -u, --map-attempts INT        perform (up to) this many mappings per read (0 for no limit) [24 paired / 64 unpaired]" << endl
    << "  -O, --max-paths INT           consider (up to) this many paths per alignment for population consistency scoring, 0 to disable [10]" << endl
    << "  -M, --max-multimaps INT       report (up to) this many mappings per read [1]" << endl
    << "  -r, --reseed-length INT       reseed SMEMs for internal MEMs if they are at least this long (0 for no reseeding) [28]" << endl
    << "  -W, --reseed-diff FLOAT       require internal MEMs to have length within this much of the SMEM's length [0.45]" << endl
    << "  -K, --clust-length INT        minimum MEM length form clusters [automatic]" << endl
    << "  -c, --hit-max INT             use at most this many hits for any MEM (0 for no limit) [1024]" << endl
    << "  -w, --approx-exp FLOAT        let the approximate likelihood miscalculate likelihood ratios by this power [10.0]" << endl
    << "  --recombination-penalty FLOAT use this log recombination penalty for GBWT haplotype scoring [20.7]" << endl
    << "  --always-check-population     always try to population-score reads, even if there is only a single mapping" << endl
    << "  --delay-population            do not apply population scoring at intermediate stages of the mapping algorithm" << endl
    << "  --force-haplotype-count INT   assume that INT haplotypes ought to run through each fixed part of the graph, if nonzero [0]" << endl
    << "  -C, --drop-subgraph FLOAT     drop alignment subgraphs whose MEMs cover this fraction less of the read than the best subgraph [0.2]" << endl
    << "  -U, --prune-exp FLOAT         prune MEM anchors if their approximate likelihood is this root less than the optimal anchors [1.25]" << endl
    << "scoring:" << endl
    << "  -q, --match INT               use this match score [1]" << endl
    << "  -z, --mismatch INT            use this mismatch penalty [4]" << endl
    << "  --score-matrix FILE           read a 5x5 integer substitution scoring matrix from a file" << endl
    << "  -o, --gap-open INT            use this gap open penalty [6]" << endl
    << "  -y, --gap-extend INT          use this gap extension penalty [1]" << endl
    << "  -L, --full-l-bonus INT        add this score to alignments that use the full length of the read [5]" << endl
    << "  -m, --remove-bonuses          remove full length alignment bonuses in reported scores" << endl
    << "computational parameters:" << endl
    << "  -t, --threads INT             number of compute threads to use" << endl
    << "  -Z, --buffer-size INT         buffer this many alignments together (per compute thread) before outputting to stdout [100]" << endl;
    
}

int main_mpmap(int argc, char** argv) {

    if (argc == 2) {
        help_mpmap(argv);
        return 1;
    }

    // initialize parameters with their default options
    #define OPT_SCORE_MATRIX 1000
    #define OPT_RECOMBINATION_PENALTY 1001
    #define OPT_ALWAYS_CHECK_POPULATION 1002
    #define OPT_DELAY_POPULATION_SCORING 1003
    #define OPT_FORCE_HAPLOTYPE_COUNT 1004
    string matrix_file_name;
    string xg_name;
    string gcsa_name;
    string gbwt_name;
    string sublinearLS_name;
    string sublinearLS_ref_path;
    string snarls_name;
    string distance_index_name;
    string fastq_name_1;
    string fastq_name_2;
    string gam_file_name;
    int match_score = default_match;
    int mismatch_score = default_mismatch;
    int gap_open_score = default_gap_open;
    int gap_extension_score = default_gap_extension;
    int full_length_bonus = default_full_length_bonus;
    bool interleaved_input = false;
    int snarl_cut_size = 5;
    int max_paired_end_map_attempts = 24;
    int max_single_end_mappings_for_rescue = 64;
    int max_single_end_map_attempts = 64;
    int max_rescue_attempts = 10;
    int population_max_paths = 10;
    // How many distinct single path alignments should we look for in a multipath, for MAPQ?
    // TODO: create an option.
    int localization_max_paths = 5;
    int max_num_mappings = 1;
    int buffer_size = 100;
    int hit_max = 1024;
    int min_mem_length = 1;
    int min_clustering_mem_length = 0;
    int reseed_length = 28;
    double reseed_diff = 0.45;
    double reseed_exp = 0.065;
    bool use_adaptive_reseed = true;
    double cluster_ratio = 0.2;
    bool use_tvs_clusterer = false;
    bool qual_adjusted = true;
    bool strip_full_length_bonus = false;
    MappingQualityMethod mapq_method = Adaptive;
    double band_padding_multiplier = 1.0;
    int max_dist_error = 12;
    int num_alt_alns = 4;
    double suboptimal_path_exponent = 1.25;
    double likelihood_approx_exp = 10.0;
    double recombination_penalty = 20.7;
    bool always_check_population = false;
    bool delay_population_scoring = false;
    size_t force_haplotype_count = 0;
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
    size_t secondary_rescue_attempts = 4;
    size_t rescue_only_min = numeric_limits<size_t>::max(); // disabling this for now
    size_t rescue_only_anchor_max = 16;
    string sample_name = "";
    string read_group = "";
    bool prefilter_redundant_hits = true;
    bool precollapse_order_length_hits = true;
    int max_sub_mem_recursion_depth = 1;
    int max_map_attempts_arg = 0;
    int secondary_rescue_subopt_diff = 35;
    int min_median_mem_coverage_for_split = 0;
    bool suppress_cluster_merging = false;
    bool dynamic_max_alt_alns = true;
    bool simplify_topologies = true;
    bool long_read_scoring = false;
    int match_score_arg = std::numeric_limits<int>::min();
    int mismatch_score_arg = std::numeric_limits<int>::min();
    int gap_open_score_arg = std::numeric_limits<int>::min();
    int gap_extension_score_arg = std::numeric_limits<int>::min();
    int full_length_bonus_arg = std::numeric_limits<int>::min();
    
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gcsa-name", required_argument, 0, 'g'},
            {"gbwt-name", required_argument, 0, 'H'},
            {"dist-name", required_argument, 0, 'd'},
            {"linear-index", required_argument, 0, 1},
            {"linear-path", required_argument, 0, 2},
            {"fastq", required_argument, 0, 'f'},
            {"gam-input", required_argument, 0, 'G'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"interleaved", no_argument, 0, 'i'},
            {"same-strand", no_argument, 0, 'e'},
            {"single-path-mode", no_argument, 0, 'S'},
            {"snarls", required_argument, 0, 's'},
            {"tvs-clusterer", no_argument, 0, 'v'},
            {"snarl-max-cut", required_argument, 0, 'X'},
            {"alt-paths", required_argument, 0, 'a'},
            {"unstranded", no_argument, 0, 'n'},
            {"frag-sample", required_argument, 0, 'b'},
            {"frag-mean", required_argument, 0, 'I'},
            {"frag-stddev", required_argument, 0, 'D'},
            {"no-calibrate", no_argument, 0, 'B'},
            {"max-p-val", required_argument, 0, 'P'},
            {"mq-max", required_argument, 0, 'Q'},
            {"padding-mult", required_argument, 0, 'p'},
            {"map-attempts", required_argument, 0, 'u'},
            {"max-paths", required_argument, 0, 'O'},
            {"max-multimaps", required_argument, 0, 'M'},
            {"reseed-length", required_argument, 0, 'r'},
            {"reseed-diff", required_argument, 0, 'W'},
            {"clustlength", required_argument, 0, 'K'},
            {"hit-max", required_argument, 0, 'c'},
            {"approx-exp", required_argument, 0, 'w'},
            {"recombination-penalty", required_argument, 0, OPT_RECOMBINATION_PENALTY},
            {"always-check-population", no_argument, 0, OPT_ALWAYS_CHECK_POPULATION},
            {"delay-population", no_argument, 0, OPT_DELAY_POPULATION_SCORING},
            {"force-haplotype-count", required_argument, 0, OPT_FORCE_HAPLOTYPE_COUNT},
            {"drop-subgraph", required_argument, 0, 'C'},
            {"prune-exp", required_argument, 0, 'U'},
            {"long-read-scoring", no_argument, 0, 'E'},
            {"match", required_argument, 0, 'q'},
            {"mismatch", required_argument, 0, 'z'},
            {"score-matrix", required_argument, 0, OPT_SCORE_MATRIX},
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
        c = getopt_long (argc, argv, "hx:g:H:d:f:G:N:R:ieSs:vX:u:O:a:nb:I:D:BP:Q:p:M:r:W:K:c:w:C:R:Eq:z:o:y:L:mAt:Z:",
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
                if (gbwt_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide GBWT index file with -H" << endl;
                    exit(1);
                }
                break;
                
            case 'd':
                distance_index_name = optarg;
                if (distance_index_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide distance index file with -d" << endl;
                    exit(1);
                }
                break;
                
            case 1: // --linear-index
                sublinearLS_name = optarg;
                break;
            
            case 2: // --linear-path
                sublinearLS_ref_path = optarg;
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
                    cerr << "error:[vg mpmap] Must provide sample name with -N." << endl;
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
                
            case 'v':
                use_tvs_clusterer = true;
                break;
                
            case 'X':
                snarl_cut_size = parse<int>(optarg);
                break;
                
            case 'a':
                num_alt_alns = parse<int>(optarg);
                break;
                
            case 'n':
                unstranded_clustering = true;
                break;
                
            case 'b':
                frag_length_sample_size = parse<int>(optarg);
                break;
                
            case 'I':
                frag_length_mean = parse<double>(optarg);
                break;
                
            case 'D':
                frag_length_stddev = parse<double>(optarg);
                break;
                
            case 'B':
                auto_calibrate_mismapping_detection = false;
                break;
                
            case 'P':
                max_mapping_p_value = parse<double>(optarg);
                break;
                
            case 'Q':
                max_mapq = parse<int>(optarg);
                break;
                
            case 'p':
                band_padding_multiplier = parse<double>(optarg);
                break;
                
            case 'u':
                max_map_attempts_arg = parse<int>(optarg);
                // let 0 be a sentinel for no limit and also a sentinel for not giving an arg
                if (max_map_attempts_arg == 0) {
                    max_map_attempts_arg = numeric_limits<int>::max();
                }
                break;
                
            case 'O':
                population_max_paths = parse<int>(optarg);
                break;
                
            case 'M':
                max_num_mappings = parse<int>(optarg);
                break;
                
            case 'r':
                reseed_length = parse<int>(optarg);
                break;
                
            case 'W':
                reseed_diff = parse<double>(optarg);
                break;
                
            case 'K':
                min_clustering_mem_length = parse<int>(optarg);
                break;
                
            case 'c':
                hit_max = parse<int>(optarg);
                break;
                
            case 'w':
                likelihood_approx_exp = parse<double>(optarg);
                break;
                
            case OPT_RECOMBINATION_PENALTY:
                recombination_penalty = parse<double>(optarg);
                break;
                
            case OPT_ALWAYS_CHECK_POPULATION:
                always_check_population = true;
                break;
                
            case OPT_DELAY_POPULATION_SCORING:
                delay_population_scoring = true;
                break;
                
            case OPT_FORCE_HAPLOTYPE_COUNT:
                force_haplotype_count = parse<size_t>(optarg);
                break;
                
            case 'C':
                cluster_ratio = parse<double>(optarg);
                break;
                
            case 'U':
                suboptimal_path_exponent = parse<double>(optarg);
                break;
                
            case 'E':
                long_read_scoring = true;
                break;
                
            case 'q':
                match_score_arg = parse<int>(optarg);
                break;
                
            case 'z':
                mismatch_score_arg = parse<int>(optarg);
                break;
                
            case OPT_SCORE_MATRIX:
                matrix_file_name = optarg;
                if (matrix_file_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide matrix file with --matrix-file." << endl;
                    exit(1);
                }
                break;
                
            case 'o':
                gap_open_score_arg = parse<int>(optarg);
                break;
                
            case 'y':
                gap_extension_score_arg = parse<int>(optarg);
                break;
                
            case 'm':
                strip_full_length_bonus = true;
                break;
                
            case 'L':
                full_length_bonus_arg = parse<int>(optarg);
                break;
                
            case 'A':
                qual_adjusted = false;
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg mpmap] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'Z':
                buffer_size = parse<int>(optarg);
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
    
    if (!distance_index_name.empty() && snarls_name.empty()) {
        cerr << "error:[vg mpmap] Snarl distance index (-d) requires a matching snarl file (-s) to also be provided." << endl;
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
        cerr << "warning:[vg mpmap] Ignoring same strand parameter (-e) because no paired end input provided." << endl;
    }
    
    if (num_alt_alns <= 0) {
        cerr << "error:[vg mpmap] Number of alternate snarl paths (-a) set to " << num_alt_alns << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (unstranded_clustering && use_tvs_clusterer) {
        cerr << "warning:[vg mpmap] Target value search clustering (-v) does not have an unstranded option (-n), ignoring unstranded option" << endl;
        unstranded_clustering = false;
    }
    else if (unstranded_clustering && !distance_index_name.empty()) {
        cerr << "warning:[vg mpmap] Snarl distance index-based clustering (-d) does not have an unstranded option (-n), ignoring unstranded option" << endl;
        unstranded_clustering = false;
    }
    
    if (frag_length_sample_size <= 0) {
        cerr << "error:[vg mpmap] Fragment length distribution sample size (-b) set to " << frag_length_sample_size << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (snarl_cut_size < 0) {
        cerr << "error:[vg mpmap] Max snarl cut size (-U) set to " << snarl_cut_size << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (max_mapping_p_value <= 0.0) {
        cerr << "error:[vg mpmap] Max mapping p-value (-P) set to " << max_mapping_p_value << ", must set to a positive number." << endl;
        exit(1);
    }
    
    if (max_mapq <= 0 && mapq_method != None) {
        cerr << "error:[vg mpmap] Maximum mapping quality (-Q) set to " << max_mapq << ", must set to a positive integer." << endl;
        exit(1);
    }
    
    if (band_padding_multiplier < 0.0) {
        cerr << "error:[vg mpmap] Band padding (-p) set to " << band_padding_multiplier << ", must set to a nonnegative number." << endl;
        exit(1);
    }
    
    if (max_map_attempts_arg < 0) {
        cerr << "error:[vg mpmap] Maximum number of mapping attempts (-u) set to " << max_map_attempts_arg << ", must set to a positive integer or 0 for no maximum." << endl;
        exit(1);
    }
    
    if (population_max_paths < 0) {
        cerr << "error:[vg mpmap] Maximum number of paths per alignment for population scoring (-O) set to " << population_max_paths << ", must set to a nonnegative integer." << endl;
        exit(1);
    }
    
    if (population_max_paths != 10 && population_max_paths != 0 && gbwt_name.empty() && sublinearLS_name.empty()) {
        // Don't allow anything but the default or the "disabled" setting without an index.
        // TODO: This restriction makes neat auto-generation of command line options for different conditions hard.
        cerr << "error:[vg mpmap] Maximum number of paths per alignment for population scoring (-O) is specified but population database (-H or --linear-index) was not provided." << endl;
        exit(1);
    }
    
    if (always_check_population && gbwt_name.empty() && sublinearLS_name.empty()) {
        cerr << "error:[vg mpmap] Cannot --always-check-population if no population database (-H or --linear-index) is provided." << endl;
        exit(1);
    }
    
    if (delay_population_scoring && gbwt_name.empty() && sublinearLS_name.empty()) {
        cerr << "warning:[vg mpmap] Cannot --delay-population scoring if no population database (-H or --linear-index) is provided. Ignoring option." << endl;
    }
    
    if (force_haplotype_count != 0 && gbwt_name.empty() && sublinearLS_name.empty()) {
        cerr << "warning:[vg mpmap] Cannot --force-haplotype-count if no population database (-H or --linear-index) is provided. Ignoring option." << endl;
    }
    
    if (!sublinearLS_name.empty() && !gbwt_name.empty()) {
        cerr << "error:[vg mpmap] GBWT index (-H) and linear haplotype index (--linear-index) both specified. Only one can be used." << endl;
        exit(1);
    }
    
    if (!sublinearLS_name.empty() && sublinearLS_ref_path.empty()) {
        cerr << "error:[vg mpmap] Linear haplotype index (--linear-index) cannot be used without a single reference path (--linear-path)." << endl;
        exit(1);
    }
    
    if (sublinearLS_name.empty() && !sublinearLS_ref_path.empty()) {
        cerr << "error:[vg mpmap] Linear haplotype ref path (--linear-path) cannot be used without an index (--linear-index)." << endl;
        exit(1);
    }
    
    // choose either the user supplied max or the default for paired/unpaired
    int max_map_attempts = max_map_attempts_arg ? max_map_attempts_arg : ((interleaved_input || !fastq_name_2.empty()) ?
                                                                          max_paired_end_map_attempts : max_single_end_map_attempts);
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
    
    if (likelihood_approx_exp < 1.0) {
        cerr << "error:[vg mpmap] Likelihood approximation exponent (-w) set to " << likelihood_approx_exp << ", must set to at least 1.0." << endl;
        exit(1);
    }
    
    if (cluster_ratio < 0.0 || cluster_ratio >= 1.0) {
        cerr << "error:[vg mpmap] Cluster drop ratio (-C) set to " << cluster_ratio << ", must set to a number between 0.0 and 1.0." << endl;
        exit(1);
    }
    
    if (use_tvs_clusterer && distance_index_name.empty()) {
        cerr << "error:[vg mpmap] The Target Value Search clusterer (-v) requires a distance index (-d)." << endl;
        exit(1);
    }
    
    if (suboptimal_path_exponent < 1.0) {
        cerr << "error:[vg mpmap] Suboptimal path likelihood root (-R) set to " << suboptimal_path_exponent << ", must set to at least 1.0." << endl;
        exit(1);
    }
    
    if ((match_score_arg != std::numeric_limits<int>::min() || mismatch_score_arg != std::numeric_limits<int>::min()) && !matrix_file_name.empty())  {
        cerr << "error:[vg mpmap] Cannot choose custom scoring matrix (--score-matrix) and custom match/mismatch score (-q/-z) simultaneously." << endl;
        exit(1);
    }
    
    if (long_read_scoring) {
        // defaults for long read scoring
        match_score = 1;
        mismatch_score = 1;
        gap_open_score = 1;
        gap_extension_score = 1;
        full_length_bonus = 0;
    }
    
    // if we indicated any other scores, apply those, possibly overriding
    if (match_score_arg != std::numeric_limits<int>::min()) {
        match_score = match_score_arg;
    }
    if (mismatch_score_arg != std::numeric_limits<int>::min()) {
        mismatch_score = mismatch_score_arg;
    }
    if (gap_open_score_arg != std::numeric_limits<int>::min()) {
        gap_open_score = gap_open_score_arg;
    }
    if (gap_extension_score_arg != std::numeric_limits<int>::min()) {
        gap_extension_score = gap_extension_score_arg;
    }
    if (full_length_bonus_arg != std::numeric_limits<int>::min()) {
        full_length_bonus = full_length_bonus_arg;
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
    
    // adjust parameters that produce irrelevant extra work or bad behavior single path mode
    
    if (single_path_alignment_mode && population_max_paths == 0) {
        // TODO: I don't like having these constants floating around in two different places, but it's not very risky, just a warning
        if (!snarls_name.empty()) {
            cerr << "warning:[vg mpmap] Snarl file (-s) is ignored in single path mode (-S) without multipath population scoring (-O)." << endl;
            // TODO: Not true!
        }
        
        if (snarl_cut_size != 5) {
            cerr << "warning:[vg mpmap] Snarl cut limit (-u) is ignored in single path mode (-S) without multipath population scoring (-O)." << endl;
        }
        
        if (num_alt_alns != 4) {
            cerr << "warning:[vg mpmap] Number of alternate alignments (-a) is ignored in single path mode (-S) without multipath population scoring (-O)." << endl;
        }
        
        num_alt_alns = 1;
    }
    
    if (single_path_alignment_mode && !long_read_scoring) {
        // we get better performance by splitting up clusters a bit more when we're forcing alignments to go to only one place
        min_median_mem_coverage_for_split = 2;
        suppress_cluster_merging = true;
    }
    
    if (single_path_alignment_mode) {
        // simplifying topologies is redundant work if we're just going to take the maximum weight path anyway
        simplify_topologies = false;
    }
    
    // ensure required parameters are provided
    
    if (xg_name.empty()) {
        cerr << "error:[vg mpmap] Multipath mapping requires an XG index, must provide XG file (-x)" << endl;
        exit(1);
    }
    
    if (gcsa_name.empty()) {
        cerr << "error:[vg mpmap] Multipath mapping requires a GCSA2 index, must provide GCSA2 file (-g)" << endl;
        exit(1);
    }
    
    // create in-memory objects
    
    ifstream xg_stream(xg_name);
    if (!xg_stream) {
        cerr << "error:[vg mpmap] Cannot open XG file " << xg_name << endl;
        exit(1);
    }
    
    ifstream gcsa_stream(gcsa_name);
    if (!gcsa_stream) {
        cerr << "error:[vg mpmap] Cannot open GCSA2 file " << gcsa_name << endl;
        exit(1);
    }
    
    string lcp_name = gcsa_name + ".lcp";
    ifstream lcp_stream(lcp_name);
    if (!lcp_stream) {
        cerr << "error:[vg mpmap] Cannot open LCP file " << lcp_name << endl;
        exit(1);
    }

    ifstream matrix_stream;
    if (!matrix_file_name.empty()) {
      matrix_stream.open(matrix_file_name);
      if (!matrix_stream) {
          cerr << "error:[vg mpmap] Cannot open scoring matrix file " << matrix_file_name << endl;
          exit(1);
      }
    }
    
    ifstream distance_index_stream;
    if (!distance_index_name.empty()) {
        distance_index_stream.open(distance_index_name);
        if (!distance_index_stream) {
            cerr << "error:[vg mpmap] Cannot open distance index file " << distance_index_name << endl;
            exit(1);
        }
    }
    
    ifstream snarl_stream;
    if (!snarls_name.empty()) {
        snarl_stream.open(snarls_name);
        if (!snarl_stream) {
            cerr << "error:[vg mpmap] Cannot open Snarls file " << snarls_name << endl;
            exit(1);
        }
    }
    
    ifstream gbwt_stream;
    ifstream ls_stream;
    if (!gbwt_name.empty()) {
        gbwt_stream.open(gbwt_name);
        if (!gbwt_stream) {
            cerr << "error:[vg mpmap] Cannot open GBWT file " << gbwt_name << endl;
            exit(1);
        }
    }
    else if (!sublinearLS_name.empty()) {
        // We want to use sublinear Li and Stephens as our haplotype scoring approach
        ls_stream.open(sublinearLS_name);
        if (!ls_stream) {
            cerr << "error:[vg mpmap] Cannot open sublinear Li & Stevens file " << sublinearLS_name << endl;
            exit(1);
        }
    }
    
    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Configure its temp directory to the system temp directory
    gcsa::TempFile::setDirectory(temp_file::get_dir());
    
    // Load required indexes
    
    xg::XG xg_index(xg_stream);
    gcsa::GCSA gcsa_index;
    gcsa_index.load(gcsa_stream);
    gcsa::LCPArray lcp_array;
    lcp_array.load(lcp_stream);
    
    // Load optional indexes
    
    gbwt::GBWT* gbwt = nullptr;
    haplo::linear_haplo_structure* sublinearLS = nullptr;
    haplo::ScoreProvider* haplo_score_provider = nullptr;
    if (!gbwt_name.empty()) {
        gbwt = new gbwt::GBWT();
        gbwt->load(gbwt_stream);
        
        // We have the GBWT available for scoring haplotypes
        haplo_score_provider = new haplo::GBWTScoreProvider<gbwt::GBWT>(*gbwt);
    } else if (!sublinearLS_name.empty()) {
        
        // TODO: we only support a single ref contig, and we use these
        // hardcoded mutation and recombination likelihoods
        
        // What is the rank of our one and only reference path
        auto xg_ref_rank = xg_index.path_rank(sublinearLS_ref_path);
        
        sublinearLS = new linear_haplo_structure(ls_stream, -9 * 2.3, -6 * 2.3, xg_index, xg_ref_rank);
        haplo_score_provider = new haplo::LinearScoreProvider(*sublinearLS);
    }
    // TODO: Allow using haplo::XGScoreProvider?
    
    SnarlManager* snarl_manager = nullptr;
    if (!snarls_name.empty()) {
        snarl_manager = new SnarlManager(snarl_stream);
    }
    
    DistanceIndex* distance_index = nullptr;
    if (!distance_index_name.empty()) {
        distance_index = new DistanceIndex(&xg_index, snarl_manager, distance_index_stream);
    }
    
    MultipathMapper multipath_mapper(&xg_index, &gcsa_index, &lcp_array, haplo_score_provider, snarl_manager, distance_index);
    
    // set alignment parameters
    multipath_mapper.set_alignment_scores(match_score, mismatch_score, gap_open_score, gap_extension_score, full_length_bonus);
    if(matrix_stream.is_open()) multipath_mapper.load_scoring_matrix(matrix_stream);
    multipath_mapper.adjust_alignments_for_base_quality = qual_adjusted;
    multipath_mapper.strip_bonuses = strip_full_length_bonus;
    multipath_mapper.band_padding_multiplier = band_padding_multiplier;
    multipath_mapper.init_band_padding_memo();
    
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
    // Use population MAPQs when we have the right option combination to make that sensible.
    multipath_mapper.use_population_mapqs = (haplo_score_provider != nullptr && population_max_paths > 0);
    multipath_mapper.population_max_paths = population_max_paths;
    multipath_mapper.recombination_penalty = recombination_penalty;
    multipath_mapper.always_check_population = always_check_population;
    multipath_mapper.delay_population_scoring = delay_population_scoring;
    multipath_mapper.force_haplotype_count = force_haplotype_count;
    
    // set pruning and clustering parameters
    multipath_mapper.use_tvs_clusterer = use_tvs_clusterer;
    multipath_mapper.max_expected_dist_approx_error = max_dist_error;
    multipath_mapper.mem_coverage_min_ratio = cluster_ratio;
    multipath_mapper.log_likelihood_approx_factor = likelihood_approx_exp;
    multipath_mapper.num_mapping_attempts = max_map_attempts;
    multipath_mapper.unstranded_clustering = unstranded_clustering;
    multipath_mapper.min_median_mem_coverage_for_split = min_median_mem_coverage_for_split;
    multipath_mapper.suppress_cluster_merging = suppress_cluster_merging;
    multipath_mapper.use_tvs_clusterer = use_tvs_clusterer;
    
    // set pair rescue parameters
    multipath_mapper.max_rescue_attempts = max_rescue_attempts;
    multipath_mapper.max_single_end_mappings_for_rescue = max(max(max_single_end_mappings_for_rescue, max_rescue_attempts), max_num_mappings);
    multipath_mapper.secondary_rescue_subopt_diff = secondary_rescue_subopt_diff;
    multipath_mapper.secondary_rescue_score_diff = secondary_rescue_score_diff;
    multipath_mapper.secondary_rescue_attempts = secondary_rescue_attempts;
    multipath_mapper.rescue_only_min = rescue_only_min;
    multipath_mapper.rescue_only_anchor_max = rescue_only_anchor_max;
    
    // set multipath alignment topology parameters
    multipath_mapper.max_snarl_cut_size = snarl_cut_size;
    multipath_mapper.num_alt_alns = num_alt_alns;
    multipath_mapper.dynamic_max_alt_alns = dynamic_max_alt_alns;
    multipath_mapper.simplify_topologies = simplify_topologies;
    multipath_mapper.max_suboptimal_path_score_ratio = suboptimal_path_exponent;
    
    // if directed to, auto calibrate the mismapping detection to the graph
    if (auto_calibrate_mismapping_detection) {
        multipath_mapper.calibrate_mismapping_detection(num_calibration_simulations, calibration_read_length);
    }
    
    // set computational paramters
    int thread_count = get_thread_count();
    multipath_mapper.set_alignment_threads(thread_count);
    
    // Establish a watchdog to find reads that take too long to map.
    // If we see any, we will issue a warning.
    unique_ptr<Watchdog> watchdog(new Watchdog(thread_count, chrono::minutes(20)));
    
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
            // For each multipath alignment, get the greedy nonoverlapping
            // single-path alignments from the top k optimal single-path
            // alignments.
            vector<Alignment> options;
            multipath_mapper.reduce_to_single_path(mp_aln, options, localization_max_paths);
            
            // There will always be at least one result. Use the optimal alignment.
            output_buf.emplace_back(std::move(options.front()));
            
            // compute the Alignment identity to make vg call happy
            output_buf.back().set_identity(identity(output_buf.back().path()));
            
            if (mp_aln.has_annotation()) {
                // Move over annotations
                output_buf.back().set_allocated_annotation(mp_aln.release_annotation());
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
            
            // Compute nonoverlapping single path alignments for each multipath alignment
            vector<Alignment> options;
            multipath_mapper.reduce_to_single_path(mp_aln_pair.first, options, localization_max_paths);
            
            // There will always be at least one result. Use the optimal alignment.
            output_buf.emplace_back(std::move(options.front()));
            
            if (mp_aln_pair.first.has_annotation()) {
                // Move over annotations
                output_buf.back().set_allocated_annotation(mp_aln_pair.first.release_annotation());
            }
            
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
            
            // Now do the second read
            options.clear();
            multipath_mapper.reduce_to_single_path(mp_aln_pair.second, options, localization_max_paths);
            output_buf.emplace_back(std::move(options.front()));
            
            if (mp_aln_pair.second.has_annotation()) {
                // Move over annotations
                output_buf.back().set_allocated_annotation(mp_aln_pair.second.release_annotation());
            }
            
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

        auto thread_num = omp_get_thread_num();

        if (watchdog) {
            watchdog->check_in(thread_num, alignment.name());
        }

        vector<MultipathAlignment> mp_alns;
        multipath_mapper.multipath_map(alignment, mp_alns, max_num_mappings);
        if (single_path_alignment_mode) {
            output_single_path_alignments(mp_alns);
        }
        else {
            output_multipath_alignments(mp_alns);
        }
        
        if (watchdog) {
            watchdog->check_out(thread_num);
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

        auto thread_num = omp_get_thread_num();

#ifdef record_read_run_times
        clock_t start = clock();
#endif

        if (watchdog) {
            watchdog->check_in(thread_num, alignment_1.name());
        }
        
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
        
        if (watchdog) {
            watchdog->check_out(thread_num);
        }
        
#ifdef record_read_run_times
        clock_t finish = clock();
#pragma omp critical
        read_time_file << alignment_1.name() << "\t" << alignment_2.name() << "\t" << double(finish - start) / CLOCKS_PER_SEC << endl;
#endif
    };
    
    // do unpaired, independent multipath alignment, and write to buffer as paired
    function<void(Alignment&, Alignment&)> do_independent_paired_alignments = [&](Alignment& alignment_1, Alignment& alignment_2) {
        // get reads on the same strand so that oriented distance estimation works correctly
        // but if we're clearing the ambiguous buffer we already RC'd these on the first pass

        auto thread_num = omp_get_thread_num();

#ifdef record_read_run_times
        clock_t start = clock();
#endif

        if (watchdog) {
            watchdog->check_in(thread_num, alignment_1.name());
        }

        if (!same_strand) {
            // TODO: the output functions undo this transformation, so we have to do it here.
        
            // remove the path so we won't try to RC it (the path may not refer to this graph)
            alignment_2.clear_path();
            reverse_complement_alignment_in_place(&alignment_2, [&](vg::id_t node_id) { return xg_index.node_length(node_id); });
        }
        
        // Align independently
        vector<MultipathAlignment> mp_alns_1, mp_alns_2;
        multipath_mapper.multipath_map(alignment_1, mp_alns_1, max_num_mappings);
        multipath_mapper.multipath_map(alignment_2, mp_alns_2, max_num_mappings);
               
        vector<pair<MultipathAlignment, MultipathAlignment>> mp_aln_pairs;
        for (size_t i = 0; i < mp_alns_1.size() && i < mp_alns_2.size(); i++) {
            // Pair arbitrarily. Stop when one side runs out of alignments.
            mp_aln_pairs.emplace_back(mp_alns_1[i], mp_alns_2[i]);
        }
        
        // TODO: Set a flag or annotation or something to say we don't really believe the pairing
       
        if (single_path_alignment_mode) {
            output_single_path_paired_alignments(mp_aln_pairs);
        }
        else {
            output_multipath_paired_alignments(mp_aln_pairs);
        }
        
        if (watchdog) {
            watchdog->check_out(thread_num);
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
                do_independent_paired_alignments(aln_pair.first, aln_pair.second);
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
    //cerr << "attempted to split " << OrientedDistanceClusterer::SPLIT_ATTEMPT_COUNTER << " of " << OrientedDistanceClusterer::PRE_SPLIT_CLUSTER_COUNTER << " clusters with " << OrientedDistanceClusterer::SUCCESSFUL_SPLIT_ATTEMPT_COUNTER << " splits successful (" << 100.0 * double(OrientedDistanceClusterer::SUCCESSFUL_SPLIT_ATTEMPT_COUNTER) / OrientedDistanceClusterer::SPLIT_ATTEMPT_COUNTER << "%) resulting in " << OrientedDistanceClusterer::POST_SPLIT_CLUSTER_COUNTER << " total clusters (" << OrientedDistanceClusterer::POST_SPLIT_CLUSTER_COUNTER - OrientedDistanceClusterer::PRE_SPLIT_CLUSTER_COUNTER << " new)" << endl;
    //cerr << "entered secondary rescue " << MultipathMapper::SECONDARY_RESCUE_TOTAL << " times with " << MultipathMapper::SECONDARY_RESCUE_COUNT << " actually attempting rescues, totaling " << MultipathMapper::SECONDARY_RESCUE_ATTEMPT << " rescues (" << double(MultipathMapper::SECONDARY_RESCUE_ATTEMPT) / MultipathMapper::SECONDARY_RESCUE_COUNT << " average per attempt)" << endl;
    
    if (snarl_manager != nullptr) {
        delete snarl_manager;
    }
   
    if (haplo_score_provider != nullptr) {
        delete haplo_score_provider;
    }
    
    if (sublinearLS != nullptr) {
        delete sublinearLS;
    }
   
    if (gbwt != nullptr) {
        delete gbwt;
    }
    
    if (distance_index != nullptr) {
        delete distance_index;
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_mpmap("mpmap", "multipath alignments of reads to a graph", main_mpmap);



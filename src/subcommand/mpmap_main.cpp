/**
 * \file mpmap_main.cpp: multipath mapping of reads to a graph
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include "subcommand.hpp"

#include "xg.hpp"
#include "gcsa.h"
#include "../utility.hpp"
#include "../alignment.hpp"
#include "../multipath_mapper.hpp"
#include "../gssw_aligner.hpp"
#include "../mem.hpp"

//#define debug_mpmap

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_mpmap(char** argv) {
    cerr
    << "usage: " << argv[0] << " mpmap [options] -x index.xg -g index.gcsa [-f reads1.fq [-f reads2.fq] | -b reads.bam | -G reads.gam] > aln.gamp" << endl
    << "Multipath align reads to a graph." << endl
    << endl
    << "basic options:" << endl
    << "graph/index:" << endl
    << "  -x, --xg-name FILE        use this xg index (required)" << endl
    << "  -g, --gcsa-name FILE      use this GCSA2/LCP index pair (required; both FILE and FILE.lcp)" << endl
    << "input:" << endl
    << "  -f, --fastq FILE          input FASTQ (possibly compressed), can be given twice for paired ends (for stdin use -)" << endl
    << "  -i, --interleaved         FASTQ is interleaved with paired ends" << endl
    << "  -b, --hts-input FILE      align reads from this htslib-compatible file (BAM/CRAM/SAM; for stdin use -)" << endl
    << "  -G, --gam-input FILE      realign .gam input (for stdin, use -)" << endl
    << "algorithm:" << endl
    << "  -S, --single-path-mode    produce single-path alignments (.gam) instead of multipath alignments (.gamp) (ignores -sua)" << endl
    << "  -s, --snarls FILE         align to alternate paths in these snarls" << endl
    << "scoring:" << endl
    << "  -A, --no-qual-adjust      do not perform base quality adjusted alignments (required if input does not have base qualities)" << endl
    << endl
    << "advanced options:" << endl
    << "algorithm:" << endl
    << "  -U, --snarl-max-cut INT   do not align to alternate paths in a snarl if an exact match is at least this long (0 for no limit) [5]" << endl
    << "  -a, --alt-paths INT       align to (up to) this many alternate paths in between MEMs or in snarls [4]" << endl
    << "  -v, --mq-method OPT       mapping quality method: 0 - none, 1 - fast approximation, 2 - exact [1]" << endl
    << "  -Q, --mq-max OPT          cap mapping quality estimates at this much [60]" << endl
    << "  -p, --band-padding INT    pad dynamic programming bands in inter-MEM alignment by this much [2]" << endl
    << "  -u, --map-attempts INT    perform (up to) this many mappings per read (0 for no limit) [32]" << endl
    << "  -M, --max-multimaps INT   report (up to) this many mappings per read [1]" << endl
    << "  -r, --reseed-length INT   reseed SMEMs for internal MEMs if they are at least this long (0 for no reseeding) [32]" << endl
    << "  -W, --reseed-diff INT     require internal MEMs to have length within tåhis much of the SMEM's length [8]" << endl
    << "  -k, --min-mem-length INT  minimum MEM length to anchor multipath alignments [1]" << endl
    << "  -c, --hit-max INT         ignore MEMs that occur greater than this many times in the graph (0 for no limit) [128]" << endl
    << "  -d, --max-dist-error INT  maximum typical deviation between distance on a reference path and distance in graph [8]" << endl
    << "  -w, --approx-exp FLOAT    let the approximate likelihood miscalculate likelihood ratios by this power [4.0]" << endl
    << "  -C, --drop-subgraph FLOAT drop alignment subgraphs whose MEMs cover this fraction less of the read than the best subgraph [0.5]" << endl
    << "  -R, --prune-ratio FLOAT   prune MEM anchors if their approximate likelihood is this ratio less than the optimal anchors [10000.0]" << endl
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
    string snarls_name;
    string fastq_name_1;
    string fastq_name_2;
    string hts_file_name;
    string gam_file_name;
    int match_score = default_match;
    int mismatch_score = default_mismatch;
    int gap_open_score = default_gap_open;
    int gap_extension_score = default_gap_extension;
    int full_length_bonus = 5;
    bool interleaved_input = false;
    int snarl_cut_size = 5;
    int max_map_attempts = 32;
    int max_num_mappings = 1;
    int buffer_size = 100;
    int hit_max = 128;
    int min_mem_length = 1;
    int reseed_length = 32;
    int reseed_diff = 8;
    double cluster_ratio = 0.5;
    bool qual_adjusted = true;
    bool strip_full_length_bonus = false;
    MappingQualityMethod mapq_method = Approx;
    int band_padding = 2;
    int max_dist_error = 8;
    int num_alt_alns = 4;
    double suboptimal_path_ratio = 10000.0;
    double likelihood_approx_exp = 4.0;
    bool single_path_alignment_mode = false;
    int max_mapq = 60;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"gcsa-name", required_argument, 0, 'g'},
            {"hts-input", required_argument, 0, 'b'},
            {"fastq", required_argument, 0, 'f'},
            {"interleaved", no_argument, 0, 'i'},
            {"gam-input", required_argument, 0, 'G'},
            {"single-path-mode", no_argument, 0, 'S'},
            {"snarls", required_argument, 0, 's'},
            {"snarl-max-cut", required_argument, 0, 'U'},
            {"alt-paths", required_argument, 0, 'a'},
            {"mq-method", required_argument, 0, 'v'},
            {"mq-max", required_argument, 0, 'Q'},
            {"band-padding", required_argument, 0, 'p'},
            {"map-attempts", required_argument, 0, 'u'},
            {"max-multimaps", required_argument, 0, 'M'},
            {"reseed-length", required_argument, 0, 'r'},
            {"reseed-diff", required_argument, 0, 'W'},
            {"min-mem-length", required_argument, 0, 'k'},
            {"hit-max", required_argument, 0, 'c'},
            {"max-dist-error", required_argument, 0, 'd'},
            {"approx-exp", required_argument, 0, 'w'},
            {"drop-subgraph", required_argument, 0, 'C'},
            {"prune-ratio", required_argument, 0, 'R'},
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
        c = getopt_long (argc, argv, "hx:g:b:f:iG:Ss:u:a:v:Q:p:M:r:W:k:c:d:w:C:R:q:z:o:y:L:mAt:Z:",
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
                
            case 'i':
                interleaved_input = true;
                break;
                
            case 'b':
                hts_file_name = optarg;
                if (hts_file_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide HTS file (SAM/BAM/CRAM) with -b." << endl;
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
                
            case 'U':
                snarl_cut_size = atoi(optarg);
                break;
                
            case 'a':
                num_alt_alns = atoi(optarg);
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
                
            case 'M':
                max_num_mappings = atoi(optarg);
                break;
                
            case 'r':
                reseed_length = atoi(optarg);
                break;
                
            case 'W':
                reseed_diff = atoi(optarg);
                break;
                
            case 'k':
                min_mem_length = atoi(optarg);
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
                
            case 'R':
                suboptimal_path_ratio = atof(optarg);
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
    
    if (num_alt_alns <= 0) {
        cerr << "error:[vg mpmap] Number of alternate snarl paths (-a) set to " << num_alt_alns << ", must set to a positive integer." << endl;
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
    
    if (reseed_diff <= 0 && reseed_length != 0) {
        cerr << "error:[vg mpmap] Reseeding length difference (-W) set to " << reseed_diff << ", must set to a positive integer if reseeding." << endl;
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
    
    if (suboptimal_path_ratio < 1.0) {
        cerr << "error:[vg mpmap] Suboptimal path likelihood ratio (-R) set to " << suboptimal_path_ratio << ", must set to at least 1.0." << endl;
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
    
    if (single_path_alignment_mode) {
        // TODO: I don't like having these constants floating around in two different places, but it's not very risky, just a warning
        if (!snarls_name.empty()) {
            cerr << "warning:[vg mpmap] Snarl file (-s) is ignored in single path mode (-S)." << endl;
        }
        if (num_alt_alns != 4) {
            cerr << "warning:[vg mpmap] Number of alternate alignments (-a) is ignored in single path mode (-S)." << endl;
        }
        if (snarl_cut_size != 5) {
            cerr << "warning:[vg mpmap] Snarl cut limit (-u) is ignored in single path mode (-S)." << endl;
        }
        snarls_name = "";
        num_alt_alns = 1;
        snarl_cut_size = 0;
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
    
    SnarlManager* snarl_manager = nullptr;
    if (!snarls_name.empty()) {
        ifstream snarl_stream(snarls_name);
        if (!snarl_stream) {
            cerr << "error:[vg mpmap] Cannot open Snarls file " << snarls_name << endl;
            exit(1);
        }
        snarl_manager = new SnarlManager(snarl_stream);
    }
    
    MultipathMapper multipath_mapper(&xg_index, &gcsa_index, &lcp_array, snarl_manager);
    
    // set alignment parameters
    multipath_mapper.set_alignment_scores(match_score, mismatch_score, gap_open_score, gap_extension_score, full_length_bonus);
    multipath_mapper.adjust_alignments_for_base_quality = qual_adjusted;
    multipath_mapper.strip_bonuses = strip_full_length_bonus;
    multipath_mapper.band_padding = band_padding;
    
    // set mem finding parameters
    multipath_mapper.hit_max = hit_max;
    multipath_mapper.min_mem_length = min_mem_length;
    multipath_mapper.mem_reseed_length = reseed_length;
    multipath_mapper.fast_reseed = true;
    multipath_mapper.fast_reseed_length_diff = reseed_diff;
    
    // set other algorithm parameters
    multipath_mapper.mapping_quality_method = mapq_method;
    multipath_mapper.max_mapping_quality = max_mapq;
    multipath_mapper.mem_coverage_min_ratio = cluster_ratio;
    multipath_mapper.max_expected_dist_approx_error = max_dist_error;
    multipath_mapper.max_snarl_cut_size = snarl_cut_size;
    multipath_mapper.num_alt_alns = num_alt_alns;
    multipath_mapper.num_mapping_attempts = max_map_attempts ? max_map_attempts : numeric_limits<int>::max();
    multipath_mapper.log_likelihood_approx_factor = likelihood_approx_exp;
    multipath_mapper.set_suboptimal_path_likelihood_ratio(suboptimal_path_ratio); // note: do this after choosing whether qual adj alignments
    
    int thread_count = get_thread_count();
    multipath_mapper.set_alignment_threads(thread_count);
    
    vector<vector<Alignment> > single_path_output_buffer(thread_count);
    vector<vector<MultipathAlignment> > multipath_output_buffer(thread_count);
    
    auto output_multipath_alignments = [&](vector<MultipathAlignment>& mp_alns) {
        auto& output_buf = multipath_output_buffer[omp_get_thread_num()];
        
        // Copy all the alignments over to the output buffer
        copy(mp_alns.begin(), mp_alns.end(), back_inserter(output_buf));
        
        stream::write_buffered(cout, output_buf, buffer_size);
    };
    
    auto output_single_path_alignments = [&](vector<MultipathAlignment>& mp_alns) {
        auto& output_buf = single_path_output_buffer[omp_get_thread_num()];
        
        // Copy all the alignments over to the output buffer
        for (const auto& mp_aln : mp_alns) {
            output_buf.emplace_back();
            optimal_alignment(mp_aln, output_buf.back());
        }
        
        stream::write_buffered(cout, output_buf, buffer_size);
    };
    
#ifdef debug_mpmap
    cerr << "[vg mpmap] created all in memory objects, beginning mapping" << endl;
#endif
    
    function<void(Alignment&)> do_unpaired_alignments = [&](Alignment& alignment) {
        vector<MultipathAlignment> mp_alns;
        multipath_mapper.multipath_map(alignment, mp_alns, max_num_mappings);
        if (single_path_alignment_mode) {
            output_single_path_alignments(mp_alns);
        }
        else {
            output_multipath_alignments(mp_alns);
        }
    };
    
    if (!fastq_name_1.empty()) {
        if (fastq_name_2.empty()) {
            fastq_unpaired_for_each_parallel(fastq_name_1, do_unpaired_alignments);
        }
        else if (interleaved_input) {
            // TODO
        }
        else {
            // TODO
        }
    }
    
    if (!hts_file_name.empty()) {
        // TODO
    }
    
    if (!gam_file_name.empty()) {
        ifstream gam_in(gam_file_name);
        if (!gam_in) {
            cerr << "error:[vg mpmap] Cannot open GAM file " << gam_file_name << endl;
            exit(1);
        }
        if (interleaved_input) {
            // TODO
        }
        else {
            stream::for_each_parallel(gam_in, do_unpaired_alignments);
        }
        gam_in.close();
    }
    
    // clear output buffers
    for (int i = 0; i < thread_count; i++) {
        vector<Alignment>& single_path_buffer = single_path_output_buffer[i];
        stream::write_buffered(cout, single_path_buffer, 0);
        
        vector<MultipathAlignment>& multipath_buffer = multipath_output_buffer[i];
        stream::write_buffered(cout, multipath_buffer, 0);
    }
    cout.flush();
    
    //cerr << "MEM cluster filtering efficiency: " << ((double) OrientedDistanceClusterer::PRUNE_COUNTER) / OrientedDistanceClusterer::CLUSTER_TOTAL << " (" << OrientedDistanceClusterer::PRUNE_COUNTER << "/" << OrientedDistanceClusterer::CLUSTER_TOTAL << ")" << endl;
    //cerr << "subgraph filtering efficiency: " << ((double) MultipathMapper::PRUNE_COUNTER) / MultipathMapper::SUBGRAPH_TOTAL << " (" << MultipathMapper::PRUNE_COUNTER << "/" << MultipathMapper::SUBGRAPH_TOTAL << ")" << endl;
    
    delete snarl_manager;
    
    return 0;
}

// Register subcommand
static Subcommand vg_mpmap("mpmap", "multipath alignments of reads to a graph", main_mpmap);



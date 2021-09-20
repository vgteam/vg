#include "subcommand.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../utility.hpp"
#include "../mapper.hpp"
#include "../hts_alignment_emitter.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

#include <unistd.h>
#include <getopt.h>
#include <chrono>

using namespace vg;
using namespace vg::subcommand;

void help_map(char** argv) {
    cerr << "usage: " << argv[0] << " map [options] -d idxbase -f in1.fq [-f in2.fq] >aln.gam" << endl
         << "Align reads to a graph." << endl
         << endl
         << "graph/index:" << endl
         << "    -d, --base-name BASE          use BASE.xg and BASE.gcsa as the input index pair" << endl
         << "    -x, --xg-name FILE            use this xg index or graph (defaults to <graph>.vg.xg)" << endl
         << "    -g, --gcsa-name FILE          use this GCSA2 index (defaults to <graph>" << gcsa::GCSA::EXTENSION << ")" << endl
         << "    -1, --gbwt-name FILE          use this GBWT haplotype index (defaults to <graph>"<<gbwt::GBWT::EXTENSION << ")" << endl
         << "algorithm:" << endl
         << "    -t, --threads N               number of compute threads to use" << endl
         << "    -k, --min-mem INT             minimum MEM length (if 0 estimate via -e) [0]" << endl
         << "    -e, --mem-chance FLOAT        set {-k} such that this fraction of {-k} length hits will by chance [5e-4]" << endl
         << "    -c, --hit-max N               ignore MEMs who have >N hits in our index (0 for no limit) [2048]" << endl
         << "    -Y, --max-mem INT             ignore mems longer than this length (unset if 0) [0]" << endl
         << "    -r, --reseed-x FLOAT          look for internal seeds inside a seed longer than FLOAT*--min-seed [1.5]" << endl
         << "    -u, --try-up-to INT           attempt to align up to the INT best candidate chains of seeds (1/2 for paired) [128]" << endl
         << "    -l, --try-at-least INT        attempt to align at least the INT best candidate chains of seeds [1]" << endl
         << "    -E, --approx-mq-cap INT       weight MQ by suffix tree based estimate when estimate less than FLOAT [0]" << endl
         << "    --id-mq-weight N              scale mapping quality by the alignment score identity to this power [2]" << endl
         << "    -W, --min-chain INT           discard a chain if seeded bases shorter than INT [0]" << endl
         << "    -C, --drop-chain FLOAT        drop chains shorter than FLOAT fraction of the longest overlapping chain [0.45]" << endl
         << "    -n, --mq-overlap FLOAT        scale MQ by count of alignments with this overlap in the query with the primary [0]" << endl
         << "    -P, --min-ident FLOAT         accept alignment only if the alignment identity is >= FLOAT [0]" << endl
         << "    -H, --max-target-x N          skip cluster subgraphs with length > N*read_length [100]" << endl
         << "    -w, --band-width INT          band width for long read alignment [256]" << endl
         << "    -O, --band-overlap INT        band overlap for long read alignment [{-w}/8]" << endl
         << "    -J, --band-jump INT           the maximum number of bands of insertion we consider in the alignment chain model [128]" << endl
         << "    -B, --band-multi INT          consider this many alignments of each band in banded alignment [16]" << endl
         << "    -Z, --band-min-mq INT         treat bands with less than this MQ as unaligned [0]" << endl
         << "    -I, --fragment STR            fragment length distribution specification STR=m:μ:σ:o:d [5000:0:0:0:1]" << endl
         << "                                  max, mean, stdev, orientation (1=same, 0=flip), direction (1=forward, 0=backward)" << endl
         << "    -U, --fixed-frag-model        don't learn the pair fragment model online, use {-I} without update" << endl
         << "    -p, --print-frag-model        suppress alignment output and print the fragment model on stdout as per {-I} format" << endl
         << "    --frag-calc INT               update the fragment model every INT perfect pairs [10]" << endl
         << "    --fragment-x FLOAT            calculate max fragment size as frag_mean+frag_sd*FLOAT [10]" << endl
         << "    --mate-rescues INT            attempt up to INT mate rescues per pair [64]" << endl
         << "    -S, --unpaired-cost INT       penalty for an unpaired read pair [17]" << endl
         << "    --no-patch-aln                do not patch banded alignments by locally aligning unaligned regions" << endl
         << "    --xdrop-alignment             use X-drop heuristic (much faster for long-read alignment)" << endl
         << "    --max-gap-length              maximum gap length allowed in each contiguous alignment (for X-drop alignment) [40]" << endl
         << "scoring:" << endl
         << "    -q, --match INT               use this match score [1]" << endl
         << "    -z, --mismatch INT            use this mismatch penalty [4]" << endl
         << "    --score-matrix FILE           read a 4x4 integer substitution scoring matrix from a file" << endl
         << "    -o, --gap-open INT            use this gap open penalty [6]" << endl
         << "    -y, --gap-extend INT          use this gap extension penalty [1]" << endl
         << "    -L, --full-l-bonus INT        the full-length alignment bonus [5]" << endl
         << "    --drop-full-l-bonus           remove the full length bonus from the score before sorting and MQ calculation" << endl
         << "    -a, --hap-exp FLOAT           the exponent for haplotype consistency likelihood in alignment score [1]" << endl
         << "    --recombination-penalty FLOAT use this log recombination penalty for GBWT haplotype scoring [20.7]" << endl
         << "    -A, --qual-adjust             perform base quality adjusted alignments (requires base quality input)" << endl
         << "preset:" << endl
         << "    -m, --alignment-model STR     use a preset alignment scoring model, either \"short\" (default) or \"long\" (for ONT/PacBio)" << endl
         << "                                  \"long\" is equivalent to `-u 2 -L 63 -q 1 -z 2 -o 2 -y 1 -w 128 -O 32`" << endl
         << "input:" << endl
         << "    -s, --sequence STR            align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -V, --seq-name STR            name the sequence using this value (for graph modification with new named paths)" << endl
         << "    -T, --reads FILE              take reads (one per line) from FILE, write alignments to stdout" << endl
         << "    -b, --hts-input FILE          align reads from htslib-compatible FILE (BAM/CRAM/SAM) stdin (-), alignments to stdout" << endl
         << "    -G, --gam-input FILE          realign GAM input" << endl
         << "    -f, --fastq FILE              input fastq or (2-line format) fasta, possibly compressed, two are allowed, one for each mate" << endl
         << "    -F, --fasta FILE              align the sequences in a FASTA file that may have multiple lines per reference sequence" << endl
         << "    -i, --interleaved             fastq or GAM is interleaved paired-ended" << endl
         << "    -N, --sample NAME             for --reads input, add this sample" << endl
         << "    -R, --read-group NAME         for --reads input, add this read group" << endl
         << "output:" << endl
         << "    -j, --output-json             output JSON rather than an alignment stream (helpful for debugging)" << endl
         << "    -%, --gaf                     output alignments in GAF format" << endl
         << "    --surject-to TYPE             surject the output into the graph's paths, writing TYPE := bam |sam | cram" << endl
         << "    --ref-paths FILE              ordered list of paths in the graph, one per line or HTSlib .dict, for HTSLib @SQ headers" << endl
         << "    --buffer-size INT             buffer this many alignments together before outputting in GAM [512]" << endl
         << "    -X, --compare                 realign GAM input (-G), writing alignment with \"correct\" field set to overlap with input" << endl
         << "    -v, --refpos-table            for efficient testing output a table of name, chr, pos, mq, score" << endl
         << "    -K, --keep-secondary          produce alignments for secondary input alignments in addition to primary ones" << endl
         << "    -M, --max-multimaps INT       produce up to INT alignments for each read [1]" << endl
         << "    -Q, --mq-max INT              cap the mapping quality at INT [60]" << endl
         << "    --exclude-unaligned           exclude reads with no alignment" << endl
         << "    -D, --debug                   print debugging information about alignment to stderr" << endl
         << "    --log-time                    print runtime to stderr" << endl;

}

int main_map(int argc, char** argv) {

    std::chrono::time_point<std::chrono::system_clock> launch = std::chrono::system_clock::now();

    if (argc == 2) {
        help_map(argv);
        return 1;
    }

    #define OPT_SCORE_MATRIX 1000
    #define OPT_RECOMBINATION_PENALTY 1001
    #define OPT_EXCLUDE_UNALIGNED 1002
    #define OPT_REF_PATHS 1003
    string matrix_file_name;
    string seq;
    string qual;
    string seq_name;
    string db_name;
    string xg_name;
    string gcsa_name;
    string gbwt_name;
    string read_file;
    string hts_file;
    string fasta_file;
    bool keep_secondary = false;
    int hit_max = 2048;
    int max_multimaps = 1;
    int thread_count = 1;
    string output_format = "GAM";
    string ref_paths_name;
    bool exclude_unaligned = false;
    bool debug = false;
    float min_score = 0;
    string sample_name;
    string read_group;
    string fastq1, fastq2;
    bool interleaved_input = false;
    int band_width = 256;
    int band_overlap = -1;
    int band_multimaps = 16;
    int max_band_jump = 128;
    bool always_rescue = false;
    bool top_pairs_only = false;
    int max_mem_length = 0;
    int min_mem_length = 0;
    int min_cluster_length = 0;
    float mem_reseed_factor = 1.5;
    int max_target_factor = 100;
    int buffer_size = 512;
    int8_t match = default_match;
    int8_t mismatch = default_mismatch;
    int8_t gap_open = default_gap_open;
    int8_t gap_extend = default_gap_extension;
    int8_t full_length_bonus = default_full_length_bonus;
    int unpaired_penalty = 17;
    double haplotype_consistency_exponent = 1;
    double recombination_penalty = 20.7;
    bool strip_bonuses = false;
    bool qual_adjust_alignments = false;
    int extra_multimaps = 128;
    int min_multimaps = 1;
    int max_mapping_quality = 60;
    double maybe_mq_threshold = 0;
    double identity_weight = 2;
    string gam_input;
    bool compare_gam = false;
    int fragment_max = 5000;
    int fragment_size = 0;
    double fragment_mean = 0;
    double fragment_stdev = 0;
    double fragment_sigma = 10;
    bool fragment_orientation = false;
    bool fragment_direction = true;
    float chance_match = 5e-4;
    bool use_fast_reseed = true;
    float drop_chain = 0.45;
    float mq_overlap = 0.0;
    int kmer_size = 0; // if we set to positive, we'd revert to the old kmer based mapper
    int kmer_stride = 0;
    int pair_window = 64; // unused
    int mate_rescues = 64;
    bool fixed_fragment_model = false;
    bool print_fragment_model = false;
    int fragment_model_update = 10;
    bool acyclic_graph = false;
    bool patch_alignments = true;
    int min_banded_mq = 0;
    int max_sub_mem_recursion_depth = 2;
    bool xdrop_alignment = false;
    uint32_t max_gap_length = 40;
    bool log_time = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"sequence", required_argument, 0, 's'},
                {"seq-name", required_argument, 0, 'V'},
                {"base-name", required_argument, 0, 'd'},
                {"xg-name", required_argument, 0, 'x'},
                {"gcsa-name", required_argument, 0, 'g'},
                {"gbwt-name", required_argument, 0, '1'},
                {"reads", required_argument, 0, 'T'},
                {"sample", required_argument, 0, 'N'},
                {"read-group", required_argument, 0, 'R'},
                {"hit-max", required_argument, 0, 'c'},
                {"max-multimaps", required_argument, 0, 'M'},
                {"threads", required_argument, 0, 't'},
                {"gam-input", required_argument, 0, 'G'},
                {"output-json", no_argument, 0, 'j'},
                {"hts-input", required_argument, 0, 'b'},
                {"keep-secondary", no_argument, 0, 'K'},
                {"exclude-unaligned", no_argument, 0, OPT_EXCLUDE_UNALIGNED},
                {"fastq", required_argument, 0, 'f'},
                {"fasta", required_argument, 0, 'F'},
                {"interleaved", no_argument, 0, 'i'},
                {"band-width", required_argument, 0, 'w'},
                {"band-overlap", required_argument, 0, 'O'},
                {"band-multi", required_argument, 0, 'B'},
                {"band-jump", required_argument, 0, 'J'},
                {"band-min-mq", required_argument, 0, 'Z'},
                {"min-ident", required_argument, 0, 'P'},
                {"debug", no_argument, 0, 'D'},
                {"min-mem", required_argument, 0, 'k'},
                {"max-mem", required_argument, 0, 'Y'},
                {"reseed-x", required_argument, 0, 'r'},
                {"min-chain", required_argument, 0, 'W'},
                {"fast-reseed", no_argument, 0, '6'},
                {"max-target-x", required_argument, 0, 'H'},
                {"buffer-size", required_argument, 0, '9'},
                {"match", required_argument, 0, 'q'},
                {"mismatch", required_argument, 0, 'z'},
                {"score-matrix", required_argument, 0, OPT_SCORE_MATRIX},
                {"gap-open", required_argument, 0, 'o'},
                {"gap-extend", required_argument, 0, 'y'},
                {"qual-adjust", no_argument, 0, 'A'},
                {"try-up-to", required_argument, 0, 'u'},
                {"compare", no_argument, 0, 'X'},
                {"fragment", required_argument, 0, 'I'},
                {"fragment-x", required_argument, 0, '3'},
                {"full-l-bonus", required_argument, 0, 'L'},
                {"hap-exp", required_argument, 0, 'a'},
                {"recombination-penalty", required_argument, 0, OPT_RECOMBINATION_PENALTY},
                {"alignment-model", required_argument, 0, 'm'},
                {"mem-chance", required_argument, 0, 'e'},
                {"drop-chain", required_argument, 0, 'C'},
                {"mq-overlap", required_argument, 0, 'n'},
                {"try-at-least", required_argument, 0, 'l'},
                {"mq-max", required_argument, 0, 'Q'},
                {"mate-rescues", required_argument, 0, '0'},
                {"approx-mq-cap", required_argument, 0, 'E'},
                {"fixed-frag-model", no_argument, 0, 'U'},
                {"print-frag-model", no_argument, 0, 'p'},
                {"frag-calc", required_argument, 0, '4'},
                {"id-mq-weight", required_argument, 0, '7'},
                {"refpos-table", no_argument, 0, 'v'},
                {"surject-to", required_argument, 0, '5'},
                {"ref-paths", required_argument, 0, OPT_REF_PATHS},
                {"no-patch-aln", no_argument, 0, '8'},
                {"drop-full-l-bonus", no_argument, 0, '2'},
                {"unpaired-cost", required_argument, 0, 'S'},
                {"max-gap-length", required_argument, 0, 1},
                {"xdrop-alignment", no_argument, 0, 2},
                {"gaf", no_argument, 0, '%'},
                {"log-time", no_argument, 0, '^'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:J:Q:d:x:g:1:T:N:R:c:M:t:G:jb:Kf:iw:P:Dk:Y:r:W:6H:Z:q:z:o:y:Au:B:I:S:l:e:C:V:O:L:a:n:E:X:UpF:m:7:v5:824:3:9:0:%^",
                         long_options, &option_index);


        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 's':
            seq = optarg;
            break;

        case 'd':
            db_name = optarg;
            break;

        case 'x':
            xg_name = optarg;
            break;

        case 'g':
            gcsa_name = optarg;
            break;
            
        case '1':
            gbwt_name = optarg;
            break;

        case 'V':
            seq_name = optarg;
            break;

        case 'c':
            hit_max = parse<int>(optarg);
            break;

        case 'M':
            max_multimaps = parse<int>(optarg);
            break;

        case '7':
            identity_weight = parse<double>(optarg);
            break;

        case 'Q':
            max_mapping_quality = parse<int>(optarg);
            break;

        case 'E':
            maybe_mq_threshold = parse<double>(optarg);
            break;

        case 'L':
            full_length_bonus = parse<int>(optarg);
            break;

        case '2':
            strip_bonuses = true;
            break;

        case 'a':
            haplotype_consistency_exponent = parse<double>(optarg);
            break;
            
        case OPT_RECOMBINATION_PENALTY:
            recombination_penalty = parse<double>(optarg);
            break;
        
        case 'm':
            if (string(optarg) == "long") {
                extra_multimaps = 2;
                full_length_bonus = 63;
                match = 1;
                mismatch = 2;
                gap_open = 2;
                gap_extend = 1;
                band_width = 128;
                band_overlap = 32;
            }
            break;

        case 'T':
            read_file = optarg;
            break;

        case 'R':
            read_group = optarg;
            break;

        case 'N':
            sample_name = optarg;
            break;

        case 'b':
            hts_file = optarg;
            break;

        case 'K':
            keep_secondary = true;
            break;

        case OPT_EXCLUDE_UNALIGNED:
            exclude_unaligned = true;
            break;

        case 'f':
            if (fastq1.empty()) fastq1 = optarg;
            else if (fastq2.empty()) fastq2 = optarg;
            else { cerr << "[vg map] error: more than two fastqs specified" << endl; exit(1); }
            break;

        case 'F':
            fasta_file = optarg;
            break;

        case 'i':
            interleaved_input = true;
            break;

        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;

        case 'D':
            debug = true;
            break;

        case 'e':
            chance_match = parse<double>(optarg);
            break;

        case 'C':
            drop_chain = parse<double>(optarg);
            break;

        case 'l':
            min_multimaps = parse<int>(optarg);
            break;

        case 'n':
            mq_overlap = parse<double>(optarg);
            break;

        case 'G':
            gam_input = optarg;
            break;

        case 'j':
            output_format = "JSON";
            break;

        case '%':
            output_format = "GAF";
            break;

        case 'w':
            band_width = parse<int>(optarg);
            band_width = band_width == 0 ? INT_MAX : band_width;
            break;

        case 'O':
            band_overlap = parse<int>(optarg);
            break;

        case 'B':
            band_multimaps = parse<int>(optarg);
            break;

        case 'J':
            max_band_jump = parse<int>(optarg);
            break;

        case 'Z':
            min_banded_mq = parse<int>(optarg);
            break;

        case 'P':
            min_score = parse<double>(optarg);
            break;

        case 'k':
            min_mem_length = parse<int>(optarg);
            break;

        case 'Y':
            max_mem_length = parse<int>(optarg);
            break;

        case 'r':
            mem_reseed_factor = parse<double>(optarg);
            break;

        case 'W':
            min_cluster_length = parse<int>(optarg);
            break;

        case 'H':
            max_target_factor = parse<int>(optarg);
            break;

        case '9':
            buffer_size = parse<int>(optarg);
            break;

        case 'q':
            match = parse<int>(optarg);
            break;

        case 'z':
            mismatch = parse<int>(optarg);
            break;

        case OPT_SCORE_MATRIX:
            matrix_file_name = optarg;
            if (matrix_file_name.empty()) {
                cerr << "error:[vg map] Must provide matrix file with --matrix-file." << endl;
                exit(1);
            }
            break;

        case 'o':
            gap_open = parse<int>(optarg);
            break;

        case 'y':
            gap_extend = parse<int>(optarg);
            break;

        case 'A':
            qual_adjust_alignments = true;
            break;

        case 'u':
            extra_multimaps = parse<int>(optarg);
            break;

        case 'X':
            compare_gam = true;
            output_format = "JSON";
            break;

        case 'v':
            output_format = "TSV";
            break;

        case '5':
            output_format = optarg;
            for (auto& c: output_format) {
                // Convert to upper case
                c = toupper(c);
            }
            if (output_format != "SAM" && output_format != "BAM" && output_format != "CRAM") {
                cerr << "error [vg map] illegal surjection type " << optarg << endl;
                return 1;
            }
            break;
            
        case OPT_REF_PATHS:
            ref_paths_name = optarg;
            break;

        case '8':
            patch_alignments = false;
            break;

        case 'I':
        {
            vector<string> parts = split_delims(string(optarg), ":");
            if (parts.size() == 1) {
                convert(parts[0], fragment_max);
            } else if (parts.size() == 5) {
                convert(parts[0], fragment_size);
                convert(parts[1], fragment_mean);
                convert(parts[2], fragment_stdev);
                convert(parts[3], fragment_orientation);
                convert(parts[4], fragment_direction);
            } else {
                cerr << "error [vg map] expected five :-delimited numbers to --fragment" << endl;
                return 1;
            }
        }
        break;

        case '3':
            fragment_sigma = parse<double>(optarg);
            break;

        case 'S':
            unpaired_penalty = parse<int>(optarg);
            break;

        case '0':
            mate_rescues = parse<int>(optarg);
            break;

        case 'U':
            fixed_fragment_model = true;
            break;

        case 'p':
            print_fragment_model = true;
            break;

        case '4':
            fragment_model_update = parse<int>(optarg);
            break;

        case 1:
            max_gap_length = atoi(optarg);     // fall through
        case 2:
            xdrop_alignment = true;
            break;

        case 'h':
        case '^':
            log_time = true;
            break;
        case '?':
            /* getopt_long already printed an error message. */
            help_map(argv);
            exit(1);
            break;


        default:
            cerr << "Unimplemented option " << (char) c << endl;
            exit(1);
        }
    }

    // Decide if we are outputting to an htslib format
    bool hts_output = (output_format == "SAM" || output_format == "BAM" || output_format == "CRAM");

    if (!ref_paths_name.empty() && !hts_output) {
        cerr << "warning:[vg map] Reference path file (--ref-paths) is only used when output format (--surject-to) is SAM, BAM, or CRAM." << endl;
        ref_paths_name = "";
    }

    if (seq.empty() && read_file.empty() && hts_file.empty() && fastq1.empty() && gam_input.empty() && fasta_file.empty()) {
        cerr << "error:[vg map] A sequence or read file is required when mapping." << endl;
        return 1;
    }

    if (!qual.empty() && (seq.length() != qual.length())) {
        cerr << "error:[vg map] Sequence and base quality string must be the same length." << endl;
        return 1;
    }

    if (qual_adjust_alignments && ((fastq1.empty() && hts_file.empty() && qual.empty() && gam_input.empty()) // must have some quality input
                                   || (!seq.empty() && qual.empty())                                         // can't provide sequence without quality
                                   || !read_file.empty()))                                                   // can't provide sequence list without qualities
    {
        cerr << "error:[vg map] Quality adjusted alignments require base quality scores for all sequences." << endl;
        return 1;
    }
    // note: still possible that hts file types don't have quality, but have to check the file to know
    
    MappingQualityMethod mapping_quality_method = Approx;

    string file_name;
    if (optind < argc) {
        file_name = get_input_file_name(optind, argc, argv);
    }

    if (gcsa_name.empty() && !file_name.empty()) {
        gcsa_name = file_name + gcsa::GCSA::EXTENSION;
    }
    
    if (gbwt_name.empty() && !file_name.empty()) {
        gbwt_name = file_name + gbwt::GBWT::EXTENSION;
    }

    if (xg_name.empty() && !file_name.empty()) {
        xg_name = file_name + ".xg";
    }

    if (!db_name.empty()) {
        xg_name = db_name + ".xg";
        gcsa_name = db_name + gcsa::GCSA::EXTENSION;
        gbwt_name = db_name + gbwt::GBWT::EXTENSION;
    }

    if (band_overlap == -1) {
        band_overlap = band_width/8;
    }

    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Configure its temp directory to the system temp directory
    gcsa::TempFile::setDirectory(temp_file::get_dir());

    // Load up our indexes.
    PathPositionHandleGraph* xgidx = nullptr;
    unique_ptr<gcsa::GCSA> gcsa;
    unique_ptr<gcsa::LCPArray> lcp;
    unique_ptr<gbwt::GBWT> gbwt;
    // Used only for memory management:
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::PathPositionVectorizableOverlayHelper overlay_helper;
    
    // One of them may be used to provide haplotype scores
    haplo::ScoreProvider* haplo_score_provider = nullptr;
    
    if(!xg_name.empty()) {
        // We have an xg index!

        // We try opening the file, and then see if it worked
        ifstream xg_stream(xg_name);
        if (!xg_stream) {
            cerr << "Error[vg map]: Unable to open xg file \"" << xg_name << "\"" << endl;
            exit(1);
        }
        xg_stream.close();
        
        // TODO: tell when the user asked for an XG vs. when we guessed one,
        // and error when the user asked for one and we can't find it.
        if(debug) {
            cerr << "Loading xg index " << xg_name << "..." << endl;
        }
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        xgidx = dynamic_cast<PathPositionHandleGraph*>(overlay_helper.apply(path_handle_graph.get()));
    }

    ifstream gcsa_stream(gcsa_name);
    if(gcsa_stream) {
        // We have a GCSA index too!
        if(debug) {
            cerr << "Loading GCSA2 index " << gcsa_name << "..." << endl;
        }
        gcsa = vg::io::VPKG::load_one<gcsa::GCSA>(gcsa_stream);
    }

    string lcp_name = gcsa_name + ".lcp";
    ifstream lcp_stream(lcp_name);
    if (lcp_stream) {
        if(debug) {
            cerr << "Loading LCP index " << lcp_name << "..." << endl;
        }
        lcp = vg::io::VPKG::load_one<gcsa::LCPArray>(lcp_stream);
    }
    
    ifstream gbwt_stream(gbwt_name);
    if(gbwt_stream) {
        // We have a GBWT index too!
        if(debug) {
            cerr << "Loading GBWT haplotype index " << gbwt_name << "..." << endl;
        }
        
        gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
        
        // We want to use this for haplotype scoring
        haplo_score_provider = new haplo::GBWTScoreProvider<gbwt::GBWT>(*gbwt);
    }

    ifstream matrix_stream;
    if (!matrix_file_name.empty()) {
      matrix_stream.open(matrix_file_name);
      if (!matrix_stream) {
          cerr << "error:[vg map] Cannot open scoring matrix file " << matrix_file_name << endl;
          exit(1);
      }
    }

    thread_count = get_thread_count();

    // TODO: We need a Mapper for every thread because the Mapper's fragment
    // length distribution isn't yet thread safe.  
    vector<Mapper*> mapper;
    mapper.resize(thread_count);
    
    // When outputting single-ended alignments, we need an empty vector to pass around
    vector<Alignment> empty_alns;
   
    // Look up all the paths we might need to surject to.
    vector<tuple<path_handle_t, size_t, size_t>> paths;
    if (hts_output) {
        paths = get_sequence_dictionary(ref_paths_name, *xgidx);
    }
    
    // Set up output to an emitter that will handle serialization and surjection
    unique_ptr<vg::io::AlignmentEmitter> alignment_emitter = get_alignment_emitter("-", output_format, paths, thread_count, xgidx);

    // We have one function to dump alignments into
    auto output_alignments = [&](vector<Alignment>& alns1, vector<Alignment>& alns2) {
        if (alns2.empty()) {
            // Single-ended read
            alignment_emitter->emit_mapped_single(std::move(alns1));
        } else {
            // Paired reads
            if (hts_output) {
                // We need a tlen limit for flags
                
                // Look up the paired end distribution stats for deciding if reads are propelry paired
                auto& stats = mapper[omp_get_thread_num()]->frag_stats;
                // Put a proper pair bound at 6 std devs.
                // If distribution hasn't been computed yet, this comes out 0 and no bound is applied.
                int64_t tlen_limit = stats.cached_fragment_length_mean + 6 * stats.cached_fragment_length_stdev;
                
                // Send the tlen limit when emitting
                alignment_emitter->emit_mapped_pair(std::move(alns1), std::move(alns2), tlen_limit);
            } else {
                // No need for a tlen limit
                alignment_emitter->emit_mapped_pair(std::move(alns1), std::move(alns2));
            }
        }
    };

    for (int i = 0; i < thread_count; ++i) {
        Mapper* m = nullptr;
        if(xgidx && gcsa.get() && lcp.get()) {
            // We have the xg and GCSA indexes, so use them
            m = new Mapper(xgidx, gcsa.get(), lcp.get(), haplo_score_provider);
        } else {
            // Can't continue with null
            throw runtime_error("Need XG, GCSA, and LCP to create a Mapper");
        }
        m->hit_max = hit_max;
        m->max_multimaps = max_multimaps;
        m->min_multimaps = max(min_multimaps, max_multimaps);
        m->band_multimaps = band_multimaps;
        m->min_banded_mq = min_banded_mq;
        m->maybe_mq_threshold = maybe_mq_threshold;
        m->exclude_unaligned = exclude_unaligned;
        m->debug = debug;
        m->min_identity = min_score;
        m->drop_chain = drop_chain;
        m->mq_overlap = mq_overlap;
        m->min_mem_length = (min_mem_length > 0 ? min_mem_length
                             : m->random_match_length(chance_match));
        m->min_cluster_length = min_cluster_length;
        m->mem_reseed_length = round(mem_reseed_factor * m->min_mem_length);
        if (debug && i == 0) {
            cerr << "[vg map] : min_mem_length = " << m->min_mem_length
                 << ", mem_reseed_length = " << m->mem_reseed_length
                 << ", min_cluster_length = " << m->min_cluster_length << endl;
        }
        m->fast_reseed = use_fast_reseed;
        m->max_sub_mem_recursion_depth = max_sub_mem_recursion_depth;
        m->max_target_factor = max_target_factor;
        if (matrix_stream.is_open()) {
            m->set_alignment_scores(matrix_stream, gap_open, gap_extend, full_length_bonus, haplotype_consistency_exponent);
            // reset the stream for the next Mapper
            matrix_stream.seekg(0);
        }
        else {
            m->set_alignment_scores(match, mismatch, gap_open, gap_extend, full_length_bonus, haplotype_consistency_exponent);
        }
        m->strip_bonuses = strip_bonuses;
        m->max_xdrop_gap_length = max_gap_length;
        m->adjust_alignments_for_base_quality = qual_adjust_alignments;
        m->extra_multimaps = extra_multimaps;
        m->mapping_quality_method = mapping_quality_method;
        m->recombination_penalty = recombination_penalty;
        m->always_rescue = always_rescue;
        m->frag_stats.fixed_fragment_model = fixed_fragment_model;
        m->frag_stats.fragment_max = fragment_max;
        m->frag_stats.fragment_sigma = fragment_sigma;
        m->unpaired_penalty = unpaired_penalty;
        if (fragment_mean) {
            m->frag_stats.fragment_size = fragment_size;
            m->frag_stats.cached_fragment_length_mean = fragment_mean;
            m->frag_stats.cached_fragment_length_stdev = fragment_stdev;
            m->frag_stats.cached_fragment_orientation_same = fragment_orientation;
            m->frag_stats.cached_fragment_direction = fragment_direction;
        }
        m->frag_stats.fragment_model_update_interval = fragment_model_update;
        m->max_mapping_quality = max_mapping_quality;
        m->mate_rescues = mate_rescues;
        m->max_band_jump = max_band_jump;
        m->identity_weight = identity_weight;
        m->assume_acyclic = acyclic_graph;
        m->patch_alignments = patch_alignments;
        mapper[i] = m;
    }
    vector<size_t> reads_mapped_by_thread(thread_count, 0);

    std::chrono::time_point<std::chrono::system_clock> init = std::chrono::system_clock::now();

    if (!seq.empty()) {
        int tid = omp_get_thread_num();

        Alignment unaligned;
        unaligned.set_sequence(seq);

        if (!qual.empty()) {
            unaligned.set_quality(qual);
        }
        
        vector<Alignment> alignments = mapper[tid]->align_multi(unaligned,
                                                                kmer_size,
                                                                kmer_stride,
                                                                max_mem_length,
                                                                band_width,
                                                                band_overlap,
                                                                xdrop_alignment);

        if(alignments.size() == 0 && !exclude_unaligned) {
            // If we didn't have any alignments, report the unaligned alignment
            alignments.push_back(unaligned);
        }


        for(auto& alignment : alignments) {
            if (!sample_name.empty()) alignment.set_sample_name(sample_name);
            if (!read_group.empty()) alignment.set_read_group(read_group);
            if (!seq_name.empty()) alignment.set_name(seq_name);
        }

        // Output the alignments in the correct format, possibly surjecting.
        output_alignments(alignments, empty_alns);
        reads_mapped_by_thread[tid] += 1;
    }

    if (!read_file.empty()) {
        ifstream in(read_file);
        bool more_data = in.good();
#pragma omp parallel shared(in)
        {
            string line;
            int tid = omp_get_thread_num();
            while (in.good()) {
                line.clear();
#pragma omp critical (readq)
                {
                    std::getline(in,line);
                }
                if (!line.empty()) {
                    // Make an alignment
                    Alignment unaligned;
                    unaligned.set_sequence(line);
                    vector<Alignment> alignments = mapper[tid]->align_multi(unaligned,
                                                                            kmer_size,
                                                                            kmer_stride,
                                                                            max_mem_length,
                                                                            band_width,
                                                                            band_overlap,
                                                                            xdrop_alignment);
                    

                    for(auto& alignment : alignments) {
                        // Set the alignment metadata
                        if (!sample_name.empty()) alignment.set_sample_name(sample_name);
                        if (!read_group.empty()) alignment.set_read_group(read_group);
                    }


                    // Output the alignments in the correct format, possibly surjecting. 
                    output_alignments(alignments, empty_alns);
                }
                reads_mapped_by_thread[tid] += 1;
            }
        }
    }

    if (!fasta_file.empty()) {
        FastaReference ref;
        ref.open(fasta_file);
        auto align_seq = [&](const string& name, const string& seq) {
            if (!seq.empty()) {
                // Make an alignment
                Alignment unaligned;
                unaligned.set_sequence(seq);
                unaligned.set_name(name);
                int tid = omp_get_thread_num();
                vector<Alignment> alignments = mapper[tid]->align_multi(unaligned,
                                                                        kmer_size,
                                                                        kmer_stride,
                                                                        max_mem_length,
                                                                        band_width,
                                                                        band_overlap,
                                                                        xdrop_alignment);
                
                for(auto& alignment : alignments) {
                    // Set the alignment metadata
                    if (!sample_name.empty()) alignment.set_sample_name(sample_name);
                    if (!read_group.empty()) alignment.set_read_group(read_group);
                }
                // Output the alignments in the correct format, possibly surjecting.
                output_alignments(alignments, empty_alns);

                reads_mapped_by_thread[tid] += 1;
            }
        };
#pragma omp parallel for
        for (size_t i = 0; i < ref.index->sequenceNames.size(); ++i) {
            auto& name = ref.index->sequenceNames[i];
            string seq = nonATGCNtoN(toUppercase(ref.getSequence(name)));
            align_seq(name, seq);
        }
    }

    if (!hts_file.empty()) {
        function<void(Alignment&)> lambda = [&](Alignment& alignment) {
            if(alignment.is_secondary() && !keep_secondary) {
                // Skip over secondary alignments in the input; we don't want several output mappings for each input *mapping*.
                return;
            }

            int tid = omp_get_thread_num();
            vector<Alignment> alignments = mapper[tid]->align_multi(alignment,
                                                                    kmer_size,
                                                                    kmer_stride,
                                                                    max_mem_length,
                                                                    band_width,
                                                                    band_overlap,
                                                                    xdrop_alignment);
                                                                    
            for(auto& alignment : alignments) {
                // Set the alignment metadata
                if (!sample_name.empty()) alignment.set_sample_name(sample_name);
                if (!read_group.empty()) alignment.set_read_group(read_group);
            }

            // Output the alignments in JSON or protobuf as appropriate.
            output_alignments(alignments, empty_alns);

            reads_mapped_by_thread[tid] += 1;
        };
        // run
        hts_for_each_parallel(hts_file, lambda);
    }

    if (!fastq1.empty()) {
        if (interleaved_input) {
            // paired interleaved
            auto output_func = [&](Alignment& aln1,
                                   Alignment& aln2,
                                   pair<vector<Alignment>, vector<Alignment>>& alnp) {
                
                if (!print_fragment_model) {
                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alnp.first, alnp.second);
                }
            };
            
            function<void(Alignment&,Alignment&)> lambda = [&](Alignment& aln1, Alignment& aln2) {
                auto our_mapper = mapper[omp_get_thread_num()];
                bool queued_resolve_later = false;
                auto alnp = our_mapper->align_paired_multi(aln1,
                                                           aln2,
                                                           queued_resolve_later,
                                                           max_mem_length,
                                                           top_pairs_only,
                                                           false,
                                                           xdrop_alignment);
                if (!queued_resolve_later) {
                    output_func(aln1, aln2, alnp);
                    // check if we should try to align the queued alignments
                    if (our_mapper->frag_stats.fragment_size != 0
                        && !our_mapper->imperfect_pairs_to_retry.empty()) {
                        int i = 0;
                        for (auto p : our_mapper->imperfect_pairs_to_retry) {
                            auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                                       queued_resolve_later,
                                                                       max_mem_length,
                                                                       top_pairs_only,
                                                                       true,
                                                                       xdrop_alignment);
                            output_func(p.first, p.second, alnp);
                        }
                        our_mapper->imperfect_pairs_to_retry.clear();
                    }
                }

                reads_mapped_by_thread[omp_get_thread_num()] += 2;
            };
            fastq_paired_interleaved_for_each_parallel(fastq1, lambda);
#pragma omp parallel
            { // clean up buffered alignments that weren't perfect
                auto our_mapper = mapper[omp_get_thread_num()];
                // if we haven't yet computed these, assume we couldn't get an estimate for fragment size
                our_mapper->frag_stats.fragment_size = fragment_max;
                for (auto p : our_mapper->imperfect_pairs_to_retry) {
                    bool queued_resolve_later = false;
                    auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                               queued_resolve_later,
                                                               max_mem_length,
                                                               top_pairs_only,
                                                               true,
                                                               xdrop_alignment);
                    output_func(p.first, p.second, alnp);
                }
                our_mapper->imperfect_pairs_to_retry.clear();
            }
        } else if (fastq2.empty()) {
            // single
            function<void(Alignment&)> lambda = [&](Alignment& alignment) {
                        int tid = omp_get_thread_num();
                        vector<Alignment> alignments = mapper[tid]->align_multi(alignment,
                                                                                kmer_size,
                                                                                kmer_stride,
                                                                                max_mem_length,
                                                                                band_width,
                                                                                band_overlap,
                                                                                xdrop_alignment);
                        //cerr << "This is just before output_alignments" << alignment.DebugString() << endl;
                        output_alignments(alignments, empty_alns);
                        reads_mapped_by_thread[tid] += 1;
                    };
            fastq_unpaired_for_each_parallel(fastq1, lambda);
        } else {
            // paired two-file
            auto output_func = [&](Alignment& aln1,
                                   Alignment& aln2,
                                   pair<vector<Alignment>, vector<Alignment>>& alnp) {
                // Make sure we have unaligned "alignments" for things that don't align.
                // Output the alignments in JSON or protobuf as appropriate.
                if (!print_fragment_model) {
                    output_alignments(alnp.first, alnp.second);
                }
            };
            function<void(Alignment&,Alignment&)> lambda = [&](Alignment& aln1, Alignment& aln2) {
                auto our_mapper = mapper[omp_get_thread_num()];
                bool queued_resolve_later = false;
                auto alnp = our_mapper->align_paired_multi(aln1,
                                                           aln2,
                                                           queued_resolve_later,
                                                           max_mem_length,
                                                           top_pairs_only,
                                                           false,
                                                           xdrop_alignment);
                if (!queued_resolve_later) {
                    output_func(aln1, aln2, alnp);
                    // check if we should try to align the queued alignments
                    if (our_mapper->frag_stats.fragment_size != 0
                        && !our_mapper->imperfect_pairs_to_retry.empty()) {
                        int i = 0;
                        for (auto p : our_mapper->imperfect_pairs_to_retry) {
                            auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                                       queued_resolve_later,
                                                                       max_mem_length,
                                                                       top_pairs_only,
                                                                       true,
                                                                       xdrop_alignment);
                            output_func(p.first, p.second, alnp);
                        }
                        our_mapper->imperfect_pairs_to_retry.clear();
                    }
                }

                reads_mapped_by_thread[omp_get_thread_num()] += 2;
            };
            fastq_paired_two_files_for_each_parallel(fastq1, fastq2, lambda);
#pragma omp parallel
            {
                auto our_mapper = mapper[omp_get_thread_num()];
                our_mapper->frag_stats.fragment_size = fragment_max;
                for (auto p : our_mapper->imperfect_pairs_to_retry) {
                    bool queued_resolve_later = false;
                    auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                               queued_resolve_later,
                                                               max_mem_length,
                                                               top_pairs_only,
                                                               true,
                                                               xdrop_alignment);
                    output_func(p.first, p.second, alnp);
                }
                our_mapper->imperfect_pairs_to_retry.clear();

                reads_mapped_by_thread[omp_get_thread_num()] += 2;
            }
        }
    }

    if (!gam_input.empty()) {
        ifstream gam_in(gam_input);
        if (interleaved_input) {
            // Paired-end GAM input
            auto output_func = [&] (Alignment& aln1,
                                    Alignment& aln2,
                                    pair<vector<Alignment>, vector<Alignment>>& alnp) {
                if (print_fragment_model) {
                    // do nothing
                } else {
                    // Output the alignments in JSON or protobuf as appropriate.
                    if (compare_gam) {
                        alnp.first.front().set_correct(overlap(aln1.path(), alnp.first.front().path()));
                        alnp.second.front().set_correct(overlap(aln2.path(), alnp.second.front().path()));
                        alignment_set_distance_to_correct(alnp.first.front(), aln1);
                        alignment_set_distance_to_correct(alnp.second.front(), aln2);
                    }
                    output_alignments(alnp.first, alnp.second);
                }
            };
            function<void(Alignment&,Alignment&)> lambda = [&](Alignment& aln1, Alignment& aln2) {
                auto our_mapper = mapper[omp_get_thread_num()];
                bool queued_resolve_later = false;
                auto alnp = our_mapper->align_paired_multi(aln1,
                                                           aln2,
                                                           queued_resolve_later,
                                                           max_mem_length,
                                                           top_pairs_only,
                                                           false,
                                                           xdrop_alignment);
                if (!queued_resolve_later) {
                    output_func(aln1, aln2, alnp);
                    // check if we should try to align the queued alignments
                    if (our_mapper->frag_stats.fragment_size != 0
                        && !our_mapper->imperfect_pairs_to_retry.empty()) {
                        int i = 0;
                        for (auto p : our_mapper->imperfect_pairs_to_retry) {
                            auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                                       queued_resolve_later,
                                                                       max_mem_length,
                                                                       top_pairs_only,
                                                                       true,
                                                                       xdrop_alignment);
                            output_func(p.first, p.second, alnp);
                        }
                        our_mapper->imperfect_pairs_to_retry.clear();
                    }
                }
                reads_mapped_by_thread[omp_get_thread_num()] += 2;
            };
            vg::io::for_each_interleaved_pair_parallel(gam_in, lambda);
#pragma omp parallel
            {
                auto our_mapper = mapper[omp_get_thread_num()];
                our_mapper->frag_stats.fragment_size = fragment_max;
                for (auto p : our_mapper->imperfect_pairs_to_retry) {
                    bool queued_resolve_later = false;
                    auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                               queued_resolve_later,
                                                               max_mem_length,
                                                               top_pairs_only,
                                                               true,
                                                               xdrop_alignment);
                    output_func(p.first, p.second, alnp);
                }
                our_mapper->imperfect_pairs_to_retry.clear();
            }
        } else {
            // Processing single-end GAM input
            function<void(Alignment&)> lambda = [&](Alignment& alignment) {
                int tid = omp_get_thread_num();
                vector<Alignment> alignments = mapper[tid]->align_multi(alignment,
                                                                        kmer_size,
                                                                        kmer_stride,
                                                                        max_mem_length,
                                                                        band_width,
                                                                        band_overlap,
                                                                        xdrop_alignment);
                if (compare_gam) {
                    // Compare against true input at mapping time
                    alignments.front().set_correct(overlap(alignment.path(), alignments.front().path()));
                    alignment_set_distance_to_correct(alignments.front(), alignment);
                }
                output_alignments(alignments, empty_alns);
                reads_mapped_by_thread[tid] += 1;
            };
            vg::io::for_each_parallel(gam_in, lambda);
        }
        gam_in.close();
    }

    if (print_fragment_model) {
        if (mapper[0]->frag_stats.fragment_size) {
            // we've calculated our fragment size, so print it and bail out
            cout << mapper[0]->frag_stats.fragment_model_str() << endl;
        } else {
            cerr << "[vg map] Error: could not calculate fragment model" << endl;
        }
    }

    if (haplo_score_provider) {
        delete haplo_score_provider;
        haplo_score_provider = nullptr;
    }
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> mapping_seconds = end - init;
    std::chrono::duration<double> index_load_seconds = init - launch;

    if (log_time){

        size_t total_reads_mapped = 0;
        for (auto& reads_mapped : reads_mapped_by_thread) {
            total_reads_mapped += reads_mapped;
        }
    
        double reads_per_second_per_thread = total_reads_mapped / (mapping_seconds.count() * thread_count);
        cerr << "Index load time: " << index_load_seconds.count() << endl;
        cerr << "Mapped " << total_reads_mapped << " reads" << endl;
        cerr << "Mapping speed: " << reads_per_second_per_thread << " reads per second per thread" << endl; 
    }
    
    cout.flush();

    // clean up our mappers
    for (uint64_t i = 0; i < mapper.size(); ++i) {
        delete mapper[i];
    }

    return 0;

}

static Subcommand vg_map("map", "MEM-based read alignment", PIPELINE, 5, main_map);

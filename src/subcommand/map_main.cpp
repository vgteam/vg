#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "../mapper.hpp"
#include "../stream.hpp"

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;

void help_map(char** argv) {
    cerr << "usage: " << argv[0] << " map [options] -d idxbase -f in1.fq [-f in2.fq] >aln.gam" << endl
         << "Align reads to a graph." << endl
         << endl
         << "graph/index:" << endl
         << "    -d, --base-name BASE    use BASE.xg and BASE.gcsa as the input index pair" << endl
         << "    -x, --xg-name FILE      use this xg index (defaults to <graph>.vg.xg)" << endl
         << "    -g, --gcsa-name FILE    use this GCSA2 index (defaults to <graph>" << gcsa::GCSA::EXTENSION << ")" << endl
         << "    -1, --gbwt-name FILE    use this GBWT haplotype index (defaults to <graph>"<<gbwt::GBWT::EXTENSION << ")" << endl
         << "algorithm:" << endl
         << "    -t, --threads N         number of compute threads to use" << endl
         << "    -k, --min-mem INT       minimum MEM length (if 0 estimate via -e) [0]" << endl
         << "    -e, --mem-chance FLOAT  set {-k} such that this fraction of {-k} length hits will by chance [5e-4]" << endl
         << "    -c, --hit-max N         ignore MEMs who have >N hits in our index (0 for no limit) [8192]" << endl
         << "    -Y, --max-mem INT       ignore mems longer than this length (unset if 0) [0]" << endl
         << "    -r, --reseed-x FLOAT    look for internal seeds inside a seed longer than FLOAT*--min-seed [1.5]" << endl
         << "    -u, --try-up-to INT     attempt to align up to the INT best candidate chains of seeds (1/2 for paired) [128]" << endl
         << "    -l, --try-at-least INT  attempt to align at least the INT best candidate chains of seeds [1]" << endl
         << "    -E, --approx-mq-cap INT weight MQ by suffix tree based estimate when estimate less than FLOAT [0]" << endl
         << "    --id-mq-weight N        scale mapping quality by the alignment score identity to this power [2]" << endl
         << "    -W, --min-chain INT     discard a chain if seeded bases shorter than INT [0]" << endl
         << "    -C, --drop-chain FLOAT  drop chains shorter than FLOAT fraction of the longest overlapping chain [0.45]" << endl
         << "    -n, --mq-overlap FLOAT  scale MQ by count of alignments with this overlap in the query with the primary [0]" << endl
         << "    -P, --min-ident FLOAT   accept alignment only if the alignment identity is >= FLOAT [0]" << endl
         << "    -H, --max-target-x N    skip cluster subgraphs with length > N*read_length [100]" << endl
         << "    -m, --acyclic-graph     improves runtime when the graph is acyclic" << endl
         << "    -w, --band-width INT    band width for long read alignment [256]" << endl
         << "    -J, --band-jump INT     the maximum jump we can see between bands (maximum length variant we can detect) [{-w}]" << endl
         << "    -B, --band-multi INT    consider this many alignments of each band in banded alignment [1]" << endl
         << "    -Z, --band-min-mq INT   treat bands with less than this MQ as unaligned [0]" << endl
         << "    -I, --fragment STR      fragment length distribution specification STR=m:μ:σ:o:d [5000:0:0:0:1]" << endl
         << "                            max, mean, stdev, orientation (1=same, 0=flip), direction (1=forward, 0=backward)" << endl
         << "    -U, --fixed-frag-model  don't learn the pair fragment model online, use {-I} without update" << endl
         << "    -p, --print-frag-model  suppress alignment output and print the fragment model on stdout as per {-I} format" << endl
         << "    --frag-calc INT         update the fragment model every INT perfect pairs [10]" << endl
         << "    --fragment-x FLOAT      calculate max fragment size as frag_mean+frag_sd*FLOAT [10]" << endl
         << "    -O, --mate-rescues INT  attempt up to INT mate rescues per pair [64]" << endl
         << "    -S, --unpaired-cost INT penalty for an unpaired read pair [17]" << endl
         << "    --no-patch-aln          do not patch banded alignments by locally aligning unaligned regions" << endl
         << "scoring:" << endl
         << "    -q, --match INT         use this match score [1]" << endl
         << "    -z, --mismatch INT      use this mismatch penalty [4]" << endl
         << "    -o, --gap-open INT      use this gap open penalty [6]" << endl
         << "    -y, --gap-extend INT    use this gap extension penalty [1]" << endl
         << "    -L, --full-l-bonus INT  the full-length alignment bonus [5]" << endl
         << "    --drop-full-l-bonus     remove the full length bonus from the score before sorting and MQ calculation" << endl
         << "    -a, --hap-exp FLOAT     the exponent for haplotype consistency likelihood in alignment score [1]" << endl
         << "    -A, --qual-adjust       perform base quality adjusted alignments (requires base quality input)" << endl
         << "input:" << endl
         << "    -s, --sequence STR      align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -V, --seq-name STR      name the sequence using this value (for graph modification with new named paths)" << endl
         << "    -T, --reads FILE        take reads (one per line) from FILE, write alignments to stdout" << endl
         << "    -b, --hts-input FILE    align reads from htslib-compatible FILE (BAM/CRAM/SAM) stdin (-), alignments to stdout" << endl
         << "    -G, --gam-input FILE    realign GAM input" << endl
         << "    -f, --fastq FILE        input fastq or (2-line format) fasta, possibly compressed, two are allowed, one for each mate" << endl
         << "    -F, --fasta FILE        align the sequences in a FASTA file that may have multiple lines per reference sequence" << endl
         << "    -i, --interleaved       fastq or GAM is interleaved paired-ended" << endl
         << "    -N, --sample NAME       for --reads input, add this sample" << endl
         << "    -R, --read-group NAME   for --reads input, add this read group" << endl
         << "output:" << endl
         << "    -j, --output-json       output JSON rather than an alignment stream (helpful for debugging)" << endl
         << "    --surject-to TYPE       surject the output into the graph's paths, writing TYPE := bam |sam | cram" << endl
         << "    --buffer-size INT       buffer this many alignments together before outputting in GAM [512]" << endl
         << "    -X, --compare           realign GAM input (-G), writing alignment with \"correct\" field set to overlap with input" << endl
         << "    -v, --refpos-table      for efficient testing output a table of name, chr, pos, mq, score" << endl
         << "    -K, --keep-secondary    produce alignments for secondary input alignments in addition to primary ones" << endl
         << "    -M, --max-multimaps INT produce up to INT alignments for each read [1]" << endl
         << "    -Q, --mq-max INT        cap the mapping quality at INT [60]" << endl
         << "    -D, --debug             print debugging information about alignment to stderr" << endl;

}

int main_map(int argc, char** argv) {

    if (argc == 2) {
        help_map(argv);
        return 1;
    }
    
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
    int hit_max = 8192;
    int max_multimaps = 1;
    int thread_count = 1;
    bool output_json = false;
    string surject_type;
    bool debug = false;
    float min_score = 0;
    string sample_name;
    string read_group;
    string fastq1, fastq2;
    bool interleaved_input = false;
    int band_width = 256;
    int band_multimaps = 1;
    int max_band_jump = -1;
    bool always_rescue = false;
    bool top_pairs_only = false;
    int max_mem_length = 0;
    int min_mem_length = 0;
    int min_cluster_length = 0;
    float mem_reseed_factor = 1.5;
    int max_target_factor = 100;
    int buffer_size = 512;
    int8_t match = 1;
    int8_t mismatch = 4;
    int8_t gap_open = 6;
    int8_t gap_extend = 1;
    int8_t full_length_bonus = 5;
    int unpaired_penalty = 17;
    double haplotype_consistency_exponent = 1;
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
    bool refpos_table = false;
    bool patch_alignments = true;
    int min_banded_mq = 0;

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
                {"fastq", required_argument, 0, 'f'},
                {"fasta", required_argument, 0, 'F'},
                {"interleaved", no_argument, 0, 'i'},
                {"band-width", required_argument, 0, 'w'},
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
                {"gap-open", required_argument, 0, 'o'},
                {"gap-extend", required_argument, 0, 'y'},
                {"qual-adjust", no_argument, 0, 'A'},
                {"try-up-to", required_argument, 0, 'u'},
                {"compare", no_argument, 0, 'X'},
                {"fragment", required_argument, 0, 'I'},
                {"fragment-x", required_argument, 0, '3'},
                {"full-l-bonus", required_argument, 0, 'L'},
                {"hap-exp", required_argument, 0, 'a'},
                {"acyclic-graph", no_argument, 0, 'm'},
                {"mem-chance", required_argument, 0, 'e'},
                {"drop-chain", required_argument, 0, 'C'},
                {"mq-overlap", required_argument, 0, 'n'},
                {"try-at-least", required_argument, 0, 'l'},
                {"mq-max", required_argument, 0, 'Q'},
                {"mate-rescues", required_argument, 0, 'O'},
                {"approx-mq-cap", required_argument, 0, 'E'},
                {"fixed-frag-model", no_argument, 0, 'U'},
                {"print-frag-model", no_argument, 0, 'p'},
                {"frag-calc", required_argument, 0, '4'},
                {"id-mq-weight", required_argument, 0, '7'},
                {"refpos-table", no_argument, 0, 'v'},
                {"surject-to", required_argument, 0, '5'},
                {"no-patch-aln", no_argument, 0, '8'},
                {"drop-full-l-bonus", no_argument, 0, '2'},
                {"unpaired-cost", required_argument, 0, 'S'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:J:Q:d:x:g:1:T:N:R:c:M:t:G:jb:Kf:iw:P:Dk:Y:r:W:6H:Z:q:z:o:y:Au:B:I:S:l:e:C:V:O:L:a:n:E:X:UpF:m7:v5:824:3:9:",
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
            hit_max = atoi(optarg);
            break;

        case 'M':
            max_multimaps = atoi(optarg);
            break;

        case '7':
            identity_weight = atof(optarg);
            break;

        case 'Q':
            max_mapping_quality = atoi(optarg);
            break;

        case 'E':
            maybe_mq_threshold = atof(optarg);
            break;

        case 'L':
            full_length_bonus = atoi(optarg);
            break;

        case '2':
            strip_bonuses = true;
            break;

        case 'a':
            haplotype_consistency_exponent = atof(optarg);
            break;
        
        case 'm':
            acyclic_graph = true;
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
            omp_set_num_threads(atoi(optarg));
            break;

        case 'D':
            debug = true;
            break;

        case 'e':
            chance_match = atof(optarg);
            break;

        case 'C':
            drop_chain = atof(optarg);
            break;

        case 'l':
            min_multimaps = atoi(optarg);
            break;

        case 'n':
            mq_overlap = atof(optarg);
            break;

        case 'G':
            gam_input = optarg;
            break;

        case 'j':
            output_json = true;
            break;

        case 'w':
            band_width = atoi(optarg);
            break;

        case 'B':
            band_multimaps = atoi(optarg);
            break;

        case 'J':
            max_band_jump = atoi(optarg);
            break;

        case 'Z':
            min_banded_mq = atoi(optarg);
            break;

        case 'P':
            min_score = atof(optarg);
            break;

        case 'k':
            min_mem_length = atoi(optarg);
            break;

        case 'Y':
            max_mem_length = atoi(optarg);
            break;

        case 'r':
            mem_reseed_factor = atof(optarg);
            break;

        case 'W':
            min_cluster_length = atoi(optarg);
            break;

        case 'H':
            max_target_factor = atoi(optarg);
            break;

        case '9':
            buffer_size = atoi(optarg);
            break;

        case 'q':
            match = atoi(optarg);
            break;

        case 'z':
            mismatch = atoi(optarg);
            break;

        case 'o':
            gap_open = atoi(optarg);
            break;

        case 'y':
            gap_extend = atoi(optarg);
            break;

        case 'A':
            qual_adjust_alignments = true;
            break;

        case 'u':
            extra_multimaps = atoi(optarg);
            break;

        case 'X':
            compare_gam = true;
            output_json = true;
            break;

        case 'v':
            refpos_table = true;
            break;

        case '5':
            surject_type = optarg;
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
            fragment_sigma = atof(optarg);
            break;

        case 'S':
            unpaired_penalty = atoi(optarg);
            break;

        case 'O':
            mate_rescues = atoi(optarg);
            break;

        case 'U':
            fixed_fragment_model = true;
            break;

        case 'p':
            print_fragment_model = true;
            break;

        case '4':
            fragment_model_update = atoi(optarg);
            break;

        case 'h':
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

    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Configure its temp directory to the system temp directory
    gcsa::TempFile::setDirectory(find_temp_dir());

    // Load up our indexes.
    xg::XG* xgidx = nullptr;
    gcsa::GCSA* gcsa = nullptr;
    gcsa::LCPArray* lcp = nullptr;
    gbwt::GBWT* gbwt = nullptr;
    
    // One of them may be used to provide haplotype scores
    haplo::ScoreProvider* haplo_score_provider = nullptr;

    // We try opening the file, and then see if it worked
    ifstream xg_stream(xg_name);

    if(xg_stream) {
        // We have an xg index!
        
        // TODO: tell when the user asked for an XG vs. when we guessed one,
        // and error when the user asked for one and we can't find it.
        if(debug) {
            cerr << "Loading xg index " << xg_name << "..." << endl;
        }
        xgidx = new xg::XG(xg_stream);
        
        // TODO: Support haplo::XGScoreProvider?
    }

    ifstream gcsa_stream(gcsa_name);
    if(gcsa_stream) {
        // We have a GCSA index too!
        if(debug) {
            cerr << "Loading GCSA2 index " << gcsa_name << "..." << endl;
        }
        gcsa = new gcsa::GCSA();
        gcsa->load(gcsa_stream);
    }

    string lcp_name = gcsa_name + ".lcp";
    ifstream lcp_stream(lcp_name);
    if (lcp_stream) {
        if(debug) {
            cerr << "Loading LCP index " << gcsa_name << "..." << endl;
        }
        lcp = new gcsa::LCPArray();
        lcp->load(lcp_stream);
    }
    
    ifstream gbwt_stream(gbwt_name);
    if(gbwt_stream) {
        // We have a GBWT index too!
        if(debug) {
            cerr << "Loading GBWT haplotype index " << gbwt_name << "..." << endl;
        }
        gbwt = new gbwt::GBWT();
        gbwt->load(gbwt_stream);
        
        // We want to use this for haplotype scoring
        haplo_score_provider = new haplo::GBWTScoreProvider<gbwt::GBWT>(*gbwt);
    }

    thread_count = get_thread_count();

    vector<Mapper*> mapper;
    mapper.resize(thread_count);
    vector<vector<Alignment> > output_buffer;
    output_buffer.resize(thread_count);
    vector<Alignment> empty_alns;

    // bam/sam/cram output
    samFile* sam_out = 0;
    int buffer_limit = 100;
    bam_hdr_t* hdr = nullptr;
    int compress_level = 9; // hard coded
    map<string, string> rg_sample;
    string sam_header;

    // if no paths were given take all of those in the index
    set<string> path_names;
    if (!surject_type.empty() && path_names.empty()) {
        for (size_t i = 1; i <= xgidx->path_count; ++i) {
            path_names.insert(xgidx->path_name(i));
        }
    }

    // for SAM header generation
    auto setup_sam_header = [&hdr, &sam_out, &surject_type, &compress_level, &xgidx, &rg_sample, &sam_header] (void) {
#pragma omp critical (hts_header)
        if (!hdr) {
            char out_mode[5];
            string out_format = "";
            strcpy(out_mode, "w");
            if (surject_type == "bam") { out_format = "b"; }
            else if (surject_type == "cram") { out_format = "c"; }
            else { out_format = ""; }
            strcat(out_mode, out_format.c_str());
            if (compress_level >= 0) {
                char tmp[2];
                tmp[0] = compress_level + '0'; tmp[1] = '\0';
                strcat(out_mode, tmp);
            }
            map<string, int64_t> path_length;
            int num_paths = xgidx->max_path_rank();
            for (int i = 1; i <= num_paths; ++i) {
                auto name = xgidx->path_name(i);
                path_length[name] = xgidx->path_length(name);
            }
            hdr = hts_string_header(sam_header, path_length, rg_sample);
            if ((sam_out = sam_open("-", out_mode)) == 0) {
                cerr << "[vg map] failed to open stdout for writing HTS output" << endl;
                exit(1);
            } else {
                // write the header
                if (sam_hdr_write(sam_out, hdr) != 0) {
                    cerr << "[vg map] error: failed to write the SAM header" << endl;
                }
            }
        }
    };

    // TODO: Refactor the surjection code out of surject_main and intto somewhere where we can just use it here!

    auto surject_alignments = [&hdr, &sam_header, &mapper, &rg_sample, &setup_sam_header, &path_names, &sam_out, &xgidx] (const vector<Alignment>& alns1, const vector<Alignment>& alns2) {
        
        if (alns1.empty()) return;
        setup_sam_header();
        vector<tuple<string, int64_t, bool, Alignment> > surjects1, surjects2;
        int tid = omp_get_thread_num();
        for (auto& aln : alns1) {
            // Surject each alignment of the first read in the pair
            string path_name;
            int64_t path_pos = -1;
            bool path_reverse = false;
            
            auto surj = mapper[tid]->surject_alignment(aln, path_names, path_name, path_pos, path_reverse);
            surjects1.push_back(make_tuple(path_name, path_pos, path_reverse, surj));
            
            // hack: if we haven't established the header, we look at the reads to guess which read groups to put in it
            if (!hdr && !surj.read_group().empty() && !surj.sample_name().empty()) {
#pragma omp critical (hts_header)
                rg_sample[surj.read_group()] = surj.sample_name();
            }
        }
        
        for (auto& aln : alns2) {
            // Surject each alignment of the second read in the pair, if any
            string path_name;
            int64_t path_pos = -1;
            bool path_reverse = false;
            
            auto surj = mapper[tid]->surject_alignment(aln, path_names, path_name, path_pos, path_reverse);
            surjects2.push_back(make_tuple(path_name, path_pos, path_reverse, surj));
            
            // Don't try and populate the header; it should have happened already
        }
        
        if (surjects2.empty()) {
            // Write out surjected single-end reads
        
            for (auto& s : surjects1) {
                auto& path_name = get<0>(s);
                auto& path_pos = get<1>(s);
                auto& path_reverse = get<2>(s);
                auto& surj = get<3>(s);
                
                size_t path_len = 0;
                if (path_name != "") {
                    path_len = xgidx->path_length(path_name);
                }
                string cigar = cigar_against_path(surj, path_reverse, path_pos, path_len, 0);
                bam1_t* b = alignment_to_bam(sam_header,
                                             surj,
                                             path_name,
                                             path_pos,
                                             path_reverse,
                                             cigar);
                int r = 0;
#pragma omp critical (cout)
                r = sam_write1(sam_out, hdr, b);
                if (r == 0) { cerr << "[vg map] error: writing to stdout failed" << endl; exit(1); }
                bam_destroy1(b);
            }
        } else {
            // Write out surjected paired-end reads
            
            // Paired-end reads come in corresponding pairs, allowing duplicate reads.
            assert(surjects1.size() == surjects2.size());
            
            for (size_t i = 0; i < surjects1.size(); i++) {
                // For each corresponding pair
                auto& s1 = surjects1[i];
                auto& s2 = surjects2[i];

                // Unpack each read
                auto& path_name1 = get<0>(s1);
                auto& path_pos1 = get<1>(s1);
                auto& path_reverse1 = get<2>(s1);
                auto& surj1 = get<3>(s1);
                
                auto& path_name2 = get<0>(s2);
                auto& path_pos2 = get<1>(s2);
                auto& path_reverse2 = get<2>(s2);
                auto& surj2 = get<3>(s2);
                
                // Compute CIGARs
                size_t path_len1, path_len2;
                if (path_name1 != "") {
                    path_len1 = xgidx->path_length(path_name1);
                }
                if (path_name2 != "") {
                    path_len2 = xgidx->path_length(path_name2);
                }
                string cigar1 = cigar_against_path(surj1, path_reverse1, path_pos1, path_len1, 0);
                string cigar2 = cigar_against_path(surj2, path_reverse2, path_pos2, path_len2, 0);
                
                // TODO: compute template length based on
                // pair distance and alignment content.
                int template_length = 0;
                
                // Make BAM records
                bam1_t* b1 = alignment_to_bam(sam_header,
                                              surj1,
                                              path_name1,
                                              path_pos1,
                                              path_reverse1,
                                              cigar1,
                                              path_name2,
                                              path_pos2,
                                              template_length);
                bam1_t* b2 = alignment_to_bam(sam_header,
                                              surj2,
                                              path_name2,
                                              path_pos2,
                                              path_reverse2,
                                              cigar2,
                                              path_name1,
                                              path_pos1,
                                              template_length);
                
                // Write the records
                int r = 0;
#pragma omp critical (cout)
                r = sam_write1(sam_out, hdr, b1);
                if (r == 0) { cerr << "[vg map] error: writing to stdout failed" << endl; exit(1); }
                bam_destroy1(b1);
                r = 0;
#pragma omp critical (cout)
                r = sam_write1(sam_out, hdr, b2);
                if (r == 0) { cerr << "[vg map] error: writing to stdout failed" << endl; exit(1); }
                bam_destroy1(b2);
            }
            
            
        }
    };

    auto write_json = [](const vector<Alignment>& alns) {
        for(auto& alignment : alns) {
            string json = pb2json(alignment);
            cout << json << "\n";
        }
    };

    auto write_refpos = [](const vector<Alignment>& alns) {
        for(auto& alignment : alns) {
            Position refpos;
            if (alignment.refpos_size()) {
                refpos = alignment.refpos(0);
            }
            cout << alignment.name() << "\t"
            << refpos.name() << "\t"
            << refpos.offset() << "\t"
            << alignment.mapping_quality() << "\t"
            << alignment.score() << "\n";
        }
    };

    // We have one function to dump alignments into
    // Make sure to flush the buffer at the end of the program!
    auto output_alignments = [&output_buffer,
                              &output_json,
                              &surject_type,
                              &surject_alignments,
                              &buffer_size,
                              &refpos_table,
                              &write_json,
                              &write_refpos](const vector<Alignment>& alns1, const vector<Alignment>& alns2) {
        if (output_json) {
            // If we want to convert to JSON, convert them all to JSON and dump them to cout.
#pragma omp critical (cout)
            {
                write_json(alns1);
                write_json(alns2);
            }
        } else if (refpos_table) {
            // keep multi alignments ordered appropriately
#pragma omp critical (cout)
            {
                write_refpos(alns1);
                write_refpos(alns2);
            }
        } else if (!surject_type.empty()) {
            // surject
            surject_alignments(alns1, alns2);
        } else {
            // Otherwise write them through the buffer for our thread
            int tid = omp_get_thread_num();
            auto& output_buf = output_buffer[tid];

            // Copy all the alignments over to the output buffer
            copy(alns1.begin(), alns1.end(), back_inserter(output_buf));
            copy(alns2.begin(), alns2.end(), back_inserter(output_buf));

            stream::write_buffered(cout, output_buf, buffer_size);
        }
    };

    for (int i = 0; i < thread_count; ++i) {
        Mapper* m = nullptr;
        if(xgidx && gcsa && lcp) {
            // We have the xg and GCSA indexes, so use them
            m = new Mapper(xgidx, gcsa, lcp, haplo_score_provider);
        } else {
            // Can't continue with null
            throw runtime_error("Need XG, GCSA, and LCP to create a Mapper");
        }
        m->hit_max = hit_max;
        m->hit_limit = max(max_multimaps, extra_multimaps);
        m->max_multimaps = max_multimaps;
        m->min_multimaps = max(min_multimaps, max_multimaps);
        m->band_multimaps = band_multimaps;
        m->min_banded_mq = min_banded_mq;
        m->maybe_mq_threshold = maybe_mq_threshold;
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
        m->max_target_factor = max_target_factor;
        m->set_alignment_scores(match, mismatch, gap_open, gap_extend, full_length_bonus, haplotype_consistency_exponent);
        m->strip_bonuses = strip_bonuses;
        m->adjust_alignments_for_base_quality = qual_adjust_alignments;
        m->extra_multimaps = extra_multimaps;
        m->mapping_quality_method = mapping_quality_method;
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
        m->max_band_jump = max_band_jump > -1 ? max_band_jump : band_width;
        m->identity_weight = identity_weight;
        m->assume_acyclic = acyclic_graph;
        m->context_depth = 3; // for surjection
        m->patch_alignments = patch_alignments;
        mapper[i] = m;
    }

    if (!seq.empty()) {
        int tid = omp_get_thread_num();

        Alignment unaligned;
        unaligned.set_sequence(seq);

        if (!qual.empty()) {
            unaligned.set_quality(qual);
        }

        vector<Alignment> alignments = mapper[tid]->align_multi(unaligned, kmer_size, kmer_stride, max_mem_length, band_width);
        if(alignments.size() == 0) {
            // If we didn't have any alignments, report the unaligned alignment
            alignments.push_back(unaligned);
        }


        for(auto& alignment : alignments) {
            if (!sample_name.empty()) alignment.set_sample_name(sample_name);
            if (!read_group.empty()) alignment.set_read_group(read_group);
            if (!seq_name.empty()) alignment.set_name(seq_name);
        }

        // Output the alignments in JSON or protobuf as appropriate.
        output_alignments(alignments, empty_alns);
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

                    vector<Alignment> alignments = mapper[tid]->align_multi(unaligned, kmer_size, kmer_stride, max_mem_length, band_width);

                    for(auto& alignment : alignments) {
                        // Set the alignment metadata
                        if (!sample_name.empty()) alignment.set_sample_name(sample_name);
                        if (!read_group.empty()) alignment.set_read_group(read_group);
                    }


                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alignments, empty_alns);
                }
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
                vector<Alignment> alignments = mapper[tid]->align_multi(unaligned, kmer_size, kmer_stride, max_mem_length, band_width);
                for(auto& alignment : alignments) {
                    // Set the alignment metadata
                    if (!sample_name.empty()) alignment.set_sample_name(sample_name);
                    if (!read_group.empty()) alignment.set_read_group(read_group);
                }
                // Output the alignments in JSON or protobuf as appropriate.
                output_alignments(alignments, empty_alns);
            }
        };
#pragma omp parallel for
        for (size_t i = 0; i < ref.index->sequenceNames.size(); ++i) {
            auto& name = ref.index->sequenceNames[i];
            string seq = nonATGCNtoN(ref.getSequence(name));
            align_seq(name, seq);
        }
    }

    if (!hts_file.empty()) {
        function<void(Alignment&)> lambda =
            [&mapper,
             &output_alignments,
             &keep_secondary,
             &kmer_size,
             &kmer_stride,
             &max_mem_length,
             &band_width,
             &empty_alns]
                (Alignment& alignment) {

                    if(alignment.is_secondary() && !keep_secondary) {
                        // Skip over secondary alignments in the input; we don't want several output mappings for each input *mapping*.
                        return;
                    }

                    int tid = omp_get_thread_num();
                    vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, max_mem_length, band_width);

                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alignments, empty_alns);
                };
        // run
        hts_for_each_parallel(hts_file, lambda);
    }

    if (!fastq1.empty()) {
        if (interleaved_input) {
            // paired interleaved
            auto output_func = [&output_alignments,
                                &compare_gam,
                                &print_fragment_model]
                (Alignment& aln1,
                 Alignment& aln2,
                 pair<vector<Alignment>, vector<Alignment>>& alnp) {
                if (!print_fragment_model) {
                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alnp.first, alnp.second);
                }
            };
            function<void(Alignment&,Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &keep_secondary,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width,
                 &pair_window,
                 &top_pairs_only,
                 &print_fragment_model,
                 &output_func](Alignment& aln1, Alignment& aln2) {
                auto our_mapper = mapper[omp_get_thread_num()];
                bool queued_resolve_later = false;
                auto alnp = our_mapper->align_paired_multi(aln1, aln2, queued_resolve_later, max_mem_length, top_pairs_only, false);
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
                                                                       true);
                            output_func(p.first, p.second, alnp);
                        }
                        our_mapper->imperfect_pairs_to_retry.clear();
                    }
                }
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
                                                               true);
                    output_func(p.first, p.second, alnp);
                }
                our_mapper->imperfect_pairs_to_retry.clear();
            }
        } else if (fastq2.empty()) {
            // single
            function<void(Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width,
                 &empty_alns]
                    (Alignment& alignment) {

                        int tid = omp_get_thread_num();
                        vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, max_mem_length, band_width);
                        //cerr << "This is just before output_alignments" << alignment.DebugString() << endl;
                        output_alignments(alignments, empty_alns);
                    };
            fastq_unpaired_for_each_parallel(fastq1, lambda);
        } else {
            // paired two-file
            auto output_func = [&output_alignments,
                                &print_fragment_model]
                (Alignment& aln1,
                 Alignment& aln2,
                 pair<vector<Alignment>, vector<Alignment>>& alnp) {
                // Make sure we have unaligned "alignments" for things that don't align.
                // Output the alignments in JSON or protobuf as appropriate.
                if (!print_fragment_model) {
                    output_alignments(alnp.first, alnp.second);
                }
            };
            function<void(Alignment&,Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &keep_secondary,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width,
                 &pair_window,
                 &top_pairs_only,
                 &print_fragment_model,
                 &output_func](Alignment& aln1, Alignment& aln2) {
                auto our_mapper = mapper[omp_get_thread_num()];
                bool queued_resolve_later = false;
                auto alnp = our_mapper->align_paired_multi(aln1, aln2, queued_resolve_later, max_mem_length, top_pairs_only, false);
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
                                                                       true);
                            output_func(p.first, p.second, alnp);
                        }
                        our_mapper->imperfect_pairs_to_retry.clear();
                    }
                }
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
                                                               true);
                    output_func(p.first, p.second, alnp);
                }
                our_mapper->imperfect_pairs_to_retry.clear();
            }
        }
    }

    if (!gam_input.empty()) {
        ifstream gam_in(gam_input);
        if (interleaved_input) {
            auto output_func = [&output_alignments,
                                &compare_gam,
                                &print_fragment_model]
                (Alignment& aln1,
                 Alignment& aln2,
                 pair<vector<Alignment>, vector<Alignment>>& alnp) {
                if (print_fragment_model) {
                    // do nothing
                } else {
                    // Output the alignments in JSON or protobuf as appropriate.
                    if (compare_gam) {
                        alnp.first.front().set_correct(overlap(aln1.path(), alnp.first.front().path()));
                        alnp.second.front().set_correct(overlap(aln2.path(), alnp.second.front().path()));
                    }
                    output_alignments(alnp.first, alnp.second);
                }
            };
            function<void(Alignment&,Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &keep_secondary,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width,
                 &compare_gam,
                 &pair_window,
                 &top_pairs_only,
                 &print_fragment_model,
                 &output_func](Alignment& aln1, Alignment& aln2) {
                auto our_mapper = mapper[omp_get_thread_num()];
                bool queued_resolve_later = false;
                auto alnp = our_mapper->align_paired_multi(aln1, aln2, queued_resolve_later, max_mem_length, top_pairs_only, false);
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
                                                                       true);
                            output_func(p.first, p.second, alnp);
                        }
                        our_mapper->imperfect_pairs_to_retry.clear();
                    }
                }
            };
            stream::for_each_interleaved_pair_parallel(gam_in, lambda);
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
                                                               true);
                    output_func(p.first, p.second, alnp);
                }
                our_mapper->imperfect_pairs_to_retry.clear();
            }
        } else {
            function<void(Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &keep_secondary,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width,
                 &compare_gam,
                 &empty_alns]
                (Alignment& alignment) {
                int tid = omp_get_thread_num();
                std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
                vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, max_mem_length, band_width);
                std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds = end-start;
                // Output the alignments in JSON or protobuf as appropriate.
                if (compare_gam) {
                    alignments.front().set_correct(overlap(alignment.path(), alignments.front().path()));
                }
                output_alignments(alignments, empty_alns);
            };
            stream::for_each_parallel(gam_in, lambda);
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

    // clean up
    for (int i = 0; i < thread_count; ++i) {
        delete mapper[i];
        auto& output_buf = output_buffer[i];
        if (!output_json && !refpos_table && surject_type.empty()) {
            stream::write_buffered(cout, output_buf, 0);
        }
    }

    // special cleanup for htslib outputs
    if (!surject_type.empty()) {
        if (hdr != nullptr) bam_hdr_destroy(hdr);
        sam_close(sam_out);
        cout.flush();
    }
    
    if (haplo_score_provider) {
        delete haplo_score_provider;
        haplo_score_provider = nullptr;
    }
    if (gbwt) {
        delete gbwt;
        gbwt = nullptr;
    }
    if (lcp) {
        delete lcp;
        lcp = nullptr;
    }
    if(gcsa) {
        delete gcsa;
        gcsa = nullptr;
    }
    if(xgidx) {
        delete xgidx;
        xgidx = nullptr;
    }

    cout.flush();

    return 0;

}

static Subcommand vg_map("map", "MEM-based read alignment", PIPELINE, 3, main_map);

#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "../mapper.hpp"
#include "../stream.hpp"

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
         << "algorithm:" << endl
         << "    -t, --threads N         number of compute threads to use" << endl
         << "    -k, --min-seed INT      minimum seed (MEM) length (set to -1 to estimate given -e) [-1]" << endl
         << "    -c, --hit-max N         ignore MEMs who have >N hits in our index [512]" << endl
         << "    -e, --seed-chance FLOAT set {-k} such that this fraction of {-k} length hits will by chance [0.0001]" << endl
         << "    -Y, --max-seed INT      ignore seeds longer than this length [0]" << endl
         << "    -r, --reseed-x FLOAT    look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]" << endl
         << "    -u, --try-up-to INT     attempt to align up to the INT best candidate chains of seeds [512]" << endl
         << "    -l, --try-at-least INT  attempt to align up to the INT best candidate chains of seeds [16]" << endl
         << "    -E, --approx-mq-cap INT weight MQ by suffix tree based estimate when estimate less than INT [60]" << endl
         << "    --id-mq-weight N        scale mapping quality by the alignment score identity to this power [2]" << endl
         << "    -W, --min-chain INT     discard a chain if seeded bases shorter than INT [0]" << endl
         << "    -C, --drop-chain FLOAT  drop chains shorter than FLOAT fraction of the longest overlapping chain [0.5]" << endl
         << "    -n, --mq-overlap FLOAT  scale MQ by count of alignments with this overlap in the query with the primary [0.5]" << endl
         << "    -P, --min-ident FLOAT   accept alignment only if the alignment identity is >= FLOAT [0]" << endl
         << "    -H, --max-target-x N    skip cluster subgraphs with length > N*read_length [100]" << endl
         << "    -v, --mq-method OPT     mapping quality method: 0 - none, 1 - fast approximation, 2 - exact [1]" << endl
         << "    -w, --band-width INT    band width for long read alignment [256]" << endl
         << "    -J, --band-jump INT     the maximum jump we can see between bands (maximum length variant we can detect) [{-w}]" << endl
         << "    -I, --fragment STR      fragment length distribution specification STR=m:μ:σ:o:d [10000:0:0:0:1]" << endl
         << "                            max, mean, stdev, orientation (1=same, 0=flip), direction (1=forward, 0=backward)" << endl
         << "    -U, --fixed-frag-model  don't learn the pair fragment model online, use {-I} without update" << endl
         << "    -p, --print-frag-model  suppress alignment output and print the fragment model on stdout as per {-I} format" << endl
         << "    -F, --frag-calc INT     update the fragment model every INT perfect pairs [10]" << endl
         << "    -S, --fragment-x FLOAT  calculate max fragment size as frag_mean+frag_sd*FLOAT [10]" << endl
         << "    -O, --mate-rescues INT  attempt up to INT mate rescues per pair [64]" << endl
         << "scoring:" << endl
         << "    -q, --match INT         use this match score [1]" << endl
         << "    -z, --mismatch INT      use this mismatch penalty [4]" << endl
         << "    -o, --gap-open INT      use this gap open penalty [6]" << endl
         << "    -y, --gap-extend INT    use this gap extension penalty [1]" << endl
         << "    -L, --full-l-bonus INT  the full-length alignment bonus [5]" << endl
         << "    -m, --include-bonuses   include bonuses in reported scores" << endl
         << "    -A, --qual-adjust       perform base quality adjusted alignments (requires base quality input)" << endl
         << "input:" << endl
         << "    -s, --sequence STR      align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -V, --seq-name STR      name the sequence using this value (for graph modification with new named paths)" << endl
         << "    -T, --reads FILE        take reads (one per line) from FILE, write alignments to stdout" << endl
         << "    -b, --hts-input FILE    align reads from htslib-compatible FILE (BAM/CRAM/SAM) stdin (-), alignments to stdout" << endl
         << "    -G, --gam-input FILE    realign GAM input" << endl
         << "    -f, --fastq FILE        input fastq (possibly compressed), two are allowed, one for each mate" << endl
         << "    -i, --interleaved       fastq is interleaved paired-ended" << endl
         << "    -N, --sample NAME       for --reads input, add this sample" << endl
         << "    -R, --read-group NAME   for --reads input, add this read group" << endl
         << "output:" << endl
         << "    -j, --output-json       output JSON rather than an alignment stream (helpful for debugging)" << endl
         << "    -Z, --buffer-size INT   buffer this many alignments together before outputting in GAM [100]" << endl
         << "    -X, --compare           realign GAM input (-G), writing alignment with \"correct\" field set to overlap with input" << endl
         << "    -K, --keep-secondary    produce alignments for secondary input alignments in addition to primary ones" << endl
         << "    -M, --max-multimaps INT produce up to INT alignments for each read [1]" << endl
         << "    -B, --band-multi INT    consider this many alignments of each band in banded alignment [4]" << endl
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
    string read_file;
    string hts_file;
    bool keep_secondary = false;
    int hit_max = 512;
    int max_multimaps = 1;
    int thread_count = 1;
    bool output_json = false;
    bool debug = false;
    float min_score = 0;
    string sample_name;
    string read_group;
    string fastq1, fastq2;
    bool interleaved_input = false;
    int band_width = 256;
    int band_multimaps = 4;
    int max_band_jump = -1;
    bool always_rescue = false;
    bool top_pairs_only = false;
    int max_mem_length = 0;
    int min_mem_length = -1;
    int min_cluster_length = 0;
    float mem_reseed_factor = 1.5;
    bool mem_chaining = true;
    int max_target_factor = 100;
    int buffer_size = 100;
    int8_t match = 1;
    int8_t mismatch = 4;
    int8_t gap_open = 6;
    int8_t gap_extend = 1;
    int8_t full_length_bonus = 5;
    bool strip_bonuses = true;
    bool qual_adjust_alignments = false;
    int extra_multimaps = 512;
    int min_multimaps = 16;
    int max_mapping_quality = 60;
    int method_code = 1;
    int maybe_mq_threshold = 60;
    double identity_weight = 2;
    string gam_input;
    bool compare_gam = false;
    int fragment_max = 1e4;
    int fragment_size = 0;
    double fragment_mean = 0;
    double fragment_stdev = 0;
    double fragment_sigma = 10;
    bool fragment_orientation = false;
    bool fragment_direction = true;
    bool use_cluster_mq = false;
    float chance_match = 0.05;
    bool use_fast_reseed = true;
    float drop_chain = 0.5;
    float mq_overlap = 0.5;
    int kmer_size = 0; // if we set to positive, we'd revert to the old kmer based mapper
    int kmer_stride = 0;
    int pair_window = 64; // unused
    int mate_rescues = 64;
    bool fixed_fragment_model = false;
    bool print_fragment_model = false;
    int fragment_model_update = 10;

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
                {"interleaved", no_argument, 0, 'i'},
                {"band-width", required_argument, 0, 'w'},
                {"band-multi", required_argument, 0, 'B'},
                {"band-jump", required_argument, 0, 'J'},
                {"min-ident", required_argument, 0, 'P'},
                {"debug", no_argument, 0, 'D'},
                {"min-seed", required_argument, 0, 'k'},
                {"max-seed", required_argument, 0, 'Y'},
                {"reseed-x", required_argument, 0, 'r'},
                {"min-chain", required_argument, 0, 'W'},
                {"fast-reseed", no_argument, 0, '6'},
                {"id-clustering", no_argument, 0, 'a'},
                {"max-target-x", required_argument, 0, 'H'},
                {"buffer-size", required_argument, 0, 'Z'},
                {"match", required_argument, 0, 'q'},
                {"mismatch", required_argument, 0, 'z'},
                {"gap-open", required_argument, 0, 'o'},
                {"gap-extend", required_argument, 0, 'y'},
                {"qual-adjust", no_argument, 0, 'A'},
                {"try-up-to", required_argument, 0, 'u'},
                {"compare", no_argument, 0, 'X'},
                {"fragment", required_argument, 0, 'I'},
                {"fragment-x", required_argument, 0, 'S'},
                {"full-l-bonus", required_argument, 0, 'L'},
                {"include-bonuses", no_argument, 0, 'm'},
                {"seed-chance", required_argument, 0, 'e'},
                {"drop-chain", required_argument, 0, 'C'},
                {"mq-overlap", required_argument, 0, 'n'},
                {"try-at-least", required_argument, 0, 'l'},
                {"mq-method", required_argument, 0, 'v'},
                {"mq-max", required_argument, 0, 'Q'},
                {"mate-rescues", required_argument, 0, 'O'},
                {"approx-mq-cap", required_argument, 0, 'E'},
                {"fixed-frag-model", no_argument, 0, 'U'},
                {"print-frag-model", no_argument, 0, 'p'},
                {"frag-calc", required_argument, 0, 'F'},
                {"id-mq-weight", required_argument, 0, '7'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:J:Q:d:x:g:T:N:R:c:M:t:G:jb:Kf:iw:P:Dk:Y:r:W:6aH:Z:q:z:o:y:Au:B:I:S:l:e:C:v:V:O:L:n:E:X:UpF:m7:",
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
            maybe_mq_threshold = atoi(optarg);
            break;

        case 'L':
            full_length_bonus = atoi(optarg);
            break;
        
        case 'm':
            strip_bonuses = false;
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

        case 'Z':
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

        case 'v':
            method_code = atoi(optarg);
            break;

        case 'X':
            compare_gam = true;
            output_json = true;
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

        case 'S':
            fragment_sigma = atof(optarg);
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

        case 'F':
            fragment_model_update = atoi(optarg);
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_map(argv);
            exit(1);
            break;


            default:
                abort ();
        }
    }

    if (seq.empty() && read_file.empty() && hts_file.empty() && fastq1.empty() && gam_input.empty()) {
        cerr << "error:[vg map] A sequence or read file is required when mapping." << endl;
        return 1;
    }

    if (!qual.empty() && (seq.length() != qual.length())) {
        cerr << "error:[vg map] Sequence and base quality string must be the same length." << endl;
        return 1;
    }

    if (qual_adjust_alignments && ((fastq1.empty() && hts_file.empty() && qual.empty()) // must have some quality input
                                   || (!seq.empty() && qual.empty())                    // can't provide sequence without quality
                                   || !read_file.empty()))                              // can't provide sequence list without qualities
    {
        cerr << "error:[vg map] Quality adjusted alignments require base quality scores for all sequences." << endl;
        return 1;
    }
    // note: still possible that hts file types don't have quality, but have to check the file to know

    MappingQualityMethod mapping_quality_method;
    if (method_code == 0) {
        mapping_quality_method = None;
    }
    else if (method_code == 1) {
        mapping_quality_method = Approx;
    }
    else if (method_code == 2) {
        mapping_quality_method = Exact;
    }
    else {
        cerr << "error:[vg map] unrecognized mapping quality method command line arg '" << method_code << "'" << endl;
        return 1;
    }


    string file_name;
    if (optind < argc) {
        file_name = get_input_file_name(optind, argc, argv);
    }

    if (gcsa_name.empty() && !file_name.empty()) {
        gcsa_name = file_name + gcsa::GCSA::EXTENSION;
    }

    if (xg_name.empty() && !file_name.empty()) {
        xg_name = file_name + ".xg";
    }

    if (!db_name.empty()) {
        xg_name = db_name + ".xg";
        gcsa_name = db_name + gcsa::GCSA::EXTENSION;
    }

    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Configure its temp directory to the system temp directory
    gcsa::TempFile::setDirectory(find_temp_dir());

    // Load up our indexes.
    xg::XG* xindex = nullptr;
    gcsa::GCSA* gcsa = nullptr;
    gcsa::LCPArray* lcp = nullptr;

    // We try opening the file, and then see if it worked
    ifstream xg_stream(xg_name);

    if(xg_stream) {
        // We have an xg index!
        if(debug) {
            cerr << "Loading xg index " << xg_name << "..." << endl;
        }
        xindex = new xg::XG(xg_stream);
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

    thread_count = get_thread_count();

    vector<Mapper*> mapper;
    mapper.resize(thread_count);
    vector<vector<Alignment> > output_buffer;
    output_buffer.resize(thread_count);

    // We have one function to dump alignments into
    // Make sure to flush the buffer at the end of the program!
    auto output_alignments = [&output_buffer, &output_json, &buffer_size](vector<Alignment>& alignments) {
        // for(auto& alignment : alignments){
        //     cerr << "This is in output_alignments" << alignment.DebugString() << endl;
        // }

        if (output_json) {
            // If we want to convert to JSON, convert them all to JSON and dump them to cout.
            for(auto& alignment : alignments) {
                string json = pb2json(alignment);
#pragma omp critical (cout)
                cout << json << "\n";
            }
        } else {
            // Otherwise write them through the buffer for our thread
            int tid = omp_get_thread_num();
            auto& output_buf = output_buffer[tid];

            // Copy all the alignments over to the output buffer
            copy(alignments.begin(), alignments.end(), back_inserter(output_buf));

            stream::write_buffered(cout, output_buf, buffer_size);
        }
    };

    for (int i = 0; i < thread_count; ++i) {
        Mapper* m = nullptr;
        if(xindex && gcsa && lcp) {
            // We have the xg and GCSA indexes, so use them
            m = new Mapper(xindex, gcsa, lcp);
        } else {
            // Can't continue with null
            throw runtime_error("Need XG, GCSA, and LCP to create a Mapper");
        }
        m->hit_max = hit_max;
        m->max_multimaps = max_multimaps;
        m->min_multimaps = min_multimaps;
        m->band_multimaps = band_multimaps;
        m->maybe_mq_threshold = maybe_mq_threshold;
        m->debug = debug;
        m->min_identity = min_score;
        m->drop_chain = drop_chain;
        m->mq_overlap = mq_overlap;
        m->min_mem_length = (min_mem_length > 0 ? min_mem_length
                             : m->random_match_length(chance_match));
        m->mem_reseed_length = round(mem_reseed_factor * m->min_mem_length);
        m->min_cluster_length = min_cluster_length;
        if (debug && i == 0) {
            cerr << "[vg map] : min_mem_length = " << m->min_mem_length
                 << ", mem_reseed_length = " << m->mem_reseed_length
                 << ", min_cluster_length = " << m->min_cluster_length << endl;
        }
        m->fast_reseed = use_fast_reseed;
        m->mem_chaining = mem_chaining;
        m->max_target_factor = max_target_factor;
        m->set_alignment_scores(match, mismatch, gap_open, gap_extend, full_length_bonus);
        m->strip_bonuses = strip_bonuses;
        m->adjust_alignments_for_base_quality = qual_adjust_alignments;
        m->extra_multimaps = extra_multimaps;
        m->mapping_quality_method = mapping_quality_method;
        m->always_rescue = always_rescue;
        m->fixed_fragment_model = fixed_fragment_model;
        m->fragment_max = fragment_max;
        m->fragment_sigma = fragment_sigma;
        if (fragment_mean) {
            m->fragment_size = fragment_size;
            m->cached_fragment_length_mean = fragment_mean;
            m->cached_fragment_length_stdev = fragment_stdev;
            m->cached_fragment_orientation = fragment_orientation;
            m->cached_fragment_direction = fragment_direction;
        }
        m->fragment_model_update_interval = fragment_model_update;
        m->max_mapping_quality = max_mapping_quality;
        m->use_cluster_mq = use_cluster_mq;
        m->mate_rescues = mate_rescues;
        m->max_band_jump = max_band_jump > -1 ? max_band_jump : band_width;
        m->identity_weight = identity_weight;
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
        output_alignments(alignments);
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
                    if(alignments.empty()) {
                        alignments.push_back(unaligned);
                    }

                    for(auto& alignment : alignments) {
                        // Set the alignment metadata
                        if (!sample_name.empty()) alignment.set_sample_name(sample_name);
                        if (!read_group.empty()) alignment.set_read_group(read_group);
                    }


                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alignments);
                }
            }
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
             &band_width]
                (Alignment& alignment) {

                    if(alignment.is_secondary() && !keep_secondary) {
                        // Skip over secondary alignments in the input; we don't want several output mappings for each input *mapping*.
                        return;
                    }

                    int tid = omp_get_thread_num();
                    vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, max_mem_length, band_width);
                    if(alignments.empty()) {
                        alignments.push_back(alignment);
                    }

                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alignments);
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
                    output_alignments(alnp.first);
                    output_alignments(alnp.second);
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
                    if (our_mapper->fragment_size != 0
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
                our_mapper->fragment_size = fragment_max;
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
                 &band_width]
                    (Alignment& alignment) {

                        int tid = omp_get_thread_num();
                        vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, max_mem_length, band_width);

                        if(alignments.empty()) {
                            // Make sure we have a "no alignment" alignment
                            alignments.push_back(alignment);
                        }

                        //cerr << "This is just before output_alignments" << alignment.DebugString() << endl;
                        output_alignments(alignments);
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
                    output_alignments(alnp.first);
                    output_alignments(alnp.second);
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
                    if (our_mapper->fragment_size != 0
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
                our_mapper->fragment_size = fragment_max;
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
                    output_alignments(alnp.first);
                    output_alignments(alnp.second);
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
                    if (our_mapper->fragment_size != 0
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
                our_mapper->fragment_size = fragment_max;
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
                 &compare_gam]
                (Alignment& alignment) {
                int tid = omp_get_thread_num();
                std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
                vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, max_mem_length, band_width);
                std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds = end-start;
                if(alignments.empty()) {
                    alignments.push_back(alignment);
                }
                // Output the alignments in JSON or protobuf as appropriate.
                if (compare_gam) {
                    alignments.front().set_correct(overlap(alignment.path(), alignments.front().path()));
                }
                output_alignments(alignments);
            };
            stream::for_each_parallel(gam_in, lambda);
        }
        gam_in.close();
    }

    if (print_fragment_model) {
        if (mapper[0]->fragment_size) {
            // we've calculated our fragment size, so print it and bail out
            cout << mapper[0]->fragment_model_str() << endl;
        } else {
            cerr << "[vg map] Error: could not calculate fragment model" << endl;
        }
    }

    // clean up
    for (int i = 0; i < thread_count; ++i) {
        delete mapper[i];
        auto& output_buf = output_buffer[i];
        if (!output_json) {
            stream::write_buffered(cout, output_buf, 0);
        }
    }

    if(gcsa) {
        delete gcsa;
        gcsa = nullptr;
    }
    if(xindex) {
        delete xindex;
        xindex = nullptr;
    }

    cout.flush();

    return 0;

}

static Subcommand vg_msga("map", "MEM-based read alignment", main_map);

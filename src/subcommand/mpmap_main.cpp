/**
 * \file explode_main.cpp: break a graph into connected components
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

#define debug_mpmap

using namespace std;
using namespace vg;
using namespace vg::subcommand;


void help_mpmap(char** argv) {
    cerr
    << "usage: " << argv[0] << " mpmap [options] -x index.xg -g index.gcsa -f reads1.fq [-f reads2.fq] > aln.gamp" << endl
    << "Multipath align reads to a graph." << endl
    << endl
    << "options:" << endl
    << "graph/index:" << endl
    << "    -x, --xg-name FILE        use this xg index" << endl
    << "    -g, --gcsa-name FILE      use this GCSA2 index (requires both .gcsa and .gcsa.lcp)" << endl
    << "input:" << endl
    << "    -f, --fastq FILE          input FASTQ (possibly compressed), can be given twice for paired ends" << endl
    << "    -i, --interleaved         FASTQ is interleaved with paired ends" << endl
    << "    -b, --hts-input FILE      align reads from htslib-compatible FILE (BAM/CRAM/SAM) stdin (-)" << endl
    << "    -G, --gam-input FILE      realign GAM input" << endl
    << "algorithm:" << endl
    << "    -s, --snarls FILE         align to alternate paths in these snarls" << endl
    << "    -M, --max-multimaps FILE  compute at most this many mappings" << endl
    << "scoring:" << endl
    << "    -q, --match INT         use this match score [1]" << endl
    << "    -z, --mismatch INT      use this mismatch penalty [4]" << endl
    << "    -o, --gap-open INT      use this gap open penalty [6]" << endl
    << "    -y, --gap-extend INT    use this gap extension penalty [1]" << endl
    << "    -L, --full-l-bonus INT  the full-length alignment bonus [5]" << endl
    << "computational parameters:" << endl
    << "    -t, --threads N           number of compute threads to use" << endl
    << "    -Z, --buffer-size INT     buffer this many alignments together before outputting in .gamp [100]" << endl;
    
}

int main_mpmap(int argc, char** argv) {

    if (argc == 2) {
        help_mpmap(argv);
        return 1;
    }

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
    int max_num_mappings = 2;
    int buffer_size = 100;
    
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
            {"snarls", required_argument, 0, 's'},
            {"max-multimaps", required_argument, 0, 'M'},
            {"match", required_argument, 0, 'q'},
            {"mismatch", required_argument, 0, 'z'},
            {"gap-open", required_argument, 0, 'o'},
            {"gap-extend", required_argument, 0, 'y'},
            {"full-l-bonus", required_argument, 0, 'L'},
            {"threads", required_argument, 0, 't'},
            {"buffer-size", required_argument, 0, 'Z'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:g:b:f:iG:s:M:q:z:o:y:L:t:Z:",
                long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'x':
                xg_name = optarg;
                break;
                
            case 'g':
                gcsa_name = optarg;
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
                if (snarls_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide HTS file with -b" << endl;
                    exit(1);
                }
                break;
                
            case 'G':
                gam_file_name = optarg;
                if (snarls_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide GAM file with -G" << endl;
                    exit(1);
                }
                break;
                
            case 's':
                snarls_name = optarg;
                if (snarls_name.empty()) {
                    cerr << "error:[vg mpmap] Must provide snarl file with -s" << endl;
                    exit(1);
                }
                break;
                
            case 'M':
                max_num_mappings = atoi(optarg);
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
                
            case 'L':
                full_length_bonus = atoi(optarg);
                break;
                
            case 't':
                omp_set_num_threads(atoi(optarg));
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

    if (xg_name.empty()) {
        cerr << "error:[vg mpmap] Multipath mapping requires an XG index, must provide XG file" << endl;
        exit(1);
    }
    
    if (gcsa_name.empty()) {
        cerr << "error:[vg mpmap] Multipath mapping requires a GCSA2 index, must provide GCSA2 file" << endl;
        exit(1);
    }
    
    if (max_num_mappings <= 0) {
        cerr << "error:[vg mpmap] Maximum number of multimappings set to " << max_num_mappings << ", set to a positive integer" << endl;
        exit(1);
    }
    
    if (match_score > std::numeric_limits<int8_t>::max() || mismatch_score > std::numeric_limits<int8_t>::max()
        || gap_open_score > std::numeric_limits<int8_t>::max() || gap_extension_score > std::numeric_limits<int8_t>::max()
        || full_length_bonus > std::numeric_limits<int8_t>::max()) {
        cerr << "error:[vg mpmap] All alignment scoring parameters must be less than " << (int) std::numeric_limits<int8_t>::max() << endl;
        exit(1);
    }
    
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
    multipath_mapper.set_alignment_scores(match_score, mismatch_score, gap_open_score,
                                          gap_extension_score, full_length_bonus);
    
    vector<vector<MultipathAlignment> > output_buffer;
    output_buffer.resize(omp_get_num_threads());
    
    auto output_alignments = [&](list<MultipathAlignment>& mp_alns) {
        auto& output_buf = output_buffer[omp_get_thread_num()];
        
        // Copy all the alignments over to the output buffer
        copy(mp_alns.begin(), mp_alns.end(), back_inserter(output_buf));
        
        stream::write_buffered(cout, output_buf, buffer_size);
    };
    
#ifdef debug_mpmap
    cerr << "[vg mpmap] created all in memory objects, beginning mapping" << endl;
#endif
    
    if (!fastq_name_1.empty()) {
        if (fastq_name_2.empty()) {
            fastq_unpaired_for_each_parallel(fastq_name_1,
                                             [&](Alignment& alignment) {
                                                 list<MultipathAlignment> mp_alns;
                                                 multipath_mapper.multipath_map(alignment, mp_alns, max_num_mappings);
                                                 output_alignments(mp_alns);
                                             });
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
        // TODO
    }
    
    delete snarl_manager;
    
    return 0;
}

// Register subcommand
static Subcommand vg_mpmap("mpmap", "multipath alignments of reads to a graph", main_mpmap);



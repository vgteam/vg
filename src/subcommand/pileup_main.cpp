/** \file pileup_main.cpp
 *
 * Defines the "vg pileup" subcommand, which makes a pileup from a GAM.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../pileup.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_pileup(char** argv) {
    cerr << "usage: " << argv[0] << " pileup [options] <graph.vg> <alignment.gam> > out.vgpu" << endl
         << "Calculate pileup for each position in graph and output in VG Pileup format (list of protobuf NodePileups)." << endl
         << endl
         << "options:" << endl
         << "    -j, --json              output in JSON" << endl
         << "    -q, --min-quality N     ignore bases with PHRED quality < N (default=10)" << endl
         << "    -m, --max-mismatches N  ignore bases with > N mismatches within window centered on read (default=1)" << endl
         << "    -w, --window-size N     size of window to apply -m option (default=0)" << endl
         << "    -d, --max-depth N       maximum depth pileup to create (further maps ignored) (default=1000)" << endl
         << "    -M, --ignore-mapq       do not combine mapping qualities with base qualities" << endl
         << "    -p, --progress          show progress" << endl
         << "    -t, --threads N         number of threads to use" << endl
         << "    -v, --verbose           print stats on bases filtered" << endl;
}

int main_pileup(int argc, char** argv) {

    if (argc <= 3) {
        help_pileup(argv);
        return 1;
    }

    bool output_json = false;
    bool show_progress = false;
    int thread_count = 1;
    int min_quality = 10;
    int max_mismatches = 1;
    int window_size = 0;
    int max_depth = 1000; // used to prevent protobuf messages getting to big
    bool verbose = false;
    bool use_mapq = true;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
            {
                {"json", required_argument, 0, 'j'},
                {"min-quality", required_argument, 0, 'q'},
                {"max-mismatches", required_argument, 0, 'm'},
                {"window-size", required_argument, 0, 'w'},
                {"progress", required_argument, 0, 'p'},
                {"max-depth", required_argument, 0, 'd'},
                {"ignore-mapq", no_argument, 0, 'M'},
                {"threads", required_argument, 0, 't'},
                {"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "jq:m:w:pd:at:v",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'j':
            output_json = true;
            break;
        case 'q':
            min_quality = atoi(optarg);
            break;
        case 'm':
            max_mismatches = atoi(optarg);
            break;
        case 'w':
            window_size = atoi(optarg);
            break;
        case 'd':
            max_depth = atoi(optarg);
            break;
        case 'M':
            use_mapq = false;
            break;
        case 'p':
            show_progress = true;
            break;
        case 't':
            thread_count = atoi(optarg);
            break;
        case 'v':
            verbose = true;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_pileup(argv);
            exit(1);
            break;
        default:
          abort ();
        }
    }
    omp_set_num_threads(thread_count);
    thread_count = get_thread_count();

    // Parse the arguments
    if (optind >= argc) {
        help_pileup(argv);
        return 1;
    }
    string graph_file_name = get_input_file_name(optind, argc, argv);
    if (optind >= argc) {
        help_pileup(argv);
        return 1;
    }
    string alignments_file_name = get_input_file_name(optind, argc, argv);
    
    if (alignments_file_name == "-" && graph_file_name == "-") {
        cerr << "error: graph and alignments can't both be from stdin." << endl;
        exit(1);
    }

    // read the graph
    if (show_progress) {
        cerr << "Reading input graph" << endl;
    }
    VG* graph;
    get_input_file(graph_file_name, [&](istream& in) {
        graph = new VG(in);
    });

    // Make Pileups makers for each thread.
    vector<Pileups> pileups(thread_count, Pileups(graph, min_quality, max_mismatches, window_size, max_depth, use_mapq));
    
    // setup alignment stream
    get_input_file(alignments_file_name, [&](istream& alignment_stream) {
        // compute the pileups.
        if (show_progress) {
            cerr << "Computing pileups" << endl;
        }
        
        function<void(Alignment&)> lambda = [&pileups, &graph](Alignment& aln) {
            int tid = omp_get_thread_num();
            pileups[tid].compute_from_alignment(aln);
        };
        stream::for_each_parallel(alignment_stream, lambda);
    });

    // single-threaded (!) merge
    if (show_progress && pileups.size() > 1) {
        cerr << "Merging pileups" << endl;
    }
    for (int i = 1; i < pileups.size(); ++i) {
        pileups[0].merge(pileups[i]);
    }

    // spit out the pileup
    if (show_progress) {
        cerr << "Writing pileups" << endl;
    }
    if (output_json == false) {
        pileups[0].write(std::cout);
    } else {
        pileups[0].to_json(std::cout);
    }

    delete graph;

    // number of bases filtered
    if (verbose) {
        cerr << "Bases filtered by min. quality: " << pileups[0]._min_quality_count << endl
             << "Bases filtered by max mismatch: " << pileups[0]._max_mismatch_count << endl
             << "Total bases:                    " << pileups[0]._bases_count << endl << endl;
    }

    return 0;
}

// Register subcommand
static Subcommand vg_pileup("pileup", "build a pileup from a set of alignments", main_pileup);


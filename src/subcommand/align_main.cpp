/** \file align_main.cpp
 *
 * Defines the "vg align" subcommand, which aligns a read against an entire
 * graph.
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_align(char** argv) {
    cerr << "usage: " << argv[0] << " align [options] <graph.vg> >alignments.gam" << endl
         << "options:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -Q, --seq-name STR    name the sequence using this value" << endl
         << "    -j, --json            output alignments in JSON format (default GAM)" << endl
         << "    -m, --match N         use this match score (default: 1)" << endl
         << "    -M, --mismatch N      use this mismatch penalty (default: 4)" << endl
         << "    -g, --gap-open N      use this gap open penalty (default: 6)" << endl
         << "    -e, --gap-extend N    use this gap extension penalty (default: 1)" << endl
         << "    -T, --full-l-bonus N  provide this bonus for alignments that are full length (default: 5)" << endl
         << "    -i, --include-bonuses include bonuses in reported scores" << endl
         << "    -b, --banded-global   use the banded global alignment algorithm" << endl
         << "    -p, --pinned          pin the (local) alignment traceback to the optimal edge of the graph" << endl
         << "    -L, --pin-left        pin the first rather than last bases of the graph and sequence" << endl
         << "    -D, --debug           print out score matrices and other debugging info" << endl
         << "options:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -Q, --seq-name STR    name the sequence using this value" << endl
         << "    -r, --reference STR   don't use an input graph--- run SSW alignment between -s and -r" << endl
         << "    -j, --json            output alignments in JSON format (default GAM)" << endl;
}

int main_align(int argc, char** argv) {

    string seq;
    string seq_name;

    if (argc == 2) {
        help_align(argv);
        return 1;
    }

    bool print_cigar = false;
    bool output_json = false;
    int match = 1;
    int mismatch = 4;
    int gap_open = 6;
    int gap_extend = 1;
    int full_length_bonus = 5;
    string ref_seq;
    bool debug = false;
    bool banded_global = false;
    bool pinned_alignment = false;
    bool pin_left = false;
    bool strip_bonuses = true;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"sequence", required_argument, 0, 's'},
            {"seq-name", no_argument, 0, 'Q'},
            {"json", no_argument, 0, 'j'},
            {"match", required_argument, 0, 'm'},
            {"mismatch", required_argument, 0, 'M'},
            {"gap-open", required_argument, 0, 'g'},
            {"gap-extend", required_argument, 0, 'e'},
            {"reference", required_argument, 0, 'r'},
            {"debug", no_argument, 0, 'D'},
            {"banded-global", no_argument, 0, 'b'},
            {"full-l-bonus", required_argument, 0, 'T'},
            {"include-bonuses", no_argument, 0, 'i'},
            {"pinned", no_argument, 0, 'p'},
            {"pin-left", no_argument, 0, 'L'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:jhQ:m:M:g:e:Dr:F:O:bT:ipL",
                long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 's':
            seq = optarg;
            break;

        case 'Q':
            seq_name = optarg;
            break;

        case 'j':
            output_json = true;
            break;

        case 'm':
            match = atoi(optarg);
            break;

        case 'M':
            mismatch = atoi(optarg);
            break;

        case 'g':
            gap_open = atoi(optarg);
            break;

        case 'e':
            gap_extend = atoi(optarg);
            break;

        case 'T':
            full_length_bonus = atoi(optarg);
            break;
            
        case 'i':
            strip_bonuses = false;
            break;

        case 'r':
            ref_seq = optarg;
            break;

        case 'D':
            debug = true;
            break;

        case 'b':
            banded_global = true;
            break;

        case 'p':
            pinned_alignment = true;
            break;

        case 'L':
            pin_left = true;
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_align(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    VG* graph = nullptr;
    if (ref_seq.empty()) {
        // Only look at a filename if we don't have an explicit reference
        // sequence.
        get_input_file(optind, argc, argv, [&](istream& in) {
            graph = new VG(in);
        });
    }
    
    Alignment alignment;
    if (!ref_seq.empty()) {
        SSWAligner ssw = SSWAligner(match, mismatch, gap_open, gap_extend);
        alignment = ssw.align(seq, ref_seq);
    } else {
        Aligner aligner = Aligner(match, mismatch, gap_open, gap_extend, full_length_bonus);
        alignment = graph->align(seq, &aligner, 0, pinned_alignment, pin_left,
            banded_global, 0, max(seq.size(), graph->length()), debug);
        
        if (strip_bonuses) {
            alignment.set_score(aligner.remove_bonuses(alignment, pinned_alignment, pin_left));
        }

    }
    
    if (!seq_name.empty()) {
        alignment.set_name(seq_name);
    }

    if (output_json) {
        cout << pb2json(alignment) << endl;
    } else {
        function<Alignment(uint64_t)> lambda =
            [&alignment] (uint64_t n) {
                return alignment;
            };
        stream::write(cout, 1, lambda);
    }

    if (graph != nullptr) {
        delete graph;
    }

    return 0;

}

// Register subcommand
static Subcommand vg_align("align", "local alignment", main_align);


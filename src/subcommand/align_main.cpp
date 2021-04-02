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

#include "../utility.hpp"
#include "../handle.hpp"
#include "../path.hpp"
#include "../split_strand_graph.hpp"
#include "../dagified_graph.hpp"
#include "../ssw_aligner.hpp"
#include "../aligner.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

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
         << "    --score-matrix FILE   read a 5x5 integer substitution scoring matrix from a file" << endl
         << "    -g, --gap-open N      use this gap open penalty (default: 6)" << endl
         << "    -e, --gap-extend N    use this gap extension penalty (default: 1)" << endl
         << "    -T, --full-l-bonus N  provide this bonus for alignments that are full length (default: 5)" << endl
         << "    -b, --banded-global   use the banded global alignment algorithm" << endl
         << "    -p, --pinned          pin the (local) alignment traceback to the optimal edge of the graph" << endl
         << "    -L, --pin-left        pin the first rather than last bases of the graph and sequence" << endl
         << "    -r, --reference STR   don't use an input graph--- run SSW alignment between -s and -r" << endl
         << "    -D, --debug           print out score matrices and other debugging info" << endl;
}

int main_align(int argc, char** argv) {

    string seq;
    string seq_name;

    if (argc == 2) {
        help_align(argv);
        return 1;
    }

    #define OPT_SCORE_MATRIX 1000
    string matrix_file_name;
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
            {"score-matrix", required_argument, 0, OPT_SCORE_MATRIX},
            {"reference", required_argument, 0, 'r'},
            {"debug", no_argument, 0, 'D'},
            {"banded-global", no_argument, 0, 'b'},
            {"full-l-bonus", required_argument, 0, 'T'},
            {"pinned", no_argument, 0, 'p'},
            {"pin-left", no_argument, 0, 'L'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:jhQ:m:M:g:e:Dr:F:O:bT:pL",
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
            match = parse<int>(optarg);
            break;

        case 'M':
            mismatch = parse<int>(optarg);
            break;

        case 'g':
            gap_open = parse<int>(optarg);
            break;

        case 'e':
            gap_extend = parse<int>(optarg);
            break;

        case 'T':
            full_length_bonus = parse<int>(optarg);
            break;

        case OPT_SCORE_MATRIX:
            matrix_file_name = optarg;
            if (matrix_file_name.empty()) {
                cerr << "error:[vg align] Must provide matrix file with --matrix-file." << endl;
                exit(1);
            }
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

    unique_ptr<PathHandleGraph> graph;
    if (ref_seq.empty()) {
        // Only look at a filename if we don't have an explicit reference
        // sequence.
        string graph_filename = get_input_file_name(optind, argc, argv);
        graph = vg::io::VPKG::load_one<PathHandleGraph>(graph_filename);
    }

    ifstream matrix_stream;
    if (!matrix_file_name.empty()) {
      matrix_stream.open(matrix_file_name);
      if (!matrix_stream) {
          cerr << "error:[vg align] Cannot open scoring matrix file " << matrix_file_name << endl;
          exit(1);
      }
    }
    

    Alignment alignment;
    if (!ref_seq.empty()) {
        if (!matrix_file_name.empty()) {
            cerr << "error:[vg align] Custom scoring matrix not supported in reference sequence mode " << endl;
            exit(1);
        }
        SSWAligner ssw = SSWAligner(match, mismatch, gap_open, gap_extend);
        alignment = ssw.align(seq, ref_seq);
    } else {
        
        // construct a score matrix
        int8_t* score_matrix;
        if (matrix_stream.is_open()) {
            score_matrix = AlignerClient::parse_matrix(matrix_stream);
        }
        else {
            score_matrix = (int8_t*) malloc(sizeof(int8_t) * 16);
            for (size_t i = 0; i < 16; ++i) {
                if (i % 5 == 0) {
                    score_matrix[i] = match;
                }
                else {
                    score_matrix[i] = -mismatch;
                }
            }
        }
        
        // initialize an aligner
        Aligner aligner = Aligner(score_matrix, gap_open, gap_extend, full_length_bonus, vg::default_gc_content);
        
        free(score_matrix);
        
        // put everything on the forward strand
        StrandSplitGraph split(&(*graph));
        
        // dagify it as far as we might ever want
        DagifiedGraph dag(&split, seq.size() + aligner.longest_detectable_gap(seq.size(), seq.size() / 2));
        
        alignment.set_sequence(seq);
        if (pinned_alignment) {
            aligner.align_pinned(alignment, dag, pin_left);
        }
        else if (banded_global) {
            aligner.align_global_banded(alignment, dag, 1, true);
        }
        else {
            aligner.align(alignment, dag, true);
        }
        
        // translate back from the overlays
        translate_oriented_node_ids(*alignment.mutable_path(), [&](vg::id_t node_id) {
            handle_t under = split.get_underlying_handle(dag.get_underlying_handle(dag.get_handle(node_id)));
            return make_pair(graph->get_id(under), graph->get_is_reverse(under));
        });
    }

    if (!seq_name.empty()) {
        alignment.set_name(seq_name);
    }

    if (output_json) {
        cout << pb2json(alignment) << endl;
    } else {
        function<Alignment(size_t)> lambda =
            [&alignment] (size_t n) {
                return alignment;
            };
        vg::io::write(cout, 1, lambda);
        vg::io::finish(cout);
    }

    return 0;

}

// Register subcommand
static Subcommand vg_align("align", "local alignment", main_align);

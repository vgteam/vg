/*
 * Define the "vg remove_duplicate" subcommand, which remove the duplicate PCRs
 *
 *
 */
#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>
#include <cstdint>

#include <subcommand.hpp>

#include "../algorithms/alignment_path_offsets.hpp"
#include "../multipath_alignment.hpp"
#include "../alignment.hpp"
#include "../vg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../stream_sorter.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

vector <pos_t> alignment_position(const Alignment &aln) {
    /// return a vector of pos_t variables which the index of the vector points to the read offset and the value points to the node pos_t
    vector <pos_t> positions;
    for (const Mapping &mapping: aln.path().mapping()) {
        long long node_offset = mapping.position().offset(); // This is the position of the first Edit
        long long node_id = mapping.position().node_id(); // The id of the node we are working with
        bool is_reverse = mapping.position().is_reverse();

        for (const Edit &edit: mapping.edit()) {

            // Match or Mismatch
            if (edit.from_length() == edit.to_length()) {

                for (size_t i = 0; i < edit.to_length(); ++i) {
                    positions.push_back(make_pos_t(node_id, is_reverse, node_offset + i));
                }
                node_offset += edit.from_length();

            } else if (edit.to_length() == 0 && edit.from_length() > edit.to_length()) { // Deletion
                node_offset += edit.from_length();

            } else if (edit.from_length() < edit.to_length()) { // insertion
                node_offset += edit.from_length();
                for (size_t i = 0; i < edit.to_length(); ++i) {
                    positions.push_back(make_pos_t(node_id, is_reverse, node_offset));
                }
            }


        }
    }
    return positions;

}

bool check_duplicate(const Alignment aln1, const Alignment aln2) {
    /// This function check if the two input alignments are duplicate or not, They are duplicate even if have one equal position on one base
    vector <pos_t> first_alignment_pos = alignment_position(aln1);
    vector <pos_t> second_alignment_pos = alignment_position(aln2);
    for (size_t i = 0; i < first_alignment_pos.size(); ++i) {
        if (id(first_alignment_pos[i]) == id(second_alignment_pos[i]) &&
            offset(first_alignment_pos[i]) == offset(second_alignment_pos[i]) &&
            (is_rev(first_alignment_pos[i]) == is_rev(second_alignment_pos[i])))
            return true;

    }
    return false;


}

void print_gam(const Alignment aln){
    cout << aln.sequence() << endl;
}

void help_rmvdup(char **argv) {
// TODO: add whatever option is needed to this list. Change long_option and getopt_long if want to add an option.
// TODO: see what input and output file formats is possible
    cerr << "usage: " << argv[0] << " rmvdup [options] inputfile.gam > output.gam " << endl
         << "Remove duplicate PCRs from the input file. A gam index file (.gam.gai) must exists." << endl
         << "  -p, --progress               Show progress." << endl
         << "    -t, --threads N            number of threads to use" << endl;

}

int main_rmvdup(int argc, char *argv[]) {
    string filename;
    bool show_progress = false;
    int threads = 1;

    int c;
    optind = 2;  // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
                {"help",     no_argument,       0, 'h'},
                {"progress", no_argument,       0, 'p'},
                {"threads",  required_argument, 0, 't'},
                {0,          0,                 0, 0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hpt",
                        long_options, &option_index);

        if (c == -1) break;

        switch (c) {

            case 'p':
                show_progress = true;
                break;

            case 't':
                threads = parse<int>(optarg);
                omp_set_num_threads(threads);
                break;

            case 'h':
            case '?':
                help_rmvdup(argv);
                exit(1);
                break;

            default:
                help_rmvdup(argv);
                abort();
        }
    }

    if (optind >= argc) {
        cerr << "[vg view] error: no filename given" << endl;
        exit(1);
    }

//
//    // This is the function that is going to find the position of each PCR in the graph

//    function<void(Alignment & )> test = [&](const Alignment &aln) {
//        vector <pos_t> temp = alignment_position(aln);
//        cout << aln.sequence() << endl;
//        for (auto it = 0; it < temp.size(); it++)
//            cout << "Read Offset " << it << "Node id " << id(temp[it]) << " Node Offset " << offset(temp[it]) << endl;
//        cout << endl;
//        cout << check_duplicate(aln, aln);
//        cout << temp.size() << temp.

//    };
//
//    string file_name = get_input_file_name(optind, argc, argv);
//    ifstream file_in(file_name);
//    vg::io::for_each_parallel(file_in, test);


//    get_input_file(optind, argc, argv, [&](istream& gam_in) {
//        vg::io::for_each_parallel(gam_in, record_truth);
//        GAMSorter gs(show_progress);
//        unique_ptr<GAMIndex> index;
//        index = unique_ptr<GAMIndex>(new GAMIndex());
//        gs.stream_sort(gam_in, cout, index.get());
//    });

    string sorted_gam_name = get_input_file_name(optind, argc, argv);
    unique_ptr<GAMIndex> gam_index;
    unique_ptr<vg::io::ProtobufIterator<Alignment>> gam_cursor;
    if (!sorted_gam_name.empty()) {
        // Load the GAM index
        gam_index = unique_ptr<GAMIndex>(new GAMIndex());
        get_input_file(sorted_gam_name + ".gai", [&](istream& in) {
            // We get it form the appropriate .gai, which must exist
            gam_index->load(in);
        });
    }
    if (gam_index.get() != nullptr) {
        // Find in sorted GAM

        get_input_file(sorted_gam_name, [&](istream& in) {
            // Make a cursor for input
            vg::io::ProtobufIterator<Alignment> cursor(in);


            // find the alignment for the specific node
            gam_index->find(cursor, 0, 100, print_gam);
        });

    }
    return 0;

}

// Register subcommand
static Subcommand vg_removeduplicate("rmvdup", "Remove duplicate PCRs from the input file", WIDGET, main_rmvdup);
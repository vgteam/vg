/** \file paths_main.cpp
 *
 * Defines the "vg paths" subcommand, which reads paths in the graph.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_paths(char** argv) {
    cerr << "usage: " << argv[0] << " paths [options] <graph.vg>" << endl
         << "options:" << endl
         << "  inspection:" << endl
         << "    -x, --extract         return (as GAM alignments) the stored paths in the graph" << endl
         << "    -L, --list            return (as a list of names, one per line) the path names" << endl
         << "    -X, --list-xg FILE    return path names (as -L) but from given xg file" << endl  
         << "    -s, --as-seqs         write each path as a sequence" << endl;
}

int main_paths(int argc, char** argv) {

    if (argc == 2) {
        help_paths(argv);
        return 1;
    }

    int max_length = 0;
    bool as_seqs = false;
    bool extract = false;
    bool list_paths = false;
    string list_paths_xg;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

        {
            {"extract", no_argument, 0, 'x'},
            {"list", no_argument, 0, 'L'},
            {"list-xg", required_argument, 0, 'X'},
            {"max-length", required_argument, 0, 'l'},
            {"as-seqs", no_argument, 0, 's'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "l:hs:xLX:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'x':
                extract = true;
                break;
                
            case 'L':
                list_paths = true;
                break;

            case 'X':
                list_paths_xg = optarg;
                break;                

            case 'l':
                max_length = atoi(optarg);
                break;

            case 's':
                as_seqs = true;
                break;

            case 'h':
            case '?':
                help_paths(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (!list_paths_xg.empty()) {
        if (optind < argc) {
            cerr << "[vg paths] paths does not accept postional arguments with -X" << endl;
            return 1;
        }
        xg::XG xindex;
        ifstream in(list_paths_xg.c_str());
        xindex.load(in);
        size_t max_path = xindex.max_path_rank();
        for (size_t i = 1; i <= max_path; ++i) {
            cout << xindex.path_name(i) << endl;
        }
        return 0;
    }

    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    if (extract) {
        vector<Alignment> alns = graph->paths_as_alignments();
        write_alignments(cout, alns);
        delete graph;
        return 0;
    }
    
    if (list_paths) {
        graph->paths.for_each_name([&](const string& name) {
            cout << name << endl;
        });
        delete graph;
        return 0;
    }

    if (max_length == 0) {
        cerr << "error:[vg paths] a --max-length is required when generating paths" << endl;
    }

    function<void(size_t,Path&)> paths_to_seqs = [graph](size_t mapping_index, Path& p) {
        string seq = graph->path_sequence(p);
#pragma omp critical(cout)
        cout << seq << endl;
    };

    function<void(size_t,Path&)> paths_to_json = [](size_t mapping_index, Path& p) {
        string json2 = pb2json(p);
#pragma omp critical(cout)
        cout<<json2<<endl;
    };

    function<void(size_t, Path&)>* callback = &paths_to_seqs;
    if (!as_seqs) {
        callback = &paths_to_json;
    }

    delete graph;

    return 0;

}

// Register subcommand
static Subcommand vg_paths("paths", "traverse paths in the graph", main_paths);


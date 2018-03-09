/** \file kmers_main.cpp
 *
 * Defines the "vg kmers" subcommand, which generates the kmers in a graph.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../vg_set.hpp"
#include "../kmer.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_kmers(char** argv) {
    cerr << "usage: " << argv[0] << " kmers [options] <graph1.vg> [graph2.vg ...] >kmers.tsv" << endl
         << "Generates kmers of the graph(s). Output is: kmer id pos" << endl
         << endl
         << "options:" << endl
         << "    -k, --kmer-size N     print kmers of size N in the graph" << endl
         << "    -e, --edge-max N      only consider paths which make edge choices at <= this many points" << endl
         << "    -j, --kmer-stride N   step distance between succesive kmers in paths (default 1)" << endl
         << "    -t, --threads N       number of threads to use" << endl
         << "    -d, --ignore-dups     filter out duplicated kmers in normal output" << endl
         << "    -n, --allow-negs      don't filter out relative negative positions of kmers in normal output" << endl
         << "    -g, --gcsa-out        output a table suitable for input to GCSA2:" << endl
         << "                          kmer, starting position, previous characters," << endl
         << "                          successive characters, successive positions." << endl
         << "                          Forward and reverse strand kmers are reported." << endl
         << "    -B, --gcsa-binary     Write the GCSA graph in binary format." << endl
         << "    -F, --forward-only    When producing GCSA2 output, don't describe the reverse strand" << endl
         << "    -P, --path-only       Only consider kmers if they occur in a path embedded in the graph" << endl
         << "    -H, --head-id N       use the specified ID for the GCSA2 head sentinel node" << endl
         << "    -T, --tail-id N       use the specified ID for the GCSA2 tail sentinel node" << endl
         << "    -p, --progress        show progress" << endl;
}

int main_kmers(int argc, char** argv) {

    if (argc == 2) {
        help_kmers(argv);
        return 1;
    }

    int kmer_size = 0;
    bool path_only = false;
    int edge_max = 0;
    int kmer_stride = 1;
    bool show_progress = false;
    bool gcsa_out = false;
    bool allow_dups = true;
    bool allow_negs = false;
    // for distributed GCSA2 kmer generation
    int64_t head_id = 0;
    int64_t tail_id = 0;
    bool forward_only = false;
    bool gcsa_binary = false;
    bool handle_alg = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

        {
            {"help", no_argument, 0, 'h'},
            {"kmer-size", required_argument, 0, 'k'},
            {"kmer-stride", required_argument, 0, 'j'},
            {"edge-max", required_argument, 0, 'e'},
            {"threads", required_argument, 0, 't'},
            {"gcsa-out", no_argument, 0, 'g'},
            {"ignore-dups", no_argument, 0, 'd'},
            {"allow-negs", no_argument, 0, 'n'},
            {"progress",  no_argument, 0, 'p'},
            {"head-id", required_argument, 0, 'H'},
            {"tail-id", required_argument, 0, 'T'},
            {"forward-only", no_argument, 0, 'F'},
            {"gcsa-binary", no_argument, 0, 'B'},
            {"path-only", no_argument, 0, 'P'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:j:pt:e:gdnH:T:FBP",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'k':
                kmer_size = atoi(optarg);
                break;

            case 'j':
                kmer_stride = atoi(optarg);
                break;

            case 'e':
                edge_max = atoi(optarg);
                break;

            case 't':
                omp_set_num_threads(atoi(optarg));
                break;

            case 'g':
                gcsa_out = true;
                break;

            case 'F':
                forward_only = true;
                break;


            case 'P':
                path_only = true;
                break;

            case 'd':
                allow_dups = false;
                break;

            case 'n':
                allow_negs = true;
                break;

            case 'p':
                show_progress = true;
                break;

            case 'H':
                head_id = atoi(optarg);
                break;

            case 'T':
                tail_id = atoi(optarg);
                break;

            case 'B':
                gcsa_binary = true;
                break;

            case 'h':
            case '?':
                help_kmers(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    vector<string> graph_file_names;
    while (optind < argc) {
        string file_name = get_input_file_name(optind, argc, argv);
        graph_file_names.push_back(file_name);
    }

    VGset graphs(graph_file_names);

    graphs.show_progress = show_progress;

    if (gcsa_out) {
        if (edge_max != 0) {
            // I have been passing this option to vg index -g for months
            // thinking it worked. But it can't work. So we should tell the user
            // they're wrong.
            cerr << "error:[vg kmers] Cannot limit edge crossing (-e) when generating GCSA kmers (-g)."
                << " Use vg mod -p to prune the graph instead." << endl;
            exit(1);
        }
        if (!gcsa_binary) {
            graphs.write_gcsa_kmers_ascii(cout, kmer_size, head_id, tail_id);
        } else {
            size_t limit = ~(size_t)0;
            graphs.write_gcsa_kmers_binary(cout, kmer_size, limit, head_id, tail_id);
        }
    } else {
        //function<void(const kmer_t& kmer)>
        auto lambda = [](const kmer_t& kmer) {
#pragma omp critical (cout)
            cout << kmer << endl;
        };
        graphs.for_each_kmer_parallel(kmer_size, lambda);
    }
    cout.flush();

    return 0;
}

// Register subcommand
static Subcommand vg_kmers("kmers", "enumerate kmers of the graph", main_kmers);


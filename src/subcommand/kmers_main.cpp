/** \file kmers_main.cpp
 *
 * Defines the "vg kmers" subcommand, which generates the kmers in a graph.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>
#include <set>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../vg_set.hpp"
#include "../kmer.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_kmers(char** argv) {
    cerr << "usage: " << argv[0] << " kmers [options] <graph1.vg> [graph2.vg ...] >kmers.tsv" << endl
         << "Generates kmers from both strands of the graph(s). Output is: kmer id pos" << endl
         << endl
         << "general options:" << endl
         << "    -k, --kmer-size N     print kmers of size N in the graph" << endl
         << "    -t, --threads N       number of threads to use" << endl
         << "    -p, --progress        show progress" << endl
         << "gcsa options:" << endl
         << "    -g, --gcsa-out        output a table suitable for input to GCSA2:" << endl
         << "                          kmer, starting position, previous characters," << endl
         << "                          successive characters, successive positions." << endl
         << "    -B, --gcsa-binary     write the GCSA graph in binary format (implies -g)" << endl
         << "    -H, --head-id N       use the specified ID for the GCSA2 head sentinel node" << endl
         << "    -T, --tail-id N       use the specified ID for the GCSA2 tail sentinel node" << endl
         << "minimizer options:" << endl
         << "    -w, --window-size N   output the smallest kmer out of N consecutive kmers" << endl
         << "    -G, --gbwt-name X     use GBWT index X (required with -w)" << endl
         << "" << endl;
}

int main_kmers(int argc, char** argv) {

    if (argc == 2) {
        help_kmers(argv);
        return 1;
    }

    // General options.
    size_t kmer_size = 0;
    bool show_progress = false;

    // GCSA options. Head and tail for distributed kmer generation.
    bool gcsa_out = false;
    bool gcsa_binary = false;
    int64_t head_id = 0;
    int64_t tail_id = 0;

    // Minimizer options.
    bool minimizer_out = false;
    size_t window_size = 0;
    std::string gbwt_name;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

        {
            // General options.
            {"kmer-size", required_argument, 0, 'k'},
            {"threads", required_argument, 0, 't'},
            {"progress",  no_argument, 0, 'p'},

            // GCSA options.
            {"gcsa-out", no_argument, 0, 'g'},
            {"gcsa-binary", no_argument, 0, 'B'},
            {"head-id", required_argument, 0, 'H'},
            {"tail-id", required_argument, 0, 'T'},

            // Minimizer options.
            {"window-size", required_argument, 0, 'w'},
            {"gbwt-name", required_argument, 0, 'G'},

            // Obsolete options.
            {"edge-max", required_argument, 0, 'e'},
            {"forward-only", no_argument, 0, 'F'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "k:t:pgBH:T:w:G:e:Fh",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            // General options.
            case 'k':
                kmer_size = parse<size_t>(optarg);
                break;
            case 't':
                omp_set_num_threads(parse<int>(optarg));
                break;
            case 'p':
                show_progress = true;
                break;

            // GCSA options.
            case 'g':
                gcsa_out = true;
                break;
            case 'B':
                gcsa_out = true;
                gcsa_binary = true;
                break;
            case 'H':
                head_id = parse<int>(optarg);
                break;
            case 'T':
                tail_id = parse<int>(optarg);
                break;

            // Minimizer options.
            case 'w':
                minimizer_out = true;
                window_size = parse<size_t>(optarg);
                break;
            case 'G':
                gbwt_name = optarg;
                break;

            // Obsolete options.
            case 'e':
                cerr << "error: [vg kmers] Option --edge-max is obsolete. Use vg prune to prune the graph instead." << endl;
                std::exit(EXIT_FAILURE);
                break;
            case 'F':
                cerr << "error: [vg kmers] Option --forward-only is obsolete" << endl;
                std::exit(EXIT_FAILURE);
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

    if (kmer_size == 0) {
        cerr << "error: [vg kmers] --kmer-size was not specified" << endl;
        std::exit(EXIT_FAILURE);
    }
    if (gcsa_out && minimizer_out) {
        cerr << "error: [vg kmers] Cannot output minimizers in GCSA format" << endl;
        std::exit(EXIT_FAILURE);
    }
    if (minimizer_out && gbwt_name.empty()) {
        cerr << "error: [vg kmers] Minimizers require a GBWT index" << endl;
        std::exit(EXIT_FAILURE);
    }

    vector<string> graph_file_names;
    while (optind < argc) {
        string file_name = get_input_file_name(optind, argc, argv);
        graph_file_names.push_back(file_name);
    }

    VGset graphs(graph_file_names);

    graphs.show_progress = show_progress;

    if (gcsa_out) {
        if (!gcsa_binary) {
            graphs.write_gcsa_kmers_ascii(cout, kmer_size, head_id, tail_id);
        } else {
            size_t limit = ~(size_t)0;
            graphs.write_gcsa_kmers_binary(cout, kmer_size, limit, head_id, tail_id);
        }
    } else if (minimizer_out) {
        // FIXME This is a placeholder.
        std::set<std::pair<std::string, pos_t>> minimizers;
        gbwt::GBWT haplotypes;
        sdsl::load_from_file(haplotypes, gbwt_name);
        auto lambda = [kmer_size, window_size, &minimizers](const GBWTTraversal& window) {
            std::string kmer = window.seq.substr(0, kmer_size);
            auto iter = window.traversal.begin();
            pos_t pos = iter->first, min_pos = iter->first;
            for (size_t i = 1; i < window_size; i++) {
                std::string temp = window.seq.substr(i, kmer_size);
                get_offset(pos)++;
                if (offset(pos) >= offset(iter->second)) {
                    ++iter;
                    pos = iter->first;
                }
                if (temp < kmer) {
                    kmer = temp;
                    min_pos = pos;
                }
            }
#pragma omp critical (minimizers)
            minimizers.insert(std::make_pair(kmer, min_pos));
        };
        graphs.for_each_kmer_parallel(haplotypes, kmer_size + window_size - 1, lambda);
        for (auto& kmer : minimizers) {
            cout << kmer.first << "\t"
                 << id(kmer.second) << (is_rev(kmer.second) ? ":-" : ":") << offset(kmer.second) << endl;
        }
    } else {
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


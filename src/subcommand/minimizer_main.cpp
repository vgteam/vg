/** \file minimizer_main.cpp
 *
 * Defines the "vg minimizer" subcommand, which builds the experimental
 * minimizer index.
 *
 * The index contains the lexicographically smallest kmer in a window of w
 * successive kmers. If the kmer contains characters other than A, C, G, and
 * T, it will not be indexed.
 *
 * By default, the index contains all minimizers in the graph. Option
 * --max-occs can be used to specify the maximum number of occurrences for
 * a kmer. Kmers more frequent than that will be removed from the index.
 *
 * At the moment, it is only possible to index haplotype-consistent kmers
 * by specifying --gbwt-name. The support for indexing all kmers in a graph
 * will be added later.
 */

#include "../gbwt_helper.hpp"
#include "../minimizer.hpp"
#include "../vg_set.hpp"
#include "subcommand.hpp"

#include <fstream>
#include <iostream>
#include <vector>

#include <getopt.h>
#include <omp.h>

using namespace vg;
using namespace vg::subcommand;


void help_minimizer(char** argv) {
    std::cerr << "usage: " << argv[0] << " minimizer [options] <graph1.vg> [graph2.vg ...]" << std::endl;
    std::cerr << "Builds a minimizer index of the graph(s)." << std::endl;
    std::cerr << "    -k, --kmer-length N    length of the kmers in the index (default: " << MinimizerIndex::KMER_LENGTH << ")" << std::endl;
    std::cerr << "    -w, --window-length N  index the smallest kmer in a window of N kmers (default: " << MinimizerIndex::WINDOW_LENGTH << ")" << std::endl;
    std::cerr << "    -m, --max-occs N       do not index minimizers with more than N occurrences" << std::endl;
    std::cerr << "    -i, --index-name X     store the index to file X (required)" << std::endl;
    std::cerr << "    -l, --load-index X     load the index from file X and insert the new kmers into it" << std::endl;
    std::cerr << "                           (overrides --kmer-length, --window-length, and --max-occs)" << std::endl;
    std::cerr << "    -g, --gbwt-name X      index only haplotype-consistent kmers using the GBWT index in file X" << std::endl;
    std::cerr << "    -p, --progress         show progress information" << std::endl;
    std::cerr << "    -t, --threads N        use N threads for index construction (default: " << omp_get_max_threads() << ")" << std::endl;
}

int main_minimizer(int argc, char** argv) {

    if (argc == 2) {
        help_minimizer(argv);
        return 1;
    }

    // Command-line options.
    size_t kmer_length = MinimizerIndex::KMER_LENGTH;
    size_t window_length = MinimizerIndex::WINDOW_LENGTH;
    size_t max_occs = MinimizerIndex::MAX_OCCS;
    std::string index_name, load_index, gbwt_name;
    bool progress = false;
    int threads;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "kmer-length", required_argument, 0, 'k' },
            { "window-length", required_argument, 0, 'w' },
            { "max-occs", required_argument, 0, 'm' },
            { "index-name", required_argument, 0, 'i' },
            { "load-index", required_argument, 0, 'l' },
            { "gbwt-name", required_argument, 0, 'g' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "k:w:m:i:l:g:pt:h", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'k':
            kmer_length = parse<size_t>(optarg);
            break;
        case 'w':
            window_length = parse<size_t>(optarg);
            break;
        case 'm':
            max_occs = parse<size_t>(optarg);
            break;
        case 'i':
            index_name = optarg;
            break;
        case 'l':
            load_index = optarg;
            break;
        case 'g':
            gbwt_name = optarg;
            break;
        case 'p':
            progress = true;
            break;
        case 't':
            threads = parse<int>(optarg);
            threads = std::min(threads, omp_get_max_threads());
            threads = std::max(threads, 1);
            omp_set_num_threads(threads);
            break;

        case 'h':
        case '?':
            help_minimizer(argv);
            return 1;
        default:
            std::abort();
        }
    }
    if (index_name.empty()) {
        std::cerr << "[vg minimizer]: option --index-name is required" << std::endl;
        return 1;
    }
    if (gbwt_name.empty()) {
        std::cerr << "[vg minimizer]: option --gbwt-name is currently required" << std::endl;
        return 1;
    }

    double start = gbwt::readTimer();

    // Input graphs.
    std::vector<std::string> file_names;
    while (optind < argc) {
        file_names.push_back(get_input_file_name(optind, argc, argv));
    }
    if (file_names.empty()) {
        help_minimizer(argv);
        return 1;
    }
    if (progress) {
        std::cerr << "Input files: " << file_names.size() << std::endl;
    }
    VGset graphs(file_names);

    // GBWT index.
    gbwt::GBWT gbwt_index;
    if (!gbwt_name.empty()) {
        if (progress) {
            std::cerr << "Loading GBWT index " << gbwt_name << std::endl;
        }
        sdsl::load_from_file(gbwt_index, gbwt_name);
    }

    // Minimizer index.
    MinimizerIndex index(kmer_length, window_length, max_occs);
    if (!load_index.empty()) {
        if (progress) {
            std::cerr << "Loading the index from " << load_index << std::endl;
        }
        std::ifstream in(load_index, std::ios_base::binary);
        if (!in) {
            std::cerr << "error: [vg minimizer] Cannot open index file " << load_index << " for reading" << std::endl;
            return 1;
        }
        index.load(in);
        in.close();
    }

    // Build the index.
    if (progress) {
        std::cerr << "Building the index" << std::endl;
    }
    auto lambda = [&index](const GBWTTraversal& window) {
        MinimizerIndex::minimizer_type minimizer = index.minimizer(window.seq.begin(), window.seq.end());
        if (minimizer.first == MinimizerIndex::NO_KEY) {
            return;
        }
        auto iter = window.traversal.begin();
        while (offset(iter->second) - offset(iter->first) <= minimizer.second) {
            minimizer.second -= offset(iter->second) - offset(iter->first);
            ++iter;
        }
        pos_t pos = iter->first;
        get_offset(pos) += minimizer.second;
#pragma omp critical (minimizer_index)
        index.insert(minimizer.first, pos);
    };
    graphs.for_each_kmer_parallel(gbwt_index, index.k() + index.w() - 1, lambda);

    // Index statistics.
    if (progress) {
        std::cerr << index.size() << " keys (" << index.unique_keys() << " unique, " << index.frequent_keys() << " too frequent)" << std::endl;
        std::cerr << "Minimizer occurrences: " << index.values() << std::endl;
        std::cerr << "Load factor: " << index.load_factor() << std::endl;
    }

    // Serialize the index.
    if (progress) {
        std::cerr << "Writing the index to " << index_name << std::endl;
    }
    std::ofstream out(index_name, std::ios_base::binary);
    if (!out) {
        std::cerr << "error: [vg minimizer] Cannot open index file " << index_name << " for writing" << std::endl;
        return 1;
    }
    index.serialize(out);
    out.close();

    double seconds = gbwt::readTimer() - start;
    if (progress) {
        std::cout << "Time usage: " << seconds << " seconds" << std::endl;
        std::cout << "Memory usage: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
    }

    return 0;
}

// FIXME change from DEVELOPMENT to TOOLKIT or PIPELINE later
// Register subcommand
static Subcommand vg_minimizer("minimizer", "build a minimizer index", DEVELOPMENT, main_minimizer);


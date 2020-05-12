/** \file minimizer_main.cpp
 *
 * Defines the "vg minimizer" subcommand, which builds the experimental
 * minimizer index.
 *
 * The index contains the lexicographically smallest kmer in a window of w
 * successive kmers and their reverse complements. If the kmer contains
 * characters other than A, C, G, and T, it will not be indexed.
 *
 * The index contains either all or haplotype-consistent minimizers. Indexing all
 * minimizers from complex graph regions can take a long time (e.g. 65 hours
 * vs 10 minutes for 1000GP), because many windows have the same minimizer.
 * As the total number of minimizers is manageable (e.g. 2.1 billion vs.
 * 1.4 billion for 1000GP), it should be possible to develop a better
 * algorithm for finding the minimizers.
 *
 * A quick idea:
 * - For each node v, extract the subgraph for the windows starting in v.
 * - Extract all k'-mers from the subgraph and use them to determine where the
 *   minimizers can start.
 */

#include "subcommand.hpp"

#include <vg/io/vpkg.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

#include <getopt.h>
#include <omp.h>

#include <gbwtgraph/index.h>

#include "../min_distance.hpp"
#include "../handle.hpp"
#include "../utility.hpp"

using namespace vg;
using namespace vg::subcommand;

// Using too many threads just wastes CPU time without speeding up the construction.
constexpr int DEFAULT_MAX_THREADS = 16;

int get_default_threads() {
    return std::min(omp_get_max_threads(), DEFAULT_MAX_THREADS);
}

void help_minimizer(char** argv) {
    std::cerr << "usage: " << argv[0] << " minimizer -g gbwt_name -i index_name [options] graph" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Builds a minimizer index of the haplotypes in the GBWT index over the input graph." << std::endl;
    std::cerr << "The graph can be any HandleGraph, but it will be transformed into a GBWTGraph. The" << std::endl;
    std::cerr << "transformation can be avoided by providing a prebuilt GBWTGraph and using option -G." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Required options:" << std::endl;
    std::cerr << "    -g, --gbwt-name X       use the GBWT index in file X" << std::endl;
    std::cerr << "    -i, --index-name X      store the index to file X" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Minimizer options:" << std::endl;
    std::cerr << "    -k, --kmer-length N     length of the kmers in the index (default " << gbwtgraph::DefaultMinimizerIndex::key_type::KMER_LENGTH << ", max " << gbwtgraph::DefaultMinimizerIndex::key_type::KMER_MAX_LENGTH << ")" << std::endl;
    std::cerr << "    -w, --window-length N   choose the minimizer from a window of N kmers (default " << gbwtgraph::DefaultMinimizerIndex::key_type::WINDOW_LENGTH << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "    -d, --distance-index X  annotate the hits with positions in this distance index" << std::endl;
    std::cerr << "    -l, --load-index X      load the index from file X and insert the new kmers into it" << std::endl;
    std::cerr << "                            (overrides --kmer-length and --window-length)" << std::endl;
    std::cerr << "    -G, --gbwt-graph        the input graph is a GBWTGraph" << std::endl;
    std::cerr << "    -p, --progress          show progress information" << std::endl;
    std::cerr << "    -t, --threads N         use N threads for index construction (default " << get_default_threads() << ")" << std::endl;
    std::cerr << "                            (using more than " << DEFAULT_MAX_THREADS << " threads rarely helps)" << std::endl;
    std::cerr << std::endl;
}

int main_minimizer(int argc, char** argv) {

    if (argc <= 6) {
        help_minimizer(argv);
        return 1;
    }

    // Command-line options.
    size_t kmer_length = gbwtgraph::DefaultMinimizerIndex::key_type::KMER_LENGTH;
    size_t window_length = gbwtgraph::DefaultMinimizerIndex::key_type::WINDOW_LENGTH;
    std::string index_name, distance_name, load_index, gbwt_name, graph_name;
    bool is_gbwt_graph = false;
    bool progress = false;
    int threads = get_default_threads();

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "gbwt-name", required_argument, 0, 'g' },
            { "index-name", required_argument, 0, 'i' },
            { "kmer-length", required_argument, 0, 'k' },
            { "window-length", required_argument, 0, 'w' },
            { "distance-index", required_argument, 0, 'd' },
            { "load-index", required_argument, 0, 'l' },
            { "gbwt-graph", no_argument, 0, 'G' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:i:k:w:d:l:Gpt:h", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'g':
            gbwt_name = optarg;
            break;
        case 'i':
            index_name = optarg;
            break;
        case 'k':
            kmer_length = parse<size_t>(optarg);
            break;
        case 'w':
            window_length = parse<size_t>(optarg);
            break;
        case 'd':
            distance_name = optarg;
            break;
        case 'l':
            load_index = optarg;
            break;
        case 'G':
            is_gbwt_graph = true;
            break;
        case 'p':
            progress = true;
            break;
        case 't':
            threads = parse<int>(optarg);
            threads = std::min(threads, omp_get_max_threads());
            threads = std::max(threads, 1);
            break;

        case 'h':
        case '?':
            help_minimizer(argv);
            return 1;
        default:
            std::abort();
        }
    }
    if (index_name.empty() || gbwt_name.empty()) {
        std::cerr << "[vg minimizer]: options --index-name and --gbwt-name are required" << std::endl;
        return 1;
    }
    if (optind + 1 != argc) {
        help_minimizer(argv);
        return 1;
    }
    graph_name = argv[optind];
    omp_set_num_threads(threads);

    double start = gbwt::readTimer();

    // GBWT index.
    if (progress) {
        std::cerr << "Loading GBWT index " << gbwt_name << std::endl;
    }
    std::unique_ptr<gbwt::GBWT> gbwt_index(vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name));

    // GBWTGraph.
    std::unique_ptr<gbwtgraph::GBWTGraph> gbwt_graph;
    if (is_gbwt_graph) {
        if (progress) {
            std::cerr << "Loading GBWTGraph " << graph_name << std::endl;
        }
        gbwt_graph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(graph_name);
        gbwt_graph->set_gbwt(*gbwt_index);
    } else {
        if (progress) {
            std::cerr << "Loading input graph " << graph_name << std::endl;
        }
        std::unique_ptr<HandleGraph> input_graph(vg::io::VPKG::load_one<HandleGraph>(graph_name));
        if (progress) {
            std::cerr << "Building GBWTGraph" << std::endl;
        }
        gbwt_graph.reset(new gbwtgraph::GBWTGraph(*gbwt_index, *input_graph));
    }

    // Minimizer index.
    std::unique_ptr<gbwtgraph::DefaultMinimizerIndex> index(new gbwtgraph::DefaultMinimizerIndex(kmer_length, window_length));
    if (!load_index.empty()) {
        if (progress) {
            std::cerr << "Loading MinimizerIndex " << load_index << std::endl;
        }
        index = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(load_index);
    }

    // Distance index.
    std::unique_ptr<MinimumDistanceIndex> distance_index;
    if (!distance_name.empty()) {
        if (progress) {
            std::cerr << "Loading MinimumDistanceIndex " << distance_name << std::endl;
        }
        distance_index = vg::io::VPKG::load_one<MinimumDistanceIndex>(distance_name);
    }

    // Build the index.
    if (progress) {
        std::cerr << "Building MinimizerIndex" << std::endl;
    }
    if (distance_name.empty()) {
        gbwtgraph::index_haplotypes(*gbwt_graph, *index, [](const pos_t&) -> gbwtgraph::payload_type {
            return MIPayload::NO_CODE;
        });
    } else {
        gbwtgraph::index_haplotypes(*gbwt_graph, *index, [&](const pos_t& pos) -> gbwtgraph::payload_type {
            return MIPayload::encode(distance_index->get_minimizer_distances(pos));
        });
    }
    gbwt_graph.reset(nullptr);
    gbwt_index.reset(nullptr);

    // Index statistics.
    if (progress) {
        std::cerr << index->size() << " keys (" << index->unique_keys() << " unique)" << std::endl;
        std::cerr << "Minimizer occurrences: " << index->values() << std::endl;
        std::cerr << "Load factor: " << index->load_factor() << std::endl;
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Construction so far: " << seconds << " seconds" << std::endl;
    }

    // Serialize the index.
    if (progress) {
        std::cerr << "Writing the index to " << index_name << std::endl;
    }
    vg::io::VPKG::save(*index, index_name);

    if (progress) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Time usage: " << seconds << " seconds" << std::endl;
        std::cerr << "Memory usage: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
    }

    return 0;
}

// FIXME change from DEVELOPMENT to TOOLKIT or PIPELINE later
// Register subcommand
static Subcommand vg_minimizer("minimizer", "build a minimizer index", DEVELOPMENT, main_minimizer);


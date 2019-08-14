/** \file minimizer_main.cpp
 *
 * Defines the "vg minimizer" subcommand, which builds the experimental
 * minimizer index.
 *
 * The index contains the lexicographically smallest kmer in a window of w
 * successive kmers and their reverse complements. If the kmer contains
 * characters other than A, C, G, and T, it will not be indexed.
 *
 * By default, the index contains all minimizers in the graph. Option
 * --max-occs can be used to specify the maximum number of occurrences for
 * a kmer. Kmers more frequent than that will be removed from the index.
 *
 * The index contains either all minimizers or haplotype-consistent minimizers
 * (option --gbwt-name). Indexing all minimizers from complex graph regions
 * can take a long time (e.g. 65 hours vs 30 minutes for 1000GP), because many
 * windows have the same minimizer. As the total number of minimizers is
 * manageable (e.g. 2.1 billion vs. 1.4 billion for 1000GP), it should be
 * possible to develop a better algorithm for finding the minimizers.
 *
 * A quick idea:
 * - For each node v, extract the subgraph for the windows starting in v.
 * - Extract all k'-mers from the subgraph and use them to determine where the
 *   minimizers can start.
 */

#include "../gapless_extender.hpp"
#include "../gbwt_helper.hpp"
#include "../minimizer.hpp"
#include "subcommand.hpp"

#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>

#include <iostream>
#include <vector>

#include <getopt.h>
#include <omp.h>

using namespace vg;
using namespace vg::subcommand;


void help_minimizer(char** argv) {
    std::cerr << "usage: " << argv[0] << " minimizer [options] graph.xg" << std::endl;
    std::cerr << "Builds a minimizer index of the graph in the XG index." << std::endl;
    std::cerr << "    -k, --kmer-length N    length of the kmers in the index (default: " << MinimizerIndex::KMER_LENGTH << ")" << std::endl;
    std::cerr << "    -w, --window-length N  index the smallest kmer in a window of N kmers (default: " << MinimizerIndex::WINDOW_LENGTH << ")" << std::endl;
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
    std::string index_name, load_index, gbwt_name, xg_name;
    bool progress = false;
    int threads = omp_get_max_threads();

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "kmer-length", required_argument, 0, 'k' },
            { "window-length", required_argument, 0, 'w' },
            { "index-name", required_argument, 0, 'i' },
            { "load-index", required_argument, 0, 'l' },
            { "gbwt-name", required_argument, 0, 'g' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "k:w:i:l:g:pt:h", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'k':
            kmer_length = parse<size_t>(optarg);
            break;
        case 'w':
            window_length = parse<size_t>(optarg);
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
    if (optind + 1 != argc) {
        help_minimizer(argv);
        return 1;
    }
    xg_name = argv[optind];

    double start = gbwt::readTimer();

    // Input graph.
    if (progress) {
        std::cerr << "Loading XG index " << xg_name << std::endl;
    }
    std::unique_ptr<PathPositionHandleGraph> xg_index;
    xg_index = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_name);

    // Minimizer index.
    std::unique_ptr<MinimizerIndex> index(new MinimizerIndex(kmer_length, window_length));
    if (!load_index.empty()) {
        if (progress) {
            std::cerr << "Loading minimizer index " << load_index << std::endl;
        }
        index = vg::io::VPKG::load_one<MinimizerIndex>(load_index);
    }

    // GBWT-backed graph.
    std::unique_ptr<gbwt::GBWT> gbwt_index;
    std::unique_ptr<GBWTGraph> gbwt_graph;
    if (!gbwt_name.empty()) {
        if (progress) {
            std::cerr << "Loading GBWT index " << gbwt_name << std::endl;
        }
        gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name);
        if (progress) {
            std::cerr << "Building GBWT-backed graph" << std::endl;
        }
        gbwt_graph.reset(new GBWTGraph(*gbwt_index, *xg_index));
        xg_index.reset(nullptr); // The XG index is no longer needed.
    }

    // Minimizer caching.
    std::vector<std::vector<std::pair<MinimizerIndex::minimizer_type, pos_t>>> cache(threads);
    constexpr size_t MINIMIZER_CACHE_SIZE = 1024;
    auto flush_cache = [&](int thread_id) {
        gbwt::removeDuplicates(cache[thread_id], false);
        #pragma omp critical (minimizer_index)
        {
            for (auto& hit : cache[thread_id]) {
                index->insert(hit.first, hit.second);
            }
        }
        cache[thread_id].clear();
    };

    // Build the index.
    if (progress) {
        std::cerr << "Building the index" << std::endl;
    }
    const HandleGraph& graph = *(gbwt_name.empty() ?
                                 static_cast<const HandleGraph*>(xg_index.get()) :
                                 static_cast<const HandleGraph*>(gbwt_graph.get()));
    auto lambda = [&](const std::vector<handle_t>& traversal, const std::string& seq) {
        std::vector<MinimizerIndex::minimizer_type> minimizers = index->minimizers(seq);
        auto iter = traversal.begin();
        size_t node_start = 0;
        int thread_id = omp_get_thread_num();
        for (MinimizerIndex::minimizer_type& minimizer : minimizers) {
            if (minimizer.empty()) {
                continue;
            }
            // Find the node covering minimizer starting position.
            size_t node_length = graph.get_length(*iter);
            while (node_start + node_length <= minimizer.offset) {
                node_start += node_length;
                ++iter;
                node_length = graph.get_length(*iter);
            }
            pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start };
            if (minimizer.is_reverse) {
                pos = reverse_base_pos(pos, node_length);
            }
            if (!MinimizerIndex::valid_offset(pos)) {
                #pragma omp critical (cerr)
                {
                    std::cerr << "error: [vg minimizer] node offset " << offset(pos) << " is too large" << std::endl;
                }
                std::exit(EXIT_FAILURE);
            }
            cache[thread_id].emplace_back(minimizer, pos);
        }
        if (cache[thread_id].size() >= MINIMIZER_CACHE_SIZE) {
            flush_cache(thread_id);
        }
    };
    // We do a lot of redundant work by traversing both orientations and finding almost the same minimizers
    // in each orientation. If we consider only the windows starting in forward (reverse) orientation,
    // we may skip windows that cross from a reverse node to a forward node (from a forward node to a
    // reverse node).
    if (gbwt_name.empty()) {
        for_each_window(*xg_index, index->k() + index->w() - 1, lambda, (threads > 1));
    } else {
        for_each_haplotype_window(*gbwt_graph, index->k() + index->w() - 1, lambda, (threads > 1));
    }
    for (int thread_id = 0; thread_id < threads; thread_id++) {
        flush_cache(thread_id);
    }
    xg_index.reset(nullptr);
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


/** \file minimizer_main.cpp
 *
 * Defines the "vg minimizer" subcommand, which builds the minimizer index.
 *
 * The index contains the lexicographically smallest kmer in a window of w
 * successive kmers and their reverse complements. If the kmer contains
 * characters other than A, C, G, and T, it will not be indexed.
 *
 * The index contains either all or haplotype-consistent minimizers. Indexing all
 * minimizers from complex graph regions can take a long time (e.g. tens of hours
 * vs 5-10 minutes for 1000GP), because many windows have the same minimizer.
 * As the total number of minimizers is manageable (e.g. 1.5x more for 1000GP)
 * it should be possible to develop a better algorithm for finding the minimizers.
 *
 * A quick idea for indexing the entire graph:
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

#include "../index_manager.hpp"

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/index.h>

#include "../min_distance.hpp"
#include "../handle.hpp"
#include "../utility.hpp"

using namespace vg;

// Using too many threads just wastes CPU time without speeding up the construction.
constexpr int DEFAULT_MAX_THREADS = 16;

int get_default_threads() {
    return std::min(omp_get_max_threads(), DEFAULT_MAX_THREADS);
}

size_t get_default_k() {
    return IndexManager::minimizer_k;
}

size_t get_default_w() {
    return IndexManager::minimizer_w;
}

size_t get_default_s() {
    return IndexManager::minimizer_s;
}

enum input_type { input_graph, input_gg, input_gbz };

void help_minimizer(char** argv) {
    std::cerr << "usage: " << argv[0] << " minimizer [options] graph" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Builds a (w, k)-minimizer index or a (k, s)-syncmer index of the threads in the GBWT" << std::endl;
    std::cerr << "index. The graph can be any HandleGraph, which will be transformed into a GBWTGraph." << std::endl;
    std::cerr << "The transformation can be avoided by providing a GBWTGraph and using option -G or -Z." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Required options:" << std::endl;
    std::cerr << "    -g, --gbwt-name X       use the GBWT index in file X (or specify -Z)" << std::endl;
    std::cerr << "    -o, --output-name X     store the index to file X" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Minimizer options:" << std::endl;
    std::cerr << "    -k, --kmer-length N     length of the kmers in the index (default " << get_default_k() << ", max " << gbwtgraph::DefaultMinimizerIndex::key_type::KMER_MAX_LENGTH << ")" << std::endl;
    std::cerr << "    -w, --window-length N   choose the minimizer from a window of N kmers (default " << get_default_w() << ")" << std::endl;
    std::cerr << "    -b, --bounded-syncmers  index bounded syncmers instead of minimizers" << std::endl;
    std::cerr << "    -s, --smer-length N     use smers of length N in bounded syncmers (default " << get_default_s() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "    -d, --distance-index X  annotate the hits with positions in this distance index" << std::endl;
    std::cerr << "    -l, --load-index X      load the index from file X and insert the new kmers into it" << std::endl;
    std::cerr << "                            (overrides -k, -w, -b, and -s)" << std::endl;
    std::cerr << "    -G, --gbwt-graph        the input graph is a GBWTGraph" << std::endl;
    std::cerr << "    -Z, --gbz-format        the input graph is in GBZ format (-g is unnecessary)" << std::endl;
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
    size_t kmer_length = get_default_k();
    size_t window_length = get_default_w();
    size_t smer_length = get_default_s();
    std::string output_name, distance_name, load_index, gbwt_name, graph_name;
    bool use_syncmers = false;
    input_type input = input_graph;
    bool progress = false;
    int threads = get_default_threads();

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "gbwt-name", required_argument, 0, 'g' },
            { "output-name", required_argument, 0, 'o' },
            { "index-name", required_argument, 0, 'i' }, // deprecated
            { "kmer-length", required_argument, 0, 'k' },
            { "window-length", required_argument, 0, 'w' },
            { "bounded-syncmers", no_argument, 0, 'b' },
            { "smer-length", required_argument, 0, 's' },
            { "distance-index", required_argument, 0, 'd' },
            { "load-index", required_argument, 0, 'l' },
            { "gbwt-graph", no_argument, 0, 'G' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:o:i:k:w:bs:d:l:GZpt:h", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'g':
            gbwt_name = optarg;
            break;
        case 'o':
            output_name = optarg;
            break;
        case 'i':
            std::cerr << "warning: [vg minimizer] --index-name is deprecated, use --output-name instead" << std::endl;
            output_name = optarg;
            break;
        case 'k':
            kmer_length = parse<size_t>(optarg);
            break;
        case 'w':
            window_length = parse<size_t>(optarg);
            break;
        case 'b':
            use_syncmers = true;
            break;
        case 's':
            smer_length = parse<size_t>(optarg);
            break;
        case 'd':
            distance_name = optarg;
            break;
        case 'l':
            load_index = optarg;
            break;
        case 'G':
            input = input_gg;
            break;
        case 'Z':
            input = input_gbz;
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
    if (output_name.empty() || (gbwt_name.empty() && input != input_gbz)) {
        std::cerr << "[vg minimizer]: options --output-name and --gbwt-name are required" << std::endl;
        return 1;
    }
    if (optind + 1 != argc) {
        help_minimizer(argv);
        return 1;
    }
    graph_name = argv[optind];
    omp_set_num_threads(threads);

    double start = gbwt::readTimer();

    // We use GBWT and GBWTGraph in this GBZ wrapper.
    gbwtgraph::GBZ gbz;
    if (input == input_graph) {
        if (progress) {
            std::cerr << "Loading GBWT index " << gbwt_name << std::endl;
        }
        std::unique_ptr<gbwt::GBWT> gbwt_index(vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name));
        if (progress) {
            std::cerr << "Loading input graph " << graph_name << std::endl;
        }
        std::unique_ptr<HandleGraph> graph(vg::io::VPKG::load_one<HandleGraph>(graph_name));
        if (progress) {
            std::cerr << "Building GBWTGraph" << std::endl;
        }
        gbz = gbwtgraph::GBZ(gbwt_index, *graph);
    } else if (input == input_gg) {
        // This is a bit ugly, because the structures are inside vg wrappers.
        if (progress) {
            std::cerr << "Loading GBWT index " << gbwt_name << std::endl;
        }
        std::unique_ptr<gbwt::GBWT> gbwt_index(vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name));
        gbz.index = std::move(*gbwt_index);
        if (progress) {
            std::cerr << "Loading GBWTGraph " << graph_name << std::endl;
        }
        std::unique_ptr<gbwtgraph::GBWTGraph> graph(vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(graph_name));
        gbz.graph = std::move(*graph);
        gbz.graph.set_gbwt(gbz.index);
    } else if (input == input_gbz) {
        if (progress) {
            std::cerr << "Loading GBZ file " << graph_name << std::endl;
        }
        sdsl::simple_sds::load_from(gbz, graph_name);
    }

    // Minimizer index.
    std::unique_ptr<gbwtgraph::DefaultMinimizerIndex> index;
    if (load_index.empty()) {
        index = std::make_unique<gbwtgraph::DefaultMinimizerIndex>(kmer_length, (use_syncmers ? smer_length : window_length), use_syncmers);
    } else {
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
        std::cerr << "Building MinimizerIndex with k = " << index->k();
        if (index->uses_syncmers()) {
            std::cerr << ", s = " << index->s();
        } else {
            std::cerr << ", w = " << index->w();
        }
        std::cerr << std::endl;
    }
    if (distance_name.empty()) {
        gbwtgraph::index_haplotypes(gbz.graph, *index, [](const pos_t&) -> gbwtgraph::payload_type {
            return MIPayload::NO_CODE;
        });
    } else {
        gbwtgraph::index_haplotypes(gbz.graph, *index, [&](const pos_t& pos) -> gbwtgraph::payload_type {
            return MIPayload::encode(distance_index->get_minimizer_distances(pos));
        });
    }

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
        std::cerr << "Writing the index to " << output_name << std::endl;
    }
    vg::io::VPKG::save(*index, output_name);

    if (progress) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Time usage: " << seconds << " seconds" << std::endl;
        std::cerr << "Memory usage: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
    }

    return 0;
}

// Register subcommand
static vg::subcommand::Subcommand vg_minimizer("minimizer", "build a minimizer index or a syncmer index", vg::subcommand::TOOLKIT, main_minimizer);

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

#include "../gbwtgraph_helper.hpp"
#include "../gbwt_helper.hpp"
#include "../index_registry.hpp"
#include "../utility.hpp"
#include "../handle.hpp"
#include "../snarl_distance_index.hpp"
#include "../zip_code.hpp"
#include "../minimizer_index_helper.hpp"
#include <gbwtgraph/index.h>

using namespace vg;

// Using too many threads just wastes CPU time without speeding up the construction.
constexpr int DEFAULT_MAX_THREADS = 16;

// For weighted minimizers.
constexpr size_t DEFAULT_THRESHOLD = 500; // This should be Giraffe hard hit cap.
constexpr size_t DEFAULT_ITERATIONS = 3;
constexpr size_t MAX_ITERATIONS = gbwtgraph::MinimizerHeader::FLAG_WEIGHT_MASK >> gbwtgraph::MinimizerHeader::FLAG_WEIGHT_OFFSET;
constexpr size_t HASH_TABLE_MIN_WIDTH = 10;
constexpr size_t HASH_TABLE_MAX_WIDTH = 36;

int get_default_threads() {
    return std::min(omp_get_max_threads(), DEFAULT_MAX_THREADS);
}

void help_minimizer(char** argv) {
    std::cerr << "usage: " << argv[0] << " minimizer [options] -d graph.dist -o graph.min graph" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Builds a (w, k)-minimizer index or a (k, s)-syncmer index of the threads in the" << std::endl;
    std::cerr << "GBWT. The graph can be any HandleGraph, which will be made into a GBWTGraph." << std::endl;
    std::cerr << "The transformation can be avoided by providing a GBWTGraph or a GBZ graph." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Required options:" << std::endl;
    std::cerr << "  -d, --distance-index FILE  annotate hits with positions in this distance index" << std::endl;
    std::cerr << "  -o, --output-name FILE     store the index in a file" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Minimizer options:" << std::endl;
    std::cerr << "  -k, --kmer-length N        length of the kmers in the index "
                                           << "[" << IndexingParameters::short_read_minimizer_k 
                                           << "] (max " << gbwtgraph::DefaultMinimizerIndex::key_type::KMER_MAX_LENGTH 
                                           << ")" << std::endl;
    std::cerr << "  -w, --window-length N      choose minimizer from a window of N kmers "
                                           << "[" << IndexingParameters::short_read_minimizer_w << "]" << std::endl;
    std::cerr << "  -c, --closed-syncmers      index closed syncmers instead of minimizers" << std::endl;
    std::cerr << "  -s, --smer-length N        use smers of length N in closed syncmers "
                                           << "[" << IndexingParameters::minimizer_s << "]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Weighted minimizers:" << std::endl;
    std::cerr << "  -W, --weighted             use weighted minimizers" << std::endl;
    std::cerr << "      --threshold N          downweight kmers with more than N hits "
                                           << "[" << DEFAULT_THRESHOLD << "]" << std::endl;
    std::cerr << "      --iterations N         downweight frequent kmers by N iterations "
                                           << "[" << DEFAULT_ITERATIONS << "]" << std::endl;
    std::cerr << "      --fast-counting        use the fast kmer counting algorithm (default)" << std::endl;
    std::cerr << "      --save-memory          use the space-efficient kmer counting algorithm" << std::endl;
    std::cerr << "      --hash-table N         use 2^N-cell hash tables for kmer counting" << std::endl;
    std::cerr << "                             (default: guess)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "  -z, --zipcode-name FILE    store the distances that are too big in afile" << std::endl;
    std::cerr << "                             if no -z, some distances may be discarded" << std::endl;
    std::cerr << "  -l, --load-index FILE      load this index and insert the new kmers into it" << std::endl;
    std::cerr << "                             (overrides minimizer / weighted minimizer options)" << std::endl;
    std::cerr << "  -g, --gbwt-name FILE       use this GBWT index (required with a non-GBZ graph)" << std::endl;
    std::cerr << "  -E, --rec-mode             use recombination-aware MinimizerIndex" << std::endl;
    std::cerr << "  -p, --progress             show progress information" << std::endl;
    std::cerr << "  -t, --threads N            use N threads for index construction "
                                           << "[" << get_default_threads() << "]" << std::endl;
    std::cerr << "                             (using more than " << DEFAULT_MAX_THREADS << " threads rarely helps)" << std::endl;
    std::cerr << "      --no-dist              build the index without distance index annotations" << std::endl;
    std::cerr << "                             (not recommended)" << std::endl;
    std::cerr << "  -h, --help                 print this help message to stderr and exit" << std::endl;
    std::cerr << std::endl;
}

template<typename PayloadType>
void construct_minimizer_dispatch(
    gbwtgraph::GBZ* gbz,
    const std::string& distance_name,
    const std::string& load_index,
    const std::string& zipcode_name,
    const std::string& output_name,
    size_t hash_table_size,
    size_t threshold,
    size_t iterations,
    double start,
    bool space_efficient_counting,
    bool weighted,
    bool use_syncmers,
    bool progress
) {
    // Minimizer index.
    using IndexType = gbwtgraph::MinimizerIndex<gbwtgraph::Key64, gbwtgraph::PositionPayload<PayloadType>>;
    IndexType index(IndexingParameters::short_read_minimizer_k, 
        (use_syncmers ? IndexingParameters::minimizer_s : IndexingParameters::short_read_minimizer_w),
        use_syncmers);
    
    // Helper for building the index.
    MinimizerIndexHelper<IndexType, PayloadType> mi_helper(gbz, index, progress);
    
    // find frequent kmers and set them in the index
    if (load_index.empty()) {
        if (weighted) {
            mi_helper.set_frequent_kmers(
                index.k(),
                threshold,
                space_efficient_counting,
                hash_table_size,
                iterations
            );
        }
    } else {
        load_minimizer(index, load_index, progress);
    }

    // Distance index 
    std::unique_ptr<SnarlDistanceIndex> distance_index;
    if (!distance_name.empty()) {
        // new distance index
        if (progress) {
            std::cerr << "Loading SnarlDistanceIndex from " << distance_name << std::endl;
        }
        distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(distance_name);
        distance_index->preload(true);
    }

    mi_helper.build(
        distance_index.get(), 
        zipcode_name,
        output_name
    );
}

int main_minimizer(int argc, char** argv) {

    if (argc <= 5) {
        help_minimizer(argv);
        return 1;
    }

    // Command-line options.
    std::string output_name, distance_name, zipcode_name, load_index, gbwt_name, graph_name;
    bool use_syncmers = false;
    bool weighted = false, space_efficient_counting = false;
    size_t threshold = DEFAULT_THRESHOLD, iterations = DEFAULT_ITERATIONS, hash_table_size = 0;
    bool progress = false;
    int threads = get_default_threads();
    bool require_distance_index = true;
    bool rec_mode = false;

    constexpr int OPT_THRESHOLD = 1001;
    constexpr int OPT_ITERATIONS = 1002;
    constexpr int OPT_FAST_COUNTING = 1003;
    constexpr int OPT_SAVE_MEMORY = 1004;
    constexpr int OPT_HASH_TABLE = 1005;
    constexpr int OPT_NO_DIST = 1100;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "gbwt-name", required_argument, 0, 'g' },
            { "rec-mode", no_argument, 0, 'E' },
            { "distance-index", required_argument, 0, 'd' },
            { "output-name", required_argument, 0, 'o' },
            { "index-name", required_argument, 0, 'i' }, // deprecated
            { "kmer-length", required_argument, 0, 'k' },
            { "window-length", required_argument, 0, 'w' },
            { "bounded-syncmers", no_argument, 0, 'b' }, // deprecated
            { "closed-syncmers", no_argument, 0, 'c' },
            { "smer-length", required_argument, 0, 's' },
            { "weighted", no_argument, 0, 'W' },
            { "threshold", required_argument, 0, OPT_THRESHOLD },
            { "iterations", required_argument, 0, OPT_ITERATIONS },
            { "fast-counting", no_argument, 0, OPT_FAST_COUNTING },
            { "save-memory", no_argument, 0, OPT_SAVE_MEMORY },
            { "hash-table", required_argument, 0, OPT_HASH_TABLE },
            { "zipcode-name", required_argument, 0, 'z' },
            { "load-index", required_argument, 0, 'l' },
            { "gbwt-graph", no_argument, 0, 'G' }, // deprecated
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { "no-dist", no_argument, 0, OPT_NO_DIST },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "g:Ed:o:i:k:w:bcs:Wz:l:Gpt:h?", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'g':
            gbwt_name = optarg;
            break;
        case 'E':
            rec_mode = true;
            std::cerr << "[vg minimizer] warning: --rec-mode is still under development" << std::endl;
            break;
        case 'd':
            distance_name = optarg;
            break;
        case 'o':
            output_name = optarg;
            break;
        case 'i':
            std::cerr << "[vg minimizer] warning: --index-name is deprecated, use --output-name instead" << std::endl;
            output_name = optarg;
            break;

        case 'k':
            IndexingParameters::short_read_minimizer_k = parse<size_t>(optarg);
            break;
        case 'w':
            IndexingParameters::short_read_minimizer_w = parse<size_t>(optarg);
            break;
        case 'b':
            std::cerr << "[vg minimizer] warning: --bounded-syncmers is deprecated, use --closed-syncmers instead" << std::endl;
            use_syncmers = true;
            break;
        case 'c':
            use_syncmers = true;
            break;
        case 's':
            IndexingParameters::minimizer_s = parse<size_t>(optarg);
            break;

        case 'W':
            weighted = true;
            break;
        case OPT_THRESHOLD:
            threshold = parse<size_t>(optarg);
            break;
        case OPT_ITERATIONS:
            iterations = parse<size_t>(optarg);
            iterations = std::max(iterations, size_t(1));
            iterations = std::min(iterations, MAX_ITERATIONS);
            break;
        case OPT_FAST_COUNTING:
            space_efficient_counting = false;
            break;
        case OPT_SAVE_MEMORY:
            space_efficient_counting = true;
            break;
        case OPT_HASH_TABLE:
            {
                size_t width = parse<size_t>(optarg);
                width = std::max(width, HASH_TABLE_MIN_WIDTH);
                width = std::min(width, HASH_TABLE_MAX_WIDTH);
                hash_table_size = size_t(1) << width;
            }
            break;

        case 'z':
            zipcode_name = optarg;
            break;
        case 'l':
            load_index = optarg;
            break;
        case 'G':
            std::cerr << "[vg minimizer] warning: --gbwt-graph is deprecated, graph format is now autodetected" << std::endl;
            break;
        case 'p':
            progress = true;
            break;
        case 't':
            threads = parse<int>(optarg);
            threads = std::min(threads, omp_get_max_threads());
            threads = std::max(threads, 1);
            break;
        case OPT_NO_DIST:
            require_distance_index = false;
            break;

        case 'h':
        case '?':
            help_minimizer(argv);
            return 1;
        default:
            std::abort();
        }
    }
    if (output_name.empty()) {
        std::cerr << "error: [vg minimizer] option --output-name is required" << std::endl;
        return 1;
    }
    if (optind + 1 != argc) {
        help_minimizer(argv);
        return 1;
    }
    graph_name = argv[optind];
    if (require_distance_index && distance_name.empty()) {
        std::cerr << "error: [vg minimizer] one of options --distance-index and --no-dist is required" << std::endl;
        return 1;
    }
    if (!load_index.empty() || use_syncmers) {
        weighted = false;
    }
    omp_set_num_threads(threads);


    double start = gbwt::readTimer();
   
    // We use GBWT and GBWTGraph in this GBZ wrapper.
    unique_ptr<gbwtgraph::GBZ> gbz;
   
    // Load whatever the graph argument is
    if (progress) {
        std::cerr << "Loading input graph from " << graph_name << std::endl;
    }
    auto input = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, gbwtgraph::GBWTGraph, HandleGraph>(graph_name);
    if (get<0>(input)) {
        // We loaded a GBZ directly
        gbz = std::move(get<0>(input));
    } else if (get<1>(input)) {
        // We loaded a GBWTGraph and need to pair it with a GBWT
        if (gbwt_name.empty()) {
            std::cerr << "error: [vg minimizer] option --gbwt-name is required when using a GBWTGraph" << std::endl;
            return 1;
        }
        gbz.reset(new gbwtgraph::GBZ());
        gbz->graph = std::move(*get<1>(input));

        // Go get the GBWT
        load_gbwt(gbz->index, gbwt_name, progress);
        // And attach them together
        gbz->graph.set_gbwt(gbz->index);
    } else if (get<2>(input)) {
        // We got a normal HandleGraph
        if (gbwt_name.empty()) {
            std::cerr << "error: [vg minimizer] option --gbwt-name is required when using a HandleGraph" << std::endl;
            return 1;
        }
        gbz.reset(new gbwtgraph::GBZ());
        
        load_gbwt(gbz->index, gbwt_name, progress);
        if (progress) {
            std::cerr << "Building GBWTGraph" << std::endl;
        }
        gbz->graph = gbwtgraph::GBWTGraph(gbz->index, *get<2>(input));
    } else {
        std::cerr << "error: [vg minimizer] input graph is not a GBZ, GBWTGraph, or HandleGraph." << std::endl;
        return 1;
    }

    if (rec_mode) {
        construct_minimizer_dispatch<gbwtgraph::PayloadXL>(
            gbz.get(), distance_name, load_index, zipcode_name, output_name,
            hash_table_size, threshold, iterations, start,
            space_efficient_counting, weighted, use_syncmers, progress
        );
    } else {
        construct_minimizer_dispatch<gbwtgraph::Payload>(
            gbz.get(), distance_name, load_index, zipcode_name, output_name,
            hash_table_size, threshold, iterations, start,
            space_efficient_counting, weighted, use_syncmers, progress
        );
    }

    return 0;
}

// Register subcommand
static vg::subcommand::Subcommand vg_minimizer("minimizer", "build a minimizer index or a syncmer index",
                                               vg::subcommand::TOOLKIT, main_minimizer);

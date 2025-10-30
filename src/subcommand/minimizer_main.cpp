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

#include <gbwtgraph/index.h>

using namespace vg;

//------------------------------------------------------------------------------

// TODO: We may need more threads with -E / --rec-mode, as it does more computation with the hits.
// Using too many threads just wastes CPU time without speeding up the construction.
constexpr int DEFAULT_MAX_THREADS = 16;

//------------------------------------------------------------------------------

int get_default_threads() {
    return std::min(omp_get_max_threads(), DEFAULT_MAX_THREADS);
}

struct MinimizerConfig {
    // File names.
    std::string output_name, distance_name, zipcode_name, graph_name;

    // Construction parameters.
    MinimizerIndexParameters params;
    bool require_distance_index = true;

    // Other options.
    // Note that params also has a 'progress' field.
    bool progress = false;
    int threads = get_default_threads();
    bool describe = false;

    MinimizerConfig(int argc, char** argv, int max_threads, const Logger& logger);
};

using code_type = gbwtgraph::KmerEncoding::code_type;
using payload_type = ZipCode::payload_type;

//------------------------------------------------------------------------------

int main_minimizer(int argc, char** argv) {
    double start = gbwt::readTimer();
    Logger logger("vg minimizer");
    MinimizerConfig config(argc, argv, get_default_threads(), logger);

    if (config.describe) {
        // Just describe the index.
        describe_minimizer_index(config.output_name, std::cout);
        return 0;
    }

    // Load the graph.
    gbwtgraph::GBZ gbz;
    load_gbz(gbz, config.graph_name, config.progress);

    // Load the distance index if requested.
    std::unique_ptr<SnarlDistanceIndex> distance_index;
    if (!config.distance_name.empty()) {
        // new distance index
        if (config.progress) {
            logger.info() << "Loading SnarlDistanceIndex from " << config.distance_name << std::endl;
        }
        distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(config.distance_name);
        distance_index->preload(true);
    }

    ZipCodeCollection oversized_zipcodes;
    gbwtgraph::DefaultMinimizerIndex index = build_minimizer_index(
        gbz,
        distance_index.get(),
        &oversized_zipcodes,
        config.params
    );

    // Serialize the index and the oversized zipcodes.
    save_minimizer(index, config.output_name);
    if (!config.zipcode_name.empty()) {
        std::ofstream zip_out(config.zipcode_name);
        oversized_zipcodes.serialize(zip_out);
        zip_out.close();

    }

    if (config.progress) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Time usage: " << seconds << " seconds" << std::endl;
        std::cerr << "Memory usage: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
    }

    return 0;
}

//------------------------------------------------------------------------------

void help_minimizer(char** argv) {
    std::cerr << "usage:" << std::endl;
    std::cerr << "    " << argv[0] << " minimizer [options] -d graph.dist -o graph.min graph.gbz" << std::endl;
    std::cerr << "    " << argv[0] << " minimizer --describe graph.min" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Builds a (w, k)-minimizer index or a (k, s)-syncmer index." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Required options:" << std::endl;
    std::cerr << "  -d, --distance-index FILE  annotate hits with positions in this distance index" << std::endl;
    std::cerr << "  -o, --output-name FILE     store the index in a file" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Minimizer options:" << std::endl;
    std::cerr << "  -k, --kmer-length N        length of the kmers in the index "
                                           << "[" << gbwtgraph::DefaultMinimizerIndex::key_type::KMER_LENGTH 
                                           << "] (max " << gbwtgraph::DefaultMinimizerIndex::key_type::KMER_MAX_LENGTH 
                                           << ")" << std::endl;
    std::cerr << "  -w, --window-length N      choose minimizer from a window of N kmers "
                                           << "[" << gbwtgraph::DefaultMinimizerIndex::key_type::WINDOW_LENGTH << "]" << std::endl;
    std::cerr << "  -c, --closed-syncmers      index closed syncmers instead of minimizers" << std::endl;
    std::cerr << "  -s, --smer-length N        use smers of length N in closed syncmers "
                                           << "[" << IndexingParameters::minimizer_s << "]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Weighted minimizers:" << std::endl;
    std::cerr << "  -W, --weighted             use weighted minimizers" << std::endl;
    std::cerr << "      --threshold N          downweight kmers with more than N hits "
                                           << "[" << MinimizerIndexParameters::DEFAULT_THRESHOLD << "]" << std::endl;
    std::cerr << "      --iterations N         downweight frequent kmers by N iterations "
                                           << "[" << MinimizerIndexParameters::DEFAULT_ITERATIONS << "]" << std::endl;
    std::cerr << "      --fast-counting        use the fast kmer counting algorithm (default)" << std::endl;
    std::cerr << "      --save-memory          use the space-efficient kmer counting algorithm" << std::endl;
    std::cerr << "      --hash-table N         use 2^N-cell hash tables for kmer counting" << std::endl;
    std::cerr << "                             (default: guess)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "  -z, --zipcode-name FILE    store the distances that are too big in a file" << std::endl;
    std::cerr << "                             if no -z, some distances may be discarded" << std::endl;
    std::cerr << "  -E, --rec-mode             build recombination-aware MinimizerIndex" << std::endl;
    std::cerr << "  -p, --progress             show progress information" << std::endl;
    std::cerr << "  -t, --threads N            use N threads for index construction "
                                           << "[" << get_default_threads() << "]" << std::endl;
    std::cerr << "                             (using more than " << DEFAULT_MAX_THREADS << " threads rarely helps)" << std::endl;
    std::cerr << "      --no-dist              build the index without distance index annotations" << std::endl;
    std::cerr << "                             (not recommended)" << std::endl;
    std::cerr << "      --describe FILE        describe the minimizer index in FILE and exit" << std::endl;
    std::cerr << "  -h, --help                 print this help message to stderr and exit" << std::endl;
    std::cerr << std::endl;
}

MinimizerConfig::MinimizerConfig(int argc, char** argv, int max_threads, const Logger& logger) {
    if (argc < 3) {
        help_minimizer(argv);
        std::exit(EXIT_FAILURE);
    }

    constexpr int OPT_THRESHOLD = 1001;
    constexpr int OPT_ITERATIONS = 1002;
    constexpr int OPT_FAST_COUNTING = 1003;
    constexpr int OPT_SAVE_MEMORY = 1004;
    constexpr int OPT_HASH_TABLE = 1005;
    constexpr int OPT_NO_DIST = 1100;
    constexpr int OPT_DESCRIBE = 1101;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "rec-mode", no_argument, 0, 'E' },
            { "distance-index", required_argument, 0, 'd' },
            { "output-name", required_argument, 0, 'o' },
            { "kmer-length", required_argument, 0, 'k' },
            { "window-length", required_argument, 0, 'w' },
            { "closed-syncmers", no_argument, 0, 'c' },
            { "smer-length", required_argument, 0, 's' },
            { "weighted", no_argument, 0, 'W' },
            { "threshold", required_argument, 0, OPT_THRESHOLD },
            { "iterations", required_argument, 0, OPT_ITERATIONS },
            { "fast-counting", no_argument, 0, OPT_FAST_COUNTING },
            { "save-memory", no_argument, 0, OPT_SAVE_MEMORY },
            { "hash-table", required_argument, 0, OPT_HASH_TABLE },
            { "zipcode-name", required_argument, 0, 'z' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { "no-dist", no_argument, 0, OPT_NO_DIST },
            { "describe", required_argument, 0, OPT_DESCRIBE },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "Ed:o:k:w:cs:Wz:pt:h?", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'E':
            this->params.paths_in_payload = true;
            logger.warn() << "--rec-mode is still under development" << std::endl;
            break;
        case 'd':
            this->distance_name = optarg;
            break;
        case 'o':
            this->output_name = optarg;
            break;

        case 'k':
            this->params.k = parse<size_t>(optarg);
            break;
        case 'w':
            this->params.w_or_s = parse<size_t>(optarg);
            break;
        case 'c':
            this->params.use_syncmers = true;
            break;
        case 's':
            this->params.w_or_s = parse<size_t>(optarg);
            break;

        case 'W':
            this->params.use_weighted_minimizers = true;
            break;
        case OPT_THRESHOLD:
            this->params.threshold = parse<size_t>(optarg);
            break;
        case OPT_ITERATIONS:
            this->params.iterations = parse<size_t>(optarg);
            break;
        case OPT_FAST_COUNTING:
            this->params.space_efficient_counting = false;
            break;
        case OPT_SAVE_MEMORY:
            this->params.space_efficient_counting = true;
            break;
        case OPT_HASH_TABLE:
            this->params.hash_table_width = parse<size_t>(optarg);
            break;

        case 'z':
            this->zipcode_name = optarg;
            break;
        case 'p':
            this->params.progress = true;
            this->progress = true;
            break;
        case 't':
            this->threads = set_thread_count(logger, optarg);
            break;
        case OPT_NO_DIST:
            this->require_distance_index = false;
            break;
        case OPT_DESCRIBE:
            this->describe = true;
            this->output_name = optarg;
            this->require_distance_index = false;
            break;

        case 'h':
        case '?':
            help_minimizer(argv);
            std::exit(EXIT_FAILURE);
        default:
            std::abort();
        }
    }
    if (this->output_name.empty()) {
        logger.error() << "option --output-name is required" << std::endl;
    }
    if (!this->describe) {
        if (optind + 1 != argc) {
            help_minimizer(argv);
            std::exit(EXIT_FAILURE);
        }
        this->graph_name = argv[optind];
    }
    if (this->require_distance_index && this->distance_name.empty()) {
        logger.error() << "one of options --distance-index and --no-dist is required" << std::endl;
    }
}

//------------------------------------------------------------------------------

// Register subcommand
static vg::subcommand::Subcommand vg_minimizer("minimizer", "build a minimizer index or a syncmer index",
                                               vg::subcommand::TOOLKIT, main_minimizer);

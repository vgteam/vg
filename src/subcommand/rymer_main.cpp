/** \file rymer_main.cpp
 *
 * Defines the "vg rymer" subcommand, which builds the RYmer index that "vg safari" uses to seed
 * ancient DNA reads. See artifacts/rymer_index.hpp for what the index holds.
 */

#include "subcommand.hpp"

#include <algorithm>
#include <iostream>
#include <string>

#include <getopt.h>
#include <omp.h>

#include <gbwt/utils.h>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>

#include "../gbwtgraph_helper.hpp"
#include "../log.hpp"
#include "../utility.hpp"
#include "../../artifacts/rymer_index.hpp"
#include "../../artifacts/rymer_seed_rescue.hpp"

using namespace vg;
using namespace vg::subcommand;

namespace {

// Using more threads than this wastes CPU time without speeding up construction.
constexpr int RYMER_MAX_DEFAULT_THREADS = 16;

constexpr size_t RYMER_DEFAULT_K = 29;
constexpr size_t RYMER_DEFAULT_W = 11;
constexpr size_t RYMER_MIN_K = 10;

int rymer_default_threads() {
    return std::min(omp_get_max_threads(), RYMER_MAX_DEFAULT_THREADS);
}

void help_rymer(char** argv) {
    std::cerr << "usage: " << argv[0] << " rymer [options] -o graph.ry graph.gbz" << std::endl
              << std::endl
              << "Builds the (w, k)-RYmer index of the haplotypes in a GBZ graph, for vg safari." << std::endl
              << "Minimizers are chosen from the DNA sequence and stored under their RYmers, the" << std::endl
              << "k-mers written with R for purines (A, G) and Y for pyrimidines (C, T), so that" << std::endl
              << "reads with deamination damage still find them." << std::endl
              << std::endl
              << "required options:" << std::endl
              << "  -o, --output-name FILE    store the index in FILE" << std::endl
              << std::endl
              << "RYmer options:" << std::endl
              << "  -k, --kmer-length N       length of the k-mers in the index [" << RYMER_DEFAULT_K << "]" << std::endl
              << "                            (" << RYMER_MIN_K << " to " << gbwtgraph::Key64::KMER_MAX_LENGTH << ")" << std::endl
              << "  -w, --window-length N     choose the minimizer from a window of N k-mers [" << RYMER_DEFAULT_W << "]" << std::endl
              << std::endl
              << "other options:" << std::endl
              << "  -p, --progress            show progress information" << std::endl
              << "  -t, --threads N           use N threads for index construction [" << rymer_default_threads() << "]" << std::endl
              << "                            (using more than " << RYMER_MAX_DEFAULT_THREADS << " threads rarely helps)" << std::endl
              << "  -h, --help                print this help message" << std::endl
              << std::endl;
    print_safari_citation(std::cerr);
    std::cerr << std::endl;
}

}

int main_rymer(int argc, char** argv) {
    Logger logger("vg rymer");

    if (argc == 2) {
        help_rymer(argv);
        return 1;
    }

    std::string output_name;
    size_t k = RYMER_DEFAULT_K;
    size_t w = RYMER_DEFAULT_W;
    bool progress = false;
    int threads = rymer_default_threads();

    constexpr int OPT_DISTANCE_INDEX = 1000;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
            {"output-name", required_argument, nullptr, 'o'},
            {"kmer-length", required_argument, nullptr, 'k'},
            {"window-length", required_argument, nullptr, 'w'},
            {"progress", no_argument, nullptr, 'p'},
            {"threads", required_argument, nullptr, 't'},
            // Accepted so that commands written for the original SAFARI still run. The RYmer
            // index stores DNA k-mers in its payload, so it has no room for distance hints.
            {"distance-index", required_argument, nullptr, OPT_DISTANCE_INDEX},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, 0, nullptr, 0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "o:k:w:pt:d:h?", long_options, &option_index);
        if (c == -1) {
            break;
        }

        switch (c) {
        case 'o':
            output_name = optarg;
            break;
        case 'k':
            k = parse<size_t>(optarg);
            break;
        case 'w':
            w = parse<size_t>(optarg);
            break;
        case 'p':
            progress = true;
            break;
        case 't':
            threads = set_thread_count(logger, optarg);
            break;
        case 'd':
        case OPT_DISTANCE_INDEX:
            logger.warn() << "the RYmer index does not use a distance index; ignoring " << optarg << std::endl;
            break;
        case 'h':
        case '?':
        default:
            help_rymer(argv);
            return 1;
        }
    }

    if (output_name.empty()) {
        logger.error() << "option --output-name is required" << std::endl;
    }
    if (optind + 1 != argc) {
        help_rymer(argv);
        return 1;
    }
    if (k < RYMER_MIN_K || k > gbwtgraph::Key64::KMER_MAX_LENGTH) {
        logger.error() << "k-mer length must be from " << RYMER_MIN_K << " to "
                       << gbwtgraph::Key64::KMER_MAX_LENGTH << std::endl;
    }
    if (w == 0) {
        logger.error() << "window length must be at least 1" << std::endl;
    }
    std::string graph_name = argv[optind];
    omp_set_num_threads(threads);

    double start = gbwt::readTimer();

    gbwtgraph::GBZ gbz;
    load_gbz(gbz, graph_name, progress);

    gbwtgraph::DefaultMinimizerIndex index(k, w, RYMER_PAYLOAD_SIZE, false);
    if (progress) {
        logger.info() << "Building RYmer index with k = " << index.k() << ", w = " << index.w() << std::endl;
    }
    index_rymers(gbz, index);

    if (progress) {
        logger.info() << index.size() << " keys (" << index.unique_keys() << " unique)" << std::endl;
        logger.info() << "RYmer occurrences: " << index.number_of_values() << std::endl;
        logger.info() << "Load factor: " << index.load_factor() << std::endl;
    }

    save_minimizer(index, output_name, progress);

    if (progress) {
        double seconds = gbwt::readTimer() - start;
        logger.info() << "Time usage: " << seconds << " seconds" << std::endl;
        logger.info() << "Memory usage: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
    }

    return 0;
}

// Register subcommand
static Subcommand vg_rymer("rymer", "build an RYmer index for ancient DNA mapping with vg safari",
                           TOOLKIT, main_rymer);

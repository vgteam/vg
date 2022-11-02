/** \file haplotypes_main.cpp
 *
 * Defines the "vg haplotypes" subcommand, which will ultimately sample haplotypes.
 *
 * This is currently highly experimental.
 */

// FIXME tests

#include "subcommand.hpp"

#include "../recombinator.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

#include <getopt.h>
#include <omp.h>

using namespace vg;

//----------------------------------------------------------------------------

void help_haplotypes(char** argv) {
    // FIXME number of jobs, GBZ output, check that the generated GBZ is a subgraph
    std::cerr << "usage: " << argv[0] << " " << argv[1] << " [options] -o output.gbwt graph.gbz" << std::endl;
    std::cerr << std::endl;
    // FIXME description
    std::cerr << "Some experiments with haplotype sampling." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Required options:" << std::endl;
    std::cerr << "    -o, --output-name X     write the output GBWT to X" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Index files:" << std::endl;
    std::cerr << "    -d, --distance-index X  use this distance index (or guess from graph name)" << std::endl;
    std::cerr << "    -r, --r-index-name X    use this r-index (or guess from graph name)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "    -p, --progress          show progress information" << std::endl;
    std::cerr << std::endl;
}

// FIXME Should these be parameters?
constexpr size_t APPROXIMATE_JOBS = 32;
constexpr size_t NUM_HAPLOTYPES = 16;
constexpr size_t TARGET_DISTANCE = 10000;

//----------------------------------------------------------------------------

bool ends_with(const std::string& str, const std::string& suffix) {
    if (str.length() < suffix.length()) {
        return false;
    }
    return (str.substr(str.length() - suffix.length()) == suffix);
}

std::string get_name(const std::string& graph_name, const std::string& extension) {
    size_t length = graph_name.length();
    if (ends_with(graph_name, gbwtgraph::GBZ::EXTENSION)) {
        length -= gbwtgraph::GBZ::EXTENSION.length();
    }
    return graph_name.substr(0, length) + extension;
}

//----------------------------------------------------------------------------

int main_haplotypes(int argc, char** argv) {
    double start = gbwt::readTimer();
    if (argc < 5) {
        help_haplotypes(argv);
        return 1;
    }
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);

    // Parse options into these.
    std::string graph_name, r_index_name, distance_name, output_name;
    bool progress = false;

    // Process the arguments.
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "output-name", required_argument, 0, 'o' },
            { "distance-index", required_argument, 0, 'd' },
            { "r_index", required_argument, 0, 'r' },
            { "progress", no_argument, 0, 'p' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "o:d:r:ph", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'o':
            output_name = optarg;
            break;

        case 'd':
            distance_name = optarg;
            break;
        case 'r':
            r_index_name = optarg;
            break;

        case 'p':
            progress = true;
            break;

        case 'h':
        case '?':
            help_haplotypes(argv);
            return 1;
        default:
            std::abort();
        }
    }

    // Determine file names.
    if (optind + 1 != argc) {
        help_haplotypes(argv);
        return 1;
    }
    graph_name = argv[optind];
    if (output_name.empty()) {
        std::cerr << "Option --output-name is required" << std::endl;
        return 1;
    }
    if (r_index_name.empty()) {
        r_index_name = get_name(graph_name, gbwt::FastLocate::EXTENSION);
        if (progress) {
            std::cerr << "Guessing that r-index is " << r_index_name << std::endl;
        }
    }
    if (distance_name.empty()) {
        distance_name = get_name(graph_name, ".dist");
        if (progress) {
            std::cerr << "Guessing that distance index is " << distance_name << std::endl;
        }
    }

    // Load the indexes.
    double checkpoint = gbwt::readTimer();
    gbwtgraph::GBZ gbz;
    load_gbz(gbz, graph_name, progress);
    gbwt::FastLocate r_index;
    load_r_index(r_index, r_index_name, progress);
    r_index.setGBWT(gbz.index);
    SnarlDistanceIndex distance_index;
    if (progress) {
        std::cerr << "Loading distance index from " << distance_name << std::endl;
    }
    distance_index.deserialize(distance_name);
    if (progress) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Loaded the indexes in " << seconds << " seconds" << std::endl;
    }

    // Create a recombinator.
    Recombinator recombinator(gbz, r_index, distance_index, progress);

    if (progress) {
        std::cerr << "Processing chains using " << omp_get_max_threads() << " threads" << std::endl;
    }
    checkpoint = gbwt::readTimer();
    std::vector<gbwt::GBWT> indexes(recombinator.chains_by_job.size());
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < recombinator.chains_by_job.size(); job++) {
        const std::vector<gbwtgraph::TopLevelChain>& chains = recombinator.chains_by_job[job];
        // FIXME buffer size?
        gbwt::GBWTBuilder builder(sdsl::bits::length(gbz.index.sigma() - 1));
        builder.index.addMetadata();
        for (size_t chain = 0; chain < chains.size(); chain++) {
            recombinator.generate_haplotypes(chains[chain], builder);
        }
        builder.finish();
        indexes[job] = builder.index;
    }
    if (progress) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Processed the chains in " << seconds << " seconds" << std::endl;
    }

    if (progress) {
        std::cerr << "Merging the partial indexes" << std::endl;
    }
    checkpoint = gbwt::readTimer();
    gbwt::GBWT merged(indexes);
    if (progress) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Merged the indexes in " << seconds << " seconds" << std::endl;
    }
    save_gbwt(merged, output_name, progress);

    if (progress) {
        double seconds = gbwt::readTimer() - start;
        double gib = gbwt::inGigabytes(gbwt::memoryUsage());
        std::cerr << "Used " << seconds << " seconds, " << gib << " GiB" << std::endl;
    }
    return 0;
}

// FIXME description
static vg::subcommand::Subcommand vg_haplotypes("haplotypes", "haplotype sampling experiments", vg::subcommand::DEVELOPMENT, main_haplotypes);

//----------------------------------------------------------------------------


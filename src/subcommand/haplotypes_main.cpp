/** \file haplotypes_main.cpp
 *
 * Defines the "vg haplotypes" subcommand, which will ultimately sample haplotypes.
 *
 * This is currently highly experimental.
 */

// FIXME tests

#include "subcommand.hpp"

#include "../recombinator.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

#include <getopt.h>
#include <omp.h>

#include <gbwtgraph/index.h>

using namespace vg;

//----------------------------------------------------------------------------

constexpr size_t DEFAULT_MAX_JOBS = 14;
constexpr size_t APPROXIMATE_JOBS = 32;
constexpr size_t NUM_HAPLOTYPES = 16;
constexpr size_t TARGET_DISTANCE = 10000;

/*
  The same parameter is used both for the number of threads in minimizer index
  construction and the number of parallel GBWT construction jobs.

  In minimizer index construction, the usual rule of thumb is that 16 threds is
  enough. With too many threads, inserting the minimizers into the hash table
  becomes a bottleneck.

  In GBWT construction, we have foreground threads for generating the paths and
  background threads for building the GBWTs. The foreground threads are
  usually idle most of the time. The structure of the graph also limits the
  number of parallel jobs. In human graphs, there is no benefit in using more
  than 14 jobs.
*/
size_t get_default_jobs() {
    size_t jobs = std::round(0.85 * omp_get_max_threads());
    jobs = std::max(jobs, size_t(1));
    return std::min(jobs, DEFAULT_MAX_JOBS);
}

void help_haplotypes(char** argv) {
    std::cerr << "usage: " << argv[0] << " " << argv[1] << " [options] -o output.gbz graph.gbz" << std::endl;
    std::cerr << std::endl;
    // FIXME description
    std::cerr << "Some experiments with haplotype sampling." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Required options:" << std::endl;
    std::cerr << "    -o, --output-name X      write the output GBZ to X" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Index files:" << std::endl;
    std::cerr << "    -d, --distance-index X   use this distance index (default: <basename>.dist)" << std::endl;
    std::cerr << "    -m, --minimizer-index X  use this minimizer index (default: build the index)" << std::endl;
    std::cerr << "    -r, --r-index X          use this r-index (default: <basename>.ri)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "    -v, --verbosity N        set verbosity level (0 = none, 1 = basic, 2 = detailed, 3 = debug; default: 0)" << std::endl;
    std::cerr << "        --parallel-jobs N    number of parallel jobs (default: " << get_default_jobs() << ")" << std::endl;
    std::cerr << "        --validate           check that the output graph is a subgraph of the input" << std::endl;
    std::cerr << std::endl;
}

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

void validate_subgraph(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph, Recombinator::Verbosity verbosity);

//----------------------------------------------------------------------------

int main_haplotypes(int argc, char** argv) {
    double start = gbwt::readTimer();
    if (argc < 5) {
        help_haplotypes(argv);
        return 1;
    }
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);

    // Parse options into these.
    std::string graph_name, output_name;
    std::string distance_name, minimizer_name, r_index_name;
    Recombinator::Verbosity verbosity = Recombinator::verbosity_silent;
    size_t parallel_jobs = get_default_jobs();
    bool validate = false;

    constexpr int OPT_PARALLEL_JOBS = 1000;
    constexpr int OPT_VALIDATE = 1001;

    static struct option long_options[] =
    {
        { "output-name", required_argument, 0, 'o' },
        { "distance-index", required_argument, 0, 'd' },
        { "minimizer-index", required_argument, 0, 'm' },
        { "r-index", required_argument, 0, 'r' },
        { "verbosity", required_argument, 0, 'v' },
        { "parallel-jobs", required_argument, 0, OPT_PARALLEL_JOBS },
        { "validate", no_argument, 0,  OPT_VALIDATE },
        { 0, 0, 0, 0 }
    };

    // Process the arguments.
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        int option_index = 0;
        c = getopt_long(argc, argv, "o:d:m:r:v:h", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'o':
            output_name = optarg;
            break;

        case 'd':
            distance_name = optarg;
            break;
        case 'm':
            minimizer_name = optarg;
            break;
        case 'r':
            r_index_name = optarg;
            break;

        case 'v':
            {
                size_t level = parse<size_t>(optarg);
                if (level > Recombinator::verbosity_debug) {
                    std::cerr << "error: [vg haplotypes] invalid verbosity level: " << level << std::endl;
                    return 1;
                }
                verbosity = static_cast<Recombinator::Verbosity>(level);
            }
            break;
        case OPT_PARALLEL_JOBS:
            {
                size_t max_threads = omp_get_max_threads();
                parallel_jobs = parse<size_t>(optarg);
                if (parallel_jobs == 0 || parallel_jobs > max_threads) {
                    std::cerr << "error: [vg haplotypes] cannot run " << parallel_jobs << " jobs in parallel on this system" << std::endl;
                    return 1;
                }
            }
            break;
        case OPT_VALIDATE:
            validate = true;
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
    if (distance_name.empty()) {
        distance_name = get_name(graph_name, ".dist");
        if (verbosity >= Recombinator::verbosity_basic) {
            std::cerr << "Guessing that distance index is " << distance_name << std::endl;
        }
    }
    if (r_index_name.empty()) {
        r_index_name = get_name(graph_name, gbwt::FastLocate::EXTENSION);
        if (verbosity >= Recombinator::verbosity_basic) {
            std::cerr << "Guessing that r-index is " << r_index_name << std::endl;
        }
    }
    omp_set_num_threads(parallel_jobs);

    // Load/build the indexes.
    double checkpoint = gbwt::readTimer();
    gbwtgraph::GBZ gbz;
    load_gbz(gbz, graph_name, verbosity >= Recombinator::verbosity_basic);
    gbwt::FastLocate r_index;
    load_r_index(r_index, r_index_name, verbosity >= Recombinator::verbosity_basic);
    r_index.setGBWT(gbz.index);

    SnarlDistanceIndex distance_index;
    if (verbosity >= Recombinator::verbosity_basic) {
        std::cerr << "Loading distance index from " << distance_name << std::endl;
    }
    distance_index.deserialize(distance_name);
    gbwtgraph::DefaultMinimizerIndex minimizer_index;
    if (minimizer_name.empty()) {
        if (verbosity >= Recombinator::verbosity_basic) {
            std::cerr << "Building minimizer index" << std::endl;
        }
        gbwtgraph::index_haplotypes(gbz.graph, minimizer_index, [&](const pos_t& pos) -> gbwtgraph::payload_type {
            return MIPayload::encode(get_minimizer_distances(distance_index, pos));
        });
    } else {
        load_minimizer(minimizer_index, minimizer_name, verbosity >= Recombinator::verbosity_basic);
    }
    if (verbosity >= Recombinator::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Prepared the indexes in " << seconds << " seconds" << std::endl;
    }

    // Create a recombinator.
    Recombinator recombinator(gbz, r_index, distance_index, minimizer_index, verbosity);

    if (verbosity >= Recombinator::verbosity_basic) {
        std::cerr << "Running " << parallel_jobs << " jobs in parallel" << std::endl;
    }
    checkpoint = gbwt::readTimer();
    std::vector<gbwt::GBWT> indexes(recombinator.chains_by_job.size());
    Recombinator::Statistics statistics;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < recombinator.chains_by_job.size(); job++) {
        const std::vector<gbwtgraph::TopLevelChain>& chains = recombinator.chains_by_job[job];
        // FIXME buffer size?
        gbwt::GBWTBuilder builder(sdsl::bits::length(gbz.index.sigma() - 1));
        builder.index.addMetadata();
        Recombinator::Statistics job_statistics;
        for (size_t chain = 0; chain < chains.size(); chain++) {
            Recombinator::Statistics chain_statistics = recombinator.generate_haplotypes(chains[chain], builder);
            job_statistics.combine(chain_statistics);
        }
        builder.finish();
        indexes[job] = builder.index;
        #pragma omp critical
        {
            if (verbosity >= Recombinator::verbosity_detailed) {
                std::cerr << "Job " << job << ": "; job_statistics.print(std::cerr) << std::endl;
            }
            statistics.combine(job_statistics);
        }
    }
    if (verbosity >= Recombinator::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Processed "; statistics.print(std::cerr) << std::endl;
        std::cerr << "Finished the jobs in " << seconds << " seconds" << std::endl;
    }

    if (verbosity >= Recombinator::verbosity_basic) {
        std::cerr << "Merging the partial indexes" << std::endl;
    }
    checkpoint = gbwt::readTimer();
    gbwt::GBWT merged(indexes);
    if (verbosity >= Recombinator::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Merged the indexes in " << seconds << " seconds" << std::endl;
    }

    if (verbosity >= Recombinator::verbosity_basic) {
        std::cerr << "Building GBWTGraph" << std::endl;
    }
    checkpoint = gbwt::readTimer();
    gbwtgraph::GBWTGraph output_graph = gbz.graph.subgraph(merged);
    if (verbosity >= Recombinator::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Built the GBWTGraph in " << seconds << " seconds" << std::endl;
    }

    if (validate) {
        validate_subgraph(gbz.graph, output_graph, verbosity);
    }
    save_gbz(merged, output_graph, output_name, verbosity >= Recombinator::verbosity_basic);

    if (verbosity >= Recombinator::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        double gib = gbwt::inGigabytes(gbwt::memoryUsage());
        std::cerr << "Used " << seconds << " seconds, " << gib << " GiB" << std::endl;
    }
    return 0;
}

// FIXME description
static vg::subcommand::Subcommand vg_haplotypes("haplotypes", "haplotype sampling experiments", vg::subcommand::DEVELOPMENT, main_haplotypes);

//----------------------------------------------------------------------------

// Returns a GBWTGraph handle as a string (id, orientation).
std::string to_string(handle_t handle) {
    gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(handle);
    return std::string("(") + std::to_string(gbwt::Node::id(node)) + std::string(", ") + std::to_string(gbwt::Node::is_reverse(node)) + std::string(")");
}

void validate_nodes(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph) {
    nid_t last_node = 0;
    bool nodes_ok = subgraph.for_each_handle([&](const handle_t& handle) -> bool {
        last_node = subgraph.get_id(handle);
        return graph.has_node(last_node);
    });
    if (!nodes_ok) {
        std::cerr << "error: [vg haplotypes] subgraph node " << last_node << " is not in the original graph" << std::endl;
    }
}

void validate_edges(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph) {
    edge_t last_edge(gbwtgraph::GBWTGraph::node_to_handle(0), gbwtgraph::GBWTGraph::node_to_handle(0));
    bool edges_ok = subgraph.for_each_edge([&](const edge_t& edge) -> bool {
        last_edge = edge;
        return graph.has_edge(edge.first, edge.second);
    });
    if (!edges_ok) {
        std::cerr << "error: [vg haplotypes] subgraph edge from " << to_string(last_edge.first) << " to " << to_string(last_edge.second) << " is not in the original graph" << std::endl;
    }
}

void validate_subgraph(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph, Recombinator::Verbosity verbosity) {
    if (verbosity >= Recombinator::verbosity_basic) {
        std::cerr << "Validating the output subgraph" << std::endl;
    }
    double start = gbwt::readTimer();

    std::thread nodes(validate_nodes, std::cref(graph), std::cref(subgraph));
    std::thread edges(validate_edges, std::cref(graph), std::cref(subgraph));
    nodes.join();
    edges.join();

    if (verbosity >= Recombinator::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Validated the subgraph in " << seconds << " seconds" << std::endl;
    }
}

//----------------------------------------------------------------------------


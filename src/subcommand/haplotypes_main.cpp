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

constexpr size_t DEFAULT_MAX_THREADS = 16;

/*
  We take a single parameter for the number of threads and adjust that for
  various tasks.

  1) In minimizer index construction, the usual rule of thumb is not using more
  than 16 threads. Beyond that, inserting the minimizers into the hash table
  often becomes a bottleneck.

  2) Each GBWT construction job nominally uses two threads: one for generating
  the paths and another for building the GBWT. Generating the paths is much
  cheaper, and those threads only run a fraction of the time. Also, in human
  graphs, there is no benefit from using more than 14 construction jobs.

  3) GBWTGraph deserialization and construction use multiple threads for
  caching information about generic and reference paths.
*/
size_t haplotypes_default_threads() {
    size_t threads = omp_get_max_threads();
    threads = std::max(threads, size_t(1));
    return std::min(threads, DEFAULT_MAX_THREADS);
}

void help_haplotypes(char** argv) {
    std::cerr << "Usage: " << argv[0] << " " << argv[1] << " [options] (-k counts.kff -g output.gbz | -H output.hapl) graph.gbz" << std::endl;
    std::cerr << std::endl;
    // FIXME description
    std::cerr << "Some experiments with haplotype sampling." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Output options:" << std::endl;
    std::cerr << "    -g, --gbz-output X        write the output GBZ to X" << std::endl;
    std::cerr << "    -H, --haplotype-output X  write haplotype information to X" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Input files:" << std::endl;
    std::cerr << "    -d, --distance-index X    use this distance index (default: <basename>.dist)" << std::endl;
    std::cerr << "    -m, --minimizer-index X   use this minimizer index (default: build the index)" << std::endl;
    std::cerr << "    -r, --r-index X           use this r-index (default: <basename>.ri)" << std::endl;
    std::cerr << "    -i, --haplotype-input X   use this haplotype information (default: generate the information)" << std::endl;
    std::cerr << "    -k, --kmer-input X        use kmer counts from this KFF file (required for --gbz-output)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "    -v, --verbosity N         set verbosity level (0 = none, 1 = basic, 2 = detailed, 3 = debug; default: 0)" << std::endl;
    std::cerr << "    -t, --threads N           approximate number of threads (default: " << haplotypes_default_threads() << ")" << std::endl;
    std::cerr << "        --validate            check that the output graph is a subgraph of the input" << std::endl;
    std::cerr << std::endl;
}

//----------------------------------------------------------------------------

size_t threads_to_jobs(size_t threads) {
    size_t jobs = std::round(0.85 * threads);
    return std::max(jobs, size_t(1));
}

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

void validate_subgraph(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph, HaplotypePartitioner::Verbosity verbosity);

//----------------------------------------------------------------------------

int main_haplotypes(int argc, char** argv) {
    double start = gbwt::readTimer();
    if (argc < 5) {
        help_haplotypes(argv);
        return 1;
    }
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    size_t max_threads = omp_get_max_threads();

    // Parse options into these.
    std::string graph_name, gbz_output, haplotype_output;
    std::string distance_name, minimizer_name, r_index_name, haplotype_input, kmer_input;
    HaplotypePartitioner::Verbosity verbosity = HaplotypePartitioner::verbosity_silent;
    size_t threads = haplotypes_default_threads();
    bool validate = false;

    constexpr int OPT_VALIDATE = 1001;

    static struct option long_options[] =
    {
        { "gbz-output", required_argument, 0, 'g' },
        { "haplotype-output", required_argument, 0, 'H' },
        { "distance-index", required_argument, 0, 'd' },
        { "minimizer-index", required_argument, 0, 'm' },
        { "r-index", required_argument, 0, 'r' },
        { "haplotype-input", required_argument, 0, 'i' },
        { "kmer-input", required_argument, 0, 'k' },
        { "verbosity", required_argument, 0, 'v' },
        { "threads", required_argument, 0, 't' },
        { "validate", no_argument, 0,  OPT_VALIDATE },
        { 0, 0, 0, 0 }
    };

    // Process the arguments.
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        int option_index = 0;
        c = getopt_long(argc, argv, "g:H:d:m:r:i:k:v:t:h", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'g':
            gbz_output = optarg;
            break;
        case 'H':
            haplotype_output = optarg;
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
        case 'i':
            haplotype_input = optarg;
            break;
        case 'k':
            kmer_input = optarg;
            break;

        case 'v':
            {
                size_t level = parse<size_t>(optarg);
                if (level > HaplotypePartitioner::verbosity_debug) {
                    std::cerr << "error: [vg haplotypes] invalid verbosity level: " << level << std::endl;
                    return 1;
                }
                verbosity = static_cast<HaplotypePartitioner::Verbosity>(level);
            }
            break;
        case 't':
            {
                threads = parse<size_t>(optarg);
                if (threads == 0 || threads > max_threads) {
                    std::cerr << "error: [vg haplotypes] cannot run " << threads << " threads in parallel on this system" << std::endl;
                    return 1;
                }
                omp_set_num_threads(threads);
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
    if (!gbz_output.empty() && kmer_input.empty()) {
        std::cerr << "error: [vg haplotypes] --gbz-output requires --kmer-input" << std::endl;
        return 1;
    }
    if (gbz_output.empty() && haplotype_output.empty()) {
        std::cerr << "error: [vg haplotypes] at least one of --gbz-output and --haplotype-output is required" << std::endl;
        return 1;
    }

    // Load the graph.
    gbwtgraph::GBZ gbz;
    load_gbz(gbz, graph_name, verbosity >= HaplotypePartitioner::verbosity_basic);

    // Generate or load haplotype information.
    Haplotypes haplotypes;
    if (haplotype_input.empty()) {
        double checkpoint = gbwt::readTimer();
        if (verbosity >= HaplotypePartitioner::verbosity_basic) {
            std::cerr << "Generating haplotype information" << std::endl;
        }

        // Distance index.
        if (distance_name.empty()) {
            distance_name = get_name(graph_name, ".dist");
            if (verbosity >= HaplotypePartitioner::verbosity_basic) {
                std::cerr << "Guessing that distance index is " << distance_name << std::endl;
            }
        }
        SnarlDistanceIndex distance_index;
        if (verbosity >= HaplotypePartitioner::verbosity_basic) {
            std::cerr << "Loading distance index from " << distance_name << std::endl;
        }
        distance_index.deserialize(distance_name);

        // Minimizer index.
        gbwtgraph::DefaultMinimizerIndex minimizer_index;
        if (minimizer_name.empty()) {
            double minimizer = gbwt::readTimer();
            if (verbosity >= HaplotypePartitioner::verbosity_basic) {
                std::cerr << "Building minimizer index" << std::endl;
            }
            gbwtgraph::index_haplotypes(gbz.graph, minimizer_index, [&](const pos_t& pos) -> gbwtgraph::payload_type {
                return MIPayload::encode(get_minimizer_distances(distance_index, pos));
            });
            if (verbosity >= HaplotypePartitioner::verbosity_basic) {
                double seconds = gbwt::readTimer() - minimizer;
                std::cerr << "Built the minimizer index in " << seconds << " seconds" << std::endl;
            }
        } else {
            load_minimizer(minimizer_index, minimizer_name, verbosity >= HaplotypePartitioner::verbosity_basic);
        }

        // R-index.
        if (r_index_name.empty()) {
            r_index_name = get_name(graph_name, gbwt::FastLocate::EXTENSION);
            if (verbosity >= HaplotypePartitioner::verbosity_basic) {
                std::cerr << "Guessing that r-index is " << r_index_name << std::endl;
            }
        }
        gbwt::FastLocate r_index;
        load_r_index(r_index, r_index_name, verbosity >= HaplotypePartitioner::verbosity_basic);
        r_index.setGBWT(gbz.index);

        // Partition the haplotypes.
        HaplotypePartitioner partitioner(gbz, r_index, distance_index, minimizer_index, verbosity);
        haplotypes = partitioner.partition_haplotypes();
        if (verbosity >= HaplotypePartitioner::verbosity_basic) {
            double seconds = gbwt::readTimer() - checkpoint;
            std::cerr << "Generated haplotype information in " << seconds << " seconds" << std::endl;
        }
    } else {
        if (verbosity >= HaplotypePartitioner::verbosity_basic) {
            std::cerr << "Loading haplotype information from " << haplotype_input << std::endl;
        }
        sdsl::simple_sds::load_from(haplotypes, haplotype_input);
    }

    // Save haplotype information if necessary.
    if (!haplotype_output.empty()) {
        if (verbosity >= HaplotypePartitioner::verbosity_basic) {
            std::cerr << "Writing haplotype information to " << haplotype_output << std::endl;
        }
        sdsl::simple_sds::serialize_to(haplotypes, haplotype_output);
    }

    if (gbz_output.empty()) {
        if (verbosity >= HaplotypePartitioner::verbosity_basic) {
            double seconds = gbwt::readTimer() - start;
            double gib = gbwt::inGigabytes(gbwt::memoryUsage());
            std::cerr << "Used " << seconds << " seconds, " << gib << " GiB" << std::endl;
        }
        return 0;
    }

    // Generate haplotypes.
    // TODO: Set number of jobs as a parameter?
    omp_set_num_threads(threads_to_jobs(threads));
    Recombinator recombinator(gbz, verbosity);
    gbwt::GBWT merged = recombinator.generate_haplotypes(haplotypes);
    omp_set_num_threads(threads); // Restore the number of threads.

    // Build GBWTGraph.
    if (verbosity >= HaplotypePartitioner::verbosity_basic) {
        std::cerr << "Building GBWTGraph" << std::endl;
    }
    double checkpoint = gbwt::readTimer();
    gbwtgraph::GBWTGraph output_graph = gbz.graph.subgraph(merged);
    if (verbosity >= HaplotypePartitioner::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Built the GBWTGraph in " << seconds << " seconds" << std::endl;
    }

    // Validate and serialize the graph.
    if (validate) {
        validate_subgraph(gbz.graph, output_graph, verbosity);
    }
    save_gbz(merged, output_graph, gbz_output, verbosity >= HaplotypePartitioner::verbosity_basic);

    if (verbosity >= HaplotypePartitioner::verbosity_basic) {
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

void validate_subgraph(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph, HaplotypePartitioner::Verbosity verbosity) {
    if (verbosity >= HaplotypePartitioner::verbosity_basic) {
        std::cerr << "Validating the output subgraph" << std::endl;
    }
    double start = gbwt::readTimer();

    std::thread nodes(validate_nodes, std::cref(graph), std::cref(subgraph));
    std::thread edges(validate_edges, std::cref(graph), std::cref(subgraph));
    nodes.join();
    edges.join();

    if (verbosity >= HaplotypePartitioner::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Validated the subgraph in " << seconds << " seconds" << std::endl;
    }
}

//----------------------------------------------------------------------------


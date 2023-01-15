/** \file haplotypes_main.cpp
 *
 * Defines the "vg haplotypes" subcommand, which will ultimately sample haplotypes.
 *
 * This is currently highly experimental.
 */

#include "subcommand.hpp"

#include "../recombinator.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <thread>
#include <unordered_map>
#include <unordered_set>
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

size_t haplotypes_default_k() {
    return Haplotypes::Header::DEFAULT_K;
}

size_t haplotypes_default_w() {
    return gbwtgraph::Key64::WINDOW_LENGTH;
}

void help_haplotypes(char** argv) {
    std::cerr << "Usage: " << argv[0] << " " << argv[1] << " [options] (-k counts.kff -g output.gbz | -H output.hapl) graph.gbz" << std::endl;
    std::cerr << std::endl;
    // FIXME description
    std::cerr << "Some experiments with haplotype sampling." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Output files:" << std::endl;
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
    // TODO: Expose other computational parameters?
    std::cerr << "Computational parameters:" << std::endl;
    std::cerr << "        --kmer-length N       kmer length for building the minimizer index (default: " << haplotypes_default_k() << ")" << std::endl;
    std::cerr << "        --window-length N     window length for building the minimizer index (default: " << haplotypes_default_w() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "    -v, --verbosity N         verbosity level (0 = none, 1 = basic, 2 = detailed, 3 = debug; default: 0)" << std::endl;
    std::cerr << "    -t, --threads N           approximate number of threads (default: " << haplotypes_default_threads() << ")" << std::endl;
    std::cerr << "        --validate            validate the generated information (may be slow)" << std::endl;
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

void validate_haplotypes(const Haplotypes& haplotypes,
                         const gbwtgraph::GBWTGraph& graph,
                         const gbwt::FastLocate& r_index,
                         const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
                         size_t expected_chains,
                         HaplotypePartitioner::Verbosity verbosity);

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
    size_t k = haplotypes_default_k(), w = haplotypes_default_w();
    HaplotypePartitioner::Verbosity verbosity = HaplotypePartitioner::verbosity_silent;
    size_t threads = haplotypes_default_threads();
    bool validate = false;

    constexpr int OPT_KMER_LENGTH = 1200;
    constexpr int OPT_WINDOW_LENGTH = 1201;
    constexpr int OPT_VALIDATE = 1300;

    static struct option long_options[] =
    {
        { "gbz-output", required_argument, 0, 'g' },
        { "haplotype-output", required_argument, 0, 'H' },
        { "distance-index", required_argument, 0, 'd' },
        { "minimizer-index", required_argument, 0, 'm' },
        { "r-index", required_argument, 0, 'r' },
        { "haplotype-input", required_argument, 0, 'i' },
        { "kmer-input", required_argument, 0, 'k' },
        { "kmer-length", required_argument, 0, OPT_KMER_LENGTH },
        { "window-length", required_argument, 0, OPT_WINDOW_LENGTH },
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

        case OPT_KMER_LENGTH:
            k = parse<size_t>(optarg);
            if (k == 0 || k > gbwtgraph::Key64::KMER_MAX_LENGTH) {
                std::cerr << "error: [vg haplotypes] kmer length must be between 1 and " << gbwtgraph::Key64::KMER_MAX_LENGTH << std::endl;
                return 1;
            }
            break;
        case OPT_WINDOW_LENGTH:
            w = parse<size_t>(optarg);
            if (w == 0) {
                std::cerr << "error: [vg haplotypes] window length cannot be 0" << std::endl;
                return 1;
            }
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
            threads = parse<size_t>(optarg);
            if (threads == 0 || threads > max_threads) {
                std::cerr << "error: [vg haplotypes] cannot run " << threads << " threads in parallel on this system" << std::endl;
                return 1;
            }
            omp_set_num_threads(threads);
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
        size_t expected_chains = 0;
        distance_index.for_each_child(distance_index.get_root(), [&](const handlegraph::net_handle_t&) {
            expected_chains++;
        });

        // Minimizer index.
        gbwtgraph::DefaultMinimizerIndex minimizer_index(k, w, false);
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
        HaplotypePartitioner::Parameters parameters;
        haplotypes = partitioner.partition_haplotypes(parameters);
        if (verbosity >= HaplotypePartitioner::verbosity_basic) {
            double seconds = gbwt::readTimer() - checkpoint;
            std::cerr << "Generated haplotype information in " << seconds << " seconds" << std::endl;
        }

        // Validate the haplotypes.
        if (validate) {
            validate_haplotypes(haplotypes, gbz.graph, r_index, minimizer_index, expected_chains, verbosity);
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
    omp_set_num_threads(threads_to_jobs(threads));
    Recombinator recombinator(gbz, verbosity);
    Recombinator::Parameters recombinator_parameters;
    gbwt::GBWT merged = recombinator.generate_haplotypes(haplotypes, recombinator_parameters);
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
        // FIXME validate haplotypes
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

void validate_error(const std::string& header, const std::string& message) {
    std::cerr << "error: [vg haplotypes]";
    if (!header.empty()) {
        std::cerr << header << ": ";
    }
    std::cerr << message << std::endl;
    std::exit(EXIT_FAILURE);
}

template<typename T>
std::string expected_got(T expected, T got) {
    return "expected " + std::to_string(expected) + ", got " + std::to_string(got);
}

template<typename T>
std::string pair_to_string(std::pair<T, T> value) {
    return "(" + std::to_string(value.first) + ", " + std::to_string(value.second) + ")";
}

template<typename T>
std::string expected_got_pair(std::pair<T, T> expected, std::pair<T, T> got) {
    return "expected " + pair_to_string(expected) + ", got " + pair_to_string(got);
}

void validate_error_chain(size_t chain_id, const std::string& message) {
    validate_error("chain " + std::to_string(chain_id), message);
}

void validate_error_subchain(size_t chain_id, size_t subchain_id, const std::string& message) {
    validate_error("chain " + std::to_string(chain_id) + ", subchain " + std::to_string(subchain_id), message);
}

void validate_error_sequence(size_t chain_id, size_t subchain_id, size_t sequence_id, const std::string& message) {
    std::string header = "chain " + std::to_string(chain_id) + ", subchain " + std::to_string(subchain_id) + ", sequence " + std::to_string(sequence_id);
    validate_error(header, message);
}

std::string validate_unary_path(const HandleGraph& graph, handle_t from, handle_t to) {
    std::unordered_set<handle_t> visited;
    handle_t curr = from;
    while (curr != to) {
        if (visited.find(curr) != visited.end()) {
            return "incoming path contains a cycle";
        }
        visited.insert(curr);
        handle_t successor = empty_gbwtgraph_handle();
        size_t successors = 0;
        graph.follow_edges(curr, false, [&](const handle_t& next) {
            successor = next;
            successors++;
        });
        if (successors != 1) {
            return "incoming path is not unary";
        }
        curr = successor;
    }
    return "";
}

// Returns true if the path from (start, offset) reaches end without revisiting start.
bool trace_path(const gbwt::GBWT& index, gbwt::node_type start, gbwt::size_type offset, gbwt::node_type end) {
    gbwt::edge_type pos(start, offset);
    while (pos.first != end) {
        pos = index.LF(pos);
        if (pos.first == gbwt::ENDMARKER || pos.first == start) {
            return false;
        }
    }
    return true;
}

// Returns the given haplotype over the given subchain.
std::string get_haplotype(const gbwtgraph::GBWTGraph& graph, Haplotypes::sequence_type sequence,
                          gbwt::node_type from, gbwt::node_type to, size_t k) {
    std::string haplotype;
    gbwt::edge_type pos;

    // Initial node with three cases (from start, suffix of a long `from`, short `from`).
    if (from == gbwt::ENDMARKER) {
        pos = graph.index->start(sequence.first);
        gbwtgraph::view_type view = graph.get_sequence_view(gbwtgraph::GBWTGraph::node_to_handle(pos.first));
        haplotype.append(view.first, view.second);
    } else {
        pos = gbwt::edge_type(from, sequence.second);
        gbwtgraph::view_type view = graph.get_sequence_view(gbwtgraph::GBWTGraph::node_to_handle(pos.first));
        if (view.second >= k) {
            haplotype.append(view.first + view.second - (k - 1), k - 1);
        } else {
            haplotype.append(view.first, view.second);
        }
    }

    while (true) {
        pos = graph.index->LF(pos);
        if (pos.first == gbwt::ENDMARKER) {
            break;
        }
        gbwtgraph::view_type view = graph.get_sequence_view(gbwtgraph::GBWTGraph::node_to_handle(pos.first));
        if (pos.first == to) {
            haplotype.append(view.first, std::min(view.second, k - 1));
        } else {
            haplotype.append(view.first, view.second);
        }
    }

    return haplotype;
}

void validate_chain(const Haplotypes::TopLevelChain& chain,
                    const gbwtgraph::GBWTGraph& graph,
                    const gbwt::FastLocate& r_index,
                    const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
                    size_t chain_id,
                    HaplotypePartitioner::Verbosity verbosity) {
    if (chain.offset != chain_id) {
        validate_error_chain(chain_id, "stored id is " + std::to_string(chain.offset));
    }
    if (chain.subchains.empty()) {
        validate_error_chain(chain_id, "the chain is empty");
    }

    const Haplotypes::Subchain* prev = nullptr;
    for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
        const Haplotypes::Subchain& subchain = chain.subchains[subchain_id];

        // Check that the subchain is of an appropriate type.
        switch (subchain.type) {
        case Haplotypes::Subchain::normal:
            break;
        case Haplotypes::Subchain::prefix:
            if (subchain_id > 0 && prev->type != Haplotypes::Subchain::suffix) {
                validate_error_subchain(chain_id, subchain_id, "a prefix inside a fragment");
            }
            break;
        case Haplotypes::Subchain::suffix:
            break;
        case Haplotypes::Subchain::full_haplotype:
            if (chain.subchains.size() != 1) {
                validate_error_subchain(chain_id, subchain_id, "full haplotypes in a nontrivial chain");
            }
            break;
        }

        // Check that the boundary nodes have been defined.
        if (subchain.has_start() && subchain.start == gbwt::ENDMARKER) {
            validate_error_subchain(chain_id, subchain_id, "missing start node");
        }
        if (subchain.has_end() && subchain.end == gbwt::ENDMARKER) {
            validate_error_subchain(chain_id, subchain_id, "missing end node");
        }

        // Check that the kmer presence bitvector is of appropriate length.
        size_t total_kmers = subchain.sequences.size() * subchain.kmers.size();
        if (subchain.kmers_present.size() != total_kmers) {
            std::string message = expected_got(total_kmers, subchain.kmers_present.size()) + " kmer occurrences";
            validate_error_subchain(chain_id, subchain_id, message);
        }

        // Check that there is a unary path from the previous subchain if the
        // appropriate boundary nodes are present.
        if (subchain_id > 0 && prev->has_end() && subchain.has_start()) {
            std::string message = validate_unary_path(graph, gbwtgraph::GBWTGraph::node_to_handle(prev->end), gbwtgraph::GBWTGraph::node_to_handle(subchain.start));
            if (!message.empty()) {
                validate_error_subchain(chain_id, subchain_id, message);
            }
        }

        // Sequences: normal subchains.
        if (subchain.type == Haplotypes::Subchain::normal) {
            std::vector<gbwt::size_type> da = r_index.decompressDA(subchain.start);
            std::unordered_set<Haplotypes::sequence_type> selected;
            for (size_t i = 0; i < da.size(); i++) {
                if (trace_path(*(graph.index), subchain.start, i, subchain.end)) {
                    selected.insert(Haplotypes::sequence_type(da[i], i));
                }
            }
            if (subchain.sequences.size() != selected.size()) {
                std::string message = expected_got(selected.size(), subchain.sequences.size()) + " sequences (normal)";
                validate_error_subchain(chain_id, subchain_id, message);
            }
            for (size_t i = 0; i < subchain.sequences.size(); i++) {
                if (selected.find(subchain.sequences[i]) == selected.end()) {
                    validate_error_sequence(chain_id, subchain_id, i, "invalid value " + pair_to_string(subchain.sequences[i]));
                }
            }
        }

        // Sequences: prefixes and suffixes.
        if (subchain.type == Haplotypes::Subchain::prefix || subchain.type == Haplotypes::Subchain::suffix) {
            gbwt::node_type node = (subchain.has_start() ? subchain.start : subchain.end);
            std::vector<gbwt::size_type> da = r_index.decompressDA(node);
            if (subchain.sequences.size() != da.size()) {
                std::string message = expected_got(da.size(), subchain.sequences.size()) + " sequences (prefix / suffix)";
                validate_error_subchain(chain_id, subchain_id, message);
            }
            for (size_t i = 0; i < subchain.sequences.size(); i++) {
                Haplotypes::sequence_type expected(da[i], i);
                if (subchain.sequences[i] != expected) {
                    std::string message = expected_got_pair(expected, subchain.sequences[i]);
                    validate_error_sequence(chain_id, subchain_id, i, message);
                }
            }
        }

        // Sequences: full haplotypes.
        if (subchain.type == Haplotypes::Subchain::full_haplotype) {
            if (subchain.sequences.empty()) {
                validate_error_subchain(chain_id, subchain_id, "full haplotypes without sequences");
            }
        }

        // Kmers.
        if (subchain.type != Haplotypes::Subchain::full_haplotype) {
            std::unordered_set<Haplotypes::Subchain::kmer_type> all_kmers;
            for (size_t i = 0; i < subchain.kmers.size(); i++) {
                all_kmers.insert(subchain.kmers[i]);
            }
            if (all_kmers.size() != subchain.kmers.size()) {
                std::string message = expected_got(subchain.kmers.size(), all_kmers.size()) + " kmers";
                validate_error_subchain(chain_id, subchain_id, message);
            }
            for (size_t i = 0; i < subchain.sequences.size(); i++) {
                std::string haplotype = get_haplotype(graph, subchain.sequences[i], subchain.start, subchain.end, minimizer_index.k());
                auto minimizers = minimizer_index.minimizers(haplotype);
                std::unordered_set<Haplotypes::Subchain::kmer_type> kmers_present;
                for (auto& minimizer : minimizers) {
                    if (minimizer_index.count(minimizer) == 1) {
                        kmers_present.insert(minimizer.key.get_key());
                    }
                }
                size_t found_kmers = 0;
                for (size_t j = 0, offset = i * subchain.kmers.size(); j < subchain.kmers.size(); j++, offset++) {
                    if (subchain.kmers_present[offset]) {
                        if (kmers_present.find(subchain.kmers[j]) == kmers_present.end()) {
                            std::string message = "kmer " + std::to_string(j) + " not present in the haplotype";
                            validate_error_sequence(chain_id, subchain_id, i, message);
                        }
                        found_kmers++;
                    } else {
                        if (kmers_present.find(subchain.kmers[j]) != kmers_present.end()) {
                            std::string message = "kmer " + std::to_string(j) + " is present in the haplotype";
                            validate_error_sequence(chain_id, subchain_id, i, message);
                        }
                    }
                }
                if (found_kmers != kmers_present.size()) {
                    std::string message = expected_got(kmers_present.size(), found_kmers) + " kmers";
                    validate_error_sequence(chain_id, subchain_id, i, message);
                }
            }
        }

        prev = &subchain;
    }
}

void validate_haplotypes(const Haplotypes& haplotypes,
                         const gbwtgraph::GBWTGraph& graph,
                         const gbwt::FastLocate& r_index,
                         const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
                         size_t expected_chains,
                         HaplotypePartitioner::Verbosity verbosity) {
    if (verbosity >= HaplotypePartitioner::verbosity_basic) {
        std::cerr << "Validating the haplotype information" << std::endl;
    }
    double start = gbwt::readTimer();

    // Header information.
    if (haplotypes.k() != minimizer_index.k()) {
        validate_error("k-mer length", expected_got(minimizer_index.k(), haplotypes.k()));
    }
    if (haplotypes.components() != expected_chains) {
        validate_error("graph components", expected_got(expected_chains, haplotypes.components()));
    }
    if (haplotypes.components() != haplotypes.chains.size()) {
        validate_error("top-level chains", expected_got(haplotypes.components(), haplotypes.chains.size()));
    }
    std::vector<size_t> chains_per_job(haplotypes.jobs(), 0);
    for (size_t chain = 0; chain < haplotypes.components(); chain++) {
        size_t job_id = haplotypes.chains[chain].job_id;
        if (job_id >= haplotypes.jobs()) {
            validate_error_chain(chain, "job id " + std::to_string(job_id) + " >= " + std::to_string(haplotypes.jobs()));
        }
        chains_per_job[job_id]++;
    }
    for (size_t job_id = 0; job_id < chains_per_job.size(); job_id++) {
        if (chains_per_job[job_id] == 0) {
            validate_error("", "job " + std::to_string(job_id) + " is empty");
        }
    }

    // Haplotype information is valid
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t chain = 0; chain < haplotypes.components(); chain++) {
        validate_chain(haplotypes.chains[chain], graph, r_index, minimizer_index, chain, verbosity);
    }

    // Kmers are globally unique.
    std::unordered_map<Haplotypes::Subchain::kmer_type, std::pair<size_t, size_t>> kmers;
    for (size_t chain_id = 0; chain_id < haplotypes.components(); chain_id++) {
        const Haplotypes::TopLevelChain& chain = haplotypes.chains[chain_id];
        for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
            const Haplotypes::Subchain& subchain = chain.subchains[subchain_id];
            for (size_t i = 0; i < subchain.kmers.size(); i++) {
                auto iter = kmers.find(subchain.kmers[i]);
                if (iter != kmers.end()) {
                    std::string message = "kmer " + std::to_string(i) + " also found in chain " + std::to_string(iter->second.first) + ", subchain " + std::to_string(iter->second.second);
                    validate_error_subchain(chain_id, subchain_id, message);
                }
                kmers[subchain.kmers[i]] = { chain_id, subchain_id };
            }
        }
    }
    kmers.clear();

    if (verbosity >= HaplotypePartitioner::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Validated the haplotype information in " << seconds << " seconds" << std::endl;
    }
}

//----------------------------------------------------------------------------

void validate_nodes(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph) {
    nid_t last_node = 0;
    bool nodes_ok = subgraph.for_each_handle([&](const handle_t& handle) -> bool {
        last_node = subgraph.get_id(handle);
        return graph.has_node(last_node);
    });
    if (!nodes_ok) {
        validate_error("", "invalid node " + std::to_string(last_node));
    }
}

void validate_edges(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph) {
    edge_t last_edge(gbwtgraph::GBWTGraph::node_to_handle(0), gbwtgraph::GBWTGraph::node_to_handle(0));
    bool edges_ok = subgraph.for_each_edge([&](const edge_t& edge) -> bool {
        last_edge = edge;
        return graph.has_edge(edge.first, edge.second);
    });
    if (!edges_ok) {
        validate_error("", "invalid edge " + to_string_gbwtgraph(last_edge.first) + " to " + to_string_gbwtgraph(last_edge.second));
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


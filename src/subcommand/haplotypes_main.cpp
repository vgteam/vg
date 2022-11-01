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

size_t get_distance(handle_t from, handle_t to, const HandleGraph& graph, const SnarlDistanceIndex& distance_index) {
    nid_t from_id = graph.get_id(from);
    bool from_is_reverse = graph.get_is_reverse(from);
    size_t from_offset = graph.get_length(from) - 1;
    nid_t to_id = graph.get_id(to);
    bool to_is_reverse = graph.get_is_reverse(to);
    // FIXME these are from the same top-level chain, so distance computation should be much faster
    size_t distance = distance_index.minimum_distance(from_id, from_is_reverse, from_offset, to_id, to_is_reverse, 0, false, &graph);
    return distance;
}

handle_t empty_handle() {
    return gbwtgraph::GBWTGraph::node_to_handle(0);
}

std::pair<handle_t, handle_t> empty_subchain() {
    return { empty_handle(), empty_handle() };
}

std::string to_string(handle_t handle) {
    gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(handle);
    return std::string("(") + std::to_string(gbwt::Node::id(node)) + std::string(", ") + std::to_string(gbwt::Node::is_reverse(node)) + std::string(")");
}

// (SA[i], i)
std::vector<std::pair<gbwt::size_type, gbwt::size_type>> get_paths(handle_t handle, const gbwt::FastLocate& r_index) {
    std::vector<gbwt::size_type> sa = r_index.decompressSA(gbwtgraph::GBWTGraph::handle_to_node(handle));
    std::vector<std::pair<gbwt::size_type, gbwt::size_type>> result;
    result.reserve(sa.size());
    for (size_t i = 0; i < sa.size(); i++) {
        result.push_back({ sa[i], i });
    }
    std::sort(result.begin(), result.end(), [&](std::pair<gbwt::size_type, gbwt::size_type> a, std::pair<gbwt::size_type, gbwt::size_type> b) -> bool {
        return (r_index.seqId(a.first) < r_index.seqId(b.first));
    });

    gbwt::size_type prev = std::numeric_limits<gbwt::size_type>::max();
    gbwt::size_type last_warning = std::numeric_limits<gbwt::size_type>::max();
    for (size_t i = 0; i < result.size(); i++) {
        gbwt::size_type seq_id = r_index.seqId(result[i].first);
        if (seq_id == prev && seq_id != last_warning) {
            #pragma omp critical
            {
                std::cerr << "warning: [vg haplotypes] GBWT sequence " << seq_id << " visits node " << to_string(handle) << " multiple times" << std::endl;
            }
            last_warning = seq_id;
        }
        prev = seq_id;
    }

    return result;
}

// FIXME handle haplotypes that visit the same node multiple times properly
// Returns (SA[i], i) at `from` for sequences that reach `to`.
std::vector<std::pair<gbwt::size_type, gbwt::size_type>> get_sequences(handle_t from, handle_t to, const gbwt::FastLocate& r_index) {
    auto from_paths = get_paths(from, r_index);
    auto to_paths = get_paths(to, r_index);

    auto from_iter = from_paths.begin();
    auto to_iter = to_paths.begin();
    std::vector<std::pair<gbwt::size_type, gbwt::size_type>> result;
    while (from_iter != from_paths.end() && to_iter != to_paths.end()) {
        gbwt::size_type from_id = r_index.seqId(from_iter->first);
        gbwt::size_type to_id = r_index.seqId(to_iter->first);
        if (from_id == to_id) {
            if (r_index.seqOffset(from_iter->first) >= r_index.seqOffset(to_iter->first)) {
                result.push_back(*from_iter);
            }
            ++from_iter; ++to_iter;
        } else if (from_id < to_id) {
            ++from_iter;
        } else if (from_id > to_id) {
            ++to_iter;
        }
    }

    return result;
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

    // FIXME are all chains in the right orientation?
    // FIXME if there are fragments without crossable subchains, take an entire haplotype
    // FIXME choose a different haplotype for a prefix/suffix?s
    if (progress) {
        std::cerr << "Processing chains using " << omp_get_max_threads() << " threads" << std::endl;
    }
    checkpoint = gbwt::readTimer();
    std::vector<gbwt::GBWT> indexes(recombinator.chains_by_job.size());
    size_t total_snarls = 0, total_skipped = 0, total_subchains = 0, total_generated = 0, total_took = 0;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < recombinator.chains_by_job.size(); job++) {
        const std::vector<gbwtgraph::TopLevelChain>& chains = recombinator.chains_by_job[job];
        // FIXME buffer size?
        gbwt::GBWTBuilder builder(sdsl::bits::length(gbz.index.sigma() - 1));
        builder.index.addMetadata();

        for (size_t chain = 0; chain < chains.size(); chain++) {
            std::vector<std::pair<handle_t, handle_t>> snarls;
            handle_t snarl_start = empty_handle();
            bool has_start = false;
            bool was_snarl = false;
            size_t step = 0;
            size_t skipped_snarls = 0;
            bool generated_haplotypes = false;

            net_handle_t curr = distance_index.get_bound(chains[chain].chain, false, true);
            net_handle_t chain_end = distance_index.get_bound(chains[chain].chain, true, false);
            while (curr != chain_end) {
                if (distance_index.is_node(curr)) {
                    handle_t handle = distance_index.get_handle(curr, &gbz.graph);
                    if (was_snarl && has_start) {
                        size_t distance = get_distance(snarl_start, handle, gbz.graph, distance_index);
                        if (distance < std::numeric_limits<size_t>::max()) {
                            snarls.push_back({ snarl_start, handle });
                        } else {
                            // Use empty subchains to mark possible haplotype breakpoints.                        
                            snarls.push_back(empty_subchain());
                            skipped_snarls++;
                        }
                    }
                    snarl_start = handle;
                    has_start = true;
                    was_snarl = false;
                } else if (distance_index.is_snarl(curr)) {
                    was_snarl = true;
                }

                net_handle_t next;
                size_t successors = 0;
                distance_index.follow_net_edges(curr, &gbz.graph, false, [&](const net_handle_t& child) {
                    successors++;
                    next = child;
                });
                if (successors != 1) {
                    #pragma omp critical
                    {
                        std::cerr << "Chain " << chains[chain].offset << ": found " << successors << " successors at step " << step << std::endl;
                    }
                    break;
                }
                curr = next;
                step++;
            }

            // Combine snarls into subchains.
            std::vector<std::pair<handle_t, handle_t>> subchains;
            size_t total_distance = 0;
            size_t max_distance = 0;
            size_t head = 0;
            while (head < snarls.size()) {
                if (snarls[head] == empty_subchain()) {
                    subchains.push_back(empty_subchain());
                    head++;
                    continue;
                }
                size_t distance = get_distance(snarls[head].first, snarls[head].second, gbz.graph, distance_index);
                size_t tail = head;
                while (tail + 1 < snarls.size()) {
                    if (snarls[tail + 1] == empty_subchain()) {
                        break;
                    }
                    size_t candidate = get_distance(snarls[head].first, snarls[tail + 1].second, gbz.graph, distance_index);
                    if (candidate <= TARGET_DISTANCE) {
                        distance = candidate;
                        tail++;
                    } else {
                        break;
                    }
                }
                total_distance += distance;
                max_distance = std::max(distance, max_distance);
                subchains.push_back({ snarls[head].first, snarls[tail].second });
                head = tail + 1;
            }

            // Generate the haplotypes.
            // FIXME include reference paths?
            std::vector<Recombinator::Haplotype> haplotypes;
            for (size_t i = 0; i < NUM_HAPLOTYPES; i++) {
                haplotypes.push_back({ chains[chain].offset, i, 0, gbwt::invalid_edge(), {} });
            }
            bool have_haplotypes = false;
            for (size_t i = 0; i < subchains.size(); i++) {
                if (subchains[i] == empty_subchain()) {
                    if (have_haplotypes) {
                        for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                            haplotypes[haplotype].finish(recombinator, builder);
                        }
                    }
                    have_haplotypes = false;
                    continue;
                }
                auto sequences = get_sequences(subchains[i].first, subchains[i].second, r_index);
                if (sequences.empty()) {
                    if (have_haplotypes) {
                        for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                            haplotypes[haplotype].finish(recombinator, builder);
                        }
                    }
                    have_haplotypes = false;
                    continue;
                }
                for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                    haplotypes[haplotype].extend(sequences[haplotype % sequences.size()], subchains[i], recombinator);
                }
                have_haplotypes = true;
                generated_haplotypes = true;
            }
            if (have_haplotypes) {
                for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                    haplotypes[haplotype].finish(recombinator, builder);
                }
                have_haplotypes = false;
            }

            // Take entire sequences if we could not generate any haplotypes.
            if (!generated_haplotypes) {
                gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(chains[chain].handle);
                auto sequences = r_index.decompressDA(node);
                for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                    haplotypes[haplotype].take(sequences[haplotype % sequences.size()], recombinator, builder);
                }
            }

            #pragma omp critical
            {
                size_t real_snarls = snarls.size() - skipped_snarls;
                size_t real_subchains = subchains.size() - skipped_snarls;
                if (progress && generated_haplotypes) {
                    std::cerr << "Chain " << chains[chain].offset << ":" << std::endl;
                    std::cerr << "  Selected " << real_snarls << " snarls and skipped " << skipped_snarls << " snarls" << std::endl;
                    std::cerr << "  Combined the snarls into " << real_subchains << " subchains" << std::endl;
                    std::cerr << "  Average distance " << (static_cast<double>(total_distance) / real_subchains) << ", max distance " << max_distance << std::endl;
                    std::cerr << "  Generated the haplotypes in " << haplotypes.front().fragment << " fragments" << std::endl;
                }
                total_snarls += real_snarls;
                total_skipped += skipped_snarls;
                total_subchains += real_subchains;
                if (generated_haplotypes) { total_generated++; }
                else { total_took++; }
            }
        }
        builder.finish();
        indexes[job] = builder.index;
    }
    if (progress) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Processed the chains in " << seconds << " seconds" << std::endl;
        std::cerr << "Combined " << total_snarls << " snarls into " << total_subchains << " subchains and skipped " << total_skipped << " snarls" << std::endl;
        std::cerr << "Generated " << NUM_HAPLOTYPES << " haplotypes in " << total_generated << " components and took existing ones in " << total_took << " components" << std::endl;
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


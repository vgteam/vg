/** \file haplotypes_main.cpp
 *
 * Defines the "vg haplotypes" subcommand, which will ultimately sample haplotypes.
 *
 * This is currently highly experimental.
 */

// FIXME tests

#include "subcommand.hpp"

#include "../gbwt_helper.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../snarl_distance_index.hpp"

#include <gbwtgraph/algorithms.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <getopt.h>
#include <omp.h>

using namespace vg;

//----------------------------------------------------------------------------

void help_haplotypes(char** argv) {
    std::cerr << "usage: " << argv[0] << " " << argv[1] << " [options] graph.gbz" << std::endl;
    std::cerr << std::endl;
    // FIXME description
    std::cerr << "Some experiments with haplotype sampling." << std::endl;
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
constexpr size_t TARGET_DISTANCE = 10000;
constexpr size_t DISTANCE_BOUND = 1000000;

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

std::pair<handle_t, handle_t> empty_interval() {
    return { empty_handle(), empty_handle() };
}

std::string to_string(handle_t handle) {
    gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(handle);
    return std::string("(") + std::to_string(gbwt::Node::id(node)) + std::string(", ") + std::to_string(gbwt::Node::is_reverse(node)) + std::string(")");
}

std::vector<gbwt::size_type> get_paths(handle_t handle, const gbwt::FastLocate& r_index) {
    std::vector<gbwt::size_type> result = r_index.decompressSA(gbwtgraph::GBWTGraph::handle_to_node(handle));
    std::sort(result.begin(), result.end(), [&](size_t first, size_t second) -> bool {
        return (r_index.seqId(first) < r_index.seqId(second));
    });

    gbwt::size_type prev = std::numeric_limits<gbwt::size_type>::max();
    gbwt::size_type last_warning = std::numeric_limits<gbwt::size_type>::max();
    for (size_t i = 0; i < result.size(); i++) {
        gbwt::size_type seq_id = r_index.seqId(result[i]);
        if (seq_id == prev && seq_id != last_warning) {
            #pragma omp critical
            {
                std::cerr << "warning: [vg haplotypes] GBWT path " << seq_id << " visits node " << to_string(handle) << " multiple times" << std::endl;
            }
            last_warning = seq_id;
        }
        prev = seq_id;
    }

    return result;
}

// FIXME handle haplotypes that visit the same node multiple times properly
size_t count_haplotypes(handle_t from, handle_t to, const gbwt::FastLocate& r_index) {
    auto from_paths = get_paths(from, r_index);
    auto to_paths = get_paths(to, r_index);

    auto from_iter = from_paths.begin();
    auto to_iter = to_paths.begin();
    size_t result = 0;
    while (from_iter != from_paths.end() && to_iter != to_paths.end()) {
        gbwt::size_type from_id = r_index.seqId(*from_iter);
        gbwt::size_type to_id = r_index.seqId(*to_iter);
        if (from_id == to_id) {
            if (r_index.seqOffset(*from_iter) >= r_index.seqOffset(*to_iter)) {
                result++;
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

bool is_unary_path(handle_t from, handle_t to, const HandleGraph& graph) {
    handle_t curr = from;
    std::unordered_set<handle_t> visited;
    while (curr != to) {
        if (visited.find(curr) != visited.end()) {
            return false;
        }
        visited.insert(curr);
        handle_t successor = empty_handle();
        size_t successors = 0;
        graph.follow_edges(curr, false, [&](const handle_t& next) {
            successor = next;
            successors++;
        });
        if (successors != 1) {
            return false;
        }
        curr = successor;
    }
    return true;
}

//----------------------------------------------------------------------------

int main_haplotypes(int argc, char** argv) {
    double start = gbwt::readTimer();
    if (argc < 3) {
        help_haplotypes(argv);
        return 1;
    }

    // Parse options into these.
    std::string graph_name, r_index_name, distance_name;
    bool progress = false;

    // Process the arguments.
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "distance-index", required_argument, 0, 'd' },
            { "r_index", required_argument, 0, 'r' },
            { "progress", no_argument, 0, 'p' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "d:r:ph", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
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

    // Determine GBWT construction jobs.
    gbwtgraph::GBZ gbz;
    load_gbz(gbz, graph_name, progress);
    double checkpoint = gbwt::readTimer();
    if (progress) {
        std::cerr << "Determining GBWT construction jobs" << std::endl;
    }
    size_t size_bound = gbz.graph.get_node_count() / APPROXIMATE_JOBS;
    gbwtgraph::ConstructionJobs jobs = gbwtgraph::gbwt_construction_jobs(gbz.graph, size_bound);
    if (progress) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Partitioned " << jobs.components << " components into " << jobs.size() << " jobs in " << seconds << " seconds" << std::endl;
    }

    // Load the distance index and partition top-level chains between jobs.
    SnarlDistanceIndex distance_index;
    if (progress) {
        std::cerr << "Loading distance index from " << distance_name << std::endl;
    }
    distance_index.deserialize(distance_name);
    std::vector<std::vector<net_handle_t>> partitioned_chains = gbwtgraph::partition_chains(distance_index, gbz.graph, jobs);
    jobs.clear(); // This saves memory.

    // Load the r-index.
    gbwt::FastLocate r_index;
    load_r_index(r_index, r_index_name, progress);
    r_index.setGBWT(gbz.index);

    // FIXME what to do with chains that do not contain selected snarls?
    // FIXME what to do with intervals not covered by any haplotype
    // FIXME what to do with chrX, chrY, and other special cases
    // FIXME are all chains in the right orientation?
    if (progress) {
        std::cerr << "Processing chains using " << omp_get_max_threads() << " threads" << std::endl;
    }
    checkpoint = gbwt::readTimer();
    size_t total_snarls = 0, total_skipped = 0, total_intervals = 0, total_haplotypes = 0, total_no_coverage = 0, total_not_connected = 0;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < partitioned_chains.size(); job++) {
        const std::vector<net_handle_t>& chains = partitioned_chains[job];
        for (size_t chain = 0; chain < chains.size(); chain++) {
            std::vector<std::pair<handle_t, handle_t>> snarls;
            handle_t snarl_start = empty_handle();
            bool has_start = false;
            bool was_snarl = false;
            size_t step = 0;
            size_t skipped_snarls = 0;

            net_handle_t curr = distance_index.get_bound(chains[chain], false, true);
            net_handle_t chain_end = distance_index.get_bound(chains[chain], true, false);
            while (curr != chain_end) {
                if (distance_index.is_node(curr)) {
                    handle_t handle = distance_index.get_handle(curr, &gbz.graph);
                    if (was_snarl && has_start) {
                        size_t distance = get_distance(snarl_start, handle, gbz.graph, distance_index);
                        if (distance <= DISTANCE_BOUND) {
                            snarls.push_back({ snarl_start, handle });
                        } else {
                            // Use empty intervals to mark possible haplotype breakpoints.                        
                            snarls.push_back(empty_interval());
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
                        std::cerr << "Job " << job << ", chain " << chain << ": found " << successors << " successors at step " << step << std::endl;
                    }
                    break;
                }
                curr = next;
                step++;
            }

            // Combine snarls into intervals.
            std::vector<std::pair<handle_t, handle_t>> intervals;
            size_t total_distance = 0;
            size_t max_distance = 0;
            size_t head = 0;
            while (head < snarls.size()) {
                if (snarls[head] == empty_interval()) {
                    intervals.push_back(empty_interval());
                    head++;
                    continue;
                }
                size_t distance = get_distance(snarls[head].first, snarls[head].second, gbz.graph, distance_index);
                size_t tail = head;
                while (tail + 1 < snarls.size()) {
                    if (snarls[tail + 1] == empty_interval()) {
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
                intervals.push_back({ snarls[head].first, snarls[tail].second });
                head = tail + 1;
            }

            // Count the number of haplotypes in each interval.
            size_t haplotypes = 0;
            std::vector<std::pair<handle_t, handle_t>> no_haplotypes;
            for (auto interval : intervals) {
                if (interval == empty_interval()) {
                    continue;
                }
                size_t intersection = count_haplotypes(interval.first, interval.second, r_index);
                haplotypes += intersection;
                if (intersection == 0) {
                    no_haplotypes.push_back(interval);
                }
            }

            // Try to connect the intervals.
            size_t not_connected = 0;
            for (size_t i = 1; i < intervals.size(); i++) {
                if (intervals[i - 1] == empty_interval() || intervals[i] == empty_interval()) {
                    continue;
                }
                if (!is_unary_path(intervals[i - 1].second, intervals[i].first, gbz.graph)) {
                    not_connected++;
                }
            }

            #pragma omp critical
            {
                size_t real_snarls = snarls.size() - skipped_snarls;
                size_t real_intervals = intervals.size() - skipped_snarls;
                if (progress && real_intervals > 0) {
                    std::cerr << "Job " << job << ", chain " << chain << ":" << std::endl;
                    std::cerr << "  Selected " << real_snarls << " snarls and skipped " << skipped_snarls << " snarls" << std::endl;
                    std::cerr << "  Found " << real_intervals << " intervals with " << (static_cast<double>(haplotypes) / real_intervals) << " haplotypes each" << std::endl;
                    for (auto interval : no_haplotypes) {
                        std::cerr << "  No haplotypes from " << to_string(interval.first)
                                << " to " << to_string(interval.second)
                                << " (distance " << get_distance(interval.first, interval.second, gbz.graph, distance_index) << ")" << std::endl;
                    }
                    if (not_connected > 0) {
                        std::cerr << "  No connetions between " << not_connected << " pairs of intervals" << std::endl;
                    }
                    std::cerr << "  Average distance " << (static_cast<double>(total_distance) / real_intervals) << ", max distance " << max_distance << std::endl;
                }
                total_snarls += real_snarls;
                total_skipped += skipped_snarls;
                total_intervals += real_intervals;
                total_haplotypes += haplotypes;
                total_no_coverage += no_haplotypes.size();
                total_not_connected += not_connected;
            }
        }
    }

    if (progress) {
        double end = gbwt::readTimer();
        std::cerr << "Processed the chains in " << (end - checkpoint) << " seconds" << std::endl;
        std::cerr << "Combined " << total_snarls << " snarls into " << total_intervals << " intervals and skipped " << total_skipped << " snarls" << std::endl;
        std::cerr << "The average interval contains " << (static_cast<double>(total_haplotypes) / total_intervals) << " haplotypes" << std::endl;
        std::cerr << "There are " << total_no_coverage << " intervals with no haplotypes" << std::endl;
        std::cerr << "No connections between " << total_not_connected << " pairs of intervals" << std::endl;
        double gib = gbwt::inGigabytes(gbwt::memoryUsage());
        std::cerr << "Used " << (end - start) << " seconds, " << gib << " GiB" << std::endl;
    }

    return 0;
}

// FIXME description
static vg::subcommand::Subcommand vg_haplotypes("haplotypes", "haplotype sampling experiments", vg::subcommand::DEVELOPMENT, main_haplotypes);

//----------------------------------------------------------------------------


#include "coverage_depth.hpp"
#include <bdsg/hash_graph.hpp>
#include "algorithms/subgraph.hpp"
#include <vg/io/stream.hpp>
#include "../path.hpp"

namespace vg {
namespace algorithms {

/// Estimate the depth of coverage of a given (sub) graph using the packer
/// Coverage is computed relative to the given path
double packed_depth(const PathHandleGraph& graph, const Packer& packer, const string& ref_path) {

    // get the path length
    path_handle_t path_handle = graph.get_path_handle(ref_path);
    size_t path_len = 0;
    for (handle_t handle : graph.scan_path(path_handle)) {
        path_len += graph.get_length(handle);
    }
    if (path_len == 0) {
        return 0;
    }

    // sum up the coverage
    size_t tot_base_coverage = 0;
    graph.for_each_handle([&] (handle_t handle) {
            Position pos;
            pos.set_node_id(graph.get_id(handle));
            size_t packer_pos = packer.position_in_basis(pos);
            size_t node_len = graph.get_length(handle);
            for (size_t offset = 0; offset < node_len; ++offset) {
                tot_base_coverage += packer.coverage_at_position(packer_pos + offset);
            }
        });

    // return average (over the path)
    return (double)tot_base_coverage / (double)path_len;
}


/// Estimate the binned coverage along a path
map<size_t, double> binned_packed_depth(const PathHandleGraph& graph, const Packer& packer, const string& ref_path,
                                        size_t step, size_t context, size_t threads) {

    // move forward along path (note: this can be sped up if we're given a PathPositionHandleGraph but I don't think
    // it matters for a couple of scans. 
    function<size_t(step_handle_t&, size_t)> advance = [&] (step_handle_t& step_handle, size_t distance) {
        size_t went = 0;
        for (; graph.has_next_step(step_handle) && went < distance; step_handle = graph.get_next_step(step_handle)) {
            went += graph.get_length(graph.get_handle_of_step(step_handle));
        }
        return went;
    };

    path_handle_t path_handle = graph.get_path_handle(ref_path);
    step_handle_t step_handle = graph.path_begin(path_handle);

    // hop along the graph, grabbing a step handle every "step" bases (or thereabouts)
    vector<pair<size_t, step_handle_t>> bin_centers;
    size_t pos = advance(step_handle, step / 2);
    if (pos >= step / 2) {
        size_t went;
        do {
            if (bin_centers.empty() || step_handle != bin_centers.back().second) {
                bin_centers.push_back(make_pair(pos, step_handle));
            }
            went = advance(step_handle, step);
            pos += went;
        } while (went >= step);
    }

    // our graph's too small to do any stepping, just use the first handle
    if (bin_centers.empty()) {
        bin_centers.push_back(make_pair(0, graph.path_begin(path_handle)));
    }

    // visit every bin center and make a subgraph to collect coverage from
    if (threads == 0) {
        threads = get_thread_count();
    }
    map<size_t, double> binned_depths;
#pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < bin_centers.size(); ++i) {
        // extract the subgraph
        bdsg::HashGraph subgraph;
        step_handle_t bin_step = bin_centers[i].second;
        handle_t bin_handle = graph.get_handle_of_step(bin_step);
        assert(graph.get_is_reverse(bin_handle) == false);
        subgraph.create_handle(graph.get_sequence(bin_handle), graph.get_id(bin_handle));
        expand_subgraph_by_steps(graph, subgraph, context);

        // sum up the coverage on the subgraph
        size_t tot_base_coverage = 0;
        size_t tot_ref_len = 0;
        subgraph.for_each_handle([&] (handle_t sub_handle) {
                // go back into the original graph because we don't have any
                // path information in the subgraph because we are unable
                // to get it without requiring the path position interface
                handle_t orig_handle = graph.get_handle(subgraph.get_id(sub_handle));
                Position pos;
                pos.set_node_id(graph.get_id(orig_handle));
                size_t packer_pos = packer.position_in_basis(pos);
                size_t node_len = graph.get_length(orig_handle);
                for (size_t offset = 0; offset < node_len; ++offset) {
                    tot_base_coverage += packer.coverage_at_position(packer_pos + offset);
                }
                // we manually test if each handle is on our reference path (again, to
                // not require path position interface)
                vector<step_handle_t> step_path_handles = graph.steps_of_handle(orig_handle);
                bool on_ref = false;
                for (size_t j = 0; j < step_path_handles.size() && !on_ref; ++j) {
                    on_ref = graph.get_path_handle_of_step(step_path_handles[j]) == path_handle;
                }
                if (on_ref) {
                    tot_ref_len += node_len;
                }
            });
        
        assert(tot_ref_len > 0);
        double avg_base_coverage = tot_base_coverage / tot_ref_len;
        
#pragma omp critical (update_binned_depth)
        binned_depths[bin_centers[i].first] = avg_base_coverage;
    }

    return binned_depths;

}

// draw (roughly) max_nodes nodes from the graph using the random seed
static unordered_map<nid_t, size_t> sample_nodes(const HandleGraph& graph, size_t max_nodes, size_t random_seed) {
    default_random_engine generator(random_seed);
    uniform_real_distribution<double> distribution(0, 1);
    double cutoff = std::min((double)1.0, (double)(max_nodes / graph.get_node_count()));
    unordered_map<nid_t, size_t> sampled_nodes;
    graph.for_each_handle([&](handle_t handle) {
        if (cutoff == 1 || cutoff < distribution(generator)) {
            sampled_nodes[graph.get_id(handle)] = 0;
        }
      });
    return sampled_nodes;
}

// update the coverage from an alignment.  only count nodes that are in the map already
static void update_sample_gam_depth(const Alignment& aln, unordered_map<nid_t, size_t>& node_coverage) {
    const Path& path = aln.path();
    for (int i = 0; i < path.mapping_size(); ++i) {
        const Mapping& mapping = path.mapping(i);
        nid_t node_id = mapping.position().node_id();
        if (node_coverage.count(node_id)) {
            ++node_coverage[node_id];
        } 
    }
}

// sum up the results from the different threads and return the average.
// if a min_coverage is given, nodes with less coverage are ignored
static double combine_and_average_node_coverages(vector<unordered_map<nid_t, size_t>>& node_coverages, size_t min_coverage) {
    for (int i = 1; i < node_coverages.size(); ++i) {
        for (const auto& node_cov : node_coverages[i]) {
            node_coverages[0][node_cov.first] += node_cov.second;
        }
    }
    size_t tot_coverage = 0;
    size_t tot_count = 0;
    for (const auto & node_cov : node_coverages[0]) {
        if (node_cov.second >= min_coverage) {
            tot_coverage += node_cov.second;
            ++tot_count;
        }
    }

    return tot_count > 0 ? (double)tot_coverage / (double)tot_count : 0;
}


double sample_gam_depth(const HandleGraph& graph, istream& gam_stream, size_t max_nodes, size_t random_seed, size_t min_coverage) {
    // one node counter per thread
    vector<unordered_map<nid_t, size_t>> node_coverages(get_thread_count(), sample_nodes(graph, max_nodes, random_seed));

    function<void(Alignment& aln)> aln_callback = [&](Alignment& aln) {
        update_sample_gam_depth(aln, node_coverages[omp_get_thread_num()]);
    };
    vg::io::for_each_parallel(gam_stream, aln_callback);
    return combine_and_average_node_coverages(node_coverages, min_coverage);
}

double sample_gam_depth(const HandleGraph& graph, const vector<Alignment>& alignments, size_t max_nodes, size_t random_seed, size_t min_coverage) {
    // one node counter per thread
    vector<unordered_map<nid_t, size_t>> node_coverages(get_thread_count(), sample_nodes(graph, max_nodes, random_seed));

#pragma omp parallel for
    for (size_t i = 0; i < alignments.size(); ++i) {
        update_sample_gam_depth(alignments[i], node_coverages[omp_get_thread_num()]);
    }
    return combine_and_average_node_coverages(node_coverages, min_coverage);
}



}




}


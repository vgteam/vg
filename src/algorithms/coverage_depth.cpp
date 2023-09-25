#include "coverage_depth.hpp"
#include <bdsg/hash_graph.hpp>
#include "algorithms/subgraph.hpp"
#include <vg/io/stream.hpp>
#include "../path.hpp"

namespace vg {
namespace algorithms {

void packed_depths(const Packer& packer, const string& path_name, size_t min_coverage, ostream& out_stream) {
    const PathHandleGraph& graph = dynamic_cast<const PathHandleGraph&>(*packer.get_graph());
    path_handle_t path_handle = graph.get_path_handle(path_name);
    step_handle_t start_step = graph.path_begin(path_handle);
    step_handle_t end_step = graph.path_end(path_handle);
    Position cur_pos;
    subrange_t subrange;
    string base_name = Paths::strip_subrange(path_name, &subrange);
    size_t path_offset = subrange == PathMetadata::NO_SUBRANGE ? 1 : 1 + subrange.first;
    
    for (step_handle_t cur_step = start_step; cur_step != end_step; cur_step = graph.get_next_step(cur_step)) {
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        nid_t cur_id = graph.get_id(cur_handle);
        size_t cur_len = graph.get_length(cur_handle);
        cur_pos.set_node_id(cur_id);
        cur_pos.set_is_reverse(graph.get_is_reverse(cur_handle));
        for (size_t i = 0; i < cur_len; ++i) {
            cur_pos.set_offset(i);
            size_t pos_coverage = packer.coverage_at_position(packer.position_in_basis(cur_pos));
            if (pos_coverage >= min_coverage) {
                out_stream << base_name << "\t" << path_offset << "\t" << pos_coverage << "\n";
            }
            ++path_offset;
        }
    }
}

pair<double, double> packed_depth_of_bin(const Packer& packer,
                                         step_handle_t start_step, step_handle_t end_plus_one_step,
                                         size_t min_coverage, bool include_deletions) {

    const PathHandleGraph& graph = dynamic_cast<const PathHandleGraph&>(*packer.get_graph());

    // coverage of each node via deletion (that's contained in the bin)
    unordered_map<nid_t, size_t> deletion_coverages;
    if (include_deletions) {
        const VectorizableHandleGraph* vec_graph = dynamic_cast<const VectorizableHandleGraph*>(packer.get_graph());
        unordered_map<handle_t, step_handle_t> deletion_candidates;
        handle_t prev_handle;
        for (step_handle_t cur_step = start_step; cur_step != end_plus_one_step; cur_step = graph.get_next_step(cur_step)) {
            handle_t cur_handle = graph.get_handle_of_step(cur_step);
            graph.follow_edges(cur_handle, true, [&] (handle_t other) {
                    if (!deletion_candidates.empty() && other!= prev_handle && deletion_candidates.count(other)) {
                        edge_t edge = graph.edge_handle(other, cur_handle);
                        size_t edge_pos = vec_graph->edge_index(edge);
                        size_t deletion_coverage = packer.edge_coverage(edge_pos);
                        // quadratic alert.  if this is too slow, can use interval tree or something
                        for (step_handle_t del_step = graph.get_next_step(deletion_candidates[other]);
                             del_step != cur_step;
                             del_step = graph.get_next_step(del_step)) {
                            handle_t del_handle = graph.get_handle_of_step(del_step);
                            nid_t del_id = graph.get_id(del_handle);
                            if (!deletion_coverages.count(del_id)) {
                                deletion_coverages[del_id] = deletion_coverage;
                            } else {
                                deletion_coverages[del_id] += deletion_coverage;
                            }
                        }
                    }
                });
            prev_handle = cur_handle;
            deletion_candidates[cur_handle] = cur_step;
        }
    }

    // compute the mean and variance of our base coverage across the bin
    size_t bin_length = 0;
    double mean = 0.0;
    double M2 = 0.0;

    for (step_handle_t cur_step = start_step; cur_step != end_plus_one_step; cur_step = graph.get_next_step(cur_step)) {
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        nid_t cur_id = graph.get_id(cur_handle);
        size_t cur_len = graph.get_length(cur_handle);
        size_t del_coverage = !include_deletions or !deletion_coverages.count(cur_id) ? 0 : deletion_coverages[cur_id];
        Position cur_pos;
        cur_pos.set_node_id(cur_id);
        cur_pos.set_is_reverse(graph.get_is_reverse(cur_handle));
        for (size_t i = 0; i < cur_len; ++i) {
            cur_pos.set_offset(i);
            size_t pos_coverage = packer.coverage_at_position(packer.position_in_basis(cur_pos)) + del_coverage;
            if (pos_coverage >= min_coverage) {
                wellford_update(bin_length, mean, M2, pos_coverage);
            }
        }
    }
    return wellford_mean_var(bin_length, mean, M2);
}

vector<tuple<size_t, size_t, double, double>> binned_packed_depth(const Packer& packer, const string& path_name, size_t bin_size,
                                                                  size_t min_coverage, bool include_deletions) {

    const PathHandleGraph& graph = dynamic_cast<const PathHandleGraph&>(*packer.get_graph());
    path_handle_t path_handle = graph.get_path_handle(path_name);
    
    // one scan of our path to collect the bins
    step_handle_t start_step = graph.path_begin(path_handle);
    step_handle_t end_step = graph.path_end(path_handle);
    vector<pair<size_t, step_handle_t>> bins; // start offset / start step of each bin
    size_t offset = 0;
    size_t cur_bin_size = bin_size;
    for (step_handle_t cur_step = start_step; cur_step != end_step; cur_step = graph.get_next_step(cur_step)) {
        if (cur_bin_size >= bin_size) {
            bins.push_back(make_pair(offset, cur_step));
            cur_bin_size = 0;
        }
        size_t node_len = graph.get_length(graph.get_handle_of_step(cur_step));
        offset += node_len;
        cur_bin_size += node_len;
    }

    // parallel scan to compute the coverages
    vector<tuple<size_t, size_t, double, double>> binned_depths(bins.size());
#pragma omp parallel for
    for (size_t i = 0; i < bins.size(); ++i) {
        step_handle_t bin_start_step = bins[i].second;
        step_handle_t bin_end_step = i < bins.size() - 1 ? bins[i+1].second : end_step;
        size_t bin_start = bins[i].first;
        size_t bin_end = i < bins.size() - 1 ? bins[i+1].first : offset;
        pair<double, double> coverage = packed_depth_of_bin(packer, bin_start_step, bin_end_step, min_coverage, include_deletions);
        binned_depths[i] = make_tuple(bin_start, bin_end, coverage.first, coverage.second);
    }

    return binned_depths;
}

BinnedDepthIndex binned_packed_depth_index(const Packer& packer,
                                           const vector<string>& path_names,
                                           size_t min_bin_size,
                                           size_t max_bin_size,
                                           double exp_growth_factor,
                                           size_t min_coverage,
                                           bool include_deletions,
                                           bool std_err) {
    const PathHandleGraph& graph = dynamic_cast<const PathHandleGraph&>(*packer.get_graph());
    
    BinnedDepthIndex depth_index;
    for (const string& path_name : path_names) {
        size_t path_max_bin = 0;
        graph.for_each_step_in_path(graph.get_path_handle(path_name), [&] (step_handle_t step_handle) {
                path_max_bin += graph.get_length(graph.get_handle_of_step(step_handle));
                return path_max_bin < max_bin_size;
            });
        path_max_bin = std::min(max_bin_size, path_max_bin);

        map<size_t, map<size_t, pair<float, float>>>& scaled_depth_map = depth_index[path_name];
        size_t prev_bin_size = 0;
        for (size_t bin_size = min_bin_size; bin_size != prev_bin_size;) {

            map<size_t, pair<float, float>>& depth_map = scaled_depth_map[bin_size];
            vector<tuple<size_t, size_t, double, double>> binned_depths = binned_packed_depth(packer, path_name, bin_size,
                                                                                              min_coverage, include_deletions);
            // todo: probably more efficent to just leave in sorted vector
            for (auto& binned_depth : binned_depths) {
                double var = get<3>(binned_depth);
                // optionally convert variance to standard error
                if (std_err) {
                    var = sqrt(var / (double)(get<1>(binned_depth) - get<0>(binned_depth)));
                }
                depth_map[get<0>(binned_depth)] = make_pair(get<2>(binned_depth), var);
            }

            prev_bin_size = bin_size;
            // todo: trim out useless last bins that are only a bit bigger than prev
            bin_size = std::min(path_max_bin, (size_t)pow(bin_size, exp_growth_factor));
        }
    }
    return depth_index;
}


pair<float, float> get_depth_from_index(const BinnedDepthIndex& depth_index, const string& path_name, size_t start_offset, size_t end_offset) {
     
    // accept backward ranges
    if (end_offset < start_offset) {
        swap(start_offset, end_offset);
    }
    size_t bin_size = 1 + end_offset - start_offset;
    // pad it out
    bin_size *= 2;

    auto ub1 = depth_index.at(path_name).upper_bound(bin_size);
    if (ub1 == depth_index.at(path_name).end()) {
        --ub1;
    }
    auto ub = ub1->second.upper_bound(start_offset);
    --ub;
    auto ub_end = ub1->second.upper_bound(end_offset);
    size_t count = 0;
    pair<float, float> total = make_pair(0, 0);
    for (auto cur = ub; cur != ub_end; ++cur, ++count) {
        total.first += cur->second.first;
        total.second += cur->second.second;
    }
    // todo: better way of combining?
    total.first /= (double)count;
    total.second /= (double)count;
    return total;
}

// draw (roughly) max_nodes nodes from the graph using the random seed
static unordered_map<nid_t, size_t> sample_nodes(const HandleGraph& graph, size_t max_nodes, size_t random_seed) {
    default_random_engine generator(random_seed);
    uniform_real_distribution<double> distribution(0, 1);
    double cutoff = std::min((double)1.0, (double)max_nodes / (double)graph.get_node_count());
    unordered_map<nid_t, size_t> sampled_nodes;
    graph.for_each_handle([&](handle_t handle) {
        if (cutoff == 1. || cutoff <= distribution(generator)) {
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
            // we add the number of bases covered
            node_coverage[node_id] += mapping_from_length(mapping);
        } 
    }
}

// sum up the results from the different threads and return the average.
// if a min_coverage is given, nodes with less coverage are ignored
static pair<double, double> combine_and_average_node_coverages(const HandleGraph& graph, vector<unordered_map<nid_t, size_t>>& node_coverages, size_t min_coverage) {
    for (int i = 1; i < node_coverages.size(); ++i) {
        for (const auto& node_cov : node_coverages[i]) {
            node_coverages[0][node_cov.first] += node_cov.second;
        }
    }
    size_t count = 0;
    double mean = 0.;
    double M2 = 0.;
    for (const auto & node_cov : node_coverages[0]) {
        if (node_cov.second >= min_coverage) {
            // we normalize the bases covered by the node length as we sum
            double node_len = graph.get_length(graph.get_handle(node_cov.first));
            wellford_update(count, mean, M2, (double)node_cov.second / node_len);
        }
    }

    return wellford_mean_var(count, mean, M2);
}


pair<double, double> sample_mapping_depth(const HandleGraph& graph, const string& input_filename, size_t max_nodes, size_t random_seed, size_t min_coverage, size_t min_mapq, const string& format) {
    // one node counter per thread
    vector<unordered_map<nid_t, size_t>> node_coverages(get_thread_count(), sample_nodes(graph, max_nodes, random_seed));

    function<void(Alignment& aln)> aln_callback = [&](Alignment& aln) {
        if (aln.mapping_quality() >= min_mapq) {
            update_sample_gam_depth(aln, node_coverages[omp_get_thread_num()]);
        }
    };
    if (format == "GAM") {
        get_input_file(input_filename, [&] (istream& gam_stream) {
                vg::io::for_each_parallel(gam_stream, aln_callback);
            });
    } else if (format == "GAF") {
        vg::io::gaf_unpaired_for_each_parallel(graph, input_filename, aln_callback);
    } else {
        throw runtime_error("vg::aglorithms::coverage_depth: Invalid format specified for sample_mapping_depth(): " +
                            format + ". Valid options are GAM and GAF.");
    }
    
    return combine_and_average_node_coverages(graph, node_coverages, min_coverage);
}



pair<double, double> sample_gam_depth(const HandleGraph& graph, const vector<Alignment>& alignments, size_t max_nodes, size_t random_seed, size_t min_coverage, size_t min_mapq) {
    // one node counter per thread
    vector<unordered_map<nid_t, size_t>> node_coverages(get_thread_count(), sample_nodes(graph, max_nodes, random_seed));

#pragma omp parallel for
    for (size_t i = 0; i < alignments.size(); ++i) {
        if (alignments[i].mapping_quality() >= min_mapq) {
            update_sample_gam_depth(alignments[i], node_coverages[omp_get_thread_num()]);
        }
    }
    return combine_and_average_node_coverages(graph, node_coverages, min_coverage);
}

void path_depths(const PathHandleGraph& graph, const string& path_name, size_t min_coverage, bool count_cycles,
                 const unordered_set<path_handle_t>& ignore_paths, ostream& out_stream) {
    assert(graph.has_path(path_name));

    path_handle_t path_handle = graph.get_path_handle(path_name);
    // big speedup
    unordered_map<path_handle_t, string> path_to_name;

    subrange_t subrange;
    string base_name = Paths::strip_subrange(path_name, &subrange);
    size_t offset = subrange == PathMetadata::NO_SUBRANGE ? 1 : 1 + subrange.first;

    graph.for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
            unordered_set<string> path_set;
            size_t step_count = 0;            
            handle_t handle = graph.get_handle_of_step(step_handle);
            graph.for_each_step_on_handle(handle, [&](step_handle_t step_handle_2) {
                    path_handle_t step_path_handle = graph.get_path_handle_of_step(step_handle_2);
                    if (ignore_paths.count(step_path_handle)) {
                        return true;
                    }
                    if (count_cycles) {
                        ++step_count;
                    } else {

                        auto it = path_to_name.find(step_path_handle);
                        if (it == path_to_name.end()) {
                            string step_path_name = graph.get_path_name(step_path_handle);
                            // disregard subpath tags when counting
                            it = path_to_name.insert(make_pair(step_path_handle, Paths::strip_subrange(step_path_name))).first;
                        }
                        path_set.insert(it->second);
                    }
                    return true;
                });
            size_t coverage = (count_cycles ? step_count : path_set.size()) - 1;
            size_t node_len = graph.get_length(handle);
            if (coverage >= min_coverage) {
                for (size_t i = 0; i < node_len; ++i) {
                    if (coverage >= min_coverage) {
                        out_stream << base_name << "\t" << (offset + i) << "\t" << coverage << "\n";
                    }
                }
            }
            offset += node_len;            
        });
}

pair<double, double> path_depth_of_bin(const PathHandleGraph& graph,
                                       step_handle_t start_step, step_handle_t end_plus_one_step,
                                       size_t min_coverage, bool count_cycles, const unordered_set<path_handle_t>& ignore_paths) {

    // compute the mean and variance of our base coverage across the bin
    size_t bin_length = 0;
    double mean = 0.0;
    double M2 = 0.0;

    // big speedup
    unordered_map<path_handle_t, string> path_to_name;

    for (step_handle_t cur_step = start_step; cur_step != end_plus_one_step; cur_step = graph.get_next_step(cur_step)) {
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        nid_t cur_id = graph.get_id(cur_handle);
        size_t cur_len = graph.get_length(cur_handle);

        unordered_set<string> path_set;
        size_t step_count = 0;
        graph.for_each_step_on_handle(cur_handle, [&](step_handle_t step_handle) {
                path_handle_t step_path_handle = graph.get_path_handle_of_step(step_handle);
                if (ignore_paths.count(step_path_handle)) {
                    return true;
                }
                if (count_cycles) {
                    ++step_count;
                } else {
                        auto it = path_to_name.find(step_path_handle);
                        if (it == path_to_name.end()) {
                            string step_path_name = graph.get_path_name(step_path_handle);
                            // disregard subpath tags when counting
                            it = path_to_name.insert(make_pair(step_path_handle, Paths::strip_subrange(step_path_name))).first;
                        }
                        path_set.insert(it->second);
                }
                return true;
            });
        size_t coverage = (count_cycles ? step_count : path_set.size()) - 1;        

        if (coverage >= min_coverage) {
            // todo: iteration here copied from packer implementation, not necessary
            for (size_t i = 0; i < cur_len; ++i) {
                wellford_update(bin_length, mean, M2, coverage);
            }
        }
    }
    return wellford_mean_var(bin_length, mean, M2);
}

vector<tuple<size_t, size_t, double, double>> binned_path_depth(const PathHandleGraph& graph,
                                                                const string& path_name,
                                                                size_t bin_size,
                                                                size_t min_coverage,
                                                                bool count_cycles,
                                                                const unordered_set<path_handle_t>& ignore_paths) {

    path_handle_t path_handle = graph.get_path_handle(path_name);
    
    // one scan of our path to collect the bins
    step_handle_t start_step = graph.path_begin(path_handle);
    step_handle_t end_step = graph.path_end(path_handle);
    vector<pair<size_t, step_handle_t>> bins; // start offset / start step of each bin
    size_t offset = 0;
    size_t cur_bin_size = bin_size;
    for (step_handle_t cur_step = start_step; cur_step != end_step; cur_step = graph.get_next_step(cur_step)) {
        if (cur_bin_size >= bin_size) {
            bins.push_back(make_pair(offset, cur_step));
            cur_bin_size = 0;
        }
        size_t node_len = graph.get_length(graph.get_handle_of_step(cur_step));
        offset += node_len;
        cur_bin_size += node_len;
    }

    // parallel scan to compute the coverages
    vector<tuple<size_t, size_t, double, double>> binned_depths(bins.size());
#pragma omp parallel for
    for (size_t i = 0; i < bins.size(); ++i) {
        step_handle_t bin_start_step = bins[i].second;
        step_handle_t bin_end_step = i < bins.size() - 1 ? bins[i+1].second : end_step;
        size_t bin_start = bins[i].first;
        size_t bin_end = i < bins.size() - 1 ? bins[i+1].first : offset;
        pair<double, double> coverage = path_depth_of_bin(graph, bin_start_step, bin_end_step, min_coverage, count_cycles,
                                                          ignore_paths);
        binned_depths[i] = make_tuple(bin_start, bin_end, coverage.first, coverage.second);
    }

    return binned_depths;
}



}




}


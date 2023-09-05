#ifndef VG_DEPTH_HPP_INCLUDED
#define VG_DEPTH_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "handle.hpp"
#include "statistics.hpp"
#include "packer.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/// print path-name offset base-coverage for every base on a path (just like samtools depth)
/// ignoring things below min_coverage.  offsets are 1-based in output stream
void packed_depths(const Packer& packer, const string& path_name, size_t min_coverage, ostream& out_stream);

/// Estimate the coverage along a given reference path interval [start_step, end_plus_one_step)
/// Coverage is obtained only from positions along the path, and variation is not counted
/// Except if "include_deletions" is true, then reference path positions covered by a deletion edge
/// (which is contained in the bin) will get the deletion edge's coverage counted.
/// Other types of events (such as SNPs) can throw off coverage in similar ways but deletions tend to be bigger
/// (and easier to find), so we hope that counting them is enough.
/// If one wants to infer deletions from the coverage, obviously this should be false, but if looking for
/// a background coverage for genotyping, then setting it to true may be helpful
pair<double, double> packed_depth_of_bin(const Packer& packer, step_handle_t start_step, step_handle_t end_plus_one_step,
                                         size_t min_coverage, bool include_deletions);

/// Use all available threads to estimate the binned packed coverage of a path using above fucntion
/// Each element is a bin's 0-based open-ended interval in the path, and its coverage mean,variance. 
vector<tuple<size_t, size_t, double, double>> binned_packed_depth(const Packer& packer, const string& path_name, size_t bin_size,
                                                                  size_t min_coverage, bool include_deletions);

/// Use the above function to retrieve the binned depths of a list of paths, and store them indexed by start
/// coordinate.  If std_err is true, store <mean, stderr> instead of <mean, variance>
/// For each path, a series of indexes is computed, for bin sizes from min_bin_size, min_bin_size^(exp_growth_factor), etc.
using BinnedDepthIndex = unordered_map<string, map<size_t, map<size_t, pair<float, float>>>>;
BinnedDepthIndex binned_packed_depth_index(const Packer& packer,
                                           const vector<string>& path_names,
                                           size_t min_bin_size,
                                           size_t max_bin_size,
                                           double exp_growth_factor,
                                           size_t min_coverage,
                                           bool include_deletions,
                                           bool std_err);

/// Query index created above
pair<float, float> get_depth_from_index(const BinnedDepthIndex& depth_index, const string& path_name, size_t start_offset, size_t end_offset);

/// Return the mean and variance of coverage of randomly sampled nodes from a mappings file
/// Nodes with less than min_coverage are ignored
/// The input_filename can be - for stdin
/// The stream is scanned in parallel with all threads
/// max_nodes is used to keep memory down
/// valid formats are "GAM" and "GAF"
pair<double, double> sample_mapping_depth(const HandleGraph& graph, const string& input_filename, size_t max_nodes, size_t random_seed, size_t min_coverage, size_t min_mapq, const string& format="GAM");

/// As above, but read a vector instead of a stream
pair<double, double> sample_mapping_depth(const HandleGraph& graph, const vector<Alignment>& alignments, size_t max_nodes, size_t random_seed, size_t min_coverage, size_t min_mapq);

/// print path-name offset base-coverage for every base on a path (just like samtools depth)
/// ignoring things below min_coverage.  offsets are 1-based in output stream
/// coverage here is the number of steps from (unique) other paths
void path_depths(const PathHandleGraph& graph, const string& path_name, size_t min_coverage, bool count_cycles, ostream& out_stream);

/// like packed_depth_of_bin (above), but use paths (as in path_depths) for measuring coverage
pair<double, double> path_depth_of_bin(const PathHandleGraph& graph, step_handle_t start_step, step_handle_t end_plus_one_step,
                                       size_t min_coverage, bool count_cycles);

vector<tuple<size_t, size_t, double, double>> binned_path_depth(const PathHandleGraph& graph, const string& path_name, size_t bin_size,
                                                                size_t min_coverage, bool count_cycles);

}
}

#endif

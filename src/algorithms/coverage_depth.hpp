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

/// Return the mean and variance of coverage of randomly sampled nodes from a GAM
/// Nodes with less than min_coverage are ignored
/// The stream is scanned in parallel with all threads
/// max_nodes is used to keep memory down
pair<double, double> sample_gam_depth(const HandleGraph& graph, istream& gam_stream, size_t max_nodes, size_t random_seed, size_t min_coverage, size_t min_mapq);

/// As above, but read a vector instead of a stream
pair<double, double> sample_gam_depth(const HandleGraph& graph, const vector<Alignment>& alignments, size_t max_nodes, size_t random_seed, size_t min_coverage, size_t min_mapq);

}
}

#endif

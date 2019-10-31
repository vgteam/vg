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

/// Estimate the depth of coverage of a given (sub) graph using the packer
/// Coverage is computed relative to the given path
double packed_depth(const PathHandleGraph& graph, const Packer& packer, const string& ref_path);

/// Estimate the binned coverage along a path using the packer
/// ref_path is scanned, and every "step" bases as subgraph is extracted using the given number of context steps
/// If threads is 0, all the threads are used
map<size_t, double> binned_packed_depth(const PathHandleGraph& graph, const Packer& packer, const string& ref_path,
                                        size_t step, size_t context, size_t threads = 0);


/// Get the depth of a bin
/// the "k_nearest" closest bins to the given position are used
/// bins with coverage below min_coverage are ignored
double get_binned_depth(const unordered_map<size_t, double>& binned_depths, size_t pos, size_t k_nearest = 3, double min_coverage = 1.0);

/// Return the average depth of coverage of randomly sampled nodes from a GAM
/// Nodes with less than min_coverage are ignored
/// The stream is scanned in parallel with all threads
/// max_nodes is used to keep memory down
double sample_gam_depth(const HandleGraph& graph, istream& gam_stream, size_t max_nodes, size_t random_seed, size_t min_coverage = 1.0);

/// As above, but read a vector instead of a stream
double sample_gam_depth(const HandleGraph& graph, const vector<Alignment>& alignments, size_t max_nodes, size_t random_seed, size_t min_coverage = 1.0);

}
}

#endif

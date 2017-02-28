//
//  vg_algorithms.hpp
//  

// A collection of utility algorithms to manipulating graphs that can be use by both VG and XG objects
//

#ifndef vg_algorithms_hpp
#define vg_algorithms_hpp

#include <stdio.h>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <vector>

#include "position.hpp"
#include "xg.hpp"
#include "vg.hpp"
#include "vg.pb.h"
#include "hash_map.hpp"
#include "json2pb.h"

namespace vg {
namespace algorithms {
    
    /// Fills Graph g with the subgraph of the VG graph vg that connects two positions. The nodes that contain
    /// the two positions will be "cut" at the position and will be tips in the returned graph. Sometimes it
    /// is necessary to duplicate nodes in order to do this, so a map is returned that translates node IDs in
    /// g to node IDs in vg. By default, the algorithm provides one and only one guarantee:
    ///   - g contains all paths between pos_1 and pos_2 under the maximum length except paths that include
    ///     a cycle involving either position
    /// The algorithm optionally provides additional guarantees at the expense of increased computational cost,
    /// but no increase in asymptotic complexity (the guarantees are described below). If no path between the
    /// two positions under the maximum length exists, g will be left empty. An error is thrown if g is not
    /// empty when passed to function.
    ///
    /// Args:
    ///  vg                         graph to extract subgraph from
    ///  g                          graph to extract into
    ///  max_len                    guarantee finding paths along which pos_1 and pos_2 are this distance apart
    ///  pos_1                      start position, subgraph paths begin from here in same orientation
    ///  pos_2                      end position, subgraph paths end here in the same orientation
    ///  include_terminal_pos       cut nodes behind the positions (so that pos_1 and pos_2 are included in g)
    ///                             rather than in front (so that pos_1 and pos_2 are removed)
    ///  detect_terminal_cycles     also find paths that include cycles involving pos_1 and/or pos_2
    ///  no_additional_tips         make the node(s) containing pos_1 and pos_2 the only tips in g
    ///  only_paths                 only extract nodes and edges if they fall on some path between pos_1 and
    ///                             pos_2 (implies no_additional_tips = true)
    ///  strict_max_len             only extract nodes and edges if they fall on some path between pos_1 and
    ///                             pos_2 that is under the maximum length (implies only_paths = true)
    unordered_map<id_t, id_t> extract_connecting_graph(VG& vg, Graph& g, int64_t max_len,
                                                       pos_t pos_1, pos_t pos_2,
                                                       bool include_terminal_positions = false,
                                                       bool detect_terminal_cycles = false,
                                                       bool no_additional_tips = false,
                                                       bool only_paths = false,
                                                       bool strict_max_len = false);
    
    
    /// Same semantics as previous, but accesses graph through an XG instead of a VG. Optionally uses
    /// an LRUCache to speed up Node queries (recommended for detecting terminal cycles in particular).
    unordered_map<id_t, id_t> extract_connecting_graph(xg::XG& xg_index, Graph& g, int64_t max_len,
                                                       pos_t pos_1, pos_t pos_2,
                                                       bool include_terminal_positions = false,
                                                       bool detect_terminal_cycles = false,
                                                       bool no_additional_tips = false,
                                                       bool only_paths = false,
                                                       bool strict_max_len = false,
                                                       LRUCache<id_t, Node>* node_cache = nullptr);
    
}
}

#endif /* vg_algorithms_hpp */

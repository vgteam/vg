//
//  vg_algorithms.hpp
//  
<<<<<<< HEAD
// A collection of utility algorithms to manipulating graphs
=======
// A collection of utility algorithms to manipulating graphs that can be use by both VG and XG objects
>>>>>>> origin/band-right
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
<<<<<<< HEAD
=======
#include "json2pb.h"
>>>>>>> origin/band-right

namespace vg {
namespace algorithms {
    
<<<<<<< HEAD
    /// Fills Graph g with the section of the graph that connects two positions. The resulting
    /// graph will contain all paths that move forward from the first position to the second position,
    /// subject to a maximum length. Nodes and Edges that are not on some path between the two positions
    /// are excluded from the graph. The nodes containing the end positions will be cut just past the
    /// position (so that the two positions themselves are NOT included in the graph). To maintain all paths
    /// between the positions while cutting, it is sometimes necessary to duplicate nodes. Accordingly,
    /// the function returns a map to translate node IDs in the g to node IDs in the original graph. If
    /// no path between the two positions is under the maximum length, g is left empty. Throws an exception
    /// if g is not empty when function is called.
    unordered_map<id_t, id_t> extract_connecting_graph(VG& vg, Graph& g, int64_t max_len,
                                                       const pos_t& pos_1, const pos_t& pos_2);
    
    /// Same semantics as previous, but accesses graph through an XG object instead of VG.
    unordered_map<id_t, id_t> extract_connecting_graph(xg::XG& xg_index, Graph& g, int64_t max_len,
                                                       const pos_t& pos_1, const pos_t& pos_2);
=======
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
>>>>>>> origin/band-right
    
}
}

#endif /* vg_algorithms_hpp */

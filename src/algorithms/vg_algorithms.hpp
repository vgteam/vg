//
//  vg_algorithms.hpp
//  
// A collection of utility algorithms to manipulating graphs
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
    
    /// Fills Graph g with the section of the graph that connects two positions. The resulting
    /// graph will contain all paths that move forward from the first position to the second position,
    /// subject to a maximum length. Nodes and Edges that are not on some path between the two positions
    /// are excluded from the graph. Will not find paths that involve a cycle containing the two positions
    /// unless specificied (at a performance cost of up to a factor of 2). The nodes containing the end
    /// positions will be cut just past the position (so that the two positions themselves are NOT included
    /// in the graph). To maintain all paths between the positions while cutting, it is sometimes necessary
    /// to duplicate nodes. Accordingly, the function returns a map to translate node IDs in the g to node
    /// IDs in the original graph. If  no path between the two positions is under the maximum length, g is
    /// left empty. Throws an exception if g is not empty when function is called.
    unordered_map<id_t, id_t> extract_connecting_graph(VG& vg, Graph& g, int64_t max_len,
                                                       pos_t pos_1, pos_t pos_2,
                                                       bool detect_terminal_cycles = false);
    
    /// Same semantics as previous, but accesses graph through an XG instead of a VG. Optionally uses
    /// an LRUCache to speed up Node queries.
    unordered_map<id_t, id_t> extract_connecting_graph(xg::XG& xg_index, Graph& g, int64_t max_len,
                                                       pos_t pos_1, pos_t pos_2,
                                                       bool detect_terminal_cycles = false,
                                                       LRUCache<id_t, Node>* node_cache = nullptr);
    
}
}

#endif /* vg_algorithms_hpp */

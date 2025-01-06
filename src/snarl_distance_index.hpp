#ifndef VG_SNARL_DISTANCE_HPP_INCLUDED
#define VG_SNARL_DISTANCE_HPP_INCLUDED

#include <bdsg/snarl_distance_index.hpp>
#include "snarls.hpp"
#include <structures/union_find.hpp>
#include "hash_map.hpp"
#include <gbwtgraph/minimizer.h>


namespace vg { 

using namespace sdsl;
using namespace handlegraph;
using namespace bdsg;

//TODO: If anyone ever remakes the distance index, it would be really helpful for the multicomponent chains to know the lengths of each component

//Minimum distance taking a pos instead of id/orientation/offset
size_t minimum_distance(const SnarlDistanceIndex& distance_index, pos_t pos1, pos_t pos2,
                        bool unoriented_distance = false, const HandleGraph* graph=nullptr); 
//Maximum distance taking a pos instead of id/orientation/offset
size_t maximum_distance(const SnarlDistanceIndex& distance_index, pos_t pos1, pos_t pos2); 

//Fill in the index
//size_limit is a limit on the number of nodes in a snarl, after which the index won't store pairwise distances
void fill_in_distance_index(SnarlDistanceIndex* distance_index, const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit = 50000, bool only_top_level_chain_distances = false, bool silence_warnings=true);

//Fill in the temporary snarl record with distances
void populate_snarl_index(SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, 
    pair<SnarlDistanceIndex::temp_record_t, size_t> snarl_index, size_t size_limit, bool only_top_level_chain_distances, const HandleGraph* graph) ;

SnarlDistanceIndex::TemporaryDistanceIndex make_temporary_distance_index(const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit, bool only_top_level_chain_distances);

//Define wang_hash for net_handle_t's so that we can use a hash_map
template<> struct wang_hash<handlegraph::net_handle_t> {
public:
    inline size_t operator()(const net_handle_t& net_handle) const {
        return wang_hash_64(as_integer(net_handle));
    }
};

//Given an alignment to a graph and a range, find the set of nodes in the
//graph for which the minimum distance from the position to any position
//in the node is within the given distance range
//If look_forward is true, then start from the start of the path forward,
//otherwise start from the end going backward
void subgraph_in_distance_range(const SnarlDistanceIndex& distance_index, const Path& path, const HandleGraph* super_graph, size_t min_distance,
                                size_t max_distance, std::unordered_set<nid_t>& subgraph, bool look_forward);
///Helper for subgraph_in_distance_range
///Given starting handles in the super graph and the distances to each handle (including the start position and
//the first position in the handle), add all nodes within the distance range, excluding nodes in seen_nodes
//traversal_start is the node that we started the search from, because we can traverse it a second time but 
//we don't want to include a loop distance to any node after it
void subgraph_in_distance_range_walk_graph(const HandleGraph* super_graph, size_t min_distance, size_t max_distance,
                        std::unordered_set<nid_t>& subgraph, vector<pair<handle_t, size_t>>& start_nodes,
                        hash_set<pair<nid_t, bool>>& seen_nodes, const pair<nid_t, bool>& traversal_start);

//Helper function for subgraph_in_distance_range
//Given a node that is a child of a chain and the distance to it, walk forward from the
//current child and add nodes to search_start_nodes
void subgraph_in_distance_range_walk_across_chain (const SnarlDistanceIndex& distance_index, 
        const HandleGraph* super_graph,std::unordered_set<nid_t>& subgraph,
        net_handle_t current_node, size_t current_distance, 
        vector<pair<handle_t, size_t>>& search_start_nodes,
        hash_set<pair<nid_t, bool>>& seen_nodes, 
        const size_t& min_distane, const size_t& max_distance, bool checked_loop=false);


//Add nodes to the subgraph if they are near the path in the snarl tree
//Walks up the snarl tree from either end of the path, then takes everything
//in between them
void subgraph_containing_path_snarls(const SnarlDistanceIndex& distance_index, const HandleGraph* graph, const Path& path, std::unordered_set<nid_t>& subgraph);


//Helper function for subgraph_containing_path_snarls
//Add all the nodes in the parent to the subgraph
void add_descendants_to_subgraph(const SnarlDistanceIndex& distance_index, const net_handle_t& parent, std::unordered_set<nid_t>& subgraph);


}

#endif

//#define debug_distance_indexing
//#define debug_snarl_traversal
//#define debug_distances
//#define debug_subgraph
//#define debug_hub_label_build
//#define debug_hub_label_storage

#include "snarl_distance_index.hpp"

using namespace std;
using namespace handlegraph;
namespace vg {

size_t minimum_distance(const SnarlDistanceIndex& distance_index, pos_t pos1, pos_t pos2,
                        bool unoriented_distance, const HandleGraph* graph) {
    return distance_index.minimum_distance( get_id(pos1), get_is_rev(pos1), get_offset(pos1),
                                            get_id(pos2), get_is_rev(pos2), get_offset(pos2),
                                            unoriented_distance, graph, nullptr); 
}

size_t minimum_nontrivial_distance(const SnarlDistanceIndex& distance_index, pos_t pos1, pos_t pos2,
                                   size_t pos2_length, const HandleGraph* graph) {
    bool shifted = false;
    if (pos1 == pos2) {
        if (pos2_length == std::numeric_limits<size_t>::max()) {
            // If we don't know the length, we can get it from the graph
            pos2_length = distance_index.minimum_length(
                distance_index.get_node_net_handle(id(pos2)));
        }
        // Must shift one position to avoid self-distance of 0
        if (offset(pos1) == pos2_length) {
            // Shift ending pos backward (not safe to shift forward)
            get_offset(pos2)--;
        } else {
            // Shift starting position forward
            get_offset(pos1)++;
        }
        
        shifted = true;
    }

    size_t distance = minimum_distance(distance_index, pos1, pos2, false, graph);
    if (shifted && distance != std::numeric_limits<size_t>::max()) {
        // This loop is possible, so add back in the shift
        distance++;
    }

    return distance;
}

size_t maximum_distance(const SnarlDistanceIndex& distance_index, pos_t pos1, pos_t pos2) {
    return distance_index.maximum_distance( get_id(pos1), get_is_rev(pos1), get_offset(pos1),
                                            get_id(pos2), get_is_rev(pos2), get_offset(pos2)); 
}

void fill_in_distance_index(SnarlDistanceIndex* distance_index, const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit, bool only_top_level_chain_distances, bool silence_warnings) {
    distance_index->set_snarl_size_limit(size_limit);
    distance_index->set_only_top_level_chain_distances(only_top_level_chain_distances);

    //Build the temporary distance index from the graph
    SnarlDistanceIndex::TemporaryDistanceIndex temp_index = make_temporary_distance_index(graph, snarl_finder, size_limit, only_top_level_chain_distances);

    if (!silence_warnings && temp_index.use_oversized_snarls) {
        cerr << "warning: distance index uses oversized snarls, (the biggest has "
             << temp_index.most_oversized_snarl_size << " nodes), which may make mapping slow" << endl;
        cerr << "\ttry increasing --snarl-limit when building the distance index" << endl;
    }

    //And fill in the permanent distance index
    vector<const SnarlDistanceIndex::TemporaryDistanceIndex*> indexes;
    indexes.emplace_back(&temp_index);
    distance_index->get_snarl_tree_records(indexes, graph);
}

void subgraph_containing_path_snarls(const SnarlDistanceIndex& distance_index, const HandleGraph* graph, const Path& path, std::unordered_set<nid_t>& subgraph) {
    //Get the start and end of the path
    pos_t start_pos = initial_position(path);
    net_handle_t start_node = distance_index.get_node_net_handle(get_id(start_pos));
    subgraph.insert(get_id(start_pos));

    pos_t end_pos = final_position(path);
    net_handle_t end_node = distance_index.get_node_net_handle(get_id(end_pos));
    subgraph.insert(get_id(end_pos));

    //Get the lowest common ancestor
    pair<net_handle_t, bool> lowest_ancestor_bool = distance_index.lowest_common_ancestor(start_node, end_node);
    net_handle_t common_ancestor = lowest_ancestor_bool.first;
    
    
    if (distance_index.is_snarl(common_ancestor) || common_ancestor == start_node) {
        //If the lowest common ancestor is a snarl, just add the entire snarl

        add_descendants_to_subgraph(distance_index, common_ancestor, subgraph);

    } else if (distance_index.is_chain(common_ancestor)) {

        //Get the ancestors of the nodes that are children of the common ancestor
        net_handle_t ancestor1 = distance_index.canonical(distance_index.get_parent(start_node));
        while (ancestor1 != common_ancestor) {
            start_node = ancestor1;
            ancestor1 = distance_index.canonical(distance_index.get_parent(start_node));
        }
        net_handle_t ancestor2 = distance_index.canonical(distance_index.get_parent(end_node));
        while (ancestor2 != common_ancestor) {
            end_node = ancestor2;
            ancestor2 = distance_index.canonical(distance_index.get_parent(end_node));
        }
#ifdef debug_distance_indexing
        assert(ancestor1 == ancestor2);
#endif


        //Walk from one ancestor to the other and add everything in the chain
        net_handle_t current_child = distance_index.canonical(distance_index.is_ordered_in_chain(start_node, end_node) ? start_node : end_node);
        net_handle_t end_child = distance_index.canonical(distance_index.is_ordered_in_chain(start_node, end_node) ? end_node : start_node);
        if (distance_index.is_reversed_in_parent(current_child)) {
            current_child = distance_index.flip(current_child);
        }
        if (distance_index.is_reversed_in_parent(end_child)) {
            end_child = distance_index.flip(end_child);
        }

        add_descendants_to_subgraph(distance_index, current_child, subgraph);
        while (current_child != end_child) {
            distance_index.follow_net_edges(current_child, graph, false, [&](const net_handle_t& next) {
                add_descendants_to_subgraph(distance_index, next, subgraph);
                current_child = next;

            });
        }

    }
    
}


//Recursively add all nodes in parent to the subgraph
void add_descendants_to_subgraph(const SnarlDistanceIndex& distance_index, const net_handle_t& parent, std::unordered_set<nid_t>& subgraph) {
    if (distance_index.is_node(parent)) {
        subgraph.insert(distance_index.node_id(parent));
    } else {
        distance_index.for_each_child(parent, [&](const net_handle_t& child) {
            add_descendants_to_subgraph(distance_index, child, subgraph);
        });
    }
}
   


}

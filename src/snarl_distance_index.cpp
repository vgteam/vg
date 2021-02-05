//#define debugIndex
//#define debugDistance
//#define debugSubgraph

#include "snarl_distance_index.hpp"

using namespace std;
namespace vg {


//SnarlDistanceIndex::SnarlDistanceIndex() {
//    //TODO: Constructor
//}
///////////////////////////////////////////////////////////////////////////////////////////////////
//Implement the SnarlDecomposition's functions for moving around the snarl tree
//


bool SnarlDistanceIndex::get_root( net_handle_t& net) const {
    // The root is the first thing in the index, the traversal is tip to tip
    net = as_net_handle(1);
    return true;
}

bool SnarlDistanceIndex::is_root(const net_handle_t& net) const {
    return get_record_offset(net) == 0;
}

bool SnarlDistanceIndex::is_snarl(const net_handle_t& net) const {
    SnarlDistanceIndex::RecordType type = snarl_tree_record_t(net).get_record_type();
    return (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL);
}

bool SnarlDistanceIndex::is_chain(const net_handle_t& net) const {
    SnarlDistanceIndex::RecordType type = snarl_tree_record_t(net).get_record_type();
    return (type == CHAIN || type == DISTANCED_CHAIN);
}

bool SnarlDistanceIndex::is_node(const net_handle_t& net) const {
    SnarlDistanceIndex::RecordType type = snarl_tree_record_t(net).get_record_type();
    return type == NODE;
}


//TODO: I'm interpreting this to mean is this handle a node that is the boundary of a snarl. 
//If its parent is a chain, then it is
bool SnarlDistanceIndex::is_sentinel(const net_handle_t& net) const {
    if (snarl_tree_record_t(net).get_record_type() == NODE) {
        SnarlDistanceIndex::RecordType type = snarl_tree_record_t(snarl_tree_record_t(net).get_parent_record_pointer()).get_record_type();
        return type == CHAIN || type == DISTANCED_CHAIN;
    } else {
        return false ;
    }
}

//TODO: I'm going to interpret this as the node
//TODO: Also I think I also need the graph to get the id of the handle?
//net_handle_t SnarlDistanceIndex::get_net(const handle_t& handle) const{
//}
//handle_t SnarlDistanceIndex::get_handle(const net_handle_t& net) const{
//}

net_handle_t SnarlDistanceIndex::get_parent(const net_handle_t& child) const {
    //Get the pointer to the parent, and keep the connectivity of the current handle
    size_t parent_pointer = snarl_tree_record_t(child).get_parent_record_pointer();

    ConnectivityType child_connectivity = get_connectivity(child);
    //TODO: I"m going into the parent record here, which could be avoided if things knew what their parents were, but I think if you're doing this you'd later go into the parent anyway so it's probably fine
    RecordType parent_type = snarl_tree_record_t(parent_pointer).get_record_type();
    ConnectivityType parent_connectivity = START_END;
    if ((child_connectivity == START_END || child_connectivity == END_START) 
        && (parent_type == CHAIN  || parent_type == DISTANCED_CHAIN)) {
        //TODO: This also needs to take into account the orientation of the child, which I might be able to get around?
        parent_connectivity = child_connectivity;
    }

    return get_net_handle(parent_pointer, parent_connectivity);
}


}

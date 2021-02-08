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
    if (get_handle_type(child) == NODE_HANDLE && 
        (parent_type == ROOT || parent_type == SNARL || parent_type == DISTANCED_SNARL || 
         parent_type == SIMPLE_SNARL || parent_type == OVERSIZED_SNARL)) {
        //If this is a node and it's parent is not a chain, we want to pretend that its 
        //parent is a chain
        return get_net_handle(parent_pointer, parent_connectivity, CHAIN_HANDLE);
    }

    return get_net_handle(parent_pointer, parent_connectivity);
}

net_handle_t SnarlDistanceIndex::get_bound(const net_handle_t& snarl, bool get_end, bool face_in) const {
    id_t id = get_end ? snarl_tree_record_t(snarl).get_end_id() : snarl_tree_record_t(snarl).get_start_id();
    bool rev_in_parent = node_record_t(id).get_rev_in_parent();
    if (get_end) {
        rev_in_parent = !rev_in_parent;
    }
    if (!face_in){
        rev_in_parent = !rev_in_parent;
    }
    ConnectivityType connectivity = rev_in_parent ? END_START : START_END;
    return get_net_handle(id, connectivity);
}

net_handle_t SnarlDistanceIndex::flip(const net_handle_t& net) const {
    ConnectivityType old_connectivity = get_connectivity(net);
    ConnectivityType new_connectivity;
    if (old_connectivity == START_END) {
        new_connectivity = END_START;
    } else if (old_connectivity == START_TIP) {
        new_connectivity = TIP_START;
    } else if (old_connectivity = END_START) {
        new_connectivity = START_END;
    } else if (old_connectivity = END_TIP) {
        new_connectivity = TIP_END;
    } else if (old_connectivity = TIP_START) {
        new_connectivity = START_TIP;
    } else if (old_connectivity = TIP_END) {
        new_connectivity = END_TIP;
    } else {
        new_connectivity = old_connectivity;
    }
    return get_net_handle(get_record_offset(net), new_connectivity);
}

net_handle_t SnarlDistanceIndex::canonical(const net_handle_t& net) const {
    snarl_tree_record_t record(net);
    ConnectivityType connectivity;
    if (record.is_start_end_connected()) {
        connectivity = START_END;
    } else if (record.is_start_tip_connected()) {
        connectivity = START_TIP;
    } else if (record.is_end_tip_connected()) {
        connectivity = END_TIP;
    } else if (record.is_start_start_connected()) {
        connectivity = START_START;
    } else if (record.is_end_end_connected()) {
        connectivity = END_END;
    } else if (record.is_tip_tip_connected()) {
        connectivity = TIP_TIP;
    } else {
        throw runtime_error("error: This node has no connectivity");
    }
    return get_net_handle(get_record_offset(net), connectivity);
}

SnarlDecomposition::endpoint_t SnarlDistanceIndex::starts_at(const net_handle_t& traversal) const {
    ConnectivityType connectivity = get_connectivity(traversal);
    if (connectivity == START_START || connectivity == START_END || connectivity == START_TIP ){
        return START;
    } else if (connectivity == END_START || connectivity == END_END || connectivity == END_TIP ){
        return END;
    } else if (connectivity == TIP_START || connectivity == TIP_END || connectivity == TIP_TIP ){
        return TIP;
    } else {
        throw runtime_error("error: This node has no connectivity");
    }
}
SnarlDecomposition::endpoint_t SnarlDistanceIndex::ends_at(const net_handle_t& traversal) const {
    ConnectivityType connectivity = get_connectivity(traversal);
    if (connectivity == START_START || connectivity == END_START || connectivity == TIP_START ){
        return START;
    } else if (connectivity == START_END || connectivity == END_END || connectivity == TIP_END ){
        return END;
    } else if (connectivity == START_TIP || connectivity == END_TIP || connectivity == TIP_TIP ){
        return TIP;
    } else {
        throw runtime_error("error: This node has no connectivity");
    }
}

//TODO: I'm also allowing this for the root
bool SnarlDistanceIndex::for_each_child_impl(const net_handle_t& traversal, const std::function<bool(const net_handle_t&)>& iteratee) const {
    //What is this according to the snarl tree
    RecordType type = snarl_tree_record_t(traversal).get_record_type();
    //What is this according to the handle 
    //(could be a trivial chain but actually a node according to the snarl tree)
    HandleType handle_type = get_handle_type(traversal);
    if (type == SNARL || type == DISTANCED_SNARL || type == OVERSIZED_SNARL) {
        snarl_record_t snarl_record(traversal);
        return snarl_record.for_each_child(iteratee);
    } else if (type == CHAIN || type == DISTANCED_CHAIN) {
        chain_record_t chain_record(traversal);
        return chain_record.for_each_child(iteratee);
    } else if (type == ROOT) {
        root_record_t root_record(traversal);
        return root_record.for_each_child(iteratee);
    } else if ((type == NODE || type == DISTANCED_NODE) && handle_type == CHAIN_HANDLE) {
        //This is actually a node but we're pretending it's a chain
        chain_record_t trivial_chain_record_t(traversal);
        return chain_record_t(traversal).for_each_child(iteratee);
    } else {
        throw runtime_error("error: Looking for children of a node");
    }
    
}

}

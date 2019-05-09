//
//  packed_graph.cpp
//

#include "packed_graph.hpp"

#include <handlegraph/util.hpp>
#include <atomic>

namespace vg {

    using namespace handlegraph;
    
    /// Define all of the static class variables
    
    const double PackedGraph::defrag_factor = .2;
    
    const size_t PackedGraph::PAGE_WIDTH = 128;
    
    const size_t PackedGraph::GRAPH_RECORD_SIZE = 2;
    const size_t PackedGraph::GRAPH_START_EDGES_OFFSET = 0;
    const size_t PackedGraph::GRAPH_END_EDGES_OFFSET = 1;
    
    const size_t PackedGraph::SEQ_START_RECORD_SIZE = 1;
    
    const size_t PackedGraph::SEQ_LENGTH_RECORD_SIZE = 1;
    
    const size_t PackedGraph::EDGE_RECORD_SIZE = 2;
    const size_t PackedGraph::EDGE_TRAV_OFFSET = 0;
    const size_t PackedGraph::EDGE_NEXT_OFFSET = 1;
    
    const size_t PackedGraph::NODE_MEMBER_RECORD_SIZE = 1;
    
    const size_t PackedGraph::MEMBERSHIP_ID_RECORD_SIZE = 1;
    const size_t PackedGraph::MEMBERSHIP_OFFSET_RECORD_SIZE = 1;
    const size_t PackedGraph::MEMBERSHIP_NEXT_RECORD_SIZE = 1;
    
    const size_t PackedGraph::STEP_RECORD_SIZE = 1;
    
    const size_t PackedGraph::PATH_RECORD_SIZE = 2;
    const size_t PackedGraph::PATH_PREV_OFFSET = 0;
    const size_t PackedGraph::PATH_NEXT_OFFSET = 1;
    
    PackedGraph::PackedGraph() :
        graph_iv(PAGE_WIDTH),
        seq_start_iv(PAGE_WIDTH),
        edge_lists_iv(PAGE_WIDTH),
        path_membership_node_iv(PAGE_WIDTH),
        path_membership_next_iv(PAGE_WIDTH),
        path_membership_offset_iv(PAGE_WIDTH),
        path_membership_id_iv(PAGE_WIDTH) {
        
    }
    
    PackedGraph::PackedGraph(istream& in) :
        graph_iv(PAGE_WIDTH),
        seq_start_iv(PAGE_WIDTH),
        edge_lists_iv(PAGE_WIDTH),
        path_membership_node_iv(PAGE_WIDTH),
        path_membership_next_iv(PAGE_WIDTH),
        path_membership_offset_iv(PAGE_WIDTH),
        path_membership_id_iv(PAGE_WIDTH) {
        
        deserialize(in);
    }
    
    PackedGraph::~PackedGraph() {
        
    }
    
    void PackedGraph::serialize(ostream& out) const {
        sdsl::write_member(max_id, out);
        sdsl::write_member(min_id, out);
        
        graph_iv.serialize(out);
        seq_start_iv.serialize(out);
        seq_length_iv.serialize(out);
        edge_lists_iv.serialize(out);
        id_to_graph_iv.serialize(out);
        seq_iv.serialize(out);
        
        path_membership_node_iv.serialize(out);
        path_membership_id_iv.serialize(out);
        path_membership_offset_iv.serialize(out);
        path_membership_next_iv.serialize(out);
        
        sdsl::write_member(paths.size(), out);
        for (const PackedPath& path : paths) {
            sdsl::write_member(path.name, out);
            sdsl::write_member(path.is_deleted, out);
            sdsl::write_member(path.is_circular, out);
            sdsl::write_member(path.head, out);
            sdsl::write_member(path.tail, out);
            sdsl::write_member(path.deleted_step_records, out);
            path.links_iv.serialize(out);
            path.steps_iv.serialize(out);
        }
        // note: path_id can be reconstructed from the paths
        sdsl::write_member(deleted_node_records, out);
        sdsl::write_member(deleted_edge_records, out);
        sdsl::write_member(deleted_membership_records, out);
    }
    
    void PackedGraph::deserialize(istream& in) {
        sdsl::read_member(max_id, in);
        sdsl::read_member(min_id, in);
        
        graph_iv.deserialize(in);
        seq_start_iv.deserialize(in);
        seq_length_iv.deserialize(in);
        edge_lists_iv.deserialize(in);
        id_to_graph_iv.deserialize(in);
        seq_iv.deserialize(in);
        
        path_membership_node_iv.deserialize(in);
        path_membership_id_iv.deserialize(in);
        path_membership_offset_iv.deserialize(in);
        path_membership_next_iv.deserialize(in);
        
        size_t num_paths;
        sdsl::read_member(num_paths, in);
        for (size_t i = 0; i < num_paths; i++) {
            string name;
            sdsl::read_member(name, in);
            paths.emplace_back(name, false); // dummy circularity here, real in a few lines
            PackedPath& path = paths.back();
            sdsl::read_member(path.is_deleted, in);
            sdsl::read_member(path.is_circular, in);
            sdsl::read_member(path.head, in);
            sdsl::read_member(path.tail, in);
            sdsl::read_member(path.deleted_step_records, in);
            path.links_iv.deserialize(in);
            path.steps_iv.deserialize(in);
        }
        
        // reconstruct the path_id mapping
        for (int64_t i = 0; i < paths.size(); i++) {
            if (!paths[i].is_deleted) {
                path_id[paths[i].name] = i;
            }
        }
        
        sdsl::read_member(deleted_node_records, in);
        sdsl::read_member(deleted_edge_records, in);
        sdsl::read_member(deleted_membership_records, in);
    }
    
    
    size_t PackedGraph::new_node_record(id_t node_id) {
        
        size_t next_g_iv_idx = graph_iv.size();
        
        // no edges yet, null pointer for linked list
        graph_iv.append(0);
        graph_iv.append(0);
        
        // record the sequence interval
        seq_start_iv.append(0);
        seq_length_iv.append(0);
        
        // initialize an empty path membership list
        path_membership_node_iv.append(0);
        
        // expand the ID vector's dimensions so it can handle the full ID interval
        if (id_to_graph_iv.empty()) {
            id_to_graph_iv.append_back(0);
        }
        else {
            for (int64_t i = node_id; i < min_id; i++) {
                id_to_graph_iv.append_front(0);
            }
            for (int64_t i = id_to_graph_iv.size(); i <= node_id - min_id; i++) {
                id_to_graph_iv.append_back(0);
            }
        }
        
        // update the min and max ID
        max_id = std::max(node_id, max_id);
        min_id = std::min(node_id, min_id);
        
        // record the mapping of the ID to the graph record
        id_to_graph_iv.set(node_id - min_id, graph_iv.size() / GRAPH_RECORD_SIZE);
        
        return next_g_iv_idx;
    }
    
    handle_t PackedGraph::create_handle(const string& sequence) {
        return create_handle(sequence, max_id + 1);
    }
    
    handle_t PackedGraph::create_handle(const string& sequence, const id_t& id) {
        
        if (id >= min_id && id < min_id + id_to_graph_iv.size()) {
            if (id_to_graph_iv.get(id - min_id) != 0) {
                cerr << "error:[PackedGraph] tried to create a node with ID " << id << ", but this ID already belongs to a different node" << endl;
                exit(1);
            }
        }
        
        size_t g_iv_idx = new_node_record(id);
        seq_start_iv.set(graph_index_to_seq_start_index(g_iv_idx), seq_iv.size());
        seq_length_iv.set(graph_index_to_seq_len_index(g_iv_idx), sequence.size());
        
        // encode the sequence interval
        for (size_t i = 0; i < sequence.size(); i++) {
            seq_iv.append(encode_nucleotide(sequence[i]));
        }
        
        return get_handle(id);
    }
    
    void PackedGraph::create_edge(const handle_t& left, const handle_t& right) {
        
        // look for the edge
        bool add_edge = follow_edges(left, false, [&](const handle_t& next) {
            return next != right;
        });
        
        // don't duplicate it
        if (!add_edge) {
            return;
        }
        
        // get the location of the edge list pointer in the graph vector
        size_t g_iv_left = graph_iv_index(left) + (get_is_reverse(left) ?
                                                   GRAPH_START_EDGES_OFFSET :
                                                   GRAPH_END_EDGES_OFFSET);
        size_t g_iv_right = graph_iv_index(right) + (get_is_reverse(right) ?
                                                     GRAPH_END_EDGES_OFFSET :
                                                     GRAPH_START_EDGES_OFFSET);
        
        // add a new linked list node pointing to the rest of the list
        edge_lists_iv.append(encode_traversal(right));
        edge_lists_iv.append(graph_iv.get(g_iv_left));
        // make this new node the head
        graph_iv.set(g_iv_left, edge_lists_iv.size() / EDGE_RECORD_SIZE);
        
        // don't double add a reversing self edge
        if (g_iv_left == g_iv_right) {
            return;
        }
        
        // add a new linked list node pointing to the rest of the list
        edge_lists_iv.append(encode_traversal(flip(left)));
        edge_lists_iv.append(graph_iv.get(g_iv_right));
        // make this new node the head
        graph_iv.set(g_iv_right, edge_lists_iv.size() / EDGE_RECORD_SIZE);
    }
    
    bool PackedGraph::has_node(id_t node_id) const {
        if (node_id < min_id || node_id - min_id >= id_to_graph_iv.size()) {
            return false;
        }
        else {
            return id_to_graph_iv.get(node_id - min_id) != 0;
        }
    }
    
    handle_t PackedGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        return handlegraph::number_bool_packing::pack(node_id, is_reverse);
    }
    
    id_t PackedGraph::get_id(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_number(handle);
    }
    
    bool PackedGraph::get_is_reverse(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_bit(handle);;
    }
    
    handle_t PackedGraph::flip(const handle_t& handle) const {
        return handlegraph::number_bool_packing::toggle_bit(handle);
    }
    
    size_t PackedGraph::get_length(const handle_t& handle) const {
        return seq_length_iv.get(graph_index_to_seq_len_index(graph_iv_index(handle)));
    }
    
    string PackedGraph::get_sequence(const handle_t& handle) const {
        size_t g_iv_index = graph_iv_index(handle);
        size_t seq_start = seq_start_iv.get(graph_index_to_seq_start_index(g_iv_index));
        size_t seq_len = seq_length_iv.get(graph_index_to_seq_len_index(g_iv_index));
        string seq(seq_len, 'N');
        for (size_t i = 0; i < seq_len; i++) {
            seq[i] = decode_nucleotide(seq_iv.get(seq_start + i));
        }
        return get_is_reverse(handle) ? reverse_complement(seq) : seq;
    }
    
    /*
    void PackedGraph::swap_handles(const handle_t& a, const handle_t& b) {
        // TODO: this doesn't actually affect the traversal order
        // We either need to encode the IDs here or use different method to keep track
        // of deleted nodes
        
        size_t g_iv_index_a = graph_iv_index(a);
        size_t g_iv_index_b = graph_iv_index(b);
        
        int64_t id_a = get_id(a);
        int64_t id_b = get_id(b);
        
        // swap all the records that occur in the same order as each other in the vectors
        for (size_t i = 0; i < GRAPH_RECORD_SIZE; i++) {
            uint64_t val = graph_iv.get(g_iv_index_a + i);
            graph_iv.set(g_iv_index_a + i, graph_iv.get(g_iv_index_b + i));
            graph_iv.set(g_iv_index_b + i, val);
        }
        for (size_t i = 0; i < SEQ_LENGTH_RECORD_SIZE; i++) {
            uint64_t val = seq_length_iv.get(graph_index_to_seq_len_index(g_iv_index_a) + i);
            seq_length_iv.set(graph_index_to_seq_len_index(g_iv_index_a) + i,
                              seq_length_iv.get(graph_index_to_seq_len_index(g_iv_index_b) + i));
            seq_length_iv.set(graph_index_to_seq_len_index(g_iv_index_b) + i, val);
        }
        for (size_t i = 0; i < SEQ_START_RECORD_SIZE; i++) {
            uint64_t val = seq_start_iv.get(graph_index_to_seq_start_index(g_iv_index_a) + i);
            seq_start_iv.set(graph_index_to_seq_start_index(g_iv_index_a) + i,
                              seq_start_iv.get(graph_index_to_seq_start_index(g_iv_index_b) + i));
            seq_start_iv.set(graph_index_to_seq_start_index(g_iv_index_b) + i, val);
        }
        for (size_t i = 0; i < NODE_MEMBER_RECORD_SIZE; i++) {
            uint64_t val = path_membership_node_iv.get(graph_index_to_node_member_index(g_iv_index_a) + i);
            path_membership_node_iv.set(graph_index_to_node_member_index(g_iv_index_a) + i,
                             path_membership_node_iv.get(graph_index_to_node_member_index(g_iv_index_b) + i));
            path_membership_node_iv.set(graph_index_to_node_member_index(g_iv_index_b) + i, val);
        }
        
        // swap the pointers from the ID vector
        id_to_graph_iv.set(id_a - min_id, g_iv_index_b / GRAPH_RECORD_SIZE + 1);
        id_to_graph_iv.set(id_b - min_id, g_iv_index_a / GRAPH_RECORD_SIZE + 1);
    }
    */
    
    bool PackedGraph::follow_edges_impl(const handle_t& handle, bool go_left,
                                        const std::function<bool(const handle_t&)>& iteratee) const {
        // toward start = true, toward end = false
        bool direction = get_is_reverse(handle) != go_left;
        // get the head of the linked list from the graph vector
        size_t edge_idx = graph_iv.get(graph_iv_index(handle)
                                       + (direction ? GRAPH_START_EDGES_OFFSET : GRAPH_END_EDGES_OFFSET));
        // traverse the linked list as long as directed
        bool keep_going = true;
        while (edge_idx && keep_going) {
            
            handle_t edge_target = decode_traversal(get_edge_target(edge_idx));
            if (go_left) {
                // match the orientation encoding
                edge_target = flip(edge_target);
            }
            
            keep_going = iteratee(edge_target);
            edge_idx = get_next_edge_index(edge_idx);
        }
        
        return keep_going;
    }
    
    size_t PackedGraph::node_size(void) const {
        return graph_iv.size() / GRAPH_RECORD_SIZE - deleted_node_records;
    }
    
    id_t PackedGraph::min_node_id(void) const {
        return min_id;
    }
    
    id_t PackedGraph::max_node_id(void) const {
        return max_id;
    }
    
    bool PackedGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee,
                                           bool parallel) const {
        
        if (parallel) {
            // TODO: would task based parallelism be better?
            atomic<bool> keep_going(true);
#pragma omp parallel for
            for (size_t i = 0; i < id_to_graph_iv.size(); i++) {
                if (keep_going && id_to_graph_iv.get(i)) {
                    if (!iteratee(get_handle(i + min_id))) {
                        keep_going = false;
                    }
                }
            }
            return keep_going;
        }
        else {
            for (size_t i = 0; i < id_to_graph_iv.size(); i++) {
                if (id_to_graph_iv.get(i)) {
                    if (!iteratee(get_handle(i + min_id))) {
                        return false;
                    }
                }
            }
            return true;
        }
        
    }
    
    char PackedGraph::get_base(const handle_t& handle, size_t index) const {
        size_t g_iv_index = graph_iv_index(handle);
        size_t seq_start = seq_start_iv.get(graph_index_to_seq_start_index(g_iv_index));
        if (get_is_reverse(handle)) {
            size_t seq_len = seq_length_iv.get(graph_index_to_seq_len_index(g_iv_index));
            return reverse_complement(decode_nucleotide(seq_iv.get(seq_start + seq_len - index - 1)));
        }
        else {
            return decode_nucleotide(seq_iv.get(seq_start + index));
        }
    }
    
    string PackedGraph::get_subsequence(const handle_t& handle, size_t index, size_t size) const {
        size_t g_iv_index = graph_iv_index(handle);
        size_t seq_start = seq_start_iv.get(graph_index_to_seq_start_index(g_iv_index));
        size_t seq_len = seq_length_iv.get(graph_index_to_seq_len_index(g_iv_index));
        
        size = min(size, seq_len - index);
        size_t subseq_start = get_is_reverse(handle) ? seq_start + seq_len - size - index : seq_start + index;
        
        string subseq(size, 'N');
        for (size_t i = 0; i < size; i++) {
            subseq[i] = decode_nucleotide(seq_iv.get(subseq_start + i));
        }
        return get_is_reverse(handle) ? reverse_complement(subseq) : subseq;
    }
    
    handle_t PackedGraph::apply_orientation(const handle_t& handle) {
        
        if (get_is_reverse(handle)) {
            size_t g_iv_idx = graph_iv_index(handle);
            
            // swap the edge lists
            size_t tmp = graph_iv.get(g_iv_idx + GRAPH_START_EDGES_OFFSET);
            graph_iv.set(g_iv_idx + GRAPH_START_EDGES_OFFSET, graph_iv.get(g_iv_idx + GRAPH_END_EDGES_OFFSET));
            graph_iv.set(g_iv_idx + GRAPH_END_EDGES_OFFSET, tmp);
            
            // reverse the orientation in the backward pointers
            for (size_t orientation : {true, false}) {
                // iterate down the entire edge list
                handle_t looking_for = orientation ? handle : flip(handle);
                size_t edge_list_idx = graph_iv.get(g_iv_idx + (orientation ? GRAPH_START_EDGES_OFFSET : GRAPH_END_EDGES_OFFSET));
                while (edge_list_idx) {
                    handle_t target = decode_traversal(get_edge_target(edge_list_idx));
                    size_t backward_edge_idx = graph_iv.get(graph_iv_index(target) + (get_is_reverse(target) ?
                                                                                      GRAPH_END_EDGES_OFFSET :
                                                                                      GRAPH_START_EDGES_OFFSET));
                    
                    while (backward_edge_idx) {
                        handle_t backward_edge_target = decode_traversal(get_edge_target(backward_edge_idx));
                        if (backward_edge_target == looking_for) {
                            set_edge_target(backward_edge_idx, flip(backward_edge_target));
                            break;
                        }
                        backward_edge_idx = get_next_edge_index(backward_edge_idx);
                    }
                    
                    edge_list_idx = get_next_edge_index(edge_list_idx);
                }
            }
            
            // reverse complement the sequence in place
            
            size_t seq_start = seq_start_iv.get(graph_index_to_seq_start_index(g_iv_idx));
            size_t seq_len = seq_length_iv.get(graph_index_to_seq_len_index(g_iv_idx));
            
            for (size_t i = 0, stop = seq_len / 2; i < stop; i++) {
                size_t j = seq_start + seq_len - i - 1;
                size_t k = seq_start + i;
                uint64_t base = seq_iv.get(k);
                seq_iv.set(k, complement_encoded_nucleotide(seq_iv.get(j)));
                seq_iv.set(j, complement_encoded_nucleotide(base));
            }
            if (seq_len % 2) {
                size_t j = seq_start + seq_len / 2;
                seq_iv.set(j, complement_encoded_nucleotide(seq_iv.get(j)));
            }
            
            // reverse the orientation of the node on all paths
            
            size_t path_membership = path_membership_node_iv.get(graph_index_to_node_member_index(g_iv_idx));
            while (path_membership) {
                
                // get the path that this membership record is on
                PackedPath& packed_path = paths[get_membership_path(path_membership)];
                
                // access and flip the step on the path
                size_t occ_idx = get_membership_step(path_membership);
                handle_t occ = decode_traversal(get_step_trav(packed_path, occ_idx));
                set_step_trav(packed_path, occ_idx, encode_traversal(flip(occ)));
                
                // move to the next membership record
                path_membership = get_next_membership(path_membership);
            }
            
            // the ID is preserved, we just need to need to return a forward version
            return flip(handle);
        }
        else {
            // it's already the way we want it
            return handle;
        }
    }
    
    std::vector<handle_t> PackedGraph::divide_handle(const handle_t& handle,
                                                     const std::vector<size_t>& offsets) {
        // put the offsets in forward orientation to simplify subsequent steps
        vector<size_t> forward_offsets = offsets;
        size_t node_length = get_length(handle);
        if (get_is_reverse(handle)) {
            for (size_t& off : forward_offsets) {
                off = node_length - off;
            }
        }
        
        // we will also build the return value in forward orientation
        handle_t forward_handle = forward(handle);
        vector<handle_t> return_val{forward_handle};
        size_t g_iv_idx = graph_iv_index(forward_handle);
        
        // offsets in the sequence vector will be measured relative to the first position of
        // the current handle
        size_t first_start = seq_start_iv.get(graph_index_to_seq_start_index(g_iv_idx));
        
        // we record the the edges out of this node so they can be transferred onto the final
        // node in the split
        size_t end_edges = graph_iv.get(g_iv_idx + GRAPH_END_EDGES_OFFSET);
        
        // init trackers for the previous iteration
        size_t last_offset = 0;
        id_t prev_id = get_id(forward_handle);
        for (const size_t& off : forward_offsets) {
            
            id_t next_id = max_id + 1;
            size_t new_g_iv_idx = new_node_record(next_id);
            // seq start
            seq_start_iv.set(graph_index_to_seq_start_index(new_g_iv_idx), first_start + off);
            
            return_val.push_back(get_handle(next_id, false));
            
            // now let's do what we still need to on the previous node
            // set the previous node's length based on the current offset
            seq_length_iv.set(graph_index_to_seq_len_index(g_iv_idx), off - last_offset);
            
            // create an edge forward onto the new node
            edge_lists_iv.append(encode_traversal(get_handle(next_id)));
            edge_lists_iv.append(0);
            // add the edge onto the previous node
            graph_iv.set(g_iv_idx + GRAPH_END_EDGES_OFFSET, edge_lists_iv.size() / EDGE_RECORD_SIZE);
            
            // create an edge backward to the previous node
            edge_lists_iv.append(encode_traversal(get_handle(prev_id, true)));
            edge_lists_iv.append(0);
            // add the edge backwards to the current node
            graph_iv.set(new_g_iv_idx + GRAPH_START_EDGES_OFFSET, edge_lists_iv.size() / EDGE_RECORD_SIZE);
            
            // switch the previous node variables onto the new node node
            g_iv_idx = new_g_iv_idx;
            prev_id = next_id;
            last_offset = off;
        }
        
        // set final node's length to the remaining sequence
        seq_length_iv.set(graph_index_to_seq_len_index(g_iv_idx), node_length - last_offset);
        
        // point the final node's end edges to the original node's end edges
        graph_iv.set(g_iv_idx + GRAPH_END_EDGES_OFFSET, end_edges);
        
        // update the back edges onto the final node in the division
        // note: we don't need to do the same for the first node since its ID stays the same
        size_t next_edge_idx = end_edges;
        handle_t looking_for = flip(handle);
        while (next_edge_idx) {
            handle_t target = decode_traversal(get_edge_target(next_edge_idx));
            size_t backward_edge_idx = graph_iv.get(graph_iv_index(target) + (get_is_reverse(target) ?
                                                                              GRAPH_END_EDGES_OFFSET :
                                                                              GRAPH_START_EDGES_OFFSET));
            
            while (backward_edge_idx) {
                handle_t backward_edge_target = decode_traversal(get_edge_target(backward_edge_idx));
                if (backward_edge_target == looking_for) {
                    set_edge_target(backward_edge_idx, flip(return_val.back()));
                    break;
                }
                backward_edge_idx = get_next_edge_index(backward_edge_idx);
            }
            
            next_edge_idx = get_next_edge_index(next_edge_idx);
        }
        
        // update all of the occurrences on paths
        
        size_t path_membership = path_membership_node_iv.get(graph_index_to_node_member_index(graph_iv_index(handle)));
        while (path_membership) {
            
            // get the path that this membership record is on
            PackedPath& packed_path = paths[get_membership_path(path_membership)];
            
            // split up the occurrence on the path
            size_t occ_idx = get_membership_step(path_membership);
            bool path_trav_rev = get_is_reverse(decode_traversal(get_step_trav(packed_path, occ_idx)));
            // make new occurrence records for the divided segments (except the first, which stays
            // in place)
            vector<size_t> divided_trav_offsets{occ_idx};
            for (size_t i = 1; i < return_val.size(); i++) {
                // the new traversals will have the same strandedness as the original occurrence
                packed_path.steps_iv.append(encode_traversal(path_trav_rev ? flip(return_val[i]) : return_val[i]));
                packed_path.links_iv.append(0);
                packed_path.links_iv.append(0);
                
                divided_trav_offsets.push_back(packed_path.steps_iv.size() / STEP_RECORD_SIZE);
                
                // record the membership of this node in this path
                size_t node_member_idx = graph_index_to_node_member_index(graph_iv_index(return_val[i]));
                path_membership_id_iv.append(get_membership_path(path_membership));
                path_membership_offset_iv.append(packed_path.steps_iv.size() / STEP_RECORD_SIZE);
                path_membership_next_iv.append(path_membership_node_iv.get(node_member_idx));
                
                // make this new membership record the head of the linked list
                path_membership_node_iv.set(node_member_idx, path_membership_next_iv.size() / MEMBERSHIP_NEXT_RECORD_SIZE);
            }
            
            if (path_trav_rev) {
                // update connection to the original node that should now be to the final divided node
                size_t other_idx = get_step_prev(packed_path, occ_idx);
                if (other_idx != 0) {
                    set_step_next(packed_path, other_idx, divided_trav_offsets.back());
                    set_step_prev(packed_path, divided_trav_offsets.back(), other_idx);
                }
                
                // add connections between the divided segments
                for (size_t i = 1; i < divided_trav_offsets.size(); i++) {
                    set_step_prev(packed_path, divided_trav_offsets[i - 1], divided_trav_offsets[i]);
                    set_step_next(packed_path, divided_trav_offsets[i], divided_trav_offsets[i - 1]);
                }
            }
            else {
                // update connection to the original node that should now be to the final divided node
                size_t other_idx = get_step_next(packed_path, occ_idx);
                if (other_idx != 0) {
                    set_step_prev(packed_path, other_idx, divided_trav_offsets.back());
                    set_step_next(packed_path, divided_trav_offsets.back(), other_idx);
                }
                
                // add connections between the divided segments
                for (size_t i = 1; i < divided_trav_offsets.size(); i++) {
                    set_step_next(packed_path, divided_trav_offsets[i - 1], divided_trav_offsets[i]);
                    set_step_prev(packed_path, divided_trav_offsets[i], divided_trav_offsets[i - 1]);
                }
            }
            
            // move to the next membership record
            path_membership = get_next_membership(path_membership);
        }
        
        if (get_is_reverse(handle)) {
            // reverse the vector to the orientation of the input handle
            reverse(return_val.begin(), return_val.end());
            for (handle_t& ret_handle : return_val) {
                ret_handle = flip(ret_handle);
            }
        }
        
        return return_val;
    }
    
    void PackedGraph::destroy_handle(const handle_t& handle) {
        
        // remove the back-references to the edges
        follow_edges(handle, false, [&](const handle_t& next) {
            remove_edge_reference(flip(next), flip(handle));
            
            // we don't actually bother removing the reference, but we will also consider
            // the edge on the deleting node to be deleted and hence count it up here
            ++deleted_edge_records;
            return true;
        });
        follow_edges(handle, true, [&](const handle_t& prev) {
            
            remove_edge_reference(prev, handle);
            
            // we don't actually bother removing the reference, but we will also consider
            // the edge on the deleting node to be deleted and hence count it up here
            ++deleted_edge_records;
            return true;
        });
        
        // TODO: should we update the path membership for the node? deleting a node on a path
        // is undefined behavior, so we don't *need* to...
        
        // remove the reference to the node
        id_to_graph_iv.set(get_id(handle) - min_id, 0);
        
        ++deleted_node_records;
        
        // maybe reallocate to address fragmentation
        defragment();
    }
    
    void PackedGraph::remove_edge_reference(const handle_t& on, const handle_t& to) {
        
        // Note: this function assumes that edge ref exists, crashes otherwise
        
        size_t g_iv_idx = graph_iv_index(on) + (get_is_reverse(on)
                                                ? GRAPH_START_EDGES_OFFSET
                                                : GRAPH_END_EDGES_OFFSET);

        size_t edge_list_idx = graph_iv.get(g_iv_idx);
        
        if (decode_traversal(get_edge_target(edge_list_idx)) == to) {
            // the edge back to the deleting node is the first in the list, so we need
            // to update the head
            graph_iv.set(g_iv_idx, get_next_edge_index(edge_list_idx));
        }
        else {
            // we need to traverse down the list and to find the edge back
            size_t prev_edge_list_idx = edge_list_idx;
            edge_list_idx = get_next_edge_index(edge_list_idx);
            while (decode_traversal(get_edge_target(edge_list_idx)) != to) {
                prev_edge_list_idx = edge_list_idx;
                edge_list_idx = get_next_edge_index(edge_list_idx);
            }
            // skip over this edge in this linked list
            edge_lists_iv.set((prev_edge_list_idx - 1) * EDGE_RECORD_SIZE + EDGE_NEXT_OFFSET,
                              get_next_edge_index(edge_list_idx));
        }
        ++deleted_edge_records;
    }
    
    void PackedGraph::destroy_edge(const handle_t& left, const handle_t& right) {
        remove_edge_reference(left, right);
        remove_edge_reference(flip(right), flip(left));
        defragment();
    }
    
    void PackedGraph::defragment_path(PackedPath& path, bool force) {
        
        // we don't want to defrag deleted paths since they have already been cleared
        if (path.is_deleted) {
            return;
        }
        
        // have we either deleted a lot of steps or forced a defrag?
        if (path.deleted_step_records > defrag_factor * (path.steps_iv.size() / PATH_RECORD_SIZE) || force) {
            
            if (path.head != 0) {
                
                // the path is non-empty, so we need to straighten it out and reallocate it
                PagedVector new_steps_iv(PAGE_WIDTH);
                PagedVector new_links_iv(PAGE_WIDTH);
                
                // we will need to record the translation between path steps so we can update memberships later
                PagedVector offset_translator(PAGE_WIDTH);
                offset_translator.resize(path.steps_iv.size() / STEP_RECORD_SIZE + 1);
                
                new_links_iv.reserve(path.links_iv.size() - path.deleted_step_records * PATH_RECORD_SIZE);
                new_steps_iv.reserve(path.steps_iv.size() - path.deleted_step_records * STEP_RECORD_SIZE);
                
                bool first_iter = true;
                size_t copying_from = path.head;
                size_t prev = 0;
                while (copying_from != 0 && (first_iter || copying_from != path.head)) {
                    // make a new record
                    new_steps_iv.append(get_step_trav(path, copying_from));
                    new_links_iv.append(prev);
                    new_links_iv.append(0);
                    
                    size_t here = new_steps_iv.size() / STEP_RECORD_SIZE;
                    
                    // record the correspondance between the old
                    offset_translator.set(copying_from, here);
                    
                    // update the point on the previous node
                    if (prev != 0) {
                        new_links_iv.set(new_links_iv.size() - 2 * PATH_RECORD_SIZE + PATH_NEXT_OFFSET, here);
                    }
                    
                    prev = here;
                    copying_from = get_step_next(path, copying_from);
                    first_iter = false;
                }
                
                // add the looping connection if this is a circular path
                if (path.is_circular) {
                    new_links_iv.set(new_links_iv.size() - PATH_RECORD_SIZE + PATH_NEXT_OFFSET, 1);
                    new_links_iv.set(PATH_PREV_OFFSET, new_links_iv.size() / PATH_RECORD_SIZE);
                }
                
                path.links_iv = move(new_links_iv);
                path.steps_iv = move(new_steps_iv);
                
                // update the head and tail of the newly allocated path
                path.head = 1;
                path.tail = prev;
                
                // retrieve the ID of this path so we can match it to membership records
                int64_t path_id_here = path_id.at(path.name);
                
                // now we need to iterate over each node on the path exactly one time to update its membership
                // records (even if the node occurs multiple times on this path), so we will use a bit deque
                // indexed by node_id - min_id to flag nodes as either translated or untranslated
                PackedDeque id_translated;
                id_translated.append_back(0);
                id_t min_translated_id = get_id(decode_traversal(get_step_trav(path, path.head)));
                
                first_iter = true;
                for (size_t here = path.head; here != 0 && (here != path.head || first_iter); here = get_step_next(path, here)) {
                    
                    handle_t handle = decode_traversal(get_step_trav(path, here));
                    id_t step_node_id = get_id(handle);
                    
                    // expand the bounds of the deque as necessary to be able to index by ID
                    if (step_node_id < min_translated_id) {
                        for (id_t i = step_node_id; i < min_translated_id; ++i) {
                            id_translated.append_front(0);
                        }
                        min_translated_id = step_node_id;
                    }
                    else if (step_node_id >= min_translated_id + id_translated.size()) {
                        for (id_t i = min_translated_id + id_translated.size(); i <= step_node_id; ++i) {
                            id_translated.append_back(0);
                        }
                    }
                    
                    // have we already translated the membership records for the path on this node?
                    // (we need to check this to avoid falsely translating pointers that we have actually
                    // already translated)
                    if (id_translated.get(step_node_id - min_translated_id) != 1) {
                        
                        size_t member_idx = path_membership_node_iv.get(graph_index_to_node_member_index(graph_iv_index(handle)));
                        while (member_idx) {
                            
                            // update the offsets for membership records on this path
                            if (get_membership_path(member_idx) == path_id_here) {
                                set_membership_step(member_idx, offset_translator.get(get_membership_step(member_idx)));
                            }
                            
                            // move to the next membership record
                            member_idx = get_next_membership(member_idx);
                        }
                        
                        // mark this node as updated so we don't re-update the offsets
                        id_translated.set(step_node_id - min_translated_id, 1);
                    }
                    
                    first_iter = false;
                }
            }
            else {
                // the path is empty, so let's make sure it's not holding onto any capacity it doesn't need
                path.links_iv = PagedVector(PAGE_WIDTH);
                path.steps_iv = PagedVector(PAGE_WIDTH);
            }
            
            path.deleted_step_records = 0;
        }
    }
    
    void PackedGraph::eject_deleted_paths() {
        
        uint64_t num_paths_deleted = 0;
        vector<uint64_t> paths_deleted_before(paths.size(), 0);
        for (size_t i = 0; i < paths.size(); i++) {
            
            paths_deleted_before[i] = num_paths_deleted;
            
            if (paths[i].is_deleted) {
                num_paths_deleted++;
                continue;
            }
            
            // move non-deleted paths into the front of the vector
            if (num_paths_deleted > 0) {
                paths[i - num_paths_deleted] = std::move(paths[i]);
            }
        }
        
        // eliminate the empty spots we left at the end of the vector
        if (num_paths_deleted > 0) {
            for (size_t i = 0; i < num_paths_deleted; ++i) {
                // we use pop back instead of resize because it doesn't require a default constructor
                paths.pop_back();
            }
            paths.shrink_to_fit();
            
            // update the path IDs
            for (size_t i = 0; i < paths.size(); ++i) {
                path_id[paths[i].name] = i;
            }
            
            // update the path IDs in the membership records
            for (size_t i = 0; i < path_membership_id_iv.size(); i += MEMBERSHIP_ID_RECORD_SIZE) {
                uint64_t current_path = path_membership_id_iv.get(i);
                path_membership_id_iv.set(i, current_path - paths_deleted_before[current_path]);
            }
        }
    }
    
    void PackedGraph::compact_ids() {
        
        // use an overlay to convert to a single stranded digraph
        StrandSplitGraph digraph(this);
        
        // get a low FAS layout using Eades-Lin-Smyth algorithm
        vector<handle_t> layout = algorithms::eades_algorithm(&digraph);
        // note: the single stranded graph will have a fully separated forward and reverse strands, so
        // we have the guarantee that every handle in this layout is forward in the strand split graph
        
        // in place, take only the handles that are forward in the source graph and convert them back
        // to the source handles
        size_t skipped = 0;
        for (size_t i = 0; i < layout.size(); ++i) {
            handle_t underlying = digraph.get_underlying_handle(layout[i]);
            if (get_is_reverse(underlying)) {
                ++skipped;
            }
            else {
                layout[i - skipped] = underlying;
            }
        }
        // remove everything we skipped
        layout.resize(layout.size() - skipped);
        
        assert(layout.size() == node_size());
        
        // use the layout to make a translator between current IDs and the IDs we will reassign
        PagedVector id_trans(PAGE_WIDTH);
        id_trans.resize(max_id - min_id + 1);
        for (size_t i = 0; i < layout.size(); ++i) {
            id_trans.set(get_id(layout[i]) - min_id, i + 1);
        }
        
        // we don't need the layout anymore, so release the memory
        layout.clear();
        
        // update the node IDs of edges
        for (size_t i = EDGE_TRAV_OFFSET; i < edge_lists_iv.size(); i += EDGE_RECORD_SIZE) {
            handle_t trav = decode_traversal(edge_lists_iv.get(i));
            // only translate edges to nodes that have not been deleted
            if (id_to_graph_iv.get(get_id(trav) - min_id)) {
                trav = get_handle(id_trans.get(get_id(trav) - min_id), get_is_reverse(trav));
                edge_lists_iv.set(i, encode_traversal(trav));
            }
        }
        
        // update the node IDs of steps on paths
        for (PackedPath& packed_path : paths){
            
            if (packed_path.is_deleted) {
                continue;
            }
            
            for (size_t i = 0; i < packed_path.steps_iv.size(); i += STEP_RECORD_SIZE) {
                handle_t trav = decode_traversal(packed_path.steps_iv.get(i));
                // only translate step records of nodes that have not been deleted
                if (id_to_graph_iv.get(get_id(trav) - min_id)) {
                    trav = get_handle(id_trans.get(get_id(trav) - min_id), get_is_reverse(trav));
                    packed_path.steps_iv.set(i, encode_traversal(trav));
                }
            }
        }
        
        // make a vector to translate the new IDs to the offset of the node
        PackedDeque new_id_to_graph_iv;
        new_id_to_graph_iv.reserve(node_size());
        for (id_t node_id = min_id; node_id <= max_id; ++node_id) {
            size_t offset = id_to_graph_iv.get(node_id - min_id);
            if (offset) {
                new_id_to_graph_iv.append_back(offset);
            }
        }
        
        id_to_graph_iv = move(new_id_to_graph_iv);
        
        if (new_id_to_graph_iv.size() > 0) {
            min_id = 1;
            max_id = id_to_graph_iv.size();
        }
        else {
            min_id = numeric_limits<id_t>::max();
            max_id = 0;
        }
    }
    
    void PackedGraph::tighten() {
        
        // remove deleted paths and force them to eject deleted material
        for (size_t i = 0; i < paths.size(); i++) {
            // force the path to defragment
            defragment_path(paths[i], true);
        }
        
        // push any paths we deleted out of the path vector
        eject_deleted_paths();
        
        // force the graph structures to reallocate in ID order and eject deleted material
        defragment(true);
        
        // make a new id_to_graph_iv of exactly the right size
        PackedDeque new_id_to_graph_iv;
        new_id_to_graph_iv.reserve(id_to_graph_iv.size());
        // transfer of the data
        for (size_t i = 0; i < id_to_graph_iv.size(); i++) {
            new_id_to_graph_iv.append_back(id_to_graph_iv.get(i));
        }
        // replace the old one
        id_to_graph_iv = std::move(new_id_to_graph_iv);
        
        // count up the total length of all non-deleted sequence (a little costly but we're
        // not tracking deleted sequence anywhere)
        size_t total_seq_len = 0;
        for (size_t i = 0; i < seq_length_iv.size(); i += SEQ_LENGTH_RECORD_SIZE) {
            total_seq_len += seq_length_iv.get(i);
        }
        
        // make a new seq_iv of exactly the right size
        PackedVector new_seq_iv;
        new_seq_iv.reserve(total_seq_len);
        static_assert(SEQ_START_RECORD_SIZE == SEQ_LENGTH_RECORD_SIZE,
                      "This loop will need to be rewritten if we change the record sizes");
        for (size_t i = 0; i < seq_start_iv.size(); i += SEQ_START_RECORD_SIZE) {
            // get the interval from the current seq_iv
            size_t begin = seq_start_iv.get(i);
            size_t end = begin + seq_length_iv.get(i);
            // switch the pointer to the new seq iv
            seq_start_iv.set(i, new_seq_iv.size());
            // transfer the actual sequence over
            for (size_t j = begin; j < end; j++) {
                new_seq_iv.append(seq_iv.get(j));
            }
        }
        // replace the old seq iv
        seq_iv = std::move(new_seq_iv);
    }
    
    void PackedGraph::defragment(bool force) {
        
        uint64_t num_nodes = graph_iv.size() / GRAPH_RECORD_SIZE - deleted_node_records;
        if (deleted_node_records > defrag_factor * (graph_iv.size() / GRAPH_RECORD_SIZE) || force) {
            
            // what's the real number of undeleted nodes in the graph?
            uint64_t num_nodes = graph_iv.size() / GRAPH_RECORD_SIZE - deleted_node_records;
            
            // adjust the start
            while (id_to_graph_iv.empty() ? false : id_to_graph_iv.get(0) == 0) {
                id_to_graph_iv.pop_front();
                min_id++;
            }
            // adjust the end
            while (id_to_graph_iv.empty() ? false : id_to_graph_iv.get(id_to_graph_iv.size() - 1) == 0) {
                id_to_graph_iv.pop_back();
            }
            max_id = min_id + id_to_graph_iv.size() - 1;
            
            // initialize new vectors to construct defragged copies in
            PagedVector new_graph_iv(PAGE_WIDTH);
            PackedVector new_seq_length_iv;
            PagedVector new_seq_start_iv(PAGE_WIDTH);
            PagedVector new_path_membership_node_iv(PAGE_WIDTH);
            
            // expand them to the size we need to avoid reallocation and get optimal compression
            new_graph_iv.reserve(num_nodes * GRAPH_RECORD_SIZE);
            new_seq_length_iv.reserve(num_nodes * SEQ_LENGTH_RECORD_SIZE);
            new_seq_start_iv.reserve(num_nodes * SEQ_START_RECORD_SIZE);
            new_path_membership_node_iv.reserve(num_nodes * NODE_MEMBER_RECORD_SIZE);
            
            for (size_t i = 0; i < id_to_graph_iv.size(); i++) {
                size_t raw_g_iv_idx = id_to_graph_iv.get(i);
                if (raw_g_iv_idx) {
                    size_t g_iv_idx = (raw_g_iv_idx - 1) * GRAPH_RECORD_SIZE;
                    // this node still exists, create a new copy
                    new_graph_iv.append(graph_iv.get(g_iv_idx + GRAPH_START_EDGES_OFFSET));
                    new_graph_iv.append(graph_iv.get(g_iv_idx + GRAPH_END_EDGES_OFFSET));
                    new_seq_length_iv.append(seq_length_iv.get(graph_index_to_seq_len_index(g_iv_idx)));
                    new_seq_start_iv.append(seq_start_iv.get(graph_index_to_seq_start_index(g_iv_idx)));
                    new_path_membership_node_iv.append(path_membership_node_iv.get(graph_index_to_node_member_index(g_iv_idx)));
                    // update the pointer into graph_iv
                    id_to_graph_iv.set(i, new_graph_iv.size() / GRAPH_RECORD_SIZE);
                }
            }
            
            // replace graph with the defragged copy
            graph_iv = std::move(new_graph_iv);
            seq_length_iv = std::move(new_seq_length_iv);
            seq_start_iv = std::move(new_seq_start_iv);
            path_membership_node_iv = std::move(new_path_membership_node_iv);
            
            deleted_node_records = 0;
        }
        
        // TODO: also defragment seq_iv?
        // for now only doing that inside tighten()
        
        if (deleted_edge_records > defrag_factor * (edge_lists_iv.size() / EDGE_RECORD_SIZE) || force) {
            
            uint64_t num_edge_records = edge_lists_iv.size() / EDGE_RECORD_SIZE - deleted_edge_records;
            
            PagedVector new_edge_lists_iv(PAGE_WIDTH);
            new_edge_lists_iv.reserve(num_edge_records * EDGE_RECORD_SIZE);
            
            for (size_t i = 0; i < id_to_graph_iv.size(); i++) {
                size_t raw_g_iv_idx = id_to_graph_iv.get(i);
                if (raw_g_iv_idx) {
                    // this node still exists
                    size_t g_iv_idx = (raw_g_iv_idx - 1) * GRAPH_RECORD_SIZE;
                    for (size_t edge_list_offset : {GRAPH_START_EDGES_OFFSET, GRAPH_END_EDGES_OFFSET}) {
                        size_t edge_list_idx = graph_iv.get(g_iv_idx + edge_list_offset);
                        if (edge_list_idx) {
                            // add a new edge record
                            new_edge_lists_iv.append(get_edge_target(edge_list_idx));
                            new_edge_lists_iv.append(0);
                            // point the graph vector at this new edge list
                            graph_iv.set(g_iv_idx + edge_list_offset, new_edge_lists_iv.size() / EDGE_RECORD_SIZE);
                            
                            edge_list_idx = get_next_edge_index(edge_list_idx);
                            while (edge_list_idx) {
                                // add a new edge record
                                new_edge_lists_iv.append(get_edge_target(edge_list_idx));
                                new_edge_lists_iv.append(0);
                                // point the previous link at this one
                                new_edge_lists_iv.set(new_edge_lists_iv.size() - 2 * EDGE_RECORD_SIZE + EDGE_NEXT_OFFSET,
                                                      new_edge_lists_iv.size() / EDGE_RECORD_SIZE);
                                
                                edge_list_idx = get_next_edge_index(edge_list_idx);
                            }
                        }
                    }
                }
            }
            
            edge_lists_iv = std::move(new_edge_lists_iv);
            
            deleted_edge_records = 0;
        }
        
        if (deleted_membership_records > defrag_factor * (path_membership_next_iv.size() / MEMBERSHIP_NEXT_RECORD_SIZE) || force) {
            
            uint64_t num_membership_records = path_membership_next_iv.size() / MEMBERSHIP_NEXT_RECORD_SIZE - deleted_membership_records;
            
            PagedVector new_path_membership_id_iv(PAGE_WIDTH);
            PagedVector new_path_membership_offset_iv(PAGE_WIDTH);
            PagedVector new_path_membership_next_iv(PAGE_WIDTH);
            
            new_path_membership_id_iv.reserve(num_membership_records * MEMBERSHIP_ID_RECORD_SIZE);
            new_path_membership_offset_iv.reserve(num_membership_records * MEMBERSHIP_OFFSET_RECORD_SIZE);
            new_path_membership_next_iv.reserve(num_membership_records * MEMBERSHIP_NEXT_RECORD_SIZE);
            
            for (size_t i = 0; i < id_to_graph_iv.size(); i++) {
                size_t raw_g_iv_idx = id_to_graph_iv.get(i);
                if (raw_g_iv_idx) {
                    // this node still exists
                    size_t g_iv_idx = (raw_g_iv_idx - 1) * GRAPH_RECORD_SIZE;
                    
                    uint64_t member_idx = path_membership_node_iv.get(graph_index_to_node_member_index(g_iv_idx));
                    if (member_idx) {
                        // make a new membership record
                        new_path_membership_id_iv.append(get_membership_path(member_idx));
                        new_path_membership_offset_iv.append(get_membership_step(member_idx));
                        new_path_membership_next_iv.append(0);
                        
                        // point the membership vector here
                        path_membership_node_iv.set(graph_index_to_node_member_index(g_iv_idx),
                                                    new_path_membership_next_iv.size() / MEMBERSHIP_NEXT_RECORD_SIZE);
                        
                        member_idx = get_next_membership(member_idx);
                        while (member_idx) {
                            // make a new membership record
                            new_path_membership_id_iv.append(get_membership_path(member_idx));
                            new_path_membership_offset_iv.append(get_membership_step(member_idx));
                            new_path_membership_next_iv.append(0);
                            // point the previous link at this one
                            new_path_membership_next_iv.set(new_path_membership_next_iv.size() - 2 * MEMBERSHIP_NEXT_RECORD_SIZE,
                                                            new_path_membership_next_iv.size() / MEMBERSHIP_NEXT_RECORD_SIZE);
                            
                            member_idx = get_next_membership(member_idx);
                        }
                    }
                }
            }
            
            path_membership_id_iv = std::move(new_path_membership_id_iv);
            path_membership_offset_iv = std::move(new_path_membership_offset_iv);
            path_membership_next_iv = std::move(new_path_membership_next_iv);
            
            deleted_membership_records = 0;
        }
    }
    
    void PackedGraph::optimize(bool allow_id_reassignment) {
        
        if (allow_id_reassignment) {
            // reassign IDs into a contiguous interval ordered by an approximate sort
            compact_ids();
        }
        
        // tighten up vector allocations and straighten out the linked lists they contain
        tighten();
    }
    
    void PackedGraph::clear(void) {
        graph_iv.clear();
        seq_start_iv.clear();
        seq_length_iv.clear();
        edge_lists_iv.clear();
        id_to_graph_iv.clear();
        seq_iv.clear();
        path_membership_node_iv.clear();
        path_membership_id_iv.clear();
        path_membership_offset_iv.clear();
        path_membership_next_iv.clear();
        paths.clear();
        path_id.clear();
        min_id = std::numeric_limits<id_t>::max();
        max_id = 0;
        deleted_edge_records = 0;
        deleted_node_records = 0;
        deleted_membership_records = 0;
    }
    
    bool PackedGraph::has_path(const std::string& path_name) const {
        return path_id.count(path_name);
    }
    
    path_handle_t PackedGraph::get_path_handle(const std::string& path_name) const {
        return as_path_handle(path_id.at(path_name));
    }
    
    string PackedGraph::get_path_name(const path_handle_t& path_handle) const {
        return paths.at(as_integer(path_handle)).name;
    }
    
    bool PackedGraph::get_is_circular(const path_handle_t& path_handle) const {
        return paths.at(as_integer(path_handle)).is_circular;
    }
    
    size_t PackedGraph::get_step_count(const path_handle_t& path_handle) const {
        // TODO: if we every allow step deletes, this will need to be adjusted
        const PackedPath& path = paths.at(as_integer(path_handle));
        return path.steps_iv.size() / STEP_RECORD_SIZE - path.deleted_step_records;
    }
    
    size_t PackedGraph::get_path_count() const {
        return path_id.size();
    }
    
    bool PackedGraph::for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const {
        for (const auto& path_id_record : path_id) {
            if (!iteratee(as_path_handle(path_id_record.second))) {
                return false;
            }
        }
        return true;
    }
    
    handle_t PackedGraph::get_handle_of_step(const step_handle_t& step_handle) const {
        const PackedPath& path = paths.at(as_integers(step_handle)[0]);
        return decode_traversal(get_step_trav(path, as_integers(step_handle)[1]));
    }
    
    step_handle_t PackedGraph::path_begin(const path_handle_t& path_handle) const {
        step_handle_t step;
        as_integers(step)[0] = as_integer(path_handle);
        as_integers(step)[1] = paths.at(as_integer(path_handle)).head;
        return step;
    }
    
    step_handle_t PackedGraph::path_end(const path_handle_t& path_handle) const {
        step_handle_t step;
        as_integers(step)[0] = as_integer(path_handle);
        as_integers(step)[1] = 0;
        return step;
    }
    
    step_handle_t PackedGraph::get_next_step(const step_handle_t& step_handle) const {
        step_handle_t next;
        as_integers(next)[0] = as_integers(step_handle)[0];
        as_integers(next)[1] = get_step_next(paths.at(as_integers(step_handle)[0]),
                                             as_integers(step_handle)[1]);
        return next;
    }
    
    step_handle_t PackedGraph::get_previous_step(const step_handle_t& step_handle) const {
        step_handle_t prev;
        as_integers(prev)[0] = as_integers(step_handle)[0];
        as_integers(prev)[1] = as_integers(step_handle)[1] != 0 ? get_step_prev(paths.at(as_integers(step_handle)[0]),
                                                                                as_integers(step_handle)[1])
                                                                : paths.at(as_integers(step_handle)[0]).tail;
        return prev;
    }
    
    path_handle_t PackedGraph::get_path_handle_of_step(const step_handle_t& step_handle) const {
        return as_path_handle(as_integers(step_handle)[0]);
    }
    
    
    bool PackedGraph::for_each_step_on_handle_impl(const handle_t& handle,
                                                   const function<bool(const step_handle_t&)>& iteratee) const {
        
        size_t path_membership = path_membership_node_iv.get(graph_index_to_node_member_index(graph_iv_index(handle)));
        while (path_membership) {
            
            // get the path that this membership record is on
            uint64_t path_id = get_membership_path(path_membership);
            const PackedPath& packed_path = paths[path_id];
            
            // get the traversal
            size_t occ_idx = get_membership_step(path_membership);
            handle_t trav = decode_traversal(get_step_trav(packed_path, occ_idx));
            
            // send along this step
            step_handle_t step_handle;
            as_integers(step_handle)[0] = path_id;
            as_integers(step_handle)[1] = occ_idx;
            if (!iteratee(step_handle)) {
                return false;
            }
            
            // move to the next membership record
            path_membership = get_next_membership(path_membership);
        }
        
        return true;
    }
    
    void PackedGraph::destroy_path(const path_handle_t& path) {
        
        PackedPath& packed_path = paths.at(as_integer(path));
        
        // remove node membership records corresponding to this path
        bool first_iter = true;
        for (uint64_t step_offset = packed_path.head;
             step_offset != 0 && (step_offset != packed_path.head || first_iter);
             step_offset = get_step_next(packed_path, step_offset)) {
            
            uint64_t trav = get_step_trav(packed_path, step_offset);
            size_t node_member_idx = graph_index_to_node_member_index(graph_iv_index(decode_traversal(trav)));
            
            // find a membership record for this path
            size_t prev = 0;
            size_t here = path_membership_node_iv.get(node_member_idx);
            while (as_path_handle(get_membership_path(here)) != path) {
                prev = here;
                here = get_next_membership(here);
                // note: we don't need to be careful about getting the exact corresponding step since this node
                // should try to delete a membership record exactly as many times as it occurs on this path -- all of
                // the records will get deleted
            }
            
            if (prev == 0) {
                // this was the first record, set following one to be the head
                path_membership_node_iv.set(node_member_idx, get_next_membership(here));
            }
            else {
                // make the link from the previous record skip over the current one
                set_next_membership(prev, get_next_membership(here));
            }
            
            ++deleted_membership_records;
            
            first_iter = false;
        }
        
        path_id.erase(packed_path.name);
        
        packed_path.is_deleted = true;
        packed_path.steps_iv.clear();
        packed_path.links_iv.clear();
        packed_path.name.clear();
        packed_path.head = 0;
        packed_path.tail = 0;
        
        defragment();
    }
    
    path_handle_t PackedGraph::create_path_handle(const string& name, bool is_circular) {
        path_id[name] = paths.size();
        path_handle_t path_handle = as_path_handle(paths.size());
        paths.emplace_back(name, is_circular);
        return path_handle;
    }
    
    step_handle_t PackedGraph::append_step(const path_handle_t& path, const handle_t& to_append) {
        
        PackedPath& packed_path = paths.at(as_integer(path));
        
        // create a new path record
        packed_path.steps_iv.append(as_integer(to_append));
        packed_path.links_iv.append(packed_path.tail);
        packed_path.links_iv.append(0);
        
        // the offset associated with the new record
        size_t step_offset = packed_path.steps_iv.size() / STEP_RECORD_SIZE;
        
        // update the pointer from the current tail
        if (packed_path.tail != 0) {
            set_step_next(packed_path, packed_path.tail, step_offset);
        }
        
        // update the head and tail of the list
        packed_path.tail = step_offset;
        if (packed_path.head == 0) {
            packed_path.head = step_offset;
        }
        
        // update the looping connection if this is a circular path
        if (packed_path.is_circular) {
            set_step_prev(packed_path, packed_path.head, packed_path.tail);
            set_step_next(packed_path, packed_path.tail, packed_path.head);
        }
        
        // record the membership of this node in this path
        size_t node_member_idx = graph_index_to_node_member_index(graph_iv_index(to_append));
        path_membership_id_iv.append(as_integer(path));
        path_membership_offset_iv.append(step_offset);
        path_membership_next_iv.append(path_membership_node_iv.get(node_member_idx));
        
        // make this new membership record the head of the linked list
        path_membership_node_iv.set(node_member_idx, path_membership_next_iv.size() / MEMBERSHIP_NEXT_RECORD_SIZE);
        
        // make and return an step handle
        step_handle_t step;
        as_integers(step)[0] = as_integer(path);
        as_integers(step)[1] = step_offset;
        return step;
    }
    
    step_handle_t PackedGraph::prepend_step(const path_handle_t& path, const handle_t& to_prepend) {
        
        PackedPath& packed_path = paths.at(as_integer(path));
        
        // create a new path record
        packed_path.steps_iv.append(as_integer(to_prepend));
        packed_path.links_iv.append(0);
        packed_path.links_iv.append(packed_path.head);
        
        // the offset associated with the new record
        size_t step_offset = packed_path.steps_iv.size() / STEP_RECORD_SIZE;
        
        // update the pointer from the current head
        if (packed_path.head != 0) {
            set_step_prev(packed_path, packed_path.head, step_offset);
        }
        
        // update the head and tail of the list
        packed_path.head = step_offset;
        if (packed_path.tail == 0) {
            packed_path.tail = step_offset;
        }
        
        // update the looping connection if this is a circular path
        if (packed_path.is_circular) {
            set_step_prev(packed_path, packed_path.head, packed_path.tail);
            set_step_next(packed_path, packed_path.tail, packed_path.head);
        }
        
        // record the membership of this node in this path
        size_t node_member_idx = graph_index_to_node_member_index(graph_iv_index(to_prepend));
        path_membership_id_iv.append(as_integer(path));
        path_membership_offset_iv.append(step_offset);
        path_membership_next_iv.append(path_membership_node_iv.get(node_member_idx));
        
        // make this new membership record the head of the linked list
        path_membership_node_iv.set(node_member_idx, path_membership_next_iv.size() / MEMBERSHIP_NEXT_RECORD_SIZE);
        
        // make and return an step handle
        step_handle_t step;
        as_integers(step)[0] = as_integer(path);
        as_integers(step)[1] = step_offset;
        return step;
    }
    
    pair<step_handle_t, step_handle_t> PackedGraph::rewrite_segment(const step_handle_t& segment_begin,
                                                                    const step_handle_t& segment_end,
                                                                    const vector<handle_t>& new_segment) {
        
        if (get_path_handle_of_step(segment_begin) != get_path_handle_of_step(segment_end)) {
            cerr << "error:[PackedGraph] attempted to rewrite a path segment delimited by steps on two different paths" << endl;
            exit(1);
        }
        
        PackedPath& packed_path =  paths.at(as_integers(segment_begin)[0]);
        
        // TODO: somewhat repetitive with a routine in destroy_path
        
        // find and erase the record of this node's membership on the path
        for (size_t step_offset = as_integers(segment_begin)[1];
             step_offset != as_integers(segment_end)[1]; ) {
            
            size_t g_iv_idx = graph_iv_index(decode_traversal(get_step_trav(packed_path, step_offset)));
            size_t node_member_idx = graph_index_to_node_member_index(g_iv_idx);
            
            // find the membership record that corresponds to this step
            size_t prev = 0;
            size_t here = path_membership_node_iv.get(node_member_idx);
            while (get_membership_path(here) != as_integers(segment_end)[0] ||
                   get_membership_step(here) != step_offset) {
                prev = here;
                here = get_next_membership(here);
            }
            
            if (prev == 0) {
                // this was the first record, set following one to be the head
                path_membership_node_iv.set(node_member_idx, get_next_membership(here));
            }
            else {
                // make the link from the previous record skip over the current one
                set_next_membership(prev, get_next_membership(here));
            }
            
            ++deleted_membership_records;
            
            // get the adjacent nodes in the path
            size_t prev_offset = get_step_prev(packed_path, step_offset);
            size_t next_offset = get_step_next(packed_path,step_offset);
            
            // make their links skip over this node
            if (prev_offset != 0) {
                set_step_next(packed_path, prev_offset, next_offset);
            }
            if (next_offset != 0) {
                set_step_prev(packed_path, next_offset, prev_offset);
            }
            
            // update the head and tail of the path if necessary
            if (step_offset == packed_path.head && step_offset == packed_path.tail) {
                // this is the last node in the path, set the head and tail to null
                packed_path.head = packed_path.tail = 0;
            }
            else if (step_offset == packed_path.head) {
                packed_path.head = next_offset;
            }
            else if (step_offset == packed_path.tail) {
                packed_path.tail = prev_offset;
            }
            
            packed_path.deleted_step_records++;
            
            // TODO: reallocating paths invalidates pointers, so we can't really do it here because
            // we need to return the range. might also be confusing to users
            
            // maybe reallocate to address fragmentation within the path
            //defragment_path(packed_path);
            
            step_offset = next_offset;
        }
        
        pair<step_handle_t, step_handle_t> new_segment_range(segment_end, segment_end);
        
        // now add in the new segment
        bool first_iter = true;
        uint64_t anchor_offset = as_integers(segment_end)[1];
        for (const handle_t& handle : new_segment) {
            
            // create a new step record
            packed_path.steps_iv.append(encode_traversal(handle));
            packed_path.links_iv.append(0);
            packed_path.links_iv.append(0);
            
            uint64_t step_offset = packed_path.steps_iv.size() / STEP_RECORD_SIZE;
            
            if (anchor_offset != 0) {
                
                // insert before
                uint64_t anchor_prev = get_step_prev(packed_path, anchor_offset);
                set_step_prev(packed_path, step_offset, get_step_prev(packed_path, anchor_offset));
                if (anchor_prev != 0) {
                    set_step_next(packed_path, anchor_prev, step_offset);
                }
                
                set_step_next(packed_path, step_offset, anchor_offset);
                set_step_prev(packed_path, anchor_offset, step_offset);
            }
            else {
                // place after the tail, since we're not putting it before anything
                if (packed_path.tail != 0) {
                    // attach to the tail
                    uint64_t tail_next = get_step_next(packed_path, packed_path.tail);
                    set_step_next(packed_path, packed_path.tail, step_offset);
                    set_step_prev(packed_path, step_offset, packed_path.tail);
                    
                    // handle the potential looping connection
                    if (tail_next != 0) {
                        set_step_next(packed_path, step_offset, tail_next);
                        set_step_prev(packed_path, tail_next, step_offset);
                    }
                }
                
                // we're placing at the end, so this is the new tail
                packed_path.tail = step_offset;
            }
            
            if (anchor_offset == packed_path.head) {
                // we're placing before the head (or the head is null), so this is new head
                packed_path.head = step_offset;
            }
            
            // put a membership record for this occurrence at the front of the membership list
            size_t node_member_idx = graph_index_to_node_member_index(graph_iv_index(handle));
            path_membership_id_iv.append(as_integers(segment_end)[0]);
            path_membership_offset_iv.append(step_offset);
            path_membership_next_iv.append(path_membership_node_iv.get(node_member_idx));
            path_membership_node_iv.set(node_member_idx, path_membership_next_iv.size() / MEMBERSHIP_NEXT_RECORD_SIZE);
            
            if (first_iter) {
                // record the start of the new range
                as_integers(new_segment_range.first)[1] = step_offset;
                first_iter = false;
            }
        }
        
        return new_segment_range;
    }
    
    void PackedGraph::set_circularity(const path_handle_t& path, bool circular) {
        PackedPath& packed_path = paths[as_integer(path)];
        // set the looping connection as appropriate
        if (circular && packed_path.head != 0) {
            set_step_prev(packed_path, packed_path.head, packed_path.tail);
            set_step_next(packed_path, packed_path.tail, packed_path.head);
        }
        else if (!circular && packed_path.head != 0) {
            set_step_prev(packed_path, packed_path.head, 0);
            set_step_next(packed_path, packed_path.tail, 0);
        }
        // set the annotation
        packed_path.is_circular = circular;
    }
}

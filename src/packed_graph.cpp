//
//  dgraph.cpp
//

#include "packed_graph.hpp"

namespace vg {
    
    const double PackedGraph::defrag_factor = .2;
    
    PackedGraph::PackedGraph() : graph_iv(PAGE_WIDTH), seq_start_iv(PAGE_WIDTH), edge_lists_iv(PAGE_WIDTH) {
        
    }
    
    PackedGraph::~PackedGraph() {
        
    }
    
    
    size_t PackedGraph::new_node_record(id_t node_id) {
        
        size_t next_g_iv_idx = graph_iv.size();
        
        // no edges yet, null pointer for linked list
        graph_iv.append(0);
        graph_iv.append(0);
        
        // record the sequence interval
        seq_start_iv.append(0);
        seq_length_iv.append(0);
        
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
    
    handle_t PackedGraph::create_handle(const std::string& sequence) {
        return create_handle(sequence, max_id + 1);
    }
    
    handle_t PackedGraph::create_handle(const std::string& sequence, const id_t& id) {
        
        // TODO: don't duplicate an existing node
        
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
        edge_lists_iv.append(encode_edge_target(right));
        edge_lists_iv.append(graph_iv.get(g_iv_left));
        // make this new node the head
        graph_iv.set(g_iv_left, edge_lists_iv.size() / EDGE_RECORD_SIZE);
        
        // don't double add a reversing self edge
        if (g_iv_left == g_iv_right) {
            return;
        }
        
        // add a new linked list node pointing to the rest of the list
        edge_lists_iv.append(encode_edge_target(flip(left)));
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
        return EasyHandlePacking::pack(node_id, is_reverse);
    }
    
    id_t PackedGraph::get_id(const handle_t& handle) const {
        return EasyHandlePacking::unpack_number(handle);
    }
    
    bool PackedGraph::get_is_reverse(const handle_t& handle) const {
        return EasyHandlePacking::unpack_bit(handle);;
    }
    
    handle_t PackedGraph::flip(const handle_t& handle) const {
        return EasyHandlePacking::toggle_bit(handle);
    }
    
    size_t PackedGraph::get_length(const handle_t& handle) const {
        return seq_length_iv.get(graph_index_to_seq_len_index(graph_iv_index(handle)));
    }
    
    std::string PackedGraph::get_sequence(const handle_t& handle) const {
        size_t g_iv_index = graph_iv_index(handle);
        size_t seq_start = seq_start_iv.get(graph_index_to_seq_start_index(g_iv_index));
        size_t seq_len = seq_length_iv.get(graph_index_to_seq_len_index(g_iv_index));
        std::string seq(seq_len, 'N');
        for (size_t i = 0; i < seq_len; i++) {
            seq[i] = decode_nucleotide(seq_iv.get(seq_start + i));
        }
        return get_is_reverse(handle) ? reverse_complement(seq) : seq;
    }
    
    void PackedGraph::swap_handles(const handle_t& a, const handle_t& b) {
        // TODO: this doesn't actually affect the traversal order
        // We either need to encode the IDs here or use different method to keep track
        // of deleted nodes
        
        size_t g_iv_index_a = graph_iv_index(a);
        size_t g_iv_index_b = graph_iv_index(b);
        
        int64_t id_a = get_id(a);
        int64_t id_b = get_id(b);
        
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
        
        id_to_graph_iv.set(id_a - min_id, g_iv_index_b / GRAPH_RECORD_SIZE + 1);
        id_to_graph_iv.set(id_b - min_id, g_iv_index_a / GRAPH_RECORD_SIZE + 1);
    }
    
    bool PackedGraph::follow_edges(const handle_t& handle, bool go_left,
                                                    const std::function<bool(const handle_t&)>& iteratee) const {
        // toward start = true, toward end = false
        bool direction = get_is_reverse(handle) != go_left;
        // get the head of the linked list from the graph vector
        size_t edge_idx = graph_iv.get(graph_iv_index(handle)
                                       + (direction ? GRAPH_START_EDGES_OFFSET : GRAPH_END_EDGES_OFFSET));
        // traverse the linked list as long as directed
        bool keep_going = true;
        while (edge_idx && keep_going) {
            
            handle_t edge_target = decode_edge_target(get_edge_target(edge_idx));
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
        return graph_iv.size() / GRAPH_RECORD_SIZE;
    }
    
    id_t PackedGraph::min_node_id(void) const {
        return min_id;
    }
    
    id_t PackedGraph::max_node_id(void) const {
        return max_id;
    }
    
    void PackedGraph::for_each_handle(const std::function<bool(const handle_t&)>& iteratee,
                                                       bool parallel) const {
        
        size_t num_id_positions = id_to_graph_iv.size();
        if (parallel) {
            // TODO: add OMP pragma back in
            // TODO: would task based parallelism be better?
//#pragma omp parallel for
            for (size_t i = 0; i < num_id_positions; i++) {
                if (id_to_graph_iv.get(i)) {
                    iteratee(get_handle(i + min_id));
                }
            }
        }
        else {
            for (size_t i = 0; i < num_id_positions; i++) {
                if (id_to_graph_iv.get(i)) {
                    iteratee(get_handle(i + min_id));
                }
            }
        }
        
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
                // iterate down the entire edge list                                 GRAPH_START_EDGES_OFFSET));
                handle_t looking_for = orientation ? handle : flip(handle);
                size_t edge_list_idx = graph_iv.get(g_iv_idx + (orientation ? GRAPH_START_EDGES_OFFSET : GRAPH_END_EDGES_OFFSET));
                while (edge_list_idx) {
                    handle_t target = decode_edge_target(get_edge_target(edge_list_idx));
                    size_t backward_edge_idx = graph_iv.get(graph_iv_index(target) + (get_is_reverse(target) ?
                                                                                      GRAPH_END_EDGES_OFFSET :
                                                                                      GRAPH_START_EDGES_OFFSET));
                    
                    while (backward_edge_idx) {
                        handle_t backward_edge_target = decode_edge_target(get_edge_target(backward_edge_idx));
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
            
            for (size_t i = 0; i < seq_len / 2; i++) {
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
        std::vector<size_t> forward_offsets = offsets;
        size_t node_length = get_length(handle);
        if (get_is_reverse(handle)) {
            for (size_t& off : forward_offsets) {
                off = node_length - off;
            }
        }
        
        // we will also build the return value in forward orientation
        handle_t forward_handle = get_is_reverse(handle) ? flip(handle) : handle;
        std::vector<handle_t> return_val{forward_handle};
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
            edge_lists_iv.append(encode_edge_target(get_handle(next_id)));
            edge_lists_iv.append(0);
            // add the edge onto the previous node
            graph_iv.set(g_iv_idx + GRAPH_END_EDGES_OFFSET, edge_lists_iv.size() / EDGE_RECORD_SIZE);
            
            // create an edge backward to the previous node
            edge_lists_iv.append(encode_edge_target(get_handle(prev_id, true)));
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
            handle_t target = decode_edge_target(get_edge_target(next_edge_idx));
            size_t backward_edge_idx = graph_iv.get(graph_iv_index(target) + (get_is_reverse(target) ?
                                                                              GRAPH_END_EDGES_OFFSET :
                                                                              GRAPH_START_EDGES_OFFSET));
            
            while (backward_edge_idx) {
                handle_t backward_edge_target = decode_edge_target(get_edge_target(backward_edge_idx));
                if (backward_edge_target == looking_for) {
                    set_edge_target(backward_edge_idx, flip(return_val.back()));
                    break;
                }
                backward_edge_idx = get_next_edge_index(backward_edge_idx);
            }
            
            next_edge_idx = get_next_edge_index(next_edge_idx);
        }
        
        if (get_is_reverse(handle)) {
            // reverse the vector to the orientation of the input handle
            std::reverse(return_val.begin(), return_val.end());
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
            // the edge on the deleting node to be deleted
            deleted_edge_records++;
            return true;
        });
        follow_edges(handle, true, [&](const handle_t& prev) {
            remove_edge_reference(prev, handle);
            // we don't actually bother removing the reference, but we will also consider
            // the edge on the deleting node to be deleted
            deleted_edge_records++;
            return true;
        });
        
        // remove the reference to the node
        id_to_graph_iv.set(get_id(handle) - min_id, 0);
        
        deleted_node_records++;
        
        // maybe reallocate to address fragmentation
        defragment();
    }
    
    void PackedGraph::remove_edge_reference(const handle_t& on, const handle_t& to) {
        
        // Note: this function assumes that edge ref exists, crashes otherwise
        
        size_t g_iv_idx = graph_iv_index(on) + (get_is_reverse(on)
                                                ? GRAPH_START_EDGES_OFFSET
                                                : GRAPH_END_EDGES_OFFSET);

        size_t edge_list_idx = graph_iv.get(g_iv_idx);
        
        if (decode_edge_target(get_edge_target(edge_list_idx)) == to) {
            // the edge back to the deleting node is the first in the list, so we need
            // to update the head
            graph_iv.set(g_iv_idx, get_next_edge_index(edge_list_idx));
        }
        else {
            // we need to traverse down the list and to find the edge back
            size_t prev_edge_list_idx = edge_list_idx;
            edge_list_idx = get_next_edge_index(edge_list_idx);
            while (decode_edge_target(get_edge_target(edge_list_idx)) != to) {
                prev_edge_list_idx = edge_list_idx;
                edge_list_idx = get_next_edge_index(edge_list_idx);
            }
            // skip over this edge in this linked list
            edge_lists_iv.set((prev_edge_list_idx - 1) * EDGE_RECORD_SIZE + EDGE_NEXT_OFFSET,
                              get_next_edge_index(edge_list_idx));
        }
        deleted_edge_records++;
    }
    
    void PackedGraph::destroy_edge(const handle_t& left, const handle_t& right) {
        remove_edge_reference(left, right);
        remove_edge_reference(flip(right), flip(left));
        defragment();
    }
    
    void PackedGraph::defragment(void) {
        if (deleted_node_records > defrag_factor * (graph_iv.size() / GRAPH_RECORD_SIZE)) {
            
            // adjust the start
            while (id_to_graph_iv.empty() ? false : id_to_graph_iv.get(0) == 0) {
                id_to_graph_iv.pop_front();
                min_id++;
            }
            // adjust the end
            while (id_to_graph_iv.empty() ? false : id_to_graph_iv.get(id_to_graph_iv.size() - 1) == 0) {
                id_to_graph_iv.pop_back();
            }
            
            PagedVector new_graph_iv(PAGE_WIDTH);
            PackedVector new_seq_length_iv;
            PagedVector new_seq_start_iv(PAGE_WIDTH);
            
            for (size_t i = 0; i < id_to_graph_iv.size(); i++) {
                size_t raw_g_iv_idx = id_to_graph_iv.get(i);
                if (raw_g_iv_idx) {
                    size_t g_iv_idx = (raw_g_iv_idx - 1) * GRAPH_RECORD_SIZE;
                    // this node still exists, create a new copy
                    new_graph_iv.append(graph_iv.get(g_iv_idx + GRAPH_START_EDGES_OFFSET));
                    new_graph_iv.append(graph_iv.get(g_iv_idx + GRAPH_END_EDGES_OFFSET));
                    new_seq_length_iv.append(graph_index_to_seq_len_index(g_iv_idx));
                    new_seq_start_iv.append(graph_index_to_seq_start_index(g_iv_idx));
                    // update the pointer into graph_iv
                    id_to_graph_iv.set(i, new_graph_iv.size() / GRAPH_RECORD_SIZE);
                }
            }
            
            // replace graph_iv with the defragged copy
            graph_iv = std::move(new_graph_iv);
            seq_length_iv = std::move(new_seq_length_iv);
            seq_start_iv = std::move(new_seq_start_iv);
            deleted_node_records = 0;
        }
        
        // TODO: defrag the seq_iv?
        
        if (deleted_edge_records > defrag_factor * (edge_lists_iv.size() / EDGE_RECORD_SIZE)) {
            
            PagedVector new_edge_lists_iv(PAGE_WIDTH);
            
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
    }
    
    void PackedGraph::clear(void) {
        graph_iv.clear();
        edge_lists_iv.clear();
        id_to_graph_iv.clear();
        seq_iv.clear();
        min_id = std::numeric_limits<id_t>::max();
        max_id = 0;
        deleted_edge_records = 0;
        deleted_node_records = 0;
    }
    
    bool PackedGraph::has_path(const std::string& path_name) const {
        return path_id.count(path_name);
    }
    
    path_handle_t PackedGraph::get_path_handle(const std::string& path_name) const {
        return as_path_handle(path_id.at(path_name));
    }
    
    string PackedGraph::get_path_name(const path_handle_t& path_handle) const {
        return paths.at(as_integer(path_handle)).first;
    }
    
    size_t PackedGraph::get_occurrence_count(const path_handle_t& path_handle) const {
        return paths.at(as_integer(path_handle)).second.size();
    }
    
    size_t PackedGraph::get_path_count() const {
        return paths.size();
    }
    
    void PackedGraph::for_each_path_handle(const std::function<void(const path_handle_t&)>& iteratee) const {
        
        for (const auto& path_record : paths) {
            iteratee(as_path_handle(path_record.first));
        }
        
    }
    
    handle_t PackedGraph::get_occurrence(const occurrence_handle_t& occurrence_handle) const {
        const PagedVector& path = paths.at(as_integers(occurrence_handle)[0]).second;
        uint64_t trav = path.get(as_integers(occurrence_handle)[1]);
        return reinterpret_cast<const handle_t&>(trav);
    }
    
    occurrence_handle_t PackedGraph::get_first_occurrence(const path_handle_t& path_handle) const {
        occurrence_handle_t occ;
        as_integers(occ)[0] = as_integer(path_handle);
        as_integers(occ)[1] = 0;
        return occ;
    }
    
    occurrence_handle_t PackedGraph::get_last_occurrence(const path_handle_t& path_handle) const {
        occurrence_handle_t occ;
        as_integers(occ)[0] = as_integer(path_handle);
        as_integers(occ)[1] = paths.at(as_integer(path_handle)).second.size() - 1;
        return occ;
    }
    
    bool PackedGraph::has_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
        return as_integers(occurrence_handle)[1] + 1 < paths.at(as_integers(occurrence_handle)[0]).second.size();
    }
    
    bool PackedGraph::has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
        return as_integers(occurrence_handle)[1] > 0;
    }
    
    occurrence_handle_t PackedGraph::get_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
        occurrence_handle_t next;
        as_integers(next)[0] = as_integers(occurrence_handle)[0];
        as_integers(next)[1] = as_integers(occurrence_handle)[0] + 1;
        return next;
    }
    
    occurrence_handle_t PackedGraph::get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
        occurrence_handle_t prev;
        as_integers(prev)[0] = as_integers(occurrence_handle)[0];
        as_integers(prev)[1] = as_integers(occurrence_handle)[0] - 1;
        return prev;
    }
    
    path_handle_t PackedGraph::get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
        return as_path_handle(as_integers(occurrence_handle)[0]);
    }
    
    size_t PackedGraph::get_ordinal_rank_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
        return as_integers(occurrence_handle)[1];
    }
    
    void PackedGraph::destroy_path(const path_handle_t& path) {
        auto& path_record = paths.at(as_integer(path));
        path_id.erase(path_record.first);
        paths.erase(as_integer(path));
    }
    
    path_handle_t PackedGraph::create_path_handle(const std::string& name) {
        path_id[name] = next_path_id;
        paths.emplace(next_path_id, pair<string, PagedVector>(name, PagedVector(PAGE_WIDTH)));
        path_handle_t path = as_path_handle(next_path_id);
        next_path_id++;
        return path;
        
    }
    
    occurrence_handle_t PackedGraph::append_occurrence(const path_handle_t& path, const handle_t& to_append) {
        PagedVector& path_vector = paths.at(as_integer(path)).second;
        
        occurrence_handle_t occ;
        as_integers(occ)[0] = as_integer(path);
        as_integers(occ)[1] = path_vector.size();
        
        path_vector.append(as_integer(to_append));
        
        return occ;
    }
}

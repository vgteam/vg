//
//  hash_graph.cpp
//

#include "hash_graph.hpp"

namespace vg {
    
    HashGraph::HashGraph() {
        
    }
    
    HashGraph::~HashGraph() {
        
    }
    
    
    handle_t HashGraph::create_handle(const string& sequence) {
        return create_handle(sequence, max_id + 1);
    }
    
    handle_t HashGraph::create_handle(const string& sequence, const id_t& id) {
        graph[id] = node_t(sequence);
        max_id = max(max_id, id);
        min_id = min(min_id, id);
        return get_handle(id, false);
    }
    
    void HashGraph::create_edge(const handle_t& left, const handle_t& right) {
        
        if (get_is_reverse(left)) {
            graph[get_id(left)].left_edges.push_back(right);
        }
        else {
            graph[get_id(left)].right_edges.push_back(right);
        }
        
        // a reversing self-edge only touches one side of one node, so we only want
        // to add it to a single edge list rather than two
        if (left != flip(right)){
            if (get_is_reverse(right)) {
                graph[get_id(right)].right_edges.push_back(flip(left));
            }
            else {
                graph[get_id(right)].left_edges.push_back(flip(left));
            }
        }
    }
    
    bool HashGraph::has_node(id_t node_id) const {
        return graph.count(node_id);
    }
    
    handle_t HashGraph::get_handle(const id_t& node_id, bool is_reverse) const {
        return EasyHandlePacking::pack(node_id, is_reverse);
    }
    
    id_t HashGraph::get_id(const handle_t& handle) const {
        return EasyHandlePacking::unpack_number(handle);
    }
    
    bool HashGraph::get_is_reverse(const handle_t& handle) const {
        return EasyHandlePacking::unpack_bit(handle);;
    }
    
    handle_t HashGraph::flip(const handle_t& handle) const {
        return EasyHandlePacking::toggle_bit(handle);
    }
    
    size_t HashGraph::get_length(const handle_t& handle) const {
        return graph.at(get_id(handle)).sequence.size();
    }
    
    string HashGraph::get_sequence(const handle_t& handle) const {
        return get_is_reverse(handle) ? reverse_complement(graph.at(get_id(handle)).sequence)
                                      : graph.at(get_id(handle)).sequence;
    }
    
    void HashGraph::swap_handles(const handle_t& a, const handle_t& b) {
        
        // TODO: implement?
    }
    
    bool HashGraph::follow_edges(const handle_t& handle, bool go_left,
                                 const std::function<bool(const handle_t&)>& iteratee) const {
        
        auto& edge_list = get_is_reverse(handle) != go_left ? graph.at(get_id(handle)).left_edges
                                                            : graph.at(get_id(handle)).right_edges;
        
        bool keep_going = true;
        for (auto it = edge_list.begin(); it != edge_list.end() && keep_going; it++) {
            keep_going = iteratee(go_left ? flip(*it) : *it);
        }
        return keep_going;
    }
    
    size_t HashGraph::node_size(void) const {
        return graph.size();
    }
    
    id_t HashGraph::min_node_id(void) const {
        return min_id;
    }
    
    id_t HashGraph::max_node_id(void) const {
        return max_id;
    }
    
    size_t HashGraph::get_degree(const handle_t& handle, bool go_left) const {
        auto& edge_list = get_is_reverse(handle) != go_left ? graph.at(get_id(handle)).left_edges
                                                            : graph.at(get_id(handle)).right_edges;
        return edge_list.size();
    }
    
    void HashGraph::for_each_handle(const std::function<bool(const handle_t&)>& iteratee,
                                    bool parallel) const {
        
        bool keep_going = true;
        if (parallel) {
#pragma omp parallel shared(keep_going)
            {
#pragma omp single
                {
                    for (auto it = graph.begin(); it != graph.end() && keep_going; it++) {
#pragma omp task
                        {
                            bool thread_keep_going = iteratee(get_handle(it->first));
#pragma omp atomic
                            keep_going = keep_going & thread_keep_going; // OMP only likes bitwise AND
                        }
                    }
                }
            }
            
        }
        else {
            for (auto it = graph.begin(); it != graph.end() && keep_going; it++) {
                keep_going = iteratee(get_handle(it->first));
            }
        }
    }
    
    handle_t HashGraph::apply_orientation(const handle_t& handle) {
        
        // don't do anything if it's already forward
        if (!get_is_reverse(handle)) {
            return handle;
        }
        
        // reverse the sequence
        node_t& node = graph[get_id(handle)];
        node.sequence = reverse_complement(node.sequence);
        
        // reverse the orientation
        for (vector<handle_t>* edge_list : {&node.left_edges, &node.right_edges}) {
            for (const handle_t& target : *edge_list) {
                node_t& other_node = graph[get_id(target)];
                auto& bwd_edge_list = get_is_reverse(target) ? other_node.right_edges : other_node.left_edges;
                for (handle_t& bwd_handle : bwd_edge_list) {
                    if (get_id(bwd_handle) == get_id(handle)) {
                        bwd_handle = flip(bwd_handle);
                        break;
                    }
                    // note: if a node has an edge to both sides of another node, we will end up flipping the
                    // same edge target twice, but this actually ends with the edge list in the desired state,
                    // so it's okay
                }
            }
        }
        
        // the edge lists switch sides
        swap(node.left_edges, node.right_edges);
        
        // update the occurrences on paths
        auto it = occurrences.find(get_id(handle));
        if (it != occurrences.end()) {
            for (auto& occurrence : it->second) {
                occurrence->handle = flip(occurrence->handle);
            }
        }
        
        // make it forward and return it
        return flip(handle);
    }
    
    vector<handle_t> HashGraph::divide_handle(const handle_t& handle, const vector<size_t>& offsets) {
        
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
        vector<handle_t> return_val;
        return_val.push_back(forward_handle);
        
        if (offsets.empty()) {
            return return_val;
        }
        
        // divvy up the sequence onto separate nodes
        for (size_t i = 0; i < forward_offsets.size(); i++) {
            size_t length = (i + 1 < forward_offsets.size() ? forward_offsets[i + 1] : node_length) - forward_offsets[i];
            return_val.push_back(create_handle(graph[get_id(handle)].sequence.substr(forward_offsets[i], length)));
        }
        graph[get_id(handle)].sequence = graph[get_id(handle)].sequence.substr(0, forward_offsets.front());
        
        // move the edges out the end of the node to the final one
        node_t& final_node = graph[get_id(return_val.back())];
        final_node.right_edges = move(graph[get_id(handle)].right_edges);
        graph[get_id(handle)].right_edges.clear();
        
        // update the backwards references back onto this node
        for (const handle_t& next : final_node.right_edges) {
            for (handle_t& bwd_target : get_is_reverse(next) ?
                                        graph[get_id(next)].right_edges :
                                        graph[get_id(next)].left_edges) {
                if (bwd_target == flip(forward_handle)) {
                    bwd_target = flip(return_val.back());
                    break;
                }
            }
        }
        
        // create edges between the segments of the original node
        for (size_t i = 1; i < return_val.size(); i++) {
            create_edge(return_val[i - 1], return_val[i]);
        }
        
        // update the paths and the occurrence records
        auto iter = occurrences.find(get_id(handle));
        if (iter != occurrences.end()) {
            for (path_mapping_t* mapping : iter->second) {
                path_t& path = paths[mapping->path_id];
                if (get_is_reverse(mapping->handle)) {
                    mapping = mapping->prev;
                    for (size_t i = return_val.size() - 1; i > 0; i--) {
                        mapping = path.insert_after(flip(return_val[i]), mapping);
                        occurrences[get_id(return_val[i])].push_back(mapping);
                    }
                }
                else {
                    for (size_t i = 1; i < return_val.size(); i++) {
                        mapping = path.insert_after(return_val[i], mapping);
                        occurrences[get_id(return_val[i])].push_back(mapping);
                    }
                }
            }
        }
        
        if (get_is_reverse(handle)) {
            // reverse the orientation of the return value to match the input
            reverse(return_val.begin(), return_val.end());
            for (handle_t& ret_handle : return_val) {
                ret_handle = flip(ret_handle);
            }
        }
        
        return return_val;
    }
    
    void HashGraph::destroy_handle(const handle_t& handle) {
        
        // remove backwards references from edges on other nodes
        node_t& node = graph[get_id(handle)];
        for (vector<handle_t>* edge_list : {&node.left_edges, &node.right_edges}) {
            for (const handle_t& next : *edge_list) {
                auto& bwd_edge_list = get_is_reverse(next) ? graph[get_id(next)].right_edges : graph[get_id(next)].left_edges;
                for (handle_t& bwd_target : bwd_edge_list) {
                    if (get_id(bwd_target) == get_id(handle)) {
                        bwd_target = bwd_edge_list.back();
                        bwd_edge_list.pop_back();
                        break;
                    }
                }
            }
        }
        
        // remove this node from the relevant indexes
        graph.erase(get_id(handle));
        occurrences.erase(get_id(handle));
    }
    
    void HashGraph::destroy_edge(const handle_t& left, const handle_t& right) {
        
        // remove this edge from left
        node_t& left_node = graph[get_id(left)];
        auto& left_edge_list = get_is_reverse(left) ? left_node.left_edges : left_node.right_edges;
        
        for (handle_t& next : left_edge_list) {
            if (next == right) {
                next = left_edge_list.back();
                left_edge_list.pop_back();
                break;
            }
        }
        
        // remove this edge from right
        node_t& right_node = graph[get_id(right)];
        auto& right_edge_list = get_is_reverse(right) ? right_node.right_edges : right_node.left_edges;
        
        for (handle_t& prev : right_edge_list) {
            if (prev == flip(left)) {
                prev = right_edge_list.back();
                right_edge_list.pop_back();
                break;
            }
        }
    }
    
    void HashGraph::clear(void) {
        max_id = 0;
        min_id = numeric_limits<id_t>::max();
        next_path_id = 1;
        graph.clear();
        path_id.clear();
        paths.clear();
        occurrences.clear();
    }
    
    bool HashGraph::has_path(const std::string& path_name) const {
        return path_id.count(path_name);
    }
    
    path_handle_t HashGraph::get_path_handle(const std::string& path_name) const {
        return as_path_handle(path_id.at(path_name));
    }
    
    string HashGraph::get_path_name(const path_handle_t& path_handle) const {
        return paths.at(as_integer(path_handle)).name;
    }
    
    size_t HashGraph::get_occurrence_count(const path_handle_t& path_handle) const {
        return paths.at(as_integer(path_handle)).count;
    }
    
    size_t HashGraph::get_path_count() const {
        return paths.size();
    }
    
    void HashGraph::for_each_path_handle(const std::function<void(const path_handle_t&)>& iteratee) const {
        for (auto it = paths.begin(); it != paths.end(); it++) {
            iteratee(as_path_handle(it->first));
        }
    }
    
    handle_t HashGraph::get_occurrence(const occurrence_handle_t& occurrence_handle) const {
        return ((path_mapping_t*) intptr_t(as_integers(occurrence_handle)[1]))->handle;
    }
    
    occurrence_handle_t HashGraph::get_first_occurrence(const path_handle_t& path_handle) const {
        occurrence_handle_t occ;
        as_integers(occ)[0] = as_integer(path_handle);
        as_integers(occ)[1] = intptr_t(paths.at(as_integer(path_handle)).head);
        return occ;
    }
    
    occurrence_handle_t HashGraph::get_last_occurrence(const path_handle_t& path_handle) const {
        occurrence_handle_t occ;
        as_integers(occ)[0] = as_integer(path_handle);
        as_integers(occ)[1] = intptr_t(paths.at(as_integer(path_handle)).tail);
        return occ;
    }
    
    bool HashGraph::has_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
        return ((path_mapping_t*) intptr_t(as_integers(occurrence_handle)[1]))->next != nullptr;
    }
    
    bool HashGraph::has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
        return ((path_mapping_t*) intptr_t(as_integers(occurrence_handle)[1]))->prev != nullptr;
    }
    
    occurrence_handle_t HashGraph::get_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
        occurrence_handle_t next;
        as_integers(next)[0] = as_integers(occurrence_handle)[0];
        as_integers(next)[1] = intptr_t(((path_mapping_t*) intptr_t(as_integers(occurrence_handle)[1]))->next);
        return next;
    }
    
    occurrence_handle_t HashGraph::get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
        occurrence_handle_t prev;
        as_integers(prev)[0] = as_integers(occurrence_handle)[0];
        as_integers(prev)[1] = intptr_t(((path_mapping_t*) intptr_t(as_integers(occurrence_handle)[1]))->prev);
        return prev;
    }
    
    path_handle_t HashGraph::get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
        return as_path_handle(as_integers(occurrence_handle)[0]);
    }
    
    vector<occurrence_handle_t> HashGraph::occurrences_of_handle(const handle_t& handle, bool match_orientation) const {
        vector<occurrence_handle_t> to_return;
        auto it = occurrences.find(get_id(handle));
        if (it != occurrences.end()) {
            for (path_mapping_t* mapping : it->second) {
                if (!match_orientation || mapping->handle == handle) {
                    to_return.emplace_back();
                    as_integers(to_return.back())[0] = mapping->path_id;
                    as_integers(to_return.back())[1] = intptr_t(mapping);
                }
            }
        }
        return to_return;
    }
    
    void HashGraph::destroy_path(const path_handle_t& path) {
        
        // remove the records of nodes occurring on this path
        for_each_occurrence_in_path(path, [&](const occurrence_handle_t& occ) {
            path_mapping_t* mapping = (path_mapping_t*) intptr_t(as_integers(occ)[1]);
            vector<path_mapping_t*>& node_occs = occurrences[get_id(mapping->handle)];
            for (size_t i = 0; i < node_occs.size(); i++) {
                if (node_occs[i] == mapping) {
                    node_occs[i] = node_occs.back();
                    node_occs.pop_back();
                    break;
                }
            }
        });
        
        // erase the path itself
        path_id.erase(paths[as_integer(path)].name);
        paths.erase(as_integer(path));
    }
    
    path_handle_t HashGraph::create_path_handle(const string& name) {
        path_id[name] = next_path_id;
        paths[next_path_id] = path_t(name, next_path_id);
        next_path_id++;
        return as_path_handle(next_path_id - 1);
    }
    
    occurrence_handle_t HashGraph::append_occurrence(const path_handle_t& path, const handle_t& to_append) {
        
        path_t& path_list = paths[as_integer(path)];
        path_mapping_t* mapping = path_list.push_back(to_append);
        occurrences[get_id(to_append)].push_back(mapping);
        
        occurrence_handle_t occ;
        as_integers(occ)[0] = as_integer(path);
        as_integers(occ)[1] = intptr_t(mapping);
        return occ;
    }
}

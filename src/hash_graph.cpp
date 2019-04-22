//
//  hash_graph.cpp
//

#include "hash_graph.hpp"

#include <handlegraph/util.hpp>

namespace vg {
    
    using namespace handlegraph;
    
    HashGraph::HashGraph() {
        
    }
    
    HashGraph::~HashGraph() {
        
    }
    
    HashGraph::HashGraph(istream& in) {
        deserialize(in);
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
        return handlegraph::number_bool_packing::pack(node_id, is_reverse);
    }
    
    id_t HashGraph::get_id(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_number(handle);
    }
    
    bool HashGraph::get_is_reverse(const handle_t& handle) const {
        return handlegraph::number_bool_packing::unpack_bit(handle);;
    }
    
    handle_t HashGraph::flip(const handle_t& handle) const {
        return handlegraph::number_bool_packing::toggle_bit(handle);
    }
    
    size_t HashGraph::get_length(const handle_t& handle) const {
        return graph.at(get_id(handle)).sequence.size();
    }
    
    string HashGraph::get_sequence(const handle_t& handle) const {
        return get_is_reverse(handle) ? reverse_complement(graph.at(get_id(handle)).sequence)
                                      : graph.at(get_id(handle)).sequence;
    }
    
    bool HashGraph::follow_edges_impl(const handle_t& handle, bool go_left,
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
    
    bool HashGraph::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee,
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
        
        return keep_going;
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
    
    size_t HashGraph::get_path_count() const {
        return paths.size();
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
    
    bool HashGraph::get_is_circular(const path_handle_t& path_handle) const {
        return paths.at(as_integer(path_handle)).is_circular;
    }
    
    size_t HashGraph::get_step_count(const path_handle_t& path_handle) const {
        return paths.at(as_integer(path_handle)).count;
    }
    
    bool HashGraph::for_each_path_handle_impl(const std::function<bool(const path_handle_t&)>& iteratee) const {
        for (auto it = paths.begin(); it != paths.end(); it++) {
            if (!iteratee(as_path_handle(it->first))) {
                return false;
            }
        }
        return true;
    }
    
    handle_t HashGraph::get_handle_of_step(const step_handle_t& step_handle) const {
        return ((path_mapping_t*) intptr_t(as_integers(step_handle)[1]))->handle;
    }
    
    step_handle_t HashGraph::path_begin(const path_handle_t& path_handle) const {
        step_handle_t step;
        as_integers(step)[0] = as_integer(path_handle);
        as_integers(step)[1] = intptr_t(paths.at(as_integer(path_handle)).head);
        return step;
    }
    
    step_handle_t HashGraph::path_end(const path_handle_t& path_handle) const {
        step_handle_t step;
        as_integers(step)[0] = as_integer(path_handle);
        as_integers(step)[1] = intptr_t(nullptr);
        return step;
    }
    
    step_handle_t HashGraph::get_next_step(const step_handle_t& step_handle) const {
        step_handle_t next;
        as_integers(next)[0] = as_integers(step_handle)[0];
        as_integers(next)[1] = intptr_t(((path_mapping_t*) intptr_t(as_integers(step_handle)[1]))->next);
        return next;
    }
    
    step_handle_t HashGraph::get_previous_step(const step_handle_t& step_handle) const {
        step_handle_t prev;
        as_integers(prev)[0] = as_integers(step_handle)[0];
        as_integers(prev)[1] = as_integers(step_handle)[1] != intptr_t(nullptr) ? intptr_t(((path_mapping_t*) intptr_t(as_integers(step_handle)[1]))->prev)
                                                                                : intptr_t(paths.at(as_integers(step_handle)[0]).tail);
        return prev;
    }
    
    path_handle_t HashGraph::get_path_handle_of_step(const step_handle_t& step_handle) const {
        return as_path_handle(as_integers(step_handle)[0]);
    }
    
    bool HashGraph::for_each_step_on_handle_impl(const handle_t& handle,
                                                 const function<bool(const step_handle_t&)>& iteratee) const {
        auto it = occurrences.find(get_id(handle));
        if (it != occurrences.end()) {
            for (path_mapping_t* mapping : it->second) {
                step_handle_t step;
                as_integers(step)[0] = mapping->path_id;
                as_integers(step)[1] = intptr_t(mapping);
                
                if (!iteratee(step)) {
                    return false;
                }
            }
        }
        return true;
    }
    
    void HashGraph::destroy_path(const path_handle_t& path) {
        
        // remove the records of nodes occurring on this path
        for_each_step_in_path(path, [&](const step_handle_t& step) {
            path_mapping_t* mapping = (path_mapping_t*) intptr_t(as_integers(step)[1]);
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
    
    path_handle_t HashGraph::create_path_handle(const string& name, bool is_circular) {
        path_id[name] = next_path_id;
        paths[next_path_id] = path_t(name, next_path_id, is_circular);
        next_path_id++;
        return as_path_handle(next_path_id - 1);
    }
    
    step_handle_t HashGraph::append_step(const path_handle_t& path, const handle_t& to_append) {
        
        path_t& path_list = paths[as_integer(path)];
        path_mapping_t* mapping = path_list.push_back(to_append);
        occurrences[get_id(to_append)].push_back(mapping);
        
        step_handle_t step;
        as_integers(step)[0] = as_integer(path);
        as_integers(step)[1] = intptr_t(mapping);
        return step;
    }
    
    void HashGraph::set_circularity(const path_handle_t& path, bool circular) {
        // set the annotation
        path_t& path_list = paths[as_integer(path)];
        path_list.is_circular = circular;
        
        // set the circular connection between the head and tail
        if (path_list.head) {
            if (circular) {
                path_list.head->prev = path_list.tail;
                path_list.tail->next = path_list.head;
            }
            else {
                path_list.head->prev = nullptr;
                path_list.tail->next = nullptr;
            }
        }
    }
    
    HashGraph::path_t::path_t() {
        
    }
    
    HashGraph::path_t::path_t(const string& name, const int64_t& path_id, bool is_circular) : name(name), path_id(path_id), is_circular(is_circular) {
        
        
    }
    
    HashGraph::path_t::path_t(path_t&& other) : head(other.head), tail(other.tail), path_id(other.path_id), count(other.count), name(move(other.name)), is_circular(other.is_circular) {
        // we grabbed the data in the initializer, now make sure the other one is in a valid state
        other.head = nullptr;
        other.tail = nullptr;
        other.path_id = 0;
        other.count = 0;
    }
    
    HashGraph::path_t& HashGraph::path_t::operator=(path_t&& other) {
        if (this != &other) {
            // free existing list
            this->~path_t();
            
            // steal other list
            head = other.head;
            tail = other.tail;
            other.head = nullptr;
            other.tail = nullptr;
            
            name = move(other.name);
            
            path_id = other.path_id;
            other.path_id = 0;
            
            count = other.count;
            other.count = 0;
        }
        return *this;
    }
    
    HashGraph::path_t::path_t(const path_t& other) : path_id(other.path_id), count(other.count), name(other.name), is_circular(other.is_circular) {
        
        path_mapping_t* prev = nullptr;
        bool first_iter = true;
        for (path_mapping_t* mapping = other.head;
             mapping != nullptr && (first_iter || mapping != other.head); // in case we loop around a circular path
             mapping = mapping->next) {
            
            path_mapping_t* copied = new path_mapping_t(mapping->handle, mapping->path_id);
            
            if (!head) {
                head = copied;
            }
            if (prev) {
                prev->next = copied;
                copied->prev = prev;
            }
            prev = copied;
            first_iter = false;
        }
        tail = prev;
        
        // we only loop over each mapping one time, so we may need to add in the last connection
        // in a circular path
        if (is_circular) {
            tail->next = head;
            head->prev = tail;
        }
    }
    
    HashGraph::path_t& HashGraph::path_t::operator=(const path_t& other) {
        if (this != &other) {
            // free existing list
            this->~path_t();
            
            // copy the other list
            path_mapping_t* prev = nullptr;
            bool first_iter = true;
            for (path_mapping_t* mapping = other.head;
                 mapping != nullptr && (first_iter || mapping != other.head); // in case we loop around a circular path
                 mapping = mapping->next) {
                
                path_mapping_t* copied = new path_mapping_t(mapping->handle, mapping->path_id);
                
                if (!head) {
                    head = copied;
                }
                if (prev) {
                    prev->next = copied;
                    copied->prev = prev;
                }
                prev = copied;
                first_iter = false;
            }
            tail = prev;
            
            // copy the rest of the info
            path_id = other.path_id;
            name = other.name;
            count = other.count;
            is_circular = other.is_circular;
            
            // we only loop over each mapping one time, so we may need to add in the last connection
            // in a circular path
            if (is_circular) {
                tail->next = head;
                head->prev = tail;
            }
        }
        return *this;
    }
    
    HashGraph::path_t::~path_t() {
        bool first_iter = true;
        for (path_mapping_t* mapping = head;
             mapping != nullptr && (first_iter || mapping != head);) { // for circular paths
             
            path_mapping_t* next = mapping->next;
            delete mapping;
            mapping = next;
            
            first_iter = false;
        }
    }
    
    HashGraph::path_mapping_t* HashGraph::path_t::push_back(const handle_t& handle) {
        return insert_after(handle, tail);
    }
    
    HashGraph::path_mapping_t* HashGraph::path_t::push_front(const handle_t& handle) {
        return insert_after(handle, nullptr);
    }
    
    void HashGraph::path_t::remove(path_mapping_t* mapping) {
        if (mapping == head) {
            head = mapping->next != mapping ? mapping->next : nullptr;
        }
        if (mapping == tail) {
            tail = mapping->prev != mapping ? mapping->prev : nullptr;
        }
        if (mapping->next) {
            mapping->next->prev = mapping->prev;
        }
        if (mapping->prev) {
            mapping->prev->next = mapping->next;
        }
        count--;
        delete mapping;
    }
    
    HashGraph::path_mapping_t* HashGraph::path_t::insert_after(const handle_t& handle, path_mapping_t* mapping) {
        
        path_mapping_t* inserting = new path_mapping_t(handle, path_id);
        
        if (mapping) {
            inserting->next = mapping->next;
            if (mapping->next) {
                mapping->next->prev = inserting;
            }
            mapping->next = inserting;
            inserting->prev = mapping;
            
            if (mapping == tail) {
                tail = inserting;
            }
        }
        else if (head) {
            
            // handle the potential circular connection
            inserting->prev = head->prev;
            if (inserting->prev) {
                inserting->prev->next = inserting;
            }
            
            inserting->next = head;
            head->prev = inserting;
            
            head = inserting;
        }
        else {
            // the list is empty so far, so initialize it
            head = tail = inserting;
            
            // make an initial circular connection
            if (is_circular) {
                inserting->next = inserting->prev = inserting;
            }
        }
        
        count++;
        return inserting;
    }
    
    void HashGraph::path_t::serialize(ostream& out) const {
        out.write((const char*) &path_id, sizeof(path_id) / sizeof(char));
        
        out.write((const char*) &is_circular, sizeof(is_circular) / sizeof(char));
        
        size_t name_size = name.size();
        out.write((const char*) &name_size, sizeof(name_size) / sizeof(char));
        out.write(name.c_str(), name.size());
        
        out.write((const char*) &count, sizeof(count) / sizeof(char));
        
        path_mapping_t* mapping = head;
        bool first_iter = true;
        while (mapping && (first_iter || mapping != head)) { // extra condition for circular paths
            
            int64_t step = as_integer(mapping->handle);
            out.write((const char*) &step, sizeof(step) / sizeof(char));
            mapping = mapping->next;
            
            first_iter = false;
        }
    }
    
    void HashGraph::path_t::deserialize(istream& in) {
        // free the current path if it exists
        this->~path_t();
        
        in.read((char*) &path_id, sizeof(path_id) / sizeof(char));
        
        in.read((char*) &is_circular, sizeof(is_circular) / sizeof(char));
        
        size_t name_size;
        in.read((char*) &name_size, sizeof(name_size) / sizeof(char));
        name.resize(name_size);
        for (size_t i = 0; i < name.size(); i++) {
            in.read((char*) &name[i], sizeof(char));
        }
        
        size_t num_mappings;
        in.read((char*) &num_mappings, sizeof(num_mappings) / sizeof(char));
        // note: count will be incremented in the push_back method
        count = 0;
        for (size_t i = 0; i < num_mappings; i++) {
            int64_t step;
            in.read((char*) &step, sizeof(step) / sizeof(char));
            push_back(as_handle(step));
        }
    }
    
    void HashGraph::node_t::serialize(ostream& out) const {
        size_t seq_size = sequence.size();
        out.write((const char*) &seq_size, sizeof(seq_size) / sizeof(char));
        out.write(sequence.c_str(), sequence.size());
        
        size_t left_edges_size = left_edges.size();
        out.write((const char*) &left_edges_size, sizeof(left_edges_size) / sizeof(char));
        for (size_t i = 0; i < left_edges.size(); i++) {
            int64_t next = as_integer(left_edges[i]);
            out.write((const char*) &next, sizeof(next) / sizeof(char));
        }
        
        size_t right_edges_size = right_edges.size();
        out.write((const char*) &right_edges_size, sizeof(right_edges_size) / sizeof(char));
        for (size_t i = 0; i < right_edges.size(); i++) {
            int64_t next = as_integer(right_edges[i]);
            out.write((const char*) &next, sizeof(next) / sizeof(char));
        }
    }
    
    void HashGraph::node_t::deserialize(istream& in) {
        
        size_t seq_size;
        in.read((char*) &seq_size, sizeof(seq_size) / sizeof(char));
        sequence.resize(seq_size);
        for (size_t i = 0; i < sequence.size(); i++) {
            in.read((char*) &sequence[i], sizeof(char));
        }
        
        size_t num_left_edges;
        in.read((char*) &num_left_edges, sizeof(num_left_edges) / sizeof(char));
        left_edges.resize(num_left_edges);
        for (size_t i = 0; i < left_edges.size(); i++) {
            int64_t next;
            in.read((char*) &next, sizeof(next) / sizeof(char));
            left_edges[i] = as_handle(next);
        }
        
        size_t num_right_edges;
        in.read((char*) &num_right_edges, sizeof(num_right_edges) / sizeof(char));
        right_edges.resize(num_right_edges);
        for (size_t i = 0; i < right_edges.size(); i++) {
            int64_t next;
            in.read((char*) &next, sizeof(next) / sizeof(char));
            right_edges[i] = as_handle(next);
        }
    }
    
    void HashGraph::serialize(ostream& out) const {
        out.write((const char*) &max_id, sizeof(max_id) / sizeof(char));
        out.write((const char*) &min_id, sizeof(min_id) / sizeof(char));
        out.write((const char*) &next_path_id, sizeof(next_path_id) / sizeof(char));
        
        size_t graph_size = graph.size();
        out.write((const char*) &graph_size, sizeof(graph_size) / sizeof(char));
        for (const pair<id_t, node_t>& node_record : graph) {
            out.write((const char*) &node_record.first, sizeof(node_record.first) / sizeof(char));
            node_record.second.serialize(out);
        }
        
        size_t paths_size = paths.size();
        out.write((const char*) &paths_size, sizeof(paths_size) / sizeof(char));
        for (const pair<int64_t, path_t>& path_record : paths) {
            path_record.second.serialize(out);
        }
    }
    
    void HashGraph::deserialize(istream& in) {
        clear();
        
        in.read((char*) &max_id, sizeof(max_id) / sizeof(char));
        in.read((char*) &min_id, sizeof(min_id) / sizeof(char));
        in.read((char*) &next_path_id, sizeof(next_path_id) / sizeof(char));
        
        size_t num_nodes;
        in.read((char*) &num_nodes, sizeof(num_nodes) / sizeof(char));
        graph.reserve(num_nodes);
        for (size_t i = 0; i < num_nodes; i++) {
            id_t node_id;
            in.read((char*) &node_id, sizeof(node_id) / sizeof(char));
            graph[node_id].deserialize(in);
        }
        
        size_t num_paths;
        in.read((char*) &num_paths, sizeof(num_paths) / sizeof(char));
        paths.reserve(num_paths);
        path_id.reserve(num_paths);
        for (size_t i = 0; i < num_paths; i++) {
            path_t path;
            path.deserialize(in);
            path_id[path.name] = path.path_id;
            paths[path.path_id] = move(path);
        }
        
        // we need to rebuild the occurrences of node mapping, which is not
        // part of the serialized format
        for (pair<const int64_t, path_t>& path_record : paths) {
            path_t& path = path_record.second;
            for (path_mapping_t* mapping = path.head; mapping != nullptr; mapping = mapping->next) {
                occurrences[get_id(mapping->handle)].push_back(mapping);
            }
        }
    }
}

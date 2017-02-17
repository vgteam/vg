//
//  snarls.cpp
//
//

#include "snarls.hpp"
#include "json2pb.h"

namespace vg {
    const vector<const Snarl*>& SnarlManager::children_of(const Snarl* snarl) {
        return children[key_form(snarl)];
    }
    
    const Snarl* SnarlManager::parent_of(const Snarl* snarl) {
        return parent[key_form(snarl)];
    }
    
    bool SnarlManager::is_leaf(const Snarl* snarl) {
        return children[key_form(snarl)].size() == 0;
    }
    
    bool SnarlManager::is_root(const Snarl* snarl) {
        return parent[key_form(snarl)] == nullptr;
    }
    
    const vector<const Snarl*>& SnarlManager::top_level_snarls() {
        return roots;
    }
    
    void SnarlManager::for_each_top_level_snarl_parallel(const function<void(const Snarl*)>& lambda) {
#pragma omp parallel for
        for (int i = 0; i < roots.size(); i++) {
            lambda(roots[i]);
        }
    }
    
    void SnarlManager::for_each_top_level_snarl(const function<void(const Snarl*)>& lambda) {
        for (const Snarl* snarl : roots) {
            lambda(snarl);
        }
    }
    
    void SnarlManager::flip(const Snarl* snarl) {
        
        // kinda ugly workaround of constness that avoids exposing non-const pointers while
        // still letting client edit Snarls in this controlled manner
        
        // get the offset of this snarl in the master list
        int64_t offset = (intptr_t) snarl - (intptr_t) snarls.data();
        // make sure this pointer is aligned to the
        if (offset % sizeof(Snarl) != 0) {
            cerr << "error:[SnarlManager] attempted to flip a Snarl with a pointer that is not owned by SnarlManager" << endl;
            assert(0);
        }
        
        Snarl& to_flip = snarls[offset / sizeof(Snarl)];
        
        // save the key used in the indices before editing the snarl
        auto old_key = key_form(snarl);
        
        // swap and reverse the start and end Visits
        int64_t start_id = to_flip.start().node_id();
        bool start_orientation = to_flip.start().backward();
        
        to_flip.mutable_start()->set_node_id(to_flip.end().node_id());
        to_flip.mutable_start()->set_backward(!to_flip.end().backward());
        
        to_flip.mutable_end()->set_node_id(start_id);
        to_flip.mutable_end()->set_backward(!start_orientation);
        
        // update parent index
        parent[key_form(snarl)] = parent[old_key];
        parent.erase(old_key);
        
        // update children index
        children[key_form(snarl)] = std::move(children[old_key]);
        children.erase(old_key);
        
    }
    
    const Snarl* SnarlManager::into_which_snarl(int64_t id, bool reverse) {
        return snarl_into.count(make_pair(id, reverse)) ? snarl_into[make_pair(id, reverse)] : nullptr;
    }
    
    // can include definition of inline function apart from forward declaration b/c only used in this file
    inline pair<pair<int64_t, bool>, pair<int64_t, bool> > SnarlManager::key_form(const Snarl* snarl) {
        return make_pair(make_pair(snarl->start().node_id(), snarl->start().backward()),
                         make_pair(snarl->end().node_id(), snarl->end().backward()));
    }
    
    void SnarlManager::build_indices() {
        
        for (Snarl& snarl : snarls) {
            // is this a top-level snarl?
            if (snarl.has_parent()) {
                // add this snarl to the parent-to-children index
                if (!children.count(key_form(&(snarl.parent()))) ) {
                    children.insert(make_pair(key_form(&snarl.parent()), vector<const Snarl*>(1, &snarl)));
                }
                else {
                    children[key_form(&(snarl.parent()))].push_back(&snarl);
                }
            }
            else {
                // record top level status
                roots.push_back(&snarl);
                parent[key_form(&snarl)] = nullptr;
            }
            
            // add the boundaries into the indices
            snarl_into[make_pair(snarl.start().node_id(), snarl.start().backward())] = &snarl;
            snarl_into[make_pair(snarl.end().node_id(), !snarl.end().backward())] = &snarl;
        }
        
        for (Snarl& snarl : snarls) {
            if (children.count(key_form(&snarl))) {
                // mark this snarl as the parent in child-to-parent map
                for (const Snarl* child : children[key_form(&snarl)]) {
                    parent[key_form(child)] = &snarl;
                }
            }
            else {
                // ensure that all snarls are in the parent-to-children map
                children.insert(make_pair(key_form(&snarl), vector<const Snarl*>()));
            }
        }
    }
    
    pair<unordered_set<Node*>, unordered_set<Edge*> > SnarlManager::shallow_contents(const Snarl* snarl, VG& graph,
                                                                                     bool include_boundary_nodes) {
        
        // construct maps that lets us "skip over" child snarls
        map<Node*, const Snarl*> child_snarl_starts;
        map<Node*, const Snarl*> child_snarl_ends;
        for (const Snarl* subsnarl : children_of(snarl)) {
            child_snarl_starts.insert(make_pair(graph.get_node(subsnarl->start().node_id()), subsnarl));
            child_snarl_ends.insert(make_pair(graph.get_node(subsnarl->end().node_id()), subsnarl));
        }
        
        pair<unordered_set<Node*>, unordered_set<Edge*> > to_return;
        
        unordered_set<Node*> already_stacked;
        
        // initialize stack for DFS traversal of site
        vector<Node*> stack;
        
        Node* start_node = graph.get_node(snarl->start().node_id());
        Node* end_node = graph.get_node(snarl->end().node_id());
        
        // mark the boundary nodes as already stacked so that paths will terminate on them
        already_stacked.insert(start_node);
        already_stacked.insert(end_node);
        
        // add boundary nodes as directed
        if (include_boundary_nodes) {
            to_return.first.insert(start_node);
            to_return.first.insert(end_node);
        }
        
        vector<Edge*> edges_of_node;
        
        // stack up the nodes one edge inside the snarl from the start
        graph.edges_of_node(start_node, edges_of_node);
        for (Edge* edge : edges_of_node) {
            
            // does the edge point into the snarl?
            if (edge->from() == snarl->start().node_id() && edge->from_start() == snarl->start().backward()) {

                Node* node = graph.get_node(edge->to());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }

                to_return.second.insert(edge);
            }
            else if (edge->to() == snarl->start().node_id() && edge->to_end() != snarl->start().backward()) {

                Node* node = graph.get_node(edge->from());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                to_return.second.insert(edge);
            }
        }
        edges_of_node.clear();
        
        // stack up the nodes one edge inside the snarl from the end
        graph.edges_of_node(end_node, edges_of_node);
        for (Edge* edge : edges_of_node) {
            // does the edge point into the snarl?
            if (edge->from() == snarl->end().node_id() && edge->from_start() != snarl->start().backward()) {
                
                Node* node = graph.get_node(edge->to());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                to_return.second.insert(edge);
            }
            else if (edge->to() == snarl->start().node_id() && edge->to_end() == snarl->start().backward()) {
                
                Node* node = graph.get_node(edge->from());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                to_return.second.insert(edge);
            }
        }
        edges_of_node.clear();
        
        // traverse the snarl with DFS, skipping over any child snarls
        // do not pay attention to valid walks since we also want to discover any tips
        while (stack.size()) {
            
            // pop the top node off the stack
            Node* node = stack.back();
            stack.pop_back();
            
            // record that this node is in the snarl
            to_return.first.insert(node);
            
            // are either the ends of the node facing into a snarl?
            bool forward_is_snarl = false;
            bool backward_is_snarl = false;
            
            if (child_snarl_starts.count(node)) {
                // this node is the start node of a snarl
                const Snarl* snarl = child_snarl_starts[node];
                if (snarl->start().backward()) {
                    backward_is_snarl = true;
                }
                else {
                    forward_is_snarl = true;
                }
                
                // stack up the node on the opposite side of the snarl
                // rather than traversing it
                Node* opposite_node = graph.get_node(snarl->end().node_id());
                if (!already_stacked.count(opposite_node)) {
                    stack.push_back(opposite_node);
                    already_stacked.insert(opposite_node);
                }
            }
            
            if (child_snarl_ends.count(node)) {
                // this node the end node of a snarl
                const Snarl* snarl = child_snarl_ends[node];
                if (snarl->end().backward()) {
                    forward_is_snarl = true;
                }
                else {
                    backward_is_snarl = true;
                }
                
                // stack up the node on the opposite side of the snarl
                // rather than traversing it
                Node* opposite_node = graph.get_node(snarl->start().node_id());
                if (!already_stacked.count(opposite_node)) {
                    stack.push_back(opposite_node);
                    already_stacked.insert(opposite_node);
                }
            }
            
            graph.edges_of_node(node, edges_of_node);
            
            for (Edge* edge : edges_of_node) {
                // which end of the edge is the current node?
                if (edge->from() == node->id()) {
                    // does this edge point forward or backward?
                    if ((edge->from_start() && !backward_is_snarl) ||
                        (!edge->from_start() && !forward_is_snarl)) {
                        
                        to_return.second.insert(edge);
                        Node* next_node = graph.get_node(edge->to());
                        
                        if (!already_stacked.count(next_node)) {
                            
                            stack.push_back(next_node);
                            already_stacked.insert(next_node);
                        }
                    }
                }
                else {
                    // does this edge point forward or backward?
                    if ((edge->to_end() && !forward_is_snarl) ||
                        (!edge->to_end() && !backward_is_snarl)) {
                        
                        to_return.second.insert(edge);
                        Node* next_node = graph.get_node(edge->from());
                        
                        if (!already_stacked.count(next_node)) {
                            
                            stack.push_back(next_node);
                            already_stacked.insert(next_node);
                        }
                    }
                }
            }
            
            edges_of_node.clear();
        }
        
        return to_return;
    }
    
    pair<unordered_set<Node*>, unordered_set<Edge*> > SnarlManager::deep_contents(const Snarl* snarl, VG& graph,
                                                                                  bool include_boundary_nodes) {
        
        pair<unordered_set<Node*>, unordered_set<Edge*> > to_return;
        
        unordered_set<Node*> already_stacked;
        
        // initialize stack for DFS traversal of site
        vector<Node*> stack;
        
        Node* start_node = graph.get_node(snarl->start().node_id());
        Node* end_node = graph.get_node(snarl->end().node_id());
        
        // mark the boundary nodes as already stacked so that paths will terminate on them
        already_stacked.insert(start_node);
        already_stacked.insert(end_node);
        
        // add boundary nodes as directed
        if (include_boundary_nodes) {
            to_return.first.insert(start_node);
            to_return.first.insert(end_node);
        }
        
        vector<Edge*> edges_of_node;
        
        // stack up the nodes one edge inside the snarl from the start
        graph.edges_of_node(start_node, edges_of_node);
        for (Edge* edge : edges_of_node) {
            // does the edge point into the snarl?
            if (edge->from() == snarl->start().node_id() && edge->from_start() == snarl->start().backward()) {
                
                Node* node = graph.get_node(edge->to());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                to_return.second.insert(edge);
            }
            else if (edge->to() == snarl->start().node_id() && edge->to_end() != snarl->start().backward()) {
                
                Node* node = graph.get_node(edge->from());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                to_return.second.insert(edge);
            }
        }
        edges_of_node.clear();
        
        // stack up the nodes one edge inside the snarl from the end
        graph.edges_of_node(end_node, edges_of_node);
        for (Edge* edge : edges_of_node) {
            // does the edge point into the snarl?
            if (edge->from() == snarl->end().node_id() && edge->from_start() != snarl->start().backward()) {
                
                Node* node = graph.get_node(edge->to());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                to_return.second.insert(edge);
            }
            else if (edge->to() == snarl->start().node_id() && edge->to_end() == snarl->start().backward()) {
                
                Node* node = graph.get_node(edge->from());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                to_return.second.insert(edge);
            }
        }
        edges_of_node.clear();
        
        // traverse the snarl with DFS, skipping over any child snarls
        // do not pay attention to valid walks since we also want to discover any tips
        while (stack.size()) {
            
            // pop the top node off the stack
            Node* node = stack.back();
            stack.pop_back();
            
            // record that this node is in the snarl
            to_return.first.insert(node);
            
            graph.edges_of_node(node, edges_of_node);
            
            for (Edge* edge : edges_of_node) {
                to_return.second.insert(edge);
                // get the other end of the edge
                Node* next_node = edge->from() == node->id() ? graph.get_node(edge->to()) :
                                                               graph.get_node(edge->from());
                if (!already_stacked.count(next_node)) {
                    stack.push_back(next_node);
                    already_stacked.insert(next_node);
                }
            }
            
            edges_of_node.clear();
        }
        
        return to_return;
    }
}







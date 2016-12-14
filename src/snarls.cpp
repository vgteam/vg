//
//  snarls.cpp
//
//

#include "snarls.hpp"

namespace vg {
    const vector<const Snarl*>& SnarlManager::children_of(const Snarl& snarl) {
        return children[key_form(snarl)];
    }
    
    const Snarl* SnarlManager::parent_of(const Snarl& snarl) {
        return parent[key_form(snarl)];
    }
    
    const vector<const Snarl*>& SnarlManager::top_level_snarls() {
        return roots;
    }
    
    void SnarlManager::for_each_top_level_snarl_parallel(const function<void(const Snarl&)>& lambda) {
#pragma omp parallel for
        for (int i = 0; i < roots.size(); i++) {
            lambda(*roots[i]);
        }
    }
    
    inline pair<pair<int64_t, bool>, pair<int64_t, bool> > SnarlManager::key_form(const Snarl& snarl) {
        return make_pair(make_pair(snarl.start().node_id(), snarl.start().backward()),
                         make_pair(snarl.end().node_id(), snarl.end().backward()));
    }
    
    void SnarlManager::build_trees() {
        
        for (Snarl& snarl : snarls) {
            
            // ensures that all snarls are in the children map
            if (!children.count(key_form(snarl))) {
                children.insert(make_pair(key_form(snarl), vector<const Snarl*>()));
            }
            
            // is this a top-level snarl?
            if (snarl.has_parent()) {
                // add parent to child-to-parent index
                parent.insert(make_pair(key_form(snarl), &snarl));
                
                // add this snarl to the parent-to-children index
                if (!children.count(key_form(snarl.parent()))) {
                    children.insert(make_pair(key_form(snarl.parent()), vector<const Snarl*>(1, &snarl)));
                }
                else {
                    children[key_form(snarl.parent())].push_back(&snarl);
                }
            }
            else {
                // add null parent to index
                parent.insert(make_pair(key_form(snarl), nullptr));
                
                // record top level status
                roots.push_back(&snarl);
            }
        }
    }
    
    pair<vector<Node*>, vector<Edge*> > SnarlManager::shallow_contents(const Snarl& snarl, VG& graph) {
        
        // construct maps that lets us "skip over" child snarls
        map<Node*, const Snarl*> child_snarl_starts;
        map<Node*, const Snarl*> child_snarl_ends;
        for (const Snarl* subsnarl : children_of(snarl)) {
            child_snarl_starts.insert(make_pair(graph.get_node(subsnarl->start().node_id()), subsnarl));
            child_snarl_ends.insert(make_pair(graph.get_node(subsnarl->end().node_id()), subsnarl));
        }
        
        unordered_set<Node*> nodes;
        unordered_set<Edge*> edges;
        
        set<Node*> already_stacked;
        
        // initialize stack for DFS traversal of site
        vector<Node*> stack;
        
        Node* start_node = graph.get_node(snarl.start().node_id());
        Node* end_node = graph.get_node(snarl.start().node_id());
        
        // mark the boundary nodes as already stacked so that paths will terminate on them
        already_stacked.insert(start_node);
        already_stacked.insert(end_node);
        
        // note: do not record the boundary nodes in the node list
        
        vector<Edge*> edges_of_node;
        
        // stack up the nodes one edge inside the snarl from the start
        graph.edges_of_node(start_node, edges_of_node);
        for (Edge* edge : edges_of_node) {
            // does the edge point into the snarl?
            if (edge->from() == snarl.start().node_id() && edge->from_start() == snarl.start().backward()) {

                Node* node = graph.get_node(edge->to());

                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }

                edges.insert(edge);
            }
            else if (edge->to() == snarl.start().node_id() && edge->to_end() != snarl.start().backward()) {

                Node* node = graph.get_node(edge->from());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                edges.insert(edge);
            }
        }
        edges_of_node.clear();
        
        // stack up the nodes one edge inside the snarl from the end
        graph.edges_of_node(end_node, edges_of_node);
        for (Edge* edge : edges_of_node) {
            // does the edge point into the snarl?
            if (edge->from() == snarl.end().node_id() && edge->from_start() != snarl.start().backward()) {
                
                Node* node = graph.get_node(edge->to());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                edges.insert(edge);
            }
            else if (edge->to() == snarl.start().node_id() && edge->to_end() == snarl.start().backward()) {
                
                Node* node = graph.get_node(edge->from());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                edges.insert(edge);
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
            nodes.insert(node);
            
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
                        
                        edges.insert(edge);
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
                        
                        edges.insert(edge);
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
        
        // put the nodes and edges that DFS traversed into the return vector
        pair<vector<Node*>, vector<Edge*> > to_return;
        
        for (Node* node : nodes) {
            to_return.first.push_back(node);
        }
        
        for (Edge* edge : edges) {
            to_return.second.push_back(edge);
        }
        
        return to_return;
    }
    
    pair<vector<Node*>, vector<Edge*> > SnarlManager::deep_contents(const Snarl& snarl, VG& graph) {
        
        unordered_set<Node*> nodes;
        unordered_set<Edge*> edges;
        
        set<Node*> already_stacked;
        
        // initialize stack for DFS traversal of site
        vector<Node*> stack;
        
        Node* start_node = graph.get_node(snarl.start().node_id());
        Node* end_node = graph.get_node(snarl.start().node_id());
        
        // mark the boundary nodes as already stacked so that paths will terminate on them
        already_stacked.insert(start_node);
        already_stacked.insert(end_node);
        
        // note: do not record the boundary nodes in the node list
        
        vector<Edge*> edges_of_node;
        
        // stack up the nodes one edge inside the snarl from the start
        graph.edges_of_node(start_node, edges_of_node);
        for (Edge* edge : edges_of_node) {
            // does the edge point into the snarl?
            if (edge->from() == snarl.start().node_id() && edge->from_start() == snarl.start().backward()) {
                
                Node* node = graph.get_node(edge->to());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                edges.insert(edge);
            }
            else if (edge->to() == snarl.start().node_id() && edge->to_end() != snarl.start().backward()) {
                
                Node* node = graph.get_node(edge->from());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                edges.insert(edge);
            }
        }
        edges_of_node.clear();
        
        // stack up the nodes one edge inside the snarl from the end
        graph.edges_of_node(end_node, edges_of_node);
        for (Edge* edge : edges_of_node) {
            // does the edge point into the snarl?
            if (edge->from() == snarl.end().node_id() && edge->from_start() != snarl.start().backward()) {
                
                Node* node = graph.get_node(edge->to());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                edges.insert(edge);
            }
            else if (edge->to() == snarl.start().node_id() && edge->to_end() == snarl.start().backward()) {
                
                Node* node = graph.get_node(edge->from());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                edges.insert(edge);
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
            nodes.insert(node);
            
            graph.edges_of_node(node, edges_of_node);
            
            for (Edge* edge : edges_of_node) {
                edges.insert(edge);
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
        
        // put the nodes and edges that DFS traversed into the return vector
        pair<vector<Node*>, vector<Edge*> > to_return;
        
        for (Node* node : nodes) {
            to_return.first.push_back(node);
        }
        
        for (Edge* edge : edges) {
            to_return.second.push_back(edge);
        }
        
        return to_return;
    }
}







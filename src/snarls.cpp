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
    
    map<NodeTraversal, const Snarl*> SnarlManager::child_boundary_index(const Snarl* snarl, VG& graph) {
        map<NodeTraversal, const Snarl*> index;
        for (const Snarl* child : children_of(snarl)) {
            index[to_node_traversal(child->start(), graph)] = child;
            index[to_rev_node_traversal(child->end(), graph)] = child;
        }
        return index;
    }
    
    map<NodeTraversal, const Snarl*> SnarlManager::child_start_index(const Snarl* snarl, VG& graph) {
        map<NodeTraversal, const Snarl*> index;
        for (const Snarl* child : children_of(snarl)) {
            index[to_node_traversal(child->start(), graph)] = child;
        }
        return index;
    }
    
    map<NodeTraversal, const Snarl*> SnarlManager::child_end_index(const Snarl* snarl, VG& graph) {
        map<NodeTraversal, const Snarl*> index;
        for (const Snarl* child : children_of(snarl)) {
            index[to_rev_node_traversal(child->end(), graph)] = child;
        }
        return index;
    }
    
    // can include definition of inline function apart from forward declaration b/c only used in this file
    inline pair<pair<int64_t, bool>, pair<int64_t, bool> > SnarlManager::key_form(const Snarl* snarl) {
        return make_pair(make_pair(snarl->start().node_id(), snarl->start().backward()),
                         make_pair(snarl->end().node_id(), snarl->end().backward()));
    }
    
    void SnarlManager::build_trees() {
        
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
        
        // construct the map that lets us skip over child snarls
        map<NodeTraversal, const Snarl*> child_index = child_boundary_index(snarl, graph);
        
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
            if (edge->from() == snarl->end().node_id() && edge->from_start() != snarl->end().backward()) {
                
                Node* node = graph.get_node(edge->to());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                to_return.second.insert(edge);
            }
            else if (edge->to() == snarl->end().node_id() && edge->to_end() == snarl->end().backward()) {
                
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
            
            if (child_index.count(NodeTraversal(node, false))) {
                // this node forward reads into a snarl
                forward_is_snarl = true;
                
                const Snarl* snarl = child_index[NodeTraversal(node, false)];
                
                // What's on the other side of the snarl?
                id_t other_id = snarl->start().node_id() == node->id() ? snarl->end().node_id(): snarl->start().node_id();
                
                // stack up the node on the opposite side of the snarl
                // rather than traversing it
                Node* opposite_node = graph.get_node(other_id);
                if (!already_stacked.count(opposite_node)) {
                    stack.push_back(opposite_node);
                    already_stacked.insert(opposite_node);
                }
            }
            
            if (child_index.count(NodeTraversal(node, true))) {
                // this node reverse reads into a snarl
                backward_is_snarl = true;
                
                const Snarl* snarl = child_index[NodeTraversal(node, true)];
                
                // What's on the other side of the snarl?
                id_t other_id = snarl->end().node_id() == node->id() ? snarl->start().node_id(): snarl->end().node_id();
                
                // stack up the node on the opposite side of the snarl
                // rather than traversing it
                Node* opposite_node = graph.get_node(other_id);
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
            if (edge->from() == snarl->end().node_id() && edge->from_start() != snarl->end().backward()) {
                
                Node* node = graph.get_node(edge->to());
                
                if (!already_stacked.count(node)) {
                    stack.push_back(node);
                    already_stacked.insert(node);
                }
                
                to_return.second.insert(edge);
            }
            else if (edge->to() == snarl->end().node_id() && edge->to_end() == snarl->end().backward()) {
                
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
    
    const Snarl* SnarlManager::manage(const Snarl& not_owned) {
        // TODO: keep the Snarls in some kind of sorted order to make lookup
        // efficient. We could also have a map<Snarl, Snarl*> but that would be
        // a tremendous waste of space.
        
        // Right now we're stuck with O(n) search and it's horrible
        for(auto& owned : snarls) {
            // Only compare start and end visits, because we may get snarls with
            // only those set from visits.
            if (owned.start() == not_owned.start() && owned.end() == not_owned.end()) {
                return &owned;
            }
        }
        // If we get here it doesn't exist.
        throw runtime_error("Unable to find snarl " +  pb2json(not_owned) + " in SnarlManager");
    }
    
    vector<Visit> visits_right(const Visit& visit, VG& graph,
        const map<NodeTraversal, const Snarl*>& child_boundary_index) {
        
#ifdef debug
        cerr << "Look right from " << visit << endl;
        cerr << "\tChild index:" << endl;
        for (auto& kv : child_boundary_index) {
            cerr << "\t\t" << kv.first << " -> " << *kv.second << endl;
        }
#endif
        
        // We'll populate this
        vector<Visit> to_return;
        
        // Find the right side of the visit we're on
        NodeSide right_side = to_right_side(visit);
        
        if (visit.node_id() == 0) {
            // We're leaving a child snarl, so we are going to need to check if
            // another child snarl shares this boundary node in the direction
            // we're going.
            
            // Make a node traversal for going towards that right side
            NodeTraversal out_of_child(graph.get_node(right_side.node), !right_side.is_end);
            
            if (child_boundary_index.count(out_of_child)) {
                // We leave the one child and immediately enter another!
                
                // Make a visit to it
                Visit child_visit;
                transfer_boundary_info(*child_boundary_index.at(out_of_child), *child_visit.mutable_snarl());
                
                if (right_side.node == child_visit.snarl().end().node_id()) {
                    // We came in its end
                    child_visit.set_backward(true);
                } else {
                    // We should have come in its start
                    assert(right_side.node == child_visit.snarl().start().node_id());
                }
                
                // Bail right now, so we don'ttry to explore inside this child snarl.
                to_return.push_back(child_visit);
                return to_return;
                
            }
            
        }
        
        for (auto attached : graph.sides_of(right_side)) {
            // For every NodeSide attached to the right side of this visit
            
            // Make it a NodeTraversal reading away from that side
            NodeTraversal attached_traversal(graph.get_node(attached.node), attached.is_end);
            
#ifdef debug
            cerr << "\tFind NodeTraversal " << attached_traversal << endl;
#endif
            
            if (child_boundary_index.count(attached_traversal)) {
                // We're reading into a child
                
                // Which child is it?
                const Snarl* child = child_boundary_index.at(attached_traversal);
                
#ifdef debug
                cerr << "\t\tGoes to Snarl " << *child << endl;
#endif
                
                if (attached.node == child->start().node_id()) {
                    // We're reading into the start of the child
                    
                    // Make a visit to the child snarl
                    Visit child_visit;
                    transfer_boundary_info(*child, *child_visit.mutable_snarl());
                    
#ifdef debug
                    cerr << "\t\tProduces Visit " << child_visit << endl;
#endif
                    
                    // Put it in in the forward orientation
                    to_return.push_back(child_visit);
                } else if (attached.node == child->end().node_id()) {
                    // We're reading into the end of the child
                    
                    // Make a visit to the child snarl
                    Visit child_visit;
                    transfer_boundary_info(*child, *child_visit.mutable_snarl());
                    child_visit.set_backward(true);
                    
#ifdef debug
                    cerr << "\t\tProduces Visit " << child_visit << endl;
#endif
                    
                    // Put it in in the reverse orientation
                    to_return.push_back(child_visit);
                } else {
                    // Should never happen
                    throw runtime_error("Read into child " + pb2json(*child) + " with non-matching traversal");
                }
            } else {
                // We just go into a normal node
                to_return.push_back(to_visit(attached_traversal));
            
#ifdef debug
                cerr << "\t\tProduces Visit " << to_return.back() << endl;
#endif
                
            }
        }
        
        return to_return;
        
    }
    
    vector<Visit> visits_left(const Visit& visit, VG& graph,
        const map<NodeTraversal, const Snarl*>& child_boundary_index) {
        
        // Reverse the visit
        Visit reversed = visit;
        reversed.set_backward(!reversed.backward());
        
        // Get everything right of the reversed version
        vector<Visit> to_return = visits_right(reversed, graph, child_boundary_index);
        
        // Un-reverse them so they are in the correct orientation to be seen
        // left of here.
        for (auto& v : to_return) {
            v = reverse(v);
        }
        
        return to_return;
        
    }
    
    bool operator==(const Visit& a, const Visit& b) {
        // IDs and orientations have to match, and nobody has a snarl or the
        // snarls match.
        return a.node_id() == b.node_id() &&
            a.backward() == b.backward() &&
            ((!a.has_snarl() && !b.has_snarl()) ||
            a.snarl() == b.snarl());
    }
    
    bool operator!=(const Visit& a, const Visit& b) {
        return !(a == b);
    }
    
    bool operator<(const Visit& a, const Visit& b) {
        if (!a.has_snarl() && !b.has_snarl()) {
            // Compare everything but the snarl
            return make_tuple(a.node_id(), a.backward()) < make_tuple(b.node_id(), b.backward());
        } else {
            // Compare including the snarl
            return make_tuple(a.node_id(), a.snarl(), a.backward()) < make_tuple(b.node_id(), b.snarl(), b.backward());
        }        
    }
    
    ostream& operator<<(ostream& out, const Visit& visit) {
        if (!visit.has_snarl()) {
            // Use the node ID
            out << visit.node_id();
        } else {
            // Use the snarl
            out << visit.snarl();
        }
        out << " " << (visit.backward() ? "rev" : "fwd");
        return out;
    }
    
    bool operator==(const SnarlTraversal& a, const SnarlTraversal& b) {
        if (a.snarl() != b.snarl()) {
            return false;
        }
        if (a.visits_size() != b.visits_size()) {
            return false;
        }
        for (size_t i = 0; i < a.visits_size(); i++) {
            if (a.visits(i) != b.visits(i)) {
                return false;
            }
        }
        // Otherwise everything we can think of matches
        return true;
    }
    
    bool operator!=(const SnarlTraversal& a, const SnarlTraversal& b) {
        return !(a == b);
    }
    
    bool operator<(const SnarlTraversal& a, const SnarlTraversal& b) {
        if (a.snarl() < b.snarl()) {
            return true;
        } else if (b.snarl() < a.snarl()) {
            return false;
        }
        for (size_t i = 0; i < b.visits_size(); i++) {
            if (i >= a.visits_size()) {
                // A has run out and B is still going
                return true;
            }
            
            if (a.visits(i) < b.visits(i)) {
                return true;
            } else if (b.visits(i) < a.visits(i)) {
                return false;
            }
        }
        
        // If we get here either they're equal or A has more visits than B
        return false;
    }
    
    bool operator==(const Snarl& a, const Snarl& b) {
        if (a.type() != b.type()) {
            return false;
        }
        if (a.start() != b.start()) {
            return false;
        }
        if (a.end() != b.end()) {
            return false;
        }
        if (a.has_parent() || b.has_parent()) {
            // Someone has a parent so we must compare them.
            return a.parent() == b.parent();
        }
        return true;
    }
    
    bool operator!=(const Snarl& a, const Snarl& b) {
        return !(a == b);
    }
    
    bool operator<(const Snarl& a, const Snarl& b) {
        if (!a.has_parent() && !b.has_parent()) {
            // Compare without parent
            return make_tuple(a.type(), a.start(), a.end()) < make_tuple(b.type(), b.start(), b.end());
        } else {
            // Compare with parent
            return make_tuple(a.type(), a.start(), a.end(), a.parent()) < make_tuple(b.type(), b.start(), b.end(), b.parent());
        }
    }
    
    ostream& operator<<(ostream& out, const Snarl& snarl) {
        return out << snarl.start() << "-" << snarl.end();
    }
}







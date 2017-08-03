//
//  snarls.cpp
//
//

//#define debug

#include "snarls.hpp"
#include "json2pb.h"

namespace vg {
    // TODO: this is duplicative with the other constructor, but protobuf won't let me make
    // a deserialization iterator to match its signature because its internal file streams
    // disallow copy constructors
    SnarlManager::SnarlManager(istream& in) {
        // add snarls to master list
        for (stream::ProtobufIterator<Snarl> iter(in); iter.has_next(); iter.get_next()) {
            snarls.push_back(*iter);
        }
        // record the tree structure and build the other indexes
        build_indexes();
    }
    
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
        // make sure this pointer is aligned to the vector
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
        
        // Update index index
        index_of[key_form(snarl)] = std::move(index_of[old_key]);
        index_of.erase(old_key);
        
        // note: snarl_into index is invariant to flipping
    }
    
    const Snarl* SnarlManager::into_which_snarl(int64_t id, bool reverse) {
        return snarl_into.count(make_pair(id, reverse)) ? snarl_into[make_pair(id, reverse)] : nullptr;
    }
    
    const Snarl* SnarlManager::into_which_snarl(const Visit& visit) {
        return visit.has_snarl() ? manage(visit.snarl()) : into_which_snarl(visit.node_id(), visit.backward());
    }
    
    unordered_map<pair<int64_t, bool>, const Snarl*> SnarlManager::snarl_boundary_index() {
        unordered_map<pair<int64_t, bool>, const Snarl*> index;
        for (Snarl& snarl : snarls) {
            index[make_pair(snarl.start().node_id(), snarl.start().backward())] = &snarl;
            index[make_pair(snarl.end().node_id(), !snarl.end().backward())] = &snarl;
        }
        return index;
    }
    
    unordered_map<pair<int64_t, bool>, const Snarl*> SnarlManager::snarl_end_index() {
        unordered_map<pair<int64_t, bool>, const Snarl*> index;
        for (Snarl& snarl : snarls) {
            index[make_pair(snarl.end().node_id(), !snarl.end().backward())] = &snarl;
        }
        return index;
    }
    
    unordered_map<pair<int64_t, bool>, const Snarl*> SnarlManager::snarl_start_index() {
        unordered_map<pair<int64_t, bool>, const Snarl*> index;
        for (Snarl& snarl : snarls) {
            index[make_pair(snarl.start().node_id(), snarl.start().backward())] = &snarl;
        }
        return index;
    }
    
    // can include definition of inline function apart from forward declaration b/c only used in this file
    inline SnarlManager::key_t SnarlManager::key_form(const Snarl* snarl) {
        return make_pair(make_pair(snarl->start().node_id(), snarl->start().backward()),
                         make_pair(snarl->end().node_id(), snarl->end().backward()));
    }
    

    void SnarlManager::build_indexes() {
        
        for (size_t i = 0; i < snarls.size(); i++) {
            Snarl& snarl = snarls[i];
            
            // Remember where each snarl is
            index_of[key_form(&snarl)] = i;
        
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
            
            const Snarl* forward_snarl = into_which_snarl(node->id(), false);
            const Snarl* backward_snarl = into_which_snarl(node->id(), true);
            if (forward_snarl) {
                // this node points into a snarl
                
                // What's on the other side of the snarl?
                id_t other_id = forward_snarl->start().node_id() == node->id() ? forward_snarl->end().node_id() : forward_snarl->start().node_id();
                
                // stack up the node on the opposite side of the snarl
                // rather than traversing it
                Node* opposite_node = graph.get_node(other_id);
                if (!already_stacked.count(opposite_node)) {
                    stack.push_back(opposite_node);
                    already_stacked.insert(opposite_node);
                }
            }
            
            if (backward_snarl) {
                // the reverse of this node points into a snarl
                
                // What's on the other side of the snarl?
                id_t other_id = backward_snarl->end().node_id() == node->id() ? backward_snarl->start().node_id(): backward_snarl->end().node_id();
                
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
                    if ((edge->from_start() && !backward_snarl) ||
                        (!edge->from_start() && !forward_snarl)) {
                        
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
                    if ((edge->to_end() && !forward_snarl) ||
                        (!edge->to_end() && !backward_snarl)) {
                        
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
        
        // Work out the key for the snarl
        key_t key = key_form(&not_owned);
        
        // Get the index of the snarl with that key.
        auto it = index_of.find(key);
        
        if (it == index_of.end()) {
            // It's not there. Someone is trying to manage a snarl we don't
            // really own. Complain.
            throw runtime_error("Unable to find snarl " +  pb2json(not_owned) + " in SnarlManager");
        }
        
        // Return the official copy of that snarl
        return &snarls.at(it->second);
    }
    
    vector<Visit> SnarlManager::visits_right(const Visit& visit, VG& graph, const Snarl* in_snarl) {
        
#ifdef debug
        cerr << "Look right from " << visit << endl;
#endif
        
        // We'll populate this
        vector<Visit> to_return;
        
        // Find the right side of the visit we're on
        NodeSide right_side = to_right_side(visit);
        
        if (visit.node_id() == 0) {
            // We're leaving a child snarl, so we are going to need to check if
            // another child snarl shares this boundary node in the direction
            // we're going.
            
            const Snarl* child = into_which_snarl(right_side.node, !right_side.is_end);
            if (child != nullptr && child != in_snarl
                && into_which_snarl(right_side.node, right_side.is_end) != in_snarl) {
                // We leave the one child and immediately enter another!
                
                // Make a visit to it
                Visit child_visit;
                transfer_boundary_info(*child, *child_visit.mutable_snarl());
                
                if (right_side.node == child->end().node_id()) {
                    // We came in its end
                    child_visit.set_backward(true);
                } else {
                    // We should have come in its start
                    assert(right_side.node == child->start().node_id());
                }
                
                // Bail right now, so we don't try to explore inside this child snarl.
                to_return.push_back(child_visit);
                return to_return;
                
            }
            
        }
        
        for (auto attached : graph.sides_of(right_side)) {
            // For every NodeSide attached to the right side of this visit
            
#ifdef debug
            cerr << "\tFind NodeSide " << attached << endl;
#endif
            
            const Snarl* child = into_which_snarl(attached.node, attached.is_end);
            if (child != nullptr && child != in_snarl
                && into_which_snarl(attached.node, !attached.is_end) != in_snarl) {
                // We're reading into a child
                
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
                to_return.emplace_back();
                Visit& next_visit = to_return.back();
                next_visit.set_node_id(attached.node);
                next_visit.set_backward(attached.is_end);
            
#ifdef debug
                cerr << "\t\tProduces Visit " << to_return.back() << endl;
#endif
                
            }
        }
        
        return to_return;
        
    }
    
    vector<Visit> SnarlManager::visits_left(const Visit& visit, VG& graph, const Snarl* in_snarl) {
        
        // Get everything right of the reversed visit
        vector<Visit> to_return = visits_right(reverse(visit), graph, in_snarl);
        
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







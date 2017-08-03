//
//  snarls.hpp
//
//  Contains object to own Snarls and keep track of their tree relationships as well as utility
//  functions that interact with snarls.
//

#ifndef snarls_hpp
#define snarls_hpp

#include <cstdint>
#include <stdio.h>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include "stream.hpp"
#include "vg.hpp"
#include "vg.pb.h"
#include "hash_map.hpp"


using namespace std;

namespace vg {
    
    /**
     * A structure to keep track of the tree relationships between Snarls and perform utility algorithms
     * on them
     */
    class SnarlManager {
    public:
        
        /// Construct a SnarlManager for the snarls returned by an iterator
        template <typename SnarlIterator>
        SnarlManager(SnarlIterator begin, SnarlIterator end);
        
        /// Construct a SnarlManager for the snarls contained in an input stream
        SnarlManager(istream& in);
        
        /// Default constructor
        SnarlManager() = default;
        
        /// Destructor
        ~SnarlManager() = default;
        
        /// Returns a vector of pointers to the children of a Snarl
        const vector<const Snarl*>& children_of(const Snarl* snarl);
        
        /// Returns a pointer to the parent of a Snarl or nullptr if there is none
        const Snarl* parent_of(const Snarl* snarl);
        
        /// Returns the Snarl that a traversal points into at either the start or end, or nullptr if
        /// the traversal does not point into any Snarl. Note that Snarls store the end Visit pointing
        /// out of rather than into the Snarl, so they must be reversed to query it.
        const Snarl* into_which_snarl(int64_t id, bool reverse);
        /// Returns the Snarl that a Visit points into. If the Visit contains a Snarl rather than a node
        /// ID, returns a pointer the managed version of that snarl.
        const Snarl* into_which_snarl(const Visit& visit);
        
        /// Returns true if snarl has no children and false otherwise
        bool is_leaf(const Snarl* snarl);
        
        /// Returns true if snarl has no parent and false otherwise
        bool is_root(const Snarl* snarl);
        
        /// Returns a reference to a vector with the roots of the Snarl trees
        const vector<const Snarl*>& top_level_snarls();
        
        /// Reverses the orientation of a snarl
        void flip(const Snarl* snarl);
        
        /// Returns the Nodes and Edges contained in this Snarl but not in any child Snarls (always includes the
        /// Nodes that form the boundaries of child Snarls, optionally includes this Snarl's own boundary Nodes)
        pair<unordered_set<Node*>, unordered_set<Edge*> > shallow_contents(const Snarl* snarl, VG& graph,
                                                                           bool include_boundary_nodes);
        
        /// Returns the Nodes and Edges contained in this Snarl, including those in child Snarls (optionally
        /// includes Snarl's own boundary Nodes)
        pair<unordered_set<Node*>, unordered_set<Edge*> > deep_contents(const Snarl* snarl, VG& graph,
                                                                        bool include_boundary_nodes);
        
        /// Look left from the given visit in the given graph and gets all the
        /// attached Visits to nodes or snarls.
        vector<Visit> visits_left(const Visit& visit, VG& graph, const Snarl* in_snarl);
        
        /// Look left from the given visit in the given graph and gets all the
        /// attached Visits to nodes or snarls.
        vector<Visit> visits_right(const Visit& visit, VG& graph, const Snarl* in_snarl);
        
        /// Returns a map from all Snarl boundaries to the Snarl they point into. Note that this means that
        /// end boundaries will be reversed.
        unordered_map<pair<int64_t, bool>, const Snarl*> snarl_boundary_index();
        
        /// Returns a map from all Snarl start boundaries to the "Snarl they point into.
        unordered_map<pair<int64_t, bool>, const Snarl*> snarl_start_index();
        
        /// Returns a map from all Snarl end boundaries to the Snarl they point into. Note that this means that
        /// end boundaries will be reversed.
        unordered_map<pair<int64_t, bool>, const Snarl*> snarl_end_index();
        
        /// Execute a function on all top level sites
        void for_each_top_level_snarl(const function<void(const Snarl*)>& lambda);
        
        /// Execute a function on all top level sites in parallel
        void for_each_top_level_snarl_parallel(const function<void(const Snarl*)>& lambda);
        
        /// Given a Snarl that we don't own (like from a Visit), find the
        /// pointer to the managed copy of that Snarl.
        const Snarl* manage(const Snarl& not_owned);
        
    private:
    
        /// Define the key type
        using key_t = pair<pair<int64_t, bool>, pair<int64_t, bool>>;
        
        /// Master list of the snarls in the graph
        vector<Snarl> snarls;
        
        /// Roots of snarl trees
        vector<const Snarl*> roots;
        
        /// Map of snarls to the child snarls they contain
        unordered_map<key_t, vector<const Snarl*>> children;
        /// Map of snarls to their parent
        unordered_map<key_t, const Snarl*> parent;
        
        /// Map of snarl keys to the indexes in the snarl array
        // TODO: should we switch to just pointers here and save an indirection?
        unordered_map<key_t, size_t> index_of;
        
        /// Map of node traversals to the snarls they point into
        unordered_map<pair<int64_t, bool>, const Snarl*> snarl_into;
        
        /// Converts Snarl to the form used as keys in internal data structures
        inline key_t key_form(const Snarl* snarl);
        
        /// Builds tree indexes after Snarls have been added
        void build_indexes();
    };
    
    /// Converts a Visit to a NodeTraversal. Throws an exception if the Visit is of a Snarl instead
    /// of a Node
    inline NodeTraversal to_node_traversal(const Visit& visit, const VG& graph);
    
    /// Converts a Visit to a NodeTraversal in the opposite orientation. Throws an exception if the
    /// Visit is of a Snarl instead of a Node
    inline NodeTraversal to_rev_node_traversal(const Visit& visit, const VG& graph);
    
    /// Converts a Visit to a node or snarl into a NodeSide for its left side.
    inline NodeSide to_left_side(const Visit& visit);
    
    /// Converts a Visit to a node or snarl into a NodeSide for its right side.
    inline NodeSide to_right_side(const Visit& visit);
    
    /// Converts a NodeTraversal to a Visit.
    inline Visit to_visit(const NodeTraversal& node_traversal);
    
    /// Converts a Mapping to a Visit. The mapping must represent a full node match.
    inline Visit to_visit(const Mapping& mapping);
    
    /// Make a Visit from a node ID and an orientation
    inline Visit to_visit(id_t node_id, bool is_reverse);
    
    /// Make a Visit from a snarl to traverse
    inline Visit to_visit(const Snarl& snarl);
    
    /// Get the reversed version of a visit
    inline Visit reverse(const Visit& visit);
    
    /// Converts a NodeTraversal to a Visit in the opposite orientation.
    inline Visit to_rev_visit(const NodeTraversal& node_traversal);
    
    /// Converts a Visit to a Mapping. Throws an exception if the Visit is of a Snarl instead
    /// of a Node. Uses a function to get node length.
    inline Mapping to_mapping(const Visit& visit, std::function<size_t(id_t)> node_length);
    
    /// Converts a Visit to a Mapping. Throws an exception if the Visit is of a Snarl instead
    /// of a Node. Uses a graph to get node length.
    inline Mapping to_mapping(const Visit& visit, VG& vg);
    
    /// Copies the boundary Visits from one Snarl into another
    inline void transfer_boundary_info(const Snarl& from, Snarl& to);
    
    // We need some Visit operators
    
    /**
     * Two Visits are equal if they represent the same traversal of the same
     * Node or Snarl.
     */
    bool operator==(const Visit& a, const Visit& b);
    /**
     * Two Visits are unequal if they are not equal.
     */
    bool operator!=(const Visit& a, const Visit& b);
    /**
     * A Visit is less than another Visit if it represents a traversal of a
     * smaller node, or it represents a traversal of a smaller snarl, or it
     * represents a traversal of the same node or snarl forward instead of
     * backward.
     */
    bool operator<(const Visit& a, const Visit& b);
    
    /**
     * A Visit can be printed.
     */
    ostream& operator<<(ostream& out, const Visit& visit);
    
    // And some operators for SnarlTraversals
    
    /**
     * Two SnarlTraversals are equal if their snarls are equal and they have the
     * same number of visits and all their visits are equal.
     */
    bool operator==(const SnarlTraversal& a, const SnarlTraversal& b);
    /**
     * Two SnarlTraversals are unequal if they are not equal.
     */
    bool operator!=(const SnarlTraversal& a, const SnarlTraversal& b);
    /**
     * A SnalTraversal is less than another if it is a traversal of a smaller
     * Snarl, or if its list of Visits has a smaller Visit first, or if its list
     * of Visits is shorter.
     */
    bool operator<(const SnarlTraversal& a, const SnarlTraversal& b);
    
    // And some operators for Snarls
    
    /**
     * Two Snarls are equal if their types are equal and their bounding Visits
     * are equal and their parents are equal.
     */
    bool operator==(const Snarl& a, const Snarl& b);
    /**
     * Two Snarls are unequal if they are not equal.
     */
    bool operator!=(const Snarl& a, const Snarl& b);
    /**
     * A Snarl is less than another Snarl if its type is smaller, or its start
     * Visit is smaller, or its end Visit is smaller, or its parent is smaller.
     */
    bool operator<(const Snarl& a, const Snarl& b);
    
    /**
     * A Snarl can be printed.
     */
    ostream& operator<<(ostream& out, const Snarl& snarl);
    
    /****
     * Template and Inlines:
     ****/
    
    template <typename SnarlIterator>
    SnarlManager::SnarlManager(SnarlIterator begin, SnarlIterator end) {
        // add snarls to master list
        for (auto iter = begin; iter != end; iter++) {
            snarls.push_back(*iter);
        }
        // record the tree structure and build the other indexes
        build_indexes();
    }
    
    inline NodeTraversal to_node_traversal(const Visit& visit, VG& graph) {
        assert(visit.node_id());
        return NodeTraversal(graph.get_node(visit.node_id()), visit.backward());
    }
    
    inline NodeTraversal to_rev_node_traversal(const Visit& visit, VG& graph) {
        assert(visit.node_id());
        return NodeTraversal(graph.get_node(visit.node_id()), !visit.backward());
    }
    
    inline NodeSide to_left_side(const Visit& visit) {
        assert(visit.node_id() || (visit.snarl().start().node_id() && visit.snarl().end().node_id()));
        if (visit.node_id()) {
            // Just report the left side of this node
            return NodeSide(visit.node_id(), visit.backward());
        } else if (visit.backward()) {
            // This is a reverse visit to a snarl, so its left side is the right
            // side of the end visit of the snarl.
            assert(visit.snarl().end().node_id());
            return to_right_side(visit.snarl().end());
        } else {
            // This is a forward visit to a snarl, so its left side is the left
            // side of the start visit of the snarl.
            assert(visit.snarl().start().node_id());
            return to_left_side(visit.snarl().start());
        }
    }
    
    inline NodeSide to_right_side(const Visit& visit) {
        assert(visit.node_id() || (visit.snarl().start().node_id() && visit.snarl().end().node_id()));
        if (visit.node_id()) {
            // Just report the right side of this node
            return NodeSide(visit.node_id(), !visit.backward());
        } else if (visit.backward()) {
            // This is a reverse visit to a snarl, so its right side is the
            // left side of the start visit of the snarl.
            assert(visit.snarl().start().node_id());
            return to_left_side(visit.snarl().start());
        } else {
            // This is a forward visit to a snarl, so its right side is the
            // right side of the end visit of the snarl.
            assert(visit.snarl().end().node_id());
            return to_right_side(visit.snarl().end());
        }
    }
    
    inline Visit to_visit(const NodeTraversal& node_traversal) {
        Visit to_return;
        to_return.set_node_id(node_traversal.node->id());
        to_return.set_backward(node_traversal.backward);
        return to_return;
    }
    
    inline Visit to_visit(const Mapping& mapping) {
        assert(mapping_is_match(mapping));
        assert(mapping.position().offset() == 0);
        Visit to_return;
        to_return.set_node_id(mapping.position().node_id());
        to_return.set_backward(mapping.position().is_reverse());
        return to_return;
    }
    
    inline Visit to_visit(id_t node_id, bool is_reverse) {
        Visit to_return;
        to_return.set_node_id(node_id);
        to_return.set_backward(is_reverse);
        return to_return;
    }
    
    inline Visit to_visit(const Snarl& snarl) {
        Visit to_return;
        // Only copy necessary fields
        *to_return.mutable_snarl()->mutable_start() = snarl.start();
        *to_return.mutable_snarl()->mutable_end() = snarl.end();
        return to_return;
    }
    
    inline Visit reverse(const Visit& visit) {
        // Copy the visit
        Visit to_return = visit;
        // And flip its orientation bit
        to_return.set_backward(!visit.backward());
        return to_return;
    }
    
    inline Visit to_rev_visit(const NodeTraversal& node_traversal) {
        Visit to_return;
        to_return.set_node_id(node_traversal.node->id());
        to_return.set_backward(!node_traversal.backward);
        return to_return;
    }
    
    inline Mapping to_mapping(const Visit& visit, std::function<size_t(id_t)> node_length) {
        // Can't have a Mapping to a snarl
        assert(visit.node_id());
    
        // Make a mapping to the right place
        Mapping mapping;
        mapping.mutable_position()->set_node_id(visit.node_id());
        mapping.mutable_position()->set_is_reverse(visit.backward());
        
        // Get the length of the node visited
        size_t length = node_length(visit.node_id());
        
        // Fill the Mapping in as a perfect match of that lenght
        Edit* match = mapping.add_edit();
        match->set_from_length(length);
        match->set_to_length(length);
        
        return mapping;
    }
    
    inline Mapping to_mapping(const Visit& visit, VG& graph) {
        return to_mapping(visit, [&](id_t id) {
            return graph.get_node(id)->sequence().size();
        });
    }
    
    inline void transfer_boundary_info(const Snarl& from, Snarl& to) {
        *to.mutable_start() = from.start();
        *to.mutable_end() = from.end();
    }
    
}

// note: this hash funtion is not used internally because we want the internal indices to ignore any
// additional information so that the Snarls stored as references map to the same place
// as the original objects
namespace std {
    /// hash function for Snarls
    template<>
    struct hash<const vg::Snarl> {
        size_t operator()(const vg::Snarl& snarl) const {
            auto hsh = hash<pair<pair<int64_t, bool>, pair<int64_t, bool> > >();
            return hsh(make_pair(make_pair(snarl.start().node_id(), snarl.start().backward()),
                                 make_pair(snarl.end().node_id(), snarl.end().backward())));
        }
    };
}

#endif /* snarls_hpp */

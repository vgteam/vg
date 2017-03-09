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
        
        /// Default constructor
        SnarlManager() = default;
        
        /// Destructor
        ~SnarlManager() = default;
        
        /// Returns a vector of pointers to the children of a Snarl
        const vector<const Snarl*>& children_of(const Snarl* snarl);
        
        /// Returns a pointer to the parent of a Snarl or nullptr if there is none
        const Snarl* parent_of(const Snarl* snarl);
        
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
        
        // Returns a map from the boundaries of the child Snarls to the respective child Snarl. End NodeTraversals
        // are reversed to point into the Snarl.
        map<NodeTraversal, const Snarl*> child_boundary_index(const Snarl* snarl, VG& graph);
        
        // Returns a map from the start boundary of the children to the respective child Snarl.
        map<NodeTraversal, const Snarl*> child_start_index(const Snarl* snarl, VG& graph);
        
        // Returns a map from the reversed end boundary of the children to the respective child Snarl.
        map<NodeTraversal, const Snarl*> child_end_index(const Snarl* snarl, VG& graph);
        
        /// Execute a function on all top level sites
        void for_each_top_level_snarl(const function<void(const Snarl*)>& lambda);
        
        /// Execute a function on all top level sites in parallel
        void for_each_top_level_snarl_parallel(const function<void(const Snarl*)>& lambda);
        
    private:
        
        /// Master list of the snarls in the graph
        vector<Snarl> snarls;
        
        /// Roots of snarl trees
        vector<const Snarl*> roots;
        
        /// Map of snarls to the child snarls they contain
        unordered_map<pair<pair<int64_t, bool>, pair<int64_t, bool> >, vector<const Snarl*> > children;
        unordered_map<pair<pair<int64_t, bool>, pair<int64_t, bool> >, const Snarl*> parent;
        
        /// Converts Snarl to the form used as keys in internal data structures
        inline pair<pair<int64_t, bool>, pair<int64_t, bool> > key_form(const Snarl* snarl);
        
        /// Builds tree indices after Snarls have been added
        void build_trees();
    };
    
    /// Converts a Visit to a NodeTraversal. Throws an exception if the Visit is of a Snarl instead
    /// of a Node
    inline NodeTraversal to_node_traversal(const Visit& visit, const VG& graph);
    
    /// Converts a Visit to a NodeTraversal in the opposite orientation. Throws an exception if the
    /// Visit is of a Snarl instead of a Node
    inline NodeTraversal to_rev_node_traversal(const Visit& visit, const VG& graph);
    
    /// Converts a Visit to a NodeSide for its left side. Throws an exception if
    /// the Visit is of a Snarl instead of a Node.
    inline NodeSide to_left_side(const Visit& visit);
    
    /// Converts a Visit to a NodeSide for its right side. Throws an exception if
    /// the Visit is of a Snarl instead of a Node.
    inline NodeSide to_right_side(const Visit& visit);
    
    /// Converts a NodeTraversal to a Visit.
    inline Visit to_visit(const NodeTraversal& node_traversal);
    
    /// Converts a Mapping to a Visit. The mapping must represent a full node match.
    inline Visit to_visit(const Mapping& mapping);
    
    /// Converts a NodeTraversal to a Visit in the opposite orientation.
    inline Visit to_rev_visit(const NodeTraversal& node_traversal);
    
    /// Converts a Visit to a Mapping. Throws an exception if the Visit is of a Snarl instead
    /// of a Node. Uses a function to get node length.
    inline Mapping to_mapping(const Visit& visit, std::function<size_t(id_t)> node_length);
    
    /// Converts a Visit to a Mapping. Throws an exception if the Visit is of a Snarl instead
    /// of a Node. Uses a graph to get node length.
    inline Mapping to_mapping(const Visit& visit, VG& vg);
    
    // Copies the boundary Visits from one Snarl into another
    inline void transfer_boundary_info(const Snarl& from, Snarl& to);
    
    /****
     * Template and Inlines:
     ****/
    
    template <typename SnarlIterator>
    SnarlManager::SnarlManager(SnarlIterator begin, SnarlIterator end) {
        // add snarls to master list
        for (auto iter = begin; iter != end; iter++) {
            snarls.push_back(*iter);
        }
        // record the tree structure
        build_trees();
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
        assert(visit.node_id());
        return NodeSide(visit.node_id(), visit.backward());
    }
    
    inline NodeSide to_right_side(const Visit& visit) {
        assert(visit.node_id());
        return NodeSide(visit.node_id(), !visit.backward());
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

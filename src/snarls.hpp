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
        
        /// Returns the Snarl that a traversal points into at either the start or end, or nullptr if
        /// the traversal does not point into any Snarl. Note that Snarls store the end Visit pointing
        /// out of rather than into the Anarl, so they must be reversed to query it.
        const Snarl* into_which_snarl(int64_t id, bool reverse);
        
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
        
        /// Map of snarls to their parent
        unordered_map<pair<pair<int64_t, bool>, pair<int64_t, bool> >, const Snarl*> parent;
        
        /// Map of node traversals to the snarls they point into
        unordered_map<pair<int64_t, bool>, const Snarl*> snarl_into;
        
        /// Converts Snarl to the form used as keys in internal data structures
        inline pair<pair<int64_t, bool>, pair<int64_t, bool> > key_form(const Snarl* snarl);
        
        /// Builds tree indices after Snarls have been added
        void build_indices();
    };
    
    /// Converts a Visit to a NodeTraversal. Throws an exception if the Visit is of a Snarl instead
    /// of a Node
    inline NodeTraversal to_node_traversal(const Visit& visit, const VG& graph);
    
    /// Converts a Visit to a NodeTraversal in the opposite orientation. Throws an exception if the
    /// Visit is of a Snarl instead of a Node
    inline NodeTraversal to_rev_node_traversal(const Visit& visit, const VG& graph);
    
    /// Converts a NodeTraversal to a Visit.
    inline Visit to_visit(const NodeTraversal& node_traversal);
    
    /// Converts a NodeTraversal to a Visit in the opposite orientation.
    inline Visit to_rev_visit(const NodeTraversal& node_traversal);
    
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
        build_indices();
    }
    
    inline NodeTraversal to_node_traversal(const Visit& visit, VG& graph) {
        assert(visit.node_id());
        return NodeTraversal(graph.get_node(visit.node_id()), visit.backward());
    }
    
    inline NodeTraversal to_rev_node_traversal(const Visit& visit, VG& graph) {
        assert(visit.node_id());
        return NodeTraversal(graph.get_node(visit.node_id()), !visit.backward());
    }
    
    inline Visit to_visit(const NodeTraversal& node_traversal) {
        Visit to_return;
        to_return.set_node_id(node_traversal.node->id());
        to_return.set_backward(node_traversal.backward);
        return to_return;
    }
    
    inline Visit to_rev_visit(const NodeTraversal& node_traversal) {
        Visit to_return;
        to_return.set_node_id(node_traversal.node->id());
        to_return.set_backward(!node_traversal.backward);
        return to_return;
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

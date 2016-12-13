//
//  snarls.hpp
//
//  Contains object to own Snarls and keep track of their tree relationships as well as utility
//  functions that interact with snarls.
//

#ifndef snarls_hpp
#define snarls_hpp

#include <stdio.h>
#include <unordered_map>
#include "vg.hpp"
#include "vg.pb.h"
#include "hash_map.hpp"

using namespace std;

namespace vg {
    
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
        const vector<const Snarl*>& children_of(const Snarl& snarl);
        
        /// Returns a pointer to the parent of a Snarl or nullptr if there is none
        const Snarl* parent_of(const Snarl& snarl);
        
        /// Returns the Nodes and Edges contained in this Snarl but not in any child Snarls (includes the
        /// Nodes that form the boundaries of child Snarls, does not include Snarl's own boundary Nodes)
        pair<vector<Node*>, vector<Edge*> > shallow_contents(const Snarl&, VG& graph);
        
        /// Returns the Nodes and Edges contained in this Snarl, including those in child Snarls (does not
        /// include Snarl's own boundary Nodes)
        pair<vector<Node*>, vector<Edge*> > deep_contents(const Snarl&, VG& graph);
        
        /// Returns a reference to a vector with the roots of the Snarl trees
        const vector<const Snarl*>& top_level_snarls();
        
        /// Execute a function on all top level sites in parallel
        void for_each_top_level_snarl_parallel(const function<void(const Snarl&)>& lambda);
        
    private:
        
        /// Master list of the snarls in the graph
        vector<Snarl> snarls;
        
        /// Roots of snarl trees
        vector<const Snarl*> roots;
        
        /// Map of snarls to the child snarls they contain
        unordered_map<pair<pair<int64_t, bool>, pair<int64_t, bool> >, vector<const Snarl*> > children;
        unordered_map<pair<pair<int64_t, bool>, pair<int64_t, bool> >, const Snarl*> parent;
        
        /// Converts Snarl to the form used as keys in internal data structures
        inline pair<pair<int64_t, bool>, pair<int64_t, bool> > key_form(const Snarl& snarl);
        
        /// Builds tree indices after Snarls have been added
        void build_trees();
        
    };
    
    /// Converts a Visit to a NodeTraversal. Throws an exception if the Visit is of a Snarl instead
    /// of a Node
    inline NodeTraversal to_node_traversal(const Visit& visit, const VG& graph);
    
    /// Converts a Visit to a NodeTraversal in the opposite orientation. Throws an exception if the
    /// Visit is of a Snarl instead of a Node
    inline NodeTraversal to_rev_node_traversal(const Visit& visit, const VG& graph);
    
    // must place template and inline function in header
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
}

#endif /* snarls_hpp */

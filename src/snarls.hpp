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
#include "handle.hpp"
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
    
    /**
     * Allow traversing a graph of nodes and child snarl chains within a snarl
     * within another HandleGraph. Uses its own internal child index because
     * it's used in the construction of snarls to feed to SnarlManagers.
     *
     * Assumes that the chains we get from Cactus are in a consistent order, so
     * the start of the first snarl is the very first thing in the chain, and
     * the end of the last snarl is the very last.
     *
     * We adapt the handle graph abstraction as follows:
     *
     * A chain becomes a single node with the ID and local forward orientation
     * of its first snarl's start.
     *
     * A chain node connects on its left to everything connected to its first
     * start and on its right to everything connected to its last end.
     *
     * A unary snarl becomes a single node, too. It is identified by its
     * boundary node's ID.
     *
     * If you're not using internal connectivity, a chain node or a unary snarl
     * node behaves just like an ordinary node.
     *
     * If you are using internal connectivity, edges are slightly faked:
     *
     * A chain node also sees out its right everything that is out its left if
     * it has a left-left connected snarl before any disconnected snarl.
     *
     * And similarly for the mirror case.
     *
     * All the edges on either side of a unary snarl node are the same.
     *
     * In this part of the code we talk about "heads" (the inward-facing base
     * graph handles used to represent child snarls/chains), and "tails" (the
     * inward-facing ending handles of child chains).
     * 
     */
    class NetGraph : public HandleGraph {
    public:
        /// Make a new NetGraph from the given chains and unary snarls in the given backing graph.
        NetGraph(const Visit& start, const Visit& end, const vector<vector<Snarl>>& child_chains,
            const vector<Snarl>& child_unary_snarls, const HandleGraph* graph, bool use_internal_connectivity = false);
    
        /// Look up the handle for the node with the given ID in the given orientation
        virtual handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
        
        /// Get the ID from a handle
        virtual id_t get_id(const handle_t& handle) const;
        
        /// Get the orientation of a handle
        virtual bool get_is_reverse(const handle_t& handle) const;
        
        /// Invert the orientation of a handle (potentially without getting its ID)
        virtual handle_t flip(const handle_t& handle) const;
        
        /// Get the length of a node
        virtual size_t get_length(const handle_t& handle) const;
        
        /// Get the sequence of a node, presented in the handle's local forward
        /// orientation.
        virtual string get_sequence(const handle_t& handle) const;
        
        /// Loop over all the handles to next/previous (right/left) nodes. Passes
        /// them to a callback which returns false to stop iterating and true to
        /// continue. Returns true if we finished and false if we stopped early.
        virtual bool follow_edges(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
        
        // Copy over the template for nice calls
        using HandleGraph::follow_edges;
        
        /// Loop over all the nodes in the graph in their local forward
        /// orientations, in their internal stored order. Stop if the iteratee returns false.
        virtual void for_each_handle(const function<bool(const handle_t&)>& iteratee) const;
        
        // Copy over the template for nice calls
        using HandleGraph::for_each_handle;
        
        /// Return the number of nodes in the graph
        /// TODO: can't be node_count because XG has a field named node_count.
        virtual size_t node_size() const;
        
    protected:
    
        // Save the backing graph
        const HandleGraph* graph;
        
        // And the start and end handles that bound the snarl we are working on.
        handle_t start;
        handle_t end;
        
        // Should we use the internal connectivity of chain nodes and unary
        // snarl nodes?
        bool use_internal_connectivity;
    
        // We keep the unary snarl boundaries, reading in with the unary snarl
        // contents to the right.
        unordered_set<handle_t> unary_boundaries;
        
        // We keep a map from handles that enter the ends of chains to the
        // reverse handles to their fronts. Whenever the backing graph tells us
        // to emit the one, we emit the other instead. This makes them look like
        // one big node.
        unordered_map<handle_t, handle_t> chain_end_rewrites;
        
        // We keep basically the reverse map, from chain start in chain forward
        // orientation to chain end in chain forward orientation. This lets us
        // find the edges off the far end of a chian.
        unordered_map<handle_t, handle_t> chain_ends_by_start;
        
        // Stores whether a chain or unary snarl, identified by the ID of its
        // start handle, is left-left, right-right, or left-right connected.
        unordered_map<id_t, tuple<bool, bool, bool>> connectivity;
        
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

///
///  \file snarls.hpp
///
///  Contains object to own Snarls and keep track of their tree relationships as well as utility
///  functions that interact with snarls.
///

#ifndef VG_SNARLS_HPP_INCLUDED
#define VG_SNARLS_HPP_INCLUDED

#include <iostream>
#include <cstdint>
#include <stdio.h>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <deque>
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/protobuf_iterator.hpp>
#include "vg.hpp"
#include "handle.hpp"
#include <vg/vg.pb.h>
#include "hash_map.hpp"
#include "cactus.hpp"

using namespace std;

namespace vg {

class SnarlManager;

/**
 * Represents a strategy for finding (nested) sites in a vg graph that can be described
 * by snarls. Polymorphic base class/interface.
 */
class SnarlFinder {
public:
    virtual ~SnarlFinder() = default;

    /**
     * Find all the snarls, and put them into a SnarlManager.
     */
    virtual SnarlManager find_snarls() = 0;

    /**
     * Find all the snarls of weakly connected components, optionally in
     * parallel. If not implemented, defaults to the single-threaded
     * implementation.
     */
    virtual SnarlManager find_snarls_parallel();
};

/**
 * Wrapper base class that can convert a bottom-up traversal of snarl
 * boundaries into a full snarl finder. Mostly worries about snarl
 * classification and connectivity information.
 */
class HandleGraphSnarlFinder : public SnarlFinder {
protected:
    /**
     * The graph we are finding snarls on. It must outlive us.
     */
    const HandleGraph* graph;
    
    /**
     * Find all the snarls, and put them into a SnarlManager, but don't finish it.
     * More snarls can be added later before it is finished.
     */
    virtual SnarlManager find_snarls_unindexed();
    
public:

    /**
     * Create a HandleGraphSnarlFinder to find snarls in the given graph.
     */
    HandleGraphSnarlFinder(const HandleGraph* graph);

    virtual ~HandleGraphSnarlFinder() = default;

    /**
     * Find all the snarls, and put them into a SnarlManager.
     */
    virtual SnarlManager find_snarls();
    
    /**
     * Visit all snarls and chains, including trivial snarls and single-node
     * empty chains.
     *
     * Calls begin_chain and end_chain when entrering and exiting chains in the
     * traversal. Within each chain, calls begin_snarl and end_snarl when
     * entering and exiting each snarl, in order. The caller is intended to
     * maintain its own stack to match up begin and end events.
     *
     * Each begin/end call receives the handle reading into/out of the snarl or
     * chain. 
     *
     * Both empty and cyclic chains have the in and out handles the same.
     * They are distinguished by context; empty chains have no shild snarls,
     * while cyclic chains do.
     *
     * Roots the decomposition at a global snarl with no bounding nodes, for
     * which begin_snarl is not called. So the first call will be begin_chain.
     *
     * Start handles are inward facing and end handles are outward facing.
     * Snarls must be oriented forward in their chains.
     */
    virtual void traverse_decomposition(const function<void(handle_t)>& begin_chain, const function<void(handle_t)>& end_chain,
        const function<void(handle_t)>& begin_snarl, const function<void(handle_t)>& end_snarl) const = 0;
};

/**
 * Snarls are defined at the Protobuf level, but here is how we define
 * chains as real objects.
 *
 * A chain is a sequence of Snarls, in either normal (false) or reverse (true)
 * orientation.
 *
 * The SnarlManager is going to have one official copy of each chain stored,
 * and it will give you a pointer to it on demand.
 */
using Chain = vector<pair<const Snarl*, bool>>;
    
/**
 * Return true if the first snarl in the given chain is backward relative to the chain.
 */
bool start_backward(const Chain& chain);
    
/**
 * Return true if the last snarl in the given chain is backward relative to the chain.
 */
bool end_backward(const Chain& chain);
    
/**
 * Get the inward-facing start Visit for a chain.
 */
Visit get_start_of(const Chain& chain);
    
/**
 * Get the outward-facing end Visit for a chain.
 */
Visit get_end_of(const Chain& chain);
    
/**
 * We want to be able to loop over a chain and get iterators to pairs of the
 * snarl and its orientation in the chain. So we define some iterators.
 */
struct ChainIterator {
    /// Advance the iterator
    ChainIterator& operator++();
    /// Get the snarl we're at and whether it is backward 
    pair<const Snarl*, bool> operator*() const;
    /// Get a pointer to the thing we get when we dereference the iterator
    const pair<const Snarl*, bool>* operator->() const;
        
    /// We need to define comparison because C++ doesn't give it to us for free.
    bool operator==(const ChainIterator& other) const;
    bool operator!=(const ChainIterator& other) const;
        
    /// Are we a reverse iterator or not?
    bool go_left;
    
    /// What position in the underlying vector are we in?
    Chain::const_iterator pos;
        
    /// What are the bounds of that underlying vector?
    Chain::const_iterator chain_start;
    Chain::const_iterator chain_end;
        
    /// Since we're using backing random access itarators to provide reverse
    /// iterators, we need a flag to see if we are rend (i.e. before the
    /// beginning)
    bool is_rend;
        
    /// When dereferencing, should we flip snarl orientations form the
    /// orientations they appear at in the chain when read left to right?
    bool complement;
        
    /// In order to dereference to a pair with -> we need a place to put the pair so we can have a pointer to it.
    /// Gets lazily set to wherever the iterator is pointing when we do ->
    mutable pair<const Snarl*, bool> scratch;
};
    
/**
 * We define free functions for getting iterators forward and backward through chains.
 */
ChainIterator chain_begin(const Chain& chain);
ChainIterator chain_end(const Chain& chain);
ChainIterator chain_rbegin(const Chain& chain);
ChainIterator chain_rend(const Chain& chain);
    
/// We also define some reverse complement iterators, which go from right to
/// left through the chains, but give us the reverse view. For ecample, if
/// all the snarls are oriented forward in the chain, we will iterate
/// through the snarls in reverse order, with each individual snarl also
/// reversed.
ChainIterator chain_rcbegin(const Chain& chain);
ChainIterator chain_rcend(const Chain& chain);
    
/// We also define a function for getting the ChainIterator (forward or reverse
/// complement) for a chain starting with a given snarl in the given inward
/// orientation. Only works for bounding snarls of the chain.
ChainIterator chain_begin_from(const Chain& chain, const Snarl* start_snarl, bool snarl_orientation);
/// And the end iterator for the chain (forward or reverse complement) viewed
/// from a given snarl in the given inward orientation. Only works for bounding
/// snarls of the chain, and should be the *same* bounding snarl as was used
/// for chain_begin_from.
ChainIterator chain_end_from(const Chain& chain, const Snarl* start_snarl, bool snarl_orientation);
    
/**
 * Allow traversing a graph of nodes and child snarl chains within a snarl
 * within another HandleGraph. Uses its own internal child index because
 * it's used in the construction of snarls to feed to SnarlManagers.
 *
 * Assumes that the snarls in the chains we get are in the order they
 * occur in the graph.
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
        
    /// Make a new NetGraph for the given snarl in the given backing graph,
    /// using the given chains as child chains. Unary snarls are stored as
    /// trivial chains just like other trivial chains.
    template<typename ChainContainer>
    NetGraph(const Visit& start, const Visit& end,
             const ChainContainer& child_chains_mixed,
             const HandleGraph* graph,
             bool use_internal_connectivity = false) : NetGraph(start, end, graph, use_internal_connectivity) {
            
        // All we need to do is index the children. They come mixed as real chains and unary snarls.
            
        for (auto& chain : child_chains_mixed) {
            if (chain.size() == 1 && chain.front().first->type() == UNARY) {
                // This is a unary snarl wrapped in a chain
                add_unary_child(chain.front().first);
            } else {
                // This is a real (but possibly singlr-snarl) chain
                add_chain_child(chain);
            }
        }
            
    }
        
    /// Make a net graph from the given chains and unary snarls (as pointers) in the given backing graph.
    template<typename ChainContainer, typename SnarlContainer>
    NetGraph(const Visit& start, const Visit& end,
             const ChainContainer& child_chains,
             const SnarlContainer& child_unary_snarls, const HandleGraph* graph,
             bool use_internal_connectivity = false) : NetGraph(start, end, graph, use_internal_connectivity) {
            
        // All we need to do is index the children.
        for (const Snarl* unary : child_unary_snarls) {
            add_unary_child(unary);
        }
            
        for (auto& chain : child_chains) {
            add_chain_child(chain);
        }
    }
            
    /// Make a net graph from the given chains and unary snarls (as raw values) in the given backing graph.
    /// Mostly for testing.
    NetGraph(const Visit& start, const Visit& end,
             const vector<vector<pair<Snarl, bool>>>& child_chains,
             const vector<Snarl>& child_unary_snarls,
             const HandleGraph* graph,
             bool use_internal_connectivity = false);

    /// Method to check if a node exists by ID
    virtual bool has_node(id_t node_id) const;
    
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
    virtual bool follow_edges_impl(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
        
    /// Loop over all the nodes in the graph in their local forward
    /// orientations, in their internal stored order. Stop if the iteratee returns false.
    virtual bool for_each_handle_impl(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
        
    /// Return the number of nodes in the graph
    virtual size_t get_node_count() const;
    
    /// Return the smallest ID used. 
    virtual id_t min_node_id() const;
    
    /// Return the largest ID used.
    virtual id_t max_node_id() const;
        
    // We also have some extra functions
        
    /// Get the inward-facing start handle for this net graph. Useful when
    /// working with traversals.
    const handle_t& get_start() const;
        
    /// Get the outward-facing end handle for this net graph. Useful when
    /// working with traversals.
    const handle_t& get_end() const;
        
    /// Returns true if the given handle represents a meta-node for a child
    /// chain or unary snarl, and false if it is a normal node actually in
    /// the net graph snarl's contents.
    bool is_child(const handle_t& handle) const;
    
    /// Get the handle in the backing graph reading into the child chain or
    /// unary snarl in the orientation represented by this handle to a node
    /// representing a child chain or unary snarl.
    handle_t get_inward_backing_handle(const handle_t& child_handle) const;
    
    /// Given a handle to a node in the backing graph that reads into a child
    /// chain or snarl (in either direction), get the handle in this graph used
    /// to represent that child chain or snarl in that orientation.
    handle_t get_handle_from_inward_backing_handle(const handle_t& backing_handle) const;
        
protected:
    
    /// Make a NetGraph without filling in any of the child indexes.
    NetGraph(const Visit& start, const Visit& end, const HandleGraph* graph, bool use_internal_connectivity = false);
    
    /// Add a unary child snarl to the indexes.
    void add_unary_child(const Snarl* unary);
        
    /// Add a chain of one or more non-unary snarls to the index.
    void add_chain_child(const Chain& chain);
    
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
    // find the edges off the far end of a chain.
    unordered_map<handle_t, handle_t> chain_ends_by_start;
        
    // Stores whether a chain or unary snarl, identified by the ID of its
    // start handle, is left-left, right-right, or left-right connected.
    unordered_map<id_t, tuple<bool, bool, bool>> connectivity;
        
};
    
/**
 * A structure to keep track of the tree relationships between Snarls and perform utility algorithms
 * on them
 */
class SnarlManager {
public:
        
    /// Construct a SnarlManager for the snarls returned by an iterator
    /// Also covers iterators of chains of snarls.
    template <typename SnarlIterator>
    SnarlManager(SnarlIterator begin, SnarlIterator end);
        
    /// Construct a SnarlManager for the snarls contained in an input stream
    SnarlManager(istream& in);
    
    /// Construct a SnarlManager from a function that calls a callback with each Snarl in turn
    SnarlManager(const function<void(const function<void(Snarl&)>&)>& for_each_snarl);
        
    /// Default constructor for an empty SnarlManager. Must call finish() once
    /// all snarls have been added with add_snarl().
    SnarlManager() = default;
        
    /// Destructor
    ~SnarlManager() = default;
        
    /// Cannot be copied because of all the internal pointer indexes
    SnarlManager(const SnarlManager& other) = delete;
    SnarlManager& operator=(const SnarlManager& other) = delete;
        // copy the SnarlManager
    /// Can be moved
    SnarlManager(SnarlManager&& other) = default;
    SnarlManager& operator=(SnarlManager&& other) = default;

    // Can be serialized
    void serialize(ostream& out) const;
    
    ///////////////////////////////////////////////////////////////////////////
    // Write API
    ///////////////////////////////////////////////////////////////////////////
    
    /// Add the given snarl to the SnarlManager. After all snarls have been
    /// added, finish() must be called to compute chains and indexes. We don't
    /// let precomputed chains be added, because we want chain orientations
    /// relative to snarls to be deterministic given an order of snarls.
    /// Returns a pointer to the managed snarl copy.
    /// Only this function may add in new Snarls.
    const Snarl* add_snarl(const Snarl& new_snarl);
    
    /// Reverses the orientation of a managed snarl.
    void flip(const Snarl* snarl);
    
    /// Reverses the order and orientation of a managed chain, leaving all the
    /// component snarls in their original orientations.
    void flip(const Chain* snarl);
        
    /// Note that we have finished calling add_snarl. Compute the snarl
    /// parent/child indexes and chains.
    void finish();
    
    ///////////////////////////////////////////////////////////////////////////
    // Read API
    ///////////////////////////////////////////////////////////////////////////

    /// Returns a vector of pointers to the children of a Snarl.
    /// If given null, returns the top-level root snarls.
    const vector<const Snarl*>& children_of(const Snarl* snarl) const;
        
    /// Returns a pointer to the parent of a Snarl or nullptr if there is none
    const Snarl* parent_of(const Snarl* snarl) const;
        
    /// Returns the Snarl that a traversal points into at either the start
    /// or end, or nullptr if the traversal does not point into any Snarl.
    /// Note that Snarls store the end Visit pointing out of rather than
    /// into the Snarl, so they must be reversed to query it.
    const Snarl* into_which_snarl(int64_t id, bool reverse) const;
        
    /// Returns the Snarl that a Visit points into. If the Visit contains a
    /// Snarl rather than a node ID, returns a pointer the managed version
    /// of that snarl.
    const Snarl* into_which_snarl(const Visit& visit) const;
    
    /// Get the Chain that the given snarl participates in. Instead of asking
    /// this class to walk the chain for you, use ChainIterators on this chain.
    /// This is always non-null.
    const Chain* chain_of(const Snarl* snarl) const;
    
    /// If the given Snarl is backward in its chain, return true. Otherwise,
    /// return false.
    bool chain_orientation_of(const Snarl* snarl) const;
    
    /// Get the rank that the given snarl appears in in its chain. If two
    /// snarls are in forward orientation in the chain, then leaving the end of
    /// the lower rank snarl will eventually reach the start of the higher rank
    /// snarl. If either or both snarls is backward, you leave/arrive at the
    /// other bounding node instead.
    ///
    /// Sorting snarls by rank will let you visit them in chain order without
    /// walking the whole chain.
    size_t chain_rank_of(const Snarl* snarl) const;
    
    /// Return true if a Snarl is part of a nontrivial chain of more than one
    /// snarl. Note that chain_of() still works for snarls in trivial chains.
    bool in_nontrivial_chain(const Snarl* here) const;
        
    /// Get all the snarls in all the chains under the given parent snarl.
    /// If the parent snarl is null, gives the top-level chains that connect and contain the top-level root snarls.
    /// Unary snarls and snarls in trivial chains will be presented as their own chains.
    /// Snarls are not necessarily oriented appropriately given their ordering in the chain.
    /// Useful for making a net graph.
    const deque<Chain>& chains_of(const Snarl* snarl) const;
        
    /// Get the net graph of the given Snarl's contents, using the given
    /// backing HandleGraph. If use_internal_connectivity is false, each
    /// chain and unary child snarl is treated as an ordinary node which is
    /// assumed to be only traversable from one side to the other.
    /// Otherwise, traversing the graph works like it would if you actually
    /// went through the internal graphs fo child snarls.
    NetGraph net_graph_of(const Snarl* snarl, const HandleGraph* graph, bool use_internal_connectivity = true) const;
        
    /// Returns true if snarl has no children and false otherwise
    bool is_leaf(const Snarl* snarl) const;
        
    /// Returns true if snarl has no parent and false otherwise
    bool is_root(const Snarl* snarl) const;

    /// Returns true if the snarl is trivial (an ultrabubble with just the
    /// start and end nodes) and false otherwise.
    /// TODO: Implement without needing the vg graph, by adding a flag to trivial snarls.
    bool is_trivial(const Snarl* snarl, const HandleGraph& graph) const;
    
    /// Returns true if the snarl lacks any nontrivial children.
    bool all_children_trivial(const Snarl* snarl, const HandleGraph& graph) const;

    /// Returns a reference to a vector with the roots of the Snarl trees
    const vector<const Snarl*>& top_level_snarls() const;
        
    /// Returns the Nodes and Edges contained in this Snarl but not in any child Snarls (always includes the
    /// Nodes that form the boundaries of child Snarls, optionally includes this Snarl's own boundary Nodes)
    pair<unordered_set<id_t>, unordered_set<edge_t> > shallow_contents(const Snarl* snarl, const HandleGraph& graph,
                                                                       bool include_boundary_nodes) const;
        
    /// Returns the Nodes and Edges contained in this Snarl, including those in child Snarls (optionally
    /// includes Snarl's own boundary Nodes)
    pair<unordered_set<id_t>, unordered_set<edge_t> > deep_contents(const Snarl* snarl, const HandleGraph& graph,
                                                                    bool include_boundary_nodes) const;
        
    /// Look left from the given visit in the given graph and gets all the
    /// attached Visits to nodes or snarls.
    vector<Visit> visits_left(const Visit& visit, const HandleGraph& graph, const Snarl* in_snarl) const;
        
    /// Look left from the given visit in the given graph and gets all the
    /// attached Visits to nodes or snarls.
    vector<Visit> visits_right(const Visit& visit, const HandleGraph& graph, const Snarl* in_snarl) const;
        
    /// Returns a map from all Snarl boundaries to the Snarl they point into. Note that this means that
    /// end boundaries will be reversed.
    unordered_map<pair<int64_t, bool>, const Snarl*> snarl_boundary_index() const;
        
    /// Returns a map from all Snarl start boundaries to the Snarl they point into.
    unordered_map<pair<int64_t, bool>, const Snarl*> snarl_start_index() const;
        
    /// Returns a map from all Snarl end boundaries to the Snarl they point into. Note that this means that
    /// end boundaries will be reversed.
    unordered_map<pair<int64_t, bool>, const Snarl*> snarl_end_index() const;
        
    /// Execute a function on all top level sites
    void for_each_top_level_snarl(const function<void(const Snarl*)>& lambda) const;
        
    /// Execute a function on all sites in a preorder traversal
    void for_each_snarl_preorder(const function<void(const Snarl*)>& lambda) const;
        
    /// Execute a function on all top level sites in parallel
    void for_each_top_level_snarl_parallel(const function<void(const Snarl*)>& lambda) const;
        
    /// Execute a function on all sites in parallel
    void for_each_snarl_parallel(const function<void(const Snarl*)>& lambda) const;

    /// Execute a function on all top level chains
    void for_each_top_level_chain(const function<void(const Chain*)>& lambda) const;

    /// Execute a function on all top level chains in parallel
    void for_each_top_level_chain_parallel(const function<void(const Chain*)>& lambda) const;

    /// Ececute a function on all chains
    void for_each_chain(const function<void(const Chain*)>& lambda) const;
    
    /// Ececute a function on all chains in parallel
    void for_each_chain_parallel(const function<void(const Chain*)>& lambda) const;

    /// Iterate over snarls as they are stored in deque<SnarlRecords>
    void for_each_snarl_unindexed(const function<void(const Snarl*)>& lambda) const;
        
    /// Given a Snarl that we don't own (like from a Visit), find the
    /// pointer to the managed copy of that Snarl.
    const Snarl* manage(const Snarl& not_owned) const;

    /// Sample snarls discrete uniformly 
    /// Returns a nullptr if no snarls are found 
    const Snarl* discrete_uniform_sample(minstd_rand0& random_engine)const;

    /// Count snarls in deque<SnarlRecords>, a master list of snarls in graph
    int num_snarls()const;

    ///Get the snarl number from the SnarlRecord* member with given snarl
    inline size_t snarl_number(const Snarl* snarl) const{
        const SnarlRecord* record = SnarlManager::record(snarl);
        return record->snarl_number;
    }
    //use the snarl number to access the Snarl*
    inline const Snarl* translate_snarl_num(size_t snarl_num){
        return unrecord(&snarls.at(snarl_num));
    }

        
private:
    
    /// To support the Snarl*-driven API, we use a struct that lays out a snarl
    /// followed by indexing metadata, one after the other in memory. We can
    /// just cast a Snarl* to a pointer to one of these to get access to all
    /// the metadata.
    struct alignas(alignof(Snarl)) SnarlRecord {
        /// With recent Protobuf, we can't inherit from Protobuf generated
        /// classes, so we rely on the first member here being at offset 0.
        /// This is achieved by making sure SnarlRecord is aligned like Snarl. 
        Snarl snarl;
        
        /// This is a vector of pointers into the master snarl container at
        /// children. We know the pointers are to valid SnarlRecords. A
        /// SnarlRecord does not own its children.
        vector<const Snarl*> children;
        
        /// This holds chains over the child snarls.
        deque<Chain> child_chains;
        
        /// This points to the parent SnarlRecord (as a snarl), or null if we
        /// are a root snarl or have not been told of our parent yet.
        const Snarl* parent = nullptr;
        
        /// This points to the chain we are in, or null if we are not in a chain.
        Chain* parent_chain = nullptr;
        /// And this is what index we are at in the chain;
        size_t parent_chain_index = 0;

        /// This holds the index of the SnarlRecord* in the deque
        /// We are doing this because a deque is not contiguous and the index lookup using a SnarlRecord* isn't easily derivable 
        size_t snarl_number;

        /// Allow assignment from a Snarl object, fluffing it up into a full SnarlRecord
        SnarlRecord& operator=(const Snarl& other) {
            // Just call the base assignment operator
            (*(Snarl*)this) = other;
            return *this;
        }
    };
    
    /// Get the const SnarlRecord for a const managed snarl
    inline const SnarlRecord* record(const Snarl* snarl) const {
        return (const SnarlRecord*) snarl;
    }
    
    /// Get the SnarlRecord for a managed snarl
    inline SnarlRecord* record(Snarl* snarl) {
        return (SnarlRecord*) snarl;
    }
    
    /// Get the const Snarl owned by a const SnarlRecord
    inline const Snarl* unrecord(const SnarlRecord* record) const {
        return (const Snarl*) record;
    }
    
    /// Get the Snarl owned by a SnarlRecord
    inline Snarl* unrecord(SnarlRecord* record) {
        return (Snarl*) record;
    }
    

    /// Master list of the snarls in the graph.
    /// Use a deque so pointers never get invalidated but we still have some locality.
    deque<SnarlRecord> snarls;
        
    /// Roots of snarl trees
    vector<const Snarl*> roots;
    /// Chains of root-level snarls. Uses a deque so Chain* pointers don't get invalidated.
    deque<Chain> root_chains;
        
    /// Map of node traversals to the snarls they point into
    unordered_map<pair<int64_t, bool>, const Snarl*> snarl_into;
        
    /// Builds tree indexes after Snarls have been added to the snarls vector
    void build_indexes();
        
    /// Actually compute chains for a set of already indexed snarls, which
    /// is important when chains were not provided. Returns the chains.
    deque<Chain> compute_chains(const vector<const Snarl*>& input_snarls);
    
    /// Modify the snarls and chains to enforce a couple of invariants:
    ///
    /// 1. The start node IDs of the snarls in a chain shall be unique.
    ///
    /// (This is needed by the distance indexing code, which identifies child
    /// snarls by their start nodes. TODO: That distance indexing code needs to
    /// also work out unary snarls abitting the ends of chains, which may be
    /// allowed eventually.)
    ///
    /// 2. Snarls will be oriented forward in their chains.
    ///
    /// 3. Snarls will be oriented in a chain to maximize the number of snarls
    /// that start with lower node IDs than they end with.
    ///
    /// Depends on the indexes from build_indexes() having been built.
    void regularize();
        
    // Chain computation uses these pseudo-chain-traversal functions, which
    // walk around based on the snarl boundary index. This basically gets
    // you chains, except for structures that look like circular chains,
    // which are actually presented as circular chains and not linear ones.
    // They also let you walk into unary snarls.
        
    /// Get a Visit to the snarl coming after the given Visit to a snarl, or
    /// a Visit with no Snarl no next snarl exists. Accounts for snarls'
    /// orientations.
    Visit next_snarl(const Visit& here) const;
        
    /// Get a Visit to the snarl coming before the given Visit to a snarl,
    /// or a Visit with no Snarl no previous snarl exists. Accounts for
    /// snarls' orientations.
    Visit prev_snarl(const Visit& here) const;
        
    /// Get the Snarl, if any, that shares this Snarl's start node as either
    /// its start or its end. Does not count this snarl, even if this snarl
    /// is unary. Basic operation used to traverse a chain. Caller must
    /// account for snarls' orientations within a chain.
    const Snarl* snarl_sharing_start(const Snarl* here) const;
        
    /// Get the Snarl, if any, that shares this Snarl's end node as either
    /// its start or its end. Does not count this snarl, even if this snarl
    /// is unary. Basic operation used to traverse a chain. Caller must
    /// account for snarls' orientations within a chain.
    const Snarl* snarl_sharing_end(const Snarl* here) const;
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
    
/// Converts a Mapping to a Visit. The mapping must represent a full node
/// match. If make_full_node_match is true, the mapping will automatically
/// be made a full node match during the conversion process.
inline Visit to_visit(const Mapping& mapping, bool make_full_node_match = false);
    
/// Make a Visit from a node ID and an orientation
inline Visit to_visit(id_t node_id, bool is_reverse);
    
/// Make a Visit from a snarl to traverse
inline Visit to_visit(const Snarl& snarl);

/// Make a Visit from a handle in a HandleGraph.
inline Visit to_visit(const handlegraph::HandleGraph& graph, const handle_t& handle);
    
/// Get the reversed version of a visit
inline Visit reverse(const Visit& visit);
    
/// Converts a NodeTraversal to a Visit in the opposite orientation.
inline Visit to_rev_visit(const NodeTraversal& node_traversal);
    
/// Converts a Visit to a Mapping. Throws an exception if the Visit is of a Snarl instead
/// of a Node. Uses a function to get node length.
inline Mapping to_mapping(const Visit& visit, std::function<size_t(id_t)> node_length);
    
/// Converts a Visit to a Mapping. Throws an exception if the Visit is of a Snarl instead
/// of a Node. Uses a graph to get node length.
inline Mapping to_mapping(const Visit& visit, const HandleGraph& vg);

/// Convert a snarl traversal into an alignment
inline Alignment to_alignment(const SnarlTraversal& trav, const HandleGraph& graph);

/// Copies the boundary Visits from one Snarl into another
inline void transfer_boundary_info(const Snarl& from, Snarl& to);

/// Make an edge_t from a pair of visits
edge_t to_edge(const handlegraph::HandleGraph& graph, const Visit& v1, const Visit& v2);

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
    for (auto iter = begin; iter != end; ++iter) {
        add_snarl(*iter);
    }
    // record the tree structure and build the other indexes
    finish();
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
    
inline Visit to_visit(const Mapping& mapping, bool make_full_node_match) {
    if (!make_full_node_match) {
        // If we're not explicitly coercing the mapping to a full node match, make sure it already is one.
        assert(mapping_is_match(mapping));
        assert(mapping.position().offset() == 0);
    }
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

inline Visit to_visit(const handlegraph::HandleGraph& graph, const handle_t& handle) {
    return to_visit(graph.get_id(handle), graph.get_is_reverse(handle));
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
    
inline Mapping to_mapping(const Visit& visit, const HandleGraph& graph) {
    return to_mapping(visit, [&](id_t id) {
            return graph.get_length(graph.get_handle(id));
        });
}

inline Alignment to_alignment(const SnarlTraversal& trav, const HandleGraph& graph) {
    Alignment aln;
    Path* path = aln.mutable_path();
    for (int i = 0; i < trav.visit_size(); ++i) {
        *path->add_mapping() = to_mapping(trav.visit(i), graph);
    }
    return aln;
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

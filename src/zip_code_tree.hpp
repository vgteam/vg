#ifndef VG_ZIP_CODE_TREE_HPP_INCLUDED

#define VG_ZIP_CODE_TREE_HPP_INCLUDED

#include "zip_code.hpp"
#include "snarl_seed_clusterer.hpp"

#include <stack>

namespace vg{
using namespace std;

/**

A ZipCodeTree represents of set of SnarlDistanceIndexCluserer::Seed's (seed alignments between a read and reference) 
as a tree structure.
The tree represents the connectivity of the seeds, based on the distance index.
Edges are labelled with distance values.
The tree can be traversed to find distances between seeds

This provides an iterator that, given a seed and a distance limit, iterates through seeds that are
reachable within the distance limit

The ZipCodeTree is constructed by the ZipCodeForest, which represents a collection of trees

*/
class ZipCodeTree {

    typedef SnarlDistanceIndexClusterer::Seed Seed;

    public:

    /// Constructor
    ZipCodeTree(const vector<Seed>* all_seeds) : seeds(all_seeds){};

    /*
      The tree will represent the seeds' placement in the snarl tree.
      Each node in the tree represents either a seed (position on the graph) or the 
      boundary of a snarl or chain.
      Edges are labelled with the distance between the two nodes

      This graph is actually represented as a vector of the nodes and edges
      Each item in the vector represents either a node (seed or boundary), an edge (distance),
      or the child count of a snarl

      A chain in the vector is bounded by a CHAIN_START and a CHAIN_END.
      The chain is comprised of alternating children (seed or snarl) and the distances between them,
      starting and ending with a child. The order would be:
        CHAIN_START, child, distance, child, distance, ..., child, CHAIN_END
      The distances represent the number of nucleotides on the minimum-length path in the variation graph
      between the structures that the zip code tree nodes represent.
      For distances terminating at a SEED, the distance includes the nucleotide the position is on.
      For distances between two SEEDs, the distance includes both of the positions.
      For two SEEDs on the same position, the distance between them would be 1.
      For distances terminating at a SNARL_START or SNARL_END, the distance reaches the inner edge
      (relative to the snarl) of the boundary node, so it includes the length of the boundary node of the snarl
      For example, given a subgraph of a chain:

                              n3
                            [GACG] ...
            n1      n2    / 
            [A] - [AGAC]
                          \   n4
                            [ACAG] ...

      for the sequence "SEED EDGE SNARL_START" representing a seed on n1 and the snarl starting at n2,
      the edge value would be 5.


      A snarl in the vector is bounded by a SNARL_START and a SNARL_END.
      A snarl is comprised of the two bounds, one or more chains, and the distances among them.
      SEEDs are always contained within a chain.
      For each element of the snarl (boundary or child chain), the distance to each element preceding
      it in the snarl is stored before the element.
      The distances are stored in reverse order of the elements that they reach.
      Immediately before the SNARL_END, there is a NODE_COUNT storing the number of children in the snarl
      A snarl would look like:
          SNARL_START, dist:start->c1, chain1, dist:c1->c2, dist:start->c2, chain2, ..., 
                ..., dist:c2->end, dist:c1->end, dist:start->end, node_count, SNARL_END


      Everything is ordered according to the order of the highest-level chain (top-level chain or child
      of a top-level snarl).
      For children of a snarl, the children are ordered according to the distance to the start of the snarl,
      and if that value is equal, in reverse order to the distance to the end of the snarl.
      In the variation graph, all chains are considered to be oriented "forward" in their parent snarl.
      However, in a start-to-end traversal of the snarl, the child chain may be traversed end-to-start.
      These chains would be considered to be reversed in the zip code tree, so the order of the children
      of the chain may be backwards relative to their order in the variation graph.
      If a snarl is the child of a chain that is traversed backwards in the zip tree, then that snarl
      and all its children are also traversed backwards.


      TODO: This is still just for DAGS
     */

    public:
    enum tree_item_type_t {SEED, SNARL_START, SNARL_END, CHAIN_START, CHAIN_END, EDGE, NODE_COUNT};
    struct tree_item_t {

        //Is this a seed, boundary, or an edge
        tree_item_type_t type;

        //For a seed, the index into seeds
        //For an edge, the distance value
        //Empty for a bound
        size_t value;

        //For seeds, is the position of the seed traversed backwards in the tree?
        bool is_reversed;
    };

private:
    /*************
     The actual data being stored
     ************/

    //The seeds that are taken as input
    const vector<Seed>* seeds;

protected:
    //The actual tree structure
    vector<tree_item_t> zip_code_tree;


public:

    ///Print the zip code tree to stderr
    /// ( and ) are used for the starts and ends of snarls
    /// [ and ] are used for the starts and ends of chains
    /// seeds are printed as their positions
    void print_self() const;

    /// Is the given node in a multicomponent chain, looping chain, or anything else that would cause
    /// it to not have exact distances?
    /// The distances are only guaranteed to be correct up to the given distance limit
    /// Cyclic snarls don't count as being invalid
    bool node_is_invalid(nid_t id, const SnarlDistanceIndex& distance_index, size_t distance_limit = std::numeric_limits<size_t>::max()) const;

    /// Is the node in a cyclic (non-dag) snarl?
    bool node_is_in_cyclic_snarl(nid_t id, const SnarlDistanceIndex& distance_index) const; 

    ///Check that the tree is correct
    void validate_zip_tree(const SnarlDistanceIndex& distance_index, size_t distance_limit = std::numeric_limits<size_t>::max()) const;

    ///Get the number of items in the tree
    size_t get_tree_size() const {return zip_code_tree.size();};

    ///Helper function to access the values in the zip_code_tree
    tree_item_t get_item_at_index(size_t index) const {return zip_code_tree[index];};

    /// Count the number of snarls involved in the tree
    /// Returns a pair of <dag count, non-dag count>
    /// Assumes that the tree has already been filled in
    std::pair<size_t, size_t> dag_and_non_dag_snarl_count(const vector<Seed>& all_seeds, const SnarlDistanceIndex& distance_index) const;

protected:

    //Helper function to get the orientation of a snarl tree node at a given depth
    //does the same thing as the zipcode decoder's get_is_reversed_in_parent, except
    //that is also considers chains that are children of irregular snarls.
    //We assume that all snarls are DAGs, so all children of snarls must only be
    //traversable in one orientation through the snarl. In a start-to-end traversal
    //of a snarl, each node will only be traversable start-to-end or end-to-start.
    //If it is traversable end-to-start, then it is considered to be oriented
    //backwards in its parent
    //TODO: Move this into the cpp file but I can't figure out how to make it const static
    const static bool seed_is_reversed_at_depth (const Seed& seed, size_t depth, const SnarlDistanceIndex& distance_index){
        if (seed.zipcode_decoder->get_is_reversed_in_parent(depth)) {
            return true;
        } else if (depth > 0 && seed.zipcode_decoder->get_code_type(depth-1) == ZipCode::IRREGULAR_SNARL) {
            //If the parent is an irregular snarl, then check the orientation of the child in the snarl
            net_handle_t snarl_handle = seed.zipcode_decoder->get_net_handle(depth-1, &distance_index);
            size_t rank = seed.zipcode_decoder->get_rank_in_snarl(depth);
            if (distance_index.distance_in_snarl(snarl_handle, 0, false, rank, false)
                        == std::numeric_limits<size_t>::max()
                &&
                distance_index.distance_in_snarl(snarl_handle, 1, false, rank, true)
                        == std::numeric_limits<size_t>::max()) {
                //If the distance from the start of the snarl to the start of the child is infinite
                //and the distance from the end of the snarl to the end of the child is infinite
                //then we assume that this child is "reversed" in the parent snarl
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }



public:
    
    /**
     * Exposed type for a reference to an orientation of a seed.
     */
    struct oriented_seed_t {
        size_t seed;
        bool is_reverse;
        
        /// Compare to other instances. TODO: Use default when we get C++20. 
        inline bool operator==(const oriented_seed_t& other) const {
            return seed == other.seed && is_reverse == other.is_reverse;
        }

        /// Compare to other instances. TODO: Use default when we get C++20. 
        inline bool operator!=(const oriented_seed_t& other) const {
            return !(*this == other);
        }
    };

    /**
     * Exposed type for a reference to an oriented seed at an associated distance.
     */
    struct seed_result_t : public oriented_seed_t {
        size_t distance;

        /// Compare to other instances. TODO: Use default when we get C++20. 
        inline bool operator==(const seed_result_t& other) const {
            return distance == other.distance && oriented_seed_t::operator==((oriented_seed_t)other);
        }

        /// Compare to other instances. TODO: Use default when we get C++20. 
        inline bool operator!=(const seed_result_t& other) const {
            return !(*this == other);
        }
    };

    /**
     * Iterator that visits all seeds right to left in the tree's in-order traversal.
     */
    class iterator {
    public:
        /// Make an iterator wrapping the given iterator, until the given end.
        iterator(vector<tree_item_t>::const_iterator begin, vector<tree_item_t>::const_iterator end);
        
        // Iterators are copyable and movable.
        iterator(const iterator& other) = default;
        iterator(iterator&& other) = default;
        iterator& operator=(const iterator& other) = default;
        iterator& operator=(iterator&& other) = default;

        /// Advance right
        iterator& operator++();

        /// Compare for equality to see if we hit end
        bool operator==(const iterator& other) const;

        /// Compare for inequality
        inline bool operator!=(const iterator& other) const {
            return !(*this == other);
        }
        
        /// Get the index and orientation of the seed we are currently at.
        oriented_seed_t operator*() const;

        /// Get the number of tree storage slots left in the iterator. We need
        /// this to make reverse iterators from forward ones.
        size_t remaining_tree() const;

    private:
        /// Where we are in the stored tree.
        vector<tree_item_t>::const_iterator it;
        /// Where the stored tree ends. We keep this to avoid needing a reference back to the ZipCodeTree.
        vector<tree_item_t>::const_iterator end;
    };

    /// Get an iterator over indexes of seeds in the tree, left to right.
    iterator begin() const;
    /// Get the end iterator for seeds in the tree, left to right.
    iterator end() const;

    /**
     * Iterator that looks left in the tree from a seed, possibly up to a maximum base distance.
     *
     * See https://github.com/benedictpaten/long_read_giraffe_chainer_prototype/blob/b590c34055474b0c901a681a1aa99f1651abb6a4/zip_tree_iterator.py.
     */
    class reverse_iterator {
    public:
        /// Make a reverse iterator wrapping the given reverse iterator, until
        /// the given rend, with the given distance limit.
        reverse_iterator(vector<tree_item_t>::const_reverse_iterator rbegin, vector<tree_item_t>::const_reverse_iterator rend, size_t distance_limit = std::numeric_limits<size_t>::max());

        // Reverse iterators need to be copyable for STL algorithms despite the relatively large stack.
        reverse_iterator(const reverse_iterator& other) = default;
        reverse_iterator(reverse_iterator&& other) = default;
        reverse_iterator& operator=(const reverse_iterator& other) = default;
        reverse_iterator& operator=(reverse_iterator&& other) = default;

        /// Move left
        reverse_iterator& operator++();

        /// Compare for equality to see if we hit end (the past-the-left position)
        bool operator==(const reverse_iterator& other) const;

        /// Compare for inequality
        inline bool operator!=(const reverse_iterator& other) const {
            return !(*this == other);
        }
        
        /// Get the index and orientation of the seed we are currently at, and the distance to it.
        seed_result_t operator*() const;

        /// Type for the state of the
        /// I-can't-believe-it's-not-a-pushdown-automaton
        enum State {
            S_START,
            S_SCAN_CHAIN,
            S_STACK_SNARL,
            S_SCAN_SNARL,
            S_SKIP_CHAIN
        };

    private:
        /// Where we are in the stored tree.
        vector<tree_item_t>::const_reverse_iterator it;
        /// Where the rend is where we have to stop
        vector<tree_item_t>::const_reverse_iterator rend;
        /// Distance limit we will go up to
        size_t distance_limit;
        /// Stack for computing distances
        std::stack<size_t> stack;

        // Now we define a mini stack language so we can do a
        // not-really-a-pushdown-automaton to parse the distance strings.
    
        /// Push a value to the stack
        void push(size_t value);

        /// Pop a value from the stack and return it
        size_t pop();

        /// Get a mutable reference to the value on top of the stack
        size_t& top();

        /// Duplicate the top item on the stack
        void dup();

        /// Check stack depth
        size_t depth() const;

        /// Reverse the top two elements of the stack
        void swap();

        /// Current state of the automaton
        State current_state;

        /// Adopt a new state.
        void state(State new_state);

        /// Stop parsing because nothing else can be below the distance limit.
        /// This moves the current iterator it.
        void halt();

        /// Tick the automaton, looking at the symbol at *it and updating the
        /// stack and current_state. Returns true to yield a value at the
        /// current symbol, or to halt, and false otherwise.
        bool tick();

    };

    /// Get a reverse iterator looking left from where a forward iterator is, up to a distance limit.
    reverse_iterator look_back(const iterator& from, size_t distance_limit = std::numeric_limits<size_t>::max()) const;
    /// Get the reverse end iterator for looking back from seeds.
    reverse_iterator rend() const;

    friend class ZipCodeForest;

}; 
/**
    A collection of ZipCodeTrees
    The ZipCodeForest takes a set of seeds and makes ZipCodeTrees
    There will be a separate tree for each connected component or slice of a chain that is
    too far from anything else on both sides, using the given distance limit
*/
class ZipCodeForest {

    typedef SnarlDistanceIndexClusterer::Seed Seed;
    typedef ZipCodeTree::tree_item_type_t tree_item_type_t;
    typedef ZipCodeTree::tree_item_t tree_item_t;

    public:

    ///The actual data, a collection of ZipCodeTrees
    vector<ZipCodeTree> trees;

    ///Constructor
    ZipCodeForest() {};

    ///Populate the zip forest
    /// If a distance limit is given, then also partition the tree into subtrees that are
    /// farther than the distance_limit from each other
    /// Otherwise, the forest will just be connected components
    /// If a distance limit is given, then distances larger than the distance limit are not
    /// guaranteed to be accurate
    void fill_in_forest(const vector<Seed>& all_seeds, const SnarlDistanceIndex& distance_index,
                      size_t distance_limit = std::numeric_limits<size_t>::max());
    private:
    //The seeds that are taken as input
    //The order of the seeds will never change, but the vector is not const because the zipcodes
    //decoders may change
    const vector<Seed>* seeds;

    public:

    void print_self() const {
        for (size_t i = 0 ; i < trees.size() ; i++) {
            const auto& tree = trees[i];
            cerr << i << ": ";
            tree.print_self();
        }
    }
    void validate_zip_forest(const SnarlDistanceIndex& distance_index, size_t distance_limit=std::numeric_limits<size_t>::max()) const {
        for (const auto& tree : trees) {
            tree.validate_zip_tree(distance_index, distance_limit);
        }
    }


    /************************
      Helper functions for construction
     ***********************/
    private:

    /// This gets used for sorting
    /// It represents one interval along zipcode_sort_order to be sorted
    /// At the relevant depth, everything in the interval will be on the same
    /// snarl tree node, and is_reversed is true if that snarl tree node
    /// is reversed relative to the top-level chain
    struct interval_and_orientation_t {
        size_t interval_start : 26; //inclusive
        size_t interval_end : 26;   //exclusive
        bool is_reversed : 1;
        ZipCode::code_type_t code_type : 5;
        size_t depth : 16;

        interval_and_orientation_t (size_t start, size_t end, size_t rev, ZipCode::code_type_t type, 
                                    size_t depth) :
            interval_start(start), interval_end(end), is_reversed(rev), code_type(type), depth(depth){}
    };

    /// Sorts the given interval (which must contain seeds on the same snarl/chain/node at the given depth)
    /// and return the intervals of the children, in the order of traversal
    /// Sorting is roughly linear along the top-level chains, in a topological-ish order in snarls
    /// Uses radix_sort_zipcodes and default_sort_zipcodes
    /// sort_root is true if sorting the root into connected components
    vector<interval_and_orientation_t> sort_one_interval(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval, 
            size_t interval_depth, const SnarlDistanceIndex& distance_index) const;

    /// Helper function to sort the seeds using radix sort
    /// Sorts the slice of seeds in the given interval of zipcode_sort_order, which is a vector of indices
    /// into seeds
    /// reverse_order is true if the order should be reversed. The interval also has an is_reversed field,
    /// which refers to the orientation in the snarl tree
    /// This should run in linear time, but it is dependent on the values being sorted on to have a small range
    void radix_sort_zipcodes(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval,
                             bool reverse_order, size_t depth, const SnarlDistanceIndex& distance_index, 
                             const std::function<size_t(const Seed& seed, size_t depth)>& get_sort_value) const; 

    /// Helper function to sort the seeds using std::sort
    /// Sorts the slice of seeds in the given interval of zipcode_sort_order, which is a vector of indices
    /// into seeds
    void default_sort_zipcodes(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval,
                             bool reverse_order, size_t depth, const SnarlDistanceIndex& distance_index, 
                             const std::function<size_t(const Seed& seed, size_t depth)>& get_sort_value) const; 

    /// Helper function to sort the seeds on a cyclic (non-dag) snarl
    /// depth is the depth of the snarl
    /// Returns the intervals on zipcode_sort_order
    /// The intervals may be duplicated and in different orientations
    vector<interval_and_orientation_t> sort_zipcodes_on_cyclic_snarl(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval,
                             size_t depth, const SnarlDistanceIndex& distance_index) const; 

    //////////////////// data structures and helper functions for building the forest

    //For children of snarls, we need to remember the siblings and start bound that came before them
    //so we can record their distances
    //This holds the indices (into zip_code_tree) of each seed or start of a chain,
    // and each start and child chain start of a snarl
    //The children are stored at the depth of their parents. For example, for a root chain,
    //the vector at index 0 would have the chain start, seeds that are on the chain, and the start
    //of snarls on the chain. Similarly, for a top-level snarl, at depth 1, the second vector would contain
    //the starts of chains at depth 2 
    //For the children of a chain, the value is the prefix sum in the chain (relative to the orientation 
    //of the top-level chain, not necessarily the chain itself)
    //For the children of a snarl, the value is the index of the CHAIN_START in zip_code_tree.
    //  The first seed in the chain will need to be found by looping through zip_code_tree
    struct child_info_t {
        ZipCodeTree::tree_item_type_t type;  //the type of the item
        size_t value;  //A value associated with the item, either the offset in a chain, index of the snarl child start
    
        //For the children of snarls, the distance to the left and right of the chain, that gets added to
        //edges in the snarl
        std::pair<size_t, size_t> distances;

        //Is the sibling reversed. 
        //This is only used for children of snarls, to indicate that the child is traversed backwards 
        bool is_reversed = false;
    };


    /// This stores information about the state of the forest as we fill it in
    struct forest_growing_state_t {

        vector<size_t> seed_sort_order;

        //Stores the previous things of the current structure at each depth
        vector<vector<child_info_t>> sibling_indices_at_depth;

        // We build a forest of trees. A new tree is formed either when a new top-level chain is found
        // (or a slice of a top-level chain if it is far enough away from the previous thing in the chain),
        // or when part of a chain in a snarl is too far from everything else in the snarl.
        // In the second case, the entire subtree is found before determining that it should be a subtree,
        // and then it is copied into a new zip_tree_t in the forest.
        // So only one tree is actively being added to at a time.
        //This keeps track of which is the active tree, as an index into trees
        size_t active_zip_tree;

        // Keep track of all open chains as an index into the current active_zip_tree of the start of the chain,
        // and a boolean that is true if the start of the chain is farther than the distance_limit from anything
        // else in the snarl tree
        // If the index is pointing to a CHAIN_START, then it includes the whole chain. If it points to a SEED,
        // then it is a slice
        // Any time something gets added to a chain or the chain is closed, check if the distance to anything
        // following is greater than the distance limit. If it is, copy everything from the start of the chain
        // or slice into a new tree in the forest.
        vector<pair<size_t, bool>> open_chains;

        //A stack of intervals representing snarl tree nodes. These are yet to be sorted and processed
        vector<interval_and_orientation_t> intervals_to_process;
    
        //Intervals that are currently open. These represent ancestors of whatever is currently being worked on
        //So the size is the depth of the snarl tree
        vector<interval_and_orientation_t> open_intervals;
    


    };

    // Open a chain that starts at the current_seed
    // If the chain is in a snarl, then add empty edges for the distances to everything before it in the snarl
    // Open the chain, and record its presence and distance-to-start in the parent snarl, if necessary
    void open_chain(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                      const size_t& distance_limit, const size_t& depth, const Seed& current_seed, 
                      bool chain_is_reversed);
    // Close a chain that ends at last_seed
    // If the chain was empty, remove it and anything relating to it in the parent snarl and sibling_indices
    // If it can be spliced out, take out a subtree
    // Otherwise, add the end of the chain and, if the chain was in a snarl, add the distances to everything
    // before it in the snarl and remember the distance to the end of the chain
    void close_chain(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                      const size_t& distance_limit, const size_t& depth, const Seed& last_seed,
                      bool chain_is_reversed);

    // Add the current seed (or snarl starting at the seed) and its distance to the previous thing in a chain
    // If the seed is far enough from the previous thing in the chain and it can be a new slice, split off 
    // a subtree
    // depth is the depth of the child of the chain (which may also be the chain depth if it is trivial)
    // seed_index is the index of the current seed in the list of seeds
    void add_child_to_chain(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                      const size_t& distance_limit, const size_t& depth, const size_t& seed_index, 
                      bool child_is_reversed, bool chain_is_reversed);

    // Start a new snarl
    void open_snarl(forest_growing_state_t& forest_state, const size_t& depth);

    // Close a snarl
    // depth is the depth of the snarl and last_seed is the last seed in the snarl
    // If the snarl has no children, then delete the whole thing
    // Otherwise, add all necessary distances and close it
    void close_snarl(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                      const size_t& depth, const Seed& last_seed, bool last_is_reversed, bool is_cyclic_snarl);

    // Add all the distances from everything in the snarl to either the last child of the snarl or,
    // if to_snarl_end is true, to the end bound of the snarl
    // depth is the depth of the snarl
    void add_snarl_distances(forest_growing_state_t& forest_state, const SnarlDistanceIndex& distance_index,
                             const size_t& depth, const Seed& seed, bool child_is_reversed, bool snarl_is_reversed, 
                             bool to_snarl_end, bool is_cyclic_snarl);

};

/// Print an item type to a stream
std::ostream& operator<<(std::ostream& out, const ZipCodeTree::tree_item_type_t& type);
/// Pritn an iterator state to a stream
std::ostream& operator<<(std::ostream& out, const ZipCodeTree::reverse_iterator::State& state);

}

namespace std {

/// Make an item type into a string
std::string to_string(const vg::ZipCodeTree::tree_item_type_t& type);
/// Make an iterator state into a string
std::string to_string(const vg::ZipCodeTree::reverse_iterator::State& state);

/// Hash functor to hash oriented_seed_t with std::hash
template <> struct hash<vg::ZipCodeTree::oriented_seed_t>
{
    /// Produce a hash of an oriented_seed_t.
    size_t operator()(const vg::ZipCodeTree::oriented_seed_t& item) const
    {
        // Hash it just as we would a pair.
        return hash<pair<size_t, bool>>()(make_pair(item.seed, item.is_reverse));
    }
};

/// Hash functor to hash seed_result_t with std::hash
template <> struct hash<vg::ZipCodeTree::seed_result_t>
{
    /// Produce a hash of a seed_result_t.
    size_t operator()(const vg::ZipCodeTree::seed_result_t& item) const
    {
        // Hash it just as we would a tuple.
        return hash<tuple<size_t, bool, size_t>>()(make_tuple(item.seed, item.is_reverse, item.distance));
    }
};

/// Explain to the STL algorithms what kind of iterator the zip code tree
/// forward iterator is.
template<>
struct iterator_traits<vg::ZipCodeTree::iterator>{
    using value_type = vg::ZipCodeTree::oriented_seed_t;   
    using iterator_category = forward_iterator_tag;
};

/// Explain to the STL algorithms what kind of iterator the zip code tree
/// reverse iterator is.
template<>
struct iterator_traits<vg::ZipCodeTree::reverse_iterator>{
    using value_type = vg::ZipCodeTree::seed_result_t;   
    using iterator_category = forward_iterator_tag;
};


}

#endif

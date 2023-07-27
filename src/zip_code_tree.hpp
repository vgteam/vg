#ifndef VG_ZIP_CODE_TREE_HPP_INCLUDED

#define VG_ZIP_CODE_TREE_HPP_INCLUDED

#include "zip_code.hpp"
#include "snarl_seed_clusterer.hpp"

#include <stack>

namespace vg{
using namespace std;

/**

A ZipCodeTree takes a set of SnarlDistanceIndexCluserer::Seed's (seed alignments between a read and reference) 
and provides an iterator that, given a seed and a distance limit, iterates through seeds that are
reachable within the distance limit

Generally, this will take a collection of seeds and build a tree structure representing the connectivity
of the seeds, based on the snarl decomposition
Edges are labelled with distance values.
The tree can be traversed to find distances between seeds
*/
class ZipCodeTree {

    typedef SnarlDistanceIndexClusterer::Seed Seed;

    public:

    /**
     * Constructor
     * The constructor creates a tree of the input seeds that is used for calculating distances
     */
    ZipCodeTree(){};
    void fill_in_tree(vector<Seed>& all_seeds, const SnarlDistanceIndex& distance_index);


    private:

    //The seeds that are taken as input
    //The order of the seeds will never change, but the vector is not const because the zipcodes
    //decoders may change
    vector<Seed>* seeds;


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
    //The actual tree structure
    vector<tree_item_t> zip_code_tree;

public:

    /// Return the sort order of the seeds
    /// Sorting is roughly linear along the top-level chains, in a topological-ish order in snarls
    /// Uses radix_sort_zipcodes and default_sort_zipcodes
    vector<size_t> sort_seeds_by_zipcode(const SnarlDistanceIndex& distance_index) const;

    /// Count the number of snarls involved in the tree
    /// Returns a pair of <dag count, non-dag count>
    /// Assumes that the tree has already been filled in
    std::pair<size_t, size_t> dag_and_non_dag_snarl_count(vector<Seed>& all_seeds, const SnarlDistanceIndex& distance_index) const;

    ///Print the zip code tree to stderr
    /// ( and ) are used for the starts and ends of snarls
    /// [ and ] are used for the starts and ends of chains
    /// seeds are printed as their positions
    void print_self() const;

    ///Helper function that returns the number of items in the zip_code_tree
    size_t get_tree_size() const {return zip_code_tree.size();};

    ///Check that the tree is correct
    void validate_zip_tree(const SnarlDistanceIndex& distance_index) const;

    ///Helper function to access the values in the zip_code_tree
    tree_item_t get_item_at_index(size_t index) const {return zip_code_tree[index];};
private:

    //Helper function to get the orientation of a snarl tree node at a given depth
    //does the same thing as the zipcode decoder's get_is_reversed_in_parent, except
    //that is also considers chains that are children of irregular snarls.
    //We assume that all snarls are DAGs, so all children of snarls must only be
    //traversable in one orientation through the snarl. In a start-to-end traversal
    //of a snarl, each node will only be traversable start-to-end or end-to-start.
    //If it is traversable end-to-start, then it is considered to be oriented
    //backwards in its parent
    bool seed_is_reversed_at_depth (const Seed& seed, size_t depth, const SnarlDistanceIndex& distance_index) const;

    /// This gets used for sorting
    /// It represents one interval along zipcode_sort_order to be sorted
    /// At the relevant depth, everything in the interval will be on the same
    /// snarl tree node, and is_reversed is true if that snarl tree node
    /// is reversed relative to the top-level chain
    struct interval_and_orientation_t {
        size_t interval_start : 32; //inclusive
        size_t interval_end : 31;   //exclusive
        bool is_reversed : 1;

        interval_and_orientation_t (size_t start, size_t end, size_t rev) :
            interval_start(start), interval_end(end), is_reversed(rev) {}
    };

    /// Helper function to sort the seeds using radix sort
    /// Sorts the slice of seeds in the given interval of zipcode_sort_order, which is a vector of indices
    /// into seeds
    /// This should run in linear time, but it is dependent on the values being sorted on to have a small range
    void radix_sort_zipcodes(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval,
                             vector<interval_and_orientation_t>& new_intervals, size_t depth,
                             const SnarlDistanceIndex& distance_index, 
                             const std::function<size_t(Seed& seed, size_t depth)>& get_sort_value) const; 

    /// Helper function to sort the seeds using std::sort
    /// Sorts the slice of seeds in the given interval of zipcode_sort_order, which is a vector of indices
    /// into seeds
    void default_sort_zipcodes(vector<size_t>& zipcode_sort_order, const interval_and_orientation_t& interval,
                             vector<interval_and_orientation_t>& new_intervals, size_t depth,
                             const SnarlDistanceIndex& distance_index, 
                             const std::function<size_t(Seed& seed, size_t depth)>& get_sort_value) const; 



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

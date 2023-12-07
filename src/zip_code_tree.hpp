#ifndef VG_ZIP_CODE_TREE_HPP_INCLUDED

#define VG_ZIP_CODE_TREE_HPP_INCLUDED

//#define DEBUG_ZIP_CODE_TREE
//#define DEBUG_ZIP_CODE_SORTING

#include "zip_code.hpp"
#include "snarl_seed_clusterer.hpp"

#include <stack>
#include <forward_list>

namespace vg{
using namespace std;

/**

A ZipCodeTree represents of set of SnarlDistanceIndexCluserer::Seed's (seed alignments between a 
read and reference) as a tree structure.
The tree represents the connectivity of the seeds, based on the distance index.
Edges are labelled with distance values.
The tree can be traversed to find distances between seeds

This provides an iterator that, given a seed and a distance limit, iterates through seeds that are
reachable within the distance limit

The ZipCodeTree is built by the ZipCodeForest, which represents a collection of trees

*/
class ZipCodeTree {

    typedef SnarlDistanceIndexClusterer::Seed Seed;

    public:

    /// Empty constructor
    /// ZipCodeTree's get filled in by ZipCodeForest's
    ZipCodeTree(){};

    /*
      The tree will represent the seeds' placement in the snarl tree.
      Each node in the tree represents either a seed (position on the graph, representing the start
      of an alignment) or the boundary of a snarl or chain.
      Edges are labelled with the distance between the two nodes

      This graph is actually represented as a vector of the nodes and edges
      Each item in the vector represents either a node (seed or boundary), an edge (distance),
      or the child count of a snarl

      A chain in the vector is bounded by a CHAIN_START and a CHAIN_END.
      The chain is comprised of alternating children (seed or snarl) and the distances between them,
      starting and ending with a child. The order would be:
        CHAIN_START, child, distance, child, distance, ..., child, CHAIN_END
      The distance from the chain start to the first child is included in the distances in the 
      chain's parent snarl, if relevant

      The distances represent the number of nucleotides on the minimum-length path in the variation 
      graph between the structures that the zip code tree nodes represent.
      Seeds represent the first nucleotide of the alignment, so when the seed is traversed forwards
      in the zip tree, the distance includes the position. If the seed is reversed in the zip tree,
      then the distance doesn't include the position
      For two SEEDs on the same position, the distance between them would be 0.
      For chain distances terminating at a SNARL_START or SNARL_END, the distance reaches the inner 
      edge (relative to the snarl) of the boundary node, so it includes the length of the boundary 
      node of the snarl

      For example, given a subgraph of a chain:

                              n3
                            [GACG] ...
            n1      n2    / 
            [A] - [AGAC]
                          \   n4
                            [ACAG] ...

      for the sequence "SEED EDGE SNARL_START" representing a seed on n1 and the snarl starting at 
      n2, the edge value would be 5.
      Within the snarl, the edge distances include the distance to the first seed in the chain.
      For a seed at position node 3 +1 (the A oriented forwards), the sequence would be
      "SNARL_START EDGE CHAIN_START SEED", and the edge value would be 1


      A snarl in the vector is bounded by a SNARL_START and a SNARL_END.
      A snarl is comprised of the two bounds, one or more chains, and the distances among them.
      SEEDs are always contained within a chain.
      For each element of the snarl (boundary or child chain), the distance to each element 
      preceding it in the snarl is stored before the element.
      The distances are stored in reverse order of the elements that they reach.
      Immediately before the SNARL_END, there is a NODE_COUNT storing the number of children in the
      snarl. A snarl would look like:
          SNARL_START, dist:start->c1, chain1, dist:c1->c2, dist:start->c2, chain2, ..., 
                ..., dist:c2->end, dist:c1->end, dist:start->end, node_count, SNARL_END

      For snarls that aren't dags (called cyclic snarls, even though they could have an inversion 
      and no cycles), the zip tree should represent all possible paths that the read could take 
      through the snarl. All seeds on the snarl are split up into "runs" of seeds on the same chain
      that are "close" to each other. The runs are sorted and orientated by their read coordinate 
      and each run is made into a separate child chain like normal. A run may occur twice, once in 
      each orientation. See get_cyclic_snarl_intervals() for details 


      Everything is ordered according to the order of the highest-level chain (top-level chain or 
      child of a top-level snarl).
      For children of a snarl, the children are ordered according to a topological sort of the 
      snarl. In the variation graph, all chains are considered to be oriented "forward" in their 
      parent snarl. However, in a start-to-end traversal of the snarl, the child chain may be 
      traversed end-to-start. These chains would be considered to be reversed in the zip code tree,
      so the order of the children of the chain may be backwards relative to their order in the 
      variation graph. If a snarl is the child of a chain that is traversed backwards in the zip 
      tree, then that snarl and all its children are also traversed backwards.

     */

    public:

    ///The type of an item in the zip code tree
    enum tree_item_type_t {SEED=0, SNARL_START, SNARL_END, CHAIN_START, CHAIN_END, EDGE, NODE_COUNT};

    /// One item in the zip code tree, representing a node or edge of the tree
    struct tree_item_t {

        private:
        //Is this a seed, boundary, or an edge
        tree_item_type_t type : 4;

        //For a seed, the index into seeds
        //For an edge, the distance value
        //Empty for a bound
        size_t value : 59;

        //For seeds, is the position of the seed traversed backwards in the tree?
        bool is_reversed;

        public:

        //Empty constructor
        tree_item_t (){};

        //Constructor so that value gets set properly
        tree_item_t ( tree_item_type_t type, size_t raw_value, bool is_reversed) 
            : type(type), is_reversed(is_reversed) {
            if (raw_value == std::numeric_limits<size_t>::max()) {
                value = ((size_t)1 << 59) - 1;
            } else {
                value = raw_value;
            }
        }
        tree_item_type_t get_type() const { return type; }
        size_t get_value() const { 
            return value == ((size_t)1 << 59) - 1
                   ? std::numeric_limits<size_t>::max()
                   : value;
        }
        bool get_is_reversed() const { return is_reversed; }
    };

    ///Get the number of items in the tree
    size_t get_tree_size() const {return zip_code_tree.size();}

    ///Access the values in the zip_code_tree
    tree_item_t get_item_at_index(size_t index) const {return zip_code_tree[index];};;

protected:
    //The actual tree structure
    vector<tree_item_t> zip_code_tree;

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
        reverse_iterator(vector<tree_item_t>::const_reverse_iterator rbegin, 
                         vector<tree_item_t>::const_reverse_iterator rend, 
                         size_t distance_limit = std::numeric_limits<size_t>::max());

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

    /// Get a reverse iterator looking left from where a forward iterator is, up to a distance limit
    reverse_iterator look_back(const iterator& from, 
                               size_t distance_limit = std::numeric_limits<size_t>::max()) const;
    /// Get the reverse end iterator for looking back from seeds.
    reverse_iterator rend() const;


public:

    /*************** Debugging functions for validating the zip tree ***********/

    ///Print the zip code tree to stderr
    /// ( and ) are used for the starts and ends of snarls
    /// [ and ] are used for the starts and ends of chains
    /// seeds are printed as their positions
    void print_self(const vector<Seed>* seeds) const;

    /// Is the given node in a multicomponent chain, looping chain, or anything else that would cause
    /// it to not have exact distances?
    /// The distances are only guaranteed to be correct up to the given distance limit
    /// Cyclic snarls don't count as being invalid
    bool node_is_invalid(nid_t id, const SnarlDistanceIndex& distance_index, 
                         size_t distance_limit = std::numeric_limits<size_t>::max()) const;

    /// Is the node in a cyclic (non-dag) snarl?
    bool node_is_in_cyclic_snarl(nid_t id, const SnarlDistanceIndex& distance_index) const; 

    ///Check that the tree is correct
    void validate_zip_tree(const SnarlDistanceIndex& distance_index, 
                           const vector<Seed>* seeds,
                           size_t distance_limit = std::numeric_limits<size_t>::max()) const;

    ///Helper function for validate_zip_tree for just a snarl
    void validate_snarl(std::vector<tree_item_t>::const_iterator zip_iterator, 
                        const SnarlDistanceIndex& distance_index, 
                        const vector<Seed>* seeds,
                        size_t distance_limit = std::numeric_limits<size_t>::max()) const;


    /// Count the number of snarls involved in the tree
    /// Returns a pair of <dag count, non-dag count>
    /// Assumes that the tree has already been filled in
    std::pair<size_t, size_t> dag_and_non_dag_snarl_count(const vector<Seed>& seeds, 
                                                    const SnarlDistanceIndex& distance_index) const;

protected:

    //Helper function to get the orientation of a snarl tree node at a given depth
    //does the same thing as the zipcode decoder's get_is_reversed_in_parent, except
    //that is also considers chains that are children of irregular snarls.
    //We assume that all snarls are DAGs, so all children of snarls must only be
    //traversable in one orientation through the snarl. In a start-to-end traversal
    //of a snarl, each node will only be traversable start-to-end or end-to-start.
    //If it is traversable end-to-start, then it is considered to be oriented
    //backwards in its parent
    static bool seed_is_reversed_at_depth (const Seed& seed, size_t depth, 
                                           const SnarlDistanceIndex& distance_index);




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
    /// The gap_distance_limit is the limit for making runs of seeds in a cyclic snarl- it 
    /// should be roughly the distance that the dynamic programming is willing to jump to 
    /// connect two consecutive minimizers 
    ///TODO: I think the distance_limit should just be the same as the gap_distance_limit
    /// If a distance_limit is given, then distances larger than the distance limit are not
    /// guaranteed to be accurate, but will be greater than the distance_limit
    template<typename Minimizer>
    void fill_in_forest(const vector<Seed>& seeds, const VectorView<Minimizer>& minimizers, 
                      const SnarlDistanceIndex& distance_index,
                      size_t gap_distance_limit,
                      size_t distance_limit = std::numeric_limits<size_t>::max());

    private:


    /***********************************************************************************************

      Data structures and helper functions for construction

      **********************************************************************************************

      Construction is done in a depth-first traversal of the snarl tree. So when each 
      snarl tree node is visited, the start of the structure is added to the zip tree, then each of
      its children is added to the zip tree along with the distances between them, then the end of
      the structure is added.
      
      The traversal of the snarl tree is accomplished by progressively sorting the seeds to identify
      the snarl tree structures that they lie on. Using the zip codes, the seeds can be sorted on 
      each snarl tree structure separately. Seeds along a chain are sorted to be ordered along a 
      chain, and seeds in a snarl are sorted by the child of the snarl that they are on. The seeds 
      get sorted using a radix-like sort on each structure at each depth of the snarl tree, starting
      with the root and moving down. 
      "Intervals" of seeds in the sort order are used to keep track of the location on the snarl 
      tree. An interval represents a range of seeds that are all on the same snarl tree structure.
      After sorting an interval at one depth, sub-intervals representing the children can be found.
      So first, the seeds are sorted into connected components and sliced into intervals 
      representing root-level snarls and chains. Each interval is then sorted to order the seeds
      along the snarl or chain, and new intervals are found representing ranges of seeds on the 
      children.

      Sorting and tree-building are done at the same time, progressively at each structure in the
      snarl tree. The order of tree-building is based on a stack of intervals. The algorithm starts 
      with an interval for each child of the root snarl. An interval is popped from the stack. Any 
      incomplete snarls or chains that the interval is not a child of must be completed. Then, the
      snarl or chain that the interval represents is started in the zip tree, and any relevant 
      distances are added. Intervals representing the children of the snarl or chain are found and
      added to the stack. This repeats until the stack is empty.

      Each snarl and chain in the zip code tree is comprised of the start and end bounds, the 
      children, and distances between children/bounds. So as each child is added, we will need
      to know what came before it in the parent snarl/chain so that we can add the distances. We 
      also need to remember the ancestors of each snarl and chain as we are building them, so that
      we can close each structure properly. All of this information is stored in a 
      forest_growing_state_t as the zip trees are being built.

     **********************************************************************************************/



    private:

    ////////////////////////////////////////////////////
    ///////////
    /////////// Data structures for building a zip tree
    //////////
    ////////////////////////////////////////////////////

    //One interval of seeds corresponding a snarl tree structure
    struct interval_state_t;

    struct sort_value_t;

    //Stores distance information about the child of a structure, so that distances can be
    //found between siblings
    struct child_info_t;

    /// This stores information about the state of the forest as we fill it in
    struct forest_growing_state_t {

        const vector<Seed>* seeds;

        const SnarlDistanceIndex* distance_index;

        vector<size_t> seed_sort_order;


        //This stores the sort value and code type of each seed 
        //This will change as forest building progresses but it will be set for the relevant seed 
        //immediately before sorting
        //The values also get used to calculate distance, as long as they have been set for the 
        //correct depth
        vector<sort_value_t> sort_values_by_seed;

        //Stores the previous things of the current structure at each depth
        //The children are stored at the depth of their parents. For example, for a root chain,
        //the vector at index 0 would have the chain start, seeds that are on the chain, and the 
        //start of snarls on the chain. Similarly, for a top-level snarl, at depth 1, the second 
        //vector would contain the starts of chains at depth 2 
        vector<vector<child_info_t>> sibling_indices_at_depth;

        // We build a forest of trees. A new tree is formed either when a new top-level chain is 
        // found (or a slice of a top-level chain if it is far enough away from the previous thing 
        // in the chain), or when part of a chain in a snarl is too far from everything else in the
        // snarl. In the second case, the entire subtree is found before determining that it should 
        // be a subtree, and then it is copied into a new zip_tree_t in the forest.
        // So only one tree is actively being added to at a time.
        // This keeps track of which is the active tree, as an index into trees
        // Note that this can't be an actual pointer to the forest because the address may move if 
        // the vectors get shifted around in memory.
        size_t active_zip_tree_i;

        // Keep track of all open chains as an index into the current active_zip_tree_i of the start 
        // of the chain, and a boolean that is true if the start of the chain is farther than the 
        // distance_limit from anything else in the snarl tree.
        // If the index is pointing to a CHAIN_START, then it includes the whole chain. If it 
        // points to a SEED, then it is a slice.
        // Any time something gets added to a chain or the chain is closed, check if the distance 
        // to anything following is greater than the distance limit. If it is, copy everything from
        // the start of the chain or slice into a new tree in the forest.
        vector<pair<size_t, bool>> open_chains;

        // A stack of intervals representing snarl tree nodes. These are yet to be sorted and added
        // to the zip tree. After an interval is popped, intervals of its children get added to 
        // intervals_to_process
        // The stack structure ensures that the snarl tree gets processed in the right order
        vector<interval_state_t> intervals_to_process;
    
        //Intervals that are currently open. These represent ancestors of whatever is currently 
        //being worked on. So the size is the depth of the snarl tree
        vector<interval_state_t> open_intervals;
    
        //For cyclic snarls, what is the limit for separating runs of seeds
        size_t gap_distance_limit;

        //The overall distance limit for splitting of new connected components
        size_t distance_limit;

        // Constructor given seeds and a distance index
        forest_growing_state_t(const vector<Seed>& seeds, const SnarlDistanceIndex& distance_index, 
                               size_t gap_distance_limit, size_t distance_limit) :
            seeds(&seeds), distance_index(&distance_index), gap_distance_limit(gap_distance_limit),
            distance_limit(distance_limit), active_zip_tree_i(std::numeric_limits<size_t>::max()) {

            //This represents the current sort order of the seeds
            seed_sort_order.assign(seeds.size(), 0);
            for (size_t i = 0 ; i < seed_sort_order.size() ; i++) {
                seed_sort_order[i] = i;
            }
            sort_values_by_seed.resize(seeds.size());
        }

    };


    /// For children of snarls, we need to remember the siblings and start bound that came before 
    /// them so we can record their distances
    /// This holds the indices (into zip_code_tree) of each seed or start of a chain,
    ///  and each start and child chain start of a snarl
    /// For the children of a chain, the value is the prefix sum in the chain (relative to the 
    /// orientation of the top-level chain, not necessarily the chain itself)
    /// For the children of a snarl, the value is the index of the CHAIN_START in zip_code_tree.
    ///  The first seed in the chain will need to be found by looping through zip_code_tree
    struct child_info_t {

        ZipCodeTree::tree_item_type_t type;  //the type of the item

        //A value associated with the item, either the offset in a chain, index of the snarl 
        //child start
        size_t value;  
    
        //For the children of snarls, the distance to the left and right of the chain, that gets 
        //added to edges in the snarl
        std::pair<size_t, size_t> distances;

        //Is the sibling reversed. 
        //This is only used for children of snarls, to indicate that the child is traversed 
        //backwards 
        bool is_reversed = false;
    };

    /// This gets used for sorting. It represents one interval along zipcode_sort_order, which 
    /// corresponds to a snarl tree node at the given depth
    struct interval_state_t {

        //Indices into zipcode_sort_order
        size_t interval_start : 26; //inclusive
        size_t interval_end : 26;   //exclusive

        // is_reversed is true if that snarl tree node is reversed relative to the 
        // top-level chain
        bool is_reversed : 1;

        //The type of the snarl tree structure.
        // For nodes on chains, all seeds on the chain that aren't nested in snarls are put in 
        // the same interval, regardless of if they are actually on the same node
        ZipCode::code_type_t code_type : 5;

        size_t depth : 14;

        //For children of cyclic snarls, an entire chain may be duplicated in the opposite 
        // orientation immediately after the first copy. In this case, when the second copy is 
        // processed, the entire interval is already in the correct order, just reversed. 
        //If this is_reverse_ordered true, then the interval is sorted in the reverse order, so it
        // needs to be flipped before processing
        bool is_reverse_ordered;

        //After flipping a reverse-ordered interval, all of the child intervals will be sorted
        //So remember if the interval doesn't need sorting
        bool is_ordered;


        interval_state_t (size_t start, size_t end, size_t rev, ZipCode::code_type_t type, 
                                    size_t depth) :
            interval_start(start), interval_end(end), is_reversed(rev), code_type(type), depth(depth){
            is_reverse_ordered = false;
            is_ordered = false;
        }
    };

    ///This is used for storing the value used for sorting seeds
    ///Since children of chains get sorted by the offset along the chain, it can also be used
    ///to find the values used for calculating distances
    struct sort_value_t {
        private:
        size_t sort_value;
        ZipCode::code_type_t code_type;

        // For chains, this is used to indicate the order of the child of a chain
        // Since the offset stored represents the space between nucleotides, two positions on different nodes
        // could have the same offset. Similarly, a snarl could have the same prefix sum as a node.
        // For example, in this graph:
        //                2
        //               [AA]
        //           1  /   \  3
        //          [AA] --- [AA]
        // The positions n1-0 and 3+0, and the snarl 1-3 all have the same offset of 2
        // To solve this, the prefix sum of a chain will always be multiplied by 3, and 1 will be added to snarls,
        // And 2 will be added to the node with an offset in the node of 0 (node 3 if the chain is traversed forward)

        size_t chain_order : 3; 

        public:
        //Constructor
        sort_value_t() : sort_value(std::numeric_limits<size_t>::max()),
                         code_type(ZipCode::EMPTY),
                         chain_order(7) {};
        sort_value_t (size_t sort_value, ZipCode::code_type_t code_type, size_t chain_order) :
            sort_value(sort_value), code_type(code_type), chain_order(chain_order) {};

        //Get the value used for sorting
        size_t get_sort_value() const {
            //The sort value for chains is actually the prefix sum*3+chain_order, 
            // to account for different nodes having the same prefix sum
            return chain_order != 7
                       ? (sort_value * 3) + chain_order
                       : sort_value;
        };

        //Get the value used for distance finding
        size_t get_distance_value() const {return sort_value;};

        //Get the code type
        ZipCode::code_type_t get_code_type() const {return code_type;};

        void set_sort_value(size_t value) {sort_value =value;};
        void set_code_type(ZipCode::code_type_t type) {code_type = type;};
        void set_chain_order(size_t order) {chain_order = order;};

    };

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////   Functions for sorting and finding intervals of seeds along the snarl tree
    /////////////////////////////////////////////////////////////////////////////////////////////


    /// Sorts the given interval (which must contain seeds on the same snarl/chain/node at the given
    /// depth)  Sorting is roughly linear along the top-level chains, in a topological-ish order in
    /// snarls. Uses radix_sort_zipcodes and default_sort_zipcodes
    void sort_one_interval(forest_growing_state_t& forest_state, 
                           const interval_state_t& interval) const;

    /// Helper function to sort the seeds using radix sort
    /// Sorts the slice of seeds in the given interval of zipcode_sort_order, which is a vector of
    /// indices into seeds
    /// reverse_order is true if the order should be reversed. The interval also has an is_reversed
    /// field, which refers to the orientation in the snarl tree
    /// This should run in linear time, but it is dependent on the values being sorted on to have a
    /// small range
    /// min_ and max_value are the minimum and maximum value being sorted on
    void radix_sort_zipcodes(vector<size_t>& zipcode_sort_order, 
                             const vector<sort_value_t>& sort_values_by_seed,
                             const interval_state_t& interval, bool reverse_order,
                             size_t min_value, size_t max_value) const; 

    /// Helper function to sort the seeds using std::sort
    /// Sorts the slice of seeds in the given interval of zipcode_sort_order, which is a vector 
    /// of indices into seeds
    void default_sort_zipcodes(vector<size_t>& zipcode_sort_order, 
                               const vector<sort_value_t>& sort_values_by_seed,
                               const interval_state_t& interval, bool reverse_order) const; 



    /// Assuming that the range of seeds in sort_values_by_seeds given by the interval is sorted,
    /// return the intervals of the children of the interval, in the order of traversal
    /// For children of chains, seeds that are on the chain itself and not nested will be put on
    /// the same interval if there are no seeds in snarls between them, even if they are not on
    /// the same node
    vector<interval_state_t> get_next_intervals(forest_growing_state_t& forest_state, 
                                                const interval_state_t& interval) const;

    /// Given intervals representing child chains on a cyclic snarl, re-partition them and return
    /// new intervals representing runs of seeds that are "close" in each chain
    /// Two seeds are close to each other if: 
    /// (1) the distance between them on the read is <= t, where t is a given distance limit, 
    /// (2) the minimum distance between them on the chain is <= t, and 
    /// (3) they are on the same strand in the read.
    /// Runs are sorted by their latest position in the read, and oriented according to the
    /// orientation of the read through the snarl. The orientation of the read in the snarl's parent
    /// chain and in the snarl children are estimated by finding the spearman correlation of the 
    /// seeds. If the orientation of a run is unclear, then it is duplicated to be oriented in each 
    /// direction 
    template<typename Minimizer>
    vector<interval_state_t> get_cyclic_snarl_intervals(forest_growing_state_t& forest_state,
            const VectorView<Minimizer>& minimizers, const interval_state_t& snarl_interval,
            const interval_state_t& parent_interval,
            const vector<interval_state_t>& child_intervals) const;

    //////////////////////////////////////////////////////
    ///////////          functions for building the trees
    /////////////////////////////////////////////////////

    // Open a chain that starts at the current_seed
    // If the chain is in a snarl, then add empty edges for the distances to everything before it 
    // in the snarl (found with sibling_indices_at_depth) 
    // Open the chain, and record its presence and distance-to-start in the parent snarl, if 
    // necessary seed_index is the index into seeds of the first seed in the chain
    void open_chain(forest_growing_state_t& forest_state, const size_t& depth, 
                    size_t seed_index, bool chain_is_reversed);

    // Close a chain that ends at last_seed
    // If the chain was empty, remove it and anything relating to it in the parent snarl and 
    // sibling_indices
    // If it can be spliced out, take out a subtree
    // Otherwise, add the end of the chain and, if the chain was in a snarl, add the distances to 
    // everything before it in the snarl and remember the distance to the end of the chain
    void close_chain(forest_growing_state_t& forest_state, const size_t& depth, 
                     const Seed& last_seed, bool chain_is_reversed);

    // Add the current seed (or snarl starting at the seed) and its distance to the previous thing 
    // in a chain
    // If the seed is far enough from the previous thing in the chain and it can be a new slice, 
    // split off a subtree
    // depth is the depth of the child of the chain (which may also be the chain depth if it 
    // is trivial)
    // seed_index is the index of the current seed in the list of seeds
    void add_child_to_chain(forest_growing_state_t& forest_state,
                      const size_t& depth, const size_t& seed_index, 
                      bool child_is_reversed, bool chain_is_reversed);

    // Start a new snarl
    void open_snarl(forest_growing_state_t& forest_state, const size_t& depth);

    // Close a snarl
    // depth is the depth of the snarl and last_seed is the last seed in the snarl
    // If the snarl has no children, then delete the whole thing
    // Otherwise, add all necessary distances and close it
    void close_snarl(forest_growing_state_t& forest_state, const size_t& depth, 
                     const Seed& last_seed, bool last_is_reversed, bool is_cyclic_snarl);

    // Add all the distances from everything in the snarl to either the last child of the snarl or,
    // if to_snarl_end is true, to the end bound of the snarl
    // depth is the depth of the snarl
    void add_snarl_distances(forest_growing_state_t& forest_state,
                             const size_t& depth, const Seed& seed, bool child_is_reversed, 
                             bool snarl_is_reversed, 
                             bool to_snarl_end, bool is_cyclic_snarl);


    /// Given a vector of value pairs, and a bool indicating if the pair is used for the correlation,
    /// return the correlation. This is the spearman correlation for now
    static double get_correlation (const vector<std::pair<size_t, size_t>>& values);



    /************ Helper functions for debugging ************/


    public:

    void print_self(const vector<Seed>* seeds) const {
        for (size_t i = 0 ; i < trees.size() ; i++) {
            const auto& tree = trees[i];
            cerr << i << ": ";
            tree.print_self(seeds);
        }
    }
    void validate_zip_forest(const SnarlDistanceIndex& distance_index, 
                             const vector<Seed>* seeds,
                             size_t distance_limit=std::numeric_limits<size_t>::max()) const;


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






















/// Implementations for the templated functions using Minimizers since the definition is in the minimizer_mapper
//TODO: This really shouldn't be in the hpp file

namespace vg {
    using namespace std;

template<typename Minimizer>
void ZipCodeForest::fill_in_forest(const vector<Seed>& seeds, const VectorView<Minimizer>& minimizers,
                                   const SnarlDistanceIndex& distance_index, size_t gap_distance_limit,
                                   size_t distance_limit) {
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Make a new forest with " << seeds.size() << " seeds with distance limit " << distance_limit << endl;
    for (auto& x : seeds) {
        cerr << x.pos << endl;
    }
    cerr << endl;
#endif
    if (seeds.size() == 0) {
        return;
    }

    /*
    The zip forest is made by sorting the seeds along chains/snarls, then adding each seed, 
    snarl/chain boundary, and distance to zip_code_tree.

    Sorting and tree-making are done at the same time, in a depth-first traversal of the snarl tree.
    Sorting is done per node in the snarl tree.
    
    Intervals representing ranges of seeds corresponding to snarl tree structures are stored in a 
    stack. The algorithm starts with an interval for each child of the root snarl. An interval is
    popped from the stack. Any incomplete snarls or chains that the interval is not a child of 
    must be completed. Then, the snarl or chain that the interval represents is added to the zip 
    tree, along with any relevant distances. Intervals representing the children of the snarl or 
    chain are found and added to the stack. This repeats until the stack is empty.

    */

    //Start by initializing the state
    //The forest state keeps track of the sort order of seeds, the intervals that need to be sorted,
    //and which intervals are open and incomplete. 
    forest_growing_state_t forest_state(seeds, distance_index, gap_distance_limit, distance_limit);

    //Start with the root as the interval over seed_sort_order containing everything
    interval_state_t first_interval (0, seeds.size(), false, ZipCode::EMPTY, 0);

    //Sort and get the intervals of the connected components
    sort_one_interval(forest_state, first_interval);
    vector<interval_state_t> new_intervals 
            = get_next_intervals(forest_state, first_interval);
    forest_state.intervals_to_process.insert(forest_state.intervals_to_process.end(),
                                             std::make_move_iterator(new_intervals.rbegin()),
                                             std::make_move_iterator(new_intervals.rend()));


    while (!forest_state.intervals_to_process.empty()) {
#ifdef DEBUG_ZIP_CODE_TREE
        print_self(&seeds);
#endif
        // For each unprocessed interval, process it
        // First, check if anything needs to be closed, which will happen if the interval's depth is 
        //   greater than or equal to that of an open interval.
        //   Distances between snarl children are added after the child is closed.
        // Get the intervals of this interval's children and add them in reverse order to the stack
        //   intervals_to_process
        // Open the current interval's snarl/chain


        //Get the interval
        interval_state_t current_interval = std::move(forest_state.intervals_to_process.back());
        forest_state.intervals_to_process.pop_back();

        /********************

         * First, check if anything needs to be closed and close it

         ************************/

#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "Process interval of type " << current_interval.code_type << " with range " 
             << current_interval.interval_start << "-" << current_interval.interval_end << endl;
        assert(current_interval.depth <= 
               seeds.at(forest_state.seed_sort_order[current_interval.interval_start]).zipcode_decoder->max_depth()+1);
        cerr << "Close anything open" << endl;
#endif
        while (!forest_state.open_intervals.empty()) {
            if (current_interval.depth <= forest_state.open_intervals.back().depth) {
                //If the current interval is not a child of the open interval
                //close the last thing in open_intervals
                //There will be an interval for every ancestor in the snarl tree, so this can just check depth

#ifdef DEBUG_ZIP_CODE_TREE
                cerr << "\tclose something at depth " << forest_state.open_intervals.size()-1 << endl;
#endif

                size_t depth = forest_state.open_intervals.size()-1;

                //The ancestor interval to close and its last seed
                const interval_state_t& ancestor_interval = forest_state.open_intervals.back();
                const Seed& last_seed = seeds.at(forest_state.seed_sort_order[ancestor_interval.interval_end-1]);

                if (ancestor_interval.code_type == ZipCode::CHAIN ||
                    ancestor_interval.code_type == ZipCode::NODE ||
                    ancestor_interval.code_type == ZipCode::ROOT_CHAIN ||
                    ancestor_interval.code_type == ZipCode::ROOT_NODE) {
                    //Close a chain

                    close_chain(forest_state, depth, 
                                last_seed, ancestor_interval.is_reversed); 
                } else {
#ifdef DEBUG_ZIP_CODE_TREE
                    assert(ancestor_interval.code_type == ZipCode::REGULAR_SNARL ||
                           ancestor_interval.code_type == ZipCode::IRREGULAR_SNARL ||
                           ancestor_interval.code_type == ZipCode::CYCLIC_SNARL ||
                           ancestor_interval.code_type == ZipCode::ROOT_SNARL);
#endif
                    //Close a snarl
                    close_snarl(forest_state, depth, last_seed, 
                                ancestor_interval.is_reversed, ancestor_interval.code_type == ZipCode::CYCLIC_SNARL); 
                }

                //Clear the list of children of the snarl tree structure at this level
                forest_state.sibling_indices_at_depth[depth].clear();

                //Take out this ancestor
                forest_state.open_intervals.pop_back();
            } else {
                //If the current interval is contained in this open interval, then it is also contained in all other
                // ancestors so break
                break;
            }
        }

        /************ 
         * Now start processing the current interval
         *
         *
         * Sort this interval and add the child intervals in reverse order to intervals_to_process  
         ***********/

        
        //For everything except non-dag snarls, sort get the intervals normally

        if (current_interval.code_type != ZipCode::NODE ) {
            //Sort the current interval and get the intervals corresponding to its children
            sort_one_interval(forest_state, current_interval);
            vector<interval_state_t> child_intervals = get_next_intervals(forest_state, current_interval);
            if (current_interval.code_type != ZipCode::CYCLIC_SNARL || current_interval.is_reverse_ordered
                    || current_interval.is_ordered){ 

                //If this is not a cyclic snarl, or it is the duplicated copy of a cyclic snarl child
                //Add the child intervals to the to_process stack, in reverse order so the first one
                //gets popped first
                //By forcing duplicated copies of a cyclic snarl child to be processed here, we 
                //prevent nested cyclic snarls from being duplicated in each copy, preventing an 
                //exponential blowup
                forest_state.intervals_to_process.insert(forest_state.intervals_to_process.end(),
                                                         std::make_move_iterator(child_intervals.rbegin()),
                                                         std::make_move_iterator(child_intervals.rend()));
            } else {
                //If this is a cyclic snarl, then we do further partitioning before adding the child intervals
                //The new intervals may include duplicates, so we want to limit how many times this happens

                vector<interval_state_t> snarl_child_intervals = get_cyclic_snarl_intervals(
                                                                              forest_state, 
                                                                              minimizers,
                                                                              current_interval,
                                                                              forest_state.open_intervals.back(),
                                                                              child_intervals);

                forest_state.intervals_to_process.insert(forest_state.intervals_to_process.end(),
                                                         std::make_move_iterator(snarl_child_intervals.rbegin()),
                                                         std::make_move_iterator(snarl_child_intervals.rend()));
            }
        }
    
    
        /**********
         *
         * Open the current interval
         * If the current interval is a snarl and a child of a chain, then add the preceding sibling seeds before the snarl
         *
         *******/

#ifdef DEBUG_ZIP_CODE_TREE
         cerr << "Open next interval or (if the interval is for nodes), add seeds" << endl;
#endif
        if (forest_state.open_intervals.size()+1 > forest_state.sibling_indices_at_depth.size()) {
            forest_state.sibling_indices_at_depth.emplace_back();
        }
        if (forest_state.open_intervals.empty()) {
            // If there is nothing open, then this is starting a new connected component
            // Just open it

#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "Start a new connected component" << endl;
            assert(current_interval.code_type == ZipCode::ROOT_NODE ||
                   current_interval.code_type == ZipCode::NODE ||
                   current_interval.code_type == ZipCode::ROOT_CHAIN ||
                   current_interval.code_type == ZipCode::ROOT_SNARL);
#endif

            if (forest_state.active_zip_tree_i == std::numeric_limits<size_t>::max() 
                || trees[forest_state.active_zip_tree_i].zip_code_tree.size() != 0) {
                trees.emplace_back();
                forest_state.active_zip_tree_i = trees.size()-1;
            }

            if (current_interval.code_type == ZipCode::ROOT_SNARL) {
                // Open the root snarl
                open_snarl(forest_state, 0);
            } else if (current_interval.code_type == ZipCode::NODE) {
                //For a root node, just add it as a chain with all the seeds

                trees[forest_state.active_zip_tree_i].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_START, 
                                                                             std::numeric_limits<size_t>::max(), 
                                                                             false);

                //Remember the start of the chain
                forest_state.sibling_indices_at_depth[0].push_back({ZipCodeTree::CHAIN_START, 0});

                //If this is a node, then the interval contains everything in it, so add the seeds and close the chain here
                for (size_t seed_i = current_interval.interval_start ; seed_i < current_interval.interval_end ; seed_i++) {

                    add_child_to_chain(forest_state, current_interval.depth, 
                                       forest_state.seed_sort_order[seed_i], current_interval.is_reversed,
                                       current_interval.is_reversed); 
                }
                close_chain(forest_state, current_interval.depth, 
                            seeds.at(forest_state.seed_sort_order[current_interval.interval_end-1]), 
                            current_interval.is_reversed); 

                
            } else {
                // Open the root chain/node
                trees[forest_state.active_zip_tree_i].zip_code_tree.emplace_back(ZipCodeTree::CHAIN_START, 
                                                                             std::numeric_limits<size_t>::max(), 
                                                                             false);

                //Remember the start of the chain
                forest_state.sibling_indices_at_depth[0].push_back({ZipCodeTree::CHAIN_START, 0});
            }                       
        } else if (forest_state.open_intervals.back().code_type == ZipCode::CHAIN || 
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_CHAIN ||
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_NODE) {
            // This is the child of a chain

            if (current_interval.code_type == ZipCode::NODE) {
                // If the type of this interval is NODE, then this is a range of seeds that are on nodes on the chain,
                // not necessarily on the same node
                // Add each seed

                bool is_trivial_chain = current_interval.depth-1 == 
                        seeds.at(forest_state.seed_sort_order[current_interval.interval_start]).zipcode_decoder->max_depth();
                for (size_t seed_i = current_interval.interval_start ; seed_i < current_interval.interval_end ; seed_i++) {


                    add_child_to_chain(forest_state, is_trivial_chain ? current_interval.depth-1 : current_interval.depth, 
                                       forest_state.seed_sort_order[seed_i], current_interval.is_reversed,
                                       forest_state.open_intervals.back().is_reversed); 
                }

            } else {
#ifdef DEBUG_ZIP_CODE_TREE
                assert(current_interval.code_type == ZipCode::REGULAR_SNARL || 
                       current_interval.code_type == ZipCode::IRREGULAR_SNARL || 
                       current_interval.code_type == ZipCode::CYCLIC_SNARL);
#endif

                //Add the snarl to the chain
                add_child_to_chain(forest_state, current_interval.depth,
                                   forest_state.seed_sort_order[current_interval.interval_start], 
                                   current_interval.is_reversed, forest_state.open_intervals.back().is_reversed);
            }
            

        } else {
        //If there is an open ancestor that isn't a chain, so the ancestor must be a snarl
#ifdef DEBUG_ZIP_CODE_TREE
            assert(forest_state.open_intervals.back().code_type == ZipCode::REGULAR_SNARL ||
                   forest_state.open_intervals.back().code_type == ZipCode::IRREGULAR_SNARL ||
                   forest_state.open_intervals.back().code_type == ZipCode::CYCLIC_SNARL ||
                   forest_state.open_intervals.back().code_type == ZipCode::ROOT_SNARL);
#endif

            //Open the child chain
            open_chain(forest_state, forest_state.open_intervals.size(), 
                       forest_state.seed_sort_order[current_interval.interval_start], current_interval.is_reversed);
            
        }

        if (current_interval.code_type != ZipCode::NODE) {
            // Add to open_intervals
            forest_state.open_intervals.emplace_back(std::move(current_interval));
        }
    }
    //Finished adding all intervals
    

    //Now close anything that remained open
    while (!forest_state.open_intervals.empty()) {
        interval_state_t& ancestor_interval = forest_state.open_intervals.back();
        const Seed& last_seed = seeds.at(forest_state.seed_sort_order[ancestor_interval.interval_end-1]);

        if (ancestor_interval.code_type == ZipCode::CHAIN ||
            ancestor_interval.code_type == ZipCode::ROOT_CHAIN ||
            ancestor_interval.code_type == ZipCode::ROOT_NODE) {
            //Close a chain

            close_chain(forest_state, forest_state.open_intervals.size()-1, 
                        last_seed, ancestor_interval.is_reversed); 
        } else {
#ifdef DEBUG_ZIP_CODE_TREE
            assert(ancestor_interval.code_type == ZipCode::REGULAR_SNARL ||
                   ancestor_interval.code_type == ZipCode::IRREGULAR_SNARL ||
                   ancestor_interval.code_type == ZipCode::CYCLIC_SNARL ||
                   ancestor_interval.code_type == ZipCode::ROOT_SNARL);
#endif
            //Close a snarl
            close_snarl(forest_state, forest_state.open_intervals.size()-1, 
                        last_seed, ancestor_interval.is_reversed, ancestor_interval.code_type == ZipCode::CYCLIC_SNARL); 
        }

        forest_state.open_intervals.pop_back();
    }

    if (trees[forest_state.active_zip_tree_i].zip_code_tree.size() == 0) {
        trees.erase(trees.begin() + forest_state.active_zip_tree_i);
    }
#ifdef DEBUG_ZIP_CODE_TREE
    print_self(&seeds);
    validate_zip_forest(distance_index, &seeds, distance_limit);
    assert(forest_state.open_chains.empty());
    assert(forest_state.open_intervals.empty());
#endif

}

template<typename Minimizer>
vector<ZipCodeForest::interval_state_t> ZipCodeForest::get_cyclic_snarl_intervals(
    forest_growing_state_t& forest_state,
    const VectorView<Minimizer>& minimizers, const interval_state_t& snarl_interval,
    const interval_state_t& parent_interval,
    const vector<interval_state_t>& child_intervals) const {

    vector<size_t>& zipcode_sort_order = forest_state.seed_sort_order;
    vector<sort_value_t>& sort_values_by_seed = forest_state.sort_values_by_seed;
    const vector<Seed>* seeds = forest_state.seeds; 
    const SnarlDistanceIndex* distance_index = forest_state.distance_index;

#ifdef DEBUG_ZIP_CODE_TREE
    assert(seeds->at(zipcode_sort_order[snarl_interval.interval_start]).zipcode_decoder->get_code_type(snarl_interval.depth) 
                == ZipCode::CYCLIC_SNARL);
    net_handle_t handle = seeds->at(zipcode_sort_order[snarl_interval.interval_start]).zipcode_decoder->get_net_handle(snarl_interval.depth, distance_index);
    cerr << "Sorting and finding intervals for cyclic snarl " << distance_index->net_handle_as_string(handle)
         << " with " << child_intervals.size() << " children" << endl;
#endif

    net_handle_t snarl_handle = seeds->at(zipcode_sort_order[snarl_interval.interval_start]).zipcode_decoder->get_net_handle(snarl_interval.depth, distance_index);


    /****** For each interval, form runs of reachable seeds 
      seeds are reachable if they are close on the read and chain (by distance to start of chain)
      and if they are on the same strand on the read                                              ***********/


    //A union find for finding runs of seeds that are reachable in the read and chain
    structures::UnionFind union_find(snarl_interval.interval_end - snarl_interval.interval_start) ;

    // Define a struct that represents a run
    // runs get merged with each other if they are close enough by checking the ranges they cover
    // in the read and chain
    struct run_t {
        // The representative seed in the union find
        // This is also an index into zipcode_sort_order if you add snarl_interval.interval_start
        size_t uf_head; 

        //The range of positions in the read spanned by the seeds in this run
        size_t read_range_start;
        size_t read_range_end;

        //The same thing but for the chain
        size_t chain_range_start;
        size_t chain_range_end;

        //The index of the original interval in child_intervals
        size_t interval_i;

        bool is_reversed_read;

        //Can this interval be traversed in both directions?
        bool can_be_reversed;
    };

    //Helper function to check if the value is close enough to a range of values
    auto is_within_range = [&] (size_t range_start, size_t range_end, size_t value) {
        if (value >= range_start && value <= range_end) {
            //If the value is inside the range
            return true;
        } else if (value < range_start && range_start - value <= forest_state.gap_distance_limit) {
            //If the value is before the range but still within the distance limit
            return true;
        } else if (value > range_end && value - range_end <= forest_state.gap_distance_limit) {
            //If the value is after the range but still within the distance limit
            return true;
        } else {
            return false;
        }
    };


    /*************

      Figure out the orientation of the read through the snarl

    ************/

    //Get pairs of read/chain offsets along the parent chain
    vector<pair<size_t, size_t>> parent_offset_values;

    //Check up to this many seeds on the parent chain
    size_t check_count = 50;
    int check_i = snarl_interval.interval_start - 1;

    //Get up to half of the values from before the snarl
    while (check_i >= parent_interval.interval_start && parent_offset_values.size() <= check_count/2) {

        if (seeds->at(zipcode_sort_order[check_i]).zipcode_decoder->max_depth() == snarl_interval.depth) {
            parent_offset_values.emplace_back(minimizers[seeds->at(zipcode_sort_order[check_i]).source].value.offset,
                                              seeds->at(zipcode_sort_order[check_i]).zipcode_decoder->get_offset_in_chain(snarl_interval.depth));
        }

        check_i--;
    }

    //Get the rest from after the snarl

    check_i = snarl_interval.interval_end;
    while (check_i < parent_interval.interval_end && parent_offset_values.size() < check_count) {

        if (seeds->at(zipcode_sort_order[check_i]).zipcode_decoder->max_depth() == snarl_interval.depth) {
            parent_offset_values.emplace_back(minimizers[seeds->at(zipcode_sort_order[check_i]).source].value.offset,
                                              seeds->at(zipcode_sort_order[check_i]).zipcode_decoder->get_offset_in_chain(snarl_interval.depth));
        }

        check_i++;
    }

    //>0 if the read flows backwards through the snarl
    double parent_correlation = get_correlation(parent_offset_values);
#ifdef DEBUG_ZIP_CODE_TREE
    cerr << "Correlation of parent chain from " << parent_offset_values.size() << " value pairs: " 
         << parent_correlation << endl;
#endif

    /*******************

      For each child of the snarl, walk through the seeds and build runs of seeds that are close
      For each seed, compare it to all other seeds found so far to see if they can be merged

      *****************/


    forward_list<run_t> all_runs;
    //For each seed, remember its offset in the read and chain to later compute the correlation
    vector<std::tuple<size_t, size_t, bool>> read_and_chain_offsets (snarl_interval.interval_end-snarl_interval.interval_start);

    for (size_t interval_i = 0 ; interval_i < child_intervals.size() ; interval_i++) {
        const auto& child_interval = child_intervals[interval_i];

        //Each interval is on one chain, but the chains aren't sorted yet so sort them
        sort_one_interval(forest_state, child_interval);

        //Check if the interval can be flipped in the snarl
        bool interval_is_reversed_in_snarl = child_interval.is_reversed != snarl_interval.is_reversed;
        bool interval_is_reversable;
        if (interval_is_reversed_in_snarl) {
            //If this interval is already going backwards in the snarl, then it is because it couldn't go forwards

#ifdef DEBUG_ZIP_CODE_TREE
            //This is how seed_is_reversed_at_depth currently works but double check this in case it changed 
            size_t rank = seeds->at(zipcode_sort_order[child_interval.interval_start]).zipcode_decoder->get_rank_in_snarl(snarl_interval.depth+1);
            assert (distance_index->distance_in_snarl(snarl_handle, 0, false, rank, false) == std::numeric_limits<size_t>::max()
                &&
                distance_index->distance_in_snarl(snarl_handle, 1, false, rank, true) == std::numeric_limits<size_t>::max());
#endif

            interval_is_reversable = false;
        } else {
            //If the interval is not reversed in the snarl, check if it can be reversed
            size_t rank = seeds->at(zipcode_sort_order[child_interval.interval_start]).zipcode_decoder->get_rank_in_snarl(snarl_interval.depth+1);
            size_t distance_start = distance_index->distance_in_snarl(snarl_handle, 0, false, rank, true);
            size_t distance_end = distance_index->distance_in_snarl(snarl_handle, 1, false, rank, false);
            interval_is_reversable = distance_start != std::numeric_limits<size_t>::max()
                               || distance_end != std::numeric_limits<size_t>::max();
        }
                                

        //Now partition the chain further

        //This is the set of runs for this particular chain
        std::forward_list<run_t> runs;


        //Go through all seeds in the chain and compare them to the open runs.
        //Add the seed to any run that it is reachable with, potentially combining runs
        for (size_t sort_i = child_interval.interval_start ; sort_i < child_interval.interval_end ; sort_i++) {
            const Seed& seed = seeds->at(zipcode_sort_order[sort_i]);
            const Minimizer& minimizer = minimizers[seed.source];

            //The relevant values for checking this seed against an existing run
            bool is_reversed_read = minimizer.value.is_reverse;
            size_t read_offset = minimizer.value.offset;
            size_t chain_offset = sort_values_by_seed[zipcode_sort_order[sort_i]].get_distance_value(); 

            //Remember the values for finding the correlation later
            std::get<0>(read_and_chain_offsets [sort_i-snarl_interval.interval_start])= read_offset;
            std::get<1>(read_and_chain_offsets [sort_i-snarl_interval.interval_start]) = 
                    sort_values_by_seed[zipcode_sort_order[sort_i]].get_sort_value();
            std::get<2>(read_and_chain_offsets [sort_i-snarl_interval.interval_start]) =
                    seed.zipcode_decoder->max_depth() <= snarl_interval.depth+2;


            //Make a new run for the seed, to be updated with anything combined with it
            run_t seed_run({sort_i - snarl_interval.interval_start,
                            read_offset, read_offset,
                            chain_offset, chain_offset,
                            interval_i,
                            is_reversed_read,
                            interval_is_reversable});

            //For each run, check if it is reachable with the seed, and remove the ones that aren't

            //To remove an element, keep track of the element (run_itr) and the previous iterator (prev_itr),
            // and remove_after the previous iterator
            auto prev_itr = runs.before_begin();
            auto run_itr = runs.begin();
            while (run_itr != runs.end()) {

                //A seed is reachable with a run if they are both on the same strand on the read,
                //the seed is close enough in the read, and if the seed is close enough in the chain 

                if (is_reversed_read == run_itr->is_reversed_read &&
                    is_within_range(run_itr->read_range_start, run_itr->read_range_end, read_offset) &&
                    is_within_range(run_itr->chain_range_start, run_itr->chain_range_end, chain_offset)) {
                    //If this run is reachable with the seed

                    //Combine the runs
                    seed_run.uf_head = union_find.union_groups(run_itr->uf_head, 
                                                               seed_run.uf_head);
                    seed_run.read_range_start = std::min(run_itr->read_range_start, 
                                                         seed_run.read_range_start);
                    seed_run.read_range_end = std::max(run_itr->read_range_end, 
                                                       seed_run.read_range_end);

                    seed_run.chain_range_start = std::min(run_itr->chain_range_start, 
                                                          seed_run.chain_range_start);
                    seed_run.chain_range_end = std::max(run_itr->chain_range_end, 
                                                        seed_run.chain_range_end);

                    //Remove this run
                    run_itr = runs.erase_after(prev_itr);
                } else {
                    //Otherwise, iterate to the new run
                    ++run_itr;
                    ++prev_itr;
                }
            }
            //Add the new run
            runs.push_front(std::move(seed_run));
            //TODO: Remove runs that are definitely too far away from anything else
        }
#ifdef DEBUG_ZIP_CODE_TREE
        cerr << "\tnew runs:" << endl;
        for (auto& run : runs) {
            auto seed_is = union_find.group(run.uf_head);
            for (size_t i : seed_is) {
                cerr << seeds->at(zipcode_sort_order[snarl_interval.interval_start+i]).pos << ", ";
            }
            cerr << "|";
        }
        cerr << endl;
#endif
        //Add this chain's runs to the overall list
        //This merging combines two sorted lists so sort first
        runs.sort([&](const run_t& a, const run_t& b) {
            if (parent_correlation < 0.0) {
               //If the read is going backwards through the snarl, then sort backwards by the first read coordinate 
                return a.read_range_start > b.read_range_start;
            } else {
                //Otherwise, sort so the last read coordinates go forwards
                return a.read_range_end < b.read_range_end;
            }
        });
        all_runs.merge(runs, [&](const run_t& a, const run_t& b) {
            if (parent_correlation < 0.0) {
               //If the read is going backwards through the snarl, then sort backwards by the first read coordinate 
                return a.read_range_start > b.read_range_start;
            } else {
                //Otherwise, sort so the last read coordinates go forwards
                return a.read_range_end < b.read_range_end;
            }
            });
    }
    //TODO: Merge consecutive runs on the same chain. This shouldn't affect correctness because separate 
    //      should be unreachable, but it would make the snarls smaller

    /******* Re-sort seeds by the new runs and make new intervals of the runs on the chains 
        The orientation of the runs is determined by the orientation of the read along the parent chain  ***********/
    

    vector<ZipCodeForest::interval_state_t> new_intervals;
    //New sort order to replace what's currently in zipcode_sort_order for this snarl 
    vector<size_t> new_sort_order;
    new_sort_order.reserve(snarl_interval.interval_end - snarl_interval.interval_start);

    for (const run_t& run : all_runs) {
        //For each run, add its seeds to the sort order
        //The seeds are already in the correct sort order for the chain in zipcode_sort_order, so
        //re-sort the run's seeds according to this order
        //Also check if the orientation of the read is backwards relative to the snarl, and if so,
        //flip the order of the run so it gets traversed backwards

        vector<size_t> run_seeds = union_find.group(run.uf_head);
        std::sort(run_seeds.begin(), run_seeds.end());

        new_intervals.emplace_back(snarl_interval.interval_start + new_sort_order.size(),
                                    snarl_interval.interval_start + new_sort_order.size() + run_seeds.size(),   
                                    child_intervals[run.interval_i].is_reversed,
                                    child_intervals[run.interval_i].code_type,
                                    child_intervals[run.interval_i].depth);

        //Figure out if the read running backwards through this run
        bool reverse_run = false;
        //Should we use both orientations?
        bool duplicate_run = false;
        
        if (run.can_be_reversed && parent_offset_values.size() > 0) {
            //If it is possible to traverse the run backwards in the chain, then check which is the correct orientation
            vector<pair<size_t, size_t>> run_values;
            run_values.reserve(run_seeds.size());
            for (size_t x : run_seeds) {
                if (std::get<2>(read_and_chain_offsets[x])){
                    run_values.emplace_back(std::get<0>(read_and_chain_offsets[x]),
                                            std::get<1>(read_and_chain_offsets[x]));
                }
            }

            double run_correlation = get_correlation(run_values);
#ifdef DEBUG_ZIP_CODE_TREE
            cerr << "Correlation of child run from " << run_values.size() << " value pairs: " 
                 << run_correlation << endl;
#endif
            if (std::abs(run_correlation) < 0.8 || std::abs(parent_correlation) < 0.6) {
                //If the correlation is too low, then just duplicate the run in both orientations
                //TODO This is very arbitrary, especially for the parent correlation
                duplicate_run = true;
            } else {

                bool snarl_is_traversed_backwards =  parent_correlation < 0.0;
                //If the parent chain is backwards, then the orientation gets flipped
                // This is necessary because the values used to get the correlation were the actual
                // prefix sums, not the order they were traversed in
                if (parent_interval.is_reversed) {
                    snarl_is_traversed_backwards = !snarl_is_traversed_backwards;
                }

                //Now decide which direction the run is traversed in
                bool run_is_traversed_backwards = run_correlation < 0.0;
                reverse_run = run_is_traversed_backwards != snarl_is_traversed_backwards;
            }

        }

        if (!reverse_run) {
            //If we can only go forwards through the run or
            //if the read is going through the snarl and partition in the same direction
            for (size_t sort_i : run_seeds) {
                new_sort_order.push_back(zipcode_sort_order[snarl_interval.interval_start+sort_i]);
            }

            //If we're also duplicating this run, add another interval for the same thing reversed
            if (duplicate_run) {
                const auto& last_interval = new_intervals.back();
                new_intervals.emplace_back(last_interval.interval_start,
                                           last_interval.interval_end,
                                           !last_interval.is_reversed,
                                           last_interval.code_type,
                                           last_interval.depth);
                //Remember to reverse the order
                new_intervals.back().is_reverse_ordered=true;
            }

        } else {
            //If the read is going through the run in the opposite direction as the snarl, then flip it
            for (int i = run_seeds.size()-1 ; i >= 0 ; --i) {
                new_sort_order.push_back(zipcode_sort_order[snarl_interval.interval_start+run_seeds[i]]);
            }
            new_intervals.back().is_reversed = !new_intervals.back().is_reversed;
        }
    }

    //Update the sort order in zipcode_sort_order
    for (size_t i = 0 ; i < new_sort_order.size() ; i++) {
        zipcode_sort_order[snarl_interval.interval_start+i] = new_sort_order[i];
    }
#ifdef DEBUG_ZIP_CODE_SORTING
    assert(new_sort_order.size() == (snarl_interval.interval_end - snarl_interval.interval_start));
    cerr << "New sort order " << endl;
    for (auto& interval : new_intervals) {
        for (size_t i = interval.interval_start ; i < interval.interval_end ; i++) {
            cerr << seeds->at(zipcode_sort_order[i]).pos << ", ";
        }
        cerr << "|";
    }
    cerr << endl;
#endif

    return new_intervals;
}
}

#endif

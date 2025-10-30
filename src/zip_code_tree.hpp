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

/*

A ZipCodeTree represents of set of SnarlDistanceIndexCluserer::Seed's
(seed alignments between a read and reference) as a tree structure.
The tree represents the connectivity of the seeds, based on the distance index.
Edges are labelled with distance values.

Two iterators are provided:
- A seed_iterator yields seeds in the tree, in left-to-right order
- A distance_iterator yields seeds within a given distance limit
reachable from a given seed, traversing the tree with a given initial direction

A ZipCodeForest, which represents a collection of trees, builds ZipCodeTree's
via its fill_in_forest(), which accepts a vector of seeds, a graph's distance
index, and an optional distance limit. This is meant to be called using the
seeds found in a read.

Main developers:
- Xian Chang (@xchang1)
- Adam Novak (@adamnovak): distance iterator
- Faith Okamoto (@faithokamoto): modifications for better cyclic snarl handing

*/
class ZipCodeTree {

    /// Convenient alias for SnarlDistanceIndexClusterer::Seed
    /// Despite the name, these are used for graph positions,
    /// so they act more like minimizers
    typedef SnarlDistanceIndexClusterer::Seed Seed;

    public:

    /// Empty constructor
    /// ZipCodeTree's get filled in by ZipCodeForest's
    ZipCodeTree(){};

    /*
      The tree will represent the seeds' placement in the snarl tree.
      Each node in the tree represents either a seed (position on the graph,
      for the start of an alignment) or the boundary of a snarl or chain.
      Edges are labelled with the distance between the two nodes

      This graph is actually represented as a vector of the nodes and edges
      Each item in the vector represents either a node (seed or boundary),
      an edge (distance), or the child count of a snarl

      ---- CHAINS ----

      A chain in the vector is bounded by a CHAIN_START and a CHAIN_END.
      A chain is comprised of alternating children (seed or snarl) and the
      distances between them, starting and ending with a child. The order is:
        CHAIN_START, child, distance, child, distance, ..., child, CHAIN_END
      The distance from the chain start to the first child is
      included in the distances in the chain's parent snarl, if relevant

      The distances represent the number of nucleotides on the minimum-length
      path in the variation graph between the structures that the zip code tree
      nodes represent. For edges within a chain, offsets are queried from the
      seeds' zipcodes.

      Seeds represent the first nucleotide of the alignment. When a seed is
      traversed forwards in the zip tree (get_is_reversed() is false), distancea
      starting from that seed includes the position. If the seed is reversed in
      the zip tree, then distances doesn't include the position.

      For two SEEDs on the same position, the distance between them would be 0.
      As a convenience, such edges are omitted; the tree will have [seed seed]
      For chain distances terminating at a snarl bound, the distance reaches the
      inner edge (relative to the snarl) of the boundary node,
      so it includes the length of the boundary node of the snarl

      For example, given a subgraph of a chain:

                              n3
                            [GACG] ...
            n1      n2    / 
            [A] - [AGAC]
                          \   n4
                            [ACAG] ...

      for the sequence "SEED EDGE SNARL_START" representing a seed on n1
      and the snarl starting at n2, the edge value would be 5
      Within the snarl, edges include the distance to the first seed in a chain.
      For a seed at position node 3 +1 (the A oriented forwards), the distance
      matrix would store an edge value of 1 for the distance between the
      SNARL_START and that chain

      ---- SNARLS ----

      A snarl in the vector is bounded by a SNARL_START and a SNARL_END.
      A snarl is comprised of the two bounds, one or more chains,
      and the distances among them. SEEDs are always contained within a chain.
      Immediately after the SNARL_START, there is a CHAIN_COUNT, then the lower
      triangle of a distance matrix made using distance_in_snarl(). Distances
      are always FROM the right side of a chain/bound TO the left side.
      A DAG snarl with two chains would look like:
          SNARL_START, CHAIN_COUNT, dist:start->c1, dist:start->c2,
                dist:c1->c2, dist:start->end, dist:c1->end, dist:c2->end,
                chain1, chain2, SNARL_END

      For snarls that aren't DAGs (called cyclic snarls, though they could have
      an inversion and no cycles), distances are stored to either end of each
      chain. Self-loops are included. Self-loops for chain bounds are stored as
      inf to keep the triangle shape but are never used.
      A cyclic snarl with one chain would look like:
          SNARL_START, CHAIN_COUNT, inf, dist:start->c1_L,
                dist:c1_L->c1_L, dist:start->c1_R, dist:c1_L->c1_R,
                dist:c1_R->c1_R, dist:start->end, dist:c1_L->end,
                dist:c1_R->end, inf, chain1, SNARL_END
    
      ---- ORDERING ----

      Everything is ordered according to the order of the highest-level chain
      (top-level chain or child of a top-level snarl).
      Children of a snarl are ordered by a topological sort of the snarl.

      In the variation graph, all chains are considered to be oriented "forward"
      in their parent snarl. However, in a start-to-end traversal of the snarl,
      the child chain may be traversed end-to-start. These chains are considered
      to be reversed in the zip code tree, so the order of the children of the
      chain may be backwards relative to their order in the variation graph.
      If a snarl is the child of a chain that is traversed backwards in the zip 
      tree, then that snarl and all its children are also traversed backwards.

     */

    public:

    /// The type of an item in the zip code tree
    enum tree_item_type_t {SEED=0, CHAIN_START, CHAIN_END, EDGE, CHAIN_COUNT, SNARL_START, SNARL_END};

    /// One item in the zip code tree, representing a node or edge of the tree
    struct tree_item_t {

        private:
        /// Is this a seed, boundary, or an edge
        tree_item_type_t type : 4;

        /// For a seed, the index into seeds
        /// For an edge, the distance value
        /// For a snarl or child chain bound, distance to the snarl start
        /// For example, (1 0 0 0 [seed]) would have 0 for the snarl start,
        /// 5 for the chain start, 7 for the chain end, and 8 for the snarl end
        size_t value : 59;

        /// For a bound, how long the internal section is
        /// e.g. a chain with one seed would be "1"
        size_t section_length : 59;

        /// For a seed, if we're walking through the tree from right to left
        /// (the default), will we traverse this position backwards?
        /// Or, for a bound is the snarl or parent snarl cyclic?
        /// Ignored for EDGE/CHAIN_COUNT and should be set to false
        bool is_reversed_or_cyclic;

        /// The maximum value of internal fields;
        /// this corresponds to std::numeric_limits<size_t>::max() outside
        static size_t internal_max() {
            return ((size_t)1 << 59) - 1;
        }

        public:

        /// Empty constructor
        tree_item_t (){};

        /// Constructor so that value gets set properly
        tree_item_t (tree_item_type_t type, size_t raw_value, bool is_reversed_or_cyclic) 
            : type(type), is_reversed_or_cyclic(is_reversed_or_cyclic) {
            set_value(raw_value);
            // Always default to max
            section_length = internal_max();
        }
        /// Constructor to set a "false" for is_reversed_or_cyclic
        tree_item_t (tree_item_type_t type, size_t raw_value) 
            : tree_item_t(type, raw_value, false) {}
        /// Constructor to set max for the value
        tree_item_t (tree_item_type_t type, bool is_reversed_or_cyclic) 
            : tree_item_t(type, std::numeric_limits<size_t>::max(), is_reversed_or_cyclic) {}
        /// Constructor for just a type
        tree_item_t (tree_item_type_t type) 
            : tree_item_t(type, std::numeric_limits<size_t>::max(), false) {}
        // Setters
        inline void set_value(size_t new_value) {
            value = (new_value == std::numeric_limits<size_t>::max()) ? internal_max()
                                                                      : new_value;
        }
        inline void set_section_length(size_t new_length) {
            section_length = (new_length == std::numeric_limits<size_t>::max()) ? internal_max()
                                                                                : new_length;
        }
        // Getters
        inline tree_item_type_t get_type() const { return type; }
        inline size_t get_value() const { 
            return value == internal_max() ? std::numeric_limits<size_t>::max()
                                           : value;
        }
        // Different getters based on context for readability
        /// Is this seed reversed in the tree? Uses is_reversed_or_cyclic
        /// Only call on a SEED
        inline bool get_is_reversed() const { 
            if (this->type != ZipCodeTree::SEED) {
                throw std::runtime_error("Can't get reversedness of a tree item that isn't a seed");
            }
            return is_reversed_or_cyclic;
        }
        /// Is this bound part of a cyclic snarl? Uses is_reversed_or_cyclic
        /// Only call on a bound
        inline bool get_is_cyclic() const {
            if (this->type != ZipCodeTree::SNARL_START && this->type != ZipCodeTree::SNARL_END
                && this->type != ZipCodeTree::CHAIN_START && this->type != ZipCodeTree::CHAIN_END) {
                throw std::runtime_error("Can't get cyclicness of a tree item that isn't a bound");
            }
            return is_reversed_or_cyclic;
        }
        inline size_t get_section_length() const { 
            return section_length == internal_max() ? std::numeric_limits<size_t>::max()
                                                    : section_length;
        }
        /// What to add or subtract to this thing's index
        /// to get index of other bound
        inline int64_t get_other_bound_offset() const {
            if (section_length == internal_max()) {
                throw std::runtime_error("Can't get other bound offset of a tree item with no section length");
            }
            if (this->type == ZipCodeTree::SNARL_START
                || this->type == ZipCodeTree::CHAIN_START) {
                return section_length + 1;
            } else if (this->type == ZipCodeTree::SNARL_END
                       || this->type == ZipCodeTree::CHAIN_END) {
                return -(section_length + 1);
            } else {
                throw std::runtime_error("Can't get other bound offset of a tree item that isn't a bound");
            }
        }
        // We should never need to set reversedness after the fact,
        // so the setter is only provided with the is_cyclic name
        inline void set_is_cyclic(bool is_cyclic) { is_reversed_or_cyclic = is_cyclic; }
    };

    /**
     * Exposed type for a reference to an orientation of a seed
     */
    struct oriented_seed_t {
        /// The index of the seed in the seeds vector
        size_t seed;
        /// Is the seed traversed backwards in the tree?
        bool is_reversed;

        /// Raw value constructor
        oriented_seed_t(size_t seed_index, bool is_reversed)
            : seed(seed_index), is_reversed(is_reversed) {}
        
        /// Compare to other instances. TODO: Use default when we get C++20. 
        inline bool operator==(const oriented_seed_t& other) const {
            return seed == other.seed && is_reversed == other.is_reversed;
        }

        /// Compare to other instances. TODO: Use default when we get C++20. 
        inline bool operator!=(const oriented_seed_t& other) const {
            return !(*this == other);
        }
    };

    /// Get the number of items in the tree
    inline size_t get_tree_size() const {return zip_code_tree.size();}

    /// Access the value at [index] in the zip_code_tree
    inline tree_item_t get_item_at_index(size_t index) const {return zip_code_tree[index];};

    /// Get all the seeds in the tree, in left-to-right order
    /// Also returns their orientations
    /// Basically seed_itr but without all the extra baggage
    inline vector<oriented_seed_t> get_all_seeds() const {
        vector<oriented_seed_t> all_seeds;
        for (const auto& item : zip_code_tree) {
            if (item.get_type() == SEED) {
                all_seeds.emplace_back(item.get_value(), item.get_is_reversed());
            }
        }
        return all_seeds;
    }

    /// Helper for add_distance_matrix()
    /// Essentially re-calculates chain.distances.first/second for seeds
    /// which are inside nested snarls
    ///
    /// If the seed is in a nested snarl, this is the distance to snarl edge
    /// Otherwise it is 0
    /// Moves i along to find the seed, and returns the offset
    /// If right_to_left is true, then search leftward for the last seed
    /// Otherwise, search rightward for the first seed
    size_t get_offset_to_seed(size_t& i, bool right_to_left) const;

    /// Add snarl or chain end of matching type
    /// and sets up their section_length values
    void add_close_bound(size_t start_index);

protected:
    /// The actual tree structure
    vector<tree_item_t> zip_code_tree;

public:

    /**
     * Exposed type for a reference to an oriented seed at a distance
     */
    struct seed_result_t : public oriented_seed_t {
        /// The distance to the seed from some reference position
        size_t distance;

        /// Raw value constructor
        seed_result_t(size_t distance, size_t seed_index, bool is_reversed)
            : oriented_seed_t{seed_index, is_reversed}, distance(distance) {}

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
     * Iterator that visits seeds left to right in the tree's in-order traversal
     * 
     * ## Iteration flow
     * 
     * Each seed is visited from left to right, in the same order as
     * get_all_seeds() would return them.
     * 
     * The first time any given seed is visited, right_to_left is true,
     * indicating that a distance iterator should be set off going in the
     * default direction of right to left.
     * 
     * Once the right-to-left step has been done, then the iterator
     * decides whether to do a second, left-to-right step, with
     * right_to_left set to false. This will happen if inside a cyclic snarl,
     * even if not necessarily the direct parent. Calling ++() will 
     * toggle right_to_left instead of moving the internal iteration index.
     * The next call to ++() would then reset right_to_left back to true,
     * and finally advance the iterator rightward to the next seed.
     * 
     * For seeds not inside cyclic snarls, ++() will simply move rightward.
     * 
     * We need to traverse seeds in cyclic snarls in both directions because
     * the chains are not necessarily traversed in a single direction, and
     * are stored in an arbitrary orientation.
     * 
     * ## How to use this iterator
     * 
     * The intent of this iterator is to hide a lot of the fiddly bits required
     * to iterate over a zip tree with cyclic and non-cyclic snarls. All a
     * caller needs to do is
     * - call ++() to move to the next seed
     * - dereference the iterator to get the seeds at the current position
     * - create the distance iterator for this position by calling
     *   find_distances(seed_itr, distance_limit)
     * 
     * Seeds from both iterators will be automatically oriented based on the
     * direction of traversal. That is, you can Just Trust Me (TM) for all seed
     * orientations which iterators yield.
     */
    class seed_iterator {
    public:
        // References from ZipCodeTree to let us look up distance matrices
        /// A reference to the ziptree vector
        const vector<tree_item_t>& zip_code_tree;

        /// Make an iterator starting from start_index
        /// until the end of the given ziptree
        seed_iterator(size_t start_index, const ZipCodeTree& ziptree);
        
        // Iterators are copyable and movable
        seed_iterator(const seed_iterator& other) = default;
        seed_iterator(seed_iterator&& other) = default;
        seed_iterator& operator=(const seed_iterator& other) = default;
        seed_iterator& operator=(seed_iterator&& other) = default;

        /// Advance right until we hit another seed or the end
        /// Automatically handle updating chain_numbers and cyclic_snarl_nestedness
        seed_iterator& operator++();

        /// Compare for equality to see if we hit end
        /// This just trusts that the two iterators are for the same tree
        inline bool operator==(const seed_iterator& other) const { return index == other.index; }

        /// Compare for inequality
        inline bool operator!=(const seed_iterator& other) const { return !(*this == other); }

        /// What item does index point to?
        inline tree_item_t current_item() const { return zip_code_tree.at(index); }
        
        /// Get the index and orientation of the seed we are currently at.
        /// Also return all other seeds on the same position
        inline vector<oriented_seed_t> operator*() const { return current_seeds; }

        // Getters
        // We need these to make reverse iterators from forward ones.
        inline size_t get_index() const { return index; }
        inline bool get_right_to_left() const { return right_to_left; }
        inline stack<size_t> get_chain_numbers() const { return chain_numbers; }

    private:
        /// Within each parent snarl, which chain are we in?
        stack<size_t> chain_numbers;
        /// Where we are in the stored tree.
        size_t index;
        /// How far inside of a cyclic snarl are we?
        size_t cyclic_snarl_nestedness;
        /// Direction a distance iterator should be set off
        /// NOT the direction this iterator moves;
        /// a seed_iterator always moves left to right
        bool right_to_left;
        /// Seeds at this position
        vector<oriented_seed_t> current_seeds;
    };

    /// Get an iterator over indexes of seeds in the tree, left to right.
    inline seed_iterator begin() const { return seed_iterator(0, *this); }
    /// Get the end iterator for seeds in the tree, left to right.
    /// (Note that the last element will never be a seed)
    inline seed_iterator end() const { return seed_iterator(zip_code_tree.size() - 1, *this); }

    /**
     * Iterator that looks sideways in the tree from a seed,
     * possibly up to a maximum base distance.
     * 
     * Iteration flow:
     *
     * Original Python implementation: 
     * https://github.com/benedictpaten/long_read_giraffe_chainer_prototype/blob/b590c34055474b0c901a681a1aa99f1651abb6a4/zip_tree_iterator.py.
     * 
     * New cyclic snarl handling flowchart: 
     * https://docs.google.com/drawings/d/1diKMXCuMteR06fF64BhSttsD5RHh3eTnd8FtyYqs17Y/edit?usp=sharing
     * 
     * In order to de-duplicate work, this iterator will not actually look at
     * _everything_ before it; if it knows that seeds from a given chain would
     * have already output transitions to its own position. This occurs in
     * cyclic snarls. For example, since previous chains in a cyclic snarl will
     * have already visited all other chains, later such chains do not have to
     * visit these previous chains. Similarly, when the traversal enters a
     * cyclic snarl, it only has to handle the chains going in the opposite
     * orientation, since the ones in the same orientation will have simply
     * already seen the traversal's starting position by exiting the snarl.
     */
    class distance_iterator {
    public:
        /// Make a reverse iterator starting from start_index, until
        /// the given rend, with the given distance limit.
        distance_iterator(size_t start_index,
                          const vector<tree_item_t>& zip_code_tree,
                          std::stack<size_t> chain_numbers = std::stack<size_t>(), bool right_to_left = true,
                          size_t distance_limit = std::numeric_limits<size_t>::max());

        /// Move in right_to_left direction until we hit another seed or the end
        distance_iterator& operator++();

        /// Move index in right_to_left direction
        inline void move_index() {
            if (right_to_left) {
                --index;
            } else {
                ++index;
            }
        }

        /// Compare for equality to see if we hit end
        /// This just trusts that the two iterators are for the same tree
        inline bool operator==(const distance_iterator& other) const { return index == other.index; }

        /// Compare for inequality
        inline bool operator!=(const distance_iterator& other) const { return !(*this == other); }

        /// Is the iteration done?
        inline bool done() const { return index == end_index; }
        
        /// Get the index and orientation of the seed we are currently at, 
        /// and the distance to it.
        /// Automatically handles orientation depending on right_to_left
        seed_result_t operator*() const;

        /// Type for the state of the
        /// I-can't-believe-it's-not-a-pushdown-automaton
        enum State {
            S_START,
            S_SCAN_CHAIN,
            S_SCAN_DAG_SNARL,
            S_SCAN_CYCLIC_SNARL
        };

    private:
        /// Where we are in the stored tree.
        size_t index;
        /// Where we have to stop
        const size_t end_index;
        /// Where we started from
        const size_t original_index;
        /// Within each parent snarl, which chain are we in?
        std::stack<size_t> chain_numbers;
        /// Distance limit we will go up to
        size_t distance_limit;
        /// Whether we are looking right to left (true) or left to right (false)
        bool right_to_left;
        /// Original direction
        const bool original_right_to_left;
        /// References to the zip code tree to let us look up distance matrices
        const vector<tree_item_t>& zip_code_tree;
        /// Stack for computing distances.
        std::stack<size_t> stack_data;

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

        // Methods to look up distances to stack

        /// Push distances relevant to a given chain onto the stack
        /// chain_num is one-indexed, so the first chain is 1, and the last is N
        /// right_side indicates if distances exit from the right or left side
        void stack_snarl_distances(size_t snarl_start_i, size_t chain_num, bool right_side);

        /// Helper for stack_snarl_distances()
        /// Stack a single value from a triangular distance matrix
        void stack_matrix_value(size_t matrix_start_i, bool has_main_diagonal, size_t row, size_t col);

        /// Helper for stack_snarl_distances()
        /// Stack a single value below the running distance
        void stack_below_top(size_t value);

        // Helper functions for the automaton's state machine

        /// Current state of the automaton
        State current_state;

        /// Adopt a new state.
        inline void state(State new_state) { current_state = new_state; }

        /// Stop parsing because nothing else can be below the distance limit.
        /// This moves the current iterator it.
        void halt();

        /// Throw a domain_error that the current state/symbol combo is unimplemented.
        void unimplemented_error();

        /// What item does index point to?
        inline tree_item_t current_item() const { return zip_code_tree.at(index); }

        /// Check if the current symbol is an entrance/exit,
        /// based on the direction the iterator is going (right_to_left)
        inline bool entered_snarl() const {
            return (right_to_left && current_item().get_type() == ZipCodeTree::SNARL_END)
                    || (!right_to_left && current_item().get_type() == ZipCodeTree::SNARL_START);
        }
        inline bool exited_snarl() const {
            return (right_to_left && current_item().get_type() == ZipCodeTree::SNARL_START)
                    || (!right_to_left && current_item().get_type() == ZipCodeTree::SNARL_END);
        }
        inline bool entered_chain() const {
            return (right_to_left && current_item().get_type() == ZipCodeTree::CHAIN_END)
                    || (!right_to_left && current_item().get_type() == ZipCodeTree::CHAIN_START);
        }
        inline bool exited_chain() const {
            return (right_to_left && current_item().get_type() == ZipCodeTree::CHAIN_START)
                    || (!right_to_left && current_item().get_type() == ZipCodeTree::CHAIN_END);
        }

        /// Skip the current chain, jumping to the matching end
        /// and then one past it to continue the snarl
        /// Bound pair index should be right below current distance
        void skip_chain();

        /// Decide what to do right after entering a new chain.
        /// This chain's distance should be on top of the stack.
        void initialize_chain();

        /// Decide what to do right after entering a new snarl.
        /// Return whether the initialization was successful,
        /// i.e. whether the snarl was a non-root snarl.
        bool initialize_snarl(size_t chain_num);

        /// Decide what to do when re-entering a snarl,
        /// having already stacked up distances
        void continue_snarl();

        /// Tick the automaton, looking at the symbol at *it and updating the
        /// stack and current_state. Returns true to yield a value at the
        /// current symbol, or to halt, and false otherwise.
        bool tick();

    };

    /// Get a iterator starting from where a forward iterator is, up to a distance limit
    distance_iterator find_distances(const seed_iterator& from,
                                     size_t distance_limit = std::numeric_limits<size_t>::max()) const;

public:

    /*************** Debugging functions for validating zip trees ***********/

    /// Print the zip code tree to stderr
    /// ( and ) are used for the starts and ends of DAG snarls
    /// { and } are used for the starts and ends of cyclic snarls
    /// [ and ] are used for the starts and ends of chains
    /// seeds are printed as their positions
    void print_self(const vector<Seed>* seeds) const;

    /// Is the given node in a multicomponent chain, looping chain,
    /// or anything else that would cause it to not have exact distances?
    /// The distances are only guaranteed to be correct up to the distance limit
    /// Cyclic snarls don't count as being invalid
    bool node_is_invalid(nid_t id, const SnarlDistanceIndex& distance_index, 
                         size_t distance_limit = std::numeric_limits<size_t>::max()) const;

    /// Check that the tree is correct:
    /// 1. All snarl/chain boundaries are closed properly
    /// 2. The order of the items is logical
    /// 3. The distances between seeds (as output by iteration) are correct
    void validate_zip_tree(const SnarlDistanceIndex& distance_index, 
                           const vector<Seed>* seeds,
                           size_t distance_limit = std::numeric_limits<size_t>::max()) const;

    /// Helper function for validate_zip_tree() to check snarl/chain boundaries
    /// Ensures that all boundaries are matched in type,
    /// and that pair indexes are set up correctly
    /// Also checks that there is at least one seed in the tree
    /// Calls validate_snarl() for each snarl in the top-level chain
    void validate_boundaries(const SnarlDistanceIndex& distance_index, 
                             const vector<Seed>* seeds,
                             size_t distance_limit = std::numeric_limits<size_t>::max()) const;

    /// Helper function for validate_zip_tree() to check for a well-formed order
    /// 1. Do seeds have logical orientations relative to each other?
    /// 2. Do chains follow a [child, dist, child, dist, ... child] order?
    /// 3. Are there CHAIN_COUNTs right after each SNARL_START?
    void validate_zip_tree_order(const SnarlDistanceIndex& distance_index, 
                                 const vector<Seed>* seeds) const;

    /// Helper function for validate_zip_tree() to check distance iteration
    /// Uses the same iterator logic that the main chaining code does
    /// and for each pair of seeds output by the distance_iterator,
    /// compares their distance to the distance index
    void validate_seed_distances(const SnarlDistanceIndex& distance_index, 
                                 const vector<Seed>* seeds,
                                 size_t distance_limit = std::numeric_limits<size_t>::max()) const;

    /// Helper function for validate_zip_tree for just a snarl
    /// zip_iterator is an iterator to the snarl start
    /// At the end of the function, zip_iterator will be set to the snarl end
    void validate_snarl(std::vector<tree_item_t>::const_iterator& zip_iterator, 
                        const SnarlDistanceIndex& distance_index, 
                        const vector<Seed>* seeds,
                        size_t distance_limit = std::numeric_limits<size_t>::max()) const;

    /// Helper function for validate_snarl for a chain
    /// zip_iterator is an iterator to the chain start
    /// At the end of the function, zip_iterator will be set to the chain end
    void validate_chain(std::vector<tree_item_t>::const_iterator& zip_iterator, 
                        const SnarlDistanceIndex& distance_index, 
                        const vector<Seed>* seeds,
                        size_t distance_limit = std::numeric_limits<size_t>::max()) const;
    
    /// Count the number of snarls involved in the tree
    /// Returns a pair of <dag count, cyclic count>
    /// Assumes that the tree has already been filled in
    std::pair<size_t, size_t> dag_and_cyclic_snarl_count() const;

protected:

    /// Helper function to get orientation of a snarl tree node at a given depth
    /// does the same thing as the zipcode decoder's get_is_reversed_in_parent,
    /// except it also considers chains that are children of irregular snarls.
    ///
    /// We assume that all snarls are DAGs, so all children of snarls must only
    /// be traversable in one orientation through the snarl. This assumption
    /// doesn't work for cyclic snarls, but as their chains are traversed in
    /// both directions, their storage orientation doesn't matter.
    ///
    /// In a start-to-end traversal of a snarl, each node will only be
    /// traversable start-to-end or end-to-start. If traversable end-to-start,
    /// then it is considered to be oriented backwards in its parent
    static bool seed_is_reversed_at_depth (const Seed& seed, size_t depth, 
                                           const SnarlDistanceIndex& distance_index);

    friend class ZipCodeForest;
}; 

/*
    A collection of ZipCodeTrees
    The ZipCodeForest takes a set of seeds and makes ZipCodeTrees
    There will be a separate tree for each connected component
    or slice of a chain that is too far from anything else on both sides,
    using the given distance limit
*/
class ZipCodeForest {

    typedef SnarlDistanceIndexClusterer::Seed Seed;
    typedef ZipCodeTree::tree_item_type_t tree_item_type_t;
    typedef ZipCodeTree::tree_item_t tree_item_t;

    public:

    /// The actual data, a collection of ZipCodeTrees
    vector<ZipCodeTree> trees;

    /// Constructor
    ZipCodeForest() {};

    /// Populate the zip forest
    /// If a distance limit is given, then partition into subtrees that are
    /// farther than distance_limit from each other
    /// Otherwise, the forest will just be connected components
    void fill_in_forest(const vector<Seed>& seeds, 
                        const SnarlDistanceIndex& distance_index,
                        size_t distance_limit = std::numeric_limits<size_t>::max());

    private:


    /***************************************************************************

      Data structures and helper functions for construction

      **************************************************************************

      Construction is done in a depth-first traversal of the snarl tree.
      When each snarl tree node is visited, the start of the structure is added
      to the zip tree, then each of its children is added to the zip tree along
      with the distances between them, then the end of the structure is added.
      
      The traversal is accomplished by progressively sorting the seeds to
      identify the snarl tree structures that they lie on. Using the zip codes,
      the seeds can be sorted on each snarl tree structure separately. Seeds
      along a chain are sorted to be ordered along a chain, and seeds in a snarl
      are sorted by the child of the snarl that they are on. The seeds get
      sorted using a radix-like sort on each structure at each depth of the
      snarl tree, starting with the root and moving down. 

      "Intervals" of seeds in the sort order are used to keep track of the
      location on the snarl tree. An interval represents a range of seeds that
      are all on the same snarl tree structure. After sorting an interval at one
      depth, sub-intervals representing the children can be found. So first, the
      seeds are sorted into connected components and sliced into intervals 
      representing root-level snarls and chains. Each interval is then sorted to
      order the seeds along the snarl or chain, and new intervals are found
      representing ranges of seeds on the children.

      Sorting and tree-building are done at the same time, progressively at each
      structure in the snarl tree. The order of tree-building is based on a
      stack of intervals. The algorithm starts with an interval for each child
      of the root snarl. An interval is popped from the stack. Any incomplete
      snarls or chains that the interval is not a child of must be completed.
      Then, the snarl or chain that the interval represents is started in the
      zip tree, and any relevant distances are added. Intervals representing the
      children of the snarl or chain are found and added to the stack.
      
      This repeats until the stack is empty.

      Each snarl/chain in the zip code tree is comprised of start/end bounds,
      its children, and distances between children/bounds. So as each child is
      added, we will need to know what came before it in the parent snarl/chain
      so that we can add the distances. We also need to remember the ancestors
      of each snarl and chain as we are building them, so that we can close each
      structure properly. This information is stored in a forest_growing_state_t
      as the zip trees are being built.

     **************************************************************************/

    private:

    ////////////////////////////////////////////////////
    ///////////
    /////////// Data structures for building a zip tree
    //////////
    ////////////////////////////////////////////////////

    /// Structs which need to be declared for forest_growing_state_t
    /// See definitions later
    struct interval_state_t;
    struct sort_value_t;
    struct child_info_t;

    /// This stores information about the state of the forest as we fill it in
    struct forest_growing_state_t {
        /// Seeds which are to be put in the forest
        const vector<Seed>* seeds;

        /// Distance index for the graph being represented
        const SnarlDistanceIndex* distance_index;

        /// Sort order for the seeds
        vector<size_t> seed_sort_order;


        /// This stores the sort value and code type of each seed 
        /// This will change as forest building progresses
        /// but it will be set for the relevant seed immediately before sorting
        /// The values also get used to calculate distance, 
        /// as long as they have been set for the correct depth
        vector<sort_value_t> sort_values_by_seed;

        /// Stores the previous things of the current structure at each depth
        /// The children are stored at the depth of their parents.
        /// For example, for a root chain, the vector at index 0 would have the
        /// chain start, seeds on the chain, and snarl starts on the chain.
        /// Similarly, for a top-level snarl, at depth 1, the second 
        /// vector would contain the starts of chains at depth 2 
        vector<vector<child_info_t>> sibling_indices_at_depth;

        /// We build a forest of trees. This is an index into trees
        /// to indicate which is actively being worked on.
        ///
        /// A new tree is formed either when a new top-level chain is found
        /// (or a slice of a top-level chain if far enough from previous thing),
        /// or when part of a chain in a snarl is too far from everything else.
        /// In the second case, the entire subtree is found before determining
        /// that it should be a subtree, and is copied into a new ZipCodeTree.
        /// So only one tree is actively being added to at a time.
        /// 
        /// Note that this can't be an actual pointer to the forest because
        /// the address may move if the vectors get shifted around in memory.
        size_t active_tree_index;

        /// If part of a chain is unreachable with the rest of the chain,
        /// then we want to split it off into a separate zipcode tree.
        /// This tracks all open chains as an index to the start of the chain
        /// in the current active tree, and a boolean for if the chain start is
        /// farther than distance_limit from anything else in the snarl tree.
        ///
        /// If the index is for a CHAIN_START, it includes the whole chain.
        /// If it points to a SEED/SNARL_START, then it is a slice.
        ///
        /// Any time something gets added to a chain or the chain is closed,
        /// check if the distance to anything following is >distance_limit.
        /// If it is, copy from the start of the chain/slice into a new tree.
        vector<pair<size_t, bool>> open_chains;

        /// A stack of intervals representing snarl tree nodes.
        /// These are yet to be sorted and added to the zip tree.
        /// After an interval is popped, child intervals are added in
        /// The stack structure ensures a depth-first processing order
        forward_list<interval_state_t> intervals_to_process;
    
        /// Intervals that are currently open.
        /// These represent ancestors of whatever is currently being worked on.
        /// So the size is the depth of the snarl tree
        vector<interval_state_t> open_intervals;


        /// The overall distance limit for splitting of new connected components
        size_t distance_limit;

        /// Constructor given seeds and a distance index
        forest_growing_state_t(const vector<Seed>& seeds, const SnarlDistanceIndex& distance_index, 
                               size_t distance_limit) :
            seeds(&seeds), distance_index(&distance_index), 
            distance_limit(distance_limit), active_tree_index(std::numeric_limits<size_t>::max()) {

            // This represents the current sort order of the seeds
            seed_sort_order.assign(seeds.size(), 0);
            for (size_t i = 0 ; i < seed_sort_order.size() ; i++) {
                seed_sort_order[i] = i;
            }
            sort_values_by_seed.resize(seeds.size());
        }

    };

    /// Store edgemost seeds in chains when creating a snarl's distance matrix
    /// We make one pass over each chain, remembering its edge seeds, and then
    /// use those seeds to calculate all distances
    struct seed_info_t {
        /// The seed at the edge of the chain (first/last)
        pos_t seed_pos;
        /// Length of the seed's node
        size_t node_length;
        /// The seed's non-reversed zipcode
        const ZipCode& zipcode;
        /// Its offset from the edge (i.e. chain.distances.first/second)
        size_t flank_offset;
        /// Whether to use the right/left side of the seed
        bool right_side;
        /// Rank of the seed's chain in the snarl
        size_t rank;
        /// If this seed is in a nested snarl, and we're calculating
        /// using minimum_distance, then we might need to subtract the
        /// offset from inner seed to inner snarl edge

        size_t nested_snarl_offset;

        /// Pass the seed's index (which is looked up from forest_state.seeds),
        /// whether its position should be reversed, and then a few raw values
        /// In addition, snarl depth & forest_state are used to look up seed info
        seed_info_t(size_t index, bool is_rev, size_t flank, bool right_side,
                    size_t nested_snarl_offset,
                    const size_t& depth, const forest_growing_state_t& forest_state);

        inline pos_t reverse_seed() const {
            return reverse(seed_pos, node_length);
        }
    };

    /// For children of snarls, we need to remember the siblings and start bound
    /// that came before them so we can record their distances
    /// This holds the indices (into zip_code_tree) of each seed or chain start,
    /// and each start and child chain start of a snarl
    /// For the children of a chain, the value is the prefix sum in the chain
    /// (relative to orientation of the top-level chain, not the chain itself)
    /// For the children of a snarl, the value is the index of the CHAIN_START
    /// First seed in the chain must be found by looping through zip_code_tree
    struct child_info_t {


        /// A value associated with the item, 
        /// offset in a chain or index of the snarl child start
        size_t value;
    
        /// For children of snarls, distance to the left and right of the chain
        /// that gets added to edges in the snarl
        /// The first item is a map of {seed index : distance},
        /// in case the first seed gets snipped off
        /// The second item is the distance to the rightmost seed in the chain
        std::pair<std::unordered_map<size_t, size_t>, size_t> distances;

        /// If the item is a child of a chain, its chain component
        size_t chain_component : 26; 

        /// Current item type
        ZipCodeTree::tree_item_type_t type : 5;


        /// Is the sibling reversed; only used for children of snarls,
        /// to indicate that the child is traversed backwards 
        bool is_reversed = false;

        /// Constructor for type/value leaving all other fields as defaults
        child_info_t(ZipCodeTree::tree_item_type_t type, size_t value) : type(type), value(value) {}
    };

    /// Used for sorting. Represents one interval along zipcode_sort_order,
    /// which corresponds to a snarl tree node at the given depth
    struct interval_state_t {

        /// Indices into zipcode_sort_order
        size_t interval_start : 26; // inclusive
        size_t interval_end : 26;   // exclusive

        /// Is the snarl tree node reversed relative to the top-level chain?
        bool is_reversed : 1;

        /// The type of the snarl tree structure
        /// For nodes on chains, all seeds on the chain not nested in snarls
        /// are put in the same interval, regardless of if on the same node
        ZipCode::code_type_t code_type : 5;

        /// Interval depth in the snarl tree
        size_t depth : 14;

        /// Constructor for raw values
        interval_state_t (size_t start, size_t end, size_t rev, ZipCode::code_type_t type, size_t depth) :
            interval_start(start), interval_end(end), is_reversed(rev), code_type(type), depth(depth) {
        }
    };

    struct sort_value_t {
        private:
        /// Value to sort on
        size_t sort_value;
        /// The type of the snarl tree structure.
        ZipCode::code_type_t code_type : 5;

        // For chains, used to indicate the order of the child of a chain
        // Since the offset represents the space between nucleotides,
        // two positions on different nodes could have the same offset.
        // Similarly, a snarl could have the same prefix sum as a node.
        // For example, in this graph:
        //                2
        //               [AA]
        //           1  /   \  3
        //          [AA] --- [AA]
        // The positions n1-0 and 3+0, and the snarl 1-3 all have an offset of 2
        // To solve this, the prefix sum of a chain will be multiplied by 3,
        // 1 will be added to snarls, and 2 will be added to the node
        // with an 0 offset in the node (node 3 if chain is traversed forward)

        /// The +0/+1/+2 offset used for chain prefix sums
        size_t chain_order : 3; 

        /// If the item is a child of a chain, its chain component
        size_t chain_component : 24;

        public:
        /// Constructor to make an empty sort value
        sort_value_t() : sort_value(std::numeric_limits<size_t>::max()),
                         code_type(ZipCode::EMPTY),
                         chain_order(7),
                         chain_component(0) {};
        /// Constructor for raw values
        sort_value_t (size_t sort_value, ZipCode::code_type_t code_type, size_t chain_order) :
            sort_value(sort_value), code_type(code_type), chain_order(chain_order), chain_component(0) {};

        /// Get the value used for sorting
        inline size_t get_sort_value() const {
            // The sort value for chains is actually prefix sum*3+chain_order, 
            // to account for different nodes having the same prefix sum
            return chain_order != 7
                       ? (sort_value * 3) + chain_order
                       : sort_value;
        };

        /// Get the value used for distance finding
        inline size_t get_distance_value() const { return sort_value; }

        // Getters
        inline ZipCode::code_type_t get_code_type() const { return code_type; }
        inline size_t get_chain_component() const { return chain_component; }

        // Setters
        inline void set_sort_value(size_t value) { sort_value = value; }
        inline void set_code_type(ZipCode::code_type_t type) { code_type = type; }
        inline void set_chain_order(size_t order) { chain_order = order; }
        inline void set_chain_component(size_t component) { chain_component = component; }

    };

    ////////////////////////////////////////////////////////////////////////////
    //  Functions for sorting/finding intervals of seeds along the snarl tree
    ////////////////////////////////////////////////////////////////////////////


    /// Sorts an interval, which must contain seeds on the same snarl/chain/node
    /// Sorting is linear-ish along top-level chains, topological-ish in snarls.
    /// Uses radix_sort_zipcodes() and default_sort_zipcodes()
    ///
    /// For chains, everything is sorted with the prefix sum value of the chain
    /// itself from the distance index, not the order in the zip code tree.
    /// Everything will be sorted in the order of the zip code tree, 
    /// but the values will be set from the distance index. Later, the values
    /// may be out of order or may need to be subtracted from the length of
    /// the chain to get the distance to the ends of the chain
    void sort_one_interval(forest_growing_state_t& forest_state, 
                           const interval_state_t& interval) const;

    /// Helper function for sort_one_interval() to sort seeds using radix sort
    /// Sorts the slice of seeds in the given interval of zipcode_sort_order,
    /// a vector of indices into seeds. Reverse order if reverse_order.
    /// The interval also has an is_reversed field, which refers to the
    /// orientation in the snarl tree. min_ and max_value are the minimum and
    /// maximum value being sorted on. If sort_by_chain_component is true, then
    //// sort on the chain component in sort_values
    ///
    /// This should run in linear time, but is dependent on the values
    /// being sorted on to have a small range
    void radix_sort_zipcodes(vector<size_t>& zipcode_sort_order, 
                             const vector<sort_value_t>& sort_values_by_seed,
                             const interval_state_t& interval, bool reverse_order,
                             size_t min_value, size_t max_value, bool sort_by_chain_component = false) const; 

    /// Helper function for sort_one_interval() to sort seeds using std::sort
    /// Sorts the slice of seeds in the given interval of zipcode_sort_order, 
    /// a vector of indices into seeds
    void default_sort_zipcodes(vector<size_t>& zipcode_sort_order, 
                               const vector<sort_value_t>& sort_values_by_seed,
                               const interval_state_t& interval, bool reverse_order) const; 

    /// Assuming that the range of seeds in forest_state.sort_values_by_seed
    /// given by interval is sorted, prepend child intervals to next_intervals.
    /// The new intervals get added in their sort order, so the start of a chain
    /// will be at the start of the list, to be popped first.
    ///
    /// For children of chains, seeds on the chain itself and not nested will
    /// be on the same interval if there are no seeds in snarls between them,
    /// even if they are not on the same node
    void get_next_intervals(forest_growing_state_t& forest_state, 
                            const interval_state_t& interval,
                            std::forward_list<interval_state_t>& next_intervals) const;

    //////////////////////////////////////////////////////
    ///////////          functions for building the trees
    /////////////////////////////////////////////////////

    /// Move a slice of a chain into a new tree
    /// Chain copied from forest_state.open_chains.back().first to end
    /// This is used when a part is too far from the rest to be in the same tree
    /// Returns whether a whole chain was moved (true) or just a slice (false)
    bool move_slice(forest_growing_state_t& forest_state, const size_t& depth);

    /// Open a chain that starts at the current_seed
    /// Also record its presence and distance-to-start in the parent snarl
    /// Assumes that the chain *has* a parent snarl, i.e. isn't top-level
    void open_chain(forest_growing_state_t& forest_state, const interval_state_t& interval);

    /// Close a chain that ends at last_seed
    /// If the chain was empty, remove it and anything relating to it
    /// If must be spliced out, take out a subtree
    /// Otherwise, add the end of the chain and, if the chain was in a snarl,
    /// remember the distance to the end of the chain
    void close_chain(forest_growing_state_t& forest_state, const size_t& depth, 
                     const Seed& last_seed, bool chain_is_reversed);

    /// Add a seed (or snarl starting at a seed) and the preceeding chain edge
    /// If the seed is far enough from the previous thing in the chain and
    /// it can be a new slice, split off a subtree
    /// depth is chain child's depth (may also be chain depth if trivial)
    /// seed_index is the index of the current seed in the list of seeds
    void add_child_to_chain(forest_growing_state_t& forest_state,
                            const size_t& depth, const size_t& seed_index, 
                            bool child_is_reversed, bool chain_is_reversed);

    /// Start a new snarl of the given type at the given depth
    void open_snarl(forest_growing_state_t& forest_state, const size_t& depth, bool is_cyclic_snarl);

    /// Close a snarl at the given depth with the given last_seed
    /// If the snarl has no children, then delete the whole thing
    /// Otherwise, add all necessary distances and close it
    void close_snarl(forest_growing_state_t& forest_state, const size_t& depth, 
                     const Seed& last_seed, bool last_is_reversed);

    /// Add a triangular distance matrix for a snarl to its front
    /// The matrix starts with a CHAIN_COUNT with the number of child chains,
    /// and then is a list of EDGEs; for each item in order, all distances to it
    /// from all previous items, possibly including self-loops
    /// Returns the size of the distance matrix added (plus the CHAIN_COUNT)
    size_t add_distance_matrix(forest_growing_state_t& forest_state, 
                               const size_t& depth, bool snarl_is_reversed);

    /// Helper for add_distance_matrix()
    /// Look up seeds for chain edges and remember seed_info_t for each
    /// Returns [c1_L, c1_R, c2_L, c2_R, ...] for each chain side in the snarl
    vector<seed_info_t> get_edge_seeds(const forest_growing_state_t& forest_state, 
                                       const size_t& depth) const;

    /// Helper for add_distance_matrix() to add the last row
    /// These are edges from everything prior to the snarl end bound
    void add_edges_to_end(vector<tree_item_t>& dist_matrix,
                          forest_growing_state_t& forest_state, 
                          const size_t& depth, const vector<seed_info_t>& edge_seeds,
                          bool snarl_is_reversed, bool is_cyclic_snarl) const;
    
    /// Helper for add_distance_matrix() to add the chains' rows
    /// These are edges from everything prior to the given chain end bound
    void add_edges_for_chains(vector<tree_item_t>& dist_matrix,
                              forest_growing_state_t& forest_state, 
                              const size_t& depth, const vector<seed_info_t>& edge_seeds,
                              bool snarl_is_reversed, bool is_cyclic_snarl) const;

    /************ Helper functions for debugging ************/

    public:

    /// Print each zip code tree in the forest to stderr
    /// Ziptrees are prefaced by their index, e.g. "0: <tree.print_self()>"
    inline void print_self(const vector<Seed>* seeds) const {
        for (size_t i = 0 ; i < trees.size() ; i++) {
            const auto& tree = trees[i];
            cerr << i << ": ";
            tree.print_self(seeds);
        }
    }

    /// Check that the forest is correct:
    /// 1. Each tree is valid in itself
    /// 2. All seed positions (i.e. ignoring duplicates) appear at least once
    void validate_zip_forest(const SnarlDistanceIndex& distance_index, 
                             const vector<Seed>* seeds,
                             size_t distance_limit=std::numeric_limits<size_t>::max()) const;
};

/// Print an item type to a stream
std::ostream& operator<<(std::ostream& out, const ZipCodeTree::tree_item_type_t& type);
/// Print an iterator state to a stream
std::ostream& operator<<(std::ostream& out, const ZipCodeTree::distance_iterator::State& state);

}

namespace std {

/// Make an item type into a string
std::string to_string(const vg::ZipCodeTree::tree_item_type_t& type);
/// Make an iterator state into a string
std::string to_string(const vg::ZipCodeTree::distance_iterator::State& state);

/// Hash functor to hash oriented_seed_t with std::hash
template <> struct hash<vg::ZipCodeTree::oriented_seed_t>
{
    /// Produce a hash of an oriented_seed_t.
    size_t operator()(const vg::ZipCodeTree::oriented_seed_t& item) const
    {
        // Hash it just as we would a pair.
        return hash<pair<size_t, bool>>()(make_pair(item.seed, item.is_reversed));
    }
};

/// Hash functor to hash seed_result_t with std::hash
template <> struct hash<vg::ZipCodeTree::seed_result_t>
{
    /// Produce a hash of a seed_result_t.
    size_t operator()(const vg::ZipCodeTree::seed_result_t& item) const
    {
        // Hash it just as we would a tuple.
        return hash<tuple<size_t, bool, size_t>>()(make_tuple(item.seed, item.is_reversed, item.distance));
    }
};

/// Explain to the STL algorithms what kind of iterator the zip code tree
/// forward iterator is.
template<>
struct iterator_traits<vg::ZipCodeTree::seed_iterator>{
    using value_type = vg::ZipCodeTree::oriented_seed_t;   
    using iterator_category = forward_iterator_tag;
};

/// Explain to the STL algorithms what kind of iterator the zip code tree
/// reverse iterator is.
template<>
struct iterator_traits<vg::ZipCodeTree::distance_iterator>{
    using value_type = vg::ZipCodeTree::seed_result_t;   
    using iterator_category = forward_iterator_tag;
};


}

#endif
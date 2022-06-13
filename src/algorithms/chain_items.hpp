#ifndef VG_ALGORITHMS_CHAIN_ITEMS_HPP_INCLUDED
#define VG_ALGORITHMS_CHAIN_ITEMS_HPP_INCLUDED

/**
 * \file
 * Algorithms for chaining subalignments into larger alignments.
 *
 * To use these algorithms, decide on the type (Item) you want to chain up.
 *
 * Then, make a ChainingSpace<Item>, or a ChainingSpace<Item, Source> if your
 * Items need to be interpreted in the context of some source object (like a
 * seed hit needs to be interpreted in the context of its source minimizer).
 *
 * Then, make a dynamic programming table: vector<int> or vector<traced_score_t>.
 *
 * Then, call chain_items_dp() to fill in the dynamic programming table and get
 * the score of the best chain.
 *
 * If you used traced_score_t, you can use chain_items_traceback() to get a
 * traceback of the chain's items in order.
 *
 * Helper entry points are find_best_chain() and score_best_chain() which set
 * up the DP for you and do the traceback if appropriate.
 */

#include "../gbwt_extender.hpp"
#include "../snarl_seed_clusterer.hpp"
#include "../seed_clusterer.hpp"
#include "../handle.hpp"

namespace vg {
namespace algorithms {

using namespace std;

#define debug_chaining

/**
 * We want to be able to chain a reordered subset of things without moving the originals, so we use this view over things stored in a vector.
 *
 * Both the backing collection and the indexes must outlive the view.
 */
template<typename Item>
struct VectorView {
    const vector<Item>& items;
    const vector<size_t>* indexes;
    
    /**
     * Make a VectorView of a whole vector. Provides an implicit conversion.
     */
    VectorView(const vector<Item>& items) : items(items), indexes(nullptr) {
        // Nothing to do!
    }
    
    /**
     * Make a VectorView of a reordered subset of a vector.
     */
    VectorView(const vector<Item>& items, const vector<size_t>& indexes) : items(items), indexes(&indexes) {
        // Nothing to do!
    }
    
    /**
     * Get an item by index.
     */
    const Item& operator[](size_t index) const {
        if (indexes) {
            return items[(*indexes)[index]];
        } else {
            return items[index];
        }
    }
    
    /**
     * Get the total number of items.
     */
    size_t size() const {
        if (indexes) {
            return indexes->size();
        } else {
            return items.size();
        }
    }
    
    /**
     * Determine if there are no items.
     */
    bool empty() const {
        if (indexes) {
            return indexes->empty();
        } else {
            return items.empty();
        }
    }
    
    /**
     * Call the given callback with a dense and properly ordered vector of the items.
     */
    void with_vector(const std::function<void(const vector<Item>&)>& callback) const {
        if (indexes) {
            // We need to reorder
            vector<Item> reordered;
            reordered.reserve(indexes->size());
            for (auto& i : *indexes) {
                reordered.emplace_back(items[i]);
            }
            callback(reordered);
        } else {
            // We already have this
            callback(items);
        }
    }
};

/// We support chaining different kinds of things, so we have a type that
/// abstracts out accessing their chaining-relevant fields and measuring
/// distance between them.
/// We support things that come from other source things (i.e. minimizers).
/// See BaseChainingSpace for the actual interface.
template<typename Item, typename Source = void>
struct ChainingSpace {
    // Nothing here; this exists to be specialized.
};

/// Base class of specialized ChainingSpace classes. Defines the ChainingSpace interface.
template<typename Item>
struct BaseChainingSpace {
    
    const Aligner& scoring;
    
    const SnarlDistanceIndex* distance_index;
    const HandleGraph* graph;
    
    BaseChainingSpace(const Aligner& scoring,
                      const SnarlDistanceIndex* distance_index = nullptr,
                      const HandleGraph* graph = nullptr) :
        scoring(scoring),
        distance_index(distance_index),
        graph(graph) {
        
        // Nothing to do!
    }
    
    /// Get the score collected by visiting the item. Should be an alignment score.
    virtual int score(const Item& item) const {
        // Default implementation assumes a perfect match
        return scoring.match * read_length(item);
    }
    
    /// Get where the item starts in the read
    virtual size_t read_start(const Item& item) const = 0;
    
    /// Get where the item past-ends in the read
    virtual size_t read_end(const Item& item) const = 0;
    
    /// Get the length of the item in the read
    virtual size_t read_length(const Item& item) const {
        return read_end(item) - read_start(item);
    }
    
    /// Find the start position of the item in the graph. Graph must be set.
    virtual pos_t graph_start(const Item& item) const = 0;
    
    /// Find the past-end position of the item in the graph. Graph must be set.
    virtual pos_t graph_end(const Item& item) const = 0;
    
    /// Get the length of the item in the graph. Graph must be set.
    virtual size_t graph_length(const Item& item) const {
        return read_length(item);
    }
    
    /// Get the number of handles in the path through the graph. Must always be 1 or more. Graph must be set.
    virtual size_t graph_path_size(const Item& item) const = 0;
    
    /// Get the handle at the given index in the path through the graph. Graph must be set.
    virtual handle_t graph_path_at(const Item& item, size_t index) const = 0;
    
    /// Get the offset on the first handle at which the item starts. Graph must be set.
    virtual size_t graph_path_offset(const Item& item) const = 0;
    
    /// Turn an item into a WFAAlignment. Graph must be set.
    virtual WFAAlignment to_wfa_alignment(const Item& item) const {
        // Default implementation assumes no mismatches.
        return {
            vector<handle_t>(graph_path_begin(item), graph_path_end(item)),
            {{WFAAlignment::match, (uint32_t)read_length(item)}},
            (uint32_t)graph_path_offset(item),
            (uint32_t)read_start(item),
            (uint32_t)read_length(item),
            score(item)
        };
    }
    
    /**
     * We use these iterators for traversing graph paths.
     */
    struct PathIterator {
        using value_type = handle_t;
        using difference_type = ptrdiff_t;
        using reference = const value_type&;
        using pointer = const value_type*;
        using iterator_category = std::input_iterator_tag;
    
        size_t index;
        int direction;
        const Item* item;
        const BaseChainingSpace<Item>* space;
        
        // We aren't necessarily iterating over something that actually stores
        // handles, so we need a place to point to for ->
        handle_t at;
        
        const handle_t& operator*() const {
            return at;
        }
        
        const handle_t* operator->() const {
            return &at;
        }
        
        PathIterator operator++(int) {
            PathIterator clone = *this;
            ++*this;
            return clone;
        }
        
        PathIterator& operator++() {
            if (index == 0 && direction == -1) {
                index = numeric_limits<size_t>::max();
            } else {
                index += direction;
                if (index < space->graph_path_size(*item)) {
                    at = space->graph_path_at(*item, index);
                }
            }
            return *this;
        }
        
        const bool operator==(const PathIterator& other) const {
            // Don't bother comparing items and spaces and ats
            return index == other.index && direction == other.direction;
        }
        
        const bool operator!=(const PathIterator& other) const {
            return !(*this == other);
        }
    };
    
    /// Get an iterator to the start of the path of handles taken by the item. Graph must be set.
    virtual PathIterator graph_path_begin(const Item& item) const {
        return {0, 1, &item, this, graph_path_at(item, 0)};
    }
    
    /// Get an iterator to the end of the path of handles taken by the item. Graph must be set.
    virtual PathIterator graph_path_end(const Item& item) const {
        return {graph_path_size(item), 1, &item, this, {}};
    }
    /// Get a reverse iterator to the back of the path of handles taken by the item. Graph must be set.
    virtual PathIterator graph_path_rbegin(const Item& item) const {
        return {graph_path_size(item) - 1, -1, &item, this, graph_path_at(item, graph_path_size(item) - 1)};
    }
    
    /// Get a reverse iterator to the past-front of the path of handles taken by the item. Graph must be set.
    virtual PathIterator graph_path_rend(const Item& item) const {
        return {numeric_limits<size_t>::max(), -1, &item, this, {}};
    }
    
    /**
     * Get the sequence in the read that the given item uses.
     */
    virtual string get_read_sequence(const Item& item, const std::string& full_read_sequence) const {
        return full_read_sequence.substr(read_start(item), read_length(item));
    }
    
    /**
     * Get the sequence in the graph that the given item uses.
     * Graph must be set.
     */
    virtual string get_graph_sequence(const Item& item) const {
        stringstream ss;
        for (size_t i = 0; i < graph_path_size(item); i++) {
            handle_t here = graph_path_at(item, i);
            ss << graph->get_sequence(here).substr((i == 0) ? graph_path_offset(item) : 0);
        }
        return ss.str();
    }
    
    /// Return true if the first item is a perfect chain that can't be beat,
    /// and false otherwise.
    virtual bool has_perfect_chain(const VectorView<Item>& items) const = 0;
    
    /**
     * Return the amount by which the end of the left item is past the position
     * before the start of the right item, in the graph. Returns 0 if they do
     * not actually overlap.
     * 
     * Doesn't actually match the whole graph paths against each other; if
     * there's a handle where the left one ends and then the right one starts,
     * there's no overlap.
     * 
     * Can return a number larger than the length of one item if it is
     * contained in the other.
     *
     * Graph must be set.
     */
    virtual size_t get_graph_overlap(const Item& left,
                                     const Item& right) const {
                                     
        if (!graph) {
            return 0;
        }
                                              
        // We need to worry about containment, or the left thing actually
        // starting after the right thing.
        
        // We also don't really care if they take different overall paths in
        // the graph. We find overlap if a start node is shared.
        
        // This is the handle the right extension starts on
        handle_t right_starts_on = graph_path_at(right, 0);
        
        // Scan back along the left path until we find the handle the right path starts on 
        auto last_right_start_in_left = graph_path_rbegin(left);
        while (last_right_start_in_left != graph_path_rend(left) && *last_right_start_in_left != right_starts_on) {
            ++last_right_start_in_left;
        }
        
        if (last_right_start_in_left != graph_path_rend(left)) {
            // The right extension starts at some handle where the left extension also is.
            
            // Work out how many bases of the left extension are before that handle.
            // We start with all the bases on the path.
            size_t left_bases_before_shared_handle = 0;
            auto it = last_right_start_in_left;
            ++it;
            if (it != graph_path_rend(left)) {
                // There are some left bases before the shared handle.
                while (it != graph_path_rend(left)) {
                    left_bases_before_shared_handle += graph->get_length(*it);
                    ++it;
                }
                // Since we have some bases, we back out the offset.
                left_bases_before_shared_handle -= graph_path_offset(left);
                
                // Work out how many bases of the left extension are on or after that handle
                size_t left_bases_remaining = graph_length(left) - left_bases_before_shared_handle;
                
                // Work out how many bases into the handle the right extension starts
                size_t right_offset = graph_path_offset(right);
                
                // Compute how many bases of the left extension the right extension
                // starts at or before and return that.
                if (left_bases_remaining <= right_offset) {
                    // They don't actually overlap
                    return 0;
                } else {
                    // They overlap by however many bases the left extension goes into
                    // this handle past where the right extension starts.
                    return left_bases_remaining - right_offset;
                }
                    
            } else {
                // The left extension actually starts on the shared handle where the right extension also starts. It might be completely before or completely after the right one actually.
                size_t left_past_end = graph_path_offset(left) + graph_length(left);
                if (left_past_end <= graph_path_offset(right)) {
                    // Actually there's no overlap, left will end before right starts.
                    return 0;
                } else {
                    // There's some overlap
                    return left_past_end - graph_path_offset(right);
                }
            }
        } else {
            // The handle the right extension starts on doesn't actually exist in the left path.
            // So the right extension either doesn't actually overlap, or else contains the left one.
            // Or else the right extension starts somewhere else and reads into the left extension's path.
            
            // This is the handle the left extension starts on
            handle_t left_starts_on = graph_path_at(left, 0);
            
            // So, find the earliest place the left extension's path could start in the right one
            auto first_left_start_in_right = graph_path_begin(right);
            // And how many bases come before it along the left extension
            size_t right_bases_before_shared_handle = 0;
            while (first_left_start_in_right != graph_path_end(right) && *first_left_start_in_right != left_starts_on) {
                right_bases_before_shared_handle += graph->get_length(*first_left_start_in_right);
                ++first_left_start_in_right;
            }
            // Since we know the right path's start handle does not appear on the
            // left path, we know we went through at least one handle and we can
            // safely back out the offset.
            right_bases_before_shared_handle -= graph_path_offset(right);
            
            if (first_left_start_in_right != graph_path_end(right)) {
                // The left path actually does have its starting handle somewhere in the right path.
                
                // We know how many right path bases are before that handle.
                // Then we add to that overlap the distance along the handle to the
                // start of the left extension, and also the length of the left
                // extension.
                return right_bases_before_shared_handle + graph_path_offset(left) + graph_length(left);
            } else {
                // Neither path has the start handle of the other.
                // We could maybe come up with some notion of overlap between points along the paths, but we can't do it for the extensions as a whole.
                // TODO: is this going to make the overall algorithm work well when these e.g. start/end on SNP nodes?
                return 0;
            }
        }
    
    }
                                    
    /**
     * Return the amount by which the end of the left item is past the
     * position before the start of the right item, in the read. Returns 0 if
     * they do not actually overlap.
     * 
     * Can return a number larger than the length of one item if it is
     * contained in the other.
     */
    virtual size_t get_read_overlap(const Item& left,
                                    const Item& right) const {
        size_t l = read_end(left);
        size_t r = read_start(right);
        if (l < r) {
            return 0;
        }
        return l - r;
    }

    /**
     * Get the minimum graph distance between the end of the left item and the
     * position before the start of the right item. Returns
     * std::numeric_limits<size_t>::max() if there is no route to take.
     *
     * Graph must be set.
     */
    virtual size_t get_graph_distance(const Item& left,
                                      const Item& right) const {
    
        if (!distance_index || !graph) {
            return numeric_limits<size_t>::max();
        }
        
        // Find where the left extension past-ends
        pos_t left_past_end = graph_end(left);
        // Find where the right one starts
        pos_t right_start = graph_start(right);
        // Get the oriented minimum distance from the index
        size_t distance = distance_index->minimum_distance(
            id(left_past_end), is_rev(left_past_end), offset(left_past_end),
            id(right_start), is_rev(right_start), offset(right_start),
            false, graph);
        // And return it
        return distance;
    }
    
    /**
     * Get the distance in the read between two items, or 0 if they abut or
     * overlap.
     */
    virtual size_t get_read_distance(const Item& left,
                                     const Item& right) const {
        size_t l = read_end(left);
        size_t r = read_start(right);
        if (r < l) {
            return 0;
        }
        return r - l;
    }
    
    /**
     * Get the highest score in points that the alignment of the region between
     * two Items could possibly get. Should assume that any bases that might
     * match actually match.
     *
     * If the items overlap in the read or the graph, returns
     * std::numeric_limits<int>::min() because one side of the alignment would
     * have negative length.
     */
    virtual int transition_score_upper_bound(const Item& left,
                                             const Item& right) const {
                                             
        if (read_end(left) == read_start(right)) {
            // They abut in the read.
            
            // How far apart are they in the graph?
            size_t graph_distance = get_graph_distance(left, right);
            if (graph_distance == 0) {
                // They also abut in the graph. No alignment to do.
#ifdef debug_chaining
                cerr << "Abutting items score 0" << endl;
#endif
                return 0;
            } else if (graph_distance == numeric_limits<size_t>::max()) {
                // They are unreachable in the graph (maybe an overlap?)
#ifdef debug_chaining
                cerr << "Items unreachable in graph" << endl;
#endif
                return numeric_limits<int>::min();
            } else {
                // There is a gap; pure deletion to 0 bp. Charge for that.
#ifdef debug_chaining
                cerr << "Items separated by " << graph_distance << " bp deletion" << endl;
#endif
                return -(scoring.gap_open + (graph_distance - 1) * scoring.gap_extension);
            }
        } else if (read_end(left) < read_start(right)) {
            // There's a gap in the read
            size_t read_distance = get_read_distance(left, right);
            
            // See how far they have to go in the graph
            size_t graph_distance = get_graph_distance(left, right);
            
            if (graph_distance == numeric_limits<size_t>::max()) {
                // These aren't actually reachable in the graph.
#ifdef debug_chaining
                cerr << "Items unreachable in graph" << endl;
#endif
                return numeric_limits<int>::min();
            } else {
                // Otherwise, see if there's any length change
                size_t length_change = (read_distance > graph_distance) ? (read_distance - graph_distance) : (graph_distance - read_distance);
                // And see how many pases could be matches
                size_t possible_matches = std::min(read_distance, graph_distance);
                
#ifdef debug_chaining
                cerr << "Items separated by " << possible_matches << " bp possible matches and at least a " << length_change << " bp indel" << endl;
#endif
                
                // The number of possible matches gives us a base score
                int jump_points = scoring.match * possible_matches;
                if (length_change > 0) {
                    // We have to charge for a gap though
                    jump_points -= (scoring.gap_open + (length_change - 1) * scoring.gap_extension);
                }
                return jump_points;
            }
        } else {
            // If there's an overlap in the read. say we can't do it.
            
#ifdef debug_chaining
            cerr << "Items overlap in read" << endl;
#endif
            
            return numeric_limits<int>::min();
        }
    }
};


/// This is how you chain up a bunch of GaplessExtension items
template<>
struct ChainingSpace<GaplessExtension, void>: public BaseChainingSpace<GaplessExtension> {
    using Item = GaplessExtension;
    
    // Keep the constructor
    using BaseChainingSpace<Item>::BaseChainingSpace;
    
    int score(const Item& item) const {
        return item.score;
    }
    
    size_t read_start(const Item& item) const {
        return item.read_interval.first;
    }
    
    size_t read_end(const Item& item) const {
        return item.read_interval.second;
    }
    
    pos_t graph_start(const Item& item) const {
        return make_pos_t(item.starting_position(*graph));
    }
    
    pos_t graph_end(const Item& item) const {
        return make_pos_t(item.tail_position(*graph));
    }
    
    size_t graph_path_size(const Item& item) const {
        return item.path.size();
    }
    
    handle_t graph_path_at(const Item& item, size_t index) const {
        return item.path[index];
    }
    
    size_t graph_path_offset(const Item& item) const {
        return item.offset;
    }
    
    virtual WFAAlignment to_wfa_alignment(const Item& item) const {
        // Handle mismatches in the GaplessExtension
        return WFAAlignment::from_extension(item);
    }
    
    bool has_perfect_chain(const VectorView<Item>& items) const {
        bool result;
        items.with_vector([&](const vector<Item>& dense_items) {
            // Once we're sure we have a dense vector we can call the
            // GaplessExtender function that answers this question, which
            // doesn't work with views.
            result = GaplessExtender::full_length_extensions(dense_items);
        });
        return result;
    }
    
};

/// This is how you chain up new seeds
template<typename Source>
struct ChainingSpace<NewSnarlSeedClusterer::Seed, Source> : public BaseChainingSpace<NewSnarlSeedClusterer::Seed> {
    using Item = NewSnarlSeedClusterer::Seed;
    
    // These seeds can't really be interpreted without their sources, which
    // they reference by index. Sources need to support a forward_offset() and a length.
    const vector<Source>& sources;
    
    ChainingSpace(const vector<Source>& sources,
                  const Aligner& scoring,
                  const SnarlDistanceIndex* distance_index,
                  const HandleGraph* graph) :
        BaseChainingSpace<Item>(scoring, distance_index, graph), sources(sources) {
        
        // Nothing to do!
    }
    
    // The seeds really only know their start positions; a single seed from a
    // minimizer occurrence in the graph might actually have multiple
    // corresponding end positions if the graph forks but retains the same
    // sequence. We need to consistently show a certain endpoint in the read
    // and the graph, so we just take the rest of the node that the seed starts
    // on.
    
    size_t read_start(const Item& item) const {
        return sources[item.source].forward_offset();
    }
    
    size_t read_end(const Item& item) const {
        return sources[item.source].forward_offset() + read_length(item);
    }
    
    size_t read_length(const Item& item) const {
        if (graph) {
            return graph_length(item);
        } else {
            // No graph so just treat the whole thing as the item.
            return sources[item.source].length;
        }
    }
    
    size_t graph_length(const Item& item) const {
        handle_t start_handle = graph_path_at(item, 0);
        // Length in graph is min of either remaining node or minimizer length
        return std::min((size_t) sources[item.source].length,
                        graph->get_length(start_handle) - graph_path_offset(item));
    }
    
    pos_t graph_start(const Item& item) const {
        return item.pos;
    }
    
    pos_t graph_end(const Item& item) const {
        pos_t start = graph_start(item);
        return make_pos_t(id(start), is_rev(start), offset(start) + graph_length(item));
    }
    
    size_t graph_path_size(const Item& item) const {
        return 1;
    }
    
    handle_t graph_path_at(const Item& item, size_t index) const {
        pos_t start = graph_start(item);
        return graph->get_handle(id(start), is_rev(start));
    }
    
    size_t graph_path_offset(const Item& item) const {
        return offset(graph_start(item));
    }
    
    bool has_perfect_chain(const VectorView<Item>& items) const {
        return false;
    }
};

/// This is how you chain up old seeds
template<typename Source>
struct ChainingSpace<SnarlSeedClusterer::Seed, Source> : public BaseChainingSpace<SnarlSeedClusterer::Seed> {
    using Item = SnarlSeedClusterer::Seed;
    
    // These seeds can't really be interpreted without their sources, which
    // they reference by index. Sources need to support a forward_offset() and a length.
    const vector<Source>& sources;
    
    ChainingSpace(const vector<Source>& sources,
                  const Aligner& scoring,
                  const HandleGraph* graph) :
        BaseChainingSpace<Item>(scoring, nullptr, graph), sources(sources) {
        
        // Nothing to do!
    }
    
    // The seeds really only know their start positions; a single seed from a
    // minimizer occurrence in the graph might actually have multiple
    // corresponding end positions if the graph forks but retains the same
    // sequence. We need to consistently show a certain endpoint in the read
    // and the graph, so we just take the rest of the node that the seed starts
    // on.
    
    size_t read_start(const Item& item) const {
        return sources[item.source].forward_offset();
    }
    
    size_t read_end(const Item& item) const {
        return sources[item.source].forward_offset() + read_length(item);
    }
    
    size_t read_length(const Item& item) const {
        if (graph) {
            return graph_length(item);
        } else {
            // No graph so just treat the whole thing as the item.
            return sources[item.source].length;
        }
    }
    
    size_t graph_length(const Item& item) const {
        handle_t start_handle = graph_path_at(item, 0);
        // Length in graph is min of either remaining node or minimizer length
        return std::min((size_t) sources[item.source].length,
                        graph->get_length(start_handle) - graph_path_offset(item));
    }
    
    pos_t graph_start(const Item& item) const {
        return item.pos;
    }
    
    pos_t graph_end(const Item& item) const {
        pos_t start = graph_start(item);
        return make_pos_t(id(start), is_rev(start), offset(start) + graph_length(item));
    }
    
    size_t graph_path_size(const Item& item) const {
        return 1;
    }
    
    handle_t graph_path_at(const Item& item, size_t index) const {
        pos_t start = graph_start(item);
        return graph->get_handle(id(start), is_rev(start));
    }
    
    size_t graph_path_offset(const Item& item) const {
        return offset(graph_start(item));
    }
    
    bool has_perfect_chain(const VectorView<Item>& items) const {
        return false;
    }
};

// For doing scores with backtracing, we use this type, which is a
// score and a number for the place we came from to get it.
// It can be maxed with std::max() as long as the destinations are all filled in right.
using traced_score_t = pair<int, size_t>;

/// Accessors for attributes of a type when used as a (possibly provenance-traced) DP table score.
template<typename Score>
struct score_traits {
    // We don't actually define the interface here, because templates don't do any sort of inheritance.
};

/// This is how to use a traced_score_t as a DP score, and all the operations
/// you can do with it.
template<>
struct score_traits<traced_score_t> {
    using Score = traced_score_t;
    
    /// What is the sentinel for an empty provenance?
    /// Use a function instead of a constant because that's easier when we're just a header.
    inline static size_t nowhere() {
        return numeric_limits<size_t>::max();
    }
    
    /// What's the default value for an empty table cell?
    /// Use a function instead of a constant because that's easier when we're just a header.
    inline static Score unset() {
        return {0, nowhere()};
    }
    
    /// Pack a score together with its provenance
    static Score annotate(int points, size_t from) {
        return {points, from};
    }

    /// Accessor to get the score
    static int& score(Score& s) {
        return s.first;
    }
    
    /// Accessor to get the score, when const
    static const int& score(const Score& s) {
        return s.first;
    }
    
    /// Accessor to get the source
    static size_t& source(Score& s) {
        return s.second;
    }
    
    /// Accessor to get the source, when const
    static const size_t& source(const Score& s) {
        return s.second;
    }
    
    /// Max in a score from a DP table. If it wins, record provenance.
    static void max_in(Score& dest, const vector<Score>& options, size_t option_number) {
        auto& option = options[option_number];
        if (score(option) > score(dest) || source(dest) == nowhere()) {
            // This is the new winner.
            score(dest) = score(option);
            source(dest) = option_number;
        }
    }
    
    /// Get a score from a table and record provenance in it.
    static Score score_from(const vector<Score>& options, size_t option_number) {
        traced_score_t got = options[option_number];
        source(got) = option_number;
        return got;
    }
    
    /// Add (or remove) points along a route to somewhere. Return a modified copy.
    static Score add_points(const Score& s, int adjustment) {
        return annotate(score(s) + adjustment, source(s));
    }
};

/// This is how to use an int as a DP score, and all the operations you can do
/// with it.
template<>
struct score_traits<int> {
    using Score = int;
    
    /// What is the sentinel for an empty provenance?
    /// Use a function instead of a constant because that's easier when we're just a header.
    inline static size_t nowhere() {
        return numeric_limits<size_t>::max();
    }
    
    /// What's the default value for an empty table cell?
    /// Use a function instead of a constant because that's easier when we're just a header.
    inline static Score unset() {
        return 0;
    }

    /// Pack a score together with its provenance
    static Score annotate(int points, size_t from) {
        return points;
    }

    /// Accessor to get the score.
    static Score& score(int& s) {
        return s;
    }
    
    /// Accessor to get the score, when const.
    static const Score& score(const int& s) {
        return s;
    }
    
    /// Accessor to get the source, when const (always nowhere)
    static size_t source(const int& s) {
        return nowhere();
    }
    // Note that source can't be set.
    
    /// Max in a score from a DP table. If it wins, record provenance.
    static void max_in(Score& dest, const vector<Score>& options, size_t option_number) {
        dest = std::max(dest, options[option_number]);
    }
    
    /// Get a score from a table.
    static Score score_from(const vector<Score>& options, size_t option_number) {
        return options[option_number];
    }
    
    /// Add (or remove) points along a route to somewhere. Return a modified copy.
    static Score add_points(const Score& s, int adjustment) {
        return s + adjustment;
    }
};

/// Print operator
ostream& operator<<(ostream& out, const traced_score_t& value);

/**
 * Fill in the given DP table for the best chain score ending with each
 * item. Returns the best observed score overall from that table,
 * with provenance to its location in the table, if tracked in the type.
 * Assumes some items exist.
 */
template<typename Score, typename Item, typename Source = void, typename Collection = VectorView<Item>>
Score chain_items_dp(vector<Score>& best_chain_score,
                     const VectorView<Item>& to_chain,
                     const ChainingSpace<Item, Source>& space);
                     
/**
 * Fill in the given DP table for the best chain score ending with each
 * item. Returns the best observed score overall from that table,
 * with provenance to its location in the table, if tracked in the type.
 * Assumes some items exist.
 *
 * Uses a finite lookback in items when checking where we can come from to
 * reach an item.
 */
template<typename Score, typename Item, typename Source = void, typename Collection = VectorView<Item>>
Score chain_items_finite_lookback_dp(vector<Score>& best_chain_score,
                                     const Collection& to_chain,
                                     const ChainingSpace<Item, Source>& space,
                                     size_t lookback = 5);

/**
 * Trace back through in the given DP table from the best chain score.
 */
template<typename Score, typename Item, typename Source = void, typename Collection = VectorView<Item>>
vector<size_t> chain_items_traceback(const vector<Score>& best_chain_score,
                                     const Collection& to_chain,
                                     const Score& best_past_ending_score_ever,
                                     const ChainingSpace<Item, Source>& space);

/**
 * Chain up the given group of items. Determines the best score and
 * traceback that can be obtained by chaining items together, using the
 * given gap open and gap extend penalties to charge for either overlaps or
 * gaps in coverage of the read.
 *
 * Overlaps are charged only gap open/extend penalties; multiple matches to
 * the same read base are scored as matches.
 *
 * Overlaps may result in one item containing another.
 *
 * Input items must be sorted by start position in the read.
 *
 * Optionally takes a distance index and a graph, and uses distances in the
 * graph alogn with distances in the read to score transitions between
 * items.
 *
 * Returns the score and the list of indexes of items visited to achieve
 * that score, in order.
 */
template<typename Item, typename Source = void, typename Collection = VectorView<Item>>
pair<int, vector<size_t>> find_best_chain(const Collection& to_chain,
                                          const algorithms::ChainingSpace<Item, Source>& space);

/**
 * Score the given group of items. Determines the best score that can be
 * obtained by chaining items together, using the given space to define gap
 * open and gap extend penalties to charge for either overlaps or gaps in
 * coverage of the read.
 *
 * Overlaps are charged only gap open/extend penalties; multiple matches to the
 * same read base are scored as matches.
 *
 * Overlaps may result in one item containing another.
 *
 * Input items must be sorted by start position.
 */
template<typename Item, typename Source = void, typename Collection = VectorView<Item>>
int score_best_chain(const Collection& to_chain,
                     const ChainingSpace<Item, Source>& space);

// --------------------------------------------------------------------------------

template<typename Score, typename Item, typename Source, typename Collection>
Score chain_items_dp(vector<Score>& best_chain_score,
                     const Collection& to_chain,
                     const ChainingSpace<Item, Source>& space) {
                    
    // Grab the traits into a short name so we can use the accessors concisely.
    using ST = score_traits<Score>;

#ifdef debug_chaining
    cerr << "Chaining group of " << to_chain.size() << " items" << endl;
#endif
    
    
    // This is a collection of one or more non-full-length extended seeds.
    
    // We use a sweep line algorithm to find relevant points along the read: item starts or ends.
    // This records the last base to be covered by the current sweep line.
    int64_t sweep_line = 0;
    // This records the first base not covered by the last sweep line.
    int64_t next_unswept = 0;
    
    // And we track the next unentered item
    size_t unentered = 0;
    // And the max observed end in read space so far
    size_t max_end = space.read_end(to_chain[unentered]); 
    
    // Extensions we are in are in this min-heap of past-end position and item number.
    using ending_at_t = pair<size_t, size_t>;
    priority_queue<ending_at_t, vector<ending_at_t>, std::greater<ending_at_t>> end_heap;
    
    // We track the best score for a chain reaching the position before this one and ending in a gap.
    // We never let it go below 0.
    // Will be 0 when there's no gap that can be open
    Score best_gap_score = ST::unset();
    
    // We track the score for the best chain ending with each item
    best_chain_score.clear();
    best_chain_score.resize(to_chain.size(), ST::unset());
    
    // And we're after the best score overall that we can reach when an item ends
    Score best_past_ending_score_ever = ST::unset();
    
    // Overlaps are more complicated.
    // We need a heap of all the items for which we have seen the
    // start and that we can thus overlap.
    // We filter things at the top of the heap if their past-end positions
    // have occurred.
    // So we store pairs of score we get backtracking to the current
    // position, and past-end position for the thing we are backtracking
    // from.
    // We can just use the standard max-heap comparator.
    priority_queue<pair<Score, size_t>> overlap_heap;
    
    // We encode the score relative to a counter that we increase by the
    // gap extend every base we go through, so we don't need to update and
    // re-sort the heap.
    int overlap_score_offset = 0;
    
    while(next_unswept <= max_end) {
        // We are processed through the position before next_unswept.
        
        // Find a place for sweep_line to go
        
        // Find the next seed start
        int64_t next_seed_start = numeric_limits<int64_t>::max();
        if (unentered < to_chain.size()) {
            next_seed_start = space.read_start(to_chain[unentered]);
        }
        
        // Find the next seed end
        int64_t next_seed_end = numeric_limits<int64_t>::max();
        if (!end_heap.empty()) {
            next_seed_end = end_heap.top().first;
        }
        
        // Whichever is closer between those points and the end, do that.
        sweep_line = min(min(next_seed_end, next_seed_start), (int64_t) max_end);
        
        // So now we're only interested in things that happen at sweep_line.
#ifdef debug_chaining
        cerr << "Sweep to " << sweep_line << endl;
#endif
        
        // Compute the distance from the previous sweep line position
        // Make sure to account for next_unswept's semantics as the next unswept base.
        int sweep_distance = sweep_line - next_unswept + 1;
        
        // We need to track the score of the best thing that past-ended here
        Score best_past_ending_score_here = ST::unset();
        
        while(!end_heap.empty() && end_heap.top().first == sweep_line) {
            // Find anything that past-ends here
            size_t past_ending = end_heap.top().second;
            
            // Mix it into the score
            ST::max_in(best_past_ending_score_here, best_chain_score, past_ending);
            
            // Remove it from the end-tracking heap
            end_heap.pop();
        }
#ifdef debug_chaining
        cerr << "Best score of an item past-ending here: " << best_past_ending_score_here << endl;
#endif
        

        // Pick between that and the best score overall
        best_past_ending_score_ever = std::max(best_past_ending_score_ever, best_past_ending_score_here);
        
        if (sweep_line == max_end) {
            // We don't need to think about gaps or backtracking anymore since everything has ended
            break;
        }
        
        // Update the overlap score offset by removing some gap extends from it.
        overlap_score_offset += sweep_distance * space.scoring.gap_extension;
        
        // The best way to backtrack to here is whatever is on top of the heap, if anything, that doesn't past-end here.
        Score best_overlap_score = ST::unset();
        while (!overlap_heap.empty()) {
            // While there is stuff on the heap
            if (overlap_heap.top().second <= sweep_line) {
                // We are already past this thing, so drop it
                overlap_heap.pop();
            } else {
                // This is at the top of the heap and we aren't past it
                // Decode and use its score offset if we only backtrack to here.
                best_overlap_score = ST::add_points(overlap_heap.top().first, overlap_score_offset);
                // Stop looking in the heap
                break;
            }
        }
#ifdef debug_chaining
        cerr << "Best score of overlapping back to here: " << best_overlap_score << endl;
#endif
        
        // The best way to end 1 before here in a gap is either:
        
        if (best_gap_score != ST::unset()) {
            // Best way to end 1 before our last sweep line position with a gap, plus distance times gap extend penalty
            ST::score(best_gap_score) -= sweep_distance * space.scoring.gap_extension;
        }
        
        // Best way to end 1 before here with an actual item, plus the gap open part of the gap open penalty.
        // (Will never be taken over an actual adjacency)
        best_gap_score = std::max(ST::unset(), std::max(best_gap_score, ST::add_points(best_past_ending_score_here, -(space.scoring.gap_open - space.scoring.gap_extension))));
#ifdef debug_chaining
        cerr << "Best score here but in a gap: " << best_gap_score << endl;
#endif
        
        while (unentered < to_chain.size() && space.read_start(to_chain[unentered]) == sweep_line) {
            // For each thing that starts here
            
            // Bump out the max end
            max_end = std::max(max_end, space.read_end(to_chain[unentered]));
            
            // We want to compute how much we should charge to come from
            // each place we could come from, and then come from the best
            // one.
            Score overlap_option = best_overlap_score;
            Score gap_option = best_gap_score;
            Score here_option = best_past_ending_score_here;
            
            if (space.distance_index && space.graph) {
                // We have a distance index, so we might charge different
                // amounts for different transitions depending on the graph
                // distance between the items.
                
                // So we adjust the options we have bgased on where we are
                // coming from and where we are going to.
                
                // We could take:
                // An overlap from source(overlap_option):
                if (ST::source(overlap_option) != ST::nowhere()) {
                    // See whether and how much they overlap in the graph
                    size_t graph_overlap = space.get_graph_overlap(to_chain[ST::source(overlap_option)], to_chain[unentered]);
                    // And also get the overlap in the read.
                    size_t read_overlap = space.get_read_overlap(to_chain[ST::source(overlap_option)], to_chain[unentered]);
                    if (graph_overlap > read_overlap) {
                        // They overlap in the graph more than they do in the read.
                        // We know read overlap is positive or we wouldn't be considering taking an overlap in the read.
                        assert(read_overlap > 0);
                        // Back out the extra gap extensions
                        overlap_option = ST::add_points(overlap_option, space.scoring.gap_extension * (graph_overlap - read_overlap));
                        // TODO: We need to not double-count the matches right?
                        // When we go to synthesize an actual alignment we can cut out the double-matched area and match it only one place.
#ifdef debug_chaining
                        cerr << "Graph overlap exceeds read overlap (" << graph_overlap << " > " << read_overlap << "), charge only difference: " << overlap_option << endl;
#endif
                    } else if (graph_overlap == read_overlap) {
                        // They overlap in the graph exactly as much as they do in the read.
                        // Back out the whole read gap
                        overlap_option = ST::add_points(overlap_option, space.scoring.gap_open + space.scoring.gap_extension * (read_overlap - 1));
                        // TODO: We need to not double-count the matches right?
#ifdef debug_chaining
                        cerr << "Graph overlap equals read overlap (" << graph_overlap << " = " << read_overlap << "), do not charge for gap: " << overlap_option << endl;
#endif
                    } else if (graph_overlap > 0) {
                        // They overlap in the graph but less than they do in the read.
                        // Back out some of the read gap
                        overlap_option = ST::add_points(overlap_option, space.scoring.gap_extension * (read_overlap - graph_overlap));
#ifdef debug_chaining
                        cerr << "Graph overlap is smaller than read overlap (" << graph_overlap << " < " << read_overlap << "), charge only difference: " << overlap_option << endl;
#endif
                    } else {
                        // They don't overlap in the graph at all; they might abut or be separated or be unreachable.
                        
                        // Compute the graph distance if they don't overlap in the graph
                        size_t graph_distance = space.get_graph_distance(to_chain[ST::source(overlap_option)], to_chain[unentered]);
                        
                        if (graph_distance == numeric_limits<size_t>::max()) {
                            // We can't actually get between the relevant graph areas to do the overlap.
                            // TODO: check if we can get between the points on either side of the overlapped area instead?
                            // Right now just say we can't make this connection.
                            overlap_option = ST::unset();
#ifdef debug_chaining
                            cerr << "Overlap option is actually unreachable in the graph: " << overlap_option << endl;
#endif
                        } else {
                            // We need to charge for this extra graph distance being deleted, in addition to what we charge for the existing overlap gap.
                            overlap_option = ST::add_points(overlap_option, -space.scoring.gap_extension * graph_distance);
#ifdef debug_chaining
                            cerr << "Graph overlap is actually a distance of " << graph_distance << " vs. read overlap of " << read_overlap << ", charge for additional distance: " << overlap_option << endl;
#endif
                        }
                    }
                }
                
                // A gap in the read from ST::source(gap_option):
                if (ST::source(gap_option) != ST::nowhere()) {
                    // Compute the graph distance
                    size_t graph_distance = space.get_graph_distance(to_chain[ST::source(gap_option)],
                                                                     to_chain[unentered]);
                    if (graph_distance == numeric_limits<size_t>::max()) {
                        // This is actually unreachable in the graph.
                        // So say we can't actually do this.
                        gap_option = ST::unset();
                    } else {
                        // Compare to the read distance in the same order
                        size_t read_distance = space.get_read_distance(to_chain[ST::source(gap_option)], to_chain[unentered]);
                        if (graph_distance == read_distance) {
                            // This is actually even length.
                            // Charge nothing for the gap by backing out the penalty.
                            gap_option = ST::add_points(gap_option, space.scoring.gap_open + space.scoring.gap_extension * (graph_distance - 1));
#ifdef debug_chaining
                            cerr << "Graph gap is equal to read gap (" << graph_distance << " = " << read_distance << "), do not charge for gap: " << gap_option << endl;
#endif
                        } else if (graph_distance > read_distance && read_distance > 0) {
                            // Graph distance is longer, but read distance is positive.
                            // Back out what we charged for the extra extensions.
                            gap_option = ST::add_points(gap_option, space.scoring.gap_extension * (graph_distance - read_distance));
#ifdef debug_chaining
                            cerr << "Graph gap exceeds read gap (" << graph_distance << " > " << read_distance << "), charge only difference: " << gap_option << endl;
#endif
                        } else if (read_distance == 0) {
                            // This isn't really a gap-in-read option, we're
                            // just looking at the gap started from the thing
                            // we come right after.
                            gap_option = ST::unset();
#ifdef debug_chaining
                            cerr << "Read gap just started and would be length 0, so prohibit it: " << gap_option << endl;
#endif
                        } else {
                            // Read distance is longer.
                            // Back out extend penalties from the part of
                            // the gap that the graph also has, leaving
                            // only the remaining surplus read gap.
                            gap_option = ST::add_points(gap_option, space.scoring.gap_extension * graph_distance);
#ifdef debug_chaining
                            cerr << "Graph gap is smaller than read gap (" << graph_distance << " < " << read_distance << "), charge only difference: " << gap_option << endl;
#endif
                        }
                    }
                }
                
                // An adjacency in the read from ST::source(here_option):
                if (ST::source(here_option) != ST::nowhere()) {
                    // Compute the graph distance
                    size_t graph_distance = space.get_graph_distance(to_chain[ST::source(here_option)],
                                                                     to_chain[unentered]);
                    if (graph_distance == numeric_limits<size_t>::max()) {
                        // This is actually unreachable in the graph.
                        // So say we can't actually do this.
                        // TODO: might be better to do this as an overlap in the graph if we can?
                        here_option = ST::unset();
#ifdef debug_chaining
                        cerr << "Adjacency in read is unreachable in graph: " << here_option << endl;
#endif
                    } else if (graph_distance > 0) { 
                        // There's a gap in the graph but not in the read.
                        // Charge for the gap
                        here_option = ST::add_points(here_option, -space.scoring.gap_open - space.scoring.gap_extension * (graph_distance - 1));
#ifdef debug_chaining
                        cerr << "Adjacency in read is a gap of " << graph_distance << " in graph, charge for gap: " << here_option << endl;
#endif
                    } else {
                        // These abut in the read and the graph, so do nothing.
#ifdef debug_chaining
                        cerr << "Adjacency in read is adjacency in graph: " << here_option << endl;
#endif
                    }
                }
                
            }
            
            // Compute its chain score, allowing us to just start here also.
            best_chain_score[unentered] = ST::add_points(
                std::max(std::max(ST::unset(),
                                  overlap_option),
                         std::max(gap_option,
                                  here_option)),
                space.score(to_chain[unentered]));
#ifdef debug_chaining
            cerr << "Best score of chain ending in item " << unentered << ": " << best_chain_score[unentered] << endl;
#endif
            
            // Compute its backtrack-to-here score and add it to the backtracking heap
            // We want how far we would have had to have backtracked to be
            // able to preceed the base we are at now, where this thing
            // starts. That's a gap open for the last base, and then an
            // extend for each base before it.
            size_t item_length = space.read_length(to_chain[unentered]);
            // We assume all items are 1 or more bases.
            Score raw_overlap_score = ST::add_points(
                ST::score_from(best_chain_score, unentered),
                -space.scoring.gap_open - space.scoring.gap_extension * (item_length - 1));
            Score encoded_overlap_score = ST::add_points(raw_overlap_score, -overlap_score_offset);
            
            // Stick it in the overlap heap
            overlap_heap.emplace(encoded_overlap_score, space.read_end(to_chain[unentered]));
#ifdef debug_chaining
            cerr << "Overlap encoded as " << encoded_overlap_score << " available until " << space.read_end(to_chain[unentered]) << endl;
#endif
            
            // Add it to the end finding heap
            end_heap.emplace(space.read_end(to_chain[unentered]), unentered);
            
            // Advance and check the next thing to start
            unentered++;
            
            if (unentered < to_chain.size()) {
                // Bump out the max end
                max_end = std::max(max_end, space.read_end(to_chain[unentered]));
            }
        }
        
        // Move next_unswept to sweep_line.
        // We need to add 1 since next_unswept is the next *un*included base
        next_unswept = sweep_line + 1;
    }
    
    return best_past_ending_score_ever; 
}

template<typename Score, typename Item, typename Source, typename Collection>
Score chain_items_finite_lookback_dp(vector<Score>& best_chain_score,
                                     const Collection& to_chain,
                                     const ChainingSpace<Item, Source>& space,
                                     size_t lookback) {
                    
    // Grab the traits into a short name so we can use the accessors concisely.
    using ST = score_traits<Score>;

#ifdef debug_chaining
    cerr << "Chaining group of " << to_chain.size() << " items with lookback of " << lookback << endl;
#endif
    
    // Make our DP table big enough
    best_chain_score.resize(to_chain.size(), ST::unset());
    
    // What's the winner so far?
    Score best_score = ST::unset();
   
    for (size_t i = 0; i < to_chain.size(); i++) {
        // For each item
        auto& here = to_chain[i];
        
        for (size_t back = 1; back < std::min(lookback, i) + 1; back++) {
            // For each possible source
            auto& source = to_chain[i - back];
            
#ifdef debug_chaining
            cerr << "Consider " << space.read_start(source) << "-" << space.read_end(source) << " to " << space.read_start(here) << "-" << space.read_end(here) << endl;
#endif
            
            // How much does it pay (+) or cost (-) to make the jump from there
            // to here?
            int jump_points = space.transition_score_upper_bound(source, here);
            
            if (jump_points != numeric_limits<int>::min()) {
                
                // Get the score we are coming from
                typename ST::Score source_score = ST::score_from(best_chain_score, i - back);
                
                // And the score with the transition
                typename ST::Score from_source_score = ST::add_points(source_score, jump_points);
                
                // Remember that we could make this jump
                best_chain_score[i] = std::max(best_chain_score[i],
                                               from_source_score);
                                               
#ifdef debug_chaining
                cerr << "We can reach " << i << " with " << from_source_score << " which is " << source_score << " and a transition of " << jump_points << endl;
#endif
            }
        }
        
#ifdef debug_chaining
        cerr << "Best way to reach " << i << " is " << best_chain_score[i] << endl;
#endif
        
        // See if this is the best overall
        ST::max_in(best_score, best_chain_score, i);
        
#ifdef debug_chaining
        cerr << "Best chain end so far: " << best_score << endl;
#endif
        
    }
    
    return best_score;
}

template<typename Score, typename Item, typename Source, typename Collection>
vector<size_t> chain_items_traceback(const vector<Score>& best_chain_score,
                               const Collection& to_chain,
                               const Score& best_past_ending_score_ever,
                               const ChainingSpace<Item, Source>& space) {
    
    // Now we need to trace back.
    vector<size_t> traceback;
    size_t here = score_traits<Score>::source(best_past_ending_score_ever);
    if (here != score_traits<Score>::nowhere()) {
#ifdef debug_chaining
        cerr << "Chain ends at #" << here << " at " << space.read_start(to_chain[here])
            << "-" << space.read_end(to_chain[here])
            << " with score " << best_past_ending_score_ever << endl;
#endif
        while(here != score_traits<Score>::nowhere()) {
            traceback.push_back(here);
#ifdef debug_chaining
            cerr << "Which gets score " << best_chain_score[here] << endl;
#endif
            here = score_traits<Score>::source(best_chain_score[here]);
#ifdef debug_chaining
            if (here != score_traits<Score>::nowhere()) {
                cerr << "And comes after #" << here
                << " at " << space.read_start(to_chain[here])
                << "-" << space.read_end(to_chain[here]) << endl;
            } else {
                cerr << "And is first" << endl;
            }
#endif
        }
        // Flip it around front-ways
        std::reverse(traceback.begin(), traceback.end());
    }
    
#ifdef debug_chaining
    cerr << "Best score of chain overall: " << best_past_ending_score_ever << endl;
#endif

    return traceback;
}

template<typename Item, typename Source, typename Collection>
pair<int, vector<size_t>> find_best_chain(const Collection& to_chain,
                                          const ChainingSpace<Item, Source>& space) {
                                                                 
    if (to_chain.empty()) {
        return std::make_pair(0, vector<size_t>());
    } else if (space.has_perfect_chain(to_chain)) {
        // These are full-length matches. We already have the score.
        return std::make_pair(space.score(to_chain[0]), vector<size_t>{0});
    } else {
        
        // We actually need to do DP
        vector<traced_score_t> best_chain_score;
        traced_score_t best_past_ending_score_ever = chain_items_finite_lookback_dp(best_chain_score,
                                                                                    to_chain,
                                                                                    space);
        // Then do the traceback and pair it up with the score.
        return std::make_pair(
            score_traits<traced_score_t>::score(best_past_ending_score_ever),
            chain_items_traceback(best_chain_score, to_chain, best_past_ending_score_ever, space));
    }
}

template<typename Item, typename Source, typename Collection>
int score_best_chain(const Collection& to_chain,
                     const ChainingSpace<Item, Source>& space) {
    
    if (to_chain.empty()) {
        return 0;
    } else if (space.has_perfect_chain(to_chain)) {
        // These are full-length matches. We already have the score.
        return space.score(to_chain[0]);
    } else {
        // Do the DP but without the traceback.
        vector<int> best_chain_score;
        return algorithms::chain_items_finite_lookback_dp(best_chain_score,
                                                          to_chain,
                                                          space);
    }
}

}
}

#endif

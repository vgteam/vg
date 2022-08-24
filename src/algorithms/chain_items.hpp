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
 *
 * There is also reseed_fallow_regions() which can orchestrate finding more
 * items to chain when there are big gaps between the existing ones.
 */
 
#include "extract_connecting_graph.hpp"

#include "../gbwt_extender.hpp"
#include "../snarl_seed_clusterer.hpp"
#include "../seed_clusterer.hpp"
#include "../handle.hpp"
#include "../explainer.hpp"
#include "../utility.hpp"

namespace vg {
namespace algorithms {

using namespace std;

// Make sure all of vg's print operators are available.
using vg::operator<<;

#define debug_chaining
#define debug_reseeding



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
            score(item),
            true
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
        return ss.str().substr(0, graph_length(item));
    }
    
    /// Return true if the first item is a perfect chain that can't be beat,
    /// and false otherwise.
    virtual bool has_perfect_chain(const VectorView<Item>& items) const {
        return false;
    }
    
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
     * Get the distance in the read between two items, or 0 if they abut.
     * Returns std::numeric_limits<size_t>::max() if they overlap.
     */
    virtual size_t get_read_distance(const Item& left,
                                     const Item& right) const {
        size_t l = read_end(left);
        size_t r = read_start(right);
        if (r < l) {
            return std::numeric_limits<size_t>::max();
        }
        return r - l;
    }
    
    /**
     * Get the score to make a transition from one item to another, or
     * std::numeric_limits<int>::min() if the transition should be disallowed.
     *
     * Default implementation: highest score in points that the alignment of
     * the region between two Items could possibly get. Should assume that any
     * bases that might match actually match.
     *
     * If the items overlap in the read or the graph, returns
     * std::numeric_limits<int>::min() because one side of the alignment would
     * have negative length.
     *
     * If a required gap would be longer than max_gap_length, returns
     * std::numeric_limits<int>::min() (because we're probably going the long
     * way around an inversion or something.
     */
    virtual int transition_score(const Item& left,
                                 const Item& right,
                                 size_t max_gap_length = std::numeric_limits<size_t>::max()) const {
                                             
        if (read_end(left) == read_start(right)) {
            // They abut in the read.
            
            if (!graph) {
                // Assume they abut in the graph too.
                return 0;
            } else {
                // Do they overlap in the graph?
                size_t graph_overlap = get_graph_overlap(left, right);
                if (graph_overlap > 0) {
#ifdef debug_chaining
                    cerr << "Abutting items in read overlap by " << graph_overlap << " bp in graph" << endl;
#endif
                    return numeric_limits<int>::min();
                }
                
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
                    if (graph_distance > max_gap_length) {
                        // This deletion is too long to allow the transition.
#ifdef debug_chaining
                        std::cerr << "Indel is too long to allow transition" << std::endl;
#endif
                        return numeric_limits<int>::min();
                    }
                    return -(scoring.gap_open + (graph_distance - 1) * scoring.gap_extension);
                }
            }
        } else if (read_end(left) < read_start(right)) {
            // There's a gap in the read
            size_t read_distance = get_read_distance(left, right);
            
            if (!graph) {
                // Assume they are the same distance in the graph.
                // Treat every base as a possible match.
                return scoring.match * read_distance;
            } else {
                // Do they overlap in the graph?
                size_t graph_overlap = get_graph_overlap(left, right);
                if (graph_overlap > 0) {
#ifdef debug_chaining
                    cerr << "Items that are " << read_distance << " bp apart in read overlap by " << graph_overlap << " bp in graph" << endl;
#endif
                    return numeric_limits<int>::min();
                }
            
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

                    if (length_change > max_gap_length) {
                        // This indel is too long to allow the transition.
#ifdef debug_chaining
                        std::cerr << "Indel is too long to allow transition" << std::endl;
#endif
                        return numeric_limits<int>::min();
                    } 
                    
                    // The number of possible matches gives us a base score
                    int jump_points = scoring.match * possible_matches;
                    if (length_change > 0) {
                        // We have to charge for a gap though
                        jump_points -= (scoring.gap_open + (length_change - 1) * scoring.gap_extension);
                    }
                    return jump_points;
                }
            }
        } else {
            // If there's an overlap in the read. say we can't do it.
            
#ifdef debug_chaining
            cerr << "Items overlap in read" << endl;
#endif
            
            return numeric_limits<int>::min();
        }
    }
    
    /**
     * Run checks to make sure the item is self-consistent.
     * Graph must be set.
     */
    virtual void validate(const Item& item, const std::string& full_read_sequence) const {
        // By default, do nothing
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

/**
 * This is how you chain seeds that come from sources.
 */
template<typename Item, typename Source>
struct SourceChainingSpace: public BaseChainingSpace<Item> {
    // These seeds can't really be interpreted without their sources, which
    // they reference by index. Source is assumed to be MinimizerIndex::Minimizer.
    const VectorView<Source> sources;
    
    SourceChainingSpace(const VectorView<Source>& sources,
                        const Aligner& scoring,
                        const SnarlDistanceIndex* distance_index,
                        const HandleGraph* graph) :
        BaseChainingSpace<Item>(scoring, distance_index, graph), sources(sources) {
        
        // Nothing to do!
    }
    
    // API for sources
    
    virtual size_t source_read_start(const Source& source) const = 0;
    
    virtual size_t source_read_end(const Source& source) const = 0;
    
    virtual size_t source_read_length(const Source& source) const = 0;
    
    // Default implementation implements everything read-related on top of that.
    
    virtual size_t read_start(const Item& item) const {
        return source_read_start(this->sources[item.source]);
    }
    
    virtual size_t read_end(const Item& item) const {
        return source_read_end(this->sources[item.source]);
    }
    
    virtual size_t read_length(const Item& item) const {
        return source_read_length(this->sources[item.source]);
    }
};

/**
 * This is how you chain minimizer-based seeds
 */
template<typename Item, typename Source>
struct MinimizerSourceChainingSpace : public SourceChainingSpace<Item, Source> {

    // Source is assumed to be MinimizerIndex::Minimizer.

    MinimizerSourceChainingSpace(const VectorView<Source>& sources,
                  const Aligner& scoring,
                  const SnarlDistanceIndex* distance_index,
                  const HandleGraph* graph) :
        SourceChainingSpace<Item, Source>(sources, scoring, distance_index, graph) {
        
        // Nothing to do!
    }
    
    /// Score items flat by minimizer length
    virtual int score(const Item& item) const {
        return this->scoring.match * this->sources[item.source].length;
    }
    
    /// Score transitions just as indel costs, with no points for implied matches/mismatches
    virtual int transition_score(const Item& left,
                                 const Item& right,
                                 size_t max_gap_length = std::numeric_limits<size_t>::max()) const {
                                 
        if (!this->graph) {
            // No graph means no known indel
            return 0;
        }
        
        size_t read_distance = this->get_read_distance(left, right);
        if (read_distance == numeric_limits<size_t>::max()) {
            // Overlap in read, so not allowed.
            return std::numeric_limits<int>::min();
        }
        
        size_t graph_distance = this->get_graph_distance(left, right);
        if (graph_distance == numeric_limits<size_t>::max()) {
            // No graph connection
            return std::numeric_limits<int>::min();
        }
        
        // Decide how much lenght changed
        size_t indel_length = (read_distance > graph_distance) ? read_distance - graph_distance : graph_distance - read_distance;
        
        // Then charge for that indel
        int score = 0;
        if (indel_length > 0) {
            score -= this->scoring.gap_open;
            if (indel_length > 1) {
                score -= this->scoring.gap_extension * (indel_length - 1);
            }
        }
        return score;
    }
    
    // The seeds know either a start or an end position, depending on
    // minimizer orientation, and a corresponding graph position.
    // So we need to work out which of the start or end is actually seeded, and
    // walk out until we use the rest of the graph node or the whole k-kmer,
    // and that's our seed. 
    
    virtual size_t source_read_start(const Source& source) const {
        if (source.value.is_reverse) {
            // If we're reverse, the start is gotten by walking from the end.
            return this->source_read_end(source) - this->source_read_length(source);
        } else {
            // If we're forward, the start is the stored offset, and we walk right.
            return source.value.offset;
        }
    };
    
    virtual size_t source_read_end(const Source& source) const {
        if (source.value.is_reverse) {
            // If we're reverse, the past-end is the stored offset + 1, and we walk left.
            return source.value.offset + 1;
        } else {
            // If we're forward, the end is gotten by walking from the start.
            return this->source_read_start(source) + this->source_read_length(source);
        }
    };
    
    virtual size_t source_read_length(const Source& source) const {
        // No graph position without an item so just treat the whole thing as the source.
        return source.length;
    };
    
    size_t read_start(const Item& item) const {
        if (this->sources[item.source].value.is_reverse) {
            // If we're reverse, the start is gotten by walking from the end.
            return this->read_end(item) - this->read_length(item);
        } else {
            // If we're forward, the start is the stored offset, and we walk right.
            return this->sources[item.source].value.offset;
        }
    }
    
    size_t read_end(const Item& item) const {
        if (this->sources[item.source].value.is_reverse) {
            // If we're reverse, the past-end is the stored offset + 1, and we walk left.
            return this->sources[item.source].value.offset + 1;
        } else {
            // If we're forward, the end is gotten by walking from the start.
            return this->read_start(item) + this->read_length(item);
        }
    }
    
    size_t read_length(const Item& item) const {
        if (this->graph) {
            return this->graph_length(item);
        } else {
            // No graph so just treat the whole thing as the item.
            return this->sources[item.source].length;
        }
    }
    
    size_t graph_length(const Item& item) const {
        if (this->sources[item.source].value.is_reverse) {
            // If we're reverse, we know the matching at the end.
            handle_t end_handle = this->graph_path_at(item, 0);
            // Length in graph is min of either used node or minimizer length
            return std::min((size_t) this->sources[item.source].length,
                            offset(item.pos) + 1);
        } else {
            // If we're forward, we know the matching at the start
            handle_t start_handle = graph_path_at(item, 0);
            // Length in graph is min of either remaining node or minimizer length
            return std::min((size_t) this->sources[item.source].length,
                            this->graph->get_length(start_handle) - offset(item.pos));
        }
        
    }
    
    pos_t graph_start(const Item& item) const {
        if (this->sources[item.source].value.is_reverse) {
            // If we're reverse, we have the end position in the graph and need to walk it back.
            pos_t end = this->graph_end(item);
            return make_pos_t(id(end), is_rev(end), offset(end) - this->graph_length(item));
        } else {
            // If we're forward, we have the start position in the graph
            return item.pos;
        }
    }
    
    pos_t graph_end(const Item& item) const {
        if (this->sources[item.source].value.is_reverse) {
            // If we're reverse, we have the end position and need to make it a past-end.
            return make_pos_t(id(item.pos), is_rev(item.pos), offset(item.pos) + 1);
        } else {
            // If we're forward, we have the start position and need to walk it forward.
            pos_t start = this->graph_start(item);
            return make_pos_t(id(start), is_rev(start), offset(start) + this->graph_length(item));
        }
    }
    
    size_t graph_path_size(const Item& item) const {
        return 1;
    }
    
    handle_t graph_path_at(const Item& item, size_t index) const {
        return this->graph->get_handle(id(item.pos), is_rev(item.pos));
    }
    
    size_t graph_path_offset(const Item& item) const {
        return offset(this->graph_start(item));
    }
    
    void validate(const Item& item, const std::string& full_read_sequence) const {
        string read_seq = this->get_read_sequence(item, full_read_sequence);
        string graph_seq = this->get_graph_sequence(item);
        string key_seq = this->sources[item.source].forward_sequence();
        if (read_seq.empty() ||
            read_seq != graph_seq || !(
                this->sources[item.source].value.is_reverse ?
                std::equal(read_seq.rbegin(), read_seq.rend(), key_seq.rbegin())
                : std::equal(read_seq.begin(), read_seq.end(), key_seq.begin())
        )) {
            // Read and graph sequence have to match, and that sequence has to
            // be a prefix or suffix of the read-orientation minimizer key, as
            // appropriate.
            throw std::runtime_error("Sequence mismatch for hit of source " + std::to_string(item.source) + "/" + std::to_string(this->sources.size()) + " on oriented key " + key_seq + " with read " + read_seq + " and graph " + graph_seq);
        }
    }
};

/// This is how you chain up new seeds
template<typename Source>
struct ChainingSpace<NewSnarlSeedClusterer::Seed, Source> : public MinimizerSourceChainingSpace<NewSnarlSeedClusterer::Seed, Source> {
    using Item = NewSnarlSeedClusterer::Seed;
    
    ChainingSpace(const VectorView<Source>& sources,
                  const Aligner& scoring,
                  const SnarlDistanceIndex* distance_index,
                  const HandleGraph* graph) :
        MinimizerSourceChainingSpace<Item, Source>(sources, scoring, distance_index, graph) {
        
        // Nothing to do!
    }
    
    
};

/// This is how you chain up old seeds
template<typename Source>
struct ChainingSpace<SnarlSeedClusterer::Seed, Source> : public MinimizerSourceChainingSpace<SnarlSeedClusterer::Seed, Source> {
    using Item = SnarlSeedClusterer::Seed;
    
    ChainingSpace(const VectorView<Source>& sources,
                  const Aligner& scoring,
                  const HandleGraph* graph) :
        MinimizerSourceChainingSpace<Item, Source>(sources, scoring, nullptr, graph) {
        
        // Nothing to do!
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

// These next functions all have to be templated on "Collection" if we want to
// ever be able to call them with vectors. An implicit conversion to
// VectorView doesn't count when type deduction is happening, so we can;t just
// use VectorView as the collection type. 

/**
 * Fill in the given DP table for the best chain score ending with each
 * item. Returns the best observed score overall from that table,
 * with provenance to its location in the table, if tracked in the type.
 * Assumes some items exist.
 *
 * Input items must be sorted by start position in the read.
 *
 * Takes the given per-item bonus for each item collected.
 *
 * Uses a finite lookback in items and in read bases when checking where we can
 * come from to reach an item. Also, once a given number of reachable items
 * have been found, stop looking back.
 *
 * Limits transitions to those involving indels of the given size or less, to
 * avoid very bad transitions still counting as "reachable".
 */
template<typename Score, typename Item, typename Source = void, typename Collection = VectorView<Item>>
Score chain_items_dp(vector<Score>& best_chain_score,
                     const Collection& to_chain,
                     const ChainingSpace<Item, Source>& space,
                     int item_bonus = 0,
                     size_t lookback_items = 500,
                     size_t lookback_bases = 1000,
                     size_t lookback_reachable_items = 5,
                     size_t max_indel_bases = 10000);

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
 * traceback that can be obtained by chaining items together.
 *
 * Input items must be sorted by start position in the read.
 *
 * Returns the score and the list of indexes of items visited to achieve
 * that score, in order.
 */
template<typename Item, typename Source = void, typename Collection = VectorView<Item>>
pair<int, vector<size_t>> find_best_chain(const Collection& to_chain,
                                          const algorithms::ChainingSpace<Item, Source>& space);

/**
 * Score the given group of items. Determines the best score that can be
 * obtained by chaining items together.
 *
 * Input items must be sorted by start position in the read.
 */
template<typename Item, typename Source = void, typename Collection = VectorView<Item>>
int score_best_chain(const Collection& to_chain,
                     const ChainingSpace<Item, Source>& space);
                     
/**
 * Look at the sources which occur between the given items in the read,
 * and locate their hits that occur in the part of the graph between the
 * given items in the graph, if any. Graph search is limited by
 * max_fallow_search_distance. 
 *
 * Takes a lazy place to put an inverse of the source score
 * order sort (which is the space they are numbered in for item
 * references), from the space's sources view.
 *
 * Also takes a function which, when given a Source and the sorted IDs of a
 * subgraph, calls back with each pos_t in the subgraph that can correspond to
 * the source to form an Item.
 *
 * Returns new, forged items in an arbitrary order.
 *
 * Forged items will not have all fields set. They will only have locations and
 * sources.
 *
 * TODO: Use a different kind of item type to avoid this forgery!
 */
template<typename Item, typename Source>
vector<Item> reseed_fallow_region(const Item& left,
                                  const Item& right,
                                  const ChainingSpace<Item, Source>& space,
                                  std::unique_ptr<VectorViewInverse>& source_sort_inverse,
                                  const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                  size_t max_fallow_search_distance = 10000);

/**
 * Reseed all fallow regions (regions in the read where the distance
 * between items is at least fallow_region_size) by locating hits of the
 * given sources that occur in the part of the graph between existing
 * items.
 *
 * Takes a lazy place to put an inverse of the source score
 * order sort (which is the space they are numbered in for item
 * references), from the space's sources view.
 *
 * Also takes a function which, when given a Source and the sorted IDs of a
 * subgraph, calls back with each pos_t in the subgraph that can correspond to
 * the source to form an Item.
 *
 * item_storage is a collection of items, and sorted_item_indexes is index
 * numbers in that collection of items, sorted by start position in the read.
 *
 * Updates item_storage with newely-found forged items (which only have
 * their position and source fields set), and updates
 * sorted_item_indexes to sort the newly expanded list of items in read
 * order.
 *
 * Returns some statistics: the number of fallow regions found, and the length
 * of the longest one.
 */
template<typename Item, typename Source>
pair<size_t, size_t> reseed_fallow_regions(vector<Item>& item_storage,
                                           vector<size_t>& sorted_item_indexes,
                                           const ChainingSpace<Item, Source>& space,
                                           std::unique_ptr<VectorViewInverse>& source_sort_inverse,
                                           const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                           size_t max_fallow_search_distance = 10000,
                                           size_t fallow_region_size = 200);

// --------------------------------------------------------------------------------

// Template implementations

template<typename Score, typename Item, typename Source, typename Collection>
Score chain_items_dp(vector<Score>& best_chain_score,
                     const Collection& to_chain,
                     const ChainingSpace<Item, Source>& space,
                     int item_bonus,
                     size_t lookback_items,
                     size_t lookback_bases,
                     size_t lookback_reachable_items,
                     size_t max_indel_bases) {
    
    DiagramExplainer diagram({{"rankdir", "LR"}});
    
    // Grab the traits into a short name so we can use the accessors concisely.
    using ST = score_traits<Score>;

#ifdef debug_chaining
    cerr << "Chaining group of " << to_chain.size() << " items" << endl;
#endif
    
    // We want to consider all the important transitions in the graph of what
    // items can come before what other items. We aren't allowing any
    // transitions between items that overlap in the read. We're going through
    // the destination items in order by read start, so we should also keep a
    // list of them in order by read end, and sweep a cursor over that, so we
    // always know the fisrt item that overlaps with or passes the current
    // destination item, in the read. Then when we look for possible
    // predecessors of the destination item, we can start just before there and
    // look left.
    vector<size_t> read_end_order = sort_permutation(to_chain.begin(), to_chain.end(), [&](const Item& a, const Item& b) {
        return space.read_end(a) < space.read_end(b);
    });
    // We use first overlapping instead of last non-overlapping because we can
    // just initialize first overlapping at the beginning and be right.
    auto first_overlapping_it = read_end_order.begin();
    
    // Make our DP table big enough
    best_chain_score.resize(to_chain.size(), ST::unset());
    
    // What's the winner so far?
    Score best_score = ST::unset();
    
    for (size_t i = 0; i < to_chain.size(); i++) {
        // For each item
        auto& here = to_chain[i];
        
        while (space.read_end(to_chain[*first_overlapping_it]) <= space.read_start(here)) {
            // Scan ahead through non-overlapping items that past-end too soon,
            // to the first overlapping item that ends earliest.
            // Ordering physics *should* constrain the iterator to not run off the end.
            ++first_overlapping_it;
            assert(first_overlapping_it != read_end_order.end());
        }
        
#ifdef debug_chaining
        cerr << "First item overlapping item " << i << " beginning at " << space.read_start(here) << " is item " << *first_overlapping_it << " past-ending at " << space.read_end(to_chain[*first_overlapping_it]) << std::endl;
#endif
        
        // How many points is it worth to collect?
        auto item_points = space.score(here) + item_bonus;
        
        std::string here_gvnode = "i" + std::to_string(i);
        
        // If we come from nowhere, we get those points.
        best_chain_score[i] = std::max(best_chain_score[i], ST::annotate(item_points, ST::nowhere()));
        
#ifdef debug_chaining
        cerr << "Look at transitions to item " << i
            << " at " << space.read_start(here) << " to " << space.read_end(here);
        if (space.graph) {
            cerr << " = " << space.graph_start(here);
        }
        cerr << endl;
#endif
        
        // Count how many places we evaluated that we could have come from.
        // We may want to limit this to prevent long gaps.
        size_t reachable_items_found = 0;
        // Count how many total source items were considered for reachability
        size_t lookback_items_tested = 0;
        auto predecessor_index_it = first_overlapping_it;
        while (predecessor_index_it != read_end_order.begin() &&
               reachable_items_found < lookback_reachable_items &&
               lookback_items_tested < lookback_items) {
            --predecessor_index_it;
            // For each source that ended before here started, in reverse order by end position...
            auto& source = to_chain[*predecessor_index_it];
            
            // Say we are considering it
            lookback_items_tested++;
            
#ifdef debug_chaining
            cerr << "\tConsider transition from item " << *predecessor_index_it 
            << " at " << space.read_start(source) << " to " << space.read_end(source); 
            if (space.graph) {
                cerr << " = " << space.graph_start(source);
            }
            cerr << " score " << best_chain_score[*predecessor_index_it]
                << " with read distance " << space.get_read_distance(source, here);
            if (space.graph) {
                cerr << " and graph distance " << space.get_graph_distance(source, here);
            }
            cerr << endl;
#endif

            if (lookback_bases != 0 && space.get_read_distance(source, here) > lookback_bases) {
                // This one is out of range in the read.
                
#ifdef debug_chaining
                cerr << "\t\tStop because the read distance is longer than the lookback threshold." << endl;
#endif
                // Read distance is monotonic because we are going by right edge.
                // So we can stop now.
                break;
            }
            
            // How much does it pay (+) or cost (-) to make the jump from there
            // to here?
            // Don't allow the transition if it seems like we're going the long
            // way around an inversion and needing a huge indel.
            int jump_points = space.transition_score(source, here, max_indel_bases);
            
            if (jump_points != numeric_limits<int>::min()) {
                // Get the score we are coming from
                typename ST::Score source_score = ST::score_from(best_chain_score, *predecessor_index_it);
                
                // And the score with the transition and the points from the item
                typename ST::Score from_source_score = ST::add_points(source_score, jump_points + item_points);
                
                // Remember that we could make this jump
                best_chain_score[i] = std::max(best_chain_score[i],
                                               from_source_score);
                                               
#ifdef debug_chaining
                cerr << "\t\tWe can reach " << i << " with " << from_source_score << " which is " << source_score << " and a transition of " << jump_points << " to collect " << item_points << endl;
#endif
                std::string source_gvnode = "i" + std::to_string(*predecessor_index_it);
                // Suggest that we have an edge, where the edges that are the best routes here are the most likely to actually show up.
                diagram.suggest_edge(source_gvnode, here_gvnode, here_gvnode, ST::score(from_source_score), {
                    {"label", std::to_string(jump_points)},
                    {"weight", std::to_string(std::max<int>(1, ST::score(from_source_score)))}
                });

                reachable_items_found++;
            }
        }
        
#ifdef debug_chaining
        cerr << "\tBest way to reach " << i << " is " << best_chain_score[i] << endl;
#endif
        
        diagram.add_node(here_gvnode, {
            {"label", std::to_string(i) + " = " + std::to_string(item_points) + "/" + std::to_string(ST::score(best_chain_score[i]))}
        });
        if (space.graph) {
            auto graph_start = space.graph_start(here);
            std::string graph_gvnode = "n" + std::to_string(id(graph_start));
            diagram.ensure_node(graph_gvnode, {
                {"label", std::to_string(id(graph_start))},
                {"shape", "box"}
            });
            diagram.add_edge(here_gvnode, graph_gvnode, {{"color", "gray"}});
        }
        
        // See if this is the best overall
        ST::max_in(best_score, best_chain_score, i);
        
#ifdef debug_chaining
        cerr << "\tBest chain end so far: " << best_score << endl;
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
        traced_score_t best_past_ending_score_ever = chain_items_dp(best_chain_score,
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
        return algorithms::chain_items_dp(best_chain_score,
                                          to_chain,
                                          space);
    }
}

template<typename Item, typename Source>
vector<Item> reseed_fallow_region(const Item& left,
                                  const Item& right,
                                  const algorithms::ChainingSpace<Item, Source>& space,
                                  std::unique_ptr<VectorViewInverse>& source_sort_inverse,
                                  const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                  size_t max_fallow_search_distance) {
    // We have two seeds here, so there aren't no sources.
    // We know there must be a collection of sources in read order, so grab that.
    const vector<Source>& sources_in_read_order = *space.sources.items;
    
    // Find the range of indexes in read order for the source hits between
    // those of the bounding items.
    size_t range_begin = space.sources.backing_index(left.source) + 1;
    size_t range_end = space.sources.backing_index(right.source);
    
    if (range_begin >= range_end) {
        // Not actually any sources between here in the read.
        return {};
    }
    
    // Query distance index to see if the items are actually plausibly reachable
    size_t graph_min_distance;
    if (space.distance_index) {
        // We have a (new) graph distance index, so use that.
        // TODO: support old graph distance index???
        graph_min_distance = space.get_graph_distance(left, right);
    } else {
        // No graph distances. Just search the limit.
        graph_min_distance = max_fallow_search_distance;
    }
    
    if (graph_min_distance == std::numeric_limits<size_t>::max()) {
        // Actually unreachable, so don't add any items.
        // TODO: detect and warn? Try different bounds?
        return {};
    }
    
    if (graph_min_distance > max_fallow_search_distance) {
        // This would be too far to search.
        // TODO: warn!
        return {};
    }
    
    // Decide how far we want to search
    size_t graph_distance = std::min(max_fallow_search_distance, graph_min_distance * 2 + 1000);
    
    // Extract the connecting graph. We need the IDs in sorted order.
    vector<nid_t> sorted_ids;
    {
        // Get the graph bounds we care about looking between
        const pos_t& left_bound = space.graph_end(left);
        const pos_t& right_bound = space.graph_start(right);
        // TODO: Add an algorithm version that doesn't bother actually extracting?
        HashGraph connecting_graph;
        unordered_map<nid_t, nid_t> extracted_to_original = algorithms::extract_connecting_graph(space.graph, &connecting_graph, graph_distance, left_bound, right_bound);
        sorted_ids.reserve(extracted_to_original.size());
        for (auto& kv : extracted_to_original) {
            // Save the original-graph node IDs
            sorted_ids.push_back(kv.second);
        }
    }
    std::sort(sorted_ids.begin(), sorted_ids.end());
    
#ifdef debug_reseeding
    cerr << "Reseeding " << (range_end - range_begin) << " sources against " << sorted_ids.size() << " graph nodes" << endl;
#endif
    
    // Find hits on these nodes, for the sources that are in the right part of the read, and forge items for them.
    vector<Item> forged_items;
    for (size_t i = range_begin; i < range_end; i++) {
        // For each source between the bounds
        const Source& m = sources_in_read_order[i];
        
        // We may see duplicates, so we want to do our own deduplication.
        unordered_set<pos_t> seen;
        
        // Find all its hits in the part of the graph between the bounds
        for_each_pos_for_source_in_subgraph(m, sorted_ids, [&](const pos_t& pos) {
            // So now we know pos corresponds to read base
            // m.value.offset, in the read's forward orientation.
            
            if (seen.count(pos)) {
                // This is a duplicate position for this source, so skip it.
                // TODO: Why do these happen?
                return;
            }
            seen.insert(pos);
            
            // Forge an item.
            forged_items.emplace_back();
            forged_items.back().pos = pos;
            // Run the read-order index it through the inverse permutation to
            // get the number everyone else will know it by.
            if (!source_sort_inverse) {
                // But lazily make sure we have the inverse permutation first.
                source_sort_inverse = std::make_unique<VectorViewInverse>(space.sources);
            }
            forged_items.back().source = (*source_sort_inverse)[i];
            
#ifdef debug_reseeding
            cerr << "Found new seed for read-order source " << i << " of read@" << space.read_start(forged_items.back()) << " = graph@" << space.graph_start(forged_items.back()) << endl;
#endif
            
        });
    }

    return forged_items;
}

template<typename Item, typename Source>
pair<size_t, size_t> reseed_fallow_regions(vector<Item>& item_storage,
                                           vector<size_t>& sorted_item_indexes,
                                           const algorithms::ChainingSpace<Item, Source>& space,
                                           std::unique_ptr<VectorViewInverse>& source_sort_inverse,
                                           const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                           size_t max_fallow_search_distance,
                                           size_t fallow_region_size) {
    pair<size_t, size_t> statistics {0, 0};
    auto& fallow_region_count = statistics.first;
    auto& longest_fallow_region_length = statistics.second;
    
    // Make a VectorView over the items
    VectorView<Item> item_view {item_storage, sorted_item_indexes};
    
    for (size_t left = 0; left + 1 < item_view.size(); left++) {
        // For each pair...
        size_t right = left + 1;
        
        // Check read distance.
        // TODO: Could we use the faster non-space way of doing this without ever touching the graph?
        size_t read_distance = space.get_read_distance(item_view[left], item_view[right]);
        
        if (read_distance > fallow_region_size) {
            // If a pair is too far apart
            
#ifdef debug_reseeding
            cerr << "Reseeding between seeds " << item_view.backing_index(left)
                << " at read index " << space.read_start(item_view[left])
                << " and " << item_view.backing_index(right)
                << " at read index " << space.read_start(item_view[right])
                << endl;
#endif
            
            // Count it
            fallow_region_count++;
            longest_fallow_region_length = std::max(longest_fallow_region_length, read_distance);
            
            // Forge some fill-in items
            vector<Item> new_items = reseed_fallow_region(item_view[left],
                                                          item_view[right],
                                                          space,
                                                          source_sort_inverse,
                                                          for_each_pos_for_source_in_subgraph,
                                                          max_fallow_search_distance);
        
            // Append the fill-in items onto the item storage vector. Safe to do
            // even when we're using a view over it.
            item_storage.reserve(item_storage.size() + new_items.size());
            std::copy(new_items.begin(), new_items.end(), std::back_inserter(item_storage));
        }
    }
    
    if (item_storage.size() != sorted_item_indexes.size()) {
        // After filling in all the fallow regions, if we actually got any new
        // items, extend and re-sort sorted_item_indexes.
        // TODO: Can we just build this in order as we go left to right along the read actually?
        sorted_item_indexes.reserve(item_storage.size());
        while(sorted_item_indexes.size() < item_storage.size()) {
            sorted_item_indexes.push_back(sorted_item_indexes.size());
        }
        
        // TODO: De-duplicate sort with the initial sort
        std::sort(sorted_item_indexes.begin(), sorted_item_indexes.end(), [&](const size_t& a, const size_t& b) -> bool {
            return space.read_start(item_storage[a]) < space.read_start(item_storage[b]);
        });
    }
    
    return statistics;
}

}
}

#endif

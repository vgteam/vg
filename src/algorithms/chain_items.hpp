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
 
#include "extract_containing_graph.hpp"

#include "../gbwt_extender.hpp"
#include "../snarl_seed_clusterer.hpp"
#include "../snarl_seed_clusterer.hpp"
#include "../handle.hpp"
#include "../explainer.hpp"
#include "../utility.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace algorithms {

using namespace std;

// Make sure all of vg's print operators are available.
using vg::operator<<;

//#define debug_chaining
//#define debug_reseeding


/// We support chaining different kinds of things, so we have a type that
/// abstracts out accessing their chaining-relevant fields and measuring
/// distance between them.
/// We support things that come from other source things (i.e. minimizers).
/// See BaseChainingSpace for the actual interface.
template<typename Item, typename Source = void>
struct ChainingSpace {
    // Nothing here; this exists to be specialized.
};

/**
 * Class to represent a scoring regime for chaining. Includes the scoring
 * parameters from an Aligner, and also a strategy for scoring jumps between
 * items.
 */
struct ChainingScorer {

    /**
     * Make a ChainingScorer off of an Aligner.
     */
    inline ChainingScorer(const Aligner& scores = Aligner()) : match(scores.match), mismatch(scores.mismatch), gap_open(scores.gap_open), gap_extension(scores.gap_extension) {
        // Nothing to do!
    }
    
    virtual ~ChainingScorer() = default;
    
    /**
     * Score a transition between two items, based on read and graph distance.
     * If read distance is std::numeric_limits<size_t>::max(), the items overlap in the read.
     * If graph distance is std::numeric_limits<size_t>::max(), the items are unreachable in the graph.
     * If graph distance is null, no graph is available.
     * Returns a score for the transition, or std::numeric_limits<int>::min() if the transition should be prohibited.
     * Prohibits any transitions with a length change over max_gap_length.
     */
    virtual int score_transition(size_t read_distance,
                                 const size_t* graph_distance,
                                 size_t max_gap_length = std::numeric_limits<size_t>::max()) const = 0;



    int match;
    int mismatch;
    int gap_open;
    int gap_extension;

};

/**
 * Scorer that assumes all non-indel bases are matches.
 */
struct MatchAssumingChainingScorer : public ChainingScorer {

    using ChainingScorer::ChainingScorer;
    
    virtual ~MatchAssumingChainingScorer() = default;
    
    /**
     * Score a transition between two items, based on read and graph distance.
     * If read distance is std::numeric_limits<size_t>::max(), the items overlap in the read.
     * If graph distance is std::numeric_limits<size_t>::max(), the items are unreachable in the graph.
     * If graph distance is null, no graph is available.
     * Returns a score for the transition, or std::numeric_limits<int>::min() if the transition should be prohibited.
     * Prohibits any transitions with a length change over max_gap_length.
     */
    virtual int score_transition(size_t read_distance,
                                 const size_t* graph_distance,
                                 size_t max_gap_length = std::numeric_limits<size_t>::max()) const;
                                 
};

/**
 * Scorer that only counts the scores of indels.
 */
struct IndelOnlyChainingScorer : public ChainingScorer {

    using ChainingScorer::ChainingScorer;
    
    virtual ~IndelOnlyChainingScorer() = default;

    /**
     * Score a transition between two items, based on read and graph distance.
     * If read distance is std::numeric_limits<size_t>::max(), the items overlap in the read.
     * If graph distance is std::numeric_limits<size_t>::max(), the items are unreachable in the graph.
     * If graph distance is null, no graph is available.
     * Returns a score for the transition, or std::numeric_limits<int>::min() if the transition should be prohibited.
     * Prohibits any transitions with a length change over max_gap_length.
     */
    virtual int score_transition(size_t read_distance,
                                 const size_t* graph_distance,
                                 size_t max_gap_length = std::numeric_limits<size_t>::max()) const;
};

/**
 * Source for distances between positions in the graph. May do caching.
 */
class DistanceSource {
public:

    virtual ~DistanceSource() = default;

    /**
     * Query the distance between two positions. Not const, because caches may
     * need to be updated.
     *
     * Returns std::numeric_limits<size_t>::max() if unreachable.
     */
    virtual size_t get_distance(const pos_t& left, const pos_t& right) = 0;
};

/**
 * Distance source that pulls straight from a SnarlDistanceIndex.
 */
class DistanceIndexDistanceSource : public DistanceSource {
public:

    DistanceIndexDistanceSource(const SnarlDistanceIndex* distance_index, const HandleGraph* graph);
    virtual ~DistanceIndexDistanceSource() = default;

    virtual size_t get_distance(const pos_t& left, const pos_t& right);
   
protected:
    const SnarlDistanceIndex* distance_index;
    const HandleGraph* graph;
};

/**
 * Distance source that pulls from a pre-computed graph of oriented distances
 * between immediate neighbors. Answers queries by doing an algorithmically bad
 * but in practice probably fast enough Dijkstra every time.
 */
class DistanceNetDistanceSource : public DistanceSource {
public:
    /**
     * Make a distance source by exploring the given graph around the given
     * relevant positions. Only queries between the given positions are
     * guaranteed to be accepted. The same position may appear multiple times.
     *
     * Distances larger than the given max distance between adjacent positions
     * may not be recorded.
     */
    DistanceNetDistanceSource(const HandleGraph* graph, const vector<pos_t>& relevant_positions, size_t max_step_distance = 200);
    
    virtual ~DistanceNetDistanceSource() = default;

    virtual size_t get_distance(const pos_t& left, const pos_t& right);
    
protected:
    std::unordered_map<std::pair<nid_t, bool>, std::unordered_map<std::pair<nid_t, bool>, size_t>> immediate_distances;
    std::unordered_map<nid_t, size_t> node_lengths;
};

/// Baser base class to let you hold chaining spaces over anything, and destroy
/// them without knowing what they are over
struct UnknownItemChainingSpace {
    virtual ~UnknownItemChainingSpace() = default;
};

/// Base class of specialized ChainingSpace classes. Defines the ChainingSpace interface.
template<typename Item>
struct BaseChainingSpace : public UnknownItemChainingSpace {
    
    const ChainingScorer& scorer;
    
    const SnarlDistanceIndex* distance_index;
    const HandleGraph* graph;
    
    mutable DistanceSource* distance_source;
    unique_ptr<DistanceSource> distance_source_storage;
    
    BaseChainingSpace(const ChainingScorer& scorer,
                      const SnarlDistanceIndex* distance_index = nullptr,
                      const HandleGraph* graph = nullptr,
                      DistanceSource* distance_source = nullptr) :
        scorer(scorer),
        distance_index(distance_index),
        graph(graph) {
        
        if (distance_source == nullptr && distance_index != nullptr && graph != nullptr) {
            // Default the distance source
            distance_source_storage.reset(new DistanceIndexDistanceSource(distance_index, graph));
            this->distance_source = distance_source_storage.get();
        } else {
            this->distance_source = distance_source;
        }
    }
    
    virtual ~BaseChainingSpace() = default;
    
    /// Replace our current distance source with this new one, which we take
    /// ownership of.
    virtual void give_distance_source(DistanceSource* new_distance_source) {
        distance_source_storage.reset(new_distance_source);
        distance_source = new_distance_source;
    }
    
    /// Get the score collected by visiting the item.
    virtual int score(const Item& item) const {
        // Default implementation assumes a perfect match
        return scorer.match * read_length(item);
    }
    
    /// Get the alignment score for the alignment represented by an item. May
    /// be different han the score used in chaining.
    virtual int alignment_score(const Item& item) const {
        // Default implementation assumes a perfect match
        return scorer.match * read_length(item);
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
            alignment_score(item),
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
    
        if (!distance_source) {
            return numeric_limits<size_t>::max();
        }
        
        // Find where the left extension past-ends
        pos_t left_past_end = graph_end(left);
        // Find where the right one starts
        pos_t right_start = graph_start(right);
        
        // Ask the distance source about it
        return distance_source->get_distance(left_past_end, right_start);
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
     * Run checks to make sure the item is self-consistent.
     * Graph must be set.
     */
    virtual void validate(const Item& item, const std::string& full_read_sequence) const {
        // By default, do nothing
    }
    
    /**
     * Turn an Item into something we can show to the user.
     */
    virtual std::string to_string(const Item& item) const {
        std::stringstream s;
        s << "{R" << this->read_start(item) << "->" << this->read_end(item);
        if (this->graph) {
            s << "=G" << this->graph_start(item);
        }
        s << "}";
        return s.str();
    }
    
    /**
     * Turn a transition between Items into something we can show to the user.
     * Describe just the transition, not the items.
     */
    virtual std::string to_string(const Item& a, const Item& b) const {
        std::stringstream s;
        s << "{Rdist: " << this->get_read_distance(a, b);
        if (this->graph) {
            s << " Gdist: " << this->get_graph_distance(a, b);
        }
        s << "}";
        return s.str();
    }
};


/// This is how you chain up a bunch of GaplessExtension items
template<>
struct ChainingSpace<GaplessExtension, void>: public BaseChainingSpace<GaplessExtension> {
    using Item = GaplessExtension;
    
    // Keep the constructor
    using BaseChainingSpace<Item>::BaseChainingSpace;
    
    virtual ~ChainingSpace() = default;
    
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
                        const ChainingScorer& scorer,
                        const SnarlDistanceIndex* distance_index,
                        const HandleGraph* graph,
                        DistanceSource* distance_source = nullptr) :
        BaseChainingSpace<Item>(scorer, distance_index, graph, distance_source), sources(sources) {
        
        // Nothing to do!
    }
    
    virtual ~SourceChainingSpace() = default;
    
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
                  const ChainingScorer& scorer,
                  const SnarlDistanceIndex* distance_index,
                  const HandleGraph* graph,
                  DistanceSource* distance_source = nullptr) :
        SourceChainingSpace<Item, Source>(sources, scorer, distance_index, graph, distance_source) {
        
        // Nothing to do!
    }
    
    virtual ~MinimizerSourceChainingSpace() = default;
    
    /// Score items flat by minimizer length
    virtual int score(const Item& item) const {
        return this->scorer.match * this->sources[item.source].length;
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
                  const ChainingScorer& scorer,
                  const SnarlDistanceIndex* distance_index,
                  const HandleGraph* graph,
                  DistanceSource* distance_source = nullptr) :
        MinimizerSourceChainingSpace<Item, Source>(sources, scorer, distance_index, graph, distance_source) {
        
        // Nothing to do!
    }
    
    virtual ~ChainingSpace() = default;
    
    
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
 * When making a chain, we need to be able to control how far back we look for
 * previous items to connect to. For controling this, we have a lookback
 * strategy.
 *
 * To control the strategy, you instantiate a subclass of this calss and configure its fields.
 *
 * To use the strategy, you call setup_problem() on it to get a
 * LookbackStrategy::Problem, which is per-chaining-problem scratch. Then you
 * use that lookback problem while solving your chaining problem, in a single
 * thread, and then throw it away.
 */
class LookbackStrategy {
public:
    /// When asked about a transition, we can check that item, skip that item
    /// and look at more before it, or stop looking.
    enum verdict_t {
        CHECK,
        SKIP,
        STOP
    };
    
    /**
     * Individual problem scratch that gets made from the strategy, like a factory.
     * Strategy must always outlive its problems.
     */
    class Problem {
    public:
        /**
         * Start working on a new destination. Destinations should be visited in
         * increasing order of read start position.
         *
         * Must be called before should_check().
         */
        virtual void advance() = 0;
        
        /**
         * Tell the problem that the given item is at the given graph position.
         * Lets the problem remember this and do a bit of its own clustering.
         */
        virtual void place_in_graph(size_t item, const pos_t& graph_pos, const SnarlDistanceIndex& distance_index, const HandleGraph& graph);
        
        /**
         * Determine if we should consider a transition between two items, given
         * their read distance and the predecessor's achieved score. Counts the
         * transition as a possibility.
         */
        virtual verdict_t should_check(size_t item_a, size_t item_b, size_t read_distance, int item_a_score) = 0;
        
        /**
         * Report that we considered a connection between two items, and found them
         * to be at the given read distance, and the given graph distance, and that
         * the transition had the given score. Distances should be
         * std::numeric_limits<size_t>::max() for things that are unreachable, and
         * the score should be std::numeric_limits<int>::min() if the transition is
         * otherwise disallowed.
         *
         * For a given destination, should be called for each source in decreasing
         * order of end position.
         */
        virtual void did_check(size_t item_a, size_t item_b, size_t read_distance, const size_t* graph_distance, int transition_score, int achieved_score) = 0;
        
        // Since we will use a pointer to this base class to own derived
        // classes, we need a virtual destructor. 
        virtual ~Problem() = default;
    };

    /**
     * Set up a new problem over the given number of items. Thread safe.
     *
     * We need to use a pointer here because the problem scratch for different
     * lookback control strategies can be of different sizes.
     */
    virtual std::unique_ptr<Problem> setup_problem(size_t item_count) const = 0;
};

/**
 * Lookback strategy that allows a certain number of total items, a certain
 * base pair distance, and a certain number of "good" items.
 */
class FlatLimitLookbackStrategy : public LookbackStrategy {
public:
    /// How many items will we look at total?
    size_t lookback_items = 500;
    /// How many good items will we look at before stopping?
    size_t lookback_good_items = 10;
    /// How many reachable items will we look at before stopping?
    size_t lookback_reachable_items = 500;
    
    /// How far back should we look before stopping.
    size_t lookback_bases = 1000;
    
    virtual std::unique_ptr<LookbackStrategy::Problem> setup_problem(size_t item_count) const;
    
    class Problem : public LookbackStrategy::Problem {
    public:
        Problem(const FlatLimitLookbackStrategy& parent);
    
        virtual void advance();
        virtual verdict_t should_check(size_t item_a, size_t item_b, size_t read_distance, int item_a_score);
        virtual void did_check(size_t item_a, size_t item_b, size_t read_distance, const size_t* graph_distance, int transition_score, int achieved_score);
        
        virtual ~Problem() = default;
        
    protected:
        const FlatLimitLookbackStrategy& strategy;
        size_t lookback_items_used = 0;
        size_t lookback_good_items_used = 0;
        size_t lookback_reachable_items_used = 0;
        int best_achieved_score = 0;
    };
};

/**
 * Lookback strategy that progressively doubles the lookback distance until
 * something with a positive overall score is found, or a hard limit is hit.
 */
class ExponentialLookbackStrategy : public LookbackStrategy {
public:
    /// How far back should we look before stopping, max?
    size_t lookback_bases = 200;
    /// How far should our initial search go?
    size_t initial_search_bases = 10;
    /// How much should we increase by?
    double scale_factor = 2.0;
    /// How many points can we lose on a transition per max distance and still have it be good?
    double min_good_transition_score_per_base = -0.1;
    
    virtual std::unique_ptr<LookbackStrategy::Problem> setup_problem(size_t item_count) const;
    
    class Problem : public LookbackStrategy::Problem {
    public:
        Problem(const ExponentialLookbackStrategy& parent);
    
        virtual void advance();
        virtual verdict_t should_check(size_t item_a, size_t item_b, size_t read_distance, int item_a_score);
        virtual void did_check(size_t item_a, size_t item_b, size_t read_distance, const size_t* graph_distance, int transition_score, int achieved_score);
        
        virtual ~Problem() = default;
        
    protected:
        const ExponentialLookbackStrategy& strategy;
        size_t limit;
        int best_transition_found;
        int best_achieved_score;
        bool good_score_found;
    };
};

/**
 * Lookback strategy that uses exponential limits and buckets items by graph reachability.
 */
class BucketLookbackStrategy : public LookbackStrategy {
public:
    /// How far back should we look before stopping, max?
    size_t lookback_bases = 200;
    /// How far should our initial search go?
    size_t initial_search_bases = 10;
    /// How much should we increase by?
    double scale_factor = 2.0;
    /// How many points can we lose on a transition per max distance and still have it be good?
    double min_good_transition_score_per_base = -0.1;
    /// How far apart should things need to be to be different buckets?
    /// TODO: What if the point at which we tick over this happens to divide the true chain?
    size_t bucket_limit = 50000;
    /// How far must an indel go in bucket coordinates to prohibit crossing it? Item lengths can count against this.
    size_t max_inferred_indel = 1100;
    /// How much should the bucket coordinates of two things differ before we decide not to check for a path between them, even when we can't actually lower-bound the distance?
    size_t max_suspicious_bucket_coordinate_difference = 10000;
    
    virtual std::unique_ptr<LookbackStrategy::Problem> setup_problem(size_t item_count) const;
    
    class Problem : public LookbackStrategy::Problem {
    public:
        Problem(const BucketLookbackStrategy& parent);
    
        virtual void advance();
        virtual void place_in_graph(size_t item, const pos_t& graph_pos, const SnarlDistanceIndex& distance_index, const HandleGraph& graph);
        virtual verdict_t should_check(size_t item_a, size_t item_b, size_t read_distance, int item_a_score);
        virtual void did_check(size_t item_a, size_t item_b, size_t read_distance, const size_t* graph_distance, int transition_score, int achieved_score);
        
        virtual ~Problem() = default;
        
    protected:
        const BucketLookbackStrategy& strategy;
        size_t limit;
        int best_transition_found;
        int best_achieved_score;
        bool good_score_found;
        std::vector<pos_t> bucket_heads;
        std::vector<size_t> bucket_sizes;
        std::unordered_map<size_t, size_t> item_to_head;
        std::unordered_map<size_t, size_t> item_to_coordinate;
    };
};

/**
 * Get rid of items that are shadowed or contained by (or are identical to) others.
 * We assume that if the start and end positions put you on the same diagonals
 * on the same nodes, where the middle of the path goes doesn't matter.
 *
 * Erases items that didn't survive from indexes, and sorts them by read start
 * position.
 */
template<typename Item, typename Source = void>
void sort_and_shadow(const std::vector<Item>& items, std::vector<size_t>& indexes, const ChainingSpace<Item, Source>& space);

/**
 * Get rid of items that are shadowed or contained by (or are identical to) others.
 * We assume that if the start and end positions put you on the same diagonals
 * on the same nodes, where the middle of the path goes doesn't matter.
 *
 * Erases items that didn't survive from items, and sorts them by read start
 * position.
 */
template<typename Item, typename Source = void>
void sort_and_shadow(std::vector<Item>& items, const ChainingSpace<Item, Source>& space);

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
 * come from to reach an item. Also, once a given number of good-looking
 * predecessor items have been found, stop looking back.
 *
 * Limits transitions to those involving indels of the given size or less, to
 * avoid very bad transitions.
 */
template<typename Score, typename Item, typename Source = void, typename Collection = VectorView<Item>>
Score chain_items_dp(vector<Score>& best_chain_score,
                     const Collection& to_chain,
                     const ChainingSpace<Item, Source>& space,
                     const LookbackStrategy& lookback_strategy = ExponentialLookbackStrategy(),
                     int item_bonus = 0,
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
 * Look at the sources which occur between the given items in the read, and
 * locate their hits that occur in the part of the graph around the given items
 * in the graph, if any. Graph search is limited by max_fallow_search_distance. 
 *
 * Takes a lazy place to put an inverse of the source score
 * order sort (which is the space they are numbered in for item
 * references), from the space's sources view.
 *
 * Also takes a function which, when given a Source and the sorted IDs of a
 * subgraph, calls back with each pos_t in the subgraph that can correspond to
 * the source to form an Item.
 *
 * Returns new, forged items in ascending order by read start position. 
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
                                  size_t read_length,
                                  const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                  size_t max_fallow_search_distance = 10000);
                                  
/**
 * Look at the sources which occur between the given collections of items in
 * the read, and locate their hits that occur in the part of the graph around
 * the given items in the graph, if any. Graph search is limited by
 * max_fallow_search_distance.
 *
 * If either collection is empty, searcches between the read start/end and the
 * items in the other collection.
 *
 * Takes a lazy place to put an inverse of the source score
 * order sort (which is the space they are numbered in for item
 * references), from the space's sources view.
 *
 * Also takes a function which, when given a Source and the sorted IDs of a
 * subgraph, calls back with each pos_t in the subgraph that can correspond to
 * the source to form an Item.
 *
 * Returns new, forged items in ascending order by read start position. 
 *
 * Forged items will not have all fields set. They will only have locations and
 * sources.
 *
 * TODO: Use a different kind of item type to avoid this forgery!
 */
template<typename Item, typename Source, typename Collection=VectorView<Item>>
vector<Item> reseed_fallow_region(const Collection& left_items,
                                  const Collection& right_items,
                                  const ChainingSpace<Item, Source>& space,
                                  std::unique_ptr<VectorViewInverse>& source_sort_inverse,
                                  size_t read_length, 
                                  const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                  size_t max_fallow_search_distance = 10000);

/**
 * Reseed all fallow regions (regions in the read where the distance
 * between items, or an item and the end of the read, is at least
 * fallow_region_size) by locating hits of the given sources that occur in the
 * part of the graph between existing items.
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
                                           size_t read_length,
                                           const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                           size_t max_fallow_search_distance = 10000,
                                           size_t fallow_region_size = 200);

// --------------------------------------------------------------------------------

// Template implementations

template<typename Item, typename Source = void>
void sort_and_shadow(const std::vector<Item>& items, std::vector<size_t>& indexes, const ChainingSpace<Item, Source>& space) {
    
    // Sort the indexes by read start ascending, and read end descending
    std::sort(indexes.begin(), indexes.end(), [&](const size_t& a, const size_t& b) {
        auto& a_item = items[a];
        auto& b_item = items[b];
        auto a_start = space.read_start(a_item);
        auto b_start = space.read_start(b_item);
        // a should be first if it starts earlier, or starts atthe same place and ends later.
        return (a_start < b_start || (a_start == b_start && space.read_end(a_item) > space.read_end(b_item)));
    });
    
    // Keep a collection of the diagonal pairs that are already represented,
    // and the read end position of the latest-ending item on those pairs that
    // we have taken. A diagonal is defined as a graph node ID, a graph strand,
    // and the difference between the graph offset and the read position. So we
    // can represent them with pos_t, and subtract the read position out of the
    // stored offset to make them.
    std::unordered_map<std::pair<pos_t, pos_t>, size_t> diagonal_progress;
    
    // Scan through and make a new collection of indexes, keeping the first on
    // any pair of diagonals, which will thus be the one with the earliest
    // start, and within those the latest end. Since we need to keep items
    // which partially overlap but don't contain each other, we also keep an
    // item if it is the new latest-ending thing we've seen for a pair of
    // diagonals.
    std::vector<size_t> kept_indexes;
    kept_indexes.reserve(indexes.size());
    for (auto i : indexes) {
        // For each item we might keep
        auto& item = items[i];
        // Fetch out the read bounds so we can reuse them
        auto item_read_start = space.read_start(item);
        auto item_read_end = space.read_end(item);
        
        // Prepare the key of the diagonals it visits
        std::pair<pos_t, pos_t> diagonal = std::make_pair(space.graph_start(item), space.graph_end(item));
        // Make the offsets store a difference between graph and read offset so
        // they really represent diagonals.
        get_offset(diagonal.first) -= item_read_start;
        get_offset(diagonal.second) -= item_read_end;
        
        auto& furthest_read_end = diagonal_progress[diagonal];
        if (furthest_read_end < item_read_end) {
            // This is the first, or latest-ending, item seen on this diagonal.
            // If there was an earlier-ending item taken, we know it started before this one, because of iteration order.
            // So take this item.
            kept_indexes.push_back(i);
            // And record that we got out this far
            furthest_read_end = item_read_end;
#ifdef debug_chaining
            std::cerr << "Keep " << space.to_string(item) << " which gets us to R" << furthest_read_end << " on diagonal (" << diagonal.first << ", " << diagonal.second << ")"  << std::endl;
#endif
        } else {
#ifdef debug_chaining
            std::cerr << "Discard " << space.to_string(item) << " as shadowed because we already got to R" << furthest_read_end << " on diagonal (" << diagonal.first << ", " << diagonal.second << ")"  << std::endl;
#endif
        }
    }
    
    // Replace the indexes with the sorted and deduplicated ones.
    indexes = std::move(kept_indexes);
}

template<typename Item, typename Source = void>
void sort_and_shadow(std::vector<Item>& items, const ChainingSpace<Item, Source>& space) {
    // Use the index-based implementation and then apply those indexes
    std::vector<size_t> indexes = range_vector(items.size());
    sort_and_shadow(items, indexes, space);
    std::vector<Item> kept_items;
    kept_items.reserve(indexes.size());
    for (auto& index : indexes) {
        kept_items.emplace_back(std::move(items[index]));
    }
    items = std::move(kept_items);
}

template<typename Score, typename Item, typename Source, typename Collection>
Score chain_items_dp(vector<Score>& best_chain_score,
                     const Collection& to_chain,
                     const ChainingSpace<Item, Source>& space,
                     const LookbackStrategy& lookback_strategy,
                     int item_bonus,
                     size_t max_indel_bases) {
    
    DiagramExplainer diagram;
    diagram.add_globals({{"rankdir", "LR"}});
    
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
    
    // Set up the lookback tracking so we know when to stop looking at
    // predecessors.
    std::unique_ptr<LookbackStrategy::Problem> lookback_problem = lookback_strategy.setup_problem(to_chain.size());
    
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
        
        // How many points is it worth to collect?
        auto item_points = space.score(here) + item_bonus;
        
        std::string here_gvnode = "i" + std::to_string(i);
        
        // If we come from nowhere, we get those points.
        best_chain_score[i] = std::max(best_chain_score[i], ST::annotate(item_points, ST::nowhere()));
        
#ifdef debug_chaining
        cerr << "Look at transitions to #" << i
            << " at " << space.to_string(here);
        cerr << endl;
#endif

#ifdef debug_chaining
        cerr << "\tFirst item overlapping #" << i << " beginning at " << space.read_start(here) << " is #" << *first_overlapping_it << " past-ending at " << space.read_end(to_chain[*first_overlapping_it]) << " so start before there." << std::endl;
#endif
        
        if (space.graph && space.distance_index) {
            // Tell the lookback problem where we are in the graph.
            lookback_problem->place_in_graph(i, space.graph_start(here), *space.distance_index, *space.graph);
        }

        // Start considering predecessors for this item.
        lookback_problem->advance();
        auto predecessor_index_it = first_overlapping_it;
        while (predecessor_index_it != read_end_order.begin()) {
            --predecessor_index_it;
            // For each source that ended before here started, in reverse order by end position...
            auto& source = to_chain[*predecessor_index_it];
            
#ifdef debug_chaining
            cerr << "\tConsider transition from #" << *predecessor_index_it << ": " << space.to_string(source) << endl;
#endif

            // How far do we go in the read?
            size_t read_distance = space.get_read_distance(source, here);

            // See if this looks like a promising source.
            LookbackStrategy::verdict_t verdict = lookback_problem->should_check(*predecessor_index_it, i, read_distance, ST::score(best_chain_score[*predecessor_index_it]));
            if (verdict == LookbackStrategy::STOP) {
#ifdef debug_chaining
                cerr << "\t\tStop because the lookback strategy says so." << endl;
#endif
                break;
            } else if (verdict == LookbackStrategy::SKIP) {
#ifdef debug_chaining
                cerr << "\t\tSkip because the lookback strategy says so." << endl;
#endif
                continue;
            }
            
            // Now it's safe to make a distance query
#ifdef debug_chaining
            cerr << "\t\tCome from score " << best_chain_score[*predecessor_index_it]
                << " across " << space.to_string(source, here) << endl;
#endif
            
            // We will actually evaluate the source.
            
            // How far do we go in the graph?
            // We use a pointer as an optional.
            size_t graph_distance_storage;
            size_t* graph_distance_ptr;
            if (space.graph) {
                graph_distance_storage = space.get_graph_distance(source, here);
                graph_distance_ptr = &graph_distance_storage;
            } else {
                graph_distance_ptr = nullptr;
            }
            
            // How much does it pay (+) or cost (-) to make the jump from there
            // to here?
            // Don't allow the transition if it seems like we're going the long
            // way around an inversion and needing a huge indel.
            int jump_points = space.scorer.score_transition(read_distance, graph_distance_ptr, max_indel_bases);
            // And how much do we end up with overall coming from there.
            int achieved_score;
            
            if (jump_points != numeric_limits<int>::min()) {
                // Get the score we are coming from
                typename ST::Score source_score = ST::score_from(best_chain_score, *predecessor_index_it);
                
                // And the score with the transition and the points from the item
                typename ST::Score from_source_score = ST::add_points(source_score, jump_points + item_points);
                
                // Remember that we could make this jump
                best_chain_score[i] = std::max(best_chain_score[i],
                                               from_source_score);
                                               
#ifdef debug_chaining
                cerr << "\t\tWe can reach #" << i << " with " << source_score << " + " << jump_points << " from transition + " << item_points << " from item = " << from_source_score << endl;
#endif
                if (ST::score(from_source_score) > 0) {
                    // Only explain edges that were actual candidates since we
                    // won't let local score go negative
                    
                    std::string source_gvnode = "i" + std::to_string(*predecessor_index_it);
                    // Suggest that we have an edge, where the edges that are the best routes here are the most likely to actually show up.
                    diagram.suggest_edge(source_gvnode, here_gvnode, here_gvnode, ST::score(from_source_score), {
                        {"label", std::to_string(jump_points)},
                        {"weight", std::to_string(std::max<int>(1, ST::score(from_source_score)))}
                    });
                }
                
                achieved_score = ST::score(from_source_score);
            } else {
#ifdef debug_chaining
                cerr << "\t\tTransition is impossible." << endl;
#endif
                achieved_score = std::numeric_limits<size_t>::min();
            }
            
            // Note that we checked out this transition and saw the observed scores and distances.
            lookback_problem->did_check(*predecessor_index_it, i, read_distance, graph_distance_ptr, jump_points, achieved_score);
        }
        
#ifdef debug_chaining
        cerr << "\tBest way to reach #" << i << " is " << best_chain_score[i] << endl;
#endif
        
        std::stringstream label_stream;
        label_stream << "#" << i << " " << space.to_string(here) << " = " << item_points << "/" << ST::score(best_chain_score[i]);
        diagram.add_node(here_gvnode, {
            {"label", label_stream.str()}
        });
        if (space.graph) {
            auto graph_start = space.graph_start(here);
            std::string graph_gvnode = "n" + std::to_string(id(graph_start)) + (is_rev(graph_start) ? "r" : "f");
            diagram.ensure_node(graph_gvnode, {
                {"label", std::to_string(id(graph_start)) + (is_rev(graph_start) ? "-" : "+")},
                {"shape", "box"}
            });
            // Show the item as connected to its source graph node
            diagram.add_edge(here_gvnode, graph_gvnode, {{"color", "gray"}});
            // Make the next graph node along the same strand
            std::string graph_gvnode2 = "n" + std::to_string(id(graph_start) + (is_rev(graph_start) ? -1 : 1)) + (is_rev(graph_start) ? "r" : "f");
            diagram.ensure_node(graph_gvnode2, {
                {"label", std::to_string(id(graph_start) + (is_rev(graph_start) ? -1 : 1)) + (is_rev(graph_start) ? "-" : "+")},
                {"shape", "box"}
            });
            // And show them as connected. 
            diagram.ensure_edge(graph_gvnode, graph_gvnode2, {{"color", "gray"}});
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
        cerr << "Chain ends at #" << here << " " << space.to_string(to_chain[here])
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
                << " " << space.to_string(to_chain[here]) << endl;
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
                                  size_t read_length,
                                  const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                  size_t max_fallow_search_distance) {
                                  
    // We can just decide the pointers are really to std::array<const Item, 1>.
    const std::array<const Item, 1>* left_items((const std::array<const Item, 1>*)&left);
    const std::array<const Item, 1>* right_items((const std::array<const Item, 1>*)&right);
    
    return reseed_fallow_region<Item, Source, std::array<const Item, 1>>(*left_items, *right_items, space, source_sort_inverse, read_length, for_each_pos_for_source_in_subgraph, max_fallow_search_distance);
}

template<typename Item, typename Source, typename Collection>
vector<Item> reseed_fallow_region(const Collection& left_items,
                                  const Collection& right_items,
                                  const algorithms::ChainingSpace<Item, Source>& space,
                                  std::unique_ptr<VectorViewInverse>& source_sort_inverse,
                                  size_t read_length,
                                  const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                  size_t max_fallow_search_distance) {
    
    // Make sure we aren't being given no seeds
    assert(!left_items.empty() || !right_items.empty());
    
    // We have two seeds here, so there aren't no sources.
    // We know there must be a collection of sources in read order, so grab that.
    const vector<Source>& sources_in_read_order = *space.sources.items;
    
    // Find the range of indexes in read order for the source hits between
    // those of the bounding items.
    size_t range_begin = left_items.empty() ? 0 : std::numeric_limits<size_t>::max();
    for (auto& left : left_items) {
        range_begin = std::min(range_begin, space.sources.backing_index(left.source) + 1);
    }
    size_t range_end = right_items.empty() ? sources_in_read_order.size() : 0;
    for (auto& right : right_items) {
        range_end = std::max(range_end, space.sources.backing_index(right.source));
    }
    
    if (range_begin >= range_end) {
        // Not actually any sources between here in the read.
        return {};
    }
    
    // Determine bounds in read we care about
    size_t read_region_start = left_items.empty() ? 0 : std::numeric_limits<size_t>::max();
    for (auto& left : left_items) {
        read_region_start = std::min(read_region_start, space.read_end(left));
    }
    size_t read_region_end =  right_items.empty() ? read_length : 0;
    for (auto& right : right_items) {
        read_region_end = std::max(read_region_end, space.read_start(right));
    }
    
    // Decide how far we want to search
    size_t graph_distance = std::min(max_fallow_search_distance, (read_region_end - read_region_start) * 2);
    
    // Collect all the seed positions to extract around.
    std::vector<pos_t> seed_positions;
    seed_positions.reserve(left_items.size() + right_items.size());
    // And the distance from them in each direction to search
    std::vector<size_t> position_forward_max_dist;
    position_forward_max_dist.reserve(seed_positions.size());
    std::vector<size_t> position_backward_max_dist;
    position_backward_max_dist.reserve(seed_positions.size());
    
    for (auto& left : left_items) {
        // We look forward from all the left items
        seed_positions.push_back(space.graph_end(left));
        position_forward_max_dist.push_back(graph_distance);
        position_backward_max_dist.push_back(0);
    }
    for (auto& right : right_items) {
        // And we look back form all the right items
        seed_positions.push_back(space.graph_start(right));
        position_forward_max_dist.push_back(0);
        position_backward_max_dist.push_back(graph_distance);
    }
    
    // Extract the containing graph. We need the IDs in sorted order.
    // Because we are not cutting any nodes, node IDs don't change.
    vector<nid_t> sorted_ids;
    {
        bdsg::HashGraph subgraph;
        algorithms::extract_containing_graph(space.graph, &subgraph, seed_positions, graph_distance);
        sorted_ids.reserve(subgraph.get_node_count());
        subgraph.for_each_handle([&](const handle_t& h) {
            sorted_ids.push_back(subgraph.get_id(h));
        });
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
            cerr << "Found new seed for read-order source " << i << " of " << space.to_string(forged_items.back()) << endl;
#endif

            // Now make sure that we don't produce any redundant items, either
            // with other ones we've produced, or with boundary items.
            
            // So, drop the item if it's not strictly between the bounding
            // items in the read, given what the space sees as its bounds.
            
            if (space.read_start(forged_items.back()) < read_region_start || space.read_end(forged_items.back()) > read_region_end) {
                // We've gone outside the bounds of the region we are supposed
                // to be working on, and might be shadowed by or shadow
                // something out there.
#ifdef debug_reseeding
                cerr << "\tNew seed is out of reseeding range " << read_region_start << "-" << read_region_end << " in read; discarding." << endl;
#endif
                forged_items.pop_back();
                return;
            }
        });
    }
    
    // Now we just need to make sure that our forged items don't shadow/duplicate each other.
    sort_and_shadow(forged_items, space);

    return forged_items;
}

template<typename Item, typename Source>
pair<size_t, size_t> reseed_fallow_regions(vector<Item>& item_storage,
                                           vector<size_t>& sorted_item_indexes,
                                           const algorithms::ChainingSpace<Item, Source>& space,
                                           std::unique_ptr<VectorViewInverse>& source_sort_inverse,
                                           size_t read_length,
                                           const std::function<void(const Source&, const std::vector<nid_t>&, const std::function<void(const pos_t&)>&)>& for_each_pos_for_source_in_subgraph,
                                           size_t max_fallow_search_distance,
                                           size_t fallow_region_size) {
    
    // You might think that we could just walk through all the items in read
    // order and do the gaps between successive items that are longer than the
    // threshold.
    //
    // But you would be wrong! Some items will be tied in read order, and place
    // the same piece of the read at different places in the graph. We need to
    // consider all participants in each tie set at each end of a break in the
    // items.
    //
    // Then we need to extract around all the bounds, in the correct
    // directions, at once, and get hits in all relevant parts of the graph. 
    
    // TODO: To really make sure we are getting the correct participants on the
    // left side, we really need to look at the items sorted by read end
    // position, so we can get everything ending a the same point. Or, we need
    // to get some latest-ending item on the left side of the gap, and
    // everything overlapping it.
    // For now we just group by read start position, which could mislead us
    // depending on where the graph node breaks are and how we define our
    // items' lengths.
    
    pair<size_t, size_t> statistics {0, 0};
    auto& fallow_region_count = statistics.first;
    auto& longest_fallow_region_length = statistics.second;
    
    // Make a VectorView over the items
    VectorView<Item> item_view {item_storage, sorted_item_indexes};
    
    // Make a set of items we started with or have found so far.
    // We need to track start position and end position, because we can have
    // reverse-strand minimizers produce hits that have different end
    // positions, but get cut off at the same start position.
    std::unordered_set<std::pair<std::pair<size_t, pos_t>, std::pair<size_t, pos_t>>> seen_items;
    for (auto& item : item_view) {
        // Make sure none of the initial items are duplicates.
        size_t read_start = space.read_start(item);
        pos_t graph_start = space.graph_start(item);
        size_t read_end = space.read_end(item);
        pos_t graph_end = space.graph_end(item);
        auto key = std::make_pair(std::make_pair(read_start, graph_start), std::make_pair(read_end, graph_end));
        auto found = seen_items.find(key);
        if (found != seen_items.end()) {
            throw std::runtime_error("Duplicate initial mapping " + space.to_string(item));
        }
        seen_items.emplace_hint(found, std::move(key));
    }
    
    // Use a function to reseed across a distance between two possibly-empty runs bounding a fallow region.
    auto reseed = [&](size_t read_distance, size_t left_run_start, size_t left_run_end, size_t right_run_start, size_t right_run_end) {
        // Count it
        fallow_region_count++;
        longest_fallow_region_length = std::max(longest_fallow_region_length, read_distance);
        
        // Find all the indexes of left items
        std::vector<size_t> left_indexes;
        left_indexes.reserve(left_run_end - left_run_start);
        for (size_t i = left_run_start; i != left_run_end; i++) {
            left_indexes.push_back(item_view.backing_index(i));
        }
        
        // And all the indexes of right items
        std::vector<size_t> right_indexes;
        right_indexes.reserve(right_run_end - right_run_start);
        for (size_t i = right_run_start; i != right_run_end; i++) {
            right_indexes.push_back(item_view.backing_index(i));
        }
    
#ifdef debug_reseeding
        cerr << "Reseeding between ";
        if (!left_indexes.empty()) {
            cerr << left_indexes.size() << " seeds near #" << item_view.backing_index(left_run_start)
                 << " " << space.to_string(item_view[left_run_start]);
        } else {
            cerr << "left read end";
        }
        cerr << " and ";
        if (!right_indexes.empty()) {
            cerr << right_indexes.size() << " seeds near #" << item_view.backing_index(right_run_start)
                << " " << space.to_string(item_view[right_run_start]);
        } else {
            cerr << "right read end";
        }
        cerr << endl;
#endif

        // Forge some fill-in items
        vector<Item> new_items = reseed_fallow_region({item_storage, left_indexes},
                                                      {item_storage, right_indexes},
                                                      space,
                                                      source_sort_inverse,
                                                      read_length,
                                                      for_each_pos_for_source_in_subgraph,
                                                      max_fallow_search_distance);
                                                      
        // The newly-found items are not going to be duplicates of
        // each other, but because reseeding subgraphs can overlap,
        // we might get the same exact hits multiple times in
        // apparently different fallow regions. So we just drop
        // duplicates now.
        //
        // TODO: Make the reseeding subgraphs not overlap so we
        // stop wasting time???
       
        size_t kept = 0;
        for (auto& item : new_items) {
            // Deduplicate the newly-found items with previous reseed queries.
            // TODO: We assume that identical bounds means identical items.
            size_t read_start = space.read_start(item);
            pos_t graph_start = space.graph_start(item);
            size_t read_end = space.read_end(item);
            pos_t graph_end = space.graph_end(item);
            auto key = std::make_pair(std::make_pair(read_start, graph_start), std::make_pair(read_end, graph_end));
            auto found = seen_items.find(key);
            if (found == seen_items.end()) {
                // Append the fill-in items onto the item storage vector. Safe to do
                // even when we're using a view over it.
                item_storage.emplace_back(std::move(item));
                seen_items.emplace_hint(found, std::move(key));
                kept++;
            }
        }
#ifdef debug_reseeding
        std::cerr << "Kept " << kept << " new seeds" << std::endl;
#endif
    };
    
    // Find the first run of things starting at the same place in the read.
    size_t current_run_start = 0;
    size_t current_run_end = 1;
    while (current_run_end < item_view.size() && space.read_start(item_view[current_run_end]) == space.read_start(item_view[current_run_start])) {
        // Scan until we find the next thing that starts after us.
        ++current_run_end;
    }
    
    if (space.read_start(item_view[current_run_start]) > fallow_region_size) {
        // Start to first run is fallow.
        // Note that all items in the run start at the same base.
        
        // Find all the indexes of current run items. Treat them as "right"
        reseed(space.read_start(item_view[current_run_start]), 0, 0, current_run_start, current_run_end);
    } else {
#ifdef debug_reseeding
        std::cerr << "No need to reseed from start of read at 0"
                  << " because run near #" << item_view.backing_index(current_run_start)
                  << " " << space.to_string(item_view[current_run_start]) 
                  << " is within " << fallow_region_size << " of read start" << std::endl;
#endif
    }
    
    while (current_run_end < item_view.size()) {
        // For each run of items sharing a start position and having something after them
        
        // Find all the next run items
        size_t next_run_start = current_run_end;
        size_t next_run_end = next_run_start + 1;
        while (next_run_end < item_view.size() && space.read_start(item_view[next_run_end]) == space.read_start(item_view[next_run_start])) {
            // Scan until we find the next thing that starts after us.
            ++next_run_end;
        }
        
        // Check read distance. TODO: Could vary depending on where the head item ends. Did we guarantee a sort by both ends?
        size_t read_distance = space.get_read_distance(item_view[current_run_start], item_view[next_run_start]);
        
        if (read_distance != std::numeric_limits<size_t>::max() && read_distance > fallow_region_size) {
            
            // Now we know that the items between current_run_start and current_run_end form a fallow region with the items between next_run_start and next_run_end.
            reseed(read_distance, current_run_start, current_run_end, next_run_start, next_run_end);
        }
        
        // Make the next run into the current run
        current_run_start = next_run_start;
        current_run_end = next_run_end;
    }
    
    // Find where the last seed ends in the read
    size_t current_run_read_end = 0;
    for (size_t i = current_run_start; i < current_run_end; i++) {
        current_run_read_end = std::max(current_run_read_end, space.read_end(item_view[i]));
    }
    if (read_length - current_run_read_end > fallow_region_size) {
        // Last run to end is fallow
        reseed(read_length - current_run_read_end, current_run_start, current_run_end, 0, 0);
    } else {
#ifdef debug_reseeding
        std::cerr << "No need to reseed to end of read at " << read_length
                  << " because run near #" << item_view.backing_index(current_run_start)
                  << " " << space.to_string(item_view[current_run_start]) << " that ends at "
                  << current_run_read_end << " is within " << fallow_region_size << " of read end" << std::endl;
#endif
    }
    
    if (item_storage.size() != sorted_item_indexes.size()) {
        // After filling in all the fallow regions, if we actually got any new
        // items, extend and re-sort sorted_item_indexes.
        // TODO: It would be nice to just do this as we go along and paste
        // together the index arrays, but the overlapping-item-run code makes
        // that fiddly because we don't just walk down the input index array.
        // So we re-sort instead.
        sorted_item_indexes.reserve(item_storage.size());
        while(sorted_item_indexes.size() < item_storage.size()) {
            sorted_item_indexes.push_back(sorted_item_indexes.size());
        }
        
        // TODO: De-duplicate sort code with the initial sort
        std::sort(sorted_item_indexes.begin(), sorted_item_indexes.end(), [&](const size_t& a, const size_t& b) -> bool {
            auto a_start = space.read_start(item_storage[a]);
            auto b_start = space.read_start(item_storage[b]);
            return a_start < b_start || (a_start == b_start && space.read_end(item_storage[a]) < space.read_end(item_storage[b]));
        });
    }
    
    return statistics;
}

}
}

#endif

#ifndef VG_FUNNEL_HPP_INCLUDED
#define VG_FUNNEL_HPP_INCLUDED

#include <string>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <vg/vg.pb.h>
#include "annotation.hpp"


/** 
 * \file funnel.hpp
 * Contains the Funnel class, for recording the history of complex multi-stage transformations of sets of results.
 */
 
namespace vg {

using namespace std;

/**
 * Represents a record of an invocation of a pipeline for an input.
 *
 * Tracks the history of "lines" of data "item" provenance through a series of
 * "stages", containing a series of "filters".
 *
 * Lines are "introduced", and "project" from earlier stages to later stages,
 * possibly "expanding" or "merging", until they "fail" a filter or reach the
 * final stage. At each stage, items occur in a linear order and are identified
 * by index.
 *
 * An item may be a "group", with a certain size.
 *
 * We also can assign "scores" or correctness/placed-ness "tags" to items at a
 * stage. Tags can cover a region of a linear read space.
 */
class Funnel {

public:
    /// Start processing the given named input.
    /// Name must not be empty.
    /// No stage or substage will be active.
    void start(const string& name);
    
    /// Stop processing the given named input.
    /// All stages and substages are stopped.
    void stop();
    
    /// Start the given stage, and end all previous stages and substages.
    /// Name must not be empty.
    /// Multiple stages with the same name will be coalesced.
    void stage(const string& name);
    
    /// Stop the current stage.
    void stage_stop();
    
    /// Start the given substage, nested insude the current stage. End all previous substages.
    /// Substages within a stage may repeat and are coalesced.
    /// Name must not be empty.
    void substage(const string& name);
    
    /// Stop the current substage.
    void substage_stop();
    
    /// Start processing the given item coming from the previous stage.
    void processing_input(size_t prev_stage_item);
    
    /// Stop processing an item from the previous stage.
    void processed_input();
    
    /// Start producing the given output item, whether it has been projected yet or not. 
    void producing_output(size_t item);
    
    /// Stop producing an output item.
    void produced_output();
    
    /// Introduce the given number of new items, starting their own lines of provenance (default 1).
    void introduce(size_t count = 1);
    
    /// Expand the given item from the previous stage into the given number of new items at this stage.
    void expand(size_t prev_stage_item, size_t count);
    
    /// Merge all the given item indexes from the previous stage into a new item at this stage.
    /// The new item will be a group, sized according to the number of previous items merged.
    template<typename Iterator>
    void merge_group(Iterator prev_stage_items_begin, Iterator prev_stage_items_end);
    
    /// Merge all the given item indexes from the previous stage into a new item at this stage.
    /// The new item will be a group, sized according to the total size of
    /// previous groups, with non-groups counting as size 1.
    template<typename Iterator>
    void merge_groups(Iterator prev_stage_items_begin, Iterator prev_stage_items_end);
    
    /// Merge all the given item indexes from the previous stage into a new item at this stage.
    /// The new item will be a single item.
    template<typename Iterator>
    void merge(Iterator prev_stage_items_begin, Iterator prev_stage_items_end);
    
    /// Record extra provenance relationships where the latest current-stage
    /// item came from the given previous-stage items. Increases the
    /// current-stage item group size by the number of previous-stage items
    /// added.
    ///
    /// Propagates tagging.
    template<typename Iterator>
    void also_merge_group(Iterator prev_stage_items_begin, Iterator prev_stage_items_end);
    
    /// Record extra provenance relationships where the latest current-stage
    /// item came from the given earlier-stage items. Increases the
    /// current-stage item group size by the number of previous-stage items
    /// added.
    ///
    /// Propagates tagging.
    ///
    /// earlier_stage_lookback determines how many stages to look back and must be
    /// 1 or more.
    template<typename Iterator>
    void also_merge_group(size_t earlier_stage_lookback, Iterator earlier_stage_items_begin, Iterator earlier_stage_items_end);
    
    /// Record an extra provenance relationship where the latest current-stage
    /// item came from the given previous-stage item, the given number of
    /// stages ago (min 1).
    ///
    /// Does not adjust group size or propagate tagging.
    void also_relevant(size_t earlier_stage_lookback, size_t earlier_stage_item); 
    
    /// Project a single item from the previous stage to a single non-group item at this stage.
    void project(size_t prev_stage_item);
    
    /// Project a single item from the previous stage to a new group item at the current stage, with the given size.
    void project_group(size_t prev_stage_item, size_t group_size);
    
    /// Fail the given item from the previous stage on the given filter and do not project it through to this stage.
    /// Items which do not fail a filter must pass the filter and be projected to something.
    /// The filter name must survive the funnel, because a pointer to it will be stored.
    /// Allows a statistic for the filtered-on value for the failing item to be recorded.
    void fail(const char* filter, size_t prev_stage_item, double statistic = nan(""));
    
    /// Pass the given item from the previous stage through the given filter at this stage.
    /// Items which do not pass a filter must fail it.
    /// All items which pass filters must do so in the same order.
    /// The filter name must survive the funnel, because a pointer to it will be stored.
    /// Allows a statistic for the filtered-on value for the passing item to be recorded.
    void pass(const char* filter, size_t prev_stage_item, double statistic = nan(""));
    
    /// Assign the given score to the given item at the current stage.
    void score(size_t item, double score);
    
    /// We can tag items as having one of these states.
    enum class State {
        NONE = 0,
        PLACED = 1,
        CORRECT = 2
    };
    
    /// Tag the given item as being in the given state at the current stage.
    /// Future items that derive from it will inherit these tags. Optionally
    /// allows specifying that the state extends over a range in read space.
    void tag(size_t item, State state, size_t tag_start = 0, size_t tag_length = std::numeric_limits<size_t>::max());

    /// Tag the given item as "correct" at the current stage. Future items that
    /// derive from it will also be tagged as correct.
    /// Optionally allows specifying that the correctness extends over a range
    /// in read space, so correctness can be tracked as a property of regions
    /// of the read, rather than the whole read.
    /// If called multiple times, with different bounds, the correct region
    /// will enclose all the correct regions provided in the different calls.
    void tag_correct(size_t item, size_t tag_start = 0, size_t tag_length = std::numeric_limits<size_t>::max());
    
    /// Return true if the given item at this stage is tagged correct, or
    /// descends from an item that was tagged correct.
    bool is_correct(size_t item) const;

    /// Return true if the given item at the previous stage is tagged correct, or
    /// descends from an item that was tagged correct.
    bool was_correct(size_t prev_stage_item) const;
    
    /// Return true if the given item at the given named previous stage is
    /// tagged correct, or descends from an item that was tagged correct.
    /// Needs a hint about what number the stage was in the order, to make
    /// lookup fast.
    bool was_correct(size_t prev_stage_index, const string& prev_stage_name, size_t prev_stage_item) const;
    
    /// Get the tagged regions with the given tag (or better) for the given item in the current stage.
    /// Regions are returned as start inclusive, end exclusive.
    /// TODO: Right now there are always 0 or 1 regions.
    std::vector<std::pair<size_t, size_t>> where_is_tagged(State state, size_t item) const;

    /// Get the name of the most recent stage that had a correct-tagged item
    /// survive into it, or "none" if no items were ever tagged correct.
    /// Optionally allows specifying a read space interval to intersect with
    /// items, so the query returns the last stage that had a correct item
    /// intersecting that range.
    string last_correct_stage(size_t tag_start = 0, size_t tag_length = std::numeric_limits<size_t>::max()) const;
    
    /// Get the name of the most recent stage that had a n item tagged with the
    /// given tag or better survive into it, or "none" if no items were ever
    /// tagged that good. Optionally allows specifying a read space interval to
    /// intersect with items, so the query returns the last stage that had an
    /// item intersecting that range and also an item witht hat tag or better.
    ///
    /// TODO: Make worse tag ranges not match queries for better tags!
    string last_tagged_stage(State tag, size_t tag_start = 0, size_t tag_length = std::numeric_limits<size_t>::max()) const;
    
    /// Get the index of the most recent item created in the current stage.
    size_t latest() const;
    
    /// Call the given callback with stage name, and vector of result item
    /// sizes at that stage, and a duration in seconds, for each stage.
    void for_each_stage(const function<void(const string&, const vector<size_t>&, const double&)>& callback) const;
    
    /// Represents the performance of a filter, for either item counts or total item sizes.
    /// Note that passing_correct and failing_correct will always be 0 if nothing is tagged correct.
    struct FilterPerformance {
        size_t passing = 0;
        size_t failing = 0;
        size_t passing_correct = 0;
        size_t failing_correct = 0;
    };
    
    /// Call the given callback with stage name, filter name, performance
    /// report for items, performance report for total size of items, values
    /// for correct items for the filter statistic, and values for incorrect
    /// (or merely not known-correct) items for the filter statistic.
    /// Runs the callback for each stage and filter, in order. Only includes
    /// filters that were actually passed or failed by any items.
    void for_each_filter(const function<void(const string&, const string&,
        const FilterPerformance&, const FilterPerformance&,
        const vector<double>&, const vector<double>&)>& callback) const;

    /// Dump information from the Funnel as a dot-format Graphviz graph to the given stream.
    /// Illustrates stages and provenance.
    void to_dot(ostream& out) const;

    /// Set an alignments annotations with the number of results at each stage
    /// if annotate_correctness is true, also annotate the alignment with the
    /// number of correct results at each stage. This assumes that we've been
    /// tracking correctness all along
    void annotate_mapped_alignment(Alignment& aln, bool annotate_correctness) const;
    
protected:
    
    /// Pick a clock to use for measuring stage duration
    using clock = std::chrono::high_resolution_clock;
    /// And a type to represent stage transition times
    using time_point = clock::time_point;
    
    /// What's the name of the funnel we start()-ed. Will be empty if nothing is running.
    string funnel_name;
    
    /// At what time did we start()
    time_point start_time;
    
    /// At what time did we stop()
    time_point stop_time;
    
    /// What's the name of the current stage? Will be empty if no stage is running.
    string stage_name;
    
    /// At what time did the stage start?
    time_point stage_start_time;
    
    /// What's the name of the current substage? Will be empty if no substage is running.
    string substage_name;
    
    /// What's the current prev-stage input we are processing?
    /// Will be numeric_limits<size_t>::max() if none.
    size_t input_in_progress = numeric_limits<size_t>::max();
    
    /// what's the current current-stage output we are generating?
    /// Will be numeric_limits<size_t>::max() if none.
    size_t output_in_progress = numeric_limits<size_t>::max();
    
    // Now members we need for provenance tracking
    
    /// Represents a flag vector over positions via a sorted interval list.
    /// Allows setting flags in a range.
    struct PaintableSpace {
        /// Mark a range as painted
        void paint(size_t start, size_t length);
        /// Check if any position in the given range is painted
        bool is_any_painted(size_t start, size_t length) const;
        
        /// Store start position and length for all painted intervals.
        std::map<size_t, size_t> regions;
    };
    
    /// Represents an Item whose provenance we track
    struct Item {
        size_t group_size = 0;
        double score = 0;
        /// Is this item tagged with a state, or a descendant of a tagged item?
        State tag = State::NONE;
        /// If the item is tagged, over what interval is it tagged?
        /// When projecting, intervals are combined by min/maxing the bounds.
        size_t tag_start = std::numeric_limits<size_t>::max();
        size_t tag_length = 0;
        /// What previous stage items were combined to make this one, if any?
        vector<size_t> prev_stage_items = {};
        /// And what items from stages before that? Recorded as (stage offset,
        /// item number) pairs; all the offsets will be >=2.
        vector<pair<size_t, size_t>> earlier_stage_items = {};
        /// What filters did the item pass at this stage, if any?
        vector<const char*> passed_filters = {};
        /// And what statistics did they have (or NaN)?
        vector<double> passed_statistics = {};
        /// What filter did the item finally fail at at this stage, if any?
        const char* failed_filter = nullptr;
        /// And what statistic did it fail with (or NaN)?
        double failed_statistic = nan("");
    };
    
    /// Represents a Stage which is a series of Items, which track their own provenance.
    struct Stage {
        string name;
        vector<Item> items;
        /// How long did the stage last, in seconds?
        float duration;
        /// How many of the items were actually projected?
        /// Needed because items may need to expand to hold information for items that have not been projected yet.
        size_t projected_count = 0;
        /// What's the best tag of anything at this stage?
        State tag = State::NONE;
        /// Where are tags applied?
        PaintableSpace tag_space;
    };
    
    /// Ensure an item with the given index exists in the current stage and return a reference to it.
    /// We need to do it this way because we might save a production duration before an item is really projected.
    /// The items of the current stage should only be modified through this.
    /// Note that you do *not* need to create an item in order to get it.
    Item& get_item(size_t index);
    
    /// Create a new item in the current stage and get its index.
    /// Advances the projected count counter.
    size_t create_item();
    
    /// Rercord all the stages, including their names and item provenance.
    /// Handles repeated stages.
    vector<Stage> stages;
};

inline std::ostream& operator<<(std::ostream& out, const Funnel::State& state) {
    switch (state) {
        case Funnel::State::NONE:
            return out << "NONE";
        case Funnel::State::PLACED:
            return out << "PLACED";
        case Funnel::State::CORRECT:
            return out << "CORRECT";
        default:
            return out << "UNKNOWN";
    }
}

template<typename Iterator>
void Funnel::merge_group(Iterator prev_stage_items_begin, Iterator prev_stage_items_end) {
    // Do a non-group merge
    merge(prev_stage_items_begin, prev_stage_items_end);
    
    // Find it
    size_t index = latest();

    // Update its size
    get_item(index).group_size = get_item(index).prev_stage_items.size();
}

template<typename Iterator>
void Funnel::merge_groups(Iterator prev_stage_items_begin, Iterator prev_stage_items_end) {
    // Do a non-group merge
    merge(prev_stage_items_begin, prev_stage_items_end);
    
    // Find it
    size_t index = latest();
    
    // Compute the total size it should have
    auto& prev_stage = stages[stages.size() - 2];
    size_t total_size = 0;
    for (auto& prev : get_item(index).prev_stage_items) {
        size_t prev_group_size = prev_stage.items[prev].group_size;
        if (prev_group_size == 0) {
            // Non-groups count as size 1
            prev_group_size = 1;
        }
        total_size += prev_group_size;
    }

    // Update its size
    get_item(index).group_size = total_size;
}

template<typename Iterator>
void Funnel::merge(Iterator prev_stage_items_begin, Iterator prev_stage_items_end) {
    // There must be a prev stage to merge from
    assert(stages.size() > 1);
    auto& prev_stage = stages[stages.size() - 2];

    // Make a new item to combine all the given items.
    size_t index = create_item();

    for (Iterator& it = prev_stage_items_begin; it != prev_stage_items_end; ++it) {
        // For each prev stage item
        size_t prev_stage_item = *it;

        // Make sure it existed
        assert(prev_stage.items.size() > prev_stage_item);

        // Record the dependency
        get_item(index).prev_stage_items.push_back(prev_stage_item);
        
        // Propagate tags
        auto& old = prev_stage.items[prev_stage_item];
        if (old.tag != State::NONE) {
            // Tag the new item if it came from something tagged.
            tag(index, old.tag, old.tag_start, old.tag_length);
        }
    }
}

template<typename Iterator>
void Funnel::also_merge_group(Iterator prev_stage_items_begin, Iterator prev_stage_items_end) {
    also_merge_group(1, prev_stage_items_begin, prev_stage_items_end);
}

template<typename Iterator>
void Funnel::also_merge_group(size_t earlier_stage_lookback, Iterator earlier_stage_items_begin, Iterator earlier_stage_items_end) {
    assert(earlier_stage_lookback > 0);
    assert(stages.size() > earlier_stage_lookback);
    auto& earlier_stage = stages[stages.size() - 1 - earlier_stage_lookback];
    auto& item = get_item(latest());
    
    for (Iterator& it = earlier_stage_items_begin; it != earlier_stage_items_end; ++it) {
        // For each earlier stage item
        size_t earlier_stage_item = *it;

        // Make sure it existed
        assert(earlier_stage.items.size() > earlier_stage_item);

        // Record the dependency
        if (earlier_stage_lookback == 1) {
            // References to the immediately preceeding stage are special
            item.prev_stage_items.push_back(earlier_stage_item);
        } else {
            // References to earlier stages include the stage offset back
            item.earlier_stage_items.emplace_back(earlier_stage_lookback, earlier_stage_item);
        }
        
        // Increase group size
        item.group_size += 1;
        
        // Propagate tags
        auto& old = earlier_stage.items[earlier_stage_item];
        if (old.tag != State::NONE) {
            // Tag the new item if it came from something tagged.
            tag(latest(), old.tag, old.tag_start, old.tag_length);
        }
    }
}

}

#endif

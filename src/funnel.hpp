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
 * We also can assign scores to items at a stage.
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

    /// Tag the given item as "correct" at the current stage. Future items that
    /// derive from it will also be tagged as correct.
    void tag_correct(size_t item);
    
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

    /// Get the name of the most recent stage that had a correct-tagged item
    /// survive into it, or "none" if no items were ever tagged correct.
    string last_correct_stage() const;
    
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
    void to_dot(ostream& out);

    /// Set an alignments annotations with the number of results at each stage
    /// if annotate_correctness is true, also annotate the alignment with the
    /// number of correct results at each stage. This assumes that we've been
    /// tracking correctness all along
    void annotate_mapped_alignment(Alignment& aln, bool annotate_correctness);
    
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
    
    /// Represents an Item whose provenance we track
    struct Item {
        size_t group_size = 0;
        double score = 0;
        /// Is this item tagged as correct, or a descendant of a tagged item?
        bool correct = false;
        /// What previous stage items were combined to make this one, if any?
        vector<size_t> prev_stage_items = {};
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
        /// Does this stage contain any items tagged as correct?
        bool has_correct = false;
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

template<typename Iterator>
void Funnel::merge_group(Iterator prev_stage_items_begin, Iterator prev_stage_items_end) {
    // There must be a prev stage to merge from
    assert(stages.size() > 1);
    auto& prev_stage = stages[stages.size() - 2];

    // Make a new item to hold all the given items.
    size_t index = create_item();

    for (Iterator& i = prev_stage_items_begin; i != prev_stage_items_end; ++i) {
        // For each prev stage item (not copying the iterator)
        size_t prev_stage_item = *i;

        // Make sure it existed
        assert(prev_stage.items.size() > prev_stage_item);

        // Record the dependency
        get_item(index).prev_stage_items.push_back(prev_stage_item);
         
        if (prev_stage.items[prev_stage_item].correct) {
            // Tag the new item correct if it came from something correct.
            tag_correct(index);
        }
    }

    // Update its size
    get_item(index).group_size = get_item(index).prev_stage_items.size();
}

}

#endif

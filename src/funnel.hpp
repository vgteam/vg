#ifndef VG_FUNNEL_HPP_INCLUDED
#define VG_FUNNEL_HPP_INCLUDED

#include <string>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <functional>
#include <iostream>
#include <limits>

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
 * "stages".
 *
 * Lines are "introduced", and "project" from earlier stages to later stages,
 * possibly "expanding" or "merging", until they are "killed" or reach the
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
    
    /// Kill the given item from the previous stage and do not project it through to this stage.
    /// Items which are not killed must be projected to something.
    void kill(size_t prev_stage_item);
    
    /// Gill all the given items from the previous stage.
    template<typename Iterator>
    void kill_all(Iterator prev_stage_items_begin, Iterator prev_stage_items_end);
    
    /// Assign the given score to the given item at the current stage.
    void score(size_t item, double score);

    /// Tag the given item as "correct" at the current stage. Future items that
    /// derive from it will also be tagged as correct.
    void tag_correct(size_t item);

    /// Get the name of the most recent stage that had a correct-tagged item
    /// survive into it, or "none" if no items were ever tagged correct.
    string last_correct_stage() const;
    
    /// Get the index of the most recent item created in the current stage.
    size_t latest() const;
    
    /// Call the given callback with stage name and number of results at that stage, for each stage.
    void for_each_stage(const function<void(const string&, size_t)>& callback) const;

    /// Dump information from the Funnel as a dot-format Graphviz graph to the given stream.
    /// Illustrates stages and provenance.
    void to_dot(ostream& out);
    
protected:
    
    /// What's the name of the funnel we start()-ed. Will be empty if nothing is running.
    string funnel_name;
    
    /// What's the name of the current stage? Will be empty if no stage is running.
    string stage_name;
    
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
    };
    
    /// Represents a Stage which is a series of Items, which track their own provenance.
    struct Stage {
        string name;
        vector<Item> items;
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

template<typename Iterator>
void Funnel::kill_all(Iterator prev_stage_items_begin, Iterator prev_stage_items_end) {
    for (; prev_stage_items_begin != prev_stage_items_end; ++prev_stage_items_begin) {
        // Kill each thing between begin and end
        kill(*prev_stage_items_begin);
    }
}


}

#endif

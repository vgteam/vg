#ifndef VG_FUNNEL_HPP_INCLUDED
#define VG_FUNNEL_HPP_INCLUDED

#include <chrono>
#include <string>
#include <vector>
#include <unordered_map>

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
 * timed "stages".
 *
 * Lines are "introduced", and "project" from earlier stages to later stages,
 * possibly "expanding" or "merging", until they are "killed" or reach the
 * final stage. At each stage, items occur in a linear order and are identified
 * by index.
 *
 * Each stage has its runtime recorded, and the runtime for substages of the
 * stage are recorded as well.
 *
 * An item may be a "group", with a certain size.
 *
 * We can also allocate time to processing input items from the previous stage,
 * or ourput items of the current stage.
 *
 * We also can assign scores to items at a stage.
 */
class Funnel {

public:
    /// Start the timer for processing the given named input.
    /// Name must not be empty.
    /// No stage or substage will be active.
    void start(const string& name);
    
    /// Stop the timer for processing the given named input.
    /// All stages and substages are stopped.
    void stop();
    
    /// Start the given stage, and end all previous stages and substages.
    /// Name must not be empty.
    /// Multiple stages with the same name will be coalesced.
    void stage(const string& name);
    
    /// Stop counting time against the current stage.
    void stage_stop();
    
    /// Start the given substage, nested insude the current stage. End all previous substages.
    /// Substages within a stage may repeat and are coalesced.
    /// Name must not be empty.
    void substage(const string& name);
    
    /// Stop counting time against the current substage.
    void substage_stop();
    
    /// Start a clock assigning runtime to processing the given item coming from the previous stage.
    void processing_input(size_t prev_stage_item);
    
    /// Stop the clock assigning runtime to processing an item from the previous stage.
    void processed_input();
    
    /// Start a clock assigning time to producing the given output item, whether it has been projected yet or not. 
    void producing_output(size_t item);
    
    /// Stop the clock assigning time to producing an output item.
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
    
    /// Get the index of the most recent item created in the current stage.
    size_t latest() const;
    
    /// Get the total number of seconds elapsed between start() and stop()
    double total_seconds() const;
    
protected:
    
    /// How do we record durations internally?
    using Duration = chrono::nanoseconds;
    /// How do we record timepoints internally?
    using Timepoint = chrono::time_point<chrono::high_resolution_clock>;
    
    // Members we need for time tracking
    
    /// What's the name of the funnel we start()-ed. Will be empty if nothing is running.
    string funnel_name;
    /// When did we start it the funnel?
    Timepoint funnel_start_time;
    /// If we finished the funnel, how long did it last?
    Duration funnel_duration;
    
    /// What's the name of the current stage? Will be empty if no stage is running.
    string stage_name;
    /// When did we start the stage?
    Timepoint stage_start_time;
    /// Records total duration of all stages by name.
    unordered_map<string, Duration> stage_durations;
    
    /// What's the name of the current substage? Will be empty if no substage is running.
    string substage_name;
    /// When did we start the substage?
    Timepoint substage_start_time;
    /// Records total duration of all substages by substage name and stage name.
    unordered_map<string, unordered_map<string, Duration>> substage_durations;
    
    /// What's the current prev-stage input we are processing?
    /// Will be numeric_limits<size_t>::max() if none.
    size_t input_in_progress = numeric_limits<size_t>::max();
    /// When did we start processing that input?
    Timepoint input_start_time;
    
    /// what's the current current-stage output we are generating?
    /// Will be numeric_limits<size_t>::max() if none.
    size_t output_in_progress = numeric_limits<size_t>::max();
    /// When did we start processing that input?
    Timepoint output_start_time;
    
    /// Produce the current time
    Timepoint now() const;
    
    // Now members we need for provenance tracking
    
    /// Represents an Item whose provenance we track
    struct Item {
        size_t group_size = 0;
        double score = 0;
        /// What previous stage items were combined to make this one, if any?
        vector<size_t> prev_stage_items = {};
        /// How long did it take to produce this item, in total?
        Duration produce_duration = 0;
        /// How long was spent processing this item as an input, in total?
        Duration process_duration = 0;
    };
    
    /// Represents a Stage which is a series of Items, which track their own provenance.
    struct Stage {
        string name;
        vector<Item> items;
        /// How many of the items were actually projected?
        /// Needed because items may need to expand to hold production times for items that have not been projected yet.
        size_t projected_count = 0;
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
    assert(!stages.empty());
    // Make a new item holding all the given items.
    size_t index = create_item();
    std::copy(prev_stage_items_begin, prev_stage_items_end, std::back_inserter(get_item(index).prev_stage_items));
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

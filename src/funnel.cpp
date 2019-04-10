#include "funnel.hpp"

#include <cassert>

/**
 * \file funnel.hpp: implementation of the Funnel class
 */
 
namespace vg {
using namespace std;

void Funnel::start(const string& name) {
    assert(!name.empty());

    // (Re)start the funnel.
    funnel_start_time = now();
    funnel_name = name;
    
    // Clear out old data
    stage_name.clear();
    substage_name.clear();
    stage_durations.clear();
    substage_durations.clear();
    stages.clear();
}

void Funnel::stop() {
    // Stop any lingering stage (which takes care of substages and processing/producing timers)
    stage_stop();
    
    // Record the ultimate duration.
    funnel_duration = now() - funnel_start_time;
}

void Funnel::stage(const string& name) {
    assert(!funnel_name.empty());
    assert(!name.empty());
    
    // Stop the previous stage if any.
    stage_stop();

    // Allocate new stage structures.
    stages.emplace_back();
    stages.back().name = name;
    
    // Save the name and start time
    stage_name = name;
    stage_start_time = now();
}

void Funnel::stage_stop() {
    if (!stage_name.empty()) {
        // A stage was running.
        
        // Stop any substage.
        substage_stop();
        // Stop any process/produce timers
        processed_input();
        produced_output();
        
        // Get the elapsed stage time
        Duration stage_duration = now() - stage_start_time;
        // Save it, coalescing by name
        stage_durations[stage_name] += stage_duration;
        
        // Say the stage is stopped 
        stage_name.clear();
    }
}

void Funnel::substage(const string& name) {
    assert(!funnel_name.empty());
    assert(!stage_name.empty());
    assert(!name.empty());

    // Stop previous substage if any.
    substage_stop();
    
    // Substages don't bound produce/process timers.
    
    // Save the name and start time
    substage_name = name;
    substage_start_time = now();
}
    
void Funnel::substage_stop() {
    if (!substage_name.empty()) {
        // A substage was running.
        
        // Substages don't bound produce/process timers.
        
        // Get the elapsed substage time
        Duration substage_duration = now() - substage_start_time;
        // Save it, coalescing by name
        substage_durations[stage_name][substage_name] += substage_duration;
        
        // Say the stage is stopped 
        substage_name.clear();
    }
}

void Funnel::processing_input(size_t prev_stage_item) {
    // We can only take input from previous stages, in a stage
    assert(!stage_name.empty());
    assert(stages.size() > 1);
    assert(prev_stage_item != numeric_limits<size_t>::max());
    assert(stages[stages.size() - 2].items.size() > prev_stage_item);

    // Stop any previous input processing timer
    processed_input();
    
    // Start this one
    input_in_progress = prev_stage_item;
    input_start_time = now();
}

void Funnel::processed_input() {
    if (input_in_progress != numeric_limits<size_t>::max()) {
        // We were processing an input
        
        // Work out how long it took
        Duration input_duration = now() - input_start_time;
        // Add it in
        stages[stages.size() - 2].items[input_in_progress].process_duration += input_duration;
        
        // Say we're done with the input.
        input_in_progress = numeric_limits<size_t>::max();
    }
}

void Funnel::producing_output(size_t item) {
    // We can only produce output in a stage
    assert(!stage_name.empty());
    assert(!stages.empty());
    assert(item != numeric_limits<size_t>::max());
    
    // Stop any previous input processing timer
    produced_output();
    
    // Start this one
    output_in_progress = item;
    output_start_time = now();
}

void Funnel::produced_output() {
    if (output_in_progress != numeric_limits<size_t>::max()) {
        // We were producing an output
        
        // Work out how long it took
        Duration output_duration = now() - output_start_time;
        // Add it in
        get_item(output_in_progress).produce_duration += output_duration;
        
        // Say we're done with the output.
        output_in_progress = numeric_limits<size_t>::max();
    }
}

void Funnel::introduce(size_t count) {
    // Create that many new items
    for (size_t i = 0; i < count; i++) {
        create_item();
    }
}

void Funnel::expand(size_t prev_stage_item, size_t count) {
    for (size_t i = 0; i < count; i++) {
        // Create the requested number of items
        project(prev_stage_item);
    }
}

void Funnel::project(size_t prev_stage_item) {
    // There must be a prev stage to project from
    assert(stages.size() > 1);
    auto& prev_stage = stages[stages.size() - 2];

    // Make one new item
    size_t index = create_item();

    // Record the ancestry
    get_item(index).prev_stage_items.push_back(prev_stage_item);

    if (prev_stage.items[prev_stage_item].correct) {
        // Tag the new item correct if it came from something correct
        tag_correct(index);
    }
}

void Funnel::project_group(size_t prev_stage_item, size_t group_size) {
    // Project the item
    project(prev_stage_item);
    // Save the group size
    get_item(latest()).group_size = group_size;
}

void Funnel::kill(size_t prev_stage_item) {
    // TODO: kill is a no-op for now.
    // We just don't project from it.
}

void Funnel::score(size_t item, double score) {
    get_item(item).score = score;
}

void Funnel::tag_correct(size_t item) {
    // Say the item is correct
    get_item(item).correct = true;
    // Say the stage has something correct.
    stages.back().has_correct = true;
}

string Funnel::last_correct_stage() const {
    // Just do a linear scan backward through stages
    for (auto it = stages.rbegin(); it != stages.rend(); ++it) {
        if (it->has_correct) {
            return it->name;
        }
    }
    return "none";
}

size_t Funnel::latest() const {
    assert(!stages.empty());
    assert(!stages.back().items.empty());
    return stages.back().items.size() - 1;
}

double Funnel::total_seconds() const {
    return chrono::duration_cast<chrono::duration<double>>(funnel_duration).count();
}

void Funnel::for_each_time(const function<void(const string&, const string&, double)>& callback) {
    // Handle overall
    callback("", "", total_seconds());

    for (auto& kv : stage_durations) {
        // Handle each stage
        callback(kv.first, "", chrono::duration_cast<chrono::duration<double>>(kv.second).count());
    }

    for (auto& kv : substage_durations) {
        for (auto& kv2 : kv.second) {
            // Handle each substage
            callback(kv.first, kv2.first, chrono::duration_cast<chrono::duration<double>>(kv2.second).count());
        }
    }
}

Funnel::Timepoint Funnel::now() const {
    return chrono::high_resolution_clock::now();
}

Funnel::Item& Funnel::get_item(size_t index) {
    assert(!stages.empty());
    if (index >= stages.back().items.size()) {
        // Allocate up through here
        stages.back().items.resize(index + 1);
    }
    return stages.back().items[index];
}

size_t Funnel::create_item() {
    assert(!stages.empty());
    
    // Work out where to put it
    size_t next_index = stages.back().projected_count;
    // Make sure the item slot exists
    get_item(next_index);
    // Record the item's creation
    stages.back().projected_count++;
    
    // Return the index used
    return next_index;
}
    



}














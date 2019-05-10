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

double Funnel::to_seconds(const Duration& time) const {
    return chrono::duration_cast<chrono::duration<double>>(time).count();
}

double Funnel::total_seconds() const {
    return to_seconds(funnel_duration);
}

void Funnel::for_each_time(const function<void(const string&, const string&, double)>& callback) const {
    // Handle overall
    callback("", "", total_seconds());

    for (auto& kv : stage_durations) {
        // Handle each stage
        callback(kv.first, "", to_seconds(kv.second));
    }

    for (auto& kv : substage_durations) {
        for (auto& kv2 : kv.second) {
            // Handle each substage
            callback(kv.first, kv2.first, to_seconds(kv2.second));
        }
    }
}

void Funnel::for_each_stage(const function<void(const string&, size_t)>& callback) const {
    for (auto& stage : stages) {
        // Report the name and item count of each stage.
        callback(stage.name, stage.items.size());
    }
}

void Funnel::to_dot(ostream& out) {
    out << "digraph graphname {" << endl;
    out << "rankdir=\"TB\";" << endl;

    // Get total time
    double total_time = total_seconds();

    for (size_t s = 0; s < stages.size(); s++) {
        // For each stage in order
        auto& stage = stages[s];

        // Compute a GraphViz ID part for the stage
        string stage_id = "s" + to_string(s);

        // Start a subgraph.
        // Prepend cluster so it draws as a box.
        out << "subgraph cluster_" << stage_id << " {" << endl;
        out << "label = \"" << stage.name << "\";" << endl;
        out << "graph[style=solid];" << endl;
        out << "rank=same;" << endl;

        for (size_t i = 0; i < stage.items.size(); i++) {
            // For each item in the stage
            auto& item = stage.items[i];

            // Compute a GraphViz ID
            string item_id = stage_id + "i" + to_string(i);

            // Emit a node
            out << item_id << "[label=\"" << i << "\" shape=circle tooltip=\"";
            if (item.group_size != 0) {
                out << "size " << item.group_size;
            }
            if (item.score != 0) {
                if (item.group_size != 0) {
                    out << ", ";
                }
                out << "score " << item.score;
            }
            out << "\"";
            if (item.correct) {
                // Make it green if it is correct
                out << " color=green";
            }
            out << "];" << endl;

            if (s > 0) {
                // There is a previous stage, so we can draw edges from it.
                for (auto& p : item.prev_stage_items) {
                    // Connect everything from the previous stage to it
                    auto& prev_item = stages[s - 1].items.at(p);

                    out << "s" << (s - 1) << "i" << p << " -> " << item_id << "[";
                    if (item.correct && prev_item.correct) {
                        // Correctness came this way
                        out << "color=green";
                    }
                    out << "];" << endl;
                }
            }

            // Assign generation time
            double item_time = to_seconds(item.produce_duration);
            if (item_time > 0 && total_time > 0) {
                // We can have a time node attached to this item
                double item_portion = item_time / total_time;

                // Make an ID for the time node
                string time_id = item_id + "t1";
                
                // Make the time node like a little bar chart
                out << time_id << "[shape=box style=filled label=\"\" color=red fixedsize=true width=\"0.5\" height=\""
                    << item_portion << "\" tooltip=\"" << item_time << " seconds" << "\"];" << endl;

                // Attach it to the item
                out << item_id << "->" << time_id << "[dir=none style=dashed constraint=false];" << endl;
            }

        }

        if (s > 1) {
            // Make time nodes in this stage for processing things from the previous stage
            for (size_t prev_stage_item = 0; prev_stage_item < stages[s - 1].items.size(); prev_stage_item++) {
                auto& item = stages[s - 1].items[prev_stage_item];

                double item_time = to_seconds(item.process_duration);
                if (item_time > 0 && total_time > 0) {
                    // We can have a time node attached to this item
                    double item_portion = item_time / total_time;

                    // Make an ID for the time node
                    string item_id = "s" + to_string(s - 1) + "i" + to_string(prev_stage_item);
                    string time_id = item_id + "t2";
                    
                    // Make the time node like a little bar chart
                    out << time_id << "[shape=box style=filled label=\"\" color=red fixedsize=true width=\"0.5\" height=\""
                        << item_portion << "\" tooltip=\"" << item_time << " seconds" << "\"];" << endl;

                    // Attach it to the item
                    out << item_id << "->" << time_id << "[dir=none style=dashed];" << endl;
                }
            }
        }

        out << "}" << endl;
    }

    out << "}" << endl;
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














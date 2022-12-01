#include "funnel.hpp"

#include <cassert>
#include <cstring>

/**
 * \file funnel.hpp: implementation of the Funnel class
 */
 
namespace vg {
using namespace std;

void Funnel::PaintableSpace::paint(size_t start, size_t length) {
    // Find the last interval starting strictly before start
    auto predecessor = regions.lower_bound(start);
    if (predecessor != regions.begin()) {
        --predecessor;
        // We have one.
        
        if (predecessor->first + predecessor->second >= start) {
            // If its length is long enough to abut or cover start
            
            if (predecessor->first + predecessor->second > start + length) {
                // It completely encompasses us, so nothing to do!
                return;
            }
            
            // Budge start back and increase length
            length += (start - predecessor->first);
            start = predecessor->first;
            
            // And remove it
            regions.erase(predecessor);
            // TODO: Can we fix it up?
        }
    }
    
    // Find the first interval starting at or after start
    auto successor = regions.upper_bound(start);
    auto range_first = regions.end();
    auto range_last = regions.end();
    while (successor != regions.end() && successor->first <= start + length) {
        // For each from there that starts at or before start + length
        // Increase length to cover up to its end
        length = std::max(successor->first + successor->second, start + length) - start;
        // And remember to remove it it
        if (range_first == regions.end()) {
            range_first = successor;
        }
        // Check the next thing
        ++successor;
        // Which provides the removal past-end
        range_last = successor;
    }
    
    // Remove the covered intervals
    regions.erase(range_first, range_last);
    
    // Add the new interval
    regions.emplace(start, length);
}

bool Funnel::PaintableSpace::is_any_painted(size_t start, size_t length) const {
    // Find the last interval starting strictly before start
    auto predecessor = regions.lower_bound(start);
    if (predecessor != regions.begin()) {
        --predecessor;
        // We have one.
        if (predecessor->first + predecessor->second > start) {
            // It covers our start, so we overlap
            return true;
        }
    }
    
    auto successor = regions.upper_bound(start);
    if (successor != regions.end()) {
        // There's something starting at or after us
        if (start + length > successor->first) {
            // And we overlap it
            return true;
        }
    }
    
    // We can't overlap anything
    return false;
}

void Funnel::start(const string& name) {
    assert(!name.empty());

    // (Re)start the funnel.
    funnel_name = name;
    start_time = clock::now();
    
    // Clear out old data
    stage_name.clear();
    substage_name.clear();
    stages.clear();
}

void Funnel::stop() {
    // Stop any lingering stage (which takes care of substages)
    stage_stop();
    // Stop the funnel overall
    stop_time = clock::now();
}

void Funnel::stage(const string& name) {
    assert(!funnel_name.empty());
    assert(!name.empty());
    
    // Stop the previous stage if any.
    stage_stop();

    // Allocate new stage structures.
    stages.emplace_back();
    stages.back().name = name;
    
    // Save the name
    stage_name = name;
    
    // Record the start time
    stage_start_time = clock::now();
}

void Funnel::stage_stop() {
    if (!stage_name.empty()) {
        // A stage was running.
        
        // Stop any substage.
        substage_stop();
        // Stop any process/produce 
        processed_input();
        produced_output();
        
        // Say the stage is stopped 
        stage_name.clear();
        
        // Record the duration in seconds
        auto stage_stop_time = clock::now();
        stages.back().duration = chrono::duration_cast<chrono::duration<double>>(stage_stop_time - stage_start_time).count();
    }
}

void Funnel::substage(const string& name) {
    assert(!funnel_name.empty());
    assert(!stage_name.empty());
    assert(!name.empty());

    // Stop previous substage if any.
    substage_stop();
    
    // Substages don't bound produce/process.
    
    // Save the name 
    substage_name = name;
}
    
void Funnel::substage_stop() {
    if (!substage_name.empty()) {
        // A substage was running.
        
        // Substages don't bound produce/process.
        
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

    // Stop any previous input processing
    processed_input();
    
    // Start this one
    input_in_progress = prev_stage_item;
}

void Funnel::processed_input() {
    if (input_in_progress != numeric_limits<size_t>::max()) {
        // We were processing an input
        
        // Say we're done with the input.
        input_in_progress = numeric_limits<size_t>::max();
    }
}

void Funnel::producing_output(size_t item) {
    // We can only produce output in a stage
    assert(!stage_name.empty());
    assert(!stages.empty());
    assert(item != numeric_limits<size_t>::max());
    
    // Stop any previous input processing
    produced_output();
    
    // Start this one
    output_in_progress = item;
}

void Funnel::produced_output() {
    if (output_in_progress != numeric_limits<size_t>::max()) {
        // We were producing an output
        
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
    
    auto& old = prev_stage.items[prev_stage_item];

    if (old.tag != State::NONE) {
        // Tag the new item if it came from something tagged
        tag(index, old.tag, old.tag_start, old.tag_length);
    }
}

void Funnel::project_group(size_t prev_stage_item, size_t group_size) {
    // Project the item
    project(prev_stage_item);
    // Save the group size
    get_item(latest()).group_size = group_size;
}

void Funnel::also_relevant(size_t earlier_stage_lookback, size_t earlier_stage_item) {
    assert(earlier_stage_lookback > 0);
    assert(stages.size() > earlier_stage_lookback);
    auto& earlier_stage = stages[stages.size() - 1 - earlier_stage_lookback];
    assert(earlier_stage.items.size() > earlier_stage_item);
    auto& item = get_item(latest());
    if (earlier_stage_lookback == 1) {
        // References to the immediately preceeding stage are special
        item.prev_stage_items.push_back(earlier_stage_item);
    } else {
        // References to earlier stages include the stage offset back
        item.earlier_stage_items.emplace_back(earlier_stage_lookback, earlier_stage_item);
    }
}

void Funnel::fail(const char* filter, size_t prev_stage_item, double statistic) {
    // There must be a prev stage to project from
    assert(stages.size() > 1);
    auto& prev_stage = stages[stages.size() - 2];

    // Record the item as having failed this filter
    prev_stage.items[prev_stage_item].failed_filter = filter;
    prev_stage.items[prev_stage_item].failed_statistic = statistic;
}

void Funnel::pass(const char* filter, size_t prev_stage_item, double statistic) {
    // There must be a prev stage to project from
    assert(stages.size() > 1);
    auto& prev_stage = stages[stages.size() - 2];

    // Record the item as having passed this filter
    prev_stage.items[prev_stage_item].passed_filters.emplace_back(filter);
    prev_stage.items[prev_stage_item].passed_statistics.emplace_back(statistic);
}

void Funnel::score(size_t item, double score) {
    get_item(item).score = score;
}

void Funnel::tag(size_t item, State state, size_t tag_start, size_t tag_length) {

#ifdef debug
    std::cerr << "Tag item " << item << " stage " << stages.back().name << " as " << state << " on " << tag_start << "-" << tag_start + tag_length << std::endl;
#endif
    
    // Say the item is tagged
    auto& to_mark = get_item(item);
    to_mark.tag = std::max(to_mark.tag, state);
    
    if (to_mark.tag_start == std::numeric_limits<size_t>::max() && to_mark.tag_length == 0) {
        // Item hasn't been tagged before, so we can just adopt the passed range.
        to_mark.tag_start = tag_start;
        to_mark.tag_length = tag_length;
    } else {
        // We need to find the enclosing range of the existing range and the new range.
        size_t correct_end = std::max(to_mark.tag_start + to_mark.tag_length, tag_start + tag_length);
        to_mark.tag_start = std::min(to_mark.tag_start, tag_start);
        to_mark.tag_length = correct_end - to_mark.tag_start;
    }
    
#ifdef debug
    std::cerr << "\tNow tagged over " << to_mark.tag_start << "-" << to_mark.tag_start + to_mark.tag_length << std::endl;
#endif
    
    // TODO: Allow different tags to cover different ranges?
    // TODO: Allow per-item gapped range tracking?
    
    // Say the stage has tag over this interval.
    stages.back().tag = std::max(stages.back().tag, state);
    stages.back().tag_space.paint(tag_start, tag_length);
}

void Funnel::tag_correct(size_t item, size_t tag_start, size_t tag_length) {
    tag(item, State::CORRECT, tag_start, tag_length);
}

bool Funnel::is_correct(size_t item) const {
    return stages.back().items[item].tag >= State::CORRECT;
}

bool Funnel::was_correct(size_t prev_stage_item) const {
    assert(stages.size() > 1);
    auto& prev_stage = stages[stages.size() - 2];
    return prev_stage.items[prev_stage_item].tag >= State::CORRECT;
}

bool Funnel::was_correct(size_t prev_stage_index, const string& prev_stage_name, size_t prev_stage_item) const {
    assert(stages.size() > prev_stage_index);
    auto& prev_stage = stages[prev_stage_index];
    assert(prev_stage.name == prev_stage_name);
    return prev_stage.items[prev_stage_item].tag >= State::CORRECT;
}

string Funnel::last_tagged_stage(State tag, size_t tag_start, size_t tag_length) const {
    // Just do a linear scan backward through stages
    for (auto it = stages.rbegin(); it != stages.rend(); ++it) {
        if (it->tag >= tag && it->tag_space.is_any_painted(tag_start, tag_length)) {
            // If we are tagged good enough and have a tag in part of that
            // area, then we are a matching stage.
            return it->name;
        }
    }
    return "none";
}

string Funnel::last_correct_stage(size_t tag_start, size_t tag_length) const {
    return last_tagged_stage(State::CORRECT, tag_start, tag_length); 
}

size_t Funnel::latest() const {
    assert(!stages.empty());
    assert(!stages.back().items.empty());
    return stages.back().items.size() - 1;
}

void Funnel::for_each_stage(const function<void(const string&, const vector<size_t>&, const double&)>& callback) const {
    for (auto& stage : stages) {
        // Make a vector of item sizes
        vector<size_t> item_sizes;
        item_sizes.reserve(stage.items.size());
        for (auto& item : stage.items) {
            item_sizes.push_back(item.group_size);
        }
        // Report the name and item count of each stage.
        callback(stage.name, item_sizes, stage.duration);
    }
}

void Funnel::for_each_filter(const function<void(const string&, const string&,
    const FilterPerformance&, const FilterPerformance&, const vector<double>&, const vector<double>&)>& callback) const {
    
    for (auto& stage : stages) {
        // Hold the names of all filters encountered
        vector<const char*> filter_names;
        // And the by-item and by-size performance stats.
        vector<pair<FilterPerformance, FilterPerformance>> filter_performances;
        // And the correct and not-known-correct filter statistic values
        vector<pair<vector<double>, vector<double>>> filter_statistics;
        
        for (auto& item : stage.items) {
            // For each item
            size_t filter_index;
            for (filter_index = 0; filter_index < item.passed_filters.size(); filter_index++) {
                // For each filter it passed
                if (filter_index >= filter_names.size()) {
                    // If it is new
                
                    // Remember its name in the list of filters
                    filter_names.push_back(item.passed_filters[filter_index]);
                    
                    // And give it an empty report
                    filter_performances.emplace_back();
                    filter_statistics.emplace_back();
                } else {
                    // Make sure the name is correct
                    // TODO: can we justy match on pointer value and not string value?
                    assert(strcmp(filter_names[filter_index], item.passed_filters[filter_index]) == 0);
                }
                
                // Record passing
                filter_performances[filter_index].first.passing++;
                filter_performances[filter_index].first.passing_correct += item.tag >= State::CORRECT;
                
                filter_performances[filter_index].second.passing += item.group_size;
                filter_performances[filter_index].second.passing_correct += item.tag >= State::CORRECT ? item.group_size : 0;
                
                if (item.tag >= State::CORRECT) {
                    // Record this statistic value as belonging to a correct item
                    filter_statistics[filter_index].first.push_back(item.passed_statistics[filter_index]);
                } else {
                    // Record this statistic value as belonging to a not necessarily correct item
                    filter_statistics[filter_index].second.push_back(item.passed_statistics[filter_index]);
                }
            }
            
            if (item.failed_filter != nullptr) {
                // For the final, failed filter, if any
                
                if (filter_index >= filter_names.size()) {
                    // If it is new
                    
                    // Remember its name in the list of filters
                    filter_names.push_back(item.failed_filter);
                    
                    // And give it an empty report
                    filter_performances.emplace_back();
                    filter_statistics.emplace_back();
                } else {
                    // Make sure the name is correct
                    // TODO: can we justy match on pointer value and not string value?
                    assert(strcmp(filter_names[filter_index], item.failed_filter) == 0);
                }
                
                // Record failing
                filter_performances[filter_index].first.failing++;
                filter_performances[filter_index].first.failing_correct += (item.tag >= State::CORRECT) ? 1 : 0;
                
                filter_performances[filter_index].second.failing += item.group_size;
                filter_performances[filter_index].second.failing_correct += (item.tag >= State::CORRECT) ? item.group_size : 0;
                
                if (item.tag >= State::CORRECT) {
                    // Record this statistic value as belonging to a correct item
                    filter_statistics[filter_index].first.push_back(item.failed_statistic);
                } else {
                    // Record this statistic value as belonging to a not necessarily correct item
                    filter_statistics[filter_index].second.push_back(item.failed_statistic);
                }
            }
        }
        
        // Now we have gone through the filters for this stage for every item.
        
        for (size_t i = 0; i < filter_names.size(); i++) {
            // For each filter
            
            // Report the results tabulated across items.
            callback(stage.name, filter_names[i],
                filter_performances[i].first, filter_performances[i].second,
                filter_statistics[i].first, filter_statistics[i].second);
        }
    }
}

void Funnel::to_dot(ostream& out) const {
    out << "digraph graphname {" << endl;
    out << "rankdir=\"TB\";" << endl;

    for (size_t s = 0; s < stages.size(); s++) {
        // For each stage in order
        auto& stage = stages[s];

        // Compute a GraphViz ID part for the stage
        string stage_id = "s" + to_string(s);

        // Start a subgraph.
        // Prepend cluster so it draws as a box.
        out << "subgraph cluster_" << stage_id << " {" << endl;
        out << "label = \"" << stage.name;
        if (stage.tag != State::NONE) {
            // Put in if it is tagged and over what area
            out << ", " << stage.tag;
            size_t range_min = std::numeric_limits<size_t>::max();
            size_t range_max = 0;
            for (auto& region : stage.tag_space.regions) {
                if (region.first != 0 || region.second != std::numeric_limits<size_t>::max()) {
                    // Expand bounds by the bounds of this region, to summarize
                    range_min = std::min(range_min, region.first);
                    range_max = std::max(range_max, region.first + region.second);
                }
            }
            if (range_min != 0 && range_max != std::numeric_limits<size_t>::max()) {
                out << " " << range_min << "-" << range_max << " in " << stage.tag_space.regions.size() << " regions";
            }
        }
        out << "\";" << endl;
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
            if (item.tag != State::NONE && (item.tag_start != std::numeric_limits<size_t>::max() || item.tag_length != 0)) {
                if (item.group_size != 0 || item.score != 0) {
                    out << ", ";
                }
                out << "tagged " << item.tag_start << " - " << item.tag_start + item.tag_length; 
            }
            out << "\"";
            if (item.tag >= State::CORRECT) {
                // Make it green if it is correct
                out << " color=green";
            } else if (item.tag >= State::PLACED) {
                out << " color=blue";
            }
            out << "];" << endl;

            if (s > 0) {
                // There is a previous stage, so we can draw edges from it.
                for (auto& p : item.prev_stage_items) {
                    // Connect everything from the previous stage to it
                    auto& prev_item = stages[s - 1].items.at(p);

                    out << "s" << (s - 1) << "i" << p << " -> " << item_id << "[";
                    if (item.tag >= State::CORRECT && prev_item.tag >= State::CORRECT) {
                        // Correctness came this way
                        out << "color=green";
                    } else if (item.tag >= State::PLACED && prev_item.tag >= State::PLACED) {
                        // Placedness came this way
                        out << "color=blue";
                    }
                    out << "];" << endl;
                }
                if (s > 1) {
                    // And there are other stages before that
                    for (auto& p : item.earlier_stage_items) {
                        // Connect everything from the earlier stages to it
                        assert(p.first > 1);
                        auto& prev_item = stages.at(s - p.first).items.at(p.second);

                        out << "s" << (s - p.first) << "i" << p.second << " -> " << item_id << "[constraint=false";
                        if (item.tag >= State::CORRECT && prev_item.tag >= State::CORRECT) {
                            // Correctness came this way
                            out << ",color=green";
                        } else if (item.tag >= State::PLACED && prev_item.tag >= State::PLACED) {
                            // Placedness came this way
                            out << ",color=blue";
                        }
                        out << "];" << endl;
                    }
                }
            }

        }

        out << "}" << endl;
    }

    out << "}" << endl;
}
void Funnel::annotate_mapped_alignment(Alignment& aln, bool annotate_correctness) const {
    // Save the total duration in the field set asside for it
    aln.set_time_used(chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count());
    
    for_each_stage([&](const string& stage, const vector<size_t>& result_sizes, const double& duration) {
        // Save the number of items
        set_annotation(aln, "stage_" + stage + "_results", (double)result_sizes.size());
        // And the per-stage duration
        set_annotation(aln, "stage_" + stage + "_time", duration);
    });
    
    set_annotation(aln, "last_placed_stage", last_tagged_stage(State::PLACED));
    for (size_t i = 0; i < aln.sequence().size(); i += 500) {
        // For each 500 bp window, annotate with the last stage that had something placed in or spanning the window.
        // TODO: This is terrible, use an array or something.
        set_annotation(aln, "last_placed_stage_" + std::to_string(i) + "bp", last_tagged_stage(State::PLACED, i, 500));
    }
    
    if (annotate_correctness) {
        // And with the last stage at which we had any descendants of the correct seed hit locations
        set_annotation(aln, "last_correct_stage", last_correct_stage());
    }
    
    // Annotate with the performances of all the filters
    // We need to track filter number
    size_t filter_num = 0;
    for_each_filter([&](const string& stage, const string& filter,
        const Funnel::FilterPerformance& by_count, const Funnel::FilterPerformance& by_size,
        const vector<double>& filter_statistics_correct, const vector<double>& filter_statistics_non_correct) {

            string filter_id = to_string(filter_num) + "_" + filter + "_" + stage;

            // Save the stats
            set_annotation(aln, "filter_" + filter_id + "_passed_count_total", (double) by_count.passing);
            set_annotation(aln, "filter_" + filter_id + "_failed_count_total", (double) by_count.failing);
            set_annotation(aln, "filter_" + filter_id + "_passed_size_total", (double) by_size.passing);
            set_annotation(aln, "filter_" + filter_id + "_failed_size_total", (double) by_size.failing);
            
            if (annotate_correctness) {
                set_annotation(aln, "filter_" + filter_id + "_passed_count_correct", (double) by_count.passing_correct);
                set_annotation(aln, "filter_" + filter_id + "_failed_count_correct", (double) by_count.failing_correct);
                set_annotation(aln, "filter_" + filter_id + "_passed_size_correct", (double) by_size.passing_correct);
                set_annotation(aln, "filter_" + filter_id + "_failed_size_correct", (double) by_size.failing_correct);
            }
            
            // Save the correct and non-correct filter statistics, even if
            // everything is non-correct because correctness isn't computed
            set_annotation(aln, "filterstats_" + filter_id + "_correct", filter_statistics_correct);
            set_annotation(aln, "filterstats_" + filter_id + "_noncorrect", filter_statistics_non_correct);
            filter_num++;
        });
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














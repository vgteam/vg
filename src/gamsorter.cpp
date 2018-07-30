#include "gamsorter.hpp"
#include "utility.hpp"
#include "json2pb.h"
#include "position.hpp"
#include "gam_index.hpp"

/**
 * \file gamsorter.cpp
 * GAMSorter: sort a gam by position and offset.
 * Store unmapped reads at node 0.
 */

using namespace std;
using namespace vg;

GAMSorter::GAMSorter(bool show_progress) {
    this->show_progress = show_progress;
}

void GAMSorter::sort(vector<Alignment>& alns) const {
    std::sort(alns.begin(), alns.end(), [&](const Alignment& a, const Alignment& b) {
        return this->less_than(a, b);
    });
}

void GAMSorter::dumb_sort(istream& gam_in, ostream& gam_out, GAMIndex* index_to) {
    std::vector<Alignment> sort_buffer;

    stream::for_each<Alignment>(gam_in, [&](Alignment &aln) {
        sort_buffer.push_back(aln);
    });

    this->sort(sort_buffer);
    
    // Write the output in non-enormous chunks, so indexing is actually useful
    vector<Alignment> out_buffer;
    
    // Make an output emitter
    stream::ProtobufEmitter<Alignment> emitter(gam_out);
    
    if (index_to != nullptr) {
        emitter.on_group([&index_to](const vector<Alignment>& group, int64_t start_vo, int64_t past_end_vo) {
            // Whenever it emits a group, index it.
            // Make sure to only capture things that will outlive the emitter
            index_to->add_group(group, start_vo, past_end_vo);
        });
    }
    
    for (auto& aln : sort_buffer) {
        // Feed in all the sorted alignments
        emitter.write(std::move(aln));
    }
    
    // Emitter destruction will terminate the file with an EOF marker
}

void GAMSorter::stream_sort(istream& gam_in, ostream& gam_out, GAMIndex* index_to) {

    // Read the input into buffers.
    std::vector<Alignment> input_buffer;
    
    // When a buffer is full, sort it and write it to a temporary file, which
    // we remember.    
    vector<string> temp_file_names;
    
    auto finish_buffer = [&]() {
        // Do a sort. TODO: do it on another thread.
        this->sort(input_buffer);
        
        // Save it to a temp file.
        string temp_name = temp_file::create();
        temp_file_names.push_back(temp_name);
        ofstream temp_stream(temp_name);
        // OK to save as one massive group here.
        // TODO: This write could also be in a thread.
        stream::write_buffered(temp_stream, input_buffer, 0);
        
        input_buffer.clear();
    };
    
    
    stream::for_each<Alignment>(gam_in, [&](const Alignment& aln) {
        // Buffer each input alignment
        input_buffer.push_back(aln);
        if (input_buffer.size() == max_buf_size) {
            // We have a full temp file's worth of data.
            finish_buffer();
        }
    });
    finish_buffer();
    
    
    // Put all the files in a priority queue based on which has an alignment that comes first.
    // We work with pointers to cursors because we don't want to be copying the actual cursors around the heap.
    // We also *reverse* the order, because priority queues put the "greatest" element forts
    using cursor_t = stream::ProtobufIterator<Alignment>;
    auto cursor_order = [&](cursor_t*& a, cursor_t*& b) {
        if (b->has_next()) {
            if(!a->has_next()) {
                // Cursors that aren't empty come first
                return true;
            }
            return less_than(*(*b), *(*a));
        }
        return false;
    };
    priority_queue<cursor_t*, vector<cursor_t*>, decltype(cursor_order)> temp_files(cursor_order);
    
    // The open files also need to live in a collection; the cursors don't own them
    // They also can't be allowed to move since we reference them
    list<ifstream> temp_ifstreams;
    // The cursors also need to live in a collection, because we don't want to be
    // moving/copying them and their internal buffers and streams.
    // And they can't move after creation either.
    list<cursor_t> temp_cursors;
    
    for (auto& name : temp_file_names) {
        // Open each file again
        temp_ifstreams.emplace_back();
        temp_ifstreams.back().open(name);
        // Make a cursor for it and put it in the heap
        temp_cursors.emplace_back(temp_ifstreams.back());
        
        // Put the cursor pointer in the queue
        temp_files.push(&temp_cursors.back());
    }
    
    // Make an output emitter
    stream::ProtobufEmitter<Alignment> emitter(gam_out);
    
    if (index_to != nullptr) {
        emitter.on_group([&index_to](const vector<Alignment>& group, int64_t start_vo, int64_t past_end_vo) {
            // Whenever it emits a group, index it.
            // Make sure to only capture things that will outlive the emitter
            index_to->add_group(group, start_vo, past_end_vo);
        });
    }
    
    while(!temp_files.empty() && temp_files.top()->has_next()) {
        // Until we have run out of data in all the temp files
        
        // Pop off the winning cursor
        cursor_t* winner = temp_files.top();
        temp_files.pop();
        
        // Grab and emit its alignment, and advance it
        emitter.write(std::move(winner->take()));
        
        // Put it back in the heap if it is not depleted
        if (winner->has_next()) {
            temp_files.push(winner);
        }
        // TODO: Maybe keep it off the heap for the next loop somehow if it still wins
    }
    
    // The output file will be flushed and finished automatically when the emitter goes away.
    
    // The temp files will get cleaned up automatically when the program ends.
}

void GAMSorter::benedict_sort(istream& gam_in, ostream& gam_out, GAMIndex* index_to) {
    if (gam_in.tellg() == -1) {
        cerr << "error:[vg gamsort]: Cannot sort an unseekable GAM" << endl;
        exit(1);
    }
    
    // Go to the end of the file
    gam_in.seekg(0, gam_in.end);
    // Get its position
    auto file_end = gam_in.tellg();
    // Go to the start
    gam_in.seekg(0);
    // Make a progress bar
    create_progress("load positions", file_end);
    
    // This will have all the item VOs and let us sort them by position
    vector<pair<pos_t, int64_t>> pos_to_vo;
    
    stream::ProtobufIterator<Alignment> cursor(gam_in);
    
    // Count reads seen so we can only update our progress bar sometimes
    size_t seen = 0;
    
    while(cursor.has_next()) {
        // Get the min position of each alignment
        pos_t min_pos = make_pos_t(get_min_position(*cursor));
        
        // Save it with the alignment-s virtual offset
        pos_to_vo.emplace_back(min_pos, cursor.tell_item());
        
        cursor.get_next();
        
        if (seen % 1000 == 0) {
            update_progress(gam_in.tellg());
        }
        seen++;
    }
    
    update_progress(gam_in.tellg());
    destroy_progress();
    create_progress("sort positions", 1);
    
    // Sort everything by pos_t key
    std::sort(pos_to_vo.begin(), pos_to_vo.end(), [&](const pair<pos_t, int64_t>& a, const pair<pos_t, int64_t>& b) {
        return this->less_than(a.first, b.first);
    });
    
    update_progress(1);
    destroy_progress();
    create_progress("reorder reads", pos_to_vo.size());
    
    // Make an output emitter
    stream::ProtobufEmitter<Alignment> emitter(gam_out);
    
    if (index_to != nullptr) {
        emitter.on_group([&index_to](const vector<Alignment>& group, int64_t start_vo, int64_t past_end_vo) {
            // Whenever it emits a group, index it.
            // Make sure to only capture things that will outlive the emitter
            index_to->add_group(group, start_vo, past_end_vo);
        });
    }
    
    // Actually do the shuffle
    for (auto& pos_and_vo : pos_to_vo) {
        // For each item in sorted order
        
        // Load it
        cursor.seek_item_and_stop(pos_and_vo.second);
        
        // Send it out
        emitter.write(std::move(cursor.take()));
        
        increment_progress();
    }
    
    destroy_progress();
}


bool GAMSorter::less_than(const Alignment &a, const Alignment &b) const {
    return less_than(get_min_position(a), get_min_position(b));
}

Position GAMSorter::get_min_position(const Alignment& aln) const {
    return get_min_position(aln.path());
}

Position GAMSorter::get_min_position(const Path& path) const {
    if (path.mapping_size() == 0) {
        // This path lives at a default Position
        return Position();
    }
    
    Position min = path.mapping(0).position();
    for(size_t i = 1; i < path.mapping_size(); i++) {
        const Position& other = path.mapping(i).position();
        if (less_than(other, min)) {
            // We found a smaller position
            min = other;
        }
    }
    
    return min;
}

bool GAMSorter::equal_to(const Position& a, const Position& b) const {
    return (a.node_id() == b.node_id() &&
            a.is_reverse() == b.is_reverse() &&
            a.offset() == b.offset());
}

bool GAMSorter::less_than(const Position& a, const Position& b) const {
    if (a.node_id() < b.node_id()) {
        return true;
    } else if (a.node_id() > b.node_id()) {
        return false;
    }
    
    if (a.is_reverse() < b.is_reverse()) {
        return true;
    } else if (a.is_reverse() > b.is_reverse()) {
        return false;
    }

    if (a.offset() < b.offset()) {
        return true;
    }
    
    return false;
}

bool GAMSorter::less_than(const pos_t& a, const pos_t& b) const {
    if (id(a) < id(b)) {
        return true;
    } else if (id(a) > id(b)) {
        return false;
    }
    
    if (is_rev(a) < is_rev(b)) {
        return true;
    } else if (is_rev(a) > is_rev(b)) {
        return false;
    }

    if (offset(a) < offset(b)) {
        return true;
    }
    
    return false;
}

bool GAMSorter::greater_than(const Position& a, const Position& b) const {
    if (a.node_id() > b.node_id()) {
        return true;
    } else if (a.node_id() < b.node_id()) {
        return false;
    }
    
    if (a.is_reverse() > b.is_reverse()) {
        return true;
    } else if (a.is_reverse() < b.is_reverse()) {
        return false;
    }

    if (a.offset() > b.offset()) {
        return true;
    }
    
    return false;
}

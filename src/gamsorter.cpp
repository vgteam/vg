#include "gamsorter.hpp"
#include "utility.hpp"
#include "json2pb.h"
#include "position.hpp"
#include "gam_index.hpp"

#include <sys/time.h>
#include <sys/resource.h>

/**
 * \file gamsorter.cpp
 * GAMSorter: sort a gam by position and offset.
 * Store unmapped reads at node 0.
 */

using namespace std;
using namespace vg;

GAMSorter::GAMSorter(bool show_progress) {
    this->show_progress = show_progress;
    
    // We would like this many FDs max, if not limited below that.
    max_fan_in = 2048;
    // We need at least this many to sort practically.
    int min_fan_in = 100;
    
    // We need this many extra FDs not used for fan-in
    int extra_fds = 10;
    
    // Work out how many FDs we are allowed
    struct rlimit fd_limit;
    if (getrlimit(RLIMIT_NOFILE, &fd_limit) != 0) {
        // We don't know; choose a conservative default.
        max_fan_in = min_fan_in;
        cerr << "warning:[vg gamsort]: Cannot determine file descriptor limits; using " << max_fan_in << " temp file fan-in" << endl;
    } else {
        // We read the limit
        if (fd_limit.rlim_cur != RLIM_INFINITY && fd_limit.rlim_cur < max_fan_in + extra_fds) {
            // Max out our FD limit
            fd_limit.rlim_cur = min<size_t>(max_fan_in + extra_fds, fd_limit.rlim_max);
            
            if (setrlimit(RLIMIT_NOFILE, &fd_limit) != 0) {
                // We asked for a value in bound sso we should have succeeded
                throw runtime_error("Error adjusting file descriptor limit to " + to_string(fd_limit.rlim_cur)
                    + " / " + to_string(fd_limit.rlim_max));
            }
        }
        
        if (fd_limit.rlim_cur != RLIM_INFINITY && fd_limit.rlim_cur < max_fan_in + extra_fds) {
            // We need to limit ourselves to under the max FD limit
            if (fd_limit.rlim_cur < extra_fds + min_fan_in) {
                // If we can't at least do a fan-in of 10 we have a big problem.
                cerr << "error:[vg gamsort]: Open file limit very low (" << fd_limit.rlim_cur << "); we need " << (extra_fds + min_fan_in) << endl;
                exit(1);
            }
            
            // Set the max fan in to be subject to the limit
            max_fan_in = min((size_t)(fd_limit.rlim_cur - extra_fds), max_fan_in);
        }
    }
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

    // We want to work out the file size, if we can.
    size_t file_size = 0;
    {
        // Save our position
        auto here = gam_in.tellg();
        // Go to the end
        gam_in.seekg(0, gam_in.end);
        // Get its position
        auto there = gam_in.tellg();
        // Go back to where we were
        gam_in.seekg(here);
            
        if (gam_in.good()) {
            // We can seek in this stream. So how far until the end?
            file_size = there - here;
        } else {
            // It's entirely possible that none of that worked. So clear the error flags and leave the size at 0.
            gam_in.clear();
        }
    }
    
    
    // Don't give an actual 0 to the progress code or it will NaN
    create_progress("break into sorted chunks", file_size == 0 ? 1 : file_size);

    // Eventually we put sorted chunks of data in temp files and put their names here
    vector<string> outstanding_temp_files;
    
    // This tracks the number of reads in each file, by file name
    unordered_map<string, size_t> reads_per_file;
    // This tracks the total reads observed on input
    size_t total_reads_read = 0;
    
    // This cursor will read in the input file.
    cursor_t input_cursor(gam_in);
    
    #pragma omp parallel shared(gam_in, input_cursor, outstanding_temp_files, reads_per_file, total_reads_read)
    {
    
        while(true) {
    
            vector<Alignment> thread_buffer;
        
            #pragma omp critical (input_cursor)
            {
                // Each thread fights for the file and the winner reads some data
                size_t buffered_message_bytes = 0;
                while (input_cursor.has_next() && buffered_message_bytes < max_buf_size) {
                    // Until we run out of input alignments or space, buffer each, recording its size.
                    buffered_message_bytes += input_cursor.get_item_size();
                    thread_buffer.emplace_back(std::move(input_cursor.take()));
                }
            
                // Update the progress bar
                update_progress(gam_in.tellg());
            }
            
            if (thread_buffer.empty()) {
                // No data was found
                break;
            }
            
            // Do a sort of the data we grabbed
            this->sort(thread_buffer);
            
            // Save it to a temp file.
            string temp_name = temp_file::create();
            ofstream temp_stream(temp_name);
            // OK to save as one massive group here.
            // TODO: This write could also be in a thread.
            stream::write_buffered(temp_stream, thread_buffer, 0);
            
            #pragma omp critical (outstanding_temp_files)
            {
                // Remember the temp file name
                outstanding_temp_files.push_back(temp_name);
                // Remember the reads in the file, for progress purposes
                reads_per_file[temp_name] = thread_buffer.size();
                // Remember how many reads we found in the total
                total_reads_read += thread_buffer.size();
            }
        }
    }
    
    // Now we know the reader threads have taken care of the input, and all the data is in temp files.
    
    destroy_progress();
    
    while (outstanding_temp_files.size() > max_fan_in) {
        // We can't merge them all at once, so merge subsets of them.
        outstanding_temp_files = streaming_merge(outstanding_temp_files, &reads_per_file);
    }
    
    // Now we can merge (and maybe index) the final layer of the tree.
    
    // Open up cursors into all the files.
    list<ifstream> temp_ifstreams;
    list<cursor_t> temp_cursors;
    open_all(outstanding_temp_files, temp_ifstreams, temp_cursors);
    
    // Make an output emitter
    emitter_t emitter(gam_out);
    
    if (index_to != nullptr) {
        emitter.on_group([&index_to](const vector<Alignment>& group, int64_t start_vo, int64_t past_end_vo) {
            // Whenever it emits a group, index it.
            // Make sure to only capture things that will outlive the emitter
            index_to->add_group(group, start_vo, past_end_vo);
        });
    }
    
    // Merge the cursors into the emitter
    streaming_merge(temp_cursors, emitter, total_reads_read);
    
    // Clean up
    temp_cursors.clear();
    temp_ifstreams.clear();
    for (auto& filename : outstanding_temp_files) {
        temp_file::remove(filename);
    }
            
}

void GAMSorter::open_all(const vector<string>& filenames, list<ifstream>& streams, list<cursor_t>& cursors) {
    // The open files need to live in a collection; the cursors don't own them.
    // They also can't be allowed to move since we reference them.
    // The cursors also need to live in a collection, because we don't want to be
    // moving/copying them and their internal buffers and streams.
    // And they can't move after creation either.
    
    // So everything lives in caller-passed lists.
    
    for (auto& filename : filenames) {
        // Open each file
        streams.emplace_back();
        streams.back().open(filename);
        // Make a cursor for it
        cursors.emplace_back(streams.back());
    }

}

void GAMSorter::streaming_merge(list<cursor_t>& cursors, emitter_t& emitter, size_t expected_reads) {

    create_progress("merge " + to_string(cursors.size()) + " files", expected_reads == 0 ? 1 : expected_reads);
    // Count the reads we actually see
    size_t observed_reads = 0;

    // Put all the files in a priority queue based on which has an alignment that comes first.
    // We work with pointers to cursors because we don't want to be copying the actual cursors around the heap.
    // We also *reverse* the order, because priority queues put the "greatest" element forts
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
    priority_queue<cursor_t*, vector<cursor_t*>, decltype(cursor_order)> cursor_queue(cursor_order);

    for (auto& cursor : cursors) {
        // Put the cursor pointers in the queue
        cursor_queue.push(&cursor);
    }
    
    while(!cursor_queue.empty() && cursor_queue.top()->has_next()) {
        // Until we have run out of data in all the temp files
        
        // Pop off the winning cursor
        cursor_t* winner = cursor_queue.top();
        cursor_queue.pop();
        
        // Grab and emit its alignment, and advance it
        emitter.write(std::move(winner->take()));
        
        // Put it back in the heap if it is not depleted
        if (winner->has_next()) {
            cursor_queue.push(winner);
        }
        // TODO: Maybe keep it off the heap for the next loop somehow if it still wins
        
        observed_reads++;
        if (expected_reads != 0) {
            update_progress(observed_reads);
        }
    }
    
    // We finished the files, so say we're done.
    // TODO: Should we warn/fail if we expected the wrong number of reads?
    update_progress(expected_reads == 0 ? 1 : expected_reads);
    destroy_progress();

}

vector<string> GAMSorter::streaming_merge(const vector<string>& temp_files_in, unordered_map<string, size_t>* reads_per_file) {
    
    // What are the names of the merged files we create?
    vector<string> temp_files_out;
    
    // We don't do this loop in parallel because the point of looping is to limit the total currently open files.
    for (size_t start_file = 0; start_file < temp_files_in.size(); start_file += max_fan_in) {
        // For each range of sufficiently few files, starting at start_file and running for file_count
        size_t file_count = min(max_fan_in, temp_files_in.size() - start_file);
    
        // Open up cursors into all the files.
        list<ifstream> temp_ifstreams;
        list<cursor_t> temp_cursors;
        open_all(vector<string>(&temp_files_out[start_file], &temp_files_out[start_file + file_count]), temp_ifstreams, temp_cursors);
        
        // Work out how many reads to expect
        size_t expected_reads = 0;
        if (reads_per_file != nullptr) {
            for (size_t i = start_file; i < start_file + file_count; i++) {
                expected_reads += reads_per_file->at(temp_files_in.at(i));
            }
        }
        
        // Open an output file
        string out_file_name = temp_file::create();
        ofstream out_stream(out_file_name);
        temp_files_out.push_back(out_file_name);
        
        // Make an output emitter
        emitter_t emitter(out_stream);
        
        // Merge the cursors into the emitter
        streaming_merge(temp_cursors, emitter, expected_reads);
        
        // The output file will be flushed and finished automatically when the emitter goes away.
        
        // Clean up the input files we used
        temp_cursors.clear();
        temp_ifstreams.clear();
        for (size_t i = start_file; i < file_count; i++) {
            temp_file::remove(temp_files_in.at(i));
        }
        
        if (reads_per_file != nullptr) {
            // Save the total reads that should be in the created file, in case we need to do another pass
            (*reads_per_file)[out_file_name] = expected_reads;
        }
    }
    
    return temp_files_out;
        
}

void GAMSorter::benedict_sort(istream& gam_in, ostream& gam_out, GAMIndex* index_to) {
    // Go to the end of the file
    gam_in.seekg(0, gam_in.end);
    // Get its position
    auto file_end = gam_in.tellg();
    // Go to the start
    gam_in.seekg(0);
    
    // This will have all the item VOs and let us sort them by position
    vector<pair<pos_t, int64_t>> pos_to_vo;
    
    stream::ProtobufIterator<Alignment> cursor(gam_in);
    
    if (cursor.tell_raw() == -1) {
        // This will catch non-blocked gzip files, as well as streaming streams.
        cerr << "error:[vg gamsort]: Cannot sort an unseekable GAM" << endl;
        exit(1);
    }
    
    // Make a progress bar
    create_progress("load positions", file_end);
    
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

/**
 * \file
 * Minimizer (sub)sampling algorithm implementation.
 */
 
#include "sample_minimal.hpp"

#include <deque>
#include <utility>
#include <iostream>

namespace vg {
namespace algorithms {

using namespace std;

//#define debug

void sample_minimal(size_t count, size_t element_length, size_t window_size, size_t sequence_length, const std::function<size_t(size_t)>& get_start, const std::function<bool(size_t, size_t)>& should_beat, const std::function<void(size_t)>& sample) {

    if (count == 0) {
        return;
    }

    // We're going to try and do the Jain et al. 2020 algorithm as a sweep line
    // algorithm. Just in case the elements aren't dense.
    // TODO: In long-read Giraffe right now the elements are dense.

    // This will hold all the elements in the sliding window of bases, except
    // that we will drop elements that are superseded by more minimal ones.
    std::deque<size_t> queue;
    // This will hold the start of the element at the front of the queue, if any
    size_t front_start;

    // This will hold the next element not in the queue.
    size_t next_element = 0;
    // This will hold the start of the next element not in the queue yet, if any.
    size_t next_start = get_start(next_element);
#ifdef debug
    std::cerr << "Element " << next_element << " starts at " << next_start << std::endl;
#endif

    // Fill the queue for the first window
    while (next_element < count && next_start + element_length <= window_size) {
#ifdef debug
        std::cerr << "Element " << next_element << " at " << next_start << " is in first window" << std::endl;
#endif
        while (!queue.empty() && should_beat(next_element, queue.back())) {
#ifdef debug
            std::cerr << "Element " << next_element << " beats element " << queue.back() << std::endl;
#endif
            queue.pop_back();
        }
        queue.push_back(next_element);
        if (queue.front() == next_element) {
            front_start = next_start;
        }
        next_element++;
        if (next_element < count) {
            next_start = get_start(next_element);
#ifdef debug
            std::cerr << "Element " << next_element << " starts at " << next_start << std::endl;
#endif
        }
    }
    if (!queue.empty()) {
        // Find the winner fo the first window
#ifdef debug
        std::cerr << "Element " << queue.front() << " is minimal in first window" << std::endl;
#endif
        sample(queue.front());
    } else {
#ifdef debug
        std::cerr << "First window is empty" << std::endl;
#endif
    }


    // This will hold our sweep-line cursor, and is the start of the last window fully entered.
    size_t cursor = 0;
    // The first thing in the queue is also already sampled.

    while (cursor + window_size < sequence_length) {
        // More windows to consider

        // Jump to the last window if nothing intervenes
        size_t sweep_to = sequence_length - window_size;
        if (next_element < count) {
            // Or to the first window the next element is in, if closer.
            sweep_to = std::min(sweep_to, next_start + element_length - window_size); 
        }
        if (!queue.empty()) {
            // Or to the first window that the first element in the queue is not in, if closer.
            sweep_to = std::min(sweep_to, front_start + 1);
        }

#ifdef debug
        std::cerr << "Sweep to window " << sweep_to << "-" << sweep_to + window_size << std::endl;
#endif

        while (!queue.empty() && sweep_to > front_start) {
            // We are going to the first window that this element is not in.
            // Drop elements from the front of the queue that were already sampled.
#ifdef debug
            std::cerr << "Going to leave element " << queue.front() << " which started at " << front_start << std::endl;
#endif
            queue.pop_front();
            if (!queue.empty()) {
                front_start = get_start(queue.front());
                if (sweep_to > front_start) {
                    // Must be another element at the same position (as we never go past the old front_start + 1)
                    // This is a tie (since it didn't beat out the one we just popped).
                    // So sample this too.
#ifdef debug
                    std::cerr << "Element " << queue.front() << " was also minimal in window " << cursor << "-" << cursor + window_size << std::endl;
#endif
                    sample(queue.front());
                }
            }
        }

        while (next_element < count && sweep_to >= next_start + element_length - window_size) {
            // We are going to the first window that the next element is in.
#ifdef debug
            std::cerr << "Element " << next_element << " at " << next_start << " is going to be visible in window " << sweep_to << "-" << sweep_to + window_size << std::endl;
#endif
            while (!queue.empty() && should_beat(next_element, queue.back())) {
#ifdef debug
                std::cerr << "Element " << next_element << " beats element " << queue.back() << " which will never be sampled" << std::endl;
#endif
                queue.pop_back();
            }
            queue.push_back(next_element);
            if (queue.front() == next_element) {
                front_start = next_start;
            }
            next_element++;
            if (next_element < count) {
                next_start = get_start(next_element);
#ifdef debug
                std::cerr << "Element " << next_element << " starts at " << next_start << std::endl;
#endif
            }
        }
        
        if (!queue.empty()) {
            // Sample the front element because either it is now minimal
            // because we removed something in the way, or it is now minimal
            // because we added it.
#ifdef debug
            std::cerr << "Element " << queue.front() << " is minimal in new window " << sweep_to << "-" << sweep_to + window_size << std::endl;
#endif
            sample(queue.front());
        }

        // Advance the sweep line since we have fully processed the next interesting window
        cursor = sweep_to;
    }

    // Now handle ties at/exiting of the last window
    if (!queue.empty()) {
        // We consider everything that started at the same place as the front element we already sampled.
        size_t tie_front_start = front_start;
#ifdef debug
        std::cerr << "Finishing last window " << cursor << "-" << cursor + window_size << std::endl;
#endif
        while (!queue.empty() && front_start == tie_front_start) {
            // Drop elements from the front of the queue that were already sampled.
#ifdef debug
            std::cerr << "Going to leave element " << queue.front() << " which started at " << front_start << std::endl;
#endif
            queue.pop_front();
            if (!queue.empty()) {
                front_start = get_start(queue.front());
                if (front_start == tie_front_start) {
                    // Another element at the same position.
                    // This is a tie (since it didn't beat out the one we just popped).
                    // So sample this too.
#ifdef debug
                    std::cerr << "Element " << queue.front() << " was also minimal in window " << cursor << "-" << cursor + window_size << std::endl;
#endif
                    sample(queue.front());
                }
            }
        }
    }
}

}
}

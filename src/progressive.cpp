#include "progressive.hpp"

#include <vg/io/stream.hpp>
#include <iostream>
#include <chrono>

namespace vg {

using namespace std;

void Progressive::with_progress(bool show_progress, const std::string& task, const std::function<void(const std::function<void(size_t, size_t)>& progress)>& callback) {
    if (!show_progress) {
        // Use the handy no-op function from libvgio.
        callback(vg::io::NO_PROGRESS);
    } else {
        // We really do need to show progress.
        Progressive progressive;
        progressive.show_progress = show_progress;
        progressive.preload_progress(task);
        bool first_progress_update = true;

        callback([&](size_t completed, size_t total) {
            if (completed != std::numeric_limits<size_t>::max() && total != std::numeric_limits<size_t>::max()) {
                // This is a real update;
                if (first_progress_update) {
                    // This is the first update. Make the bar.
                    progressive.create_progress(total);
                    first_progress_update = false;
                }
                // Tell the bar how big to be.
                progressive.update_progress(completed);
            }
        });

        progressive.destroy_progress();
    }
}


void Progressive::create_progress(const string& message, long count) {
    if (show_progress) {
        progress_message = message;
        create_progress(count);
    }
}

void Progressive::create_progress(long count) {
    if (show_progress) {
        progress_count = count;
        last_progress = 0;
        progress_seen = 0;
        if (progress) {
            // Get rid of the old one.
            delete progress;
        }
        progress_message.resize(30, ' ');
        progress = new ProgressBar(progress_count, progress_message.c_str());
        progress->Progressed(0);
    }
}

void Progressive::preload_progress(const string& message) {
    if (show_progress && !progress) {
        // Set the message. We can't change it if a progress bar exists
        // currently because that progress bar has a pointer to the old
        // message's data.
        progress_message = message;
    }
}

void Progressive::update_progress(long i) {
    if (show_progress && progress) {
        if ((i <= progress_count
             && (long double) (i - last_progress) / (long double) progress_count >= 0.001)
            || i == progress_count) {
#pragma omp critical (progress)
            {
                progress->Progressed(i);
                last_progress = i;
            }
        }
#pragma omp critical (progress_seen)
        {
            // Remember the amout of progress we have seen.
            progress_seen = i;
        }
    }
}

void Progressive::increment_progress() {
#pragma omp critical (increment_progress)
    {
        // Only one increment can happen at a time, so we don't get lower values
        // coming after higher values because the increments and the displays
        // interleaved.

        // What should we set the progress bar to?
        size_t new_progress;

#pragma omp critical (progress_seen)
        {   
            // We need to synchronize with update calls.
            
            // Bump up whatever we have to 1 more than that.
            new_progress = progress_seen + 1;
        }
        
        // TODO: interleaving updates and increments from different threads will
        // be weird, and progress may go backward. But that makes no real sense
        // anyway, so we just have to not explode/deadlock when it happens.
    
        update_progress(new_progress);
    }
}

void Progressive::destroy_progress(void) {
    if (show_progress && progress) {
        update_progress(progress_count);
        cerr << endl;
        progress_message = "progress";
        progress_count = 0;
    }
    if (progress) {
        delete progress;
        progress = nullptr;
    }
}

}

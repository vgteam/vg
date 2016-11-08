#include "progressive.hpp"

#include <iostream>

namespace vg {

using namespace std;

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
    }
}

void Progressive::increment_progress() {
#pragma omp critical (progress)
    {
        update_progress(last_progress + 1);
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

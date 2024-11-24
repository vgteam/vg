#ifndef VG_PROGRESSIVE_HPP_INCLUDED
#define VG_PROGRESSIVE_HPP_INCLUDED

// progressive.hpp: defines a Progressive mixin class that gives any object a
// progress bar that can be turned on and off.

#include <string>
#include <functional>

#include "progress_bar.hpp"

namespace vg {

using namespace std;

/**
 * Inherit form this class to give your class create_progress(),
 * update_progress(), and destroy_progress() methods, and a public show_progress
 * field that can be toggled on and off.
 *
 * Must not be destroyed while a progress bar is active.
 */
class Progressive {

public:

    /**
     * Static callback-based progress system for places where we can't inherit from the class.
     *
     * Calls the callback with a progress function that either updates a
     * progress bar on a reasonable schedule or doesn't, depending on
     * show_progress.
     */
    static void with_progress(bool show_progress, const std::string& task, const std::function<void(const std::function<void(size_t, size_t)>& progress)>& callback);

    // Should progress bars be shown when the progress methods are called?
    bool show_progress = false;
    
    /**
     * If no progress bar is currently displayed, set the message to use for
     * the next progress bar to be created. Does nothing if show_progress is
     * false or when a progress bar is displayed.
     *
     * Public so that users of a class can provide descriptive messages for
     * generic progress operations (like VG's for_each_kmer_parallel).
     */
    void preload_progress(const string& message);
    /**
     * Create a progress bar showing the given message, with the given number of
     * items to process. Does nothing if show_progress is false. Replaces any
     * existing progress bar.
     */
    void create_progress(const string& message, long count);
    /**
     * Create a progress bar with the given number of items to process, using
     * either a default message, or the message passed to the last
     * preload_progress call since a progress bar was destroyed. Does nothing if
     * show_progress is false. Replaces any existing progress bar.
     */
    void create_progress(long count);
    /**
     * Update the progress bar, noting that the given number of items have been
     * processed. Does nothing if no progress bar is displayed.
     */
    void update_progress(long i);
    /**
     * Update the progress bar, noting that one additional item has been
     * processed. Does nothing if no progress bar is displayed.
     */
    void increment_progress();
    /**
     * Destroy the current progress bar, if it exists.
     */
    void destroy_progress(void);
    
private:
    string progress_message = "progress";
    // How many total ticks of progress are there?
    long progress_count;
    // What's the last progress value we displayed?
    long last_progress;
    // What's the last progress value we've actually seen, either through an
    // explicit update or an increment?
    long progress_seen;
    // What's the actual progress bar renderer we're using?
    ProgressBar* progress = nullptr;
};

}

#endif

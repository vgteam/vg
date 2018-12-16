#include "watchdog.hpp"

#include <iostream>
#include <cassert>

#include "memusage.hpp"

namespace vg {

using namespace std;

Watchdog::Watchdog(size_t thread_count, const duration& timeout) :
    stop_watcher(false),
    state(thread_count),
    timeout(timeout),
    memory_high_water_kb(get_max_rss_kb()),
    watcher(&Watchdog::watcher_loop, this) {

    // Nothing to do
}

Watchdog::~Watchdog() {
    // We need to shut down the watcher thread before we can clean up the object.
    stop_watcher = true;
    watcher.join();
}

void Watchdog::check_in(size_t thread, const string& task) {
    // Find the state for the thread we are talking about
    auto& t = state.at(thread); 
    
    // Lock this thread's data
    lock_guard<mutex> lock(t.access_mutex);
    
    if (t.is_checked_in) {
        // Don't let anyone check in twice
        throw runtime_error("vg::Watchdog: Thread " + to_string(thread) +
            " is checked in already on " + t.task_name +
            " but is trying to check in for " + task + "!");
    }
    
    // Record the checkin
    t.is_checked_in = true;
    t.timed_out = false;
    t.last_checkin = clock::now();
    t.task_name = task;
    t.checkin_high_water_kb = memory_high_water_kb;
}

void Watchdog::check_out(size_t thread) {
    // Find the state for the thread we are talking about
    auto& t = state.at(thread); 

    // Lock this thread's data
    lock_guard<mutex> lock(t.access_mutex);
    
    if (!t.is_checked_in) {
        // Don't let anyone check out twice
        throw runtime_error("vg::Watchdog: Thread " + to_string(thread) +
            " is checked out from last task " + t.task_name +
            " but is trying to check out again!");
    }
    
    if (t.timed_out) {
        // The thread already hit the timeout and we reported a warning. We should follow up.
        
        // How long was the thread checked in?
        auto checked_in_duration = clock::now() - t.last_checkin;
        auto checked_in_seconds = chrono::duration_cast<chrono::seconds>(checked_in_duration);
        
        // While it was checked in, how much did the high water memory usage mark rise
        auto high_water_increase_kb = memory_high_water_kb - t.checkin_high_water_kb;
        
        #pragma omp critical (cerr)
        cerr << "warning[vg::Watchdog]: Thread " << thread << " finally checked out after "
            << checked_in_seconds.count() << " seconds and " << high_water_increase_kb
            << " kb memory growth processing: " << t.task_name << endl;
    }
    
    // Record the checkout
    t.is_checked_in = false;
}

void Watchdog::watcher_loop() {
    while (!stop_watcher) {
        // Keep looping until we're asked to shut down
        
        // Update our memory usage estimate
        memory_high_water_kb = get_max_rss_kb(); 
        
        for (size_t i = 0; i < state.size(); i++) {
            // For each thread we are watching
            
            // Find the state for the thread we are talking about
            auto& t = state.at(i); 
            
            {
                // Lock its data
                lock_guard<mutex> lock(t.access_mutex);
                
                if (t.is_checked_in && !t.timed_out) {
                    // How long has it been checked in?
                    auto checked_in_duration = clock::now() - t.last_checkin;
                    
                    if (checked_in_duration > timeout) {
                        // The thread has been checked in too long! Yell about it.
                        auto checked_in_seconds = chrono::duration_cast<chrono::seconds>(checked_in_duration);
                        
                        // While it was checked in, how much did the high water memory usage mark rise
                        auto high_water_increase_kb = memory_high_water_kb - t.checkin_high_water_kb;
                        
                        #pragma omp critical (cerr)
                        cerr << "warning[vg::Watchdog]: Thread " << i << " has been checked in for "
                            << checked_in_seconds.count() << " seconds processing: " << t.task_name << endl;
                            
                        // Only report once per task
                        t.timed_out = true;
                    }
                }
            }
            
            // See if we were asked to stop
            if (stop_watcher) {
                break;
            }
        }
        
        // Work out how long to sleep: 1 second or our timeout, whichever is shorter.
        auto sleep_time = chrono::milliseconds(1000);
        if (timeout < sleep_time) {
            sleep_time = chrono::duration_cast<chrono::milliseconds>(timeout);
        }
        // Sleep and check for timeouts again
        this_thread::sleep_for(sleep_time);
    }
}

}

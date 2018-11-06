#include "watchdog.hpp"

#include <iostream>
#include <cassert>

namespace vg {

using namespace std;

Watchdog::Watchdog(size_t thread_count, const duration& timeout) :
    stop_watcher(false),
    thread_locks(thread_count),
    is_checked_in(thread_count, false),
    last_checkin(thread_count),
    task_name(thread_count),
    timed_out(thread_count, false),
    timeout(timeout),
    watcher(&Watchdog::watcher_loop, this) {

    // Nothing to do
}

Watchdog::~Watchdog() {
    // We need to shut down the watcher thread before we can clean up the object.
    stop_watcher = true;
    watcher.join();
}

void Watchdog::check_in(size_t thread, const string& task) {
    // Lock this thread's data
    lock_guard<mutex> lock(thread_locks.at(thread));
    
    // The thread must be checked out to check in
    assert(is_checked_in[thread] == false);
    
    // Record the checkin
    is_checked_in[thread] = true;
    timed_out[thread] = false;
    last_checkin[thread] = clock::now();
    task_name[thread] = task;
}

void Watchdog::check_out(size_t thread) {
    // Lock this thread's data
    lock_guard<mutex> lock(thread_locks.at(thread));
    
    // The thread must be checked in to check out
    assert(is_checked_in[thread] == true);
    
    // Record the checkout
    is_checked_in[thread] = false;
}

void Watchdog::watcher_loop() {
    while (!stop_watcher) {
        // Keep looping until we're asked to shut down
        for (size_t i = 0; i < thread_locks.size(); i++) {
            // For each thread we are watching
            
            {
                // Lock its data
                lock_guard<mutex> lock(thread_locks[i]);
                
                if (is_checked_in[i] && !timed_out[i]) {
                    // How long has it been checked in?
                    auto checked_in_duration = clock::now() - last_checkin[i];
                    
                    if (checked_in_duration > timeout) {
                        // The thread has been checked in too long! Yell about it.
                        auto checked_in_seconds = chrono::duration_cast<chrono::seconds>(checked_in_duration);
                        cerr << "error[vg::Watchdog]: Thread " << i << " has been checked in for "
                            << checked_in_seconds.count() << " seconds processing: " << task_name[i] << endl;
                            
                        // Only report once per task
                        timed_out[i] = true;
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

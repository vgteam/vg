#ifndef VG_WATCHDOG_HPP_INCLUDED
#define VG_WATCHDOG_HPP_INCLUDED

/**
 * \file watchdog.hpp
 * Defines a watchdog timer that lets us detect reads that are excessively difficult to map.
 */

#include <chrono>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>

namespace vg {

using namespace std;

/**
 * Represents a watchdog timer. Each instance owns its own watching thread.
 * Other threads will check in and check out as they start and complete tasks,
 * and the watchdog thread will complain and possibly terminate the program if
 * a thread stays checked in for too long.
 *
 * All synchronization is managed internally. Threads are responsible for
 * knowing their ID numbers, and we can only handle a certain number of
 * threads.
 */
class Watchdog {

public:
    
    // What time types will we use?
    using clock = chrono::steady_clock;
    using duration = clock::duration;
    using time_point = clock::time_point;
    
    /**
     * Create a Watchdog monitoring the specified number of threads.
     * Complain if any thread is checked in longer than timeout.
     * Automatically starts the monitoring thread.
     */
    Watchdog(size_t thread_count, const duration& timeout);
    
    /**
     * Destroy the Watchdog. Automatically stops and joins on the monitoring thread.
     */
    ~Watchdog();
    
    /**
     * Check the given thread in, to do the given task.
     */
    void check_in(size_t thread, const string& task);
    
    /**
     * Check the given thread out of the task it is checked in for.
     */
    void check_out(size_t thread);
    
private:
    // Since we are accessed by the watcher thread, we can't be copied or moved
    
    Watchdog(const Watchdog& other) = delete;
    Watchdog(Watchdog&& other) = delete;
    
    Watchdog& operator=(const Watchdog& other) = delete;
    Watchdog& operator=(Watchdog&& other) = delete;
    
protected:

    /// Have we been asked to stop the watcher thread?
    atomic<bool> stop_watcher;

    /// Each thread we track has a mutex to control access to its data between it and the watcher thread
    vector<mutex> thread_locks;
    
    /// Is each thread checked in (true) or not?
    vector<bool> is_checked_in;
    /// When did each thread last check in?
    vector<chrono::steady_clock::time_point> last_checkin;
    /// What task did each thread last check into?
    vector<string> task_name;
    /// Have we reported a watchdog timeout since the last checkin for each thread?
    vector<bool> timed_out;
    
    /// How long should we give a task to be checked in before complaining?
    duration timeout;
    
    /// We have a thread that does our watchdogging
    thread watcher;
    
    /// Function run in the watcher thread.
    void watcher_loop();
    
    
   

};

}

#endif

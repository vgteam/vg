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

    // We can't use any vector<bool> data structures protected by these mutexes
    // because they are not guaranteed to allow concurrent access to different
    // elements because of the packing. 
    
    /// Wrap all the per-thread state in a struct to avoid vectors of bools and
    /// simplify things.
    struct thread_state_t {
        /// Lock this before reading or writing any fields.
        /// Regulates access between the thread itself and the watcher thread.
        mutex access_mutex;
        /// Is the thread checked in (true) or not?
        bool is_checked_in = false;
        /// Have we reported a watchdog timeout since the last checkin for the thread?
        bool timed_out = false;
        /// When did the thread last check in?
        chrono::steady_clock::time_point last_checkin;
        /// what was our memory high water mark at last checkin?
        size_t checkin_high_water_kb;
        /// What task did the thread last check into?
        string task_name;
    };
    
    /// Holds the state of each thread, along with its mutex.
    vector<thread_state_t> state;
    
    /// How long should we give a task to be checked in before complaining?
    duration timeout;
    
    /// What's the most recent process memory high water mark estimate?
    /// We report on this when tasks take a long time in case they also are using a lot of memory.
    atomic<size_t> memory_high_water_kb;
    
    /// We have a thread that does our watchdogging
    thread watcher;
    
    /// Function run in the watcher thread.
    void watcher_loop();
    
};

}

#endif

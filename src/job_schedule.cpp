/**
 * \file job_schedule.cpp: implements JobSchedule
 */
#include "job_schedule.hpp"

#include <thread>
#include <chrono>
#include <mutex>
#include <atomic>

#include "utility.hpp"

namespace vg {

using namespace std;

JobSchedule::JobSchedule(const vector<pair<int64_t, int64_t>>& job_requirements,
                         const function<void(int64_t)>& job_func)
    : job_func(job_func)
{
    for (int64_t i = 0; i < job_requirements.size(); ++i) {
        queue.emplace_back(job_requirements[i].second, i);
    }
    // sort in decreasing order by time required
    queue.sort([&](const pair<int64_t, int64_t>& a,
                   const pair<int64_t, int64_t>& b) {
        return job_requirements[a.second].first > job_requirements[b.second].first;
    });
}

void JobSchedule::execute(int64_t target_memory_usage) {
    
    atomic<int64_t> est_memory_usage(0);
    mutex queue_lock;
    int num_threads = get_thread_count();
    vector<thread> workers;
    for (int i = 0; i < num_threads; ++i) {
        workers.emplace_back([&]() {
            while (true) {
                
                int64_t job_memory = -1, job_idx = -1;
                queue_lock.lock();
                if (queue.empty()) {
                    // the queue emptied out while we were waiting
                    queue_lock.unlock();
                    break;
                }
                if (est_memory_usage.load() == 0) {
                    // even if we don't have the memory budget to do this job, we're
                    // going to have to at some point and the memory situation will
                    // never get any better than this
                    tie(job_memory, job_idx) = queue.front();
                    queue.pop_front();
                    est_memory_usage.fetch_add(job_memory);
                }
                else {
                    // find the longest-running job that can be done with the available
                    // memory budget
                    for (auto it = queue.begin(); it != queue.end(); ++it) {
                        if (it->first + est_memory_usage.load() <= target_memory_usage) {
                            tie(job_memory, job_idx) = *it;
                            queue.erase(it);
                            est_memory_usage.fetch_add(job_memory);
                            break;
                        }
                    }
                }
                queue_lock.unlock();
                
                if (job_idx == -1) {
                    // there's nothing we can do right now, so back off a second
                    // before trying again
                    this_thread::sleep_for(chrono::seconds(1));
                }
                else {
                    // we think we have enough memory available to attempt this job
                    job_func(job_idx);
                    est_memory_usage.fetch_sub(job_memory);
                }
            }
        });
    }
    // barrier sync
    for (auto& worker : workers) {
        worker.join();
    }
}
}


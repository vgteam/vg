#ifndef VG_JOB_SECHEDULE_HPP_INCLUDED
#define VG_JOB_SECHEDULE_HPP_INCLUDED

/** \file
 * job_schedule.hpp: defines JobSchedule
 */
#include <vector>
#include <utility>
#include <functional>
#include <queue>
#include <cstdint>
#include <list>

namespace vg {

using namespace std;

/*
 * A parallel job scheduler that tries to (if possible) respect a
 * cap on memory usage. Works best with a moderate number of
 * relatively large jobs.
 */
class JobSchedule {
public:
    
    // job requirements are given in pairs of (time estimate, memory estimate)
    // with the memory estimate being in bytes, and the time estimate in
    // arbitrary units
    // the job function should execute the i-th job when called
    JobSchedule(const vector<pair<int64_t, int64_t>>& job_requirements,
                const function<void(int64_t)>& job_func);
    ~JobSchedule() = default;
    
    // execute the job schedule with a target maximum memory usage
    void execute(int64_t target_memory_usage);
    
private:
    
    function<void(int64_t)> job_func;
    list<pair<int64_t, int64_t>> queue;
    
};

}

#endif

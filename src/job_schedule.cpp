/**
 * \file job_schedule.cpp: implements JobSchedule
 */
#include "job_schedule.hpp"

#include <omp.h>
#include <thread>
#include <chrono>


namespace vg {

using namespace std;

JobSchedule::JobSchedule(const vector<pair<int64_t, int64_t>>& job_requirements,
                         const function<void(int64_t)>& job_func)
    : job_func(job_func)
{
    for (int64_t i = 0; i < job_requirements.size(); ++i) {
        queue.emplace(job_requirements[i].first,
                      job_requirements[i].second,
                      i);
        
    }
    
}

void JobSchedule::execute(int64_t target_memory_usage) {
    
    int64_t jobs_ongoing = 0;
    int64_t est_memory_usage = 0;
    
    auto do_job = [&](tuple<int64_t, int64_t, int64_t> job) {
#pragma omp atomic update
        ++jobs_ongoing;
#pragma omp atomic update
        est_memory_usage += get<1>(job);
        this->job_func(get<2>(job));
#pragma omp atomic update
        est_memory_usage -= get<1>(job);
#pragma omp atomic update
        --jobs_ongoing;
    };
    
#pragma omp parallel default(none) shared(target_memory_usage, jobs_ongoing, est_memory_usage, do_job)
#pragma omp single
    {
        
        while (!queue.empty()) {
            int64_t curr_mem, curr_jobs;
#pragma omp atomic read
            curr_mem = est_memory_usage;
#pragma omp atomic read
            curr_jobs = jobs_ongoing;
            
            // get the remaining job with the highest estimated run time
            auto job = queue.top();
            if (omp_get_num_threads() == 1) {
                // single threaded, just do the job in the scheduler thread
                queue.pop();
                do_job(job);
            }
            else if (curr_jobs == 0 ||
                     (curr_mem + get<1>(queue.top()) < target_memory_usage
                      && curr_jobs < omp_get_num_threads())) {
                // we have memory and threads available to do the next job
                queue.pop();
#pragma omp task default(none) firstprivate(job) shared(jobs_ongoing, est_memory_usage, do_job)
                do_job(job);
            }
            else {
                // TODO: it would be nice to do useful work with this thread, maybe
                // with the smallest job?
                // but performance could crater if this thread isn't available
                // when other threads finish
                
                // wait 1 sec and try again
                this_thread::sleep_for(chrono::seconds(1));
            }
        }
    }
    
}

}


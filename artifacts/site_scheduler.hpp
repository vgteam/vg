#ifndef VG_SITE_SCHEDULER_HPP_INCLUDED
#define VG_SITE_SCHEDULER_HPP_INCLUDED

/** \file
 * The order and the parallel jobs in which a caller calls the sites of a site tree.
 */

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

#include <omp.h>

#include "../src/handle.hpp"
#include "../src/utility.hpp"

namespace vg {

using namespace std;

/**
 * Calls the sites of a tree in parallel, a level at a time.
 *
 * First it calls every top-level site. Each call can queue sites to call after it, usually the
 * site's children. Then it calls the queued sites, then the sites those calls queued, and so on
 * until none are queued.
 *
 * A site's key is the lower of its two boundary node IDs. Sorting by key puts sites with nearby
 * node IDs next to each other.
 *
 * With a batch window of 0, each top-level site is its own parallel job. With a batch window of w
 * node IDs, the top-level sites are grouped into batches instead. A batch holds the sites whose
 * key falls in one window [k*w, (k+1)*w), and is one job, which calls its sites in key order.
 * Queued sites are always one job each, started in key order when batching.
 *
 * `Site` names a site, and is copied into queues and batches.
 */
template<typename Site>
class SiteScheduler {
public:
    /// Lists sites, by calling the function it is given on each.
    using SiteLister = function<void(const function<void(const Site&)>&)>;
    /// Calls one site. `queue` is the calling thread's list of sites to call next, to which the
    /// call adds the sites it wants called.
    using SiteCall = function<void(const Site& site, vector<Site>& queue)>;

    /// Seconds a batch runs before it hands its remaining sites out as tasks.
    ///
    /// Threads that have run out of batches take the tasks, so one slow window does not leave the
    /// other threads idle.
    static constexpr double SPLIT_AFTER_SECONDS = 5.0;
    /// Number of sites in each task a slow batch hands out.
    static constexpr size_t SPLIT_SITES = 16;

    /// Make a scheduler that calls each site with `call`.
    ///
    /// `key` gives a site's key, and is needed only when `batch_window` is not 0.
    SiteScheduler(size_t batch_window, function<nid_t(const Site&)> key, SiteCall call) :
        batch_window(batch_window), key(std::move(key)), call(std::move(call)) {
        // Nothing to do
    }

    /// Call the top-level sites.
    ///
    /// Without batching they come from `for_each_parallel`, which calls its function on several
    /// threads, one job per site. With batching they come from `for_each`.
    void call_top_level(const SiteLister& for_each, const SiteLister& for_each_parallel) {
        if (batch_window == 0) {
            for_each_parallel([this](const Site& site) {
                call(site, queues[omp_get_thread_num()]);
            });
            return;
        }
        vector<Site> roots;
        for_each([&roots](const Site& site) {
            roots.push_back(site);
        });
        sort_by_key(roots);
        const vector<pair<size_t, size_t>> batches = batch_ranges(roots);

#pragma omp parallel for schedule(dynamic, 1)
        for (int b = 0; b < (int)batches.size(); ++b) {
            call_batch(roots, batches[b].first, batches[b].second);
        }
    }

    /// Call the queued sites, a level at a time, until none are queued.
    void call_queued() {
        while (!std::ranges::all_of(queues, [](const vector<Site>& site_vec) {return site_vec.empty();})) {
            vector<Site> cur_queue;
            for (vector<Site>& thread_queue : queues) {
                cur_queue.reserve(cur_queue.size() + thread_queue.size());
                std::ranges::move(thread_queue, std::back_inserter(cur_queue));
                thread_queue.clear();
            }

            if (batch_window > 0) {
                sort_by_key(cur_queue);
            }

#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < cur_queue.size(); ++i) {
                call(cur_queue[i], queues[omp_get_thread_num()]);
            }
        }
    }

private:
    /// Sort sites by key.
    ///
    /// Each site's key is read once, beforehand, because reading a key can follow a pointer into a
    /// scattered record. The sort compares the same keys in the same sequence as sorting the sites
    /// themselves by key would, so it puts sites with equal keys in the same order too.
    void sort_by_key(vector<Site>& sites) const {
        vector<pair<nid_t, Site>> keyed(sites.size());
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < sites.size(); ++i) {
            keyed[i] = make_pair(key(sites[i]), sites[i]);
        }
        // std::ranges::sort may order equal keys differently, which would change the output.
        std::sort(keyed.begin(), keyed.end(),
                  [](const pair<nid_t, Site>& a, const pair<nid_t, Site>& b) {
                      return a.first < b.first;
                  });
        std::ranges::transform(keyed, sites.begin(),
                               [](const pair<nid_t, Site>& keyed_site) { return keyed_site.second; });
    }

    /// Split sites sorted by key into batches, as [begin, end) index ranges, one per window.
    vector<pair<size_t, size_t>> batch_ranges(const vector<Site>& sorted_sites) const {
        vector<pair<size_t, size_t>> batches;
        size_t begin = 0;
        while (begin < sorted_sites.size()) {
            const auto window = (size_t)(key(sorted_sites[begin]) / (nid_t)batch_window);
            size_t end = begin + 1;
            while (end < sorted_sites.size() &&
                   (size_t)(key(sorted_sites[end]) / (nid_t)batch_window) == window) {
                ++end;
            }
            batches.emplace_back(begin, end);
            begin = end;
        }
        return batches;
    }

    /// Call the sites in [begin, end) of `sites` in order, as one job.
    ///
    /// If the job runs longer than `SPLIT_AFTER_SECONDS`, hand its remaining sites out as tasks of
    /// `SPLIT_SITES` sites, and wait for them. This thread runs tasks while it waits.
    void call_batch(const vector<Site>& sites, size_t begin, size_t end) {
        const auto batch_start = std::chrono::steady_clock::now();
        size_t i = begin;
        while (i < end) {
            call(sites[i], queues[omp_get_thread_num()]);
            ++i;
            if (i < end && std::chrono::duration<double>(std::chrono::steady_clock::now()
                                                         - batch_start).count()
                               > SPLIT_AFTER_SECONDS) {
                break;
            }
        }
        for (size_t first = i; first < end; first += SPLIT_SITES) {
            const size_t last = std::min(end, first + SPLIT_SITES);
#pragma omp task firstprivate(first, last) shared(sites)
            for (size_t k = first; k < last; ++k) {
                call(sites[k], queues[omp_get_thread_num()]);
            }
        }
#pragma omp taskwait
    }

    size_t batch_window;
    function<nid_t(const Site&)> key;
    SiteCall call;
    /// One queue per thread, of the sites to call at the next level.
    vector<vector<Site>> queues = vector<vector<Site>>(get_thread_count());
};

/**
 * A tally of the sites a walk has called, which reports its progress to stderr.
 *
 * It counts top-level sites until `finish_top_level()`, and nested sites after. It reports every
 * 100000th site of each kind, and each total when that kind is finished. It is safe to count from
 * several threads.
 */
class SiteCallTally {
public:
    /// Count one site called.
    void count() {
        if (top_level) {
            ++top_snarl_count;
            if (top_snarl_count % 100000 == 0) {
#pragma omp critical (cerr)
                cerr << "[vg call]: Processed " << top_snarl_count << " top-level snarls" << endl;
            }
        } else {
            std::int64_t counted = ++nested_snarl_count;
            if (counted % 100000 == 0) {
#pragma omp critical (cerr)
                cerr << "[vg call]: Processed " << counted << " nested snarls" << endl;
            }
        }
    }

    /// Report the top-level total, and count later sites as nested.
    ///
    /// No thread may be counting while this runs.
    void finish_top_level() {
        cerr << "[vg call]: Finished processing " << top_snarl_count << " top-level snarls" << endl;
        top_level = false;
    }

    /// Report the nested total, if any nested sites were counted.
    void finish_nested() const {
        if (nested_snarl_count > 0) {
            cerr << "[vg call]: Finished processing " << nested_snarl_count << " nested snarls" << endl;
        }
    }

private:
    std::atomic<std::int64_t> top_snarl_count{0};
    std::atomic<std::int64_t> nested_snarl_count{0};
    /// Whether sites counted now are top-level sites.
    bool top_level = true;
};

/**
 * Which children of the sites in a forest a walk goes on to call.
 */
enum class RecurseType {
    /// Call the children of each site whose call failed.
    OnFail,
    /// Call the children of every site.
    Always,
    /// Call only the top-level sites.
    Never
};

/**
 * A view of the nested sites of a graph, which lists the top-level sites and the children of each
 * site.
 *
 * Each top-level site is the root of one tree of sites.
 *
 * `Site` names a site, and is copied into queues and batches.
 */
template<typename Site>
class SiteForest {
public:
    virtual ~SiteForest() = default;

    /// Return the key of a site, the lower of its two boundary node IDs.
    virtual nid_t key(const Site& site) const = 0;

    /// Call `iteratee` on each top-level site, on this thread.
    virtual void for_each_top(const function<void(const Site&)>& iteratee) const = 0;

    /// Call `iteratee` on each top-level site, spread over the OpenMP threads.
    virtual void for_each_top_parallel(const function<void(const Site&)>& iteratee) const = 0;

    /// Return true if a site has nothing to call.
    ///
    /// A walk neither calls nor counts such a site, and does not go on to its children.
    virtual bool is_trivial(const Site& site) const = 0;

    /// Call `iteratee` on each child of a site.
    virtual void for_each_child(const Site& site, const function<void(const Site&)>& iteratee) const = 0;
};

/**
 * Call the sites of a forest from the top level down.
 *
 * Call every top-level site in parallel, then the children `recurse_type` asks for, and so on
 * down. Run the calls in the jobs that a `SiteScheduler` with `batch_window` gives.
 *
 * @param call Calls a site, and returns false if it failed.
 * @param show_progress Report the sites called to stderr.
 */
template<typename Site>
void walk_site_tree(const SiteForest<Site>& forest, RecurseType recurse_type, size_t batch_window,
                    bool show_progress, const function<bool(const Site&)>& call) {
    SiteCallTally tally;

    auto process_site = [&forest, &call, recurse_type, show_progress, &tally](
                            const Site& site, vector<Site>& thread_queue) {
        if (forest.is_trivial(site)) {
            return;
        }
        if (bool was_called = call(site);
            recurse_type == RecurseType::Always || (!was_called && recurse_type == RecurseType::OnFail)) {
            forest.for_each_child(site, [&thread_queue](const Site& child) {
                thread_queue.push_back(child);
            });
        }
        if (show_progress) {
            tally.count();
        }
    };
    SiteScheduler<Site> scheduler(batch_window, [&forest](const Site& site) { return forest.key(site); },
                                  process_site);

    scheduler.call_top_level(
        [&forest](const function<void(const Site&)>& iteratee) { forest.for_each_top(iteratee); },
        [&forest](const function<void(const Site&)>& iteratee) { forest.for_each_top_parallel(iteratee); });
    if (show_progress) {
        tally.finish_top_level();
    }

    scheduler.call_queued();
    if (show_progress) {
        tally.finish_nested();
    }
}

}

#endif

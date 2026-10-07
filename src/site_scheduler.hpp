#ifndef VG_SITE_SCHEDULER_HPP_INCLUDED
#define VG_SITE_SCHEDULER_HPP_INCLUDED

/** \file
 * The order and the parallel jobs in which a caller calls the sites of a site tree.
 */

#include <algorithm>
#include <functional>
#include <iterator>
#include <utility>
#include <vector>

#include <omp.h>

#include "handle.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

/**
 * Calls the sites of a tree in parallel: first every top-level site, then the sites those calls
 * queued, then the sites those calls queued, and so on until none are queued. A call queues the
 * sites it wants called after it, usually its children.
 *
 * By default each site is its own parallel job. With a batch window of w node IDs, the top-level
 * sites are grouped into batches instead: a batch holds the sites whose key falls in one window of
 * w IDs, [k*w, (k+1)*w) for some k, each batch is one job, and a job calls its sites in key order.
 * Queued sites are still one job each, started in key order. A site's key is the lower of its two
 * boundary node IDs, so that key order puts sites with nearby node IDs next to each other.
 *
 * `Site` names a site, and is cheap to copy.
 */
template<typename Site>
class SiteScheduler {
public:
    /// Lists sites, by calling the function it is given on each.
    using SiteLister = function<void(const function<void(const Site&)>&)>;
    /// Calls one site. `queue` is the calling thread's list of sites to call next, to which the
    /// call adds the sites it wants called.
    using SiteCall = function<void(const Site& site, vector<Site>& queue)>;

    /// `key` gives a site's key, and is needed only when `batch_window` is not 0.
    SiteScheduler(size_t batch_window, function<nid_t(const Site&)> key, SiteCall call) :
        batch_window(batch_window), key(std::move(key)), call(std::move(call)),
        queues(get_thread_count()) {
        // Nothing to do
    }

    /// Call the top-level sites. Without batching they come from `for_each_parallel`, which calls
    /// its function on several threads, one job per site; with batching they come from `for_each`.
    void call_top_level(const SiteLister& for_each, const SiteLister& for_each_parallel) {
        auto call_site = [&](const Site& site) {
            call(site, queues[omp_get_thread_num()]);
        };
        if (batch_window == 0) {
            for_each_parallel(call_site);
            return;
        }
        // Call in node-ID order, batched by window, so that sites with nearby node IDs are
        // called together.
        vector<Site> roots;
        for_each([&](const Site& site) {
            roots.push_back(site);
        });
        sort_by_key(roots);

        // Split the sorted sites into runs that share a window, with one parallel job per run
        // rather than per site.
        vector<pair<size_t, size_t>> windows;
        size_t begin = 0;
        while (begin < roots.size()) {
            size_t window = (size_t)(key(roots[begin]) / (nid_t)batch_window);
            size_t end = begin + 1;
            while (end < roots.size() &&
                   (size_t)(key(roots[end]) / (nid_t)batch_window) == window) {
                ++end;
            }
            windows.emplace_back(begin, end);
            begin = end;
        }

#pragma omp parallel for schedule(dynamic, 1)
        for (int w = 0; w < (int)windows.size(); ++w) {
            for (size_t i = windows[w].first; i < windows[w].second; ++i) {
                call_site(roots[i]);
            }
        }
    }

    /// Call the queued sites, a level at a time, until none are queued.
    void call_queued() {
        while (!std::all_of(queues.begin(), queues.end(),
                            [](const vector<Site>& site_vec) {return site_vec.empty();})) {
            vector<Site> cur_queue;
            for (vector<Site>& thread_queue : queues) {
                cur_queue.reserve(cur_queue.size() + thread_queue.size());
                std::move(thread_queue.begin(), thread_queue.end(), std::back_inserter(cur_queue));
                thread_queue.clear();
            }

            if (batch_window > 0) {
                // Keep queued sites in node-ID order too.
                sort_by_key(cur_queue);
            }

#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < cur_queue.size(); ++i) {
                call(cur_queue[i], queues[omp_get_thread_num()]);
            }
        }
    }

private:
    /// Sort sites by key, so that sites with nearby node IDs are called one after another.
    ///
    /// Each site's key is read once, beforehand: reading it at every comparison followed a pointer
    /// into each site's scattered record, millions of times for a whole genome's top-level sites.
    /// The sort compares the same keys in the same sequence, so the result, including the order of
    /// sites with equal keys, is the one sorting the sites themselves gives.
    void sort_by_key(vector<Site>& sites) const {
        vector<pair<nid_t, Site>> keyed(sites.size());
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < sites.size(); ++i) {
            keyed[i] = make_pair(key(sites[i]), sites[i]);
        }
        std::sort(keyed.begin(), keyed.end(),
                  [](const pair<nid_t, Site>& a, const pair<nid_t, Site>& b) {
                      return a.first < b.first;
                  });
        for (size_t i = 0; i < sites.size(); ++i) {
            sites[i] = keyed[i].second;
        }
    }

    size_t batch_window;
    function<nid_t(const Site&)> key;
    SiteCall call;
    /// One queue per thread, of the sites to call at the next level.
    vector<vector<Site>> queues;
};

}

#endif

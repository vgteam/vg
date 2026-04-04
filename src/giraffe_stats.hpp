#ifndef VG_GIRAFFE_STATS_HPP_INCLUDED
#define VG_GIRAFFE_STATS_HPP_INCLUDED

/**
 * \file giraffe_stats.hpp
 * Aggregate per-stage timing statistics across reads for the Giraffe mapper.
 */

#include "funnel.hpp"

#include <string>
#include <vector>
#include <unordered_map>
#include <ostream>

namespace vg {

using namespace std;

/**
 * Thread-safe aggregate statistics collector for the Giraffe mapper.
 *
 * Each thread pushes data into its own ThreadData slot (no locking on the hot
 * path). Summary statistics are computed and printed after all reads are
 * mapped.
 *
 * Also logs a per-stage breakdown to stderr for any individual read whose
 * total mapping time exceeds a configurable threshold.
 */
class GiraffeStats {
public:
    GiraffeStats(size_t thread_count, double slow_threshold_s);

    /**
     * Record timing/item data from one read's Funnel.
     * Must be called after funnel.stop().
     * Uses omp_get_thread_num() to select the per-thread slot.
     */
    void record_read(const Funnel& funnel, const std::string& read_name);

    /**
     * Print a summary table of per-stage timing percentiles to out.
     * Thread-safe to call after all record_read() calls are done.
     */
    void print_summary(std::ostream& out) const;

private:
    struct ThreadData {
        // Per-stage accumulated durations and item counts.
        // Stage names appear in insertion order the first time they're seen.
        std::vector<std::string> stage_order;
        std::unordered_map<std::string, std::vector<double>> stage_durations;
        std::unordered_map<std::string, std::vector<size_t>> stage_item_counts;

        // Per-substage accumulated durations, keyed as "stage/substage".
        std::unordered_map<std::string, std::vector<double>> substage_durations;

        // Total per-read duration.
        std::vector<double> read_durations;

        size_t slow_read_count = 0;
    };

    std::vector<ThreadData> per_thread;
    double slow_threshold_s;

    // Compute percentile p in [0,1] from a sorted vector.
    static double percentile(const std::vector<double>& sorted, double p);
};

} // namespace vg

#endif // VG_GIRAFFE_STATS_HPP_INCLUDED

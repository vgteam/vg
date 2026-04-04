/**
 * \file giraffe_stats.cpp
 */

#include "giraffe_stats.hpp"

#include <omp.h>
#include <algorithm>
#include <format>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <iostream>

namespace vg {

GiraffeStats::GiraffeStats(size_t thread_count, double slow_threshold_s)
    : per_thread(thread_count), slow_threshold_s(slow_threshold_s) {}

void GiraffeStats::record_read(const Funnel& funnel, const std::string& read_name) {
    int tid = omp_get_thread_num();
    ThreadData& td = per_thread.at(tid);

    double total = funnel.total_seconds();
    td.read_durations.push_back(total);

    bool is_slow = slow_threshold_s > 0.0 && total > slow_threshold_s;
    if (is_slow) {
        td.slow_read_count++;
    }

    // Accumulate per-stage data and optionally build the slow-read message.
    std::ostringstream slow_msg;
    if (is_slow) {
        slow_msg << std::format("warning[vg::Giraffe]: Slow read \"{}\" took {:.3f}s:\n", read_name, total);
    }

    funnel.for_each_stage([&](const std::string& stage,
                               const std::vector<size_t>& result_sizes,
                               const std::vector<double>& /*correct_scores*/,
                               const std::vector<double>& /*noncorrect_scores*/,
                               const double& duration,
                               const std::unordered_map<std::string, double>& sub_durations) {
        // Record into per-thread storage.
        if (!td.stage_durations.count(stage)) {
            td.stage_order.push_back(stage);
        }
        td.stage_durations[stage].push_back(duration);
        td.stage_item_counts[stage].push_back(result_sizes.size());

        for (auto& kv : sub_durations) {
            std::string key = stage + "/" + kv.first;
            td.substage_durations[key].push_back(kv.second);
        }

        if (is_slow) {
            slow_msg << std::format("  {}: {:.3f}s ({} items)\n", stage, duration, result_sizes.size());
            for (auto& kv : sub_durations) {
                slow_msg << "    " << kv.first << ": "
                         << std::fixed << std::setprecision(3) << kv.second << "s\n";
            }
        }
    });

    if (is_slow) {
        #pragma omp critical (cerr)
        std::cerr << slow_msg.str() << std::flush;
    }
}

double GiraffeStats::percentile(const std::vector<double>& sorted, double p) {
    if (sorted.empty()) return 0.0;
    if (sorted.size() == 1) return sorted[0];
    double idx = p * (sorted.size() - 1);
    size_t lo = (size_t)idx;
    size_t hi = lo + 1;
    if (hi >= sorted.size()) return sorted.back();
    double frac = idx - lo;
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

void GiraffeStats::print_summary(std::ostream& out) const {
    // Merge all thread data.
    // Use the stage order from thread 0 (or whichever thread saw stages first),
    // then append any stages seen by other threads.
    std::vector<std::string> stage_order;
    std::unordered_map<std::string, std::vector<double>, StringHash, StringEqual> stage_durations;
    std::unordered_map<std::string, std::vector<size_t>, StringHash, StringEqual> stage_item_counts;
    std::unordered_map<std::string, std::vector<double>, StringHash, StringEqual> substage_durations;
    std::vector<double> read_durations;
    size_t slow_read_count = 0;
    size_t total_reads = 0;

    for (auto& td : per_thread) {
        // Merge stage order (preserve first-seen ordering).
        for (auto& s : td.stage_order) {
            if (!stage_durations.count(s)) {
                stage_order.push_back(s);
            }
        }
        for (auto& kv : td.stage_durations) {
            auto& dst = stage_durations[kv.first];
            dst.insert(dst.end(), kv.second.begin(), kv.second.end());
        }
        for (auto& kv : td.stage_item_counts) {
            auto& dst = stage_item_counts[kv.first];
            dst.insert(dst.end(), kv.second.begin(), kv.second.end());
        }
        for (auto& kv : td.substage_durations) {
            auto& dst = substage_durations[kv.first];
            dst.insert(dst.end(), kv.second.begin(), kv.second.end());
        }
        read_durations.insert(read_durations.end(),
                              td.read_durations.begin(), td.read_durations.end());
        slow_read_count += td.slow_read_count;
        total_reads += td.read_durations.size();
    }

    if (total_reads == 0) return;

    // Sort per-read durations for percentiles.
    std::vector<double> sorted_reads = read_durations;
    std::sort(sorted_reads.begin(), sorted_reads.end());

    // Print header.
    out << "\n=== Giraffe Per-Stage Timing (" << total_reads << " reads) ===\n";
    out << std::left
        << std::setw(28) << "Stage"
        << std::right
        << std::setw(10) << "Mean(ms)"
        << std::setw(10) << "P50(ms)"
        << std::setw(10) << "P95(ms)"
        << std::setw(10) << "P99(ms)"
        << std::setw(10) << "Max(ms)"
        << std::setw(12) << "MeanItems"
        << std::setw(10) << "MaxItems"
        << "\n";

    // Helper to print one row.
    auto print_row = [&](const std::string& label, std::vector<double>& durs,
                         const std::vector<size_t>* items) {
        std::sort(durs.begin(), durs.end());
        double mean_ms = 1000.0 * std::accumulate(durs.begin(), durs.end(), 0.0) / durs.size();
        double p50_ms  = 1000.0 * percentile(durs, 0.50);
        double p95_ms  = 1000.0 * percentile(durs, 0.95);
        double p99_ms  = 1000.0 * percentile(durs, 0.99);
        double max_ms  = 1000.0 * durs.back();

        out << std::left << std::setw(28) << label << std::right
            << std::fixed << std::setprecision(2)
            << std::setw(10) << mean_ms
            << std::setw(10) << p50_ms
            << std::setw(10) << p95_ms
            << std::setw(10) << p99_ms
            << std::setw(10) << max_ms;

        if (items && !items->empty()) {
            double mean_items = (double)std::accumulate(items->begin(), items->end(), (size_t)0) / items->size();
            size_t max_items  = *std::max_element(items->begin(), items->end());
            out << std::setw(12) << std::setprecision(1) << mean_items
                << std::setw(10) << max_items;
        }
        out << "\n";
    };

    for (auto& stage : stage_order) {
        print_row(stage, stage_durations[stage], &stage_item_counts[stage]);

        // Print any substages for this stage, indented.
        for (auto& kv : substage_durations) {
            // Key format is "stage/substage".
            size_t slash = kv.first.find('/');
            if (slash == std::string::npos) continue;
            if (kv.first.substr(0, slash) != stage) continue;
            std::string substage_label = "  " + kv.first.substr(slash + 1);
            print_row(substage_label, kv.second, nullptr);
        }
    }

    // Total read timing.
    double mean_ms = 1000.0 * std::accumulate(read_durations.begin(), read_durations.end(), 0.0) / total_reads;
    double p99_ms  = 1000.0 * percentile(sorted_reads, 0.99);
    double max_ms  = 1000.0 * sorted_reads.back();

    out << "\nTotal per-read: mean=" << std::fixed << std::setprecision(2) << mean_ms
        << "ms, P99=" << p99_ms << "ms, max=" << max_ms << "ms\n";

    if (slow_threshold_s > 0.0) {
        out << "Slow reads (>" << slow_threshold_s << "s): "
            << slow_read_count << " / " << total_reads
            << " (" << std::fixed << std::setprecision(3)
            << 100.0 * slow_read_count / total_reads << "%)\n";
    }
    out << std::flush;
}

} // namespace vg

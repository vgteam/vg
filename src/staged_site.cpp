#include "flow_caller.hpp"

namespace vg {

/// Total over the per-thread queues.
template <typename Queues>
static size_t total_queued(const Queues& queues) {
    size_t n = 0;
    for (const auto& queue : queues) {
        n += queue.size();
    }
    return n;
}

size_t FlowCaller::pending_record_count() const {
    return total_queued(pending_records);
}

size_t FlowCaller::render_record_count() const {
    return total_queued(render_records);
}

const vector<int>& FlowCaller::cached_panel_alleles(PendingRecord& rec) {
    if (!rec.panel_cached) {
        rec.panel_cache = panel_lookup.alleles(rec.travs);
        rec.panel_cached = true;
    }
    return rec.panel_cache;
}

vector<FlowCaller::PendingRecord*> FlowCaller::records_for_render(bool for_phasing) {
    vector<PendingRecord*> out;
    out.reserve(render_record_count() + deferred_pending.size());
    for (auto& queue : render_records) {
        for (PendingRecord& rec : queue) {
            out.push_back(&rec);
        }
    }
    for (PendingRecord& rec : deferred_pending) {
        // The same records the hand-off holds back, so that this matches what is rendered: a
        // dropped chain, which the parent's chosen genotype does not carry, and a
        // `reported_inline` one, which an enclosing block's ALT spells. Tested each time, since both
        // can change between linkage passes.
        if (rec.dropped || rec.reported_inline) {
            continue;
        }
        // A chain with no reference path has no REF or POS, so it is not rendered, but it is
        // genotyped, gets anchors, and has a meaningful strand, so phasing includes it.
        if (rec.no_reference && !for_phasing) {
            continue;
        }
        out.push_back(&rec);
    }
    return out;
}

}

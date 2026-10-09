#include <array>
#include <unordered_map>

#include "flow_caller.hpp"

namespace vg {

unordered_map<size_t, array<int, 3>> FlowCaller::chosen_snapshot() {
    // Each record's chosen pair and ploidy, keyed by record key.
    unordered_map<size_t, array<int, 3>> out;
    if (!linker.enabled()) {
        return out;
    }
    // Looked up on several threads, then filed in record order, so that the map is built exactly
    // as one loop over the records builds it.
    const vector<StagedSite*> records = staged_sites.in_order();
    vector<size_t> keys(records.size());
    for (size_t i = 0; i < records.size(); ++i) {
        keys[i] = records[i]->record_key;
    }
    vector<array<int, 3>> chosen;
    vector<char> found;
    linker.collector()->chosen_traversals_for(keys, chosen, found);
    for (size_t i = 0; i < keys.size(); ++i) {
        if (found[i]) {
            out[keys[i]] = chosen[i];
        }
    }
    return out;
}

size_t FlowCaller::snapshot_digest(const unordered_map<size_t, array<int, 3>>& snap) {
    // Independent of order, since the snapshot is a hash map: each record's contribution is
    // combined with a commutative mix.
    size_t acc = snap.size() * 1000003ULL;
    for (const auto& kv : snap) {
        size_t h = kv.first;
        h = h * 1000003ULL + (size_t)(kv.second[0] + 3);
        h = h * 1000003ULL + (size_t)(kv.second[1] + 3);
        h = h * 1000003ULL + (size_t)kv.second[2];
        acc ^= h + 0x9e3779b97f4a7c15ULL + (acc << 6) + (acc >> 2);
    }
    return acc;
}

size_t FlowCaller::chosen_changed(const unordered_map<size_t, array<int, 3>>& before,
                                  const unordered_map<size_t, array<int, 3>>& after) {
    size_t moved = 0;
    for (const auto& kv : after) {
        auto found = before.find(kv.first);
        if (found == before.end() || found->second != kv.second) {
            ++moved;
        }
    }
    // A record that had a chosen answer and now has none has changed too.
    for (const auto& kv : before) {
        if (after.count(kv.first) == 0) {
            ++moved;
        }
    }
    return moved;
}

}

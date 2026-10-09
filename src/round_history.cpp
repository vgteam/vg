#include <array>
#include <unordered_map>

#include "round_history.hpp"

namespace vg {

RoundHistory::State RoundHistory::state(StagedSiteTable& sites, const LinkageCollector* model) {
    State out;
    if (model == nullptr) {
        return out;
    }
    // Looked up on several threads, then filed in record order, so that the map is built exactly
    // as one loop over the records builds it.
    const vector<StagedSite*> records = sites.in_order();
    vector<size_t> keys(records.size());
    for (size_t i = 0; i < records.size(); ++i) {
        keys[i] = records[i]->record_key;
    }
    vector<array<int, 3>> chosen;
    vector<char> found;
    model->chosen_traversals_for(keys, chosen, found);
    for (size_t i = 0; i < keys.size(); ++i) {
        if (found[i]) {
            out[keys[i]] = chosen[i];
        }
    }
    return out;
}

size_t RoundHistory::digest(const State& snap) {
    // Each site is mixed in in the map's iteration order, so two equal states give the same digest
    // when they were built alike. `state` builds every state the same way, in record order.
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

size_t RoundHistory::changed(const State& before, const State& after) {
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

void RoundHistory::remember(const State& state) {
    digests.push_back(digest(state));
}

size_t RoundHistory::first_round_with(const State& state) const {
    const size_t d = digest(state);
    for (size_t i = 0; i < digests.size(); ++i) {
        if (digests[i] == d) {
            return i + 1;
        }
    }
    return 0;
}

}

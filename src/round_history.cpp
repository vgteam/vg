#include <array>
#include <unordered_map>

#include "round_history.hpp"

namespace vg::multipass {

RoundHistory::State RoundHistory::state(StagedSiteTable& sites, const LinkageCollector* model) {
    State out;
    if (model == nullptr) {
        return out;
    }
    // Look up the choices on several threads, then file them in record order, so that every state
    // is built by the same sequence of insertions and iterates in the same order.
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
    // Mix in the sites in the map's iteration order. Equal states then have equal digests only if
    // they were built by the same insertions, which `state()` guarantees.
    size_t acc = snap.size() * 1000003ULL;
    for (const auto& [key, choice] : snap) {
        size_t h = key;
        h = h * 1000003ULL + (size_t)(choice[0] + 3);
        h = h * 1000003ULL + (size_t)(choice[1] + 3);
        h = h * 1000003ULL + (size_t)choice[2];
        acc ^= h + 0x9e3779b97f4a7c15ULL + (acc << 6) + (acc >> 2);
    }
    return acc;
}

size_t RoundHistory::changed(const State& before, const State& after) {
    size_t moved = 0;
    for (const auto& [key, choice] : after) {
        auto found = before.find(key);
        if (found == before.end() || found->second != choice) {
            ++moved;
        }
    }
    for (const auto& [key, choice] : before) {
        if (!after.contains(key)) {
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

#include <algorithm>
#include <array>
#include <unordered_map>

#include "round_history.hpp"

namespace vg::multipass {

RoundHistory::State RoundHistory::state(StagedSiteTable& sites, const LinkageCollector* model) {
    State out;
    if (model == nullptr) {
        return out;
    }
    // Look up the choices on several threads, then file them.
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

RoundHistory::SortedState RoundHistory::sorted(const State& state) {
    SortedState out(state.begin(), state.end());
    std::ranges::sort(out, {}, &SortedState::value_type::first);
    return out;
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
    remembered.push_back(sorted(state));
}

size_t RoundHistory::first_round_with(const State& state) const {
    const SortedState wanted = sorted(state);
    for (size_t i = 0; i < remembered.size(); ++i) {
        if (remembered[i] == wanted) {
            return i + 1;
        }
    }
    return 0;
}

}

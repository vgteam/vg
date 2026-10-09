#include <algorithm>
#include <iterator>

#include <omp.h>

#include "staged_site.hpp"

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

void StagedSite::set_call(unique_ptr<SnarlCaller::CallInfo> info, SiteScore* typed) {
    call_info = std::move(info);
    score = typed;
}

void StagedSite::set_score(unique_ptr<SiteScore> typed) {
    score = typed.get();
    call_info = std::move(typed);
}

unique_ptr<SiteScore> StagedSite::take_score() {
    unique_ptr<SiteScore> out(score);
    // `out` owns it now.
    call_info.release();
    score = nullptr;
    return out;
}

const vector<int>& StagedSite::panel_alleles(const PanelLookup& lookup) {
    if (!panel_cached) {
        panel_cache = lookup.alleles(travs);
        panel_cached = true;
    }
    return panel_cache;
}

void StagedSiteTable::start(size_t threads) {
    // resize, not assign: StagedSite holds a unique_ptr, so it cannot be copied.
    nested_lists.clear();
    nested_lists.resize(max(threads, (size_t)1));
    queues.clear();
    queues.resize(max(threads, (size_t)1));
}

void StagedSiteTable::add_top_level(StagedSite&& site) {
    queues[omp_get_thread_num()].push_back(std::move(site));
}

void StagedSiteTable::add_nested(StagedSite&& site) {
    nested_lists[omp_get_thread_num()].push_back(std::move(site));
}

vector<StagedSite>& StagedSiteTable::gather_nested() {
    const size_t before = nested_sites.size();
    nested_sites.reserve(nested_sites.size() + total_queued(nested_lists));
    for (auto& list : nested_lists) {
        std::move(list.begin(), list.end(), std::back_inserter(nested_sites));
        list.clear();
    }
    if (nested_sites.size() != before) {
        // Built from every gathered site, so that dropping a chain can drop everything under it.
        children.clear();
        children.reserve(nested_sites.size() * 2);
        for (size_t i = 0; i < nested_sites.size(); ++i) {
            children[nested_sites[i].parent_record_key].push_back(i);
        }
    }
    return nested_sites;
}

const vector<size_t>* StagedSiteTable::children_of(size_t parent_key) const {
    auto found = children.find(parent_key);
    return found == children.end() ? nullptr : &found->second;
}

unordered_map<size_t, StagedSite*> StagedSiteTable::by_key() {
    unordered_map<size_t, StagedSite*> out;
    out.reserve((nested_sites.size() + queued_count()) * 2);
    for (StagedSite& site : nested_sites) {
        out[site.record_key] = &site;
    }
    for (auto& queue : queues) {
        for (StagedSite& site : queue) {
            out[site.record_key] = &site;
        }
    }
    return out;
}

void StagedSiteTable::for_each_parent_top_down(
    const unordered_map<size_t, StagedSite*>& by_key,
    const function<void(const StagedSite& parent, const vector<size_t>& children)>& visit) const {
    vector<pair<uint8_t, size_t>> parents;
    parents.reserve(children.size());
    for (const auto& kv : children) {
        auto parent = by_key.find(kv.first);
        if (parent != by_key.end()) {
            parents.emplace_back(parent->second->level, kv.first);
        }
    }
    sort(parents.begin(), parents.end());
    for (const pair<uint8_t, size_t>& level_key : parents) {
        visit(*by_key.at(level_key.second), children.at(level_key.second));
    }
}

vector<StagedSite*> StagedSiteTable::in_order(bool with_off_reference) {
    vector<StagedSite*> out;
    out.reserve(queued_count() + nested_sites.size());
    for (auto& queue : queues) {
        for (StagedSite& site : queue) {
            out.push_back(&site);
        }
    }
    for (StagedSite& site : nested_sites) {
        if (site.dropped || site.reported_inline) {
            continue;
        }
        if (site.no_reference && !with_off_reference) {
            continue;
        }
        out.push_back(&site);
    }
    return out;
}

void StagedSiteTable::for_each(const function<void(StagedSite&)>& visit) {
    for (auto& queue : queues) {
        for (StagedSite& site : queue) {
            visit(site);
        }
    }
    for (StagedSite& site : nested_sites) {
        visit(site);
    }
}

StagedSiteTable::HandOff StagedSiteTable::hand_off() {
    HandOff out;
    size_t next_queue = 0;
    for (StagedSite& site : nested_sites) {
        if (site.dropped) {
            continue;
        }
        if (site.reported_inline) {
            ++out.inline_unrendered;
            continue;
        }
        if (site.no_reference) {
            ++out.no_ref_unrendered;
            continue;
        }
        queues[next_queue % queues.size()].push_back(std::move(site));
        ++next_queue;
    }
    nested_sites.clear();
    children.clear();
    return out;
}

size_t StagedSiteTable::queued_count() const {
    return total_queued(queues);
}

}

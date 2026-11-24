#include <algorithm>
#include <iterator>
#include <map>
#include <tuple>

#include <omp.h>

#include "staged_site.hpp"
#include "vcf_record.hpp"

namespace vg::multipass {

/// Count the sites in all of the per-thread queues.
template <typename Queues>
static size_t total_queued(const Queues& queues) {
    size_t n = 0;
    for (const auto& queue : queues) {
        n += queue.size();
    }
    return n;
}

/// Turn a path name and a position on it into a locus, clamping the position at 0.
static SiteLocus locus_of(pair<string, int64_t> pos_info) {
    if (const string locus = PathMetadata::parse_locus_name(pos_info.first);
        locus != PathMetadata::NO_LOCUS_NAME) {
        pos_info.first = locus;
    }
    return SiteLocus{pos_info.first, (size_t)max((int64_t)0, pos_info.second)};
}

SiteLocus ReferencePlacement::locus(const PathPositionHandleGraph& graph,
                                    const EditSnarlDecompositionOverlay& sites,
                                    const net_handle_t& net) const {
    if (no_reference) {
        // `get_ref_position` would fail an assertion on a site the path does not visit.
        return locus_of(make_pair(path_name, position_from_parent));
    }
    // Place a site inside a node at its own base, so that the sites inside one node are ordered
    // and spaced. Its bounds read the node the way the reference path does.
    pair<string, int64_t> position = get_ref_position(graph, sites.get_site_start(net),
                                                      sites.get_site_end(net), path_name, offset);
    if (sites.is_edit(net)) {
        position.second += (int64_t)sites.get_site_start_offset(net);
    }
    return locus_of(position);
}

void StagedSite::set_call(unique_ptr<SnarlCaller::CallInfo> new_call, SiteScore* typed) {
    call_info = std::move(new_call);
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
    if (!panel_cache) {
        panel_cache = lookup.alleles(alleles);
    }
    return *panel_cache;
}

void StagedSiteTable::start(size_t threads) {
    // Use resize, not assign, because StagedSite holds a unique_ptr, so it cannot be copied.
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
        std::ranges::move(list, std::back_inserter(nested_sites));
        list.clear();
    }
    if (nested_sites.size() != before) {
        // Index every gathered site, so that dropping a chain can drop everything under it.
        children.clear();
        children.reserve(nested_sites.size() * 2);
        for (size_t i = 0; i < nested_sites.size(); ++i) {
            children[nested_sites[i].placement.parent_record_key].push_back(i);
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
    vector<pair<size_t, size_t>> parents;
    parents.reserve(children.size());
    for (const auto& [parent_key, kids] : children) {
        if (auto parent = by_key.find(parent_key); parent != by_key.end()) {
            parents.emplace_back(parent->second->placement.level, parent_key);
        }
    }
    std::ranges::sort(parents);
    for (const auto& [level, parent_key] : parents) {
        visit(*by_key.at(parent_key), children.at(parent_key));
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
        if (site.dropped || site.placement.reported_inline) {
            continue;
        }
        if (site.reference.no_reference && !with_off_reference) {
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
        if (site.placement.reported_inline || site.reference.no_reference) {
            if (site.placement.reported_inline) {
                ++out.inline_unrendered;
            } else {
                ++out.no_ref_unrendered;
            }
            // Keep the site for its parent's record, which spells its chosen allele, and free its
            // scores and its kept panel alleles.
            site.set_call(nullptr, nullptr);
            site.panel_cache.reset();
            held_back.push_back(std::move(site));
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

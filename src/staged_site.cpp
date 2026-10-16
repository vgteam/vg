#include <algorithm>
#include <iterator>
#include <map>
#include <tuple>

#include <omp.h>

#include "staged_site.hpp"

namespace vg {
namespace multipass {

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
        panel_cache = lookup.alleles(travs, &sequences);
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
        if (site.reported_inline || site.no_reference) {
            if (site.reported_inline) {
                ++out.inline_unrendered;
            } else {
                ++out.no_ref_unrendered;
            }
            // Kept for its parent's record, which spells its chosen allele, without its scores.
            site.set_call(nullptr, nullptr);
            site.panel_cache.clear();
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

StagedSiteTree::StagedSiteTree(StagedSiteTable& sites, const HandleGraph& graph,
                               const function<string(const SiteBounds&)>& name) {
    // A site named by its bounds is one node, however many staged sites it encloses and whether
    // or not it was staged itself, so it is found by its bounds whichever way round they are read.
    using BoundsKey = tuple<uint64_t, uint64_t, bool>;
    auto key_of = [&](const SiteBounds& bounds) {
        const BoundsKey forward(handlegraph::as_integer(bounds.start),
                                handlegraph::as_integer(bounds.end), bounds.inside_node);
        const BoundsKey turned(handlegraph::as_integer(graph.flip(bounds.end)),
                               handlegraph::as_integer(graph.flip(bounds.start)),
                               bounds.inside_node);
        return min(forward, turned);
    };
    auto ends_of_bounds = [&](const SiteBounds& bounds) {
        return SiteEnds{graph.get_id(bounds.start), graph.get_is_reverse(bounds.start),
                        graph.get_id(bounds.end), graph.get_is_reverse(bounds.end)};
    };
    map<BoundsKey, const Node*> by_bounds;
    auto node_for = [&](const SiteBounds& bounds, const Node* parent) {
        auto found = by_bounds.find(key_of(bounds));
        if (found != by_bounds.end()) {
            return found->second;
        }
        nodes.push_back(Node{ends_of_bounds(bounds), string(), parent});
        by_bounds.emplace(key_of(bounds), &nodes.back());
        return (const Node*)&nodes.back();
    };
    sites.for_each([&](StagedSite& site) {
        // The enclosing sites, outermost first, so that each is added after its own parent.
        const Node* parent = nullptr;
        for (auto it = site.enclosing.rbegin(); it != site.enclosing.rend(); ++it) {
            parent = node_for(*it, parent);
        }
        if (site.id == name(site.bounds)) {
            node_for(site.bounds, parent);
        } else {
            nodes.push_back(Node{ends_of_bounds(site.bounds), site.id, parent});
        }
    });
}

void StagedSiteTree::for_each_site(const function<void(site_t)>& visit, bool in_preorder) const {
    if (in_preorder) {
        for (const Node& node : nodes) {
            visit(&node);
        }
        return;
    }
#pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < nodes.size(); ++i) {
        visit(&nodes[i]);
    }
}

SiteTree::site_t StagedSiteTree::parent_of(site_t site) const {
    return static_cast<const Node*>(site)->parent;
}

SiteEnds StagedSiteTree::ends_of(site_t site) const {
    return static_cast<const Node*>(site)->ends;
}

const string* StagedSiteTree::id_of(site_t site) const {
    const Node* node = static_cast<const Node*>(site);
    return node->id.empty() ? nullptr : &node->id;
}

}
}

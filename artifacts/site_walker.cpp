#include "site_walker.hpp"

#include <algorithm>

#include "decomposition_sites.hpp"
#include "site_scheduler.hpp"

namespace vg::multipass {

SiteWalker::SiteWalker(const SnarlDecomposition& decomposition, const HandleGraph& graph) :
    decomposition(decomposition), graph(graph) {
}

namespace {

/**
 * The forest of the sites in a snarl decomposition.
 *
 * No site is trivial, because `top_level_sites()` and `child_sites()` leave out trivial snarls.
 */
class DecompositionSiteForest : public SiteForest<net_handle_t> {
public:
    /// Show the sites of `decomposition` over `graph`. Neither is owned.
    DecompositionSiteForest(const SnarlDecomposition& decomposition, const HandleGraph& graph) :
        decomposition(decomposition), graph(graph), top_level(top_level_sites(decomposition)) {
    }

    nid_t key(const net_handle_t& site) const override {
        const auto [start, end] = oriented_ends(decomposition, graph, site);
        return std::min(graph.get_id(start), graph.get_id(end));
    }

    void for_each_top(const function<void(const net_handle_t&)>& iteratee) const override {
        for (const net_handle_t& site : top_level) {
            iteratee(site);
        }
    }

    void for_each_top_parallel(const function<void(const net_handle_t&)>& iteratee) const override {
#pragma omp parallel
        {
#pragma omp single
            {
                for (size_t i = 0; i < top_level.size(); i++) {
#pragma omp task firstprivate(i)
                    {
                        iteratee(top_level[i]);
                    }
                }
            }
        }
    }

    bool is_trivial(const net_handle_t&) const override {
        return false;
    }

    void for_each_child(const net_handle_t& site,
                        const function<void(const net_handle_t&)>& iteratee) const override {
        for (const net_handle_t& chain : child_chains(decomposition, site)) {
            for (const net_handle_t& child : child_sites(decomposition, chain)) {
                iteratee(child);
            }
        }
    }

private:
    const SnarlDecomposition& decomposition;
    const HandleGraph& graph;
    /// The top-level sites, in chain order.
    const vector<net_handle_t> top_level;
};

}

void SiteWalker::walk(GraphCaller::RecurseType recurse_type, size_t batch_window,
                      bool show_progress,
                      const function<bool(const net_handle_t&)>& call_site) const {
    walk_site_tree<net_handle_t>(DecompositionSiteForest(decomposition, graph), recurse_type,
                                 batch_window, show_progress, call_site);
}

}

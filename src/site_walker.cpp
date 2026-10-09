#include "site_walker.hpp"

#include <algorithm>

#include "decomposition_sites.hpp"
#include "site_scheduler.hpp"

namespace vg {

SiteWalker::SiteWalker(const SnarlDecomposition& decomposition, const HandleGraph& graph) :
    decomposition(decomposition), graph(graph) {
}

void SiteWalker::walk(GraphCaller::RecurseType recurse_type, size_t batch_window,
                      bool show_progress,
                      const function<bool(const SiteView&)>& call_site) const {
    const vector<net_handle_t> top_level = top_level_sites(decomposition);
    walk_site_tree<net_handle_t>(
        batch_window,
        // The lower of the two boundary node IDs, as for a snarl.
        [&](const net_handle_t& site) {
            const SiteBounds bounds = oriented_bounds(decomposition, graph, site);
            return std::min(graph.get_id(bounds.start), graph.get_id(bounds.end));
        },
        recurse_type == GraphCaller::RecurseAlways, recurse_type == GraphCaller::RecurseOnFail,
        show_progress,
        [&](const function<void(const net_handle_t&)>& lambda) {
            for (const net_handle_t& site : top_level) {
                lambda(site);
            }
        },
        [&](const function<void(const net_handle_t&)>& lambda) {
#pragma omp parallel
            {
#pragma omp single
                {
                    for (size_t i = 0; i < top_level.size(); i++) {
#pragma omp task firstprivate(i)
                        {
                            lambda(top_level[i]);
                        }
                    }
                }
            }
        },
        // The decomposition shows no site that holds nothing to call.
        [](const net_handle_t&) { return false; },
        [&](const net_handle_t& site) {
            return call_site(SiteView{site, oriented_bounds(decomposition, graph, site),
                                      enclosing_sites(decomposition, graph, site)});
        },
        [&](const net_handle_t& site, vector<net_handle_t>& queue) {
            for (const net_handle_t& chain : child_chains(decomposition, site)) {
                for (const net_handle_t& child : child_sites(decomposition, chain)) {
                    queue.push_back(child);
                }
            }
        });
}

}

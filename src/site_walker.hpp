#ifndef VG_SITE_WALKER_HPP_INCLUDED
#define VG_SITE_WALKER_HPP_INCLUDED

#include <functional>

#include <handlegraph/snarl_decomposition.hpp>

#include "graph_caller.hpp"
#include "handle.hpp"
#include "site_values.hpp"

namespace vg {

using namespace std;

/**
 * Walks the sites of a snarl decomposition for the read-likelihood caller, as
 * `GraphCaller::call_top_level_snarls` walks a SnarlManager's snarls: every top-level site in
 * parallel, then the child sites of the sites the recursion asks for, a level at a time (see
 * `walk_site_tree`). It hands each site to its callback as a `SiteView`, with the bounds the
 * decomposition orients it by and the bounds of the sites enclosing it.
 */
class SiteWalker {
public:
    /// Walk the sites of `decomposition` over `graph`. Neither is owned.
    SiteWalker(const SnarlDecomposition& decomposition, const HandleGraph& graph);

    /// Call `call_site` on every top-level site, then on the children of every site
    /// (`RecurseAlways`), of every site `call_site` returns false for (`RecurseOnFail`), or of
    /// none (`RecurseNever`), and so on down. `batch_window` batches the top-level sites by
    /// node-ID window, as `GraphCaller::set_snarl_batching` does; 0 makes each site its own job.
    /// With `show_progress`, reports the sites called to stderr.
    void walk(GraphCaller::RecurseType recurse_type, size_t batch_window, bool show_progress,
              const function<bool(const SiteView&)>& call_site) const;

private:
    const SnarlDecomposition& decomposition;
    const HandleGraph& graph;
};

}

#endif

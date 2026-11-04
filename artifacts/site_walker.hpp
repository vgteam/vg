#ifndef VG_SITE_WALKER_HPP_INCLUDED
#define VG_SITE_WALKER_HPP_INCLUDED

#include <functional>

#include <handlegraph/snarl_decomposition.hpp>

#include "../src/graph_caller.hpp"
#include "../src/handle.hpp"
#include "site_values.hpp"

namespace vg::multipass {

using namespace std;

/**
 * A walker over the sites of a snarl decomposition.
 *
 * It calls every top-level site in parallel, then the child sites that the recursion asks for, a
 * level at a time, as `walk_site_tree` does. It hands each site to its callback as its net handle.
 */
class SiteWalker {
public:
    /// Walk the sites of `decomposition` over `graph`. Neither is owned.
    SiteWalker(const SnarlDecomposition& decomposition, const HandleGraph& graph);

    /**
     * Call `call_site` on every site that `recurse_type` reaches.
     *
     * Start with the top-level sites. Then go on to the children of every site (`RecurseAlways`),
     * of every site that `call_site` returns false for (`RecurseOnFail`), or of none
     * (`RecurseNever`), and so on down.
     *
     * @param batch_window Batch the top-level sites into one parallel job per window of this many
     * node IDs. 0 makes each site its own job.
     * @param show_progress Report the sites called to stderr.
     */
    void walk(GraphCaller::RecurseType recurse_type, size_t batch_window, bool show_progress,
              const function<bool(const net_handle_t&)>& call_site) const;

private:
    const SnarlDecomposition& decomposition;
    const HandleGraph& graph;
};

}

#endif

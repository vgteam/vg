#ifndef VG_SNARL_GRAPH_HPP_INCLUDED
#define VG_SNARL_GRAPH_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "handle.hpp"
#include "snarls.hpp"
#include "traversal_finder.hpp"
#include "snarl_caller.hpp"
#include "region.hpp"
#include "zstdutil.hpp"
#include "vg/io/alignment_emitter.hpp"
#include "gref.hpp"

namespace vg {

using namespace std;
using vg::io::AlignmentEmitter;

/** Simplification of a NetGraph that ignores chains.  It is designed only for
    traversal finding.  Todo: generalize NestedFlowCaller to the point where we 
    can remove this and use NetGraph instead */
class SnarlGraph : virtual public HandleGraph {
public:
    // note: can only deal with one snarl "level" at a time
    SnarlGraph(const HandleGraph* backing_graph, SnarlManager& snarl_manager, vector<const Snarl*> snarls);

    // go from node to snarl (first val false if not a snarl)
    pair<bool, handle_t> node_to_snarl(handle_t handle) const;

    // go from edge to snarl (first val false if not a virtual edge)
    tuple<bool, handle_t, edge_t> edge_to_snarl_edge(edge_t edge) const;

    // replace a snarl node with an actual snarl in the traversal
    void embed_snarl(Visit& visit);
    void embed_snarls(SnarlTraversal& traversal);

    // replace a refpath through the snarl with the actual snarl in the traversal
    // todo: this is a bed of a hack
    void embed_ref_path_snarls(SnarlTraversal& traversal);

    ////////////////////////////////////////////////////////////////////////////
    // Handle-based interface (which is all identical to backing graph)
    ////////////////////////////////////////////////////////////////////////////
    bool has_node(nid_t node_id) const;
    handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const;
    nid_t get_id(const handle_t& handle) const;
    bool get_is_reverse(const handle_t& handle) const;
    handle_t flip(const handle_t& handle) const;
    size_t get_length(const handle_t& handle) const;
    std::string get_sequence(const handle_t& handle) const;    
    size_t get_node_count() const;
    nid_t min_node_id() const;
    nid_t max_node_id() const;
    
protected:

    bool for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
    
    /// this is the only function that's changed to do anything different from the backing graph:
    /// it is changed to "pass through" snarls by pretending there are edges from into snarl starts out of ends and
    /// vice versa.
    bool follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const;    

    /// the backing graph
    const HandleGraph* backing_graph;

    /// the snarl manager
    SnarlManager& snarl_manager;

    /// the snarls (indexed both ways).  flag is true for original orientation
    unordered_map<handle_t, pair<handle_t, bool>> snarls;
};

}

#endif

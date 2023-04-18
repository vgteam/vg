#ifndef VG_CLIP_HPP_INCLUDED
#define VG_CLIP_HPP_INCLUDED

/**
 * \file clip.hpp
 *
 * Clip regions out of a graph
 */

#include "handle.hpp"
#include "snarls.hpp"
#include "region.hpp"

namespace vg {

using namespace std;

/**
 * Visit each snarl if it is fully contained in at least one region from the input set.
 * Only the top-most snarl is visited.
 * The parameters to visit_fn are: 
 *    <the snarl, start_step, end_step, steps_reversed, the containing input region>
 */
void visit_contained_snarls(const PathPositionHandleGraph* graph, const vector<Region>& regions, SnarlManager& snarl_manager,
                            bool include_endpoints,
                            function<void(const Snarl*, step_handle_t, step_handle_t, int64_t, int64_t, bool, const Region*)> visit_fn);

/*
 * Cut nodes out of a graph, and chop up any paths that contain them, using (and resolving) supbath
 * naming conventions from Paths class in path.hpp
 * If a chopped path has a fragment with length < min_fragment_len, don't bother writing the new path
 * The fragments_per_path map is optional, and will collect some stats if present
 */
void delete_nodes_and_chop_paths(MutablePathMutableHandleGraph* graph,
                                 const unordered_set<nid_t>& nodes_to_delete,
                                 const unordered_set<edge_t>& edges_to_delete,
                                 int64_t min_fragment_len,
                                 unordered_map<string, size_t>* fragments_per_path = nullptr);


/**
 * If a given bed region spans a snarl (overlaps its end nodes, and forms a traversal)
 * then clip out all other nodes (ie nodes that don't lie on the traversal)
 *
 * IMPORTANT: for any given snarl, the first region that contains it is used.  
 *
 * Update: now accepts some snarl complexity thresholds to ignore simple enough snarls
 */
void clip_contained_snarls(MutablePathMutableHandleGraph* graph, PathPositionHandleGraph* pp_graph, const vector<Region>& regions,
                           SnarlManager& snarl_manager, bool include_endpoints, int64_t min_fragment_len,
                           size_t max_nodes, size_t max_edges, size_t max_nodes_shallow, size_t max_edges_shallow,
                           double max_avg_degree, double max_reflen_prop, size_t max_reflen,
                           bool out_bed, bool verbose);


/**
 * Clip out nodes that don't pass depth threshold (depth < min_depth).  
 * "depth" is the number of paths that step on the node. 
 * Nodes on path with given prefix ignored (todo: should really switch to regex or something)
 * iterate_handles is a hack to generalize this function to whole graphs or snarls
 */
void clip_low_depth_nodes_and_edges_generic(MutablePathMutableHandleGraph* graph,
                                            function<void(function<void(handle_t, const Region*)>)> iterate_handles,
                                            function<void(function<void(edge_t, const Region*)>)> iterate_edges,
                                            int64_t min_depth, const vector<string>& ref_prefixes,
                                            int64_t min_fragment_len, bool verbose);

/**
 * Run above function on graph
 */
void clip_low_depth_nodes_and_edges(MutablePathMutableHandleGraph* graph, int64_t min_depth, const vector<string>& ref_prefixes,
                          int64_t min_fragment_len, bool verbose);

/**
 * Or on contained snarls
 */
void clip_contained_low_depth_nodes_and_edges(MutablePathMutableHandleGraph* graph, PathPositionHandleGraph* pp_graph, const vector<Region>& regions,
                                    SnarlManager& snarl_manager, bool include_endpoints, int64_t min_depth, int64_t min_fragment_len, bool verbose);

/**
 * clip out deletion edges 
 */
void clip_deletion_edges(MutablePathMutableHandleGraph* graph, int64_t max_deletion, int64_t context_steps,
                         const vector<string>& ref_prefixes, int64_t min_fragment_len, bool verbose);

/**
* clip out stubs
*/
void clip_stubs(MutablePathMutableHandleGraph* graph, const vector<string>& ref_prefixes, int64_t min_fragment_len, bool verbose);

void clip_contained_stubs(MutablePathMutableHandleGraph* graph, PathPositionHandleGraph* pp_graph, const vector<Region>& regions,
                          SnarlManager& snarl_manager, bool include_endpoints, int64_t min_fragment_len, bool verbose);


}


#endif

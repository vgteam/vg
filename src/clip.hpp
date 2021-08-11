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
 *    <the snarl, start_step, end_step, the containing input region>
 */
void visit_contained_snarls(PathPositionHandleGraph* graph, const vector<Region>& regions, SnarlManager& snarl_manager,
                            bool include_endpoints,
                            function<void(const Snarl*, const step_handle_t&, const step_handle_t&, const Region*)> visit_fn);


/**
 * If a given bed region spans a snarl (overlaps its end nodes, and forms a traversal)
 * then clip out all other nodes (ie nodes that don't lie on the traversal)
 *
 * IMPORTANT: for any given snarl, the first region that contains it is used.  
 */
void clip_contained_snarls(MutablePathMutableHandleGraph* graph, PathPositionHandleGraph* pp_graph, const vector<Region>& regions,
                           SnarlManager& snarl_manager, bool include_endpoints, bool verbose);



                                                      
}


#endif

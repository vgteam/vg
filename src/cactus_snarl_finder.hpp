///
///  \file cactus_snarl_finder.hpp
///
///  Defines a widget for finding snarls using the pinchesAndCacti library.
///

#ifndef VG_CACTUS_SNARL_FINDER_HPP_INCLUDED
#define VG_CACTUS_SNARL_FINDER_HPP_INCLUDED

#include "snarls.hpp"

namespace vg {

using namespace std;


/**
 * Class for finding all snarls using the base-level Cactus snarl decomposition
 * interface.
 */
class CactusSnarlFinder : public HandleGraphSnarlFinder {

protected:
    /// Holds the names of reference path hints
    unordered_set<string> hint_paths;
    
    /// Call begin and end functions on this snarl, and all child chains and snarls.
    void recursively_emit_snarls(const Visit& start, const Visit& end,
        const Visit& parent_start, const Visit& parent_end,
        stList* chains_list, stList* unary_snarls_list,
        const function<void(handle_t)>& begin_chain, const function<void(handle_t)>& end_chain,
        const function<void(handle_t)>& begin_snarl, const function<void(handle_t)>& end_snarl);

public:
    /**
     * Make a new CactusSnarlFinder to find snarls in the given graph.
     * We can't filter trivial bubbles because that would break our chains.
     *
     * Skip breaking into connected components if "known_single_component" is true.
     *
     * Optionally takes a hint path name.
     */
    CactusSnarlFinder(const PathHandleGraph& graph, const string& hint_path = "", bool known_single_component = false);
        
    /**
     * Find all the snarls of weakly connected components in parallel.
     * Even single-threaded, this may be worth using as it will use less
     * memory by only considering each component in the context of itself.
     */
    virtual SnarlManager find_snarls_parallel();
    
};

}

#endif

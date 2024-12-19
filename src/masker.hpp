#ifndef VG_MASKER_HPP_INCLUDED
#define VG_MASKER_HPP_INCLUDED

#include <memory>

#include <gbwtgraph/gbwtgraph.h>

#include "handle.hpp"
#include "snarls.hpp"
#include "region.hpp"

namespace vg {

/*
 * Class to mask out graph regions with N's
 */
class Masker {
    
public:
    
    // construct own snarl decomposition
    // use external snarl decomposition
    Masker(gbwtgraph::GBWTGraph& graph, const SnarlManager* snarl_manager = nullptr);
    Masker(MutablePathDeletableHandleGraph& graph, const SnarlManager* snarl_manager = nullptr);
    Masker() = default;
    ~Masker() = default;
    
    // move
    Masker(Masker&& other) = default;
    Masker& operator=(Masker&& other) = default;
    
    // replace regions with N's, with regions given as tuples of (path name, region begin, region end).
    // regions indexes are assumed to be 0-based and end-exclusive
    void mask_sequences(const std::vector<std::tuple<std::string, size_t, size_t>>& regions) const;
    
private:
    
    // make our own snarl manager
    void init_snarls();
    
    // takes tuples of (first step, offset in first step, last step, end in last step)
    // uses the lambda function provided to apply mask sequences to the graph
    void apply_mask_sequences(const std::vector<std::tuple<step_handle_t, size_t, step_handle_t, size_t>>& mask_intervals,
                              const std::function<void(handle_t, const std::string&)>& apply_mask) const;
    
    PathHandleGraph* graph = nullptr;
    const SnarlManager* snarl_manager = nullptr;
    
    // if not given a snarl manager, we construct one and make the snarl_manager
    // point to it
    std::unique_ptr<SnarlManager> own_snarl_manager;
    
};

}

#endif

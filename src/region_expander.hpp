#ifndef VG_REGION_EXPANDER_HPP_INCLUDED
#define VG_REGION_EXPANDER_HPP_INCLUDED

#include "handle.hpp"
#include "snarls.hpp"
#include "gff_reader.hpp"

namespace vg {

    class RegionExpander {
        
    public:
        RegionExpander(PathPositionHandleGraph* xg_index, const SnarlManager* snarl_manager);
        ~RegionExpander() = default;
        
        map<pair<id_t, bool>, pair<uint64_t,uint64_t >> expanded_subgraph(const GFFRecord& gff_record);
        
    private:
        
        PathPositionHandleGraph* xg_index = nullptr;
        const SnarlManager* snarl_manager = nullptr;
        
    };

}

#endif

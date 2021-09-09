#ifndef VG_SNARL_DISTANCE_HPP_INCLUDED
#define VG_SNARL_DISTANCE_HPP_INCLUDED

#include <bdsg/snarl_distance_index.hpp>
#include "snarls.hpp"
#include <structures/union_find.hpp>
#include "hash_map.hpp"


namespace vg { 

using namespace sdsl;
using namespace handlegraph;
using namespace bdsg;
//Fill in the index
void fill_in_distance_index(SnarlDistanceIndex* distance_index, const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit = 500);

//Fill in the temporary snarl record with distances
void populate_snarl_index(SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, 
    pair<SnarlDistanceIndex::temp_record_t, size_t> snarl_index, size_t size_limit, const HandleGraph* graph) ;

SnarlDistanceIndex::TemporaryDistanceIndex make_temporary_distance_index(const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit);

//Define wang_hash for net_handle_t's so that we can use a hash_map
template<> struct wang_hash<handlegraph::net_handle_t> {
public:
    inline size_t operator()(const net_handle_t& net_handle) const {
        return wang_hash_64(reinterpret_cast<const size_t>(as_integer(net_handle)));
    }
};
}
#endif


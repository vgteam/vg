#ifndef VG_SNARL_DISTANCE_HPP_INCLUDED
#define VG_SNARL_DISTANCE_HPP_INCLUDED

//#define debug_indexing

#include <structures/union_find.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "position.hpp"
#include "snarls.hpp"


using namespace sdsl;
using namespace handlegraph;
namespace vg { 

//Fill in the index
make_distance_index(bdsg::SnarlDistanceIndex* distance_index, const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit = 500);

//Fill in the temporary snarl record with distances
void populate_snarl_index(TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, 
    pair<temp_record_t, size_t> snarl_index, size_t size_limit, const HandleGraph* graph) ;

bdsg::SnarlDistanceIndex::TemporaryDistanceIndex make_temporary_distance_index(const HandleGraph* graph, const HandleGraphSnarlFinder* snarl_finder, size_t size_limit);

#endif

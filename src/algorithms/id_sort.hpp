#ifndef VG_ALGORITHMS_ID_SORT_HPP_INCLUDED
#define VG_ALGORITHMS_ID_SORT_HPP_INCLUDED

/**
 * \file id_sort.hpp
 *
 * Defines a by-ID sort algorithm for handle graphs.
 */

#include <unordered_map>

#include "../handle.hpp"


namespace vg {
namespace algorithms {

using namespace std;


/**
 * Order all the handles in the graph in ID order. All orientations are forward.
 */
vector<handle_t> id_order(const HandleGraph* g);
                                                      
}
}

#endif

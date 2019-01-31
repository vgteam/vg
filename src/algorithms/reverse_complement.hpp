#ifndef VG_ALGORITHMS_REVERSE_COMPLEMENT_HPP_INCLUDED
#define VG_ALGORITHMS_REVERSE_COMPLEMENT_HPP_INCLUDED

/**
 * \file reverse_complement.hpp
 *
 * Defines algorithm for reverse complementing the sequence in a graph
 */

#include "../handle.hpp"
#include "../utility.hpp"

namespace vg {
namespace algorithms {

using namespace std;

    /// Fills a MutableHandleGraph 'into' with a graph that has the same sequence and path
    /// space as 'source', but the forward strand of every node is flipped to the reverse
    /// strand. Reports an error and exits if 'into' is not empty.
    unordered_map<id_t, pair<id_t, bool>> reverse_complement_graph(const HandleGraph* source,
                                                                   MutableHandleGraph* into);

}
}

#endif

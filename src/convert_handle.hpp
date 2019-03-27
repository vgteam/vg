//
//  convert_handle.hpp
//

#ifndef VG_CONVERT_HANDLE_HPP_INCLUDED
#define VG_CONVERT_HANDLE_HPP_INCLUDED

#include "handle.hpp"
#include "vg.hpp"

namespace vg {
    using namespace std;
    // Takes in a pointer to a HandleGraph and converts it to a MutableHandleGraph graph.
    void convert_handle_graph(const HandleGraph* converting, MutableHandleGraph* converted);
    
    // Change paths to a mutable path.
    void convert_path_handle_graph(const PathHandleGraph* converting, MutablePathDeletableHandleGraph* converted);

}

#endif

//
//  convert_handle.hpp
//

#ifndef VG_CONVERT_HANDLE_HPP_INCLUDED
#define VG_CONVERT_HANDLE_HPP_INCLUDED

#include "handle.hpp"
#include "vg.hpp"

namespace vg {
    using namespace std;
    /// Takes in a pointer to a HandleGraph and copies it into an empty MutableHandleGraph graph.
    void convert_handle_graph(const HandleGraph* converting, MutableHandleGraph* converted);
    
    /// Takes in a pointer to a PathHandleGraph and copies it, inclding paths, into an empty MutablePathMutableHandleGraph graph.
    void convert_path_handle_graph(const PathHandleGraph* converting, MutablePathMutableHandleGraph* converted);

}

#endif

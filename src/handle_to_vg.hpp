#ifndef VG_HANDLE_TO_VG_HPP_INCLUDED
#define VG_HANDLE_TO_VG_HPP_INCLUDED

#include "handle.hpp"
#include "vg.hpp"

namespace vg {
	using namespace std;
	/// Takes in a pointer to a HandleGraph and converts it to a VG graph.
	VG handle_to_vg(const HandleGraph* g);
}

#endif

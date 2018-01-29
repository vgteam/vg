#ifndef VG_CONVERTER_HPP_INCLUDED
#define VG_CONVERTER_HPP_INCLUDED

#include "handle.hpp"
#include "xg.hpp"
#include "vg.hpp"

namespace vg {
	using namespace std;
	/// Takes in a pointer to a HandleGraph and converts it to a VG graph.
	VG converter(const HandleGraph* g);
}

#endif
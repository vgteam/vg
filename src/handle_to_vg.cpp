#include "handle_to_vg.hpp"

namespace vg {
	using namespace std;
	
	VG handle_to_vg(const HandleGraph* xg) {
		// If xg is a null pointer, throw a runtime error
		if (xg == nullptr) {
			throw runtime_error("There is no xg to convert"); 
		} 
		// Initialize the VG graph
		VG vg;
		// Iterate through each handle in xg and create the same handle in vg
		xg->for_each_handle([&](const handle_t& here) {
			// Get the id of the xg handle
			id_t xg_id = xg->get_id(here);
			// Get the sequence of the xg handle
			string xg_seq = xg->get_sequence(here);
			// Create a handle in vg using the xg id and sequence
			vg.create_handle(xg_seq,xg_id);
		});
		// Iterate through each handle in xg 
		xg->for_each_handle([&](const handle_t& handle) {
			id_t id = xg->get_id(handle);
			bool rev = xg->get_is_reverse(handle);
			// Return a vg handle using the xg handle's id and orientation 
			handle_t current = vg.get_handle(id,rev);
			// Follow the right edges of the xg handle
			xg->follow_edges(handle, false, [&](const handle_t& r) {
				id_t id_r = xg->get_id(r);
				bool rev_r = xg->get_is_reverse(r);
				// Return a vg handle using the xg handle's id and orientation
				handle_t next = vg.get_handle(id_r, rev_r);
				// Create an edge in vg using the handles 
				vg.create_edge(current,next);
			});
			// Follow the left edges of the xg handle
			xg->follow_edges(handle, true, [&](const handle_t& l) {
				id_t id_l = xg->get_id(l);
				bool rev_l = xg->get_is_reverse(l);
				// Return a vg handle using the xg handle's id and orientation
				handle_t prev = vg.get_handle(id_l, rev_l);
				// Use the handles created from following the xg edges to create a vg edge
				vg.create_edge(prev,current); //error here
			});
		});
		return vg;
	}

}

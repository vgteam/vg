#ifndef VG_IO_SAVE_HANDLE_GRAPH_IO_HPP_INCLUDED
#define VG_IO_REGISTER_LIBVG_IO_HPP_INCLUDED

/**
 * \save_handle_graph.hpp
 * Use vpkg to serialize a HandleGraph object
 */

#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"
#include "vg.hpp"
#include "xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

namespace vg {

namespace io {

using namespace std;


/**
 * Save a handle graph using the VPKG::save() function. 
 * Todo: should this be somewhere else (ie in vgio with new types registered?)
 */
inline void save_handle_graph(HandleGraph* graph, ostream& os) {
    if (dynamic_cast<VG*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<VG*>(graph), os);
    } else if (dynamic_cast<bdsg::HashGraph*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<bdsg::HashGraph*>(graph), os);
    } else if (dynamic_cast<bdsg::PackedGraph*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<bdsg::PackedGraph*>(graph), os);
    } else if (dynamic_cast<bdsg::ODGI*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<bdsg::ODGI*>(graph), os);
    } else if (dynamic_cast<xg::XG*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<xg::XG*>(graph), os);
    } else {
        throw runtime_error("Internal error: unable to serialize graph");
    }
}

inline void save_handle_graph(HandleGraph* graph, const string& dest_path) {
    if (dynamic_cast<VG*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<VG*>(graph), dest_path);
    } else if (dynamic_cast<bdsg::HashGraph*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<bdsg::HashGraph*>(graph), dest_path);
    } else if (dynamic_cast<bdsg::PackedGraph*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<bdsg::PackedGraph*>(graph), dest_path);
    } else if (dynamic_cast<bdsg::ODGI*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<bdsg::ODGI*>(graph), dest_path);
    } else if (dynamic_cast<xg::XG*>(graph) != nullptr) {
        vg::io::VPKG::save(*dynamic_cast<xg::XG*>(graph), dest_path);
    } else {
        throw runtime_error("Internal error: unable to serialize graph");
    }    
}
    
}

}

#endif

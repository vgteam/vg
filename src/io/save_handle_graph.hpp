#ifndef VG_IO_SAVE_HANDLE_GRAPH_HPP_INCLUDED
#define VG_IO_SAVE_HANDLE_GRAPH_HPP_INCLUDED

/**
 * \file save_handle_graph.hpp
 * Use vpkg to serialize a HandleGraph object
 */

#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"
#include "vg.hpp"
#include "xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

#include <memory>

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

// Check that output format specifier is a valid graph type
inline bool valid_output_format(const string& fmt_string) {
    return fmt_string == "vg" || fmt_string == "pg" || fmt_string == "hg";
}

// Create a new graph (of handle graph type T) where the implementation is chosen using the format string
template<class T>
unique_ptr<T> new_output_graph(const string& fmt_string) {
    if (fmt_string == "vg") {
        return make_unique<VG>();
    } else if (fmt_string == "pg") {
        return make_unique<bdsg::PackedGraph>();
    } else if (fmt_string == "hg") {
        return make_unique<bdsg::HashGraph>();
    } else {
        return unique_ptr<T>();
    }
}

}

}

#endif

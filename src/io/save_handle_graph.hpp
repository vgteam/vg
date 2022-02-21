#ifndef VG_IO_SAVE_HANDLE_GRAPH_HPP_INCLUDED
#define VG_IO_SAVE_HANDLE_GRAPH_HPP_INCLUDED

/**
 * \file save_handle_graph.hpp
 * Use vpkg to serialize a HandleGraph object
 */

#include <handlegraph/serializable_handle_graph.hpp>
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"
#include "vg.hpp"
#include "xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "gfa.hpp"
#include <memory>

// TODO: move GFAIDMapInfo out of here
#include "../algorithms/gfa_to_handle.hpp"

namespace vg {

/// Use to load in GFAs and remember where they came from
class GFAHandleGraph : public bdsg::PackedGraph {
public:
   GFAHandleGraph() : bdsg::PackedGraph() {}
   virtual ~GFAHandleGraph() = default;
   
   /// Store the translation from graph ID space (as initially read in) back to
   /// GFA ID space.
   /// Won't be useful if the graph is modified.
   vg::algorithms::GFAIDMapInfo gfa_id_space;
};

namespace io {

using namespace std;


/**
 * Save a handle graph. 
 * Todo: should this be somewhere else (ie in vgio with new types registered?)
 */
inline void save_handle_graph(HandleGraph* graph, ostream& os) {

    if (dynamic_cast<GFAHandleGraph*>(graph) != nullptr) {
        // We loaded a GFA into a handle graph, want to write back to GFA
        graph_to_gfa(dynamic_cast<GFAHandleGraph*>(graph), os);
    } else if (dynamic_cast<SerializableHandleGraph*>(graph) != nullptr) {
        // SerializableHandleGraphs are all serialized bare, without VPKG framing, for libbdsg compatibility.
        dynamic_cast<SerializableHandleGraph*>(graph)->serialize(os);
    } else if (dynamic_cast<VG*>(graph) != nullptr) {
        // vg::VG doesn't use a magic number and isn't a SerializableHandleGraph
        vg::io::VPKG::save(*dynamic_cast<VG*>(graph), os);
    } else {
        throw runtime_error("Internal error: unable to serialize graph");
    }
}

inline void save_handle_graph(HandleGraph* graph, const string& dest_path) {
    if (dynamic_cast<GFAHandleGraph*>(graph) != nullptr) {
        // We loaded a GFA into a handle graph, want to write back to GFA
        ofstream os(dest_path);
        if (!os) {
            throw runtime_error("error[save_handle_graph]: Unable to write to: " + dest_path);
        }
        graph_to_gfa(dynamic_cast<GFAHandleGraph*>(graph), os);
    } else if (dynamic_cast<SerializableHandleGraph*>(graph) != nullptr) {
        // SerializableHandleGraphs are all serialized bare, without VPKG framing, for libbdsg compatibility.
        dynamic_cast<SerializableHandleGraph*>(graph)->serialize(dest_path);
    } else if (dynamic_cast<VG*>(graph) != nullptr) {
        // vg::VG doesn't use a magic number and isn't a SerializableHandleGraph
        vg::io::VPKG::save(*dynamic_cast<VG*>(graph), dest_path);
    } else {
        throw runtime_error("Internal error: unable to serialize graph");
    }    
}

// Check that output format specifier is a valid graph type
inline bool valid_output_format(const string& fmt_string) {
    return fmt_string == "vg" || fmt_string == "pg" || fmt_string == "hg" || fmt_string == "gfa";
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
    } else if (fmt_string == "gfa") {
        return make_unique<GFAHandleGraph>();
    } else {
        return unique_ptr<T>();
    }
}

}

}

#endif

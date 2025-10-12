#ifndef VG_GBWTGRAPH_HELPER_HPP_INCLUDED
#define VG_GBWTGRAPH_HELPER_HPP_INCLUDED

/** \file 
 * Utility classes and functions for working with GBWTGraph.
 */

#include <gbwtgraph/gfa.h>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>
#include "position.hpp"
#include <unordered_map>
#include <vector>

namespace vg {

//------------------------------------------------------------------------------

/**
 * Get the best configuration to use for the GBWTGraph library GFA parser, to
 * best matcch the behavior of vg's GFA parser.
 */
gbwtgraph::GFAParsingParameters get_best_gbwtgraph_gfa_parsing_parameters();

/*
    These are the proper ways of saving and loading GBWTGraph structures. In
    case of a failure, these will print an error message and exit the program.

    The vg::io::VPKG interface is effectively the same, but it does not handle
    errors in a consistent way.
*/

/// Load GBWTGraph from the file.
/// NOTE: Call `graph.set_gbwt()` afterwards with the appropriate GBWT index.
void load_gbwtgraph(gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Load GBZ from the file.
void load_gbz(gbwtgraph::GBZ& gbz, const std::string& filename, bool show_progress = false);

/// Load GBZ from separate GBWT / GBWTGraph files.
void load_gbz(gbwtgraph::GBZ& gbz, const std::string& gbwt_name, const std::string& graph_name, bool show_progress = false);

/// Load GBWT and GBWTGraph from the GBZ file.
void load_gbz(gbwt::GBWT& index, gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Load a minimizer index from the file.
void load_minimizer(gbwtgraph::DefaultMinimizerIndex& index, const std::string& filename, bool show_progress = false);

/// Load a minimizerXL index from the file.
void load_minimizer(gbwtgraph::MinimizerIndexXL& index, const std::string& filename, bool show_progress = false);


/// Save GBWTGraph to the file.
void save_gbwtgraph(const gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Save GBZ to the file.
void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& filename, bool show_progress = false);

/// Save GBWT and GBWTGraph to the GBZ file.
void save_gbz(const gbwt::GBWT& index, gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Save GBZ to separate GBWT / GBWTGraph files.
void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& gbwt_name, const std::string& graph_name, bool show_progress = false);

/// Save a minimizer index to the file.
template<typename IndexType>
void save_minimizer(const IndexType& index, const std::string& filename, bool show_progress = false) {
    if (show_progress) {
        std::cerr << "Saving MinimizerIndex to " << filename << std::endl;
    }

    try {
        std::ofstream out(filename, std::ios_base::binary);
        if (!out) {
            throw sdsl::simple_sds::CannotOpenFile(filename, true);
        }
        out.exceptions(std::ofstream::badbit | std::ofstream::failbit);
        index.serialize(out);
        out.close();
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [save_minimizer()] cannot save MinimizerIndex to " << filename << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

//------------------------------------------------------------------------------

/// Return a mapping of the original segment ids to a list of chopped node ids
std::unordered_map<std::string, std::vector<nid_t>> load_translation_map(const gbwtgraph::GBWTGraph& graph);

/// Return a backwards mapping of chopped node to original segment position (id,offset pair)
std::unordered_map<nid_t, std::pair<std::string, size_t>> load_translation_back_map(const gbwtgraph::GBWTGraph& graph);

//------------------------------------------------------------------------------

/// Returns an empty GBWTGraph handle corresponding to the GBWT endmarker.
inline handle_t empty_gbwtgraph_handle() {
    return gbwtgraph::GBWTGraph::node_to_handle(0);
}

/// Returns a string representation of a GBWTGraph handle.
std::string to_string_gbwtgraph(handle_t handle);

/// Returns a string representation of a GBWTGraph node.
std::string to_string_gbwtgraph(gbwt::node_type node);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GBWTGRAPH_HELPER_HPP_INCLUDED

#ifndef VG_GBWTGRAPH_HELPER_HPP_INCLUDED
#define VG_GBWTGRAPH_HELPER_HPP_INCLUDED

/** \file 
 * Utility classes and functions for working with GBWTGraph.
 */

#include <gbwtgraph/gfa.h>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>

#include "position.hpp"
#include "snarl_distance_index.hpp"
#include "zip_code.hpp"

#include <vector>
#include <unordered_map>

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
/// require_payload() can be used afterwards to ensure the payload type is correct.
void load_minimizer(gbwtgraph::DefaultMinimizerIndex& index, const std::string& filename, bool show_progress = false);

/// Save GBWTGraph to the file.
void save_gbwtgraph(const gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Save GBZ to the file.
void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& filename, bool show_progress = false);

/// Save GBWT and GBWTGraph to the GBZ file.
/// NOTE: GBZ tags will be empty, apart from the source tag.
void save_gbz(const gbwt::GBWT& index, gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Save GBZ to separate GBWT / GBWTGraph files.
void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& gbwt_name, const std::string& graph_name, bool show_progress = false);

/// Save a minimizer index to the file.
void save_minimizer(const gbwtgraph::DefaultMinimizerIndex& index, const std::string& filename, bool show_progress = false);

//------------------------------------------------------------------------------

enum GraphCompatibilityFlags {
    GRAPH_COMPATIBILITY_DEFAULT = 0x00,
    // Require both graphs to have names.
    GRAPH_COMPATIBILITY_STRICT = 0x01,
    // Allow the first graph to be a subgraph of the second graph.
    GRAPH_COMPATIBILITY_SUBGRAPH = 0x02,
    // Print extended information about the relationship between the graphs.
    GRAPH_COMPATIBILITY_VERBOSE = 0x04
};

/// Implementation of require_compatible_graphs().
void require_compatible_graphs_impl(
    const gbwtgraph::GraphName& first_name, const std::string& first_decription,
    const gbwtgraph::GraphName& second_name, const std::string& second_description,
    GraphCompatibilityFlags flags
);

/**
 * Checks that the objects are compatible with each other, according to the
 * gbwtgraph::GraphName information possibly stored in the tags. Prints an error
 * message and exits on failure.
 *
 * Both objects must have a graph_name() method returning gbwtgraph::GraphName.
 * If either object does not have a name, the check will succeed, unless strict
 * mode is enabled. By default, the names should be the same (the corresponding
 * graphs are identical). If subgraph is true, the first object may be for a
 * subgraph of the second object. If verbose is true, extended information about
 * the possible relationship between the graphs is printed.
 */
template <class T1, class T2>
void require_compatible_graphs(
    const T1& first, const std::string& first_decription,
    const T2& second, const std::string& second_description,
    GraphCompatibilityFlags flags = GRAPH_COMPATIBILITY_DEFAULT
) {
    gbwtgraph::GraphName first_name = first.graph_name();
    gbwtgraph::GraphName second_name = second.graph_name();
    require_compatible_graphs_impl(first_name, first_decription, second_name, second_description, flags);
}

//------------------------------------------------------------------------------

/// Minimizer index construction parameters.
struct MinimizerIndexParameters {
    /// Default for `threshold`. Should be the same as Giraffe hard hit cap.
    constexpr static size_t DEFAULT_THRESHOLD = 500;

    /// Default for `iterations`.
    constexpr static size_t DEFAULT_ITERATIONS = 3;

    /// Maximum allowed value for `iterations`.
    constexpr static size_t MAX_ITERATIONS = gbwtgraph::MinimizerHeader::FLAG_WEIGHT_MASK >> gbwtgraph::MinimizerHeader::FLAG_WEIGHT_OFFSET;

    /// Lower bound for `hash_table_width`.
    constexpr static size_t HASH_TABLE_MIN_WIDTH = 10;

    /// Upper bound for `hash_table_width`.
    constexpr static size_t HASH_TABLE_MAX_WIDTH = 36;

    /// Number of words used for a zipcode payload.
    constexpr static size_t ZIPCODE_PAYLOAD_SIZE = sizeof(ZipCode::payload_type) / sizeof(gbwtgraph::KmerEncoding::code_type); 

    enum PayloadType {
        /// No payload.
        PAYLOAD_NONE,
        /// Zipcode payload.
        PAYLOAD_ZIPCODES,
        /// Zipcode payload with path information.
        PAYLOAD_ZIPCODES_WITH_PATHS
    };

    /// Tag used to indicate the payload type in the minimizer index.
    const static std::string PAYLOAD_KEY; // "payload"

    /// k-mer length.
    size_t k = gbwtgraph::DefaultMinimizerIndex::key_type::KMER_LENGTH;

    /// Window length (with minimizers) or s-mer length (with syncmers).
    size_t w_or_s = gbwtgraph::DefaultMinimizerIndex::key_type::WINDOW_LENGTH;

    /// Whether to use syncmers instead of minimizers.
    bool use_syncmers = false;

    /// Whether to include path information in the payload (for recombination-aware mapping).
    /// Ignored if there is no zipcode payload.
    bool paths_in_payload = false;

    /// Whether to use weighted minimizers (cannot be used with syncmers).
    bool use_weighted_minimizers = false;

    /// Downweight kmers with more than this many hits.
    size_t threshold = DEFAULT_THRESHOLD;

    /// Number of iterations for weighted minimizer calculation.
    size_t iterations = DEFAULT_ITERATIONS;

    /// Whether to use the space-efficient k-mer counting algorithm with weighted minimizers.
    bool space_efficient_counting = false;

    /// Initial hash table width (in bits) for kmer counting (0 = guess).
    size_t hash_table_width = 0;

    /// Print progress information during construction.
    bool progress = false;

    /// Sets minimizer parameters.
    MinimizerIndexParameters& minimizers(size_t k, size_t w) {
        this->k = k;
        this->w_or_s = w;
        this->use_syncmers = false;
        return *this;
    }

    /// Sets syncmer parameters.
    MinimizerIndexParameters& syncmers(size_t k, size_t s) {
        this->k = k;
        this->w_or_s = s;
        this->use_syncmers = true;
        return *this;
    }

    /// Includes path information in the payload.
    MinimizerIndexParameters& with_paths(bool paths_in_payload = true) {
        this->paths_in_payload = paths_in_payload;
        return *this;
    }

    /// Sets weighted minimizer parameters.
    MinimizerIndexParameters& weighted(
        bool use_weighted_minimizers,
        size_t threshold = DEFAULT_THRESHOLD, size_t iterations = DEFAULT_ITERATIONS
    ) {
        this->use_weighted_minimizers = use_weighted_minimizers;
        this->threshold = threshold;
        this->iterations = iterations;
        return *this;
    }

    /// Sets k-mer counting parameters.
    MinimizerIndexParameters& kmer_counting(bool space_efficient_counting = false, size_t hash_table_width = 0) {
        this->space_efficient_counting = space_efficient_counting;
        this->hash_table_width = hash_table_width;
        return *this;
    }

    /// Sets progress printing.
    MinimizerIndexParameters& verbose(bool progress = true) {
        this->progress = progress;
        return *this;
    }

    /// Returns an error message if the parameters are invalid.
    std::string validate() const;

    /// Returns a string representation of the payload type.
    /// This will be used as a tag value in the minimizer index.
    static std::string payload_str(PayloadType type);
};

/// Builds a new minimizer index. If a distance index is provided, zipcodes are
/// stored as payload. If a zipcode collection is also provided, zipcodes that
/// do not fit in the payload are stored there.
///
/// Prints an error message and exits on failure. Uses OpenMP for multithreading.
gbwtgraph::DefaultMinimizerIndex build_minimizer_index(
    const gbwtgraph::GBZ& gbz,
    const SnarlDistanceIndex* distance_index,
    ZipCodeCollection* oversized_zipcodes,
    const MinimizerIndexParameters& params
);

/// Checks that the minimizer index has the expected payload type.
///
/// The check is based on tag "payload" stored in the index.
/// Prints an error message and exits on failure.
void require_payload(const gbwtgraph::DefaultMinimizerIndex& index, MinimizerIndexParameters::PayloadType expected_payload);

/// Checks that the minimizer index has one of the expected payload types.
///
/// The check is based on tag "payload" stored in the index.
/// Prints an error message and exits on failure.
void require_payload(
    const gbwtgraph::DefaultMinimizerIndex& index,
    const std::vector<MinimizerIndexParameters::PayloadType>& expected_payloads
);

/// Returns true if the minimizer index has the given payload, and false otherwise.
///
/// The check is based on tag "payload" stored in the index.
bool has_payload(const gbwtgraph::DefaultMinimizerIndex& index, MinimizerIndexParameters::PayloadType payload);

/// Returns true if the minimizer index has any of the given payloads, and false otherwise.
///
/// The check is based on tag "payload" stored in the index.
bool has_payload(
    const gbwtgraph::DefaultMinimizerIndex& index,
    const std::vector<MinimizerIndexParameters::PayloadType>& payloads
);

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

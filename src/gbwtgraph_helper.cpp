#include "gbwtgraph_helper.hpp"
#include "gbwt_helper.hpp"
#include "gbzgraph.hpp"

#include <gbwtgraph/index.h>
#include <vg/io/alignment_io.hpp>

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t MinimizerIndexParameters::DEFAULT_THRESHOLD;
constexpr size_t MinimizerIndexParameters::DEFAULT_ITERATIONS;
constexpr size_t MinimizerIndexParameters::MAX_ITERATIONS;
constexpr size_t MinimizerIndexParameters::HASH_TABLE_MIN_WIDTH;
constexpr size_t MinimizerIndexParameters::HASH_TABLE_MAX_WIDTH;
constexpr size_t MinimizerIndexParameters::ZIPCODE_PAYLOAD_SIZE;

// Other static members.
const std::string MinimizerIndexParameters::PAYLOAD_KEY = "payload";

//------------------------------------------------------------------------------

gbwtgraph::GFAParsingParameters get_best_gbwtgraph_gfa_parsing_parameters() {
    gbwtgraph::GFAParsingParameters parameters;
    // Configure GBWTGraph GFA parsing to be as close to the vg GFA parser as we can get.
    // TODO: Make it closer.
    parameters.path_name_formats.clear();

    // Parse panSN with a fragment after it (e.g. HG002#1#chr3#1).
    parameters.path_name_formats.emplace_back(
        gbwtgraph::GFAParsingParameters::PAN_SN_REGEX + "#([0-9][0-9]*)",
        gbwtgraph::GFAParsingParameters::PAN_SN_FIELDS + "F",
        gbwtgraph::GFAParsingParameters::PAN_SN_SENSE
    );

    // Parse panSN with a range after it as a normal but with a fragment based
    // on start position (e.g. HG002#1#chr3[1566235] or HG002#1#chr3[1566235-2397571]).
    parameters.path_name_formats.emplace_back(
        gbwtgraph::GFAParsingParameters::PAN_SN_REGEX + "\\[([0-9][0-9]*)(-[0-9]*)?\\]",
        gbwtgraph::GFAParsingParameters::PAN_SN_FIELDS + "F",
        gbwtgraph::GFAParsingParameters::PAN_SN_SENSE
    );

    // Parse standard panSN as what we think that is (e.g. HG002#1#chr3).
    parameters.path_name_formats.emplace_back(
        gbwtgraph::GFAParsingParameters::PAN_SN_REGEX,
        gbwtgraph::GFAParsingParameters::PAN_SN_FIELDS,
        gbwtgraph::GFAParsingParameters::PAN_SN_SENSE
    );

    // Parse path names with a sample and a contig (e.g. GRCh38#chr3).
    parameters.path_name_formats.emplace_back(
        "(.*)#(.*)",
        "XSC",
        PathSense::HAPLOTYPE
    );

    // Parse paths with just a name and a range as generic paths with a contig
    // and a fragment. Sample for generic paths gets provided automatically.
    parameters.path_name_formats.emplace_back(
        "(.*)\\[([0-9][0-9]*)(-[0-9]*)?\\]",
        "XCF",
        PathSense::GENERIC
    );

    // Parse paths with nothing to distinguish them the default way (as generic named paths)
    parameters.path_name_formats.emplace_back(
        gbwtgraph::GFAParsingParameters::DEFAULT_REGEX,
        gbwtgraph::GFAParsingParameters::DEFAULT_FIELDS,
        gbwtgraph::GFAParsingParameters::DEFAULT_SENSE
    );

    return parameters;
}

void load_gbwtgraph(gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading GBWTGraph from " << filename << std::endl;
    }

    // This mimics Simple-SDS serialization.
    try {
        std::ifstream in(filename, std::ios_base::binary);
        if (!in) {
            throw sdsl::simple_sds::CannotOpenFile(filename, false);
        }
        in.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
        graph.deserialize(in);
        in.close();
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [load_gbwtgraph()] cannot load GBWTGraph " << filename << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void load_gbz(gbwtgraph::GBZ& gbz, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading GBZ from " << filename << std::endl;
    }
    try {
        sdsl::simple_sds::load_from(gbz, filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [load_gbz()] cannot load GBZ " << filename << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void load_gbz(gbwt::GBWT& index, gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading GBWT and GBWTGraph from " << filename << std::endl;
    }
    gbwtgraph::GBZ loaded;
    try {
        sdsl::simple_sds::load_from(loaded, filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [load_gbz()] cannot load GBZ " << filename << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(loaded.index);
    graph = std::move(loaded.graph);
    graph.set_gbwt_address(index); // We moved the GBWT out from the GBZ, so we have to update the pointer.
}

void load_gbz(gbwtgraph::GBZ& gbz, const std::string& gbwt_name, const std::string& graph_name, bool show_progress) {
    gbz = gbwtgraph::GBZ();
    load_gbwt(gbz.index, gbwt_name, show_progress);
    load_gbwtgraph(gbz.graph, graph_name, show_progress);
    gbz.graph.set_gbwt(gbz.index); // We know the GBWT index corresponding to the graph.
}

void load_minimizer(gbwtgraph::DefaultMinimizerIndex& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading MinimizerIndex from " << filename << std::endl;
    }

    // This mimics Simple-SDS serialization.
    try {
        std::ifstream in(filename, std::ios_base::binary);
        if (!in) {
            throw sdsl::simple_sds::CannotOpenFile(filename, false);
        }
        in.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
        index.deserialize(in);
        in.close();
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [load_minimizer()] cannot load MinimizerIndex " << filename << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void save_gbwtgraph(const gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving GBWTGraph to " << filename << std::endl;
    }

    // This mimics Simple-SDS serialization.
    try {
        std::ofstream out(filename, std::ios_base::binary);
        if (!out) {
            throw sdsl::simple_sds::CannotOpenFile(filename, true);
        }
        out.exceptions(std::ofstream::badbit | std::ofstream::failbit);
        graph.serialize(out);
        out.close();
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [save_gbwtgraph()] cannot save GBWTGraph to " << filename << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving GBZ to " << filename << std::endl;
    }
    try {
        sdsl::simple_sds::serialize_to(gbz, filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [save_gbz()] cannot save GBZ to " << filename << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void save_gbz(const gbwt::GBWT& index, gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving GBWT and GBWTGraph to " << filename << std::endl;
    }

    // This mimics Simple-SDS serialization.
    try {
        std::ofstream out(filename, std::ios_base::binary);
        if (!out) {
            throw sdsl::simple_sds::CannotOpenFile(filename, true);
        }
        out.exceptions(std::ofstream::badbit | std::ofstream::failbit);
        gbwtgraph::GBZ::simple_sds_serialize(index, graph, out);
        out.close();
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [save_gbz()] cannot save GBWT and GBWTGraph to " << filename << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& gbwt_name, const std::string& graph_name, bool show_progress) {
    save_gbwt(gbz.index, gbwt_name, show_progress);
    save_gbwtgraph(gbz.graph, graph_name, show_progress);
}

void save_minimizer(const gbwtgraph::DefaultMinimizerIndex& index, const std::string& filename, bool show_progress) {
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

GraphCompatibilityFlags operator|(GraphCompatibilityFlags a, GraphCompatibilityFlags b) {
    return static_cast<GraphCompatibilityFlags>(static_cast<int>(a) | static_cast<int>(b));
}

GraphCompatibilityFlags& operator|=(GraphCompatibilityFlags& a, GraphCompatibilityFlags b) {
    a = a | b;
    return a;
}

void require_compatible_graphs_impl(
    const gbwtgraph::GraphName& first_name, const std::string& first_decription,
    const gbwtgraph::GraphName& second_name, const std::string& second_description,
    GraphCompatibilityFlags flags
) {
    if (!(flags & GRAPH_COMPATIBILITY_STRICT)) {
        if (!first_name.has_name() || !second_name.has_name()) {
            return;
        }
    }

    if (first_name.same(second_name)) {
        return;
    }
    if (flags & GRAPH_COMPATIBILITY_SUBGRAPH) {
        if (first_name.subgraph_of(second_name)) {
            return;
        }
    }

    std::cerr << "error: \"" << first_decription << "\" and \"" << second_description << "\" are not compatible" << std::endl;
    std::string relationship = first_name.describe_relationship(second_name, first_decription, second_description);
    std::cerr << relationship << std::endl;
    std::exit(EXIT_FAILURE);
}

void require_compatible_reference(
    const std::string& gaf_filename,
    const HandleGraph* handle_graph, const gbwtgraph::GBZ* gbz,
    bool strict
) {
    if (handle_graph == nullptr && gbz == nullptr) {
        return;
    }

    gbwtgraph::GraphName graph_name;
    if (gbz != nullptr) {
        graph_name = gbz->graph_name();
    } else {
        // GBZGraph is the only HandleGraph implementation that currently supports graph names.
        const GBZGraph* gbz_graph = dynamic_cast<const GBZGraph*>(handle_graph);
        if (gbz_graph == nullptr) {
            return;
        }
        graph_name = gbz_graph->gbz.graph_name();
    }

    // This will exit if the file does not exist and do nothing if the file is stdin ("-").
    std::vector<std::string> gaf_header_lines = vg::io::read_gaf_header_lines(gaf_filename);
    gbwtgraph::GraphName gaf_name(gaf_header_lines);

    GraphCompatibilityFlags flags = GRAPH_COMPATIBILITY_SUBGRAPH;
    if (strict) {
        flags |= GRAPH_COMPATIBILITY_STRICT;
    }
    require_compatible_graphs_impl(gaf_name, "GAF file", graph_name, "reference graph", flags);
}

//------------------------------------------------------------------------------

std::string MinimizerIndexParameters::validate() const {
    if (this->k < 1 || this->k > gbwtgraph::DefaultMinimizerIndex::key_type::KMER_MAX_LENGTH) {
        return "k-mer length must be between 1 and " + std::to_string(gbwtgraph::DefaultMinimizerIndex::key_type::KMER_MAX_LENGTH);
    }

    if (this->use_syncmers) {
        if (this->w_or_s < 1 || this->w_or_s > this->k) {
            return "s-mer length must be between 1 and k-mer length when using syncmers";
        }
        if (this->use_weighted_minimizers) {
            return "weighted minimizers cannot be used with syncmers";
        }
    } else {
        if (this->w_or_s < 1) {
            return "window length must be at least 1 when using minimizers";
        }
    }

    if (this->use_weighted_minimizers) {
        if (this->iterations > MAX_ITERATIONS) {
            return "number of iterations must be at most " + std::to_string(MAX_ITERATIONS) + " when using weighted minimizers";
        }
        if (this->hash_table_width != 0 &&
            (this->hash_table_width < HASH_TABLE_MIN_WIDTH || this->hash_table_width > HASH_TABLE_MAX_WIDTH)) {
            return "hash table width must be between " + std::to_string(HASH_TABLE_MIN_WIDTH) +
                   " and " + std::to_string(HASH_TABLE_MAX_WIDTH) + " when using weighted minimizers";
        }
    }

    return "";
}

std::string MinimizerIndexParameters::payload_str(PayloadType type) {
    switch (type) {
        case PAYLOAD_NONE:
            return "none";
        case PAYLOAD_ZIPCODES:
            return "zipcodes";
        case PAYLOAD_ZIPCODES_WITH_PATHS:
            return "zipcodes with paths";
        default:
            return "unknown";
    }
}

size_t trailing_zeros(size_t value) {
    size_t result = 0;
    if (value == 0) {
        return result;
    }
    while ((value & 1) == 0) {
        value >>= 1;
        result++;
    }
    return result;
}

size_t estimate_hash_table_size(const gbwtgraph::GBZ& gbz, bool progress) {
    if (progress) {
        std::cerr << "Estimating genome size" << std::endl;
    }
    size_t genome_size = 0;

    if (gbz.graph.get_path_count() > 0) {
        gbz.graph.for_each_path_of_sense({PathSense::REFERENCE, PathSense::GENERIC}, [&](const path_handle_t& path_handle) {
            std::string path_name = gbz.graph.get_path_name(path_handle);
            if (!Paths::is_alt(path_name)) {
                gbz.graph.for_each_step_in_path(path_handle, [&](const step_handle_t& step_handle) {
                    handle_t handle = gbz.graph.get_handle_of_step(step_handle);
                    genome_size += gbz.graph.get_length(handle);
                });
            }
        });
        if (progress) {
            std::cerr << "Estimated size based on reference / non-alt generic paths: " << genome_size << std::endl;
        }
    }

    if (genome_size == 0) {
        gbz.graph.for_each_handle([&](const handle_t& handle) {
            genome_size += gbz.graph.get_length(handle);
        });
        if (progress) {
            std::cerr << "Estimated size based on total sequence length: " << genome_size << std::endl;
        }
    }

    // Genome size / 2 should be a reasonably tight upper bound for the number of kmers
    // with any specific base in the middle position.
    size_t hash_table_size = gbwtgraph::KmerIndex<gbwtgraph::Key64>::minimum_size(genome_size / 2);
    if (progress) {
        std::cerr << "Estimated hash table size: 2^" << trailing_zeros(hash_table_size) << std::endl; 
    }

    return hash_table_size;
}

using key_type = gbwtgraph::DefaultMinimizerIndex::key_type;
using code_type = gbwtgraph::KmerEncoding::code_type;
using payload_t = ZipCode::payload_type;

std::vector<key_type> find_frequent_kmers(const gbwtgraph::GBZ& gbz, const MinimizerIndexParameters& params) {
    std::vector<key_type> frequent_kmers;
    if (!params.use_weighted_minimizers) {
        return frequent_kmers;
    }

    double start = gbwt::readTimer();
    if (params.progress) {
        std::string algorithm = (params.space_efficient_counting ? "space-efficient" : "fast");
        std::cerr << "Finding frequent kmers using the " << algorithm << " algorithm" << std::endl;
    }
    size_t hash_table_size = 0;
    if (params.hash_table_width == 0) {
        hash_table_size = estimate_hash_table_size(gbz, params.progress);
    } else {
        hash_table_size = size_t(1) << params.hash_table_width;
    }
    frequent_kmers = gbwtgraph::frequent_kmers<gbwtgraph::Key64>(
        gbz.graph, params.k, params.threshold, params.space_efficient_counting, hash_table_size
    );
    if (params.progress) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Found " << frequent_kmers.size() << " kmers with more than " << params.threshold << " hits in " << seconds << " seconds" << std::endl;
    }

    return frequent_kmers;
}

// TODO: Should we multithread this?
void cache_payloads(
    const gbwtgraph::GBZ& gbz,
    const SnarlDistanceIndex& distance_index,
    hash_map<nid_t, payload_t>& node_id_to_payload,
    ZipCodeCollection* oversized_zipcodes,
    bool progress
) {
    double start = gbwt::readTimer();
    if (progress) {
        std::cerr << "Caching payloads" << std::endl;
    }

    const handlegraph::HandleGraph* graph_ptr = (const handlegraph::HandleGraph*) &gbz.graph;

    gbz.graph.for_each_handle([&](const handle_t& handle) {
        nid_t node_id = gbz.graph.get_id(handle);
        ZipCode zipcode;
        pos_t pos = make_pos_t(node_id, false, 0);
        zipcode.fill_in_zipcode_from_pos(distance_index, pos, true, graph_ptr);
        payload_t payload = zipcode.get_payload_from_zip();
        if (payload == MIPayload::NO_CODE && oversized_zipcodes != nullptr) {
            // The zipcode is too large for the payload field.
            // Add it to the oversized zipcode list.
            zipcode.fill_in_full_decoder();
            size_t offset = oversized_zipcodes->size();
            oversized_zipcodes->emplace_back(zipcode);
            payload = { 0, offset };
        }
        node_id_to_payload.emplace(node_id, payload);
    });

    if (progress) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Cached payloads in " << seconds << " seconds" << std::endl;
    }
}

gbwtgraph::DefaultMinimizerIndex build_minimizer_index(
    const gbwtgraph::GBZ& gbz,
    const SnarlDistanceIndex* distance_index,
    ZipCodeCollection* oversized_zipcodes,
    const MinimizerIndexParameters& params
) {
    double start = gbwt::readTimer();

    // Minimal validation to avoid inconsistent parameters that would
    // otherwise require manual handling.
    std::string err = params.validate();
    if (!err.empty()) {
        std::cerr << "error: [build_minimizer_index()] " << err << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Determine payload size and type.
    size_t payload_size = 0;
    MinimizerIndexParameters::PayloadType payload_type = MinimizerIndexParameters::PAYLOAD_NONE;
    if (distance_index != nullptr) {
        payload_size = MinimizerIndexParameters::ZIPCODE_PAYLOAD_SIZE;
        if (params.paths_in_payload) {
            payload_size++;
            payload_type = MinimizerIndexParameters::PAYLOAD_ZIPCODES_WITH_PATHS;
        } else {
            payload_type = MinimizerIndexParameters::PAYLOAD_ZIPCODES;
        }
    }
    std::string payload_str = MinimizerIndexParameters::payload_str(payload_type);

    // Create an empty minimizer index.
    std::vector<key_type> frequent_kmers = find_frequent_kmers(gbz, params);
    gbwtgraph::DefaultMinimizerIndex index(params.k, params.w_or_s, payload_size, params.use_syncmers);
    if (params.use_weighted_minimizers && !frequent_kmers.empty()) {
        index.add_frequent_kmers(frequent_kmers, params.iterations);
    }
    index.set_tag(MinimizerIndexParameters::PAYLOAD_KEY, payload_str);

    // Build the index.
    if (params.progress) {
        std::cerr << "Building MinimizerIndex with k = " << index.k();
        if (index.uses_syncmers()) {
            std::cerr << ", s = " << index.s();
        } else {
            std::cerr << ", w = " << index.w();
        }
        std::cerr << ", payload = " << payload_str << std::endl;
    }

    if (distance_index == nullptr) {
        gbwtgraph::index_haplotypes(gbz, index, [](const pos_t&) { return nullptr; });
    } else {
        // Cache payloads before building the index.
        // A zipcode only depends on the node id.
        hash_map<nid_t, payload_t> node_id_to_payload;
        node_id_to_payload.reserve(gbz.graph.max_node_id() - gbz.graph.min_node_id());
        cache_payloads(gbz, *distance_index, node_id_to_payload, oversized_zipcodes, params.progress);

        auto get_payload = [&](const pos_t& pos) -> const code_type* {
            auto iter = node_id_to_payload.find(id(pos));
            if (iter != node_id_to_payload.end()) {
                return reinterpret_cast<const code_type*>(&iter->second);
            } else {
                return reinterpret_cast<const code_type*>(&MIPayload::NO_CODE);
            }
        };
        if (params.paths_in_payload) {
            gbwtgraph::index_haplotypes_with_paths(gbz, index, get_payload);
        } else {
            gbwtgraph::index_haplotypes(gbz, index, get_payload);
        }
    }

    // Index statistics.
    if (params.progress) {
        std::cerr << index.size() << " keys (" << index.unique_keys() << " unique)" << std::endl;
        std::cerr << "Minimizer occurrences: " << index.number_of_values() << std::endl;
        std::cerr << "Load factor: " << index.load_factor() << std::endl;
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Construction time: " << seconds << " seconds" << std::endl;
    }

    return index;
}

void require_payload(const gbwtgraph::DefaultMinimizerIndex& index, MinimizerIndexParameters::PayloadType expected_payload) {
    std::string found = index.get_tag(MinimizerIndexParameters::PAYLOAD_KEY);
    std::string expected_str = MinimizerIndexParameters::payload_str(expected_payload);
    if (found != expected_str) {
        std::cerr << "error: expected a minimizer index with payload type \""
            << expected_str << "\" but found \"" << found << "\"" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void require_payload(
    const gbwtgraph::DefaultMinimizerIndex& index,
    const std::vector<MinimizerIndexParameters::PayloadType>& expected_payloads
) {
    std::string found = index.get_tag(MinimizerIndexParameters::PAYLOAD_KEY);
    for (MinimizerIndexParameters::PayloadType type : expected_payloads) {
        std::string expected_str = MinimizerIndexParameters::payload_str(type);
        if (found == expected_str) {
            return;
        }
    }
    std::cerr << "error: expected a minimizer index with one of payload types {";
    for (size_t i = 0; i < expected_payloads.size(); ++i) {
        std::cerr << " \"" << MinimizerIndexParameters::payload_str(expected_payloads[i]) << "\"";
        if (i + 1 < expected_payloads.size()) {
            std::cerr << ",";
        }
    }
    std::cerr << " } but found \"" << found << "\"" << std::endl;
    std::exit(EXIT_FAILURE);
}

// TODO: I couldn't figure out a way to factor out a checker/stringifier
// function and share code between require_payload() and has_payload() without
// some cost for the abstraction.

bool has_payload(const gbwtgraph::DefaultMinimizerIndex& index, MinimizerIndexParameters::PayloadType payload) {
    std::string found = index.get_tag(MinimizerIndexParameters::PAYLOAD_KEY);
    std::string expected_str = MinimizerIndexParameters::payload_str(payload);
    return found == expected_str;
}

bool has_payload(
    const gbwtgraph::DefaultMinimizerIndex& index,
    const std::vector<MinimizerIndexParameters::PayloadType>& payloads
) {
    std::string found = index.get_tag(MinimizerIndexParameters::PAYLOAD_KEY);
    for (MinimizerIndexParameters::PayloadType type : payloads) {
        std::string expected_str = MinimizerIndexParameters::payload_str(type);
        if (found == expected_str) {
            return true;
        }
    }
    return false;
}

//------------------------------------------------------------------------------

/// Return a mapping of the original segment ids to a list of chopped node ids
/// (mimicking logic and interface from function of same name in gbwt_helper.cpp)
unordered_map<string, vector<nid_t>> load_translation_map(const gbwtgraph::GBWTGraph& graph) {
    unordered_map<string, vector<nid_t>> translation_map;
    graph.for_each_segment([&](const std::string& segment_id, std::pair<nid_t, nid_t> nodes) -> bool {
            vector<nid_t>& val = translation_map[segment_id];
            val.reserve(nodes.second - nodes.first);
            for (nid_t node_id = nodes.first; node_id < nodes.second; ++node_id) {
                val.push_back(node_id);
            }
            return true;
        });
    return translation_map;
}

/// Return a backwards mapping of chopped node to original segment position (id,offset pair)
/// (mimicking logic and interface from function of same name in gbwt_helper.cpp)
unordered_map<nid_t, pair<string, size_t>> load_translation_back_map(const gbwtgraph::GBWTGraph& graph) {
    unordered_map<nid_t, pair<string, size_t>> translation_back_map;
    graph.for_each_segment([&](const std::string& segment_id, std::pair<nid_t, nid_t> nodes) -> bool {
            size_t offset = 0;
            for (nid_t node_id = nodes.first; node_id < nodes.second; ++node_id) {
                translation_back_map[node_id] = make_pair(segment_id, offset);
                offset += graph.get_length(graph.get_handle(node_id));
            }
            return true;
        });
                
    return translation_back_map;
}

//------------------------------------------------------------------------------

std::string to_string_gbwtgraph(handle_t handle) {
    return to_string_gbwtgraph(gbwtgraph::GBWTGraph::handle_to_node(handle));
}

std::string to_string_gbwtgraph(gbwt::node_type node) {
    return std::string("(") + std::to_string(gbwt::Node::id(node)) + std::string(", ") + std::to_string(gbwt::Node::is_reverse(node)) + std::string(")");
}

//------------------------------------------------------------------------------

} // namespace vg

#include "gbwtgraph_helper.hpp"
#include "gbwt_helper.hpp"

#include <vg/io/vpkg.hpp>

namespace vg {

//------------------------------------------------------------------------------

gbwtgraph::GFAParsingParameters get_best_gbwtgraph_gfa_parsing_parameters() {
    gbwtgraph::GFAParsingParameters parameters;
    // Configure GBWTGraph GFA parsing to be as close to the vg GFA parser as we can get.
    // TODO: Make it closer.
    parameters.path_name_formats.clear();
    // Parse panSN with a fragment after it as a haplotype.
    parameters.path_name_formats.emplace_back(
        gbwtgraph::GFAParsingParameters::PAN_SN_REGEX + "#([0-9][0-9]*)",
        gbwtgraph::GFAParsingParameters::PAN_SN_FIELDS + "F",
        PathSense::HAPLOTYPE
    );
    // Parse panSN with a range after it as a normal but with a fragment based
    // on start position.
    parameters.path_name_formats.emplace_back(
        gbwtgraph::GFAParsingParameters::PAN_SN_REGEX + "\\[([0-9][0-9]*)(-[0-9]*)?\\]",
        gbwtgraph::GFAParsingParameters::PAN_SN_FIELDS + "F",
        gbwtgraph::GFAParsingParameters::PAN_SN_SENSE
    );
    // Parse standard panSN as what we think that is (reference)
    parameters.path_name_formats.emplace_back(
        gbwtgraph::GFAParsingParameters::PAN_SN_REGEX,
        gbwtgraph::GFAParsingParameters::PAN_SN_FIELDS,
        gbwtgraph::GFAParsingParameters::PAN_SN_SENSE
    );
    // Parse paths with just a name and a range as generic paths with a contig and a fragment
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
    std::unique_ptr<gbwtgraph::GBWTGraph> loaded = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_gbwtgraph()] cannot load GBWTGraph " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    graph = std::move(*loaded);
}

void load_gbz(gbwtgraph::GBZ& gbz, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading GBZ from " << filename << std::endl;
    }
    std::unique_ptr<gbwtgraph::GBZ> loaded = vg::io::VPKG::load_one<gbwtgraph::GBZ>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_gbz()] cannot load GBZ " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    gbz = std::move(*loaded);
}

void load_gbz(gbwt::GBWT& index, gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading GBWT and GBWTGraph from " << filename << std::endl;
    }
    std::unique_ptr<gbwtgraph::GBZ> loaded = vg::io::VPKG::load_one<gbwtgraph::GBZ>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_gbz()] cannot load GBZ " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(loaded->index);
    graph = std::move(loaded->graph);
    graph.set_gbwt(index); // We moved the GBWT out from the GBZ, so we have to update the pointer.
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
    std::unique_ptr<gbwtgraph::DefaultMinimizerIndex> loaded = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_minimizer()] cannot load MinimizerIndex " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void save_gbwtgraph(const gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving GBWTGraph to " << filename << std::endl;
    }
    graph.serialize(filename);
}

void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving GBZ to " << filename << std::endl;
    }
    sdsl::simple_sds::serialize_to(gbz, filename);
}

void save_gbz(const gbwt::GBWT& index, gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving GBWT and GBWTGraph to " << filename << std::endl;
    }
    std::ofstream out(filename, std::ios_base::binary);
    if (!out) {
        std::cerr << "error: [save_gbz()] cannot open file " << filename << " for writing" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    gbwtgraph::GBZ::simple_sds_serialize(index, graph, out);
    out.close();
}

void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& gbwt_name, const std::string& graph_name, bool show_progress) {
    save_gbwt(gbz.index, gbwt_name, show_progress);
    save_gbwtgraph(gbz.graph, graph_name, show_progress);
}

void save_minimizer(const gbwtgraph::DefaultMinimizerIndex& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving MinimizerIndex to " << filename << std::endl;
    }
    std::ofstream out(filename, std::ios_base::binary);
    if (!out) {
        std::cerr << "error: [save_minimizer()] cannot open file " << filename << " for writing" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index.serialize(out);
    out.close();
}

//------------------------------------------------------------------------------

} // namespace vg

#include "gbwtgraph_helper.hpp"
#include "gbwt_helper.hpp"

#include <vg/io/vpkg.hpp>

namespace vg {

//------------------------------------------------------------------------------

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
}

void load_gbz(gbwtgraph::GBZ& gbz, const std::string& gbwt_name, const std::string& graph_name, bool show_progress) {
    gbz = gbwtgraph::GBZ();
    load_gbwt(gbz.index, gbwt_name, show_progress);
    load_gbwtgraph(gbz.graph, graph_name, show_progress);
    gbz.graph.set_gbwt(gbz.index);
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

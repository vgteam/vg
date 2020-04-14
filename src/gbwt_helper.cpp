#include "gbwt_helper.hpp"
#include "utility.hpp"

#include <sstream>

namespace vg {

std::vector<std::string> parseGenotypes(const std::string& vcf_line, size_t num_samples) {
    std::vector<std::string> result;

    // The 9th tab-separated field should start with "GT".
    size_t offset = 0;
    for (int i = 0; i < 8; i++) {
        size_t pos = vcf_line.find('\t', offset);
        if (pos == std::string::npos) {
            std::cerr << "error: [vg index] VCF line does not contain genotype information" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        offset = pos + 1;
    }
    if (vcf_line.substr(offset, 2) != "GT") {
        std::cerr << "error: [vg index] VCF line does not contain genotype information" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Genotype strings are the first colon-separated fields in the 10th+ tab-separated fields.
    offset = vcf_line.find('\t', offset);
    while (offset != std::string::npos && offset + 1 < vcf_line.length()) {
        offset++;
        size_t pos = vcf_line.find_first_of("\t:", offset);
        if (pos == std::string::npos) {
            pos = vcf_line.length();
        }
        result.emplace_back(vcf_line.substr(offset, pos - offset));
        offset = vcf_line.find('\t', offset);
    }

    if (result.size() != num_samples) {
        std::cerr << "error: [vg index] expected " << num_samples << " samples, got " << result.size() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return result;
}

//------------------------------------------------------------------------------

std::string insert_gbwt_path(MutablePathHandleGraph& graph, const gbwt::GBWT& gbwt_index, gbwt::size_type id) {

    gbwt::size_type sequence_id = gbwt::Path::encode(id, false);
    if (sequence_id >= gbwt_index.sequences()) {
        std::cerr << "error: [insert_gbwt_path()] invalid path id: " << id << std::endl;
        return "";
    }

    std::string path_name = thread_name(gbwt_index, id);
    if (path_name.empty()) {
        path_name = std::to_string(id);
    }
    if (graph.has_path(path_name)) {
        std::cerr << "error: [insert_gbwt_path()] path name already exists: " << path_name << std::endl;
        return "";
    }

    path_handle_t handle = graph.create_path_handle(path_name);
    gbwt::edge_type pos = gbwt_index.start(sequence_id);
    while (pos.first != gbwt::ENDMARKER) {
        graph.append_step(handle, gbwt_to_handle(graph, pos.first));
        pos = gbwt_index.LF(pos);
    }

    return path_name;
}

//------------------------------------------------------------------------------

std::string thread_name(const gbwt::GBWT& gbwt_index, size_t id) {
    if (!gbwt_index.hasMetadata() || !gbwt_index.metadata.hasPathNames() || id >= gbwt_index.metadata.paths()) {
        return "";
    }

    const gbwt::PathName& path = gbwt_index.metadata.path(id);
    std::stringstream stream;
    stream << "_thread_";
    if (gbwt_index.metadata.hasSampleNames()) {
        stream << gbwt_index.metadata.sample(path.sample);
    } else {
        stream << path.sample;
    }
    stream << "_";
    if (gbwt_index.metadata.hasContigNames()) {
        stream << gbwt_index.metadata.contig(path.contig);
    } else {
        stream << path.contig;
    }
    stream << "_" << path.phase << "_" << path.count;
    return stream.str();
}

//------------------------------------------------------------------------------

gbwt::GBWT get_gbwt(const std::vector<gbwt::vector_type>& paths) {
    gbwt::size_type node_width = 1, total_length = 0;
    for (auto& path : paths) {
        for (auto node : path) {
            node_width = std::max(node_width, gbwt::bit_length(gbwt::Node::encode(node, true)));
        }
        total_length += 2 * (path.size() + 1);
    }

    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder builder(node_width, total_length);
    for (auto& path : paths) {
        builder.insert(path, true);
    }
    builder.finish();

    return gbwt::GBWT(builder.index);
}

//------------------------------------------------------------------------------

} // namespace vg

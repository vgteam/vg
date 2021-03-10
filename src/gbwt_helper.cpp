#include "gbwt_helper.hpp"
#include "utility.hpp"

#include <vg/io/vpkg.hpp>

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

gbwt::vector_type extract_as_gbwt_path(const PathHandleGraph& graph, const std::string& path_name) {
    gbwt::vector_type result;
    if (!graph.has_path(path_name)) {
        return result;
    }
    path_handle_t path_handle = graph.get_path_handle(path_name);
    result.reserve(graph.get_step_count(path_handle));
    for (handle_t handle : graph.scan_path(path_handle)) {
        result.push_back(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
    }
    return result;
}

gbwt::vector_type path_predecessors(const PathHandleGraph& graph, const std::string& path_name) {
    gbwt::vector_type result;
    if (!graph.has_path(path_name)) {
        return result;
    }
    path_handle_t path_handle = graph.get_path_handle(path_name);
    if (graph.get_step_count(path_handle) == 0) {
        return result;
    }
    step_handle_t step = graph.path_begin(path_handle);
    handle_t handle = graph.get_handle_of_step(step);
    graph.follow_edges(handle, true, [&] (handle_t prev) {
        if (prev != handle) {
            result.push_back(gbwt::Node::encode(graph.get_id(prev), graph.get_is_reverse(prev)));
        }
    });
    return result;
}

//------------------------------------------------------------------------------

gbwt::size_type gbwt_node_width(const HandleGraph& graph) {
    return gbwt::bit_length(gbwt::Node::encode(graph.max_node_id(), true));
}

void finish_gbwt_constuction(gbwt::GBWTBuilder& builder,
    const std::vector<std::string>& sample_names,
    const std::vector<std::string>& contig_names,
    size_t haplotype_count, bool print_metadata,
    const std::string& header) {

    builder.finish();
    builder.index.metadata.setSamples(sample_names);
    builder.index.metadata.setHaplotypes(haplotype_count);
    builder.index.metadata.setContigs(contig_names);
    if (print_metadata) {
        #pragma omp critical
        {
            std::cerr << header << ": ";
            gbwt::operator<<(std::cerr, builder.index.metadata);
            std::cerr << std::endl;
        }
    }
}

//------------------------------------------------------------------------------

void load_gbwt(const std::string& filename, gbwt::GBWT& index, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading compressed GBWT from " << filename << std::endl;
    }
    std::unique_ptr<gbwt::GBWT> loaded = vg::io::VPKG::load_one<gbwt::GBWT>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_gbwt()] could not load compressed GBWT " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void load_gbwt(const std::string& filename, gbwt::DynamicGBWT& index, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading dynamic GBWT from " << filename << std::endl;
    }
    std::unique_ptr<gbwt::DynamicGBWT> loaded = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_gbwt()] could not load dynamic GBWT " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void GBWTHandler::use_compressed() {
    if (this->in_use == index_compressed) {
        return;
    } else if (this->in_use == index_dynamic) {
        if (this->show_progress) {
            std::cerr << "Converting dynamic GBWT into compressed GBWT" << std::endl;
        }
        this->compressed = gbwt::GBWT(this->dynamic);
        this->dynamic = gbwt::DynamicGBWT();
        this->in_use = index_compressed;
    } else {
        load_gbwt(this->filename, this->compressed, this->show_progress);
        this->in_use = index_compressed;
    }
}

void GBWTHandler::use_dynamic() {
    if (this->in_use == index_dynamic) {
        return;
    } else if (this->in_use == index_compressed) {
        if (this->show_progress) {
            std::cerr << "Converting compressed GBWT into dynamic GBWT" << std::endl;
        }
        this->dynamic = gbwt::DynamicGBWT(this->compressed);
        this->compressed = gbwt::GBWT();
        this->in_use = index_dynamic;
    } else {
        load_gbwt(this->filename, this->dynamic, this->show_progress);
        this->in_use = index_dynamic;
    }
}

void GBWTHandler::use(gbwt::GBWT& new_index) {
    this->clear();
    this->compressed.swap(new_index);
    this->in_use = index_compressed;
}

void GBWTHandler::use(gbwt::DynamicGBWT& new_index) {
    this->clear();
    this->dynamic.swap(new_index);
    this->in_use = index_dynamic;
}

void GBWTHandler::unbacked() {
    this->filename = std::string();
}

void GBWTHandler::serialize(const std::string& new_filename) {
    if (this->show_progress) {
        std::cerr << "Serializing the GBWT to " << new_filename << std::endl;
    }
    if (this->in_use == index_none) {
        std::cerr << "warning: [GBWTHandler] no GBWT to serialize" << std::endl;
        return;
    } else if (this->in_use == index_compressed) {
        vg::io::VPKG::save(this->compressed, new_filename);
    } else {
        vg::io::VPKG::save(this->dynamic, new_filename);
    }
    this->filename = new_filename;
}

void GBWTHandler::clear() {
    this->compressed = gbwt::GBWT();
    this->dynamic = gbwt::DynamicGBWT();
    this->in_use = index_none;
}

//------------------------------------------------------------------------------

std::vector<gbwt::size_type> threads_for_sample(const gbwt::GBWT& gbwt_index, const std::string& sample_name) {
    if (gbwt_index.hasMetadata() && gbwt_index.metadata.hasSampleNames() && gbwt_index.metadata.hasPathNames()) {
        gbwt::size_type sample_id = gbwt_index.metadata.sample(sample_name);
        if (sample_id < gbwt_index.metadata.samples()) {
            return gbwt_index.metadata.pathsForSample(sample_id);
        }
    }
    return std::vector<gbwt::size_type>();
}

std::vector<gbwt::size_type> threads_for_contig(const gbwt::GBWT& gbwt_index, const std::string& contig_name) {
    if (gbwt_index.hasMetadata() && gbwt_index.metadata.hasContigNames() && gbwt_index.metadata.hasPathNames()) {
        gbwt::size_type contig_id = gbwt_index.metadata.contig(contig_name);
        if (contig_id < gbwt_index.metadata.contigs()) {
            return gbwt_index.metadata.pathsForContig(contig_id);
        }
    }
    return std::vector<gbwt::size_type>();
}

std::string insert_gbwt_path(MutablePathHandleGraph& graph, const gbwt::GBWT& gbwt_index, gbwt::size_type id, std::string path_name) {

    gbwt::size_type sequence_id = gbwt::Path::encode(id, false);
    if (sequence_id >= gbwt_index.sequences()) {
        std::cerr << "error: [insert_gbwt_path()] invalid path id: " << id << std::endl;
        return "";
    }

    // If the path name was not specified, try first using the default name generated from GBWT metadata.
    // If that fails, simply use the id.
    if (path_name.empty()) { path_name = thread_name(gbwt_index, id); }
    if (path_name.empty()) { path_name = std::to_string(id); }
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

Path extract_gbwt_path(const HandleGraph& graph, const gbwt::GBWT& gbwt_index, gbwt::size_type id) {

    Path result;
    gbwt::size_type sequence_id = gbwt::Path::encode(id, false);
    if (sequence_id >= gbwt_index.sequences()) {
        std::cerr << "error: [insert_gbwt_path()] invalid path id: " << id << std::endl;
        return result;
    }

    std::string path_name = thread_name(gbwt_index, id);
    if (path_name.empty()) {
        path_name = std::to_string(id);
    }
    result.set_name(path_name);

    gbwt::edge_type pos = gbwt_index.start(sequence_id);
    size_t rank = 1;
    while (pos.first != gbwt::ENDMARKER) {
        Mapping* m = result.add_mapping();
        Position* p = m->mutable_position();
        p->set_node_id(gbwt::Node::id(pos.first));
        p->set_is_reverse(gbwt::Node::is_reverse(pos.first));
        Edit* e = m->add_edit();
        size_t len = graph.get_length(gbwt_to_handle(graph, pos.first));
        e->set_to_length(len);
        e->set_from_length(len);
        m->set_rank(rank);
        pos = gbwt_index.LF(pos);
        rank++;
    }

    return result;
}

std::string thread_name(const gbwt::GBWT& gbwt_index, gbwt::size_type id) {
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

std::string thread_sample(const gbwt::GBWT& gbwt_index, gbwt::size_type id) {
    if (!gbwt_index.hasMetadata() || !gbwt_index.metadata.hasPathNames() || id >= gbwt_index.metadata.paths()) {
        return "";
    }
    
    const gbwt::PathName& path = gbwt_index.metadata.path(id);
    std::stringstream stream;
    if (gbwt_index.metadata.hasSampleNames()) {
        stream << gbwt_index.metadata.sample(path.sample);
    } else {
        stream << path.sample;
    }
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

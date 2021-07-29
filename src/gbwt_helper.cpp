#include "gbwt_helper.hpp"
#include "utility.hpp"

#include <vg/io/vpkg.hpp>

#include <sstream>
#include <unordered_map>

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

void load_gbwt(gbwt::GBWT& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading compressed GBWT from " << filename << std::endl;
    }
    std::unique_ptr<gbwt::GBWT> loaded = vg::io::VPKG::load_one<gbwt::GBWT>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_gbwt()] cannot load compressed GBWT " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void load_gbwt(gbwt::DynamicGBWT& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading dynamic GBWT from " << filename << std::endl;
    }
    std::unique_ptr<gbwt::DynamicGBWT> loaded = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_gbwt()] cannot load dynamic GBWT " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void load_r_index(gbwt::FastLocate& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading r-index from " << filename << std::endl;
    }
    std::unique_ptr<gbwt::FastLocate> loaded = vg::io::VPKG::load_one<gbwt::FastLocate>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_r_index()] cannot load r-index " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void save_gbwt(const gbwt::GBWT& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving compressed GBWT to " << filename << std::endl;
    }
    sdsl::simple_sds::serialize_to(index, filename);
}

void save_gbwt(const gbwt::DynamicGBWT& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving dynamic GBWT to " << filename << std::endl;
    }
    sdsl::simple_sds::serialize_to(index, filename);
}

void save_r_index(const gbwt::FastLocate& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving r-index to " << filename << std::endl;
    }
    if (!sdsl::store_to_file(index, filename)) {
        std::cerr << "error: [save_r_index()] cannot write r-index to " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

//------------------------------------------------------------------------------

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
        load_gbwt(this->compressed, this->filename, this->show_progress);
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
        load_gbwt(this->dynamic, this->filename, this->show_progress);
        this->in_use = index_dynamic;
    }
}

void GBWTHandler::use(gbwt::GBWT& new_index) {
    this->clear();
    this->unbacked();
    this->compressed.swap(new_index);
    this->in_use = index_compressed;
}

void GBWTHandler::use(gbwt::DynamicGBWT& new_index) {
    this->clear();
    this->unbacked();
    this->dynamic.swap(new_index);
    this->in_use = index_dynamic;
}

void GBWTHandler::unbacked() {
    this->filename = std::string();
}

void GBWTHandler::serialize(const std::string& new_filename) {
    if (this->in_use == index_none) {
        std::cerr << "warning: [GBWTHandler] no GBWT to serialize" << std::endl;
        return;
    } else if (this->in_use == index_compressed) {
        save_gbwt(this->compressed, new_filename, this->show_progress);
    } else {
        save_gbwt(this->dynamic, new_filename, this->show_progress);
    }
    this->filename = new_filename;
}

void GBWTHandler::clear() {
    this->compressed = gbwt::GBWT();
    this->dynamic = gbwt::DynamicGBWT();
    this->in_use = index_none;
}

//------------------------------------------------------------------------------

gbwt::GBWT rebuild_gbwt(const gbwt::GBWT& gbwt_index, const std::vector<std::pair<gbwt::vector_type, gbwt::vector_type>>& mappings) {
    if (gbwt_index.empty() || mappings.empty()) {
        return gbwt_index;
    }

    // Partition the mappings by the first node and determine node width.
    gbwt::size_type node_width = sdsl::bits::length(gbwt_index.sigma() - 1);
    std::unordered_map<gbwt::node_type, std::vector<std::pair<gbwt::vector_type, gbwt::vector_type>>> mappings_by_first_node;
    for (auto& mapping : mappings) {
        if (mapping.first.empty()) {
            continue;
        }
        mappings_by_first_node[mapping.first.front()].push_back(mapping);
        for (auto node : mapping.second) {
            node_width = std::max(node_width, static_cast<gbwt::size_type>(sdsl::bits::length(node)));
        }
    }

    // Build the new GBWT.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::GBWTBuilder builder(node_width);
    for (gbwt::size_type id = 0; id < gbwt_index.sequences(); id += 2) {
        gbwt::vector_type path = gbwt_index.extract(id);
        gbwt::vector_type mapped;
        size_t i = 0;
        while (i < path.size()) {
            auto iter = mappings_by_first_node.find(path[i]);
            bool found = false;
            if (iter != mappings_by_first_node.end()) {
                for (auto& mapping : iter->second) {
                    size_t j = 1;
                    while (j < mapping.first.size() && i + j < path.size() && mapping.first[j] == path[i + j]) {
                        j++;
                    }
                    if (j >= mapping.first.size()) {
                        mapped.insert(mapped.end(), mapping.second.begin(), mapping.second.end());
                        i += j;
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                mapped.push_back(path[i]);
                i++;
            }
        }
        builder.insert(mapped, true);
    }
    builder.finish();

    // Copy the metadata.
    if (gbwt_index.hasMetadata()) {
        builder.index.addMetadata();
        builder.index.metadata = gbwt_index.metadata;
    }

    return gbwt::GBWT(builder.index);
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

std::string thread_name(const gbwt::GBWT& gbwt_index, gbwt::size_type id, bool short_name) {
    if (!gbwt_index.hasMetadata() || !gbwt_index.metadata.hasPathNames() || id >= gbwt_index.metadata.paths()) {
        return "";
    }

    const gbwt::PathName& path = gbwt_index.metadata.path(id);
    std::stringstream stream;
    if (!short_name) {
        stream << "_thread_";
    }
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
    if (!short_name) {
        stream << "_" << path.phase << "_" << path.count;
    }
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

int thread_phase(const gbwt::GBWT& gbwt_index, gbwt::size_type id) {
    if (!gbwt_index.hasMetadata() || !gbwt_index.metadata.hasPathNames() || id >= gbwt_index.metadata.paths()) {
        return -1;
    }
    
    const gbwt::PathName& path = gbwt_index.metadata.path(id);
    return path.phase;
}

gbwt::size_type thread_count(const gbwt::GBWT& gbwt_index, gbwt::size_type id) {
    if (!gbwt_index.hasMetadata() || !gbwt_index.metadata.hasPathNames() || id >= gbwt_index.metadata.paths()) {
        return 0;
    }
    
    const gbwt::PathName& path = gbwt_index.metadata.path(id);
    return path.count;
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

unordered_map<nid_t, vector<nid_t>> load_translation_map(ifstream& input_stream) {
    string buffer;
    size_t line = 1;
    unordered_map<nid_t, vector<nid_t>> translation_map;
    try {
        while (getline(input_stream, buffer)) {
            vector<string> toks = split_delims(buffer, "\t");
            if (toks.size() == 3 && toks[0] == "T") {
                vector<string> toks2 = split_delims(toks[2], ",");
                nid_t segment_id = parse<nid_t>(toks[1]);
                vector<nid_t> node_ids;
                node_ids.reserve(toks2.size());
                for (string& node_id_str : toks2) {
                    node_ids.push_back(parse<nid_t>(node_id_str));
                }
                vector<nid_t>& val = translation_map[segment_id];
                if (!val.empty()) {
                    throw runtime_error("Segment " + toks[0] + " already in map");
                }
                translation_map[segment_id] = node_ids;
            } else {
                throw runtime_error("Invalid columns");
            }
            ++line;
        }
    } catch (const std::exception& e) {
        throw runtime_error("[load_translation_map] error: unable to parse line " + to_string(line) +
                            " of translation map: " + e.what());
    }
    return translation_map;
}

unordered_map<nid_t, pair<nid_t, size_t>> load_translation_back_map(HandleGraph& graph, ifstream& input_stream) {
    string buffer;
    size_t line = 1;
    unordered_map<nid_t, pair<nid_t, size_t>> translation_back_map;
    try {
        while (getline(input_stream, buffer)) {
            vector<string> toks = split_delims(buffer, "\t");
            if (toks.size() == 3 && toks[0] == "T") {
                vector<string> toks2 = split_delims(toks[2], ",");
                nid_t segment_id = stol(toks[1]);
                size_t offset = 0;
                for (string& node_id_str : toks2) {
                    nid_t node_id = stol(node_id_str);
                    if (translation_back_map.count(node_id)) {
                        throw runtime_error("Node ID " + node_id_str + " already in map");
                    }
                    translation_back_map[node_id] = make_pair(segment_id, offset);
                    offset += graph.get_length(graph.get_handle(node_id));
                }
            } else {
                throw runtime_error("Invalid columns");
            }
            ++line;
        }
    } catch (const std::exception& e) {
        throw runtime_error("[load_translation_back_map] error: unable to parse line " + to_string(line) +
                            " of translation map: " + e.what());
    }
    return translation_back_map;
}


} // namespace vg

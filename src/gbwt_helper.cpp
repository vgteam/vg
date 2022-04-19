#include "gbwt_helper.hpp"
#include "utility.hpp"

#include <vg/io/vpkg.hpp>
#include <handle.hpp>
#include <gbwtgraph/utils.h>

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

// Partition the GBWT seqeuences between jobs by the first node.
std::vector<std::vector<gbwt::size_type>> partition_gbwt_sequences(const gbwt::GBWT& gbwt_index, const std::unordered_map<nid_t, size_t>& node_to_job, size_t num_jobs) {
    std::vector<std::vector<gbwt::size_type>> result(num_jobs);
    for (gbwt::size_type sequence = 0; sequence < gbwt_index.sequences(); sequence += 2) {
        gbwt::edge_type start = gbwt_index.start(sequence);
        if (start != gbwt::invalid_edge()) {
            nid_t node = gbwt::Node::id(start.first);
            auto iter = node_to_job.find(node);
            if (iter != node_to_job.end()) {
                result[iter->second].push_back(sequence);
            } else if (start.first == gbwt::ENDMARKER) {
                result[0].push_back(sequence);
            }
        }
    }
    return result;
}

// Build a GBWT by inserting the specified sequences and applying the specified mappings.
gbwt::GBWT rebuild_gbwt_job(const gbwt::GBWT& gbwt_index, const RebuildJob& job, const std::vector<gbwt::size_type>& sequences, const RebuildParameters& parameters) {

    // Partition the mappings by the first node and determine node width.
    gbwt::size_type node_width = sdsl::bits::length(gbwt_index.sigma() - 1);
    std::unordered_map<gbwt::node_type, std::vector<RebuildJob::mapping_type>> mappings_by_first_node;
    for (auto& mapping : job.mappings) {
        if (mapping.first.empty() || mapping.first == mapping.second) {
            continue;
        }
        mappings_by_first_node[mapping.first.front()].push_back(mapping);
        std::pair<gbwt::vector_type, gbwt::vector_type> reverse;
        gbwt::reversePath(mapping.first, reverse.first);
        gbwt::reversePath(mapping.second, reverse.second);
        mappings_by_first_node[reverse.first.front()].push_back(reverse);
        for (auto node : mapping.second) {
            node_width = std::max(node_width, static_cast<gbwt::size_type>(sdsl::bits::length(node)));
        }
    }

    // Insert the sequences from the original GBWT and apply the mappings.
    gbwt::GBWTBuilder builder(node_width, parameters.batch_size, parameters.sample_interval);
    for (gbwt::size_type sequence : sequences) {
        gbwt::vector_type path = gbwt_index.extract(sequence);
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
                        // Leave the last node unprocessed if it does not change.
                        if (mapping.first.size() > 1 && mapping.second.size() > 0 && mapping.first.back() == mapping.second.back()) {
                            mapped.insert(mapped.end(), mapping.second.begin(), mapping.second.end() - 1);
                            i += mapping.first.size() - 1;
                        } else {
                            mapped.insert(mapped.end(), mapping.second.begin(), mapping.second.end());
                            i += mapping.first.size();
                        }
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

    return gbwt::GBWT(builder.index);
}

// Copy metadata from source to target, but reorder path names according to the merging order.
void copy_metadata(const gbwt::GBWT& source, gbwt::GBWT& target, const std::vector<std::vector<gbwt::size_type>>& jobs, const std::vector<size_t>& job_order) {
    if (!source.hasMetadata()) {
        return;
    }
    target.addMetadata();
    target.metadata = source.metadata;
    if (source.metadata.hasPathNames()) {
        target.metadata.clearPathNames();
        for (size_t job : job_order) {
            for (gbwt::size_type sequence : jobs[job]) {
                target.metadata.addPath(source.metadata.path(gbwt::Path::id(sequence)));
            }
        }
    }
}

gbwt::GBWT rebuild_gbwt(const gbwt::GBWT& gbwt_index,
                        const std::vector<RebuildJob>& jobs,
                        const std::unordered_map<nid_t, size_t>& node_to_job,
                        const RebuildParameters& parameters) {

    if (gbwt_index.empty() || jobs.empty()) {
        return gbwt_index;
    }
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);

    // Sort the jobs in descending order by size.
    std::vector<size_t> jobs_by_size(jobs.size());
    for (size_t i = 0; i < jobs_by_size.size(); i++) {
        jobs_by_size[i] = i;
    }
    std::sort(jobs_by_size.begin(), jobs_by_size.end(), [&](size_t a, size_t b) -> bool {
        return (jobs[a].total_size > jobs[b].total_size);
    });

    // Build indexes in parallel.
    if (parameters.show_progress) {
        std::cerr << "rebuild_gbwt(): Building " << jobs.size() << " partial GBWTs using up to " << parameters.num_jobs << " parallel jobs" << std::endl;
    }
    std::vector<gbwt::GBWT> indexes(jobs.size());
    std::vector<std::vector<gbwt::size_type>> sequences_by_job = partition_gbwt_sequences(gbwt_index, node_to_job, jobs.size());
    int old_max_threads = omp_get_max_threads();
    omp_set_num_threads(parameters.num_jobs);
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < jobs.size(); i++) {
        size_t job = jobs_by_size[i];
        if (parameters.show_progress) {
            #pragma omp critical
            {
                std::cerr << "rebuild_gbwt(): Starting job " << job << std::endl;
            }
        }
        indexes[job] = rebuild_gbwt_job(gbwt_index, jobs[job], sequences_by_job[job], parameters);
        if (parameters.show_progress) {
            #pragma omp critical
            {
                std::cerr << "rebuild_gbwt(): Inserted " << sequences_by_job[job].size() << " threads in job " << job << std::endl;
            }
        }
    }
    omp_set_num_threads(old_max_threads);

    // We can avoid merging if we had only one job.
    if (indexes.size() == 1) {
        gbwt::GBWT result(std::move(indexes.front()));
        if (gbwt_index.hasMetadata()) {
            result.addMetadata();
            result.metadata = gbwt_index.metadata;
        }
        return result;
    }

    // Merge the partial GBWTs and copy the metadata.
    if (parameters.show_progress) {
        std::cerr << "rebuild_gbwt(): Merging the partial GBWTs" << std::endl;
    }
    gbwt::GBWT merged(indexes);
    indexes.clear();
    copy_metadata(gbwt_index, merged, sequences_by_job, jobs_by_size);

    return merged;
}

gbwt::GBWT rebuild_gbwt(const gbwt::GBWT& gbwt_index, const std::vector<RebuildJob::mapping_type>& mappings) {
    std::vector<RebuildJob> jobs {
        { mappings, 0 }
    };
    std::unordered_map<nid_t, size_t> node_to_job;
    for (gbwt::size_type i = 0; i < gbwt_index.sequences(); i += 2) {
        gbwt::edge_type start = gbwt_index.start(i);
        if (start != gbwt::invalid_edge()) {
            node_to_job[gbwt::Node::id(start.first)] = 0;
        }
    }
    RebuildParameters parameters;
    return rebuild_gbwt(gbwt_index, jobs, node_to_job, parameters);
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

    gbwt::size_type sequence_id = gbwt_index.bidirectional() ? gbwt::Path::encode(id, false) : id;

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
    gbwt::size_type sequence_id = gbwt_index.bidirectional() ? gbwt::Path::encode(id, false) : id;
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

    auto& metadata = gbwt_index.metadata;
    const gbwt::PathName& path = metadata.path(id);

    if (short_name) {
        // We want a name with just sample and contig.
        // Spit out a name in reference sense format, which should suffice.
        return PathMetadata::create_path_name(PathSense::REFERENCE,
                                              metadata.hasSampleNames() ? metadata.sample(path.sample) : std::to_string(path.sample),
                                              metadata.hasContigNames() ? metadata.contig(path.contig) : std::to_string(path.contig),
                                              PathMetadata::NO_HAPLOTYPE,
                                              PathMetadata::NO_PHASE_BLOCK,
                                              PathMetadata::NO_SUBRANGE);
    }

    // Try and work out the sense of the path just from the GBWT.
    // TODO: needs to be kept in sync with GBWTGraph.
    PathSense path_sense;
    if (metadata.hasSampleNames()) {
        auto sample_name = metadata.sample(path.sample);
        if (sample_name.size() <= gbwtgraph::NAMED_PATH_SAMPLE_PREFIX.size() &&
            std::equal(gbwtgraph::NAMED_PATH_SAMPLE_PREFIX.begin(), gbwtgraph::NAMED_PATH_SAMPLE_PREFIX.end(), sample_name.begin())) {
            // It starts with the right prefix

            if (sample_name.size() == gbwtgraph::NAMED_PATH_SAMPLE_PREFIX.size()) {
                // Actually only the prefix is there

                // If we assign the path to the special reference sample it's a generic
                // path with just the one locus name.
                path_sense = PathSense::GENERIC;
            } else {
                // If we assign it to some other sample that starts with that string, it's a
                // reference path and not a generic one.
                path_sense = PathSense::REFERENCE;
            }
        } else {
            // Otherwise it's a haplotype thread.
            path_sense = PathSense::HAPLOTYPE;
        }
    } else {
        // Otherwise it's a haplotype thread.
        path_sense = PathSense::HAPLOTYPE;
    }

    // Compose a name based on the sense.
    switch(path_sense) {
    case PathSense::GENERIC:
        if (metadata.hasContigNames()) {
            // The contig name is the exposed path name.
            return metadata.contig(path.contig);
        } else {
            // We must have sample names but no contig names. This probably shouldn't happen.
            // Just use the contig number.
            return std::to_string(path.contig);
        }
        break;
    case PathSense::REFERENCE: // Fall-through
    case PathSense::HAPLOTYPE:
        // The path name must be composed.
        // If contig names are missing, make them up.
        return PathMetadata::create_path_name(path_sense,
                                              metadata.sample(path.sample),
                                              metadata.hasContigNames() ? metadata.contig(path.contig) : std::to_string(path.contig),
                                              path.phase,
                                              path_sense == PathSense::HAPLOTYPE ? path.count : PathMetadata::NO_PHASE_BLOCK,
                                              PathMetadata::NO_SUBRANGE);
        break;
    default:
        throw std::runtime_error("Unimplemented sense!");
    }
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

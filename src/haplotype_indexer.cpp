/**
 * \file haplotype_indexer.cpp: implementations of haplotype indexing with the GBWT
 */

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include <vg/io/stream.hpp>

#include "gbwt_helper.hpp"

#include "haplotype_indexer.hpp"

#include "path.hpp"
#include "alignment.hpp"

using namespace std;

namespace vg {

HaplotypeIndexer::HaplotypeIndexer() {
    // Build the GBWTs silently.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
}

std::vector<std::string> HaplotypeIndexer::parse_vcf(const std::string& filename, const PathHandleGraph& graph, const std::string& job_name) const {

    // Parse all non-alt paths.
    std::vector<path_handle_t> path_handles;
    graph.for_each_path_handle([&](path_handle_t path_handle) {
        std::string path_name = graph.get_path_name(path_handle);
        if (!Paths::is_alt(path_name)) {
            path_handles.push_back(path_handle);
        }
    });

    return this->parse_vcf(filename, graph, path_handles, job_name);
}

std::vector<std::string> HaplotypeIndexer::parse_vcf(const std::string& filename, const PathHandleGraph& graph, const std::vector<path_handle_t>& paths, const std::string& job_name) const {

    // Open the VCF file.
    vcflib::VariantCallFile variant_file;
    variant_file.parseSamples = false; // vcflib parsing is very slow if there are many samples.
    std::string temp_filename = filename;
    variant_file.open(temp_filename);
    if (!variant_file.is_open()) {
        std::cerr << "error: [HaplotypeIndexer::parse_vcf] could not open " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // How many samples are there?
    size_t num_samples = variant_file.sampleNames.size();
    if (num_samples == 0) {
        std::cerr << "error: [HaplotypeIndexer::parse_vcf] variant file '" << filename << "' does not contain phasings" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Determine the samples we want to index.
    std::pair<size_t, size_t> sample_range = this->sample_range;
    sample_range.second = std::min(sample_range.second, num_samples);
    std::vector<std::string> sample_names(variant_file.sampleNames.begin() + sample_range.first, variant_file.sampleNames.begin() + sample_range.second);
    if (this->show_progress) {
        #pragma omp critical
        {
            std::cerr << job_name << ": Parsing VCF file " << filename << " with options";
            if (!this->phase_homozygous) {
                std::cerr << " --actual-phasing";
            }
            if (this->force_phasing) {
                std::cerr << " --force-phasing";
            }
            if (this->discard_overlaps) {
                std::cerr << " --discard-overlaps";
            }
            if (!this->rename_variants) {
                std::cerr << "  --vcf-variants";
            }
            std::cerr << std::endl;
            std::cerr << job_name << ": Samples " << sample_range.first << " to " << (sample_range.second - 1) << ", batch size " << samples_in_batch << std::endl;
        }
    }

    // Parse the contigs we are interested in.
    std::vector<std::string> result;
    size_t total_variants_processed = 0;
    std::mt19937 rng(0xDEADBEEF);
    std::uniform_int_distribution<std::mt19937::result_type> random_bit(0, 1);
    size_t found_missing_variants = 0;
    for (size_t path_id = 0; path_id < paths.size(); path_id++) {
        std::string path_name = graph.get_path_name(paths[path_id]);
        std::string vcf_contig_name = (this->path_to_vcf.count(path_name) > 0 ? this->path_to_vcf.at(path_name) : path_name);

        // Set the VCF region or process the entire contig.
        if (this->regions.count(vcf_contig_name)) {
            std::pair<size_t, size_t> region = this->regions.at(vcf_contig_name);
            variant_file.setRegion(vcf_contig_name, region.first, region.second);
        } else {
            variant_file.setRegion(vcf_contig_name);
        }

        // Check that the VCF file contains this contig.
        vcflib::Variant var(variant_file);
        if (!(variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == vcf_contig_name)) {
            std::cerr << "warning: [HaplotypeIndexer::parse_vcf] contig " << vcf_contig_name << " not present in file " << filename << std::endl;
            continue;
        }
        if (this->show_progress) {
            #pragma omp critical
            {
                std::cerr << job_name << ": Path " << path_name << " matches VCF contig " << vcf_contig_name;
                if (this->regions.count(vcf_contig_name)) {
                    std::pair<size_t, size_t> region = this->regions.at(vcf_contig_name);
                    std::cerr << ", region " << region.first << " to " << region.second;
                }
                std::cerr << std::endl;
            }
        }

        // Structures to parse the VCF file into.
        std::string parse_file = (this->batch_file_prefix.empty() ? gbwt::TempFile::getName("parse") : this->batch_file_prefix + '_' + vcf_contig_name);
        gbwt::VariantPaths variants(graph.get_step_count(paths[path_id]));
        variants.setSampleNames(sample_names);
        variants.setContigName(path_name);
        std::vector<gbwt::PhasingInformation> phasings;

        // Add the reference to VariantPaths.
        for (handle_t handle : graph.scan_path(paths[path_id])) {
            variants.appendToReference(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
        }
        variants.indexReference();

        // Create a PhasingInformation file for each batch.
        for (size_t batch_start = sample_range.first; batch_start < sample_range.second; batch_start += samples_in_batch) {
            size_t batch_size = std::min(samples_in_batch, sample_range.second - batch_start);
            if (!this->batch_file_prefix.empty()) {
                // Use a permanent file.
                phasings.emplace_back(parse_file, batch_start, batch_size);
            } else {
                // Use a temporary file that persists until the program exits.
                phasings.emplace_back(batch_start, batch_size);
                phasings.back().makePersistent();
            }
            variants.addFile(phasings.back().name(), phasings.back().offset(), phasings.back().size());
        }

        // Parse the variants and the phasings.
        size_t variants_processed = 0;
        std::vector<bool> was_diploid(sample_range.second, true); // Was the sample diploid at the previous site?
        do {
            // Skip variants with non-DNA sequence, as they are not included in the graph.
            bool isDNA = allATGC(var.ref);
            for (std::vector<std::string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                if (!allATGC(*a)) isDNA = false;
            }
            if (!isDNA) {
                continue;
            }
            
            if (this->rename_variants) {
                // We need to move the variant over to the contig name
                // used in the graph, in order to get the right id for
                // it in the graph.
                var.sequenceName = path_name;
            }

            // Determine the reference nodes for the current variant and create a variant site.
            // If the variant is not an insertion, there should be a path for the ref allele.
            std::string var_name = make_variant_id(var);
            std::string ref_path_name = "_alt_" + var_name + "_0";
            gbwt::vector_type ref_path = extract_as_gbwt_path(graph, ref_path_name);
            size_t ref_pos = variants.invalid_position();
            if (!ref_path.empty()) {
                ref_pos = variants.firstOccurrence(ref_path.front());
                if (ref_pos == variants.invalid_position()) {
                    #pragma omp critical
                    {
                        std::cerr << "warning: [HaplotypeIndexer::parse_vcf] invalid ref path for " << var_name << " at "
                            << var.sequenceName << ":" << var.position << std::endl;
                    }
                    continue;
                }
            } else {
                // Try using the alternate alleles instead.
                bool found = false;
                for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                    std::string alt_path_name = "_alt_" + var_name + "_" + std::to_string(alt_index);
                    size_t candidate_pos = 0;
                    bool candidate_found = false;
                    gbwt::vector_type pred_nodes = path_predecessors(graph, alt_path_name);
                    if (!pred_nodes.empty()) {
                        for (auto node : pred_nodes) {
                            size_t pred_pos = variants.firstOccurrence(node);
                            if (pred_pos != variants.invalid_position()) {
                                candidate_pos = std::max(candidate_pos, pred_pos + 1);
                                candidate_found = true;
                                found = true;
                            }
                        }
                        // For each alternate allele, find the rightmost reference node among
                        // its predecessors. If multiple alleles have candidates for the
                        // reference position, choose the leftmost one.
                        if (candidate_found) {
                            ref_pos = std::min(ref_pos, candidate_pos);
                        }
                    }
                }
                if (!found) {
                    // This variant from the VCF is just not in the graph, so skip it.
                    found_missing_variants++;
                    if (this->warn_on_missing_variants && found_missing_variants <= this->max_missing_variant_warnings) {
                        #pragma omp critical
                        {
                            // The user might not know it. Warn them in case they mixed up their VCFs.
                            std::cerr << "warning: [HaplotypeIndexer::parse_vcf] alt and ref paths for " << var_name
                                    << " at " << var.sequenceName << ":" << var.position
                                    << " missing/empty! Was the variant skipped during construction?" << std::endl;
                            if (found_missing_variants == this->max_missing_variant_warnings) {
                                std::cerr << "warning: [HaplotypeIndexer::parse_vcf] suppressing further missing variant warnings" << std::endl;
                            }
                        }
                    }
                    continue;
                }
            }
            variants.addSite(ref_pos, ref_pos + ref_path.size());

            // Add alternate alleles to the site.
            for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                std::string alt_path_name = "_alt_" + var_name + "_" + std::to_string(alt_index);
                if (graph.has_path(alt_path_name)) {
                    variants.addAllele(extract_as_gbwt_path(graph, alt_path_name));
                } else {
                    variants.addAllele(ref_path);
                }
            }

            // Store the phasings in PhasingInformation structures.
            std::vector<std::string> genotypes = parseGenotypes(var.originalLine, num_samples);
            for (size_t batch = 0; batch < phasings.size(); batch++) {
                std::vector<gbwt::Phasing> current_phasings;
                for (size_t sample = phasings[batch].offset(); sample < phasings[batch].limit(); sample++) {
                    current_phasings.emplace_back(genotypes[sample], was_diploid[sample], this->phase_homozygous);
                    was_diploid[sample] = current_phasings.back().diploid;
                    if(this->force_phasing) {
                        current_phasings.back().forcePhased([&]() {
                            return random_bit(rng);
                        });
                    }
                }
                phasings[batch].append(current_phasings);
            }
            variants_processed++;
        }
        while (variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == vcf_contig_name); // End of variants.
        if (this->show_progress) {
            size_t phasing_bytes = 0;
            for (size_t batch = 0; batch < phasings.size(); batch++) {
                phasing_bytes += phasings[batch].bytes();
            }
            #pragma omp critical
            {
                std::cerr << job_name << ": Processed " << variants_processed << " variants on path " << path_name << ", " << gbwt::inMegabytes(phasing_bytes) << " MiB phasing information" << std::endl;
                std::cerr << job_name << ": Saving the VCF parse for path " << path_name << " to " << parse_file << std::endl;
            }
        }

        // Save the VCF parse.
        if (!sdsl::store_to_file(variants, parse_file)) {
            std::cerr << "error: [HaplotypeIndexer::parse_vcf] cannot write parse file " << parse_file << std::endl;
            std::exit(EXIT_FAILURE);
        }
        result.push_back(parse_file);

        // End of haplotype generation for the current contig.
        total_variants_processed += variants_processed;
    } // End of contigs.
        
    if (this->warn_on_missing_variants && found_missing_variants > 0) {
        #pragma omp critical
        {
            std::cerr << "warning: [HaplotypeIndexer::parse_vcf] Found " << found_missing_variants << "/" << total_variants_processed
                << " variants in phasing VCF but not in graph! Do your graph and VCF match?" << std::endl;
        }
    }

    return result;
}

std::unique_ptr<gbwt::DynamicGBWT> HaplotypeIndexer::build_gbwt(const std::vector<std::string>& vcf_parse_files, const std::string& job_name) const {

    // GBWT metadata.
    std::vector<std::string> sample_names, contig_names;
    std::set<gbwt::range_type> haplotypes;

    // GBWT index.
    std::unique_ptr<gbwt::DynamicGBWT> index(new gbwt::DynamicGBWT());
    index->addMetadata();
    if (vcf_parse_files.empty()) {
        return index;
    }

    // Construction for each contig.
    for (const std::string& filename : vcf_parse_files) {
        gbwt::VariantPaths variants;
        if (!sdsl::load_from_file(variants, filename)) {
            std::cerr << "error: [HaplotypeIndexer::build_gbwt] cannot load VCF parse from " << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (!variants.hasContigName() || !variants.hasSampleNames()) {
            std::cerr << "error: [HaplotypeIndexer::build_gbwt] VCF parse file " << filename << " does not contain sample/contig names" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (sample_names.empty()) {
            sample_names = variants.getSampleNames();
        } else if (sample_names != variants.getSampleNames()) {
            std::cerr << "error: [HaplotypeIndexer::build_gbwt] invalid sample names in VCF parse file " << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }
        contig_names.emplace_back(variants.getContigName());
        if (this->show_progress) {
            #pragma omp critical
            {
                std::cerr << job_name << ": Generating haplotypes for path " << variants.getContigName() << " from file " << filename << std::endl;
            }
        }
        gbwt::GBWTBuilder builder(variants.nodeWidth(true), this->gbwt_buffer_size * gbwt::MILLION, this->id_interval);
        builder.swapIndex(*index);
        gbwt::generateHaplotypes(variants, std::set<std::string>(),
        [&](gbwt::size_type sample_id) -> bool {
            return (this->excluded_samples.find(sample_names[sample_id]) == this->excluded_samples.end());
        }, [&](const gbwt::Haplotype& haplotype) {
            builder.insert(haplotype.path, true); // Insert in both orientations.
            builder.index.metadata.addPath(haplotype.sample, contig_names.size() - 1, haplotype.phase, haplotype.count);
            haplotypes.insert(gbwt::range_type(haplotype.sample, haplotype.phase));
        }, [&](gbwt::size_type, gbwt::size_type) -> bool {
            // For each overlap, discard it if our global flag is set.
            return this->discard_overlaps;
        });
        builder.finish();
        builder.swapIndex(*index);
    }

    // Finish the construction.
    index->metadata.setSamples(sample_names);
    index->metadata.setContigs(contig_names);
    index->metadata.setHaplotypes(haplotypes.size());
    if (this->show_progress) {
        std::cerr << job_name << ": ";
        gbwt::operator<<(std::cerr, index->metadata);
        std::cerr << std::endl;
    }
    return index;
}

std::unique_ptr<gbwt::DynamicGBWT> HaplotypeIndexer::build_gbwt(const PathHandleGraph& graph) const {

    // GBWT metadata.
    std::vector<std::string> sample_names, contig_names;
    size_t haplotype_count = 0;
    if (this->paths_as_samples) {
        contig_names.push_back("0"); // An artificial contig.
    }
    else {
        sample_names.push_back("ref"); // An artificial sample.
        haplotype_count = 1;
    }

    // GBWT construction.
    gbwt::GBWTBuilder builder(gbwt_node_width(graph), this->gbwt_buffer_size * gbwt::MILLION, this->id_interval);
    builder.index.addMetadata();

    // Actual work.
    if (show_progress) {
        #pragma omp critical
        {
            std::cerr << "Indexing embedded paths" << std::endl;
        }
    }
    graph.for_each_path_handle([&](path_handle_t path_handle) {
        std::string path_name = graph.get_path_name(path_handle);
        if (graph.is_empty(path_handle) || Paths::is_alt(path_name)) {
            return;
        }
        gbwt::vector_type buffer;
        buffer.reserve(graph.get_step_count(path_handle));
        for (handle_t handle : graph.scan_path(path_handle)) {
            buffer.push_back(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
        }
        builder.insert(buffer, true); // Insert in both orientations.
        if (this->paths_as_samples) {
            builder.index.metadata.addPath(sample_names.size(), 0, 0, 0);
            sample_names.push_back(path_name);
            haplotype_count++;
        } else {
            builder.index.metadata.addPath(0, contig_names.size(), 0, 0);
            contig_names.push_back(path_name);
        }
    });
        
    // Finish the construction and extract the index.
    finish_gbwt_constuction(builder, sample_names, contig_names, haplotype_count, this->show_progress);
    std::unique_ptr<gbwt::DynamicGBWT> built(new gbwt::DynamicGBWT());
    builder.swapIndex(*built);
    return built;
}

std::unique_ptr<gbwt::DynamicGBWT> HaplotypeIndexer::build_gbwt(const PathHandleGraph& graph,
    const std::vector<std::string>& aln_filenames, const std::string& aln_format) const {

    // GBWT metadata.
    std::vector<std::string> sample_names, contig_names;
    std::map<std::string, std::pair<size_t, size_t>> sample_info; // name -> (id, count)
    contig_names.push_back("0"); // An artificial contig.
    size_t haplotype_count = 0;

    // GBWT construction.
    gbwt::GBWTBuilder builder(gbwt_node_width(graph), this->gbwt_buffer_size * gbwt::MILLION, this->id_interval);
    builder.index.addMetadata();

    // Actual work.
    if (this->show_progress) {
        #pragma omp critical
        {
            std::cerr << "Converting " << aln_format << " to threads" << std::endl;
        }
    }
    std::function<void(Alignment&)> lambda = [&](Alignment& aln) {
        gbwt::vector_type buffer;
        for (auto& m : aln.path().mapping()) {
            buffer.push_back(mapping_to_gbwt(m));
        }
        builder.insert(buffer, true); // Insert in both orientations.
        size_t sample_id = 0, sample_count = 0;
        auto iter = sample_info.find(aln.name());
        if (iter == sample_info.end()) {
            sample_id = sample_names.size();
            sample_names.push_back(aln.name());
            sample_info[aln.name()] = std::pair<size_t, size_t>(sample_id, sample_count);
            haplotype_count++;
        } else {
            sample_id = iter->second.first;
            sample_count = iter->second.second;
            iter->second.second++;
        }
        builder.index.metadata.addPath(sample_id, 0, 0, sample_count);
    };
    for (auto& file_name : aln_filenames) {
        if (aln_format == "GAM") {
            get_input_file(file_name, [&](istream& in) {
                vg::io::for_each(in, lambda);
            });
        } else {
            assert(aln_format == "GAF");
            vg::io::gaf_unpaired_for_each(graph, file_name, lambda);
        }
    }
        
    // Finish the construction and extract the index.
    finish_gbwt_constuction(builder, sample_names, contig_names, haplotype_count, this->show_progress);
    std::unique_ptr<gbwt::DynamicGBWT> built(new gbwt::DynamicGBWT());
    builder.swapIndex(*built);
    return built;
}

}

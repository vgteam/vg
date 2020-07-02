/**
 * \file haplotype_indexer.cpp: implementations of haplotype indexing with the GBWT
 */

#include <iostream>
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
    // Use the same temp directory as VG for GBWT temp files.
    gbwt::TempFile::setDirectory(temp_file::get_dir());

    // Build the GBWTs silently.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
}

size_t HaplotypeIndexer::parse_vcf(const PathHandleGraph* graph, const std::vector<path_handle_t>& contigs,
    vcflib::VariantCallFile& variant_file, std::vector<std::string>& sample_names,
    const function<void(size_t, const gbwt::VariantPaths&, gbwt::PhasingInformation&)>& handle_contig_haplotype_batch) const {

    size_t total_variants_processed = 0;
    std::mt19937 rng(0xDEADBEEF);
    std::uniform_int_distribution<std::mt19937::result_type> random_bit(0, 1);

    // How many samples are there?
    size_t num_samples = variant_file.sampleNames.size();
    if (num_samples == 0) {
        std::cerr << "error: [HaplotypeIndexer::parse_vcf] the variant file does not contain phasings" << std::endl;
        std::exit(1);
    }

    // Determine the samples we want to index.
    std::pair<size_t, size_t> sample_range = this->sample_range;
    sample_range.second = std::min(sample_range.second, num_samples);
    sample_names.clear();
    sample_names.insert(sample_names.end(),
                        variant_file.sampleNames.begin() + sample_range.first,
                        variant_file.sampleNames.begin() + sample_range.second);
    size_t haplotype_count = 2 * (sample_range.second - sample_range.first); // Assuming a diploid genome
    if (show_progress) {
        cerr << "Haplotype generation parameters:" << endl;
        cerr << "- Samples " << sample_range.first << " to " << (sample_range.second - 1) << endl;
        cerr << "- Batch size " << samples_in_batch << endl;
        if (this->phase_homozygous) {
            cerr << "- Phase homozygous genotypes" << endl;
        }
        if (this->force_phasing) {
            cerr << "- Force phasing" << endl;
        }
        if (this->discard_overlaps) {
            cerr << "- Discard overlaps" << endl;
        }
    }

    // Parse the contigs we are interested in.
    size_t found_missing_variants = 0;
    for (size_t contig = 0; contig < contigs.size(); contig++) {
        std::string path_name = graph->get_path_name(contigs[contig]);
        std::string vcf_contig_name = (this->path_to_vcf.count(path_name) > 0 ? this->path_to_vcf.at(path_name) : path_name);
        if (this->show_progress) {
            cerr << "Processing path " << path_name << " as VCF contig " << vcf_contig_name << endl;
        }
        std::string parse_file = (this->batch_file_prefix.empty() ? "" : this->batch_file_prefix + '_' + vcf_contig_name);

        // Structures to parse the VCF file into.
        gbwt::VariantPaths variants(graph->get_step_count(contigs[contig]));
        variants.setSampleNames(sample_names);
        variants.setContigName(path_name);
        std::vector<gbwt::PhasingInformation> phasings;

        // Add the reference to VariantPaths.
        for (handle_t handle : graph->scan_path(contigs[contig])) {
            variants.appendToReference(gbwt::Node::encode(graph->get_id(handle), graph->get_is_reverse(handle)));
        }
        variants.indexReference();

        // Create a PhasingInformation file for each batch.
        for (size_t batch_start = sample_range.first; batch_start < sample_range.second; batch_start += samples_in_batch) {
            size_t batch_size = std::min(samples_in_batch, sample_range.second - batch_start);
            if (!this->batch_file_prefix.empty()) {
                // Use a permanent file.
                phasings.emplace_back(parse_file, batch_start, batch_size);
                variants.addFile(phasings.back().name(), phasings.back().offset(), phasings.back().size());
            } else {
                // Use a temporary file.
                phasings.emplace_back(batch_start, batch_size);
            }
        }

        // Set the VCF region or process the entire contig.
        if (this->regions.count(vcf_contig_name)) {
            std::pair<size_t, size_t> region = this->regions.at(vcf_contig_name);
            if (show_progress) {
                std::cerr << "- Setting region " << region.first << " to " << region.second << std::endl;
            }
            variant_file.setRegion(vcf_contig_name, region.first, region.second);
        } else {
            variant_file.setRegion(vcf_contig_name);
        }
        
        if (this->rename_variants && this->show_progress) {
            std::cerr << "- Moving variants from " << vcf_contig_name << " to " << path_name << std::endl;
        }

        // Parse the variants and the phasings.
        vcflib::Variant var(variant_file);
        size_t variants_processed = 0;
        std::vector<bool> was_diploid(sample_range.second, true); // Was the sample diploid at the previous site?
        while (variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == vcf_contig_name) {
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
            gbwt::vector_type ref_path = extract_as_gbwt_path(*graph, ref_path_name);
            size_t ref_pos = variants.invalid_position();
            if (!ref_path.empty()) {
                ref_pos = variants.firstOccurrence(ref_path.front());
                if (ref_pos == variants.invalid_position()) {
                    std::cerr << "warning: [HaplotypeIndexer::parse_vcf] invalid ref path for " << var_name << " at "
                         << var.sequenceName << ":" << var.position << std::endl;
                    continue;
                }
            } else {
                // Try using the alternate alleles instead.
                bool found = false;
                for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                    std::string alt_path_name = "_alt_" + var_name + "_" + std::to_string(alt_index);
                    size_t candidate_pos = 0;
                    bool candidate_found = false;
                    gbwt::vector_type pred_nodes = path_predecessors(*graph, alt_path_name);
                    if (!pred_nodes.empty()) {
                        for (auto node : pred_nodes) {
                            size_t pred_pos = variants.firstOccurrence(node);
                            if (pred_pos != variants.invalid_position()) {
#ifdef debug
                                cerr << "Found predecessor node " << gbwt::Node::id(node) << " " << gbwt::Node::is_reverse(node)
                                     << " occurring at valid pos " << pred_pos << endl;
#endif
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
                        // The user might not know it. Warn them in case they mixed up their VCFs.
                        std::cerr << "warning: [HaplotypeIndexer::parse_vcf] alt and ref paths for " << var_name
                                << " at " << var.sequenceName << ":" << var.position
                                << " missing/empty! Was the variant skipped during construction?" << std::endl;
                        if (found_missing_variants == this->max_missing_variant_warnings) {
                            std::cerr << "warning: [HaplotypeIndexer::parse_vcf] suppressing further missing variant warnings" << std::endl;
                        }
                    }
                    continue;
                }
            }
            variants.addSite(ref_pos, ref_pos + ref_path.size());

            // Add alternate alleles to the site.
            for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                std::string alt_path_name = "_alt_" + var_name + "_" + std::to_string(alt_index);
                if (graph->has_path(alt_path_name)) {
                    variants.addAllele(extract_as_gbwt_path(*graph, alt_path_name));
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
        } // End of variants.
        if (this->show_progress) {
            std::cerr << "- Parsed " << variants_processed << " variants" << std::endl;
            size_t phasing_bytes = 0;
            for (size_t batch = 0; batch < phasings.size(); batch++) {
                phasing_bytes += phasings[batch].bytes();
            }
            std::cerr << "- Phasing information: " << gbwt::inMegabytes(phasing_bytes) << " MB" << std::endl;
        }

        // Save memory by closing the phasings files.
        for (size_t batch = 0; batch < phasings.size(); batch++) {
            phasings[batch].close();
        }

        // Save the VCF parse
        if (!this->batch_file_prefix.empty()) {
            if (!sdsl::store_to_file(variants, parse_file)) {
                std::cerr << "error: [HaplotypeIndexer::parse_vcf] cannot write parse file " << parse_file << std::endl;
                std::exit(1);
            }
        }

        for (size_t batch = 0; batch < phasings.size(); batch++) {
            // Send all the batches to our callback.
            handle_contig_haplotype_batch(contig, variants, phasings[batch]);
        }

        // End of haplotype generation for the current contig.
        // Record the number of variants we saw on this contig
        total_variants_processed += variants_processed;
    
    } // End of contigs.
        
    if (this->warn_on_missing_variants && found_missing_variants > 0) {
        std::cerr << "warning: [HaplotypeIndexer::parse_vcf] Found " << found_missing_variants << "/" << total_variants_processed
             << " variants in phasing VCF but not in graph! Do your graph and VCF match?" << std::endl;
    }
    
    return haplotype_count;
}

std::unique_ptr<gbwt::DynamicGBWT> HaplotypeIndexer::build_gbwt(const PathHandleGraph* graph, std::string vcf_filename) const {

    // GBWT metadata.
    std::vector<std::string> sample_names, contig_names;
    size_t haplotype_count = 0;

    // GBWT construction.
    gbwt::GBWTBuilder builder(gbwt_node_width(*graph), this->gbwt_buffer_size * gbwt::MILLION, this->id_interval);
    builder.index.addMetadata();

    // Open the VCF file.
    size_t total_variants_processed = 0;
    vcflib::VariantCallFile variant_file;
    variant_file.parseSamples = false; // vcflib parsing is very slow if there are many samples.
    variant_file.open(vcf_filename);
    if (!variant_file.is_open()) {
        std::cerr << "error: [HaplotypeIndexer::build_gbwt] could not open " << vcf_filename << std::endl;
        throw runtime_error("Could not open " + vcf_filename);
    } else if (this->show_progress) {
        std::cerr << "Opened variant file " << vcf_filename << std::endl;
    }
    
    // Determine the non-alt paths and use their names as contig names.
    std::vector<path_handle_t> path_handles;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        std::string path_name = graph->get_path_name(path_handle);
        if (!Paths::is_alt(path_name)) {
            path_handles.push_back(path_handle);
            contig_names.push_back(path_name);
        }
    });

    // Actual work.
    unordered_set<gbwt::size_type> skipped_sample_numbers;
    size_t parsed_haplotypes = this->parse_vcf(graph, path_handles, variant_file, sample_names,
        [&](size_t contig, const gbwt::VariantPaths& variants, gbwt::PhasingInformation& phasings_batch) {
        gbwt::generateHaplotypes(variants, phasings_batch, [&](gbwt::size_type sample) -> bool {
            // Decide if we should process this sample or not.
            if (excluded_samples.find(variant_file.sampleNames[sample]) == excluded_samples.end()) {
                return true;
            } else {
                skipped_sample_numbers.insert(sample);
                return false;
            }
        }, [&](const gbwt::Haplotype& haplotype) {
            builder.insert(haplotype.path, true); // Insert in both orientations.
            builder.index.metadata.addPath({
                static_cast<gbwt::PathName::path_name_type>(haplotype.sample),
                static_cast<gbwt::PathName::path_name_type>(contig),
                static_cast<gbwt::PathName::path_name_type>(haplotype.phase),
                static_cast<gbwt::PathName::path_name_type>(haplotype.count)
            });
        }, [&](gbwt::size_type, gbwt::size_type) -> bool {
            // For each overlap, discard it if our global flag is set.
            return this->discard_overlaps;
        });
        if (this->show_progress) {
            std::cerr << "- Processed samples " << phasings_batch.offset() << " to " << (phasings_batch.offset() + phasings_batch.size() - 1) << std::endl;
        }
    });
    haplotype_count += parsed_haplotypes - skipped_sample_numbers.size() * 2;
        
    // Finish the construction and extract the index.
    finish_gbwt_constuction(builder, sample_names, contig_names, haplotype_count, this->show_progress);
    std::unique_ptr<gbwt::DynamicGBWT> built(new gbwt::DynamicGBWT());
    builder.swapIndex(*built);
    return built;
}

std::unique_ptr<gbwt::DynamicGBWT> HaplotypeIndexer::build_gbwt(const PathHandleGraph* graph) const {

    // GBWT metadata.
    std::vector<std::string> sample_names, contig_names;
    sample_names.push_back("ref"); // An artificial sample.
    size_t haplotype_count = 1;

    // GBWT construction.
    gbwt::GBWTBuilder builder(gbwt_node_width(*graph), this->gbwt_buffer_size * gbwt::MILLION, this->id_interval);
    builder.index.addMetadata();

    // Actual work.
    if (show_progress) {
        cerr << "Converting paths to threads..." << endl;
    }
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        std::string path_name = graph->get_path_name(path_handle);
        if (graph->is_empty(path_handle) || Paths::is_alt(path_name)) {
            return;
        }
        gbwt::vector_type buffer;
        buffer.reserve(graph->get_step_count(path_handle));
        for (handle_t handle : graph->scan_path(path_handle)) {
            buffer.push_back(gbwt::Node::encode(graph->get_id(handle), graph->get_is_reverse(handle)));
        }
        builder.insert(buffer, true); // Insert in both orientations.
        builder.index.metadata.addPath({
            static_cast<gbwt::PathName::path_name_type>(0),
            static_cast<gbwt::PathName::path_name_type>(contig_names.size()),
            static_cast<gbwt::PathName::path_name_type>(0),
            static_cast<gbwt::PathName::path_name_type>(0)
        });
        contig_names.push_back(path_name);
    });
        
    // Finish the construction and extract the index.
    finish_gbwt_constuction(builder, sample_names, contig_names, haplotype_count, this->show_progress);
    std::unique_ptr<gbwt::DynamicGBWT> built(new gbwt::DynamicGBWT());
    builder.swapIndex(*built);
    return built;
}

std::unique_ptr<gbwt::DynamicGBWT> HaplotypeIndexer::build_gbwt(const PathHandleGraph* graph,
    const std::vector<std::string>& aln_filenames, const std::string& aln_format) const {

    // GBWT metadata.
    std::vector<std::string> sample_names, contig_names;
    contig_names.push_back("0"); // An artificial contig.
    size_t haplotype_count = 0;

    // GBWT construction.
    gbwt::GBWTBuilder builder(gbwt_node_width(*graph), this->gbwt_buffer_size * gbwt::MILLION, this->id_interval);
    builder.index.addMetadata();

    // Actual work.
    if (this->show_progress) {
        std::cerr << "Converting " << aln_format << " to threads..." << std::endl;
    }
    std::function<void(Alignment&)> lambda = [&](Alignment& aln) {
        gbwt::vector_type buffer;
        for (auto& m : aln.path().mapping()) {
            buffer.push_back(mapping_to_gbwt(m));
        }
        builder.insert(buffer, true); // Insert in both orientations.
        builder.index.metadata.addPath({
            static_cast<gbwt::PathName::path_name_type>(sample_names.size()),
            static_cast<gbwt::PathName::path_name_type>(0),
            static_cast<gbwt::PathName::path_name_type>(0),
            static_cast<gbwt::PathName::path_name_type>(0)
        });
        sample_names.push_back(aln.name());
        haplotype_count++;
    };
    for (auto& file_name : aln_filenames) {
        if (aln_format == "GAM") {
            get_input_file(file_name, [&](istream& in) {
                vg::io::for_each(in, lambda);
            });
        } else {
            assert(aln_format == "GAF");
            gaf_unpaired_for_each(*graph, file_name, lambda);
        }
    }
        
    // Finish the construction and extract the index.
    finish_gbwt_constuction(builder, sample_names, contig_names, haplotype_count, this->show_progress);
    std::unique_ptr<gbwt::DynamicGBWT> built(new gbwt::DynamicGBWT());
    builder.swapIndex(*built);
    return built;
}

}

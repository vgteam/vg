/**
 * \file index_manager.cpp: implementations of common indexing functionality
 */

#include <iostream>
#include <vector>
#include <string>

#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>

#include <gbwtgraph/minimizer.h>

#include "index_manager.hpp"
#include "gbwt_helper.hpp"

using namespace std;

namespace vg {


size_t HaplotypeIndexer::parse_vcf(const PathHandleGraph* graph, map<string, Path>& alt_paths, const vector<path_handle_t>& contigs, vcflib::VariantCallFile& variant_file, std::vector<std::string>& sample_names, const function<void(size_t, const gbwt::VariantPaths&, const gbwt::PhasingInformation&)>& handle_contig_haplotype_batch) {

    // Use the same temp directory as VG.
    gbwt::TempFile::setDirectory(temp_file::get_dir());

    if (alt_paths.empty()) {
        // Populate alt paths from the graph if not already available.
        // TODO: avoid copying to Path objects and use straight from the graph.
        graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            if (Paths::is_alt(path_name)) {
                alt_paths[path_name] = path_from_path_handle(*graph, path_handle);
            }
        });
    }

    size_t total_variants_processed = 0;
    std::mt19937 rng(0xDEADBEEF);
    std::uniform_int_distribution<std::mt19937::result_type> random_bit(0, 1);

    // How many samples are there?
    size_t num_samples = variant_file.sampleNames.size();
    if (num_samples == 0) {
        cerr << "error: [HaplotypeIndexer::parse_vcf] the variant file does not contain phasings" << endl;
        exit(1);
    }

    // Determine the samples we want to index.
    // TODO: this modifies the object.
    sample_range.second = std::min(sample_range.second, num_samples);
    
    sample_names.insert(sample_names.end(),
                        variant_file.sampleNames.begin() + sample_range.first,
                        variant_file.sampleNames.begin() + sample_range.second);
    size_t haplotype_count = 2 * (sample_range.second - sample_range.first);  // Assuming a diploid genome
    if (show_progress) {
        cerr << "Haplotype generation parameters:" << endl;
        cerr << "- Samples " << sample_range.first << " to " << (sample_range.second - 1) << endl;
        cerr << "- Batch size " << samples_in_batch << endl;
        if (phase_homozygous) {
            cerr << "- Phase homozygous genotypes" << endl;
        }
        if (force_phasing) {
            cerr << "- Force phasing" << endl;
        }
        if (discard_overlaps) {
            cerr << "- Discard overlaps" << endl;
        }
    }

    size_t max_path_rank = contigs.size();
    for (size_t path_rank = 1; path_rank <= max_path_rank; path_rank++) {
        string path_name = graph->get_path_name(contigs[path_rank - 1]);
        if (Paths::is_alt(path_name)) {
            continue;
        }
        string vcf_contig_name = path_to_vcf.count(path_name) ? path_to_vcf[path_name] : path_name;
        if (show_progress) {
            cerr << "Processing path " << path_name << " as VCF contig " << vcf_contig_name << endl;
        }
        string parse_file = batch_file_prefix.empty() ? "" : batch_file_prefix + '_' + vcf_contig_name;

        // Structures to parse the VCF file into.
        gbwt::VariantPaths variants(graph->get_step_count(graph->get_path_handle(path_name)));
        variants.setSampleNames(sample_names);
        variants.setContigName(path_name);
        std::vector<gbwt::PhasingInformation> phasings;

        // Add the reference to VariantPaths.
        for (handle_t handle : graph->scan_path(contigs[path_rank - 1])) {
            variants.appendToReference(gbwt::Node::encode(graph->get_id(handle), graph->get_is_reverse(handle)));
        }
        variants.indexReference();

        // Create a PhasingInformation for each batch.
        for (size_t batch_start = sample_range.first; batch_start < sample_range.second; batch_start += samples_in_batch) {
            if (!batch_file_prefix.empty()) {
                // Use a permanent file.
                phasings.emplace_back(parse_file, batch_start, std::min(samples_in_batch, sample_range.second - batch_start));
                variants.addFile(phasings.back().name(), phasings.back().offset(), phasings.back().size());
            } else {
                // Use a temporary file.
                phasings.emplace_back(batch_start, std::min(samples_in_batch, sample_range.second - batch_start));
            }
        }

        // Set the VCF region or process the entire contig.
        if (regions.count(vcf_contig_name)) {
            auto region = regions[vcf_contig_name];
            if (show_progress) {
                cerr << "- Setting region " << region.first << " to " << region.second << endl;
            }
            variant_file.setRegion(vcf_contig_name, region.first, region.second);
        } else {
            variant_file.setRegion(vcf_contig_name);
        }
        
        if (rename_variants && show_progress) {
            cerr << "- Moving variants from " << vcf_contig_name << " to " << path_name << endl;
        }

        // Parse the variants and the phasings.
        vcflib::Variant var(variant_file);
        size_t variants_processed = 0;
        std::vector<bool> was_diploid(sample_range.second, true); // Was the sample diploid at the previous site?
        while (variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == vcf_contig_name) {
            // Skip variants with non-DNA sequence, as they are not included in the graph.
            bool isDNA = allATGC(var.ref);
            for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                if (!allATGC(*a)) isDNA = false;
            }
            if (!isDNA) {
                continue;
            }
            
            if (rename_variants) {
                // We need to move the variant over to the contig name
                // used in the graph, in order to get the right id for
                // it in the graph.
                var.sequenceName = path_name;
            }

            // Determine the reference nodes for the current variant and create a variant site.
            // If the variant is not an insertion, there should be a path for the ref allele.
            
            std::string var_name = make_variant_id(var);
            std::string ref_path_name = "_alt_" + var_name + "_0";
            auto ref_path_iter = alt_paths.find(ref_path_name);
            gbwt::vector_type ref_path;
            size_t ref_pos = variants.invalid_position();
            if (ref_path_iter != alt_paths.end() && ref_path_iter->second.mapping_size() != 0) {
                ref_path = path_to_gbwt(ref_path_iter->second);
                ref_pos = variants.firstOccurrence(ref_path.front());
                if (ref_pos == variants.invalid_position()) {
                    cerr << "warning: [HaplotypeIndexer::parse_vcf] invalid ref path for " << var_name << " at "
                         << var.sequenceName << ":" << var.position << endl;
                    continue;
                }
            } else {
                // Try using alt paths instead.
                bool found = false;
                for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                    std::string alt_path_name = "_alt_" + var_name + "_" + to_string(alt_index);
                    size_t candidate_pos = 0;
                    bool candidate_found = false;
                    auto alt_path_iter = alt_paths.find(alt_path_name);
                    if (alt_path_iter != alt_paths.end()) {
                        gbwt::vector_type pred_nodes = path_predecessors(*graph, alt_path_iter->second);
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
                    // This variant from the VCF is just not in the graph

                    found_missing_variants++;

                    if (warn_on_missing_variants) {
                        if (found_missing_variants <= max_missing_variant_warnings) {
                            // The user might not know it. Warn them in case they mixed up their VCFs.
                            cerr << "warning: [HaplotypeIndexer::parse_vcf] alt and ref paths for " << var_name
                                 << " at " << var.sequenceName << ":" << var.position
                                 << " missing/empty! Was the variant skipped during construction?" << endl;
                            if (found_missing_variants == max_missing_variant_warnings) {
                                cerr << "warning: [HaplotypeIndexer::parse_vcf] suppressing further missing variant warnings" << endl;
                            }
                        }
                    }

                    // Skip this variant and move on to the next as if it never appeared.
                    continue;
                }
            }
            variants.addSite(ref_pos, ref_pos + ref_path.size());

            // Add alternate alleles to the site.
            for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                std::string alt_path_name = "_alt_" + var_name + "_" + to_string(alt_index);
                auto alt_path_iter = alt_paths.find(alt_path_name);
                if (alt_path_iter != alt_paths.end()) {
                    variants.addAllele(path_to_gbwt(alt_path_iter->second));
                } else {
                    variants.addAllele(ref_path);
                }
            }

            // Store the phasings in PhasingInformation structures.
            std::vector<std::string> genotypes = parseGenotypes(var.originalLine, num_samples);
            for (size_t batch = 0; batch < phasings.size(); batch++) {
                std::vector<gbwt::Phasing> current_phasings;
                for (size_t sample = phasings[batch].offset(); sample < phasings[batch].limit(); sample++) {
                    current_phasings.emplace_back(genotypes[sample], was_diploid[sample], phase_homozygous);
                    was_diploid[sample] = current_phasings.back().diploid;
                    if(force_phasing) {
                        current_phasings.back().forcePhased([&]() {
                                return random_bit(rng);
                            });
                    }
                }
                phasings[batch].append(current_phasings);
            }
            variants_processed++;
        } // End of variants.
        if (show_progress) {
            cerr << "- Parsed " << variants_processed << " variants" << endl;
            size_t phasing_bytes = 0;
            for (size_t batch = 0; batch < phasings.size(); batch++) {
                phasing_bytes += phasings[batch].bytes();
            }
            cerr << "- Phasing information: " << gbwt::inMegabytes(phasing_bytes) << " MB" << endl;
        }

        // Save memory:
        // - Delete the alt paths if we no longer need them.
        // - Close the phasings files.
        if (path_rank == max_path_rank) {
            alt_paths.clear();
            // TODO: also delete the graph if we don't need it anymore.
        }
        for (size_t batch = 0; batch < phasings.size(); batch++) {
            phasings[batch].close();
        }

        // Save the VCF parse
        if (!batch_file_prefix.empty()) {
            if (!sdsl::store_to_file(variants, parse_file)) {
                cerr << "error: [HaplotypeIndexer::parse_vcf] cannot write parse file " << parse_file << endl;
                return 1;
            }
        }
        
        for (size_t batch = 0; batch < phasings.size(); batch++) {
            // Send all the batches to our callback.
            handle_contig_haplotype_batch(path_rank - 1, variants, phasings[batch]);
        }
        
        // End of haplotype generation for the current contig.
    
        // Record the number of variants we saw on this contig
        total_variants_processed += variants_processed;
    
    } // End of contigs.
        
    if (warn_on_missing_variants && found_missing_variants > 0) {
        cerr << "warning: [HaplotypeIndexer::parse_vcf] Found " << found_missing_variants << "/" << total_variants_processed
             << " variants in phasing VCF but not in graph! Do your graph and VCF match?" << endl;
    }
    
    return haplotype_count;
}

unique_ptr<gbwt::GBWT> HaplotypeIndexer::make_gbwt(const PathHandleGraph* graph, bool index_paths, const string& vcf_filename, const vector<string>& gam_filenames) {
    unique_ptr<gbwt::GBWT> to_return;
    return to_return;
}

}

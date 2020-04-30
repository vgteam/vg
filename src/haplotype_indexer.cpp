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

using namespace std;

namespace vg {

HaplotypeIndexer::HaplotypeIndexer() {
    // Use the same temp directory as VG for GBWT temp files.
    gbwt::TempFile::setDirectory(temp_file::get_dir());
}

size_t HaplotypeIndexer::parse_vcf(const PathHandleGraph* graph, map<string, Path>& alt_paths, const vector<path_handle_t>& contigs, vcflib::VariantCallFile& variant_file, std::vector<std::string>& sample_names, const function<void(size_t, const gbwt::VariantPaths&, gbwt::PhasingInformation&)>& handle_contig_haplotype_batch) {

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

tuple<vector<string>, size_t, vector<string>> HaplotypeIndexer::generate_threads(const PathHandleGraph* graph, map<string, Path>& alt_paths,
    bool index_paths, const string& vcf_filename, const vector<string>& gam_filenames,
    const function<void(size_t)>& bit_width_ready, const function<void(const gbwt::vector_type&, const gbwt::size_type (&)[4])>& each_thread) {
    
    // GBWT metadata.
    std::vector<std::string> sample_names, contig_names;
    size_t haplotype_count = 0;
    size_t true_sample_offset = 0; // Id of the first VCF sample.

    // Determine node id width.
    size_t id_width;
    if (gam_filenames.empty()) {
        nid_t max_id = graph->max_node_id();
        if (show_progress) {
            cerr << "Maximum node id in graph: " << max_id << endl;
        }
        id_width = gbwt::bit_length(gbwt::Node::encode(max_id, true));
    } else { // indexing a GAM
        if (show_progress) {
            cerr << "Finding maximum node id in GAM..." << endl;
        }
        nid_t max_id = 0;
        size_t alignments_in_gam = 0;
        function<void(Alignment&)> lambda = [&](Alignment& aln) {
            gbwt::vector_type buffer;
            for (auto& m : aln.path().mapping()) {
                max_id = max((nid_t) m.position().node_id(), max_id);
            }
            alignments_in_gam++;
        };
        for (auto& file_name : gam_filenames) {
            get_input_file(file_name, [&](istream& in) {
                vg::io::for_each(in, lambda);
            });
        }
        id_width = gbwt::bit_length(gbwt::Node::encode(max_id, true));
        sample_names.reserve(alignments_in_gam); // We store alignment names as sample names.
    }
    
    if (show_progress) {
        cerr << "Node id width: " << id_width << endl;
    }
    
    // Announce it to the caller. They will need it to prepare a GBWT.
    // TODO: Why not always just use the max ID in the graph???
    bit_width_ready(id_width);

    // Store contig names.
    if (index_paths || !vcf_filename.empty()) {
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                string path_name = graph->get_path_name(path_handle);
                if (!Paths::is_alt(path_name)) {
                    contig_names.push_back(path_name);
                }
            });
    }
    // Convert paths to threads
    if (index_paths) {
        if (show_progress) {
            cerr << "Converting paths to threads..." << endl;
        }
        size_t path_rank = 0;
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                ++path_rank;
                if (graph->is_empty(path_handle) || Paths::is_alt(graph->get_path_name(path_handle))) {
                    return;
                }
                gbwt::vector_type buffer;
                buffer.reserve(graph->get_step_count(path_handle));
                for (handle_t handle : graph->scan_path(path_handle)) {
                    buffer.push_back(gbwt::Node::encode(graph->get_id(handle), graph->get_is_reverse(handle)));
                }
                each_thread(buffer, {true_sample_offset, path_rank - 1, 0, 0});
            });

        // GBWT metadata: We assume that the XG index contains the reference paths.
        sample_names.emplace_back("ref");
        haplotype_count++;
        true_sample_offset++;
    }

    // Index GAM, using alignment names as sample names.
    if (!gam_filenames.empty()) {
        if (show_progress) {
            cerr << "Converting GAM to threads..." << endl;
        }
        function<void(Alignment&)> lambda = [&](Alignment& aln) {
            gbwt::vector_type buffer;
            for (auto& m : aln.path().mapping()) {
                buffer.push_back(mapping_to_gbwt(m));
            }
            each_thread(buffer, {sample_names.size(), 0, 0, 0});
            sample_names.emplace_back(aln.name());
            haplotype_count++;
            true_sample_offset++;
        };
        for (auto& file_name : gam_filenames) {
            get_input_file(file_name, [&](istream& in) {
                vg::io::for_each(in, lambda);
            });
        }
    }

    // Generate haplotypes from VCF input
    if (!vcf_filename.empty()) {

        size_t total_variants_processed = 0;
        vcflib::VariantCallFile variant_file;
        variant_file.parseSamples = false; // vcflib parsing is very slow if there are many samples.
        // TODO: vcflib needs a non-const string to open a file
        string mutable_filename = vcf_filename;
        variant_file.open(mutable_filename);
        if (!variant_file.is_open()) {
            cerr << "error: [HaplotypeIndexer::generate_threads] could not open " << vcf_filename << endl;
            throw runtime_error("Could not open " + vcf_filename);
        } else if (show_progress) {
            cerr << "Opened variant file " << vcf_filename << endl;
        }
        
        // Process each VCF contig corresponding to an XG path.
        vector<path_handle_t> path_handles;
        // 1st pass: scan for all non-alt paths (they are handled separately)
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                if (!alt_paths.count(graph->get_path_name(path_handle))) {
                    path_handles.push_back(path_handle);
                }
            });
        size_t max_path_rank = path_handles.size();

        // Track all the skipped samples
        unordered_set<gbwt::size_type> skipped_sample_numbers;

        size_t parsed_hapolotypes = this->parse_vcf(graph, alt_paths, path_handles, variant_file, sample_names,
            [&](size_t contig, const gbwt::VariantPaths& variants, gbwt::PhasingInformation& phasings_batch) {
            
            // For each (modifiable) batch of phasing info for a contig (in serial)
                
            // We need to generte haplotypes from the parsed VCF.
            gbwt::generateHaplotypes(variants, phasings_batch, [&](gbwt::size_type sample) -> bool {
                // Decide if we should process this sample or not.
                if (excluded_samples.find(variant_file.sampleNames[sample]) == excluded_samples.end()) {
                    return true;
                } else {
                    // This sample is to be excluded. Remember that we
                    // saw it, so we can properly count haplotypes
                    // actually done.
                    skipped_sample_numbers.insert(sample);
                    return false;
                }
            }, [&](const gbwt::Haplotype& haplotype) {
                // Store each haplotype in the GBWT
                each_thread(haplotype.path, {haplotype.sample + true_sample_offset - sample_range.first,
                                             contig,
                                             haplotype.phase,
                                             haplotype.count});
            }, [&](gbwt::size_type, gbwt::size_type) -> bool {
                // For each overlap, discard it if our global flag is set.
                return discard_overlaps;
            });
            
            if (show_progress) {
                cerr << "- Processed samples " << phasings_batch.offset() << " to " << (phasings_batch.offset() + phasings_batch.size() - 1) << endl;
            }
        });

        // Assume all the skipped samples were diploid and back them out of the number of haplotypes.
        parsed_hapolotypes -= skipped_sample_numbers.size() * 2;

        // And add into thew total haplotype count, together with other sources
        haplotype_count += parsed_hapolotypes;
            
    } // End of haplotypes.
        
    // Clear out alt_paths since it is no longer needed.
    alt_paths.clear();
    
    // Return the info needed for GBWT metadata
    return make_tuple(sample_names, haplotype_count, contig_names);
}

unique_ptr<gbwt::DynamicGBWT> HaplotypeIndexer::build_gbwt(const PathHandleGraph* graph, map<string, Path>& alt_paths,
    bool index_paths, const string& vcf_filename, const vector<string>& gam_filenames) {
    
    // GBWT metadata.
    std::vector<std::string> sample_names, contig_names;
    size_t haplotype_count = 0;
    
    // GBWT builder object
    unique_ptr<gbwt::GBWTBuilder> gbwt_builder;
    
    std::tie(sample_names, haplotype_count, contig_names) = this->generate_threads(graph, alt_paths,
        index_paths, vcf_filename, gam_filenames, [&](size_t id_width) {
        
        // Start the GBWT construction.
        if (show_progress) {
            cerr << "GBWT parameters: buffer size " << gbwt_buffer_size << ", id interval " << id_interval << endl;
        }
        gbwt::Verbosity::set(gbwt::Verbosity::SILENT);  // Make the construction thread silent.
        gbwt_builder.reset(new gbwt::GBWTBuilder(id_width, gbwt_buffer_size * gbwt::MILLION, id_interval));
        gbwt_builder->index.addMetadata();
        
    }, [&](const gbwt::vector_type& thread, const gbwt::size_type (&thread_name)[4]) {
        // We got a thread!
        // Save the thread itself
        gbwt_builder->insert(thread, true); // Insert in both orientations.
        // Save thread metadata
        gbwt_builder->index.metadata.addPath({
            static_cast<gbwt::PathName::path_name_type>(thread_name[0]),
            static_cast<gbwt::PathName::path_name_type>(thread_name[1]),
            static_cast<gbwt::PathName::path_name_type>(thread_name[2]),
            static_cast<gbwt::PathName::path_name_type>(thread_name[3])
        });
    });
    
    // Finish up writing
    gbwt_builder->finish();
    gbwt_builder->index.metadata.setSamples(sample_names);
    gbwt_builder->index.metadata.setHaplotypes(haplotype_count);
    gbwt_builder->index.metadata.setContigs(contig_names);
    if (show_progress) {
        cerr << "GBWT metadata: "; gbwt::operator<<(cerr, gbwt_builder->index.metadata); cerr << endl;
    }
        
    // Extract the built index
    unique_ptr<gbwt::DynamicGBWT> built(new gbwt::DynamicGBWT());
    gbwt_builder->swapIndex(*built);
    
    // Clear the builder
    gbwt_builder.reset();
    
    // Return the finished index
    return built;
}

}

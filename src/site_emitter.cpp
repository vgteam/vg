#include "site_emitter.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>

namespace vg {

SiteEmitter::SiteEmitter(SortedVCFWriter& output) : output(output) {
    // Every record has AT.
    output.add_info_line(
        "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">");
}

bool SiteEmitter::same_as_reference(const CalledSite&, int) {
    return false;
}

int64_t SiteEmitter::phase(const CalledSite&, const vector<int>&, const map<int, int>&, string&) {
    return -1;
}

int SiteEmitter::write_blocks(const PathPositionHandleGraph&, const CalledSite&,
                              const SiteRecord&) {
    return -1;
}

void SiteEmitter::finish_record(const CalledSite&, vcflib::Variant&) {
    // The record is written as emit built it.
}

void SiteEmitter::on_site_done(const CalledSite&, const map<int, int>&, bool) {
    // Only an emitter that tracks how its sites were written needs this.
}

const pair<string, vector<NodeVisit>>& SiteEmitter::spelling(const allele_source_t& alleles,
                                                             Spellings& spelled, size_t trav) {
    auto found = spelled.find(trav);
    if (found == spelled.end()) {
        found = spelled.try_emplace(trav, alleles(trav)).first;
    }
    return found->second;
}

bool SiteEmitter::emit(const PathPositionHandleGraph& graph, const CalledSite& site,
                       const allele_source_t& alleles, const field_writer_t& write_fields) {
    SiteRecord record = build_record(graph, site, alleles, write_fields);
    return file_record(graph, site, record);
}

SiteRecord SiteEmitter::build_record(const PathPositionHandleGraph& graph, const CalledSite& site,
                                     const allele_source_t& alleles,
                                     const field_writer_t& write_fields) {
    SiteRecord record;
    Spellings spelled;
    record.variant.ref = spelling(alleles, spelled, site.ref_trav_idx).first;

    // The VCF allele number of each allele sequence.
    map<string, int> allele_to_gt;
    // The traversal each VCF allele is spelled from, or -1 for "*".
    vector<int> site_trav;
    number_called_alleles(site, alleles, spelled, record, allele_to_gt, site_trav);
    if (genotype_snarls && site_trav.size() <= 1) {
        add_uncalled_alleles(site, alleles, spelled, allele_to_gt, site_trav);
    }

    // The visits of each VCF allele, with none for "*".
    vector<vector<NodeVisit>> allele_visits;
    for (int trav : site_trav) {
        allele_visits.push_back(trav >= 0 ? spelling(alleles, spelled, trav).second
                                          : vector<NodeVisit>());
    }

    record.variant.alt.resize(allele_to_gt.size() - 1);
    record.variant.alleles.resize(allele_to_gt.size());
    record.variant.info["AT"].resize(allele_to_gt.size());
    for (const auto& [allele_string, allele_number] : allele_to_gt) {
        if (allele_number > 0) {
            record.variant.alt[allele_number - 1] = allele_string;
        }
        record.variant.alleles[allele_number] = allele_string;
        add_allele_path_to_info(record.variant, allele_number, allele_visits.at(allele_number),
                                false, output.translation());
    }

    place_record(graph, site, output.translation(), record);
    record.variant.filter = "PASS";
    record.variant.updateAlleleIndexes();

    const string& sample = output.sample_name();
    record.variant.format.emplace_back("GT");
    record.variant.samples[sample]["GT"].push_back(
        !site.genotype.empty() ? unphased_gt_text(record.genotype) : no_call_gt_text(site.ploidy));

    // Phase the genotype here, where `trav_to_allele` maps traversals to the alleles of this
    // record.
    int64_t phase_set = -1;
    if (!site.genotype.empty()) {
        phase_set = phase(site, record.genotype, record.trav_to_allele,
                          record.variant.samples[sample]["GT"][0]);
    }

    write_fields(site_trav, record.genotype, record.variant);

    // Put PS last in FORMAT.
    if (phase_set >= 0) {
        record.variant.format.emplace_back("PS");
        record.variant.samples[sample]["PS"].push_back(std::to_string(phase_set));
    }

    // When genotyping snarls, trim the alleles only as far as the boundary nodes, so that every
    // sample gets the same alleles.
    size_t start_trim = 0;
    size_t end_trim = 0;
    if (genotype_snarls) {
        start_trim = graph.get_length(site.start);
        end_trim = graph.get_length(site.end);
    }
    flatten_common_allele_ends(record.variant, true, end_trim);
    flatten_common_allele_ends(record.variant, false, start_trim);

    // Merge near-identical called ALTs, turning 1/2 into 1/1.
    record.alleles_merged = merge_similar_alleles(graph, allele_visits, record.genotype, sample,
                                                  record.variant, gl_layout,
                                                  allele_merge_threshold, allele_merge_min_len);

    add_star_for_missing(record.variant, record.genotype, sample);

    return record;
}

void SiteEmitter::number_called_alleles(const CalledSite& site, const allele_source_t& alleles,
                                        Spellings& spelled, SiteRecord& record,
                                        map<string, int>& allele_to_gt, vector<int>& site_trav) {
    site_trav = {site.ref_trav_idx};
    allele_to_gt[record.variant.ref] = 0;
    record.trav_to_allele[site.ref_trav_idx] = 0;
    // The number of the "*" allele, once a genotype entry needs it.
    int star_allele_idx = -1;
    for (int trav : site.genotype) {
        if (trav == STAR_ALLELE_MARKER) {
            if (star_allele_idx < 0) {
                star_allele_idx = (int)allele_to_gt.size();
                allele_to_gt["*"] = star_allele_idx;
                site_trav.push_back(-1);
            }
            record.genotype.push_back(star_allele_idx);
        } else if (trav == MISSING_ALLELE_MARKER) {
            record.genotype.push_back(MISSING_ALLELE_MARKER);
        } else if (trav == site.ref_trav_idx || same_as_reference(site, trav)) {
            // A traversal that should be treated as the reference allele gets allele 0.
            record.genotype.push_back(0);
            record.trav_to_allele[trav] = 0;
        } else {
            if (const string& allele_string = spelling(alleles, spelled, trav).first;
                allele_to_gt.contains(allele_string)) {
                record.genotype.push_back(allele_to_gt[allele_string]);
            } else {
                site_trav.push_back(trav);
                record.genotype.push_back((int)allele_to_gt.size());
                allele_to_gt[allele_string] = record.genotype.back();
            }
            record.trav_to_allele[trav] = record.genotype.back();
        }
    }
}

void SiteEmitter::add_uncalled_alleles(const CalledSite& site, const allele_source_t& alleles,
                                       Spellings& spelled, map<string, int>& allele_to_gt,
                                       vector<int>& site_trav) {
    map<string, int> allele_map;
    for (size_t i = 0; i < site.traversal_count; ++i) {
        allele_map.try_emplace(spelling(alleles, spelled, i).first, (int)i);
    }
    // Consider the first `max_uncalled_alleles` sequences, including any already numbered.
    size_t considered = 0;
    for (const auto& [allele_string, trav] : allele_map) {
        if (considered == max_uncalled_alleles) {
            break;
        }
        ++considered;
        if (!allele_to_gt.contains(allele_string)) {
            allele_to_gt[allele_string] = (int)allele_to_gt.size();
            site_trav.push_back(trav);
        }
    }
}

bool SiteEmitter::file_record(const PathPositionHandleGraph& graph, const CalledSite& site,
                              SiteRecord& record) {
    // Let the subclass split the site into several records. A result of -1 means it did not, and
    // the site gets one record.
    if (int block_lines = write_blocks(graph, site, record); block_lines >= 0) {
        // Each of the records numbers its own alleles, so there is no allele map to report.
        on_site_done(site, map<int, int>(), block_lines > 0);
        return block_lines > 0;
    }

    vcflib::Variant& out_variant = record.variant;
    // A record with no ALT is hidden, unless `genotype_snarls` asks for reference calls to be
    // written too.
    bool hidden = !genotype_snarls && out_variant.alt.empty();
    if (!hidden) {
        finish_record(site, out_variant);
    }
    bool added = output.add_variant(out_variant, 0, hidden);
    on_site_done(site, record.trav_to_allele, added && !hidden);
    if (!added) {
        stringstream ss;
        ss << out_variant;
        cerr << "Warning [vg call]: Skipping variant at " << out_variant.sequenceName << ":"
             << out_variant.position << " with ID=" << out_variant.id
             << " because its line length of " << ss.str().length()
             << " exceeds vg's limit of " << SortedVCFWriter::max_vcf_line_length << endl;
    }
    return added;
}

}

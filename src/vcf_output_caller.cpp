#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "vcf_output_caller.hpp"
#include "graph_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"
#include "snarl_traversal_text.hpp"
#include "utility.hpp"

//#define debug

namespace vg {

bool VCFOutputCaller::BufferedRecordKey::operator<(const BufferedRecordKey& other) const {
    return std::tie(contig, position, id, block)
           < std::tie(other.contig, other.position, other.id, other.block);
}

void VCFOutputCaller::add_format_header(string line) {
    format_header_lines.push_back(std::move(line));
}

void VCFOutputCaller::add_info_header(string line) {
    info_header_lines.push_back(std::move(line));
}

const unordered_map<nid_t, pair<string, size_t>>* VCFOutputCaller::get_translation() const {
    return translation;
}

Region record_ref_region(const vcflib::Variant& record) {
    return Region{record.sequenceName, record.position,
                  record.position + (int64_t)record.ref.length() - 1};
}

VCFOutputCaller::VCFOutputCaller(const string& sample_name) : sample_name(sample_name), translation(nullptr), include_nested(false)
{
    output_variants.resize(get_thread_count());
    suppressed_ref_info.resize(get_thread_count());
}

VCFOutputCaller::~VCFOutputCaller() {
}

string VCFOutputCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides) const {
    stringstream ss;
    ss << "##fileformat=VCFv4.2" << endl;    
    for (int i = 0; i < contigs.size(); ++i) {
        const string& contig = contigs[i];
        size_t length;
        if (i < contig_length_overrides.size()) {
            // length override provided
            length = contig_length_overrides[i];
        } else {
            length = 0;
            for (handle_t handle : graph.scan_path(graph.get_path_handle(contig))) {
                length += graph.get_length(handle);
            }
        }
        ss << "##contig=<ID=" << contig << ",length=" << length << ">" << endl;
    }
    if (include_nested) {
        ss << nesting_info_headers();
    }
    for (const string& line : format_header_lines) {
        ss << line << endl;
    }
    ss << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << endl;
    for (const string& line : info_header_lines) {
        ss << line << endl;
    }
    if (allele_merge_threshold < 1.0) {
        ss << "##INFO=<ID=MAT,Number=.,Type=String,Description=\"Merged Allele Traversal: "
           << "ALT alleles merged after genotyping by -L/--cluster, as OLD>NEW:SIMILARITY using "
           << "pre-merge allele numbers. AD and GL are folded onto the surviving allele and MAD is "
           << "recomputed; DP, QUAL, GQ, GP and FILTER are as computed over the pre-merge allele set. "
           << "In a nested run this record gives the collapsed view of the site and its child "
           << "records the precise one, so they disagree by design.\">"
           << endl;
    }
    return ss.str();
}


bool VCFOutputCaller::add_variant(vcflib::Variant& var, size_t block) const {
    var.setVariantCallFile(output_vcf);
    stringstream ss;
    ss << var;
    string dest;
    if (ss.str().length() > VCFOutputCaller::max_vcf_line_length) {
        return false;
    }         
    int ret = zstdutil::CompressString(ss.str(), dest);
    assert(ret == 0);
    // the Variant object is too big to keep in memory when there are many genotypes, so we
    // store it in a zstd-compressed string
    output_variants[omp_get_thread_num()].emplace_back(
        BufferedRecordKey{var.sequenceName, (size_t)var.position, var.id, block}, dest);
    return true;
}

void VCFOutputCaller::write_variants(ostream& out_stream, const SnarlManager* snarl_manager) {
    assert(include_nested == false || snarl_manager != nullptr);
    if (snarl_manager == nullptr) {
        write_variants(out_stream, (const SiteTree*)nullptr);
    } else {
        const SnarlManagerSiteTree sites(*snarl_manager);
        write_variants(out_stream, &sites);
    }
}

void VCFOutputCaller::write_variants(ostream& out_stream, const SiteTree* sites) {
    assert(include_nested == false || sites != nullptr);
    if (include_nested) {
        update_nesting_info_tags(*sites);
    }
    vector<pair<BufferedRecordKey, string>> all_variants;
    // Reserve room for every record at once, rather than growing for each buffer.
    size_t total_variants = 0;
    for (const auto& buf : output_variants) {
        total_variants += buf.size();
    }
    all_variants.reserve(total_variants);
    // Move the records out of each buffer and free it. `buf` is not const, because moving from a
    // const range copies, and the whole VCF is in memory here.
    for (auto& buf : output_variants) {
        std::ranges::move(buf, std::back_inserter(all_variants));
        buf.clear();
        buf.shrink_to_fit();
    }
    std::sort(all_variants.begin(), all_variants.end(),
              [](const pair<BufferedRecordKey, string>& v1,
                 const pair<BufferedRecordKey, string>& v2) {
                  return v1.first < v2.first;
              });

    // Decompress the records in parallel, a batch at a time, and write each batch in order. Only
    // one batch of text is in memory at once.
    const size_t batch_records = 1 << 16;
    vector<string> lines;
    for (size_t batch_start = 0; batch_start < all_variants.size(); batch_start += batch_records) {
        const size_t batch_end = min(all_variants.size(), batch_start + batch_records);
        lines.assign(batch_end - batch_start, string());
#pragma omp parallel for schedule(dynamic, 256)
        for (size_t record_i = batch_start; record_i < batch_end; ++record_i) {
            const string& compressed = all_variants[record_i].second;
            string& dest = lines[record_i - batch_start];
            int ret = zstdutil::DecompressString(compressed, dest);
            assert(ret == 0);
        }
        for (const string& line : lines) {
            // endl would flush, making one write per record.
            out_stream << line << '\n';
        }
    }
}



static int countAlts(vcflib::Variant& var, int alleleIndex) {
    int alts = 0;
    for (auto s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                if (g->first == alleleIndex) {
                    alts += g->second;
                }
            }
        }
    }
    return alts;
}

static int countAlleles(vcflib::Variant& var) {
    int alleles = 0;
    for (auto s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
		if (g->first != vcflib::NULL_ALLELE) {
		    alleles += g->second;
		}
            }
        }
    }
    return alleles;
}

// this isn't from vcflib, but seems to make more sense than just returning the number of samples in the file again and again
static int countSamplesWithData(vcflib::Variant& var) {
    int samples_with_data = 0;
    for (auto s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        bool has_data = false;
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
		if (g->first != vcflib::NULL_ALLELE) {
                    has_data = true;
                    break;
		}
            }
        }
        if (has_data) {
            ++samples_with_data;
        }
    }
    return samples_with_data;
}

void VCFOutputCaller::vcf_fixup(vcflib::Variant& var) const {
    // copied from https://github.com/vgteam/vcflib/blob/master/src/vcffixup.cpp
    
    stringstream ns;
    ns << countSamplesWithData(var);
    var.info["NS"].clear();
    var.info["NS"].push_back(ns.str());

    var.info["AC"].clear();
    var.info["AF"].clear();
    var.info["AN"].clear();

    int allelecount = countAlleles(var);
    stringstream an;
    an << allelecount;
    var.info["AN"].push_back(an.str());

    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
        string& allele = *a;
        int altcount = countAlts(var, var.getAltAlleleIndex(allele) + 1);
        stringstream ac;
        ac << altcount;
        var.info["AC"].push_back(ac.str());
        stringstream af;
        double faf = (double) altcount / (double) allelecount;
        if(faf != faf) faf = 0;
        af << faf;
        var.info["AF"].push_back(af.str());
    }
}

void VCFOutputCaller::set_translation(const unordered_map<nid_t, pair<string, size_t>>* translation) {
    this->translation = translation;
}

void VCFOutputCaller::set_nested(bool nested) {
    include_nested = nested;
}

void VCFOutputCaller::set_gref_levels(map<string, int> levels) {
    this->gref_levels = std::move(levels);
}

void VCFOutputCaller::set_allele_merge(double threshold, int64_t min_len) {
    allele_merge_threshold = threshold;
    allele_merge_min_len = min_len;
}

unordered_set<string> VCFOutputCaller::get_output_contigs() const {
    unordered_set<string> contigs;
    // The sort key is (sequenceName, position) (see add_variant), so the contig is right
    // there and nothing has to be decompressed.
    for (const auto& thread_buf : output_variants) {
        for (const auto& [key, record] : thread_buf) {
            contigs.insert(key.contig);
        }
    }
    return contigs;
}

string VCFOutputCaller::prune_header_contigs(const string& header,
                                             const unordered_set<string>& keep) const {
    static const string contig_prefix = "##contig=<ID=";
    stringstream pruned;
    vector<string> lines = split_delims(header, "\n");
    for (const string& line : lines) {
        if (line.compare(0, contig_prefix.size(), contig_prefix) == 0) {
            // Parse the ID back out the same way it was written, rather than scanning for a
            // delimiter: contig names are path names and nothing stops one containing ',' or
            // '>'.  Both producers emit exactly ##contig=<ID=NAME,length=N> -- vcf_header()
            // above and Deconstructor::add_contigs_to_vcf_header().
            static const string contig_suffix = ",length=";
            size_t id_start = contig_prefix.size();
            size_t id_end = line.rfind(contig_suffix);
            if (id_end == string::npos || id_end < id_start) {
                // not a shape we wrote; leave it alone rather than guess
                pruned << line << "\n";
                continue;
            }
            string id = line.substr(id_start, id_end - id_start);
            if (!keep.count(id)) {
                continue;
            }
        }
        pruned << line << "\n";
    }
    string result = pruned.str();
    if (!header.empty() && header.back() != '\n' && !result.empty()) {
        // input had no trailing newline, so don't invent one
        result.pop_back();
    }
    return result;
}




/// The writer of a site called by a SnarlCaller, whose traversals are SnarlTraversals and whose
/// INFO and FORMAT fields the SnarlCaller writes.
class VCFOutputCaller::SnarlTraversalEmitter : public VCFOutputCaller::SiteEmitter {
public:
    /// Make an emitter for the record of `snarl`, which reads all the arguments by reference, so
    /// they must outlive it.
    ///
    /// @param trav_to_string Spells an allele. When null, an allele is spelled as the sequence of
    /// its traversal.
    SnarlTraversalEmitter(VCFOutputCaller& output, const HandleGraph& graph, SnarlCaller& snarl_caller,
                          const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                          const vector<int>& genotype, int ref_trav_idx,
                          const unique_ptr<SnarlCaller::CallInfo>& call_info,
                          bool genotype_snarls,
                          const function<string(const vector<SnarlTraversal>&, const vector<int>&,
                                                int, int, int)>& trav_to_string) :
        SiteEmitter(output, genotype_snarls), graph(graph), snarl_caller(snarl_caller), snarl(snarl),
        called_traversals(called_traversals), genotype(genotype), ref_trav_idx(ref_trav_idx),
        call_info(call_info), trav_to_string(trav_to_string) {
    }

protected:
    string spell(int trav, int genotype_index) override {
        if (trav_to_string == nullptr) {
            return vg::trav_string(graph, called_traversals[trav], output.translation);
        }
        return trav_to_string(called_traversals, genotype, trav, genotype_index, ref_trav_idx);
    }

    vector<NodeVisit> visits(int trav) override {
        return visits_of(called_traversals[trav]);
    }

    void fill_info(const vector<int>& site_trav, const vector<int>& site_genotype,
                   vcflib::Variant& variant) override {
        // A negative traversal number is the "*" allele, which gets an empty traversal.
        vector<SnarlTraversal> site_traversals;
        site_traversals.reserve(site_trav.size());
        for (int trav : site_trav) {
            site_traversals.push_back(trav >= 0 ? called_traversals[trav] : SnarlTraversal());
        }
        snarl_caller.update_vcf_info(snarl, site_traversals, site_genotype, call_info,
                                     output.sample_name, variant);
    }

private:
    const HandleGraph& graph;
    SnarlCaller& snarl_caller;
    const Snarl& snarl;
    const vector<SnarlTraversal>& called_traversals;
    const vector<int>& genotype;
    int ref_trav_idx;
    const unique_ptr<SnarlCaller::CallInfo>& call_info;
    const function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)>&
        trav_to_string;
};

bool VCFOutputCaller::emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                                   const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                                   const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                                   const string& ref_path_name, int64_t ref_offset, bool genotype_snarls, int ploidy,
                                   const function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)>& trav_to_string) {
    SnarlTraversalEmitter emitter(*this, graph, snarl_caller, snarl, called_traversals, genotype,
                                  ref_trav_idx, call_info, genotype_snarls, trav_to_string);
    CalledSite site{
        .start = graph.get_handle(snarl.start().node_id(), snarl.start().backward()),
        .end = graph.get_handle(snarl.end().node_id(), snarl.end().backward()),
        .ref_path_name = ref_path_name,
        .ref_offset = ref_offset,
        .genotype = genotype,
        .ref_trav_idx = ref_trav_idx,
        .traversal_count = called_traversals.size(),
        .ploidy = ploidy,
    };
    return emitter.emit(graph, site);
}

VCFOutputCaller::SiteEmitter::SiteEmitter(VCFOutputCaller& output, bool genotype_snarls) :
    output(output), genotype_snarls(genotype_snarls) {
}

bool VCFOutputCaller::SiteEmitter::same_as_reference(int, int) {
    return false;
}

void VCFOutputCaller::SiteEmitter::count_site(const vector<int>&, int) {
    // Only a caller that keeps counts of its sites needs this.
}

int64_t VCFOutputCaller::SiteEmitter::phase(const vector<int>&, const map<int, int>&, string&) {
    return -1;
}

GLLayout VCFOutputCaller::SiteEmitter::gl_layout() {
    return GLLayout::IMajor;
}

int VCFOutputCaller::SiteEmitter::write_blocks(const PathPositionHandleGraph&, const vector<int>&,
                                               int, const SiteRecord&, GLLayout) {
    return -1;
}

void VCFOutputCaller::SiteEmitter::finish_record(vcflib::Variant&) {
    // The record is finished as build_record left it.
}

void VCFOutputCaller::SiteEmitter::fill_info(const vector<int>&, const vector<int>&,
                                             vcflib::Variant&) {
    // The record gets no INFO or FORMAT fields beyond the ones build_record writes.
}

void VCFOutputCaller::SiteEmitter::site_filed(const map<int, int>&, size_t, bool) {
    // Only a caller that tracks how its sites were written needs this.
}

bool VCFOutputCaller::SiteEmitter::emit(const PathPositionHandleGraph& graph,
                                        const CalledSite& site) {
    count_site(site.genotype, site.ref_trav_idx);
    GLLayout layout = gl_layout();
    SiteRecord record = build_record(graph, site, layout);
    return file_record(graph, record, site, layout);
}

SiteRecord VCFOutputCaller::SiteEmitter::build_record(const PathPositionHandleGraph& graph,
                                                      const CalledSite& site, GLLayout layout) {
    SiteRecord record;

    // Spell REF at a genotype entry that holds the reference traversal, if there is one, since
    // spelling can depend on the entry.
    auto ref_gt_it = std::ranges::find(site.genotype, site.ref_trav_idx);
    record.variant.ref = spell(site.ref_trav_idx, ref_gt_it != site.genotype.end()
                                                      ? ref_gt_it - site.genotype.begin()
                                                      : 0);

    // The VCF allele number of each allele sequence.
    map<string, int> allele_to_gt;
    // The traversal each VCF allele is spelled from, or -1 for "*".
    vector<int> site_trav;
    number_called_alleles(site, record, allele_to_gt, site_trav);
    if (genotype_snarls && site_trav.size() <= 1) {
        add_uncalled_alleles(site, allele_to_gt, site_trav);
    }

    // The visits of each VCF allele, with none for "*".
    vector<vector<NodeVisit>> allele_visits;
    for (int trav : site_trav) {
        allele_visits.push_back(trav >= 0 ? visits(trav) : vector<NodeVisit>());
    }

    record.variant.alt.resize(allele_to_gt.size() - 1);
    record.variant.alleles.resize(allele_to_gt.size());
    record.variant.info["AT"].resize(allele_to_gt.size());
    for (const auto& [allele_string, allele_number] : allele_to_gt) {
        if (allele_number > 0) {
            record.variant.alt[allele_number - 1] = allele_string;
        }
        record.variant.alleles[allele_number] = allele_string;
        vg::add_allele_path_to_info(record.variant, allele_number, allele_visits.at(allele_number),
                                    false, output.translation);
    }

    place_record(graph, site, output.translation, record);
    record.variant.filter = "PASS";
    record.variant.updateAlleleIndexes();

    record.variant.format.emplace_back("GT");
    record.variant.samples[output.sample_name]["GT"].push_back(
        !site.genotype.empty() ? unphased_gt_text(record.genotype) : no_call_gt_text(site.ploidy));

    // Phase the genotype here, where `trav_to_allele` maps traversals to the alleles of this
    // record.
    int64_t phase_set = -1;
    if (!site.genotype.empty()) {
        phase_set = phase(record.genotype, record.trav_to_allele,
                          record.variant.samples[output.sample_name]["GT"][0]);
    }

    fill_info(site_trav, record.genotype, record.variant);

    // Put PS last in FORMAT.
    if (phase_set >= 0) {
        record.variant.format.emplace_back("PS");
        record.variant.samples[output.sample_name]["PS"].push_back(std::to_string(phase_set));
    }

    // When genotyping snarls, trim the alleles only as far as the boundary nodes, so that every
    // sample gets the same alleles.
    size_t start_trim = 0;
    size_t end_trim = 0;
    if (genotype_snarls) {
        start_trim = graph.get_length(site.start);
        end_trim = graph.get_length(site.end);
    }
    vg::flatten_common_allele_ends(record.variant, true, end_trim);
    vg::flatten_common_allele_ends(record.variant, false, start_trim);

    // Merge near-identical called ALTs, turning 1/2 into 1/1.
    record.alleles_merged = merge_similar_alleles(
        graph, allele_visits, record.genotype, output.sample_name, record.variant, layout,
        output.allele_merge_threshold, output.allele_merge_min_len);

    add_star_for_missing(record.variant, record.genotype, output.sample_name);

    return record;
}

void VCFOutputCaller::SiteEmitter::number_called_alleles(const CalledSite& site, SiteRecord& record,
                                                         map<string, int>& allele_to_gt,
                                                         vector<int>& site_trav) {
    site_trav = {site.ref_trav_idx};
    allele_to_gt[record.variant.ref] = 0;
    record.trav_to_allele[site.ref_trav_idx] = 0;
    // The number of the "*" allele, once a genotype entry needs it.
    int star_allele_idx = -1;
    for (int i = 0; i < site.genotype.size(); ++i) {
        if (site.genotype[i] == STAR_ALLELE_MARKER) {
            if (star_allele_idx < 0) {
                star_allele_idx = allele_to_gt.size();
                allele_to_gt["*"] = star_allele_idx;
                site_trav.push_back(-1);
            }
            record.genotype.push_back(star_allele_idx);
        } else if (site.genotype[i] == MISSING_ALLELE_MARKER) {
            record.genotype.push_back(MISSING_ALLELE_MARKER);
        } else if (site.genotype[i] == site.ref_trav_idx
                   || same_as_reference(site.genotype[i], site.ref_trav_idx)) {
            // A traversal that should be treated as the reference allele gets allele 0.
            record.genotype.push_back(0);
            record.trav_to_allele[site.genotype[i]] = 0;
        } else {
            if (string allele_string = spell(site.genotype[i], i);
                allele_to_gt.contains(allele_string)) {
                record.genotype.push_back(allele_to_gt[allele_string]);
            } else {
                site_trav.push_back(site.genotype[i]);
                record.genotype.push_back(allele_to_gt.size());
                allele_to_gt[allele_string] = record.genotype.back();
            }
            record.trav_to_allele[site.genotype[i]] = record.genotype.back();
        }
    }
}

void VCFOutputCaller::SiteEmitter::add_uncalled_alleles(const CalledSite& site,
                                                        map<string, int>& allele_to_gt,
                                                        vector<int>& site_trav) {
    map<string, int> allele_map;
    for (int i = 0; i < site.traversal_count; ++i) {
        // TODO: Is this the right genotype index? It only affects alleles that were not called.
        string allele_string = spell(i, max(0, (int)site.genotype.size() - 1));
        if (!allele_map.contains(allele_string)) {
            allele_map[allele_string] = i;
        }
    }
    // Consider the first `max_uncalled_alleles` sequences, including any already numbered.
    size_t considered = 0;
    for (const auto& [allele_string, trav] : allele_map) {
        if (considered == output.max_uncalled_alleles) {
            break;
        }
        ++considered;
        if (!allele_to_gt.contains(allele_string)) {
            allele_to_gt[allele_string] = allele_to_gt.size();
            site_trav.push_back(trav);
        }
    }
}

bool VCFOutputCaller::SiteEmitter::file_record(const PathPositionHandleGraph& graph,
                                               SiteRecord& record, const CalledSite& site,
                                               GLLayout gl_layout) {
    vcflib::Variant& out_variant = record.variant;

    // Let the subclass split the site into several records. A result of -1 means it did not, and
    // the site gets one record.
    if (const int block_lines = write_blocks(graph, site.genotype, site.ref_trav_idx, record,
                                             gl_layout);
        block_lines >= 0) {
        // Each of the records numbers its own alleles, so there is no allele map to report.
        site_filed(map<int, int>(), 0, block_lines > 0);
        return block_lines > 0;
    }

    // Whether the site gets a line. A record with no ALT gets one only if `genotype_snarls` asks
    // for reference calls to be written too.
    const bool wants_line = genotype_snarls || !out_variant.alt.empty();
    bool added = false;
    if (wants_line) {
        finish_record(out_variant);
        added = output.add_variant(out_variant);
    } else if (output.include_nested) {
        // Keep the reference interval of a site with no line, for the RC, RS and RD of its
        // children.
        output.suppressed_ref_info[omp_get_thread_num()][out_variant.id] =
            record_ref_region(out_variant);
    }
    site_filed(record.trav_to_allele, site.traversal_count, added);
    if (wants_line && !added) {
        stringstream ss;
        ss << out_variant;
        cerr << "Warning [vg call]: Skipping variant at " << out_variant.sequenceName << ":" << out_variant.position
             << " with ID=" << out_variant.id << " because its line length of " << ss.str().length() << " exceeds vg's limit of "
             << VCFOutputCaller::max_vcf_line_length << endl;
    }
    // Return false only when add_variant refused a line the site wanted.
    return wants_line ? added : true;
}

string VCFOutputCaller::snarl_caller_vcf_header(const PathHandleGraph& graph,
                                                const vector<string>& contigs,
                                                const vector<size_t>& contig_length_overrides,
                                                const SnarlCaller& snarl_caller) const {
    return genotyper_vcf_header(graph, contigs, contig_length_overrides,
                                [&](string& header) { snarl_caller.update_vcf_header(header); });
}

string VCFOutputCaller::genotyper_vcf_header(const PathHandleGraph& graph,
                                             const vector<string>& contigs,
                                             const vector<size_t>& contig_length_overrides,
                                             const function<void(string& header)>& add_genotyper_lines) const {
    string header = VCFOutputCaller::vcf_header(graph, contigs, contig_length_overrides);
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    add_genotyper_lines(header);
    header += "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    header += "##SAMPLE=<ID=" + sample_name + ">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name;
    assert(output_vcf.openForOutput(header));
    header += "\n";
    return header;
}

string VCFOutputCaller::nesting_info_headers() {
    stringstream ss;
    ss << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree counting only ancestors whose record is on this record's own reference contig (0=top level for this contig)\">" << endl;
    ss << "##INFO=<ID=CH,Number=1,Type=Integer,Description=\"Nesting steps between VCF reference contigs: how many coordinate-system changes separate this record from a linear reference. Counted as the greater of the in-VCF ancestor hops and the record's own gref contig level, because counting only ancestors that happened to emit a record made a record on a gref fragment whose parent produced no line indistinguishable from one on the linear reference -- 29,843 of 41,669 off-reference records on a gref-covered chr20. So CH >= 1 no longer implies an in-VCF parent, and therefore no longer implies PS\">" << endl;
    ss << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << endl;
    ss << "##INFO=<ID=RC,Number=1,Type=String,Description=\"CHROM of the topmost ancestor record in this VCF, or this record's own CHROM when it has none. On a gref fragment, where that own CHROM would be no use, the enclosing site is named even if it produced no record of its own; the tags are absent when there is no such site either\">" << endl;
    ss << "##INFO=<ID=RS,Number=1,Type=Integer,Description=\"Start of the site named by RC: the POS of its record, or where the site begins when it produced none. A position on that contig, not a span of the snarl, so it can precede the sequence this record describes\">" << endl;
    ss << "##INFO=<ID=RD,Number=1,Type=Integer,Description=\"End of the site named by RC: RS plus the length of that site's REF allele\">" << endl;
    return ss.str();
}


namespace {

/// Table from the boundary visits of a site to the name of its record in the VCF.
using NameOfEnds = unordered_map<SiteEnds, const string*>;

/// Return the name of a site's record: its own ID, or else its boundary visits as `site_name`
/// prints them.
string record_name_of_site(const SiteTree& sites, SiteTree::site_t site,
                           const NodeTranslation* translation) {
    if (const string* id = sites.id_of(site)) {
        return *id;
    }
    const SiteEnds e = sites.ends_of(site);
    return site_name(e.start_id, e.start_backward, e.end_id, e.end_backward, translation, false);
}

/// Return the name the record of a site has when the site is read in the other orientation.
///
/// A site with an ID has that name in both orientations.
string flipped_record_name_of_site(const SiteTree& sites, SiteTree::site_t site,
                                   const NodeTranslation* translation) {
    if (const string* id = sites.id_of(site)) {
        return *id;
    }
    const SiteEnds e = sites.ends_of(site);
    return site_name(e.end_id, !e.end_backward, e.start_id, !e.start_backward, translation, false);
}

/// Read the boundary visits of a site back from an untranslated record name like ">12<17".
///
/// Return false for a name that printing the visits read from it would not give back, so that a
/// name and a pair of visits correspond one to one.
bool read_record_name_ends(const string& name, SiteEnds& ends) {
    if (name.size() < 4 || (name[0] != '<' && name[0] != '>')) {
        return false;
    }
    const size_t middle = name.find_first_of("<>", 1);
    if (middle == string::npos || middle < 2 || middle + 1 >= name.size()) {
        return false;
    }
    const char* text = name.data();
    auto start = std::from_chars(text + 1, text + middle, ends.start_id);
    auto end = std::from_chars(text + middle + 1, text + name.size(), ends.end_id);
    if (start.ec != std::errc() || start.ptr != text + middle
        || end.ec != std::errc() || end.ptr != text + name.size()) {
        return false;
    }
    ends.start_backward = name[0] == '<';
    ends.end_backward = name[middle] == '<';
    return site_name(ends.start_id, ends.start_backward, ends.end_id, ends.end_backward, nullptr,
                     false) == name;
}

/// Call `match` on each record name in `chrom_of_name` that an untranslated site has.
///
/// A site with an ID matches by its ID. Any other site matches by its boundary visits, first in
/// its own orientation and then flipped, since a caller can write a site either way round.
template<typename Match>
void for_each_record_name_of_site(const SiteTree& sites, SiteTree::site_t site,
                                  const unordered_map<string, uint32_t>& chrom_of_name,
                                  const NameOfEnds& name_of_ends, Match&& match) {
    if (const string* id = sites.id_of(site)) {
        if (auto named = chrom_of_name.find(*id); named != chrom_of_name.end()) {
            match(named->first);
        }
        return;
    }
    const SiteEnds e = sites.ends_of(site);
    if (auto own = name_of_ends.find(e); own != name_of_ends.end()) {
        match(*own->second);
    }
    if (auto flipped = name_of_ends.find(SiteEnds{e.end_id, !e.end_backward, e.start_id,
                                                  !e.start_backward});
        flipped != name_of_ends.end()) {
        match(*flipped->second);
    }
}

/// Map each record name in `chrom_of_name` to the untranslated site that has it.
///
/// Each name is read back into boundary visits once, and the sites are matched by their visits,
/// which is much faster than printing both names of every site. Where two sites have one name,
/// the name goes to the one that the preorder visits last.
unordered_map<string, SiteTree::site_t>
index_untranslated_sites_by_record_name(const SiteTree& sites,
                                        const unordered_map<string, uint32_t>& chrom_of_name) {
    unordered_map<string, SiteTree::site_t> name_to_snarl;
    name_to_snarl.reserve(chrom_of_name.size());
    NameOfEnds name_of_ends;
    name_of_ends.reserve(chrom_of_name.size());
    for (const auto& [name, chrom] : chrom_of_name) {
        if (SiteEnds ends; read_record_name_ends(name, ends)) {
            name_of_ends.try_emplace(ends, &name);
        }
    }
    // Match the sites on several threads, each into a list of its own. The matches depend on the
    // order of the sites only if two different sites match one name, and then we match again in
    // preorder.
    vector<vector<pair<const string*, SiteTree::site_t>>> found(max(1, omp_get_max_threads()));
    sites.for_each_site([&](SiteTree::site_t snarl) {
        auto& mine = found[omp_get_thread_num()];
        for_each_record_name_of_site(sites, snarl, chrom_of_name, name_of_ends,
                                     [&](const string& name) {
            mine.emplace_back(&name, snarl);
        });
    }, false);
    bool ambiguous = false;
    for (const auto& thread_found : found) {
        for (const auto& [name, snarl] : thread_found) {
            auto [placed, inserted] = name_to_snarl.try_emplace(*name, snarl);
            if (!inserted && placed->second != snarl) {
                ambiguous = true;
            }
        }
    }
    if (ambiguous) {
        name_to_snarl.clear();
        sites.for_each_site([&](SiteTree::site_t snarl) {
            for_each_record_name_of_site(sites, snarl, chrom_of_name, name_of_ends,
                                         [&](const string& name) {
                name_to_snarl[name] = snarl;
            });
        }, true);
    }
    return name_to_snarl;
}

/// Map each record name in `chrom_of_name` to the site that has it, under a node translation.
///
/// Translated names are not node IDs, so each site's names are printed and looked up.
unordered_map<string, SiteTree::site_t>
index_translated_sites_by_record_name(const SiteTree& sites,
                                      const unordered_map<string, uint32_t>& chrom_of_name,
                                      const NodeTranslation* translation) {
    unordered_map<string, SiteTree::site_t> name_to_snarl;
    name_to_snarl.reserve(chrom_of_name.size());
    sites.for_each_site([&](SiteTree::site_t snarl) {
        if (string snarl_name = record_name_of_site(sites, snarl, translation);
            chrom_of_name.contains(snarl_name)) {
            name_to_snarl[std::move(snarl_name)] = snarl;
        }
        if (string flipped_name = flipped_record_name_of_site(sites, snarl, translation);
            chrom_of_name.contains(flipped_name)) {
            name_to_snarl[std::move(flipped_name)] = snarl;
        }
    }, true);
    return name_to_snarl;
}

/// Read the name of the site of each buffered record, and the stretch of the reference its REF
/// covers, decompressing the buffers in parallel. The name is the same for every record of a site
/// split into several.
///
/// The result has one list per buffer, in the order of the buffer.
vector<vector<pair<string, Region>>> read_record_regions(
    const vector<vector<pair<VCFOutputCaller::BufferedRecordKey, string>>>& buffers) {
    vector<vector<pair<string, Region>>> placed(buffers.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < buffers.size(); ++b) {
        placed[b].reserve(buffers[b].size());
        string output_variant_string;
        for (const auto& [key, compressed] : buffers[b]) {
            output_variant_string.clear();
            int ret = zstdutil::DecompressString(compressed, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 5);
            int64_t position = stoll(toks[1]);
            placed[b].emplace_back(block_site_name(toks[2]),
                                   Region{std::move(toks[0]), position,
                                          position + (int64_t)toks[3].length() - 1});
        }
    }
    return placed;
}

}

void VCFOutputCaller::update_nesting_info_tags(const SiteTree& sites) {

    auto name_of = [&](SiteTree::site_t site) {
        return record_name_of_site(sites, site, translation);
    };
    auto flipped_name_of = [&](SiteTree::site_t site) {
        return flipped_record_name_of_site(sites, site, translation);
    };

    // Merge the per-thread suppressed-site intervals collected during calling.  These are sites
    // that never reached the VCF, so pass 1 below cannot see them, but a record nested under one
    // has no other way to name a reference position.
    unordered_map<string, Region> suppressed_ref;
    for (auto& buf : suppressed_ref_info) {
        for (auto& kv : buf) {
            suppressed_ref.emplace(kv.first, std::move(kv.second));
        }
        buf.clear();
        buf.rehash(0);
    }

    // pass 1) index sites in vcf
    // (todo: this could be done more quickly upstream)
    //
    // One index, not two: presence in chrom_of_name IS "this snarl name is in the VCF", and
    // the value is which reference contig its record landed on.  Keeping a separate
    // names_in_vcf set alongside would store all 400k snarl-ID strings twice, which measured
    // as +70 MB of peak RSS on chr22 -- the keys, not the values, are what costs.
    //
    // Contig names are interned rather than stored per record for the same reason: there are
    // at most a few thousand distinct ones, and all we ever ask is whether two are the same.
    unordered_map<string, uint32_t> chrom_index;
    // Whether each interned contig is a synthetic gref fragment, by the same index.  Stored as a
    // bit per contig rather than looked up by name later, so the names are still stored once.
    // is_gref_name(), not is_gref_derived(): a gref copy of a real reference contig is a perfectly
    // good coordinate system -- it is what the whole VCF is deconstructed against -- and only the
    // "_<N>_alt" fragments are positions a reader cannot look up.
    vector<bool> chrom_is_gref_fragment;
    auto intern_chrom = [&](const string& chrom) -> uint32_t {
        auto result = chrom_index.emplace(chrom, (uint32_t)chrom_index.size());
        if (result.second) {
            chrom_is_gref_fragment.push_back(GrefCover::is_gref_name(chrom));
        }
        return result.first->second;
    };
    // One entry per snarl name.  A snarl ID is not unique -- a cyclic reference path that
    // traverses the same snarl twice emits two records with the same ID (see
    // nesting/cyclic_ref_multiple_variants.gfa) -- but both occurrences are traversals of one
    // path, so they are on the same contig and it does not matter which one wins here.
    unordered_map<string, uint32_t> chrom_of_name;
    // The site name and reference stretch of each record, decompressed once, for the indexes and
    // the second pass.
    vector<vector<pair<string, Region>>> placed = read_record_regions(output_variants);
    // Fill the indexes in record order.
    for (const vector<pair<string, Region>>& buffer : placed) {
        for (const auto& [site_name, region] : buffer) {
            chrom_of_name.emplace(site_name, intern_chrom(region.seq));
        }
    }

    // Index the sites by the names of their records. Only names in the VCF are looked up, so only
    // they are indexed.
    unordered_map<string, SiteTree::site_t> name_to_snarl =
        translation == nullptr
            ? index_untranslated_sites_by_record_name(sites, chrom_of_name)
            : index_translated_sites_by_record_name(sites, chrom_of_name, translation);

    // pass 2) identify top-level snarls (those with no ancestors in VCF)
    // and store reference info only for them
    // Keyed by snarl name, then by (chrom, pos), because a snarl ID can carry more than one
    // record: a cyclic reference emits two, both with the same ID (see
    // nesting/cyclic_ref_multiple_variants.gfa, which gives two <5<1 records at POS 20 and 44).
    // A plain name -> interval map was last-write-wins, so both records were handed the
    // surviving one's interval and the record at POS 20 reported RS=44.  The inner map is
    // ordered so that picking begin() is deterministic regardless of thread scheduling.
    //
    // The inner map gives the end of each record's reference stretch by its contig and start.
    unordered_map<string, map<pair<string, int64_t>, int64_t>> top_level_ref_info;

    // Helper to check if a snarl is top-level (no ancestors in VCF)
    auto is_top_level = [&](const string& name) {
        auto it = name_to_snarl.find(name);
        if (it == name_to_snarl.end()) return true; // not found, treat as top-level
        SiteTree::site_t snarl = it->second;
        while ((snarl = sites.parent_of(snarl))) {
            string cur_name = name_of(snarl);
            string flipped_name = flipped_name_of(snarl);
            if (chrom_of_name.contains(cur_name) || chrom_of_name.contains(flipped_name)) {
                return false; // has ancestor in VCF
            }
        }
        return true; // no ancestors in VCF
    };

    // Find the records whose sites are top level, in parallel. Then store their reference
    // intervals in record order.
    // Whether the site of each record is top level, by buffer.
    vector<vector<bool>> top_level(placed.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < placed.size(); ++b) {
        for (const auto& [site_name, region] : placed[b]) {
            top_level[b].push_back(is_top_level(site_name));
        }
    }
    for (size_t b = 0; b < placed.size(); ++b) {
        for (size_t i = 0; i < placed[b].size(); ++i) {
            if (top_level[b][i]) {
                const auto& [site_name, region] = placed[b][i];
                top_level_ref_info[site_name][make_pair(region.seq, region.start)] = region.end;
            }
        }
    }
    vector<vector<pair<string, Region>>>().swap(placed);

    // determine the tags from the index
    //
    // There are exactly two ways a snarl can nest inside its parent's record, and they need
    // to be told apart.  A site inside a *deletion* is covered by its parent contig's own
    // reference allele, so it has coordinates on that contig and its record's CHROM is the
    // same.  A site inside an *insertion* has no path of the parent's contig through it at
    // all, so it is only callable once some other reference (a gref fragment) covers the
    // inserted allele -- and its record's CHROM is therefore different.  So:
    //
    //   contig_level    ancestors whose record is on this record's own CHROM, i.e. how deep
    //                   the site is in its own coordinate system
    //   contig_hops     steps in the chain where CHROM changed, i.e. how many insertions deep
    //                   the site is
    //
    // Returns: (contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
    //           suppressed_name)
    // ref_chrom_name is the topmost ancestor in the VCF that sits on a reference contig rather
    // than a gref one, and suppressed_name the topmost ancestor that was dropped for having no
    // variant.  Both feed the RC/RS/RD choice below; neither affects LV/CH/PS.
    function<tuple<size_t, size_t, string, string, string, string>(const string&, const string&)> get_nesting_tags =
        [&](const string& name, const string& my_chrom) {
        string parent_name;
        string ref_chrom_name;
        string suppressed_name;
        string top_level_name = name;  // default to self (for the top-level case)
        size_t contig_level = 0;
        size_t contig_hops = 0;
        // Our own contig, and the contig of the previously visited link in the chain.
        // chrom_index is complete after pass 1, so this lookup always hits.
        uint32_t my_chrom_id = chrom_index.at(my_chrom);
        uint32_t prev_chrom_id = my_chrom_id;
        SiteTree::site_t snarl = name_to_snarl.at(name);

        assert(snarl != nullptr);
        // walk up the snarl tree
        while ((snarl = sites.parent_of(snarl))) {
            string cur_name = name_of(snarl);

            // Since it is possible that the snarl is actually flipped in the vcf, check for the
            // flipped version too
            string flipped_name = flipped_name_of(snarl);
            const string* hit = nullptr;
            if (chrom_of_name.contains(cur_name)) {
                // only count snarls that are in the vcf
                hit = &cur_name;
            } else if (chrom_of_name.contains(flipped_name)) {
                // snarl is in vcf under flipped orientation
                hit = &flipped_name;
            }
            if (hit == nullptr) {
                // Not in the VCF.  If it was dropped for having no variant we still know where it
                // sits, and it may be the only ancestor that can give a reference position.
                auto sup_it = suppressed_ref.find(cur_name);
                if (sup_it == suppressed_ref.end()) {
                    sup_it = suppressed_ref.find(flipped_name);
                }
                if (sup_it != suppressed_ref.end()) {
                    suppressed_name = sup_it->first;
                }
                continue;
            }

            auto chrom_it = chrom_of_name.find(*hit);
            uint32_t anc_chrom_id = chrom_it == chrom_of_name.end() ? my_chrom_id
                                                                    : chrom_it->second;
            if (anc_chrom_id == my_chrom_id) {
                ++contig_level;
            }
            if (anc_chrom_id != prev_chrom_id) {
                ++contig_hops;
            }
            prev_chrom_id = anc_chrom_id;

            if (parent_name.empty()) {
                // remember the first parent
                parent_name = *hit;
            }
            // keep updating top_level to find the topmost ancestor in VCF
            top_level_name = *hit;
            // ...and, separately, the topmost one actually on a reference contig.  An ancestor on
            // another gref contig can name a position, but not one a reader can look up in the
            // reference, so it is the weaker answer of the two.
            if (!chrom_is_gref_fragment[anc_chrom_id]) {
                ref_chrom_name = *hit;
            }
        }
        return make_tuple(contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
                          suppressed_name);
    };

    // pass 3) add the LV, PS, RC, RS, RD tags
#pragma omp parallel for
    for (uint64_t i = 0; i < output_variants.size(); ++i) {
        auto& thread_buf = output_variants[i];
        for (auto& output_variant_record : thread_buf) {
            string output_variant_string;
            int ret = zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 9);
            // Look the record up by the name of its site, so that a block record gets the tags of
            // its site.
            const string name = block_site_name(toks[2]);

            auto [contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
                  suppressed_name] = get_nesting_tags(name, toks[0]);
            // LV is the level within this record's own reference contig, so that a gRef fragment's
            // records start at level 0 on their own contig.
            //
            // CH counts the ancestors with a record, which is 0 on a gRef fragment whose enclosing
            // site has no line. So we use the gRef level of the contig as a floor for CH. Then
            // CH >= 1 does not imply a parent record, or PS.
            size_t gref_level = 0;
            {
                auto it = gref_levels.find(toks[0]);
                if (it != gref_levels.end() && it->second > 0) {
                    gref_level = (size_t)it->second;
                }
            }
            string nesting_tags = ";LV=" + std::to_string(contig_level);
            nesting_tags += ";CH=" + std::to_string(max(contig_hops, gref_level));
            if (!parent_name.empty()) {
                // Not "if (lv != 0)": those were equivalent only while LV was the absolute
                // count.  A record can now legitimately be at LV=0 and still have a parent on
                // another contig, and it must keep PS -- vcfbub's rescue of the children of
                // popped bubbles is keyed on it.
                nesting_tags += ";PS=" + parent_name;
            }

            // Add RC, RS, RD tags: where to look this record up in the reference.
            //
            // Prefer, in order, the topmost ancestor in the VCF that is on a reference contig;
            // then the topmost ancestor in the VCF at all; then the topmost ancestor that was
            // dropped for having no variant.  The last is what rescues a gref fragment whose
            // parent snarl only the reference and its own gref copy span: the site is real and
            // has a reference interval, it just had nothing to report.
            //
            // If none of those exist there is no reference position to give, and the tags are
            // left off.  They used to fall back to this record's own contig and position, which
            // is not a reference coordinate at all -- on a gref contig it is a self-reference
            // that a reader cannot tell apart from the genuine case.
            const string* ref_source = nullptr;
            if (!ref_chrom_name.empty()) {
                ref_source = &ref_chrom_name;
            } else if (top_level_name != name) {
                ref_source = &top_level_name;
            }
            bool have_ref = true;
            // The stretch of the reference to report, 1-based and end-inclusive.
            Region top_ref;
            if (ref_source == nullptr) {
                if (!GrefCover::is_gref_name(toks[0])) {
                    // Not on a gref fragment, so our own interval is already a position a reader
                    // can look up, and it is the narrower answer of the two.  Keep it rather than
                    // reach for an enclosing site that produced no record: doing that would
                    // repoint every such record on a reference contig at a site LV and CH say it
                    // has no ancestor in.  The records that need the reach are the fragments
                    // below, which have no usable coordinate of their own.
                    int64_t position = stoll(toks[1]);
                    top_ref = {toks[0], position, position + (int64_t)toks[3].length() - 1};
                } else {
                    auto sup_it = suppressed_name.empty() ? suppressed_ref.end()
                                                          : suppressed_ref.find(suppressed_name);
                    if (sup_it != suppressed_ref.end()) {
                        top_ref = sup_it->second;
                    } else {
                        have_ref = false;
                    }
                }
            } else {
                const auto& candidates = top_level_ref_info.at(*ref_source);
                // If the ancestor produced several records, prefer one on our own contig;
                // failing that take the smallest (chrom, pos).  Which one is "right" is
                // genuinely ambiguous, so pick deterministically rather than by chance.
                auto chosen = candidates.begin();
                for (auto it = candidates.begin(); it != candidates.end(); ++it) {
                    if (it->first.first == toks[0]) {
                        chosen = it;
                        break;
                    }
                }
                top_ref = {chosen->first.first, chosen->first.second, chosen->second};
            }
            if (have_ref) {
                nesting_tags += ";RC=" + top_ref.seq;
                nesting_tags += ";RS=" + std::to_string(top_ref.start);
                nesting_tags += ";RD=" + std::to_string(top_ref.end + 1);
            }

            // rewrite the output string using the updated info toks
            output_variant_string.clear();
            for (size_t i = 0; i < toks.size(); ++i) {
                output_variant_string += toks[i];
                if (i == 7) {
                    output_variant_string += nesting_tags;
                }
                if (i != toks.size() - 1) {
                    output_variant_string += "\t";
                }
            }
            output_variant_record.second.clear();
            ret = zstdutil::CompressString(output_variant_string, output_variant_record.second);
            assert(ret == 0);
        }
    }
}
}


#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "vcf_output_caller.hpp"
#include "graph_caller.hpp"
#include "symbolic_allele.hpp"
#include "read_likelihood_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"
#include "utility.hpp"

//#define debug

namespace vg {

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
    if (emit_phasing) {
        // FORMAT/PS is the VCF phase set, which phasing tools read. It is unrelated to INFO/PS
        // above, vg's parent-snarl field; the two are in different namespaces, so both are legal,
        // and their descriptions say which is which.
        ss << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set: the phase of a "
           << "genotype is comparable only with others carrying the same PS. One phase set per "
           << "chain, so blocks are chromosome-scale -- much longer than a read-based phaser "
           << "gives, because the phase comes from the haplotype panel rather than from reads "
           << "spanning consecutive sites. Not the INFO/PS emitted under -A, which is a parent "
           << "snarl pointer\">" << endl;
    }
    ss << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << endl;
    if (block_records.is_enabled()) {
        ss << "##INFO=<ID=SB,Number=2,Type=Integer,Description=\"Index and count of this "
           << "difference block within its snarl. A snarl is written as one record per difference "
           << "block where the reference and the called haplotypes differ from each other in more "
           << "than one place inside it, or where its own record would repeat a child snarl's, so "
           << "the count can be 1. A block record's ID is the snarl's ID with _ and the index "
           << "appended. DOUBLE COUNTING: the per-sample evidence is the SNARL's, repeated on every "
           << "block, not apportioned between them -- AD, GL, GQ, GQI, GP and QUAL are identical "
           << "across the set, because the genotype likelihood was computed over whole-snarl "
           << "traversals and has no per-block decomposition. DP, DR and BL are per-site read "
           << "counts and are site-level by definition. So any consumer that sums, averages or "
           << "otherwise aggregates evidence across records must group by the snarl's ID first and "
           << "count each snarl once. Records without SB are unaffected: they are the only record "
           << "their snarl emitted.\">" << endl;
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

void VCFOutputCaller::set_linkage(LinkageCollector* collector, const gbwt::GBWT* gbwt,
                                  const vector<size_t>* sequence_to_haplotype) {
    this->linkage_collector = collector;
    this->panel_lookup = PanelLookup(gbwt, sequence_to_haplotype,
                                     collector != nullptr ? collector->panel_size() : 0);
}

bool VCFOutputCaller::buffered_record_key_less(const BufferedRecordKey& a, const BufferedRecordKey& b) {
    if (a.contig != b.contig) {
        return a.contig < b.contig;
    }
    if (a.position != b.position) {
        return a.position < b.position;
    }
    if (a.id != b.id) {
        return a.id < b.id;
    }
    return a.block < b.block;
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
    output_variants[omp_get_thread_num()].push_back(
        make_pair(BufferedRecordKey{var.sequenceName, (size_t)var.position, var.id, block}, dest));
    return true;
}

void VCFOutputCaller::resolve_linkage() {
    if (linkage_resolved) {
        return;
    }
    if (linkage_collector == nullptr) {
        resolve_linkage_level(0, true);
        return;
    }
    // Resolve every level, since chain construction skips entries of later levels
    // than the one being resolved. `max_level()` is read again on each pass, since a pass can
    // add a chain at a deeper level.
    for (size_t gen = 0;; ++gen) {
        const size_t deepest = linkage_collector->max_level();
        resolve_linkage_level(gen, gen >= deepest);
        if (gen >= deepest) {
            break;
        }
    }
}

/// The ID of the site a record belongs to: a block record's ID without the "_<index>" that
/// tells the site's block records apart, and any other record's ID unchanged.
static string block_site_name(const string& id) {
    size_t underscore = id.rfind('_');
    return underscore == string::npos ? id : id.substr(0, underscore);
}

size_t VCFOutputCaller::record_key_of(const Snarl& snarl) const {
    return vg::record_key_of(print_snarl(snarl, false));
}

// Each read's strand log-odds for the render, used by the anchors. Built here rather than taken
// from re-genotyping, which may not have run and whose table is built before `phase_sites` is
// final.
size_t VCFOutputCaller::phase_set_id(const string& contig, size_t phase_set) {
    return phase_set_ids.emplace(make_pair(contig, phase_set), phase_set_ids.size()).first->second;
}

void VCFOutputCaller::build_render_lambda() {
    render_lambda.clear();
    render_lambda_site.clear();
    render_lambda_phase_set.clear();
    render_lambda_temper = 0.0;
    render_lambda_ceiling = 1.0;
    if (phase_sites.empty()) {
        return;
    }
    RegenotypeCounters scratch;
    accumulate_lambda(phase_sites, phase_flips, render_lambda, scratch);
    for (const PhaseSite& site : phase_sites) {
        render_lambda_site[site.record_key] = &site;
    }
    // The last PhaseCall written winning, as in `build_render_phases`.
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        render_lambda_phase_set[pc.record_key] = phase_set_id(pc.contig, pc.phase_set);
    }
    // The summed strand log-odds overstate how sure the strand is, so they are tempered. Use the
    // temper re-genotyping fitted, where it ran; otherwise fit one here.
    if (regenotype_counters.fitted_temper > 0.0) {
        render_lambda_temper = regenotype_counters.fitted_temper;
        render_lambda_ceiling = regenotype_counters.fitted_ceiling;
    } else {
        double temper = -1.0;
        double ceiling = regenotype_params.ceiling < 0.0 ? 1.0 : regenotype_params.ceiling;
        RegenotypeCounters fit_scratch;
        fit_calibration(phase_sites, phase_flips, render_lambda, regenotype_params, temper, ceiling,
                        fit_scratch);
        if (fit_scratch.fitted_temper > 0.0) {
            render_lambda_temper = fit_scratch.fitted_temper;
            render_lambda_ceiling = fit_scratch.fitted_ceiling;
        }
    }
}

double VCFOutputCaller::read_strand_log_odds(size_t record_key, std::string_view read_name) const {
    if (render_lambda.empty() || render_lambda_temper <= 0.0) {
        return 0.0;
    }
    const uint64_t key = (uint64_t)std::hash<std::string_view>{}(read_name);
    const auto found = render_lambda.find(key);
    const auto ps = render_lambda_phase_set.find(record_key);
    const size_t phase_set = ps != render_lambda_phase_set.end() ? ps->second : NO_PHASE_SET;
    if (found == render_lambda.end()) {
        // The read reached no phased site.
        return 0.0;
    }
    if (!read_strand_usable(found->second, phase_set)) {
        // The read has a strand, but in another phase set, whose strands do not correspond to
        // this site's. NaN rather than 0, because a split homozygous site drops such a read but
        // places one with no strand by a coin (see `build_site_anchors`).
        return std::numeric_limits<double>::quiet_NaN();
    }
    double value = found->second.lambda;
    size_t sites = found->second.sites;
    // Subtract this record's own contribution, so that a site is not judged by its own evidence;
    // if it was the only one, there is nothing left.
    const auto site = render_lambda_site.find(record_key);
    if (site != render_lambda_site.end()) {
        unordered_map<uint64_t, double> own;
        site_own_log_odds(*site->second, phase_flips.count(record_key) != 0, own);
        const auto mine = own.find(key);
        if (mine != own.end()) {
            value -= mine->second;
            if (sites > 0) {
                --sites;
            }
        }
    }
    if (sites == 0) {
        return 0.0;
    }
    return calibrated_log_odds(value, render_lambda_temper, render_lambda_ceiling);
}

void VCFOutputCaller::build_render_phases() {
    // Built from the phasing the linkage pass accumulated, as read phasing left it. Sites with no line
    // are included; they are simply never looked up.
    render_phases.clear();
    if (!emit_phasing) {
        return;
    }
    render_phases.reserve(linkage_phased.size() * 2);
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        // Where a site has more than one PhaseCall, the last one written wins.
        render_phases[pc.record_key] = pc;
    }
}

void VCFOutputCaller::finalise_linkage_outputs() {
    // Built after every record has been rendered, since the mosaic needs to know which sites have
    // a line, which is not known while genotypes are being resolved.
    if (linkage_collector == nullptr) {
        return;
    }
    // Read from the collector, since each PhaseCall's `emitted` was copied before any line was
    // written.
    const std::unordered_set<size_t> emitted_records = linkage_collector->emitted_records();
    size_t unexplained = 0;
    size_t order_arbitrary = 0;
    // Count the phased sites, separating those that became records from those that did not.
    size_t phased_unwritten = 0;
    for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
        if (emitted_records.count(pc.record_key) == 0) {
            // Phased, since its children take their strand from it, but not a record, so it is kept
            // out of the mosaic and the record counts.
            ++phased_unwritten;
            continue;
        }
        // Count only the strands a site has. A haploid site has one strand and a wildcard, and the
        // wildcard can be in either slot: a haploid contig fills the first slot, while a nested
        // site on its parent's second strand fills the second.
        unexplained += (pc.ploidy == 1)
                       ? (pc.hap_first == LinkageModel::WILDCARD
                          && pc.hap_second == LinkageModel::WILDCARD)
                       : (pc.hap_first == LinkageModel::WILDCARD
                          || pc.hap_second == LinkageModel::WILDCARD);
        order_arbitrary += pc.order_arbitrary;
    }
    cerr << "[vg call] linkage: " << linkage_collector->num_sites() << " sites, "
         << (linkage_collector->bytes() / (1024.0 * 1024.0)) << " MB retained, "
         << linkage_changed << " genotypes moved by linkage, " << linkage_seconds << " s" << endl;
    if (linkage_collector->num_duplicate_live_keys() > 0) {
        // Duplicate keys need not change the output, but `retract` cannot handle those sites, since
        // it retracts only the first live entry.
        cerr << "[vg call] linkage: " << linkage_collector->num_duplicate_live_keys()
             << " sites recorded onto a key that already had a live entry; the retract path cannot"
             << " address these" << endl;
    }
    if (linkage_collector->model_params().hp_prior > 0.0) {
        cerr << "[vg call] linkage: " << linkage_collector->num_site_prior_entries()
             << " live entries decoded at a run-length site's own frequency exponent (--hp-prior)"
             << endl;
    }
    if (emit_phasing) {
        // At sites where a strand is on the wildcard, no panel haplotype names it, so the phase
        // across them rests on the transitions alone.
        cerr << "[vg call] phasing: " << (linkage_phased.size() - phased_unwritten)
             << " sites phased, " << unexplained
             << " with a strand the panel does not explain" << endl;
        if (phased_unwritten > 0) {
            // Sites that wrote no VCF line but are phased. A parent whose alleles differ only inside
            // its children is written as the reference and has no line, and its children still need
            // to know which of its strands carries the chain.
            cerr << "[vg call] phasing: " << phased_unwritten
                 << " collapsed sites phased with no line of their own, so their children can"
                 << " inherit a strand" << endl;
        }
        if (order_arbitrary > 0) {
            // Heterozygous sites where no panel haplotype on either strand carries either called
            // allele. The record is still phased and in the phase set, but its order came from
            // sorting the pair, so it is arbitrary.
            cerr << "[vg call] phasing: " << order_arbitrary
                 << " heterozygous sites carry an allele order the panel does not determine"
                 << endl;
        }
    }
    if (mosaic_writer.is_enabled()) {
        // Records only: the mosaic's segments are runs over sites of the call set, and it accounts
        // for exactly the written records.
        vector<LinkageCollector::PhaseCall> written;
        written.reserve(linkage_phased.size());
        for (const LinkageCollector::PhaseCall& pc : linkage_phased) {
            if (emitted_records.count(pc.record_key) != 0) {
                written.push_back(pc);
            }
        }
        mosaic_writer.write(written, panel_lookup, sample_name);
    }
}

void VCFOutputCaller::resolve_linkage_level(size_t level, bool last) {
    linkage_resolved = true;
    if (linkage_collector == nullptr) {
        return;
    }
    // Time the pass and report the collector's size.
    auto start = std::chrono::steady_clock::now();
    // `linkage_phased` accumulates across levels, since the model needs the earlier ones: a
    // nested site's strand is read from its parent's PhaseCall, and a clamped site's phase is
    // pinned to its chosen pair.
    const size_t moved =
        linkage_collector->resolve_level(level, last,
                                              emit_phasing ? &linkage_phased : nullptr);
    double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    linkage_seconds += seconds;
    // How many sites the model moved off the genotype the reads alone chose.
    linkage_changed += moved;
    if (!last) {
        // One line per level except the last: its site count, how many of its genotypes the
        // linkage model moved, and the seconds it took.
        cerr << "[vg call] linkage level " << level << ": "
             << linkage_collector->num_sites_at(level) << " sites, "
             << moved << " genotypes moved by linkage, " << seconds << " s" << endl;
        return;
    }

}

void VCFOutputCaller::write_variants(ostream& out_stream, const SnarlManager* snarl_manager) {
    assert(include_nested == false || snarl_manager != nullptr);
    if (include_nested) {
        update_nesting_info_tags(SnarlManagerSiteTree(*snarl_manager));
    }
    vector<pair<BufferedRecordKey, string>> all_variants;
    // Reserve once: doing it inside the loop below reallocates per thread buffer.
    size_t total_variants = 0;
    for (const auto& buf : output_variants) {
        total_variants += buf.size();
    }
    all_variants.reserve(total_variants);
    // `buf` must not be const, since std::move() over const iterators copies, and the whole VCF is
    // in memory here. Each buffer is freed as it is moved. This makes write_variants() usable only
    // once.
    for (auto& buf : output_variants) {
        std::move(buf.begin(), buf.end(), std::back_inserter(all_variants));
        buf.clear();
        buf.shrink_to_fit();
    }
    std::sort(all_variants.begin(), all_variants.end(),
              [](const pair<BufferedRecordKey, string>& v1,
                 const pair<BufferedRecordKey, string>& v2) {
                  return buffered_record_key_less(v1.first, v2.first);
              });
    // Resolve the linkage model, if it has not been resolved, before the records are written.
    resolve_linkage();
    finalise_linkage_outputs();


    // Each record is decompressed and given its linkage qualities on its own, so the records are
    // finished on several threads, a batch at a time, and each batch is written in order. Only
    // one batch of text is held at once.
    const size_t batch_records = 1 << 16;
    vector<string> lines;
    for (size_t batch_start = 0; batch_start < all_variants.size(); batch_start += batch_records) {
        const size_t batch_end = min(all_variants.size(), batch_start + batch_records);
        lines.assign(batch_end - batch_start, string());
#pragma omp parallel for schedule(dynamic, 256)
        for (size_t record_i = batch_start; record_i < batch_end; ++record_i) {
            const auto& v = all_variants[record_i];
            string& dest = lines[record_i - batch_start];
            int ret = zstdutil::DecompressString(v.second, dest);
            assert(ret == 0);
            // The record key is the hash of the site's ID, as `record_key_of` computes it, so the line
            // itself gives the site's identity; a block record's ID carries it before its suffix.
            // Computed once, when first needed; several records can share a (contig, position), and
            // each must get its own site's values.
            size_t line_key = 0;
            bool have_line_key = false;
            auto id_key = [&]() -> size_t {
                if (!have_line_key) {
                    size_t a = dest.find('\t');
                    size_t b = a == string::npos ? string::npos : dest.find('\t', a + 1);
                    size_t c = b == string::npos ? string::npos : dest.find('\t', b + 1);
                    if (c != string::npos) {
                        line_key = std::hash<string>{}(block_site_name(dest.substr(b + 1, c - b - 1)));
                    }
                    have_line_key = true;
                }
                return line_key;
            };
            if (linkage_collector != nullptr) {
                // Quality first, then phasing. The line already carries the chosen genotype, since it
                // was built from it.
                const auto& quality = linkage_collector->moved_quality();
                if (!quality.empty()) {
                    auto found = quality.find(id_key());
                    if (found != quality.end()) {
                        if (!ReadLikelihoodSnarlCaller::rewrite_quality_for_chosen_genotype(
                                dest, found->second, linkage_min_confidence)) {
                            ++quality_declined;
                        }
                    }
                }
            }
        }
        for (const string& line : lines) {
            // Not endl: flushing after every record made one write per record, millions on a whole
            // genome. The stream is flushed before vg exits.
            out_stream << line << '\n';
        }
    }
    if (phase_declined.load() > 0 || quality_declined.load() > 0) {
        cerr << "[vg call] linkage: " << phase_declined.load()
             << " phases refused by the record they were rendered onto, and "
             << quality_declined.load() << " quality rewrites refused" << endl;
    }
    // Reported after the records are rendered, since block emission happens as they are.
    block_records.report();
}


vector<int> VCFOutputCaller::phase_ordered_genotype(size_t record_key,
                                                    const vector<int>& genotype) const {
    vector<int> ordered = genotype;
    if (!emit_phasing || ordered.size() != 2) {
        return ordered;
    }
    const auto found = render_phases.find(record_key);
    // Only on an exact reversal. A PhaseCall that is not a permutation of the chosen pair is left
    // alone, as `emit_variant` refuses to apply one. For a homozygote the swap changes nothing.
    if (found != render_phases.end() && found->second.ploidy == 2
        && found->second.trav_first == ordered[1]
        && found->second.trav_second == ordered[0]) {
        std::swap(ordered[0], ordered[1]);
    }
    return ordered;
}

int VCFOutputCaller::phase_haploid_slot(size_t record_key, const vector<int>& genotype) const {
    if (!emit_phasing || genotype.size() != 1) {
        return 0;
    }
    const auto found = render_phases.find(record_key);
    if (found == render_phases.end() || found->second.ploidy != 1
        || found->second.nested_strand < 0) {
        // No nested strand means a haploid locus, such as chrY or a haploid --ploidy-bed region,
        // where slot 1 means nothing.
        return 0;
    }
    // Only where the phase names the allele chosen for this site, as `phase_ordered_genotype` and
    // `emit_variant` require.
    if (found->second.trav_first != genotype[0]) {
        return 0;
    }
    return (int)found->second.nested_strand;
}

void VCFOutputCaller::collect_anchors_for(const Snarl& snarl, const vector<int>& genotype,
                                          int haploid_slot,
                                          const unique_ptr<SnarlCaller::CallInfo>& call_info,
                                          bool is_leaf, double gqn, size_t record_key) {
    if (anchor_path.empty() || anchor_writer == nullptr || call_info == nullptr) {
        return;
    }
    if (anchor_params.leaf_only && !is_leaf) {
        return;
    }
    const auto* info =
        dynamic_cast<const ReadLikelihoodSnarlCaller::ReadLikelihoodCallInfo*>(call_info.get());
    if (info == nullptr || info->anchor_evidence == nullptr) {
        // A genotype derived from a parent rather than scored here, or a run whose caller is not the
        // read-likelihood one. There are no per-read responsibilities to partition on.
        return;
    }
    vector<AnchorWriter::Anchor> anchors;
    // Each read's strand log-odds, leaving out this record.
    vector<double> read_strand;
    // Built only where `build_site_anchors` reads it: at a diploid homozygote that may be split, or
    // at a heterozygous site under --anchors-phase-hets or --anchors-strict-hets. The test must
    // match its gate, which also checks the vector's length against `evidence.reads`.
    const bool splittable_hom = genotype.size() == 2 && genotype[0] == genotype[1];
    const bool tiltable_het = (anchor_params.phase_hets || anchor_params.strict_hets)
                              && genotype.size() == 2
                              && genotype[0] != genotype[1];
    if ((anchor_params.hom_split && splittable_hom) || tiltable_het) {
        read_strand.reserve(info->anchor_evidence->reads.size());
        for (const AnchorRead& read : info->anchor_evidence->reads) {
            read_strand.push_back(read_strand_log_odds(record_key, read_names().name(read.read)));
        }
    }
    build_site_anchors(*info->anchor_evidence, genotype, print_snarl(snarl),
                       gqn,
                       info->explained_share, haploid_slot, anchor_params, *anchor_params.counters,
                       anchors,
                       (anchor_params.hom_split || anchor_params.phase_hets
                        || anchor_params.strict_hets) ? &read_strand
                                                                             : nullptr);
    // A check for --anchors-hom-split, reported per run: at heterozygous sites, whose alleles show
    // which strand each read is on, how often the read's strand log-odds agree. The log-odds leave
    // the site out. Computed only when splitting is on, and not when the heterozygous placement
    // itself uses the strand log-odds, since the check would then compare the strand with
    // itself.
    if (anchor_params.hom_split && !anchor_params.phase_hets && !anchor_params.strict_hets
        && anchors.size() >= 2) {
        int slot_of_allele[2] = {-1, -1};
        int allele_of_slot[2] = {-1, -1};
        for (const AnchorWriter::Anchor& anchor : anchors) {
            if (anchor.slot >= 0 && anchor.slot < 2) {
                allele_of_slot[anchor.slot] = anchor.allele;
            }
        }
        if (allele_of_slot[0] >= 0 && allele_of_slot[1] >= 0
            && allele_of_slot[0] != allele_of_slot[1]) {
            (void)slot_of_allele;
            unordered_set<uint32_t> counted;
            for (const AnchorWriter::Anchor& anchor : anchors) {
                if (anchor.slot < 0 || anchor.slot > 1) {
                    continue;
                }
                for (const AnchorWriter::ReadRow& row : anchor.reads) {
                    if (!counted.insert(row.read).second) {
                        continue;   // both pins carry the same partition; count each read once
                    }
                    const double lo = read_strand_log_odds(record_key, read_names().name(row.read));
                    if (std::isnan(lo) || lo == 0.0) {
                        anchor_params.counters->phase_no_opinion.fetch_add(1);
                        continue;
                    }
                    const int phase_slot = lo > 0.0 ? 0 : 1;
                    const bool agree = phase_slot == anchor.slot;
                    anchor_params.counters->phase_checked.fetch_add(1);
                    if (agree) {
                        anchor_params.counters->phase_agree.fetch_add(1);
                    }
                    // The same threshold the split uses, --split-min-q, in natural-log units.
                    if (std::abs(lo) >= anchor_params.phase_min) {
                        anchor_params.counters->phase_confident.fetch_add(1);
                        if (agree) {
                            anchor_params.counters->phase_confident_agree.fetch_add(1);
                        }
                    }
                }
            }
        }
    }
    for (AnchorWriter::Anchor& anchor : anchors) {
        anchor_writer->add(std::move(anchor));
    }
}

void VCFOutputCaller::write_anchors() {
    if (anchor_path.empty() || anchor_writer == nullptr) {
        return;
    }
    size_t anchors = anchor_writer->anchor_count();
    size_t rows = anchor_writer->read_row_count();
    if (!anchor_writer->write(anchor_path, anchor_graph_name, sample_name, anchor_reads_source,
                              anchor_mismap_min, anchor_params)) {
        return;
    }
    cerr << "[vg call] anchors: " << anchors << " written over " << rows
         << " read placements to " << anchor_path << endl;
    if (anchor_params.counters != nullptr) {
        anchor_params.counters->report(cerr);
    }
}


static int countAlts(vcflib::Variant& var, int alleleIndex) {
    int alts = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
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
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
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

// this isn't from vcflib, but seems to make more sense than just returning the number of samples in
// the file again and again
static int countSamplesWithData(vcflib::Variant& var) {
    int samples_with_data = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
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
        for (const auto& output_variant_record : thread_buf) {
            contigs.insert(output_variant_record.first.contig);
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

void VCFOutputCaller::add_allele_path_to_info(const HandleGraph* graph, vcflib::Variant& v, int allele, const Traversal& trav,
                                              bool reversed, bool one_based) const {
    vector<NodeVisit> visits;
    visits.reserve(trav.size());
    for (const handle_t& handle : trav) {
        visits.emplace_back(graph->get_id(handle), graph->get_is_reverse(handle));
    }
    vg::add_allele_path_to_info(v, allele, visits, reversed, translation);
}

void VCFOutputCaller::add_allele_path_to_info(vcflib::Variant& v, int allele, const SnarlTraversal& trav,
                                              bool reversed, bool one_based) const {
    vg::add_allele_path_to_info(v, allele, visits_of(trav), reversed, translation);
}

string VCFOutputCaller::trav_string(const HandleGraph& graph, const SnarlTraversal& trav) const {
    string seq;
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        if (visit.node_id() > 0) {
            seq += graph.get_sequence(graph.get_handle(visit.node_id(), visit.backward()));
        } else {
            seq += print_snarl(visit.snarl(), true);
        }
    }
    return seq;    
}

thread_local VCFOutputCaller::NestedContext VCFOutputCaller::nested_context;
thread_local size_t VCFOutputCaller::current_level = 0;

bool VCFOutputCaller::is_symbolically_reference(const vector<SnarlTraversal>& called_traversals,
                                                int trav_idx, int ref_trav_idx,
                                                const Snarl& snarl) const {
    // Only when symbolic collapsing is on.
    if (symbolic_manager == nullptr || ref_trav_idx < 0 || trav_idx < 0 ||
        ref_trav_idx >= (int)called_traversals.size() ||
        trav_idx >= (int)called_traversals.size()) {
        return false;
    }
    return symbolically_equal(called_traversals[trav_idx], called_traversals[ref_trav_idx],
                              snarl, *symbolic_manager);
}

void VCFOutputCaller::set_symbolic_collapsing(const SnarlManager* manager) {
    this->symbolic_manager = manager;
    block_records.set_manager(manager);
    if (manager == nullptr) {
        record_steps.same_as_reference = nullptr;
        record_steps.count_site = nullptr;
        return;
    }
    record_steps.same_as_reference = [this](const Snarl& site, const vector<SnarlTraversal>& travs,
                                            int trav, int ref_trav_idx) {
        return is_symbolically_reference(travs, trav, ref_trav_idx, site);
    };
    record_steps.count_site = [this](const PathPositionHandleGraph& graph, const Snarl& site,
                                     const vector<SnarlTraversal>& travs,
                                     const vector<int>& genotype, int ref_trav_idx) {
        block_records.count_site(site, travs, ref_trav_idx);
    };
}

RecordOptions VCFOutputCaller::record_options() const {
    return RecordOptions{
        .sample_name = sample_name,
        .translation = translation,
        .max_uncalled_alleles = max_uncalled_alleles,
        .allele_merge_threshold = allele_merge_threshold,
        .allele_merge_min_len = allele_merge_min_len,
    };
}

int64_t VCFOutputCaller::phase_record_genotype(const Snarl& site, const vector<int>& site_genotype,
                                               const map<int, int>& trav_to_allele,
                                               string& gt) const {
    if (!emit_phasing || render_phases.empty()) {
        return -1;
    }
    auto found = render_phases.find(record_key_of(site));
    if (found == render_phases.end()) {
        return -1;
    }
    const LinkageCollector::PhaseCall& phase = found->second;
    // `find`, since `operator[]` would insert a default 0 on a miss, and the map's size is not
    // a bound on traversal indices.
    const auto found_a = trav_to_allele.find(phase.trav_first);
    const auto found_b = trav_to_allele.find(phase.trav_second);
    const int a = (phase.trav_first >= 0 && found_a != trav_to_allele.end())
                      ? found_a->second : -1;
    const int b = (phase.trav_second >= 0 && found_b != trav_to_allele.end())
                      ? found_b->second : -1;
    // The phased genotype must be a permutation of the one this record carries, so that
    // phasing cannot change a genotype.
    bool same = false;
    if (phase.ploidy == 1 && site_genotype.size() == 1) {
        same = (a >= 0 && a == site_genotype[0]);
    } else if (phase.ploidy == 2 && site_genotype.size() == 2) {
        same = (a >= 0 && b >= 0)
               && ((a == site_genotype[0] && b == site_genotype[1])
                   || (a == site_genotype[1] && b == site_genotype[0]));
    }

    if (!same) {
        ++phase_declined;
        return -1;
    }
    if (phase.ploidy == 1 && phase.nested_strand >= 0) {
        // A nested ploidy-1 site is one strand of a diploid locus, since the parent's
        // other allele deletes the chain. Written as a phased pair with "." on the other
        // strand, which is how the VCF records which strand carries the allele.
        gt = nested_strand_genotype(a, phase.nested_strand);
    } else if (phase.ploidy == 1) {
        // A haploid locus: one allele and no order; PS labels its phase set. "a|a"
        // would claim a homozygous diploid call.
        gt = std::to_string(a);
    } else {
        gt = std::to_string(a) + "|" + std::to_string(b);
    }
    return (int64_t)phase.phase_set;
}

bool VCFOutputCaller::emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                                   const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                                   const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                                   const string& ref_path_name, int ref_offset, bool genotype_snarls, int ploidy,
                                   function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_string) {
    
#ifdef debug
    cerr << "emitting variant for " << pb2json(snarl) << endl;
    for (int i = 0; i < called_traversals.size(); ++i) {
        if (i == ref_trav_idx) {
            cerr << "*";
        }
        cerr << "ct[" << i << "]=" << pb2json(called_traversals[i]) << endl;
    }
    for (int i = 0; i < genotype.size(); ++i) {
        cerr << "gt[" << i << "]=" << genotype[i] << endl;
    }
#endif

    if (trav_to_string == nullptr) {
        trav_to_string = [&](const vector<SnarlTraversal>& travs, const vector<int>& travs_genotype, int trav_allele, int genotype_allele, int ref_trav_idx) {
            return trav_string(graph, travs[trav_allele]);    
        };
    }

    if (record_steps.count_site) {
        record_steps.count_site(graph, snarl, called_traversals, genotype, ref_trav_idx);
    }

    SiteAlleles alleles;
    alleles.spell = [&](int trav, int genotype_index) {
        return trav_to_string(called_traversals, genotype, trav, genotype_index, ref_trav_idx);
    };
    alleles.visits = [&](int trav) {
        return visits_of(called_traversals[trav]);
    };

    SiteHooks hooks;
    if (record_steps.same_as_reference) {
        hooks.same_as_reference = [&](int trav) {
            return record_steps.same_as_reference(snarl, called_traversals, trav, ref_trav_idx);
        };
    }
    if (record_steps.phase) {
        hooks.phase = [&](const vector<int>& site_genotype, const map<int, int>& trav_to_allele,
                          string& gt) {
            return record_steps.phase(snarl, site_genotype, trav_to_allele, gt);
        };
    }
    hooks.fill_info = [&](const vector<int>& site_trav, const vector<int>& site_genotype,
                          vcflib::Variant& variant) {
        // The "*" placeholder is an empty traversal.
        vector<SnarlTraversal> site_traversals;
        site_traversals.reserve(site_trav.size());
        for (int trav : site_trav) {
            site_traversals.push_back(trav >= 0 ? called_traversals[trav] : SnarlTraversal());
        }
        snarl_caller.update_vcf_info(snarl, site_traversals, site_genotype, call_info, sample_name,
                                     variant);
    };
    hooks.gl_layout = record_steps.gl_layout ? record_steps.gl_layout(call_info.get())
                                             : GLLayout::IMajor;

    const SiteToWrite site{
        .start = graph.get_handle(snarl.start().node_id(), snarl.start().backward()),
        .end = graph.get_handle(snarl.end().node_id(), snarl.end().backward()),
        .ref_path_name = ref_path_name,
        .ref_offset = ref_offset,
        .genotype = genotype,
        .ref_trav_idx = ref_trav_idx,
        .traversal_count = called_traversals.size(),
        .ploidy = ploidy,
        .genotype_snarls = genotype_snarls,
    };
    SiteRecord record = build_site_record(graph, site, alleles, hooks, record_options());
    vcflib::Variant& out_variant = record.variant;

    // One record per difference block, where that changes the output. The site record above is
    // finished, so the blocks take every field they do not redefine from it. -1 means the site was
    // declined, and the site record below is written as it is.
    if (record_steps.write_blocks) {
        const int block_lines = record_steps.write_blocks(graph, snarl, called_traversals, genotype,
                                                          ref_trav_idx, record, hooks.gl_layout,
                                                          genotype_snarls);
        if (block_lines >= 0) {
            // There is no single line for this snarl, and each block numbers its own alleles.
            if (record_steps.site_filed) {
                record_steps.site_filed(snarl, map<int, int>(), 0, block_lines > 0);
            }
            return block_lines > 0;
        }
    }

    // Whether this site wants a line. A pair of traversals differing from the reference only inside
    // child chains is written as allele 0 and leaves `alt` empty; such a site has no line, but its
    // children need it.
    const bool wants_line = genotype_snarls || !out_variant.alt.empty();
    bool added = false;
    if (wants_line) {
        added = add_variant(out_variant);
    } else if (include_nested) {
        // A site with nothing to report still knows where its children sit, so its reference
        // interval is kept for their RC, RS and RD, as Deconstructor::deconstruct_site does.
        suppressed_ref_info[omp_get_thread_num()][out_variant.id] =
            {out_variant.sequenceName, static_cast<size_t>(out_variant.position),
             out_variant.ref.length()};
    }
    if (record_steps.site_filed) {
        record_steps.site_filed(snarl, record.trav_to_allele, called_traversals.size(), added);
    }
    if (wants_line && !added) {
        stringstream ss;
        ss << out_variant;
        cerr << "Warning [vg call]: Skipping variant at " << out_variant.sequenceName << ":" << out_variant.position
             << " with ID=" << out_variant.id << " because its line length of " << ss.str().length() << " exceeds vg's limit of "
             << VCFOutputCaller::max_vcf_line_length << endl;
    }
    // True when the record has nothing to write, false only when add_variant refused a line; the
    // linkage pass depends on the difference.
    return wants_line ? added : true;
}

tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> VCFOutputCaller::get_ref_interval(
    const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name) const {
    return vg::get_ref_interval(graph, graph.get_handle(snarl.start().node_id(), snarl.start().backward()),
                                graph.get_handle(snarl.end().node_id(), snarl.end().backward()),
                                ref_path_name);
}

pair<string, int64_t> VCFOutputCaller::get_ref_position(const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name,
                                                        int64_t ref_path_offset) const {
    return vg::get_ref_position(graph, graph.get_handle(snarl.start().node_id(), snarl.start().backward()),
                                graph.get_handle(snarl.end().node_id(), snarl.end().backward()),
                                ref_path_name, ref_path_offset);
}

void VCFOutputCaller::flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override) const {
    vg::flatten_common_allele_ends(variant, backward, len_override);
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

string VCFOutputCaller::print_snarl(const HandleGraph* graph, const handle_t& snarl_start,
                                    const handle_t& snarl_end, bool in_brackets) const {
    return print_snarl(graph->get_id(snarl_start), graph->get_is_reverse(snarl_start),
                       graph->get_id(snarl_end), graph->get_is_reverse(snarl_end), in_brackets);
}
string VCFOutputCaller::print_snarl(const Snarl& snarl, bool in_brackets) const {
    return print_snarl(snarl.start().node_id(), snarl.start().backward(), snarl.end().node_id(),
                       snarl.end().backward(), in_brackets);
}
string VCFOutputCaller::print_flipped_snarl(const Snarl& snarl, bool in_brackets) const {
    return print_snarl(snarl.end().node_id(), !snarl.end().backward(), snarl.start().node_id(),
                       !snarl.start().backward(), in_brackets);
}
string VCFOutputCaller::print_snarl(nid_t start_node_id, bool start_backward, nid_t end_node_id,
                                    bool end_backward, bool in_brackets) const {
    return site_name(start_node_id, start_backward, end_node_id, end_backward, translation,
                     in_brackets);
}

void VCFOutputCaller::scan_snarl(const string& allele_string, function<void(const string&, Snarl&)> callback) const {
    int left = -1;
    int last = 0;
    Snarl snarl;
    string frag;
    for (int i = 0; i < allele_string.length(); ++i) {
        if (allele_string[i] == '(') {
            assert(left == -1);
            if (last < i) {
                frag = allele_string.substr(last, i-last);
                callback(frag, snarl);
            }
            left = i;
        } else if (allele_string[i] == ')') {
            assert(left >= 0 && i > left + 3);
            frag = allele_string.substr(left + 1, i - left - 1);
            auto toks = split_delims(frag, "><");
            assert(toks.size() == 2);
            assert(frag[0] == '<' || frag[0] == '>');
            int64_t start = std::stoi(toks[0]);
            snarl.mutable_start()->set_node_id(start);
            snarl.mutable_start()->set_backward(frag[0] == '<');
            assert(frag[toks[0].size() + 1] == '<' || frag[toks[0].size() + 1] == '>');
            int64_t end = std::stoi(toks[1]);
            snarl.mutable_end()->set_node_id(abs(end));
            snarl.mutable_end()->set_backward(frag[toks[0].size() + 1] == '<');
            callback("", snarl);
            left = -1;
            last = i + 1;
        }
    }
    if (last == 0) {
        callback(allele_string, snarl);
    } else {
        frag = allele_string.substr(last);
        callback(frag, snarl);
    }
}

void VCFOutputCaller::update_nesting_info_tags(const SiteTree& sites) {

    // A site's name as print_snarl spells it, and the name it has when read the other way.
    auto name_of = [&](SiteTree::site_t site) {
        const SiteEnds e = sites.ends_of(site);
        return print_snarl(e.start_id, e.start_backward, e.end_id, e.end_backward, false);
    };
    auto flipped_name_of = [&](SiteTree::site_t site) {
        const SiteEnds e = sites.ends_of(site);
        return print_snarl(e.end_id, !e.end_backward, e.start_id, !e.start_backward, false);
    };

    // Merge the per-thread suppressed-site intervals collected during calling.  These are sites
    // that never reached the VCF, so pass 1 below cannot see them, but a record nested under one
    // has no other way to name a reference position.
    unordered_map<string, SuppressedRef> suppressed_ref;
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
    // What passes 1 and 2 read from each record: its site's name, CHROM, POS and REF length. The
    // records are decompressed once, in parallel, and the indexes are then filled in record
    // order, as reading the records one by one fills them.
    struct RecordFields {
        string name;
        string chrom;
        string pos;
        size_t ref_len = 0;
        bool top_level = false;
    };
    vector<vector<RecordFields>> record_fields(output_variants.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < output_variants.size(); ++b) {
        vector<RecordFields>& fields = record_fields[b];
        fields.reserve(output_variants[b].size());
        string output_variant_string;
        for (auto& output_variant_record : output_variants[b]) {
            output_variant_string.clear();
            int ret = zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 5);
            RecordFields f;
            f.name = block_site_name(toks[2]);
            f.ref_len = toks[3].length();
            f.chrom = std::move(toks[0]);
            f.pos = std::move(toks[1]);
            fields.push_back(std::move(f));
        }
    }
    for (const vector<RecordFields>& fields : record_fields) {
        for (const RecordFields& f : fields) {
            chrom_of_name.emplace(f.name, intern_chrom(f.chrom));
        }
    }

    // index the snarl tree by name
    //
    // Only the names of sites in the VCF are ever looked up, so only they are indexed. The
    // snarls are visited in the same order as for an index of every name, so a name that two
    // snarls print goes to the same one.
    unordered_map<string, SiteTree::site_t> name_to_snarl;
    name_to_snarl.reserve(chrom_of_name.size());
    if (translation == nullptr) {
        // A name is the snarl's two boundary visits, so each VCF name is read back into its visits
        // once, and every snarl of the graph is matched by its visits instead of by printing both
        // its names, which meant tens of millions of names and a string-table lookup for each. A
        // name is read back only if printing what was read gives the name again, so a name and a
        // pair of visits correspond one to one, and a snarl matches a name exactly when it prints
        // that name.
        struct Ends {
            nid_t start_id;
            nid_t end_id;
            bool start_backward;
            bool end_backward;
            bool operator==(const Ends& other) const {
                return start_id == other.start_id && end_id == other.end_id
                       && start_backward == other.start_backward
                       && end_backward == other.end_backward;
            }
        };
        struct EndsHash {
            size_t operator()(const Ends& e) const {
                size_t h = std::hash<nid_t>()(e.start_id);
                h ^= std::hash<nid_t>()(e.end_id) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                return h ^ ((size_t)e.start_backward << 1) ^ (size_t)e.end_backward;
            }
        };
        auto read_ends = [&](const string& name, Ends& ends) -> bool {
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
            return print_snarl(ends.start_id, ends.start_backward, ends.end_id, ends.end_backward,
                               false) == name;
        };
        unordered_map<Ends, const string*, EndsHash> name_of_ends;
        name_of_ends.reserve(chrom_of_name.size());
        for (const auto& kv : chrom_of_name) {
            Ends ends;
            if (read_ends(kv.first, ends)) {
                name_of_ends.emplace(ends, &kv.first);
            }
        }
        // The VCF names a snarl matches: its own, then its flipped one (as call sometimes messes
        // with orientation).
        auto for_each_match = [&](SiteTree::site_t snarl, const function<void(const string&)>& match) {
            const SiteEnds e = sites.ends_of(snarl);
            auto own = name_of_ends.find(Ends{e.start_id, e.end_id, e.start_backward, e.end_backward});
            if (own != name_of_ends.end()) {
                match(*own->second);
            }
            auto flipped = name_of_ends.find(Ends{e.end_id, e.start_id, !e.end_backward,
                                                  !e.start_backward});
            if (flipped != name_of_ends.end()) {
                match(*flipped->second);
            }
        };
        // The snarls are matched on several threads, each into a list of its own. Only a name
        // that two different snarls match could depend on the order the snarls are visited in;
        // if there is one, the matches are made again in preorder, as they always were.
        vector<vector<pair<const string*, SiteTree::site_t>>> found(max(1, omp_get_max_threads()));
        sites.for_each_site([&](SiteTree::site_t snarl) {
            auto& mine = found[omp_get_thread_num()];
            for_each_match(snarl, [&](const string& name) {
                mine.emplace_back(&name, snarl);
            });
        }, false);
        bool ambiguous = false;
        for (const auto& thread_found : found) {
            for (const auto& name_and_snarl : thread_found) {
                auto placed = name_to_snarl.emplace(*name_and_snarl.first, name_and_snarl.second);
                if (!placed.second && placed.first->second != name_and_snarl.second) {
                    ambiguous = true;
                }
            }
        }
        if (ambiguous) {
            name_to_snarl.clear();
            sites.for_each_site([&](SiteTree::site_t snarl) {
                for_each_match(snarl, [&](const string& name) {
                    name_to_snarl[name] = snarl;
                });
            }, true);
        }
    } else {
        // Translated names are not node IDs, so they are printed and compared.
        sites.for_each_site([&](SiteTree::site_t snarl) {
                string snarl_name = name_of(snarl);
                if (chrom_of_name.count(snarl_name) != 0) {
                    name_to_snarl[std::move(snarl_name)] = snarl;
                }
                // also add a map from the flipped snarl (as call sometimes messes with orientation)
                string flipped_name = flipped_name_of(snarl);
                if (chrom_of_name.count(flipped_name) != 0) {
                    name_to_snarl[std::move(flipped_name)] = snarl;
                }
            }, true);
    }

    // pass 2) identify top-level snarls (those with no ancestors in VCF)
    // and store reference info only for them
    struct RefInfo {
        string chrom;
        size_t pos;
        size_t ref_len;
    };
    // Keyed by snarl name, then by (chrom, pos), because a snarl ID can carry more than one
    // record: a cyclic reference emits two, both with the same ID (see
    // nesting/cyclic_ref_multiple_variants.gfa, which gives two <5<1 records at POS 20 and 44).
    // A plain name -> RefInfo map was last-write-wins, so both records were handed the
    // surviving one's interval and the record at POS 20 reported RS=44.  The inner map is
    // ordered so that picking begin() is deterministic regardless of thread scheduling.
    unordered_map<string, map<pair<string, size_t>, size_t>> top_level_ref_info;

    // Helper to check if a snarl is top-level (no ancestors in VCF)
    auto is_top_level = [&](const string& name) -> bool {
        auto it = name_to_snarl.find(name);
        if (it == name_to_snarl.end()) return true; // not found, treat as top-level
        SiteTree::site_t snarl = it->second;
        while ((snarl = sites.parent_of(snarl))) {
            string cur_name = name_of(snarl);
            string flipped_name = flipped_name_of(snarl);
            if (chrom_of_name.count(cur_name) || chrom_of_name.count(flipped_name)) {
                return false; // has ancestor in VCF
            }
        }
        return true; // no ancestors in VCF
    };

    // Second pass through variants to extract ref info only for top-level snarls. Whether a
    // record's site is top level depends only on the indexes above, so the records are tested in
    // parallel; the ref info is then stored in record order.
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < record_fields.size(); ++b) {
        for (RecordFields& f : record_fields[b]) {
            f.top_level = is_top_level(f.name);
        }
    }
    for (const vector<RecordFields>& fields : record_fields) {
        for (const RecordFields& f : fields) {
            if (f.top_level) {
                top_level_ref_info[f.name][make_pair(f.chrom, static_cast<size_t>(stoul(f.pos)))] =
                    f.ref_len;
            }
        }
    }
    vector<vector<RecordFields>>().swap(record_fields);

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
            if (chrom_of_name.count(cur_name)) {
                // only count snarls that are in the vcf
                hit = &cur_name;
            } else if (chrom_of_name.count(flipped_name)) {
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
            //string& output_variant_string = output_variant_record.second;
            vector<string> toks = split_delims(output_variant_string, "\t", 9);
            // Keyed by the site, so that a block record gets its site's tags.
            const string name = block_site_name(toks[2]);

            auto [contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
                  suppressed_name] = get_nesting_tags(name, toks[0]);
            // LV is the level within this record's own reference contig, so that a gRef fragment's
            // records start at level 0 on their own contig.
            //
            // CH counts the ancestors that have a record here, so it would be 0 for a record on a
            // gRef fragment whose enclosing site wrote no line. The contig's gRef level, which
            // equals the CH of every record on a fragment, is used as a floor. So CH >= 1 does not
            // imply a parent record in the VCF, or PS.
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
            RefInfo top_ref;
            if (ref_source == nullptr) {
                if (!GrefCover::is_gref_name(toks[0])) {
                    // Not on a gref fragment, so our own interval is already a position a reader
                    // can look up, and it is the narrower answer of the two.  Keep it rather than
                    // reach for an enclosing site that produced no record: doing that would
                    // repoint every such record on a reference contig at a site LV and CH say it
                    // has no ancestor in.  The records that need the reach are the fragments
                    // below, which have no usable coordinate of their own.
                    top_ref = {toks[0], static_cast<size_t>(stoul(toks[1])), toks[3].length()};
                } else {
                    auto sup_it = suppressed_name.empty() ? suppressed_ref.end()
                                                          : suppressed_ref.find(suppressed_name);
                    if (sup_it != suppressed_ref.end()) {
                        top_ref = {sup_it->second.chrom, sup_it->second.pos, sup_it->second.ref_len};
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
                nesting_tags += ";RC=" + top_ref.chrom;
                nesting_tags += ";RS=" + std::to_string(top_ref.pos);
                nesting_tags += ";RD=" + std::to_string(top_ref.pos + top_ref.ref_len);
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


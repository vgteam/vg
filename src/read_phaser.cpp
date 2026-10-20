#include <omp.h>

#include "read_phaser.hpp"
#include "read_likelihood_caller.hpp"

namespace vg::multipass {

void ReadPhaser::configure(bool enable, const ReadPhasingParams& settings) {
    on = enable;
    params = settings;
}

/**
 * Reduce one record's read evidence to the heterozygous pair its phase names.
 *
 * Slot 0 of the result holds the phase's first traversal, so slot 0 is strand 0, as for the first
 * field of GT and the slot column of the anchor file.
 *
 * @return Whether `site` was filled. It is not filled for a site with no two distinct strands to
 * order, no read evidence, or no read that fits either allele.
 */
static bool pair_evidence(const StagedSite& rec, const LinkageCollector::PhaseCall& pc,
                          PhaseSite& site) {
    if (pc.ploidy != 2 || pc.trav_first < 0 || pc.trav_second < 0
        || pc.trav_first == pc.trav_second) {
        return false;
    }
    const SiteScore* score = rec.score;
    if (score == nullptr) {
        return false;
    }
    PhaseReadEvidence converted;
    const PhaseReadEvidence* pe = score->read_phasing_evidence(converted);
    if (pe == nullptr) {
        return false;
    }
    const auto a0 = (size_t)pc.trav_first;
    const auto a1 = (size_t)pc.trav_second;
    if (a0 >= pe->n_alleles || a1 >= pe->n_alleles) {
        // The phase names a traversal the evidence has no column for, so any column we read
        // would belong to another allele.
        return false;
    }
    site = reduce_to_pair(*pe, a0, a1);
    if (site.read_key.empty()) {
        return false;
    }
    site.record_key = rec.record_key;
    site.position = pc.position;
    site.panel_ordered = !pc.order_arbitrary;
    return true;
}

void ReadPhaser::collect_sites(
    StagedSiteTable& staged, const PhaseTable& phases,
    const unordered_map<size_t, size_t>& phase_index, vector<PhaseSite>& sites,
    const function<size_t(const string& contig, size_t phase_set)>& phase_set_id) const {
    const vector<LinkageCollector::PhaseCall>& calls = phases.calls();
    // Build the sites on several threads, a block of records at a time into each block's own
    // list. Then gather them in record order, numbering the phase sets as they are first seen.
    const vector<StagedSite*> records = staged.in_order(true);
    const size_t block_records = 4096;
    const size_t n_blocks = (records.size() + block_records - 1) / block_records;
    vector<vector<PhaseSite>> block_sites(n_blocks);
    // For each site in a block, the index in `calls` of its phase.
    vector<vector<size_t>> block_calls(n_blocks);
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < n_blocks; ++b) {
        const size_t end = min(records.size(), (b + 1) * block_records);
        for (size_t r = b * block_records; r < end; ++r) {
            const StagedSite& rec = *records[r];
            const auto found = phase_index.find(rec.record_key);
            if (found == phase_index.end()) {
                continue;
            }
            PhaseSite site;
            if (pair_evidence(rec, calls[found->second], site)) {
                block_sites[b].push_back(std::move(site));
                block_calls[b].push_back(found->second);
            }
        }
    }
    size_t total_sites = 0;
    for (const vector<PhaseSite>& block : block_sites) {
        total_sites += block.size();
    }
    sites.reserve(total_sites);
    for (size_t b = 0; b < n_blocks; ++b) {
        for (size_t i = 0; i < block_sites[b].size(); ++i) {
            const LinkageCollector::PhaseCall& pc = calls[block_calls[b][i]];
            block_sites[b][i].phase_set = phase_set_id(pc.contig, pc.phase_set);
            sites.push_back(std::move(block_sites[b][i]));
        }
        vector<PhaseSite>().swap(block_sites[b]);
    }
}

/**
 * List the phased sites of `staged` that are not dropped, with their parents and levels.
 *
 * A chain that an enclosing record spells inline (`reported_inline`) is included. It still has
 * anchors read from its strand, and the strands of its children depend on its own. A dropped
 * chain is left out, since the sample carries neither it nor anything inside it.
 */
static vector<PhaseTable::NestedLink> nested_links(
    StagedSiteTable& staged, const unordered_map<size_t, size_t>& phase_index) {
    vector<PhaseTable::NestedLink> links;
    staged.for_each([&](const StagedSite& rec) {
        if (!rec.dropped && phase_index.contains(rec.record_key)) {
            links.emplace_back(rec.record_key, rec.placement.parent_record_key,
                               rec.placement.level);
        }
    });
    return links;
}

void ReadPhaser::phase(
    StagedSiteTable& staged, PhaseTable& phases, ReadStrandTable& strands,
    const function<size_t(const string& contig, size_t phase_set)>& phase_set_id) {
    if (!on || phases.calls().empty()) {
        return;
    }
    // Re-genotyping calls this again on new genotypes, and the report describes only the phase
    // the output carries.
    counters = ReadPhasingCounters();
    const std::unordered_map<size_t, size_t> phase_index = phases.index();

    // Keep the sites in `strands`, since re-genotyping uses them.
    vector<PhaseSite>& sites = strands.sites();
    sites.clear();
    collect_sites(staged, phases, phase_index, sites, phase_set_id);
    if (sites.empty()) {
        return;
    }

    strands.read_ordered().clear();
    strands.flips() = read_phase_flips(sites, params, counters, &strands.read_ordered());

    // Block records spell the phase in their ALTs, so swapping a nested site can change the
    // allele numbers in its GT.
    counters.strands_rederived += phases.swap_strands(strands.flips(),
                                                      nested_links(staged, phase_index));
    report();
}

void ReadPhaser::report() const {
    const ReadPhasingCounters& c = counters;
    cerr << "[vg call] read phasing: " << c.sites << " het sites, " << c.reliable
         << " reliable, " << c.chains << " blocks, " << c.breaks << " chain breaks ("
         << c.breaks_no_reads << " with no spanning read), " << c.hung
         << " sites hung off the chain (" << c.hung_no_reads << " with no read), " << c.flipped
         << " re-phased against the panel, " << c.strands_rederived
         << " nested strands carried with their parent" << endl;
    if (c.demoted_incoherent > 0) {
        cerr << "[vg call] read phasing: " << c.demoted_incoherent
             << " sites demoted for low phase coherence over " << c.coherence_rounds_run
             << " rounds";
        if (c.coherence_unconverged) {
            cerr << " -- " << c.coherence_unconverged
                 << " chains were STILL demoting at the round cap, so their chain is not"
                    " a coherent fixed point";
        } else {
            cerr << " (every chain reached a coherent fixed point)";
        }
        cerr << endl;
    }
}

}

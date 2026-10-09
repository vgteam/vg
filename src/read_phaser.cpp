#include <omp.h>

#include "read_phaser.hpp"
#include "read_likelihood_caller.hpp"

namespace vg {

void ReadPhaser::configure(bool on, const ReadPhasingParams& params) {
    this->on = on;
    this->params = params;
}

void ReadPhaser::phase(
    StagedSiteTable& staged, PhaseTable& phases, ReadStrandTable& strands,
    const function<size_t(const string& contig, size_t phase_set)>& phase_set_id) {
    if (!on || phases.calls().empty()) {
        return;
    }
    // Reset, since re-genotyping calls this again on the new genotypes, and the report should
    // describe the phase the output carries.
    counters = ReadPhasingCounters();
    // Index the phasing by record key, the last one written winning.
    const std::unordered_map<size_t, size_t> phase_index = phases.index();
    const vector<LinkageCollector::PhaseCall>& calls = phases.calls();

    // Kept in `strands`, since re-genotyping uses these sites.
    vector<PhaseSite>& sites = strands.sites();
    sites.clear();
    // Each record's site is built from its own evidence alone, so the sites are built on several
    // threads, a block of records at a time into the block's own list, and then gathered in record
    // order. Phase sets are numbered in the order they are first seen, so a site's number is given
    // in that gathering pass, in record order.
    const vector<StagedSite*> records = staged.in_order(true);
    const size_t block_records = 4096;
    const size_t n_blocks = (records.size() + block_records - 1) / block_records;
    vector<vector<PhaseSite>> block_sites(n_blocks);
    // For each site in a block, its PhaseCall's index in `calls`.
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
            const LinkageCollector::PhaseCall& pc = calls[found->second];
            if (pc.ploidy != 2 || pc.trav_first < 0 || pc.trav_second < 0
                || pc.trav_first == pc.trav_second) {
                // Homozygous, haploid, or unplaced: no two strands to order.
                continue;
            }
            const SiteScore* info = rec.score;
            if (info == nullptr) {
                continue;
            }
            PhaseReadEvidence converted;
            const PhaseReadEvidence* pe = info->read_phasing_evidence(converted);
            if (pe == nullptr) {
                continue;
            }
            const size_t a0 = (size_t)pc.trav_first, a1 = (size_t)pc.trav_second;
            if (a0 >= pe->n_alleles || a1 >= pe->n_alleles) {
                // The PhaseCall names a traversal this site's matrix does not have. Skipped, since
                // reading another column would take the phase from the wrong allele.
                continue;
            }
            // Slot order is the PhaseCall's order, so slot 0 is strand 0, as for GT's first field and
            // the anchor file's slot column.
            PhaseSite site = reduce_to_pair(*pe, a0, a1);
            if (site.read_key.empty()) {
                continue;
            }
            site.record_key = rec.record_key;
            site.position = pc.position;
            block_sites[b].push_back(std::move(site));
            block_calls[b].push_back(found->second);
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
    if (sites.empty()) {
        return;
    }

    strands.flips() = read_phase_flips(sites, params, counters);

    // Apply by swapping the chosen pair's order, and carry the swaps down the nesting tree. Nested
    // sites are reordered too. Under -A, block records spell the phase in their ALTs, so
    // reordering a nested site can change its GT's allele numbers. Every recorded chain is linked,
    // including one whose line an enclosing block's ALT spells (`reported_inline`): it still has
    // anchors, read from its strand, and its children's strands depend on its own. A dropped
    // chain is left out, since the sample does not carry it or anything inside it.
    vector<PhaseTable::NestedLink> links;
    staged.for_each([&](const StagedSite& rec) {
        if (!rec.dropped && phase_index.count(rec.record_key) != 0) {
            links.push_back({rec.record_key, rec.parent_record_key, rec.level});
        }
    });
    counters.strands_rederived += phases.swap_strands(strands.flips(), std::move(links));

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
             << " rounds"
             << (c.coherence_unconverged
                     ? " -- " + std::to_string(c.coherence_unconverged)
                           + " chains were STILL demoting at the round cap, so their chain is not"
                             " a coherent fixed point"
                     : " (every chain reached a coherent fixed point)")
             << endl;
    }
}

}

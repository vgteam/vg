#include <algorithm>
#include <limits>

#include "block_records.hpp"
#include "utility.hpp"

//#define debug

namespace vg {

// The names of the AtomizeRefusal reasons, in order. The initializer sets the
// size, so that the check below fails when a name is missing as well as when one is extra.
static const char* const g_atomize_refuse_name[] = {
    "the genotyper returned no genotype: ploidy 0, or no read the matrix could place",
    "no reference traversal",
    "the snarl does not resolve",
    "the reference projection is empty",
    "the alignment's reference-to-alt step map has the wrong length",
    "no difference blocks: every called haplotype takes the reference route here",
    "an indel needs an anchor base and the snarl has none to its left",
    "the visit left of the anchor has no sequence",
    "every block spelled the reference's own bases: a route difference with no sequence difference",
    "one block, saying what the site record already says",
    "-L merged the called alleles, which blocks would spell apart again",
    "one block that two strands' routes spell differently within one site allele, so the site cannot give it a genotype",
    "one block, where a chain crossed more than once may leave its record short of what the site record says",
};
static_assert(sizeof(g_atomize_refuse_name) / sizeof(g_atomize_refuse_name[0])
                  == (size_t)AtomizeRefusal::Count,
              "each AtomizeRefusal must have a name in g_atomize_refuse_name, "
              "and each name must have a reason");

void BlockRecordWriter::report() const {
    size_t unresolvable = counters.site_unresolvable.load();
    // Keyed on whether block emission ran at all, not on any refusal counter, so that the line
    // below is written whenever block emission ran.
    if (counters.sites.load() == 0) {
        return;
    }


    // `unresolvable` should be zero on the ordinary path, where every site is a managed snarl and a
    // reversed one resolves through its reversed boundaries. Under -I/--chains it need not be: a
    // chain piece is a constructed snarl the manager does not know. The second number counts sites
    // that resolved only through their reversed boundaries.
    cerr << "[vg call] atomize: " << unresolvable
         << " sites where projection is inert because the snarl does not resolve, "
         << counters.site_reversed.load()
         << " resolved as the reversal flip_snarl produces" << endl;


    if (counters.child_inlined.load() > 0) {
        cerr << "[vg call] atomize: " << counters.child_inlined.load()
             << " child chains left without a line because a block ALT already spells them" << endl;
    }
    {
        // One line, listing only the reasons that occurred.
        const size_t reasons = (size_t)AtomizeRefusal::Count;
        size_t total = 0;
        for (size_t i = 0; i < reasons; ++i) {
            total += counters.refuse[i].load();
        }
        if (total > 0) {
            cerr << "[vg call] atomize: " << total << " sites declined block emission, so the site"
                 << " record stands:";
            bool first = true;
            for (size_t i = 0; i < reasons; ++i) {
                size_t n = counters.refuse[i].load();
                if (n > 0) {
                    cerr << (first ? " " : "; ") << n << " " << g_atomize_refuse_name[i];
                    first = false;
                }
            }
            cerr << endl;
        }
    }
    if (counters.split_sites.load() > 0) {
        cerr << "[vg call] atomize: " << counters.split_sites.load()
             << " sites written as their difference blocks rather than their site record, "
             << counters.split_lines.load()
             << " lines" << endl;
    }
}

BlockRecordWriter::ChainInlineContext BlockRecordWriter::chain_inline_context(
    const HandleGraph& graph, const SiteChildren& children, const vector<Traversal>& travs,
    const vector<int>& genotype, int ref_trav_idx) const {
    ChainInlineContext ctx;
    // Only under block emission: with one record per snarl, no chain is inside a block.
    if (!enabled || !nested) {
        return ctx;
    }
    if (ref_trav_idx < 0 || (size_t)ref_trav_idx >= travs.size() || genotype.empty()) {
        return ctx;
    }
    // A snarl whose projection has no symbols cannot answer: every child would read as not
    // reported and be dropped.
    if (!children.known) {
        return ctx;
    }
    // A genotype with the reference allele matches every reference step, including the chain, so
    // the answer is false for every child.
    for (int allele : genotype) {
        if (allele == ref_trav_idx) {
            return ctx;
        }
    }

    ctx.sref = symbolic_allele(graph, travs[ref_trav_idx], children);

    for (int allele : genotype) {
        if (allele < 0 || (size_t)allele >= travs.size()) {
            continue;
        }
        ChainInlineContext::Alt alt;
        alt.salt = symbolic_allele(graph, travs[allele], children);
        alt.blocks = symbolic_diff(ctx.sref, alt.salt);
        ctx.alts.push_back(std::move(alt));
    }
    ctx.usable = true;
    return ctx;
}

bool BlockRecordWriter::chain_reported_inline(const HandleGraph& graph,
                                              const ChainInlineContext& ctx,
                                              const ChildChain& chain) const {
    if (!ctx.usable) {
        return false;
    }
    const pair<nid_t, nid_t> bounds(graph.get_id(chain.start), graph.get_id(chain.end));

    // Where the chain sits in the reference projection. If it is not there, the reference does not
    // cross it, which the caller handles.
    bool in_reference = false;
    for (size_t i = 0; i < ctx.sref.size(); ++i) {
        if (ctx.sref[i].is_chain() && ctx.sref[i].id == bounds.first &&
            ctx.sref[i].end_id == bounds.second) {
            in_reference = true;
            break;
        }
    }
    if (!in_reference) {
        return false;
    }

    bool any_crossing = false;
    for (const ChainInlineContext::Alt& alt : ctx.alts) {
        // This strand's own crossings of the chain, found in its own projection. A strand that
        // deletes the chain has none, so it has no copy for a block to report, even when the
        // reference's step for the chain lies inside a block.
        for (size_t j = 0; j < alt.salt.size(); ++j) {
            if (!alt.salt[j].is_chain() || alt.salt[j].id != bounds.first ||
                alt.salt[j].end_id != bounds.second) {
                continue;
            }
            any_crossing = true;
            bool inside = false;
            for (const DiffBlock& b : alt.blocks) {
                if ((size_t)b.alt_begin <= j && j < (size_t)b.alt_end) {
                    inside = true;
                    break;
                }
            }
            if (!inside) {
                return false;   // matched on this haplotype: its own record reports it
            }
        }
    }

    if (!any_crossing) {
        // No called strand crosses this chain, so no block ALT spells it.
        return false;
    }

    // Every crossing by every called strand falls inside a difference block, whose ALT spells the
    // route through the chain.
    ++counters.child_inlined;
    return true;
}

void BlockRecordWriter::count_site(const SiteChildren& children, const vector<Traversal>& travs,
                                   int ref_trav_idx) const {
    if (!nested || ref_trav_idx < 0 || (size_t)ref_trav_idx >= travs.size()) {
        return;
    }
    ++counters.sites;
    const bool site_reversed = children.reversed;
    if (!children.known) {
        // The projection would be a bare node list here, so the site is counted and skipped.
        ++counters.site_unresolvable;
        return;
    }
    if (site_reversed) {
        // Counted here, once per record, rather than in the resolver, which runs once per
        // projection.
        ++counters.site_reversed;
    }

}

int BlockRecordWriter::write(const PathPositionHandleGraph& graph, const SiteChildren& children,
                             const vector<Traversal>& called_traversals,
                             const vector<int>& genotype, int ref_trav_idx,
                             const string& sample_name, const NodeTranslation* translation,
                             const SiteRecord& record, GLLayout gl_layout, bool genotype_snarls,
                             const function<bool(vcflib::Variant&, size_t)>& add_line) const {
    const vcflib::Variant& site = record.variant;
    const map<int, int>& trav_to_allele = record.trav_to_allele;
    const int64_t site_position = record.unflattened_position;
    // Every refusal below returns -1, meaning the site record is written as it is. Block emission
    // being off is not a refusal, so it is not counted.
    if (!enabled || !nested || genotype_snarls) {
        return -1;
    }
    if (genotype.empty()) {
        ++counters.refuse[(size_t)AtomizeRefusal::NoGenotype];
        return -1;
    }
    if (record.alleles_merged) {
        ++counters.refuse[(size_t)AtomizeRefusal::AllelesMerged];
        return -1;
    }
    if (ref_trav_idx < 0 || (size_t)ref_trav_idx >= called_traversals.size()) {
        ++counters.refuse[(size_t)AtomizeRefusal::NoReferenceTraversal];
        return -1;
    }
    if (!children.known) {
        // The projection would see no child chains here.
        ++counters.refuse[(size_t)AtomizeRefusal::Unresolvable];
        return -1;
    }

    const Traversal& ref_trav = called_traversals[ref_trav_idx];
    vector<pair<int, int>> ref_ranges;
    SymbolicAllele sref = symbolic_allele(graph, ref_trav, children, &ref_ranges);
    const size_t m = sref.size();
    if (m == 0 || ref_ranges.size() != m) {
        ++counters.refuse[(size_t)AtomizeRefusal::EmptyReferenceProjection];
        return -1;
    }

    // Base offset of every visit boundary of the reference traversal from the snarl's first base.
    // The reference traversal is consecutive reference-path steps, so the running sum of node
    // lengths is the offset.
    vector<size_t> ref_visit_off(ref_trav.size() + 1, 0);
    for (size_t v = 0; v < ref_trav.size(); ++v) {
        ref_visit_off[v + 1] = ref_visit_off[v] + graph.get_length(ref_trav[v]);
    }

    auto visit_of_step = [](const vector<pair<int, int>>& ranges, size_t step,
                            const Traversal& t) -> int {
        // The ranges partition the walk's handles contiguously, so the handle index at step
        // boundary k is ranges[k].first, and one past the end is the walk's length.
        return step < ranges.size() ? ranges[step].first : (int)t.size();
    };
    // `max(vb, 0)`, so that the helper never reads t[-1], even though callers already refuse
    // vb <= 0.
    auto seq_of = [&](const Traversal& t, int vb, int ve) -> string {
        string s;
        for (int v = std::max(vb, 0); v < ve && v < (int)t.size(); ++v) {
            s += graph.get_sequence(t[v]);
        }
        return s;
    };

    struct HapAlign {
        int trav = -1;
        SymbolicAllele sym;
        vector<pair<int, int>> ranges;
        vector<DiffBlock> blocks;
        vector<int> alt_before_ref;
    };
    vector<HapAlign> haps(genotype.size());
    for (size_t s = 0; s < genotype.size(); ++s) {
        haps[s].trav = genotype[s];
        if (genotype[s] < 0 || (size_t)genotype[s] >= called_traversals.size()) {
            continue;   // a star or missing allele: nothing to align
        }
        if (genotype[s] == ref_trav_idx) {
            continue;   // the reference itself: every step matches, so no blocks
        }
        haps[s].sym = symbolic_allele(graph, called_traversals[genotype[s]], children,
                                      &haps[s].ranges);
        haps[s].blocks = symbolic_diff(sref, haps[s].sym, &haps[s].alt_before_ref);
        if (haps[s].alt_before_ref.size() != m + 1) {
            ++counters.refuse[(size_t)AtomizeRefusal::StepMapLength];
            return -1;
        }
    }

    // Cluster all strands' blocks by overlap of their reference step ranges, counting touching
    // ranges as overlapping, so that a deletion on one strand next to an insertion on the other is
    // one record: as two records, the alleles would disagree about the same reference span.
    vector<pair<int, int>> ivs;
    for (const HapAlign& h : haps) {
        for (const DiffBlock& b : h.blocks) {
            ivs.emplace_back(b.ref_begin, b.ref_end);
        }
    }
    if (ivs.empty()) {
        ++counters.refuse[(size_t)AtomizeRefusal::NoDifferenceBlocks];
        return -1;
    }
    sort(ivs.begin(), ivs.end());
    vector<pair<int, int>> clusters;
    clusters.push_back(ivs[0]);
    for (size_t k = 1; k < ivs.size(); ++k) {
        if (ivs[k].first <= clusters.back().second) {
            clusters.back().second = std::max(clusters.back().second, ivs[k].second);
        } else {
            clusters.push_back(ivs[k]);
        }
    }

    // Build every record before writing any, so that the decision to split is made on the finished
    // set.
    vector<vcflib::Variant> built;
    built.reserve(clusters.size());
    // For each record, whether no two of its alleles stand for one site allele.
    vector<bool> built_one_to_one;

    for (const pair<int, int>& cluster : clusters) {
        const size_t rb = (size_t)cluster.first;
        const size_t re = (size_t)cluster.second;
        int vb = visit_of_step(ref_ranges, rb, ref_trav);
        int ve = visit_of_step(ref_ranges, re, ref_trav);
        string ref_str = seq_of(ref_trav, vb, ve);

        // Each haplotype's allele over this same reference span, as the visits that spell it: its
        // own visits inside its difference blocks, and the reference's over the steps it matches.
        // A matched chain step is only the same chain, which the haplotype may cross by another
        // route, and that route is the chain's own record to report.
        vector<Traversal> slot_span(genotype.size());
        vector<string> slot_str(genotype.size());
        vector<bool> slot_marker(genotype.size(), false);
        auto append_visits = [](Traversal& span, const Traversal& t, int from, int to) {
            for (int v = std::max(from, 0); v < to && v < (int)t.size(); ++v) {
                span.push_back(t[v]);
            }
        };
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (haps[s].trav < 0) {
                slot_marker[s] = true;
                continue;
            }
            Traversal& span = slot_span[s];
            if (haps[s].trav == ref_trav_idx) {
                append_visits(span, ref_trav, vb, ve);
                slot_str[s] = ref_str;
                continue;
            }
            const Traversal& t = called_traversals[haps[s].trav];
            // Every block that overlaps or touches the cluster lies inside it, since the clusters
            // are the unions of all blocks, so `next` walks the cluster's reference steps in order.
            size_t next = rb;
            for (const DiffBlock& b : haps[s].blocks) {
                if ((size_t)b.ref_begin <= re && rb <= (size_t)b.ref_end) {
                    append_visits(span, ref_trav, visit_of_step(ref_ranges, next, ref_trav),
                                  visit_of_step(ref_ranges, (size_t)b.ref_begin, ref_trav));
                    append_visits(span, t, visit_of_step(haps[s].ranges, (size_t)b.alt_begin, t),
                                  visit_of_step(haps[s].ranges, (size_t)b.alt_end, t));
                    next = (size_t)b.ref_end;
                }
            }
            append_visits(span, ref_trav, visit_of_step(ref_ranges, next, ref_trav), ve);
            slot_str[s] = seq_of(span, 0, (int)span.size());
        }

        // VCF has no empty allele, so an indel takes the base before it, as flatten_common_allele_ends
        // leaves it.
        bool needs_anchor = ref_str.empty();
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (!slot_marker[s] && slot_str[s].empty()) {
                needs_anchor = true;
            }
        }
        int64_t pos = site_position + (int64_t)ref_visit_off[vb];
        if (needs_anchor) {
            if (vb <= 0) {
                ++counters.refuse[(size_t)AtomizeRefusal::NoAnchorBase];
                // Also stops seq_of(ref_trav, -1, 0) below from reading ref_trav[-1], which
                // can happen when the snarl's start node appears twice in the reference traversal.
                return -1;
            }
            string left = seq_of(ref_trav, vb - 1, vb);
            if (left.empty()) {
                ++counters.refuse[(size_t)AtomizeRefusal::AnchorWithoutSequence];
                return -1;
            }
            string base(1, left.back());
            ref_str = base + ref_str;
            for (size_t s = 0; s < genotype.size(); ++s) {
                if (!slot_marker[s]) {
                    slot_str[s] = base + slot_str[s];
                }
            }
            pos -= 1;
        }

        // Merge alleles with the same sequence, as the site record does, so two strands taking
        // different routes to the same sequence are homozygous.
        map<string, int> allele_to_gt;
        vector<string> alleles;
        vector<int> site_of_block;
        allele_to_gt[ref_str] = 0;
        alleles.push_back(ref_str);
        site_of_block.push_back(0);
        vector<int> block_gt(genotype.size(), MISSING_ALLELE_MARKER);
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (slot_marker[s]) {
                continue;
            }
            auto found = allele_to_gt.find(slot_str[s]);
            if (found != allele_to_gt.end()) {
                block_gt[s] = found->second;
                continue;
            }
            int a = (int)alleles.size();
            allele_to_gt[slot_str[s]] = a;
            alleles.push_back(slot_str[s]);
            // The site allele this block allele takes its evidence from: the one the same strand
            // carries at the site.
            auto sa = trav_to_allele.find(haps[s].trav);
            site_of_block.push_back(sa != trav_to_allele.end() ? sa->second : 0);
            block_gt[s] = a;
        }
        if (alleles.size() < 2) {
            continue;   // every haplotype is the reference here; nothing to report
        }

        vcflib::Variant b_var;
        b_var = site;   // inherit INFO, FORMAT and every site-level field
        b_var.position = pos;
        b_var.ref = alleles[0];
        b_var.alt.assign(alleles.begin() + 1, alleles.end());
        b_var.alleles = alleles;
        b_var.info["AT"].clear();
        b_var.info["AT"].resize(alleles.size());
        // The second field, the number of records, is known once every cluster is built.
        b_var.info["SB"] = {std::to_string(built.size()), ""};

        // AT per block allele, over the visit range this record actually spells.
        {
            Traversal ref_span;
            for (int v = vb; v < ve && v < (int)ref_trav.size(); ++v) {
                ref_span.push_back(ref_trav[v]);
            }
            add_allele_path_to_info(b_var, 0, visits_of(graph, ref_span), false, translation);
            for (size_t a = 1; a < alleles.size(); ++a) {
                Traversal span;
                for (size_t s = 0; s < genotype.size(); ++s) {
                    if (!slot_marker[s] && block_gt[s] == (int)a) {
                        span = slot_span[s];
                        break;
                    }
                }
                add_allele_path_to_info(b_var, a, visits_of(graph, span), false, translation);
            }
        }

        // GT, keeping the site's phase. The site's GT is in site allele space and this record's
        // slots are in genotyper order, so each slot is mapped to the site allele its strand
        // carries. Three forms carry a phase set:
        //
        //   - "a|b", a phased diploid pair, whose order carries over by slot.
        //   - "a|." or ".|a", a nested chain at ploidy 1, on one strand of its parent. A block of
        //     it is part of the same allele on the same strand, so the strand carries over.
        //   - "a", a haploid locus (chrY, or chrX outside the pseudoautosomal regions), with no
        //     order, but PS still labels its phase set.
        //
        // Only a slash-separated GT is unphased, and only there is PS removed.
        bool keep_phase_set = false;
        int nested_strand = -1;   // haploid record: which side of "a|." this allele sits on
        vector<int> slot_order(block_gt.size());
        for (size_t s = 0; s < block_gt.size(); ++s) {
            slot_order[s] = (int)s;
        }
        const string* site_gt_text = nullptr;
        {
            auto site_gt = site.samples.find(sample_name);
            if (site_gt != site.samples.end()) {
                auto gt_field = site_gt->second.find("GT");
                if (gt_field != site_gt->second.end() && !gt_field->second.empty()) {
                    site_gt_text = &gt_field->second[0];
                }
            }
        }
        if (site_gt_text != nullptr && site_gt_text->find('/') == string::npos) {
            const size_t bar = site_gt_text->find('|');
            if (bar == string::npos) {
                keep_phase_set = block_gt.size() == 1 && block_gt[0] >= 0;
            } else {
                const string left = site_gt_text->substr(0, bar);
                const string right = site_gt_text->substr(bar + 1);
                if (block_gt.size() == 1 && block_gt[0] >= 0 && (left == "." || right == ".")) {
                    keep_phase_set = true;
                    nested_strand = left == "." ? 1 : 0;
                } else if (block_gt.size() == 2 && block_gt[0] >= 0 && block_gt[1] >= 0) {
                    int a = -1;
                    int b2 = -1;
                    try {
                        a = std::stoi(left);
                        b2 = std::stoi(right);
                    } catch (const std::exception&) {
                        a = -1;
                    }
                    auto site_allele_of_slot = [&](size_t sl) {
                        auto found = trav_to_allele.find(haps[sl].trav);
                        return found != trav_to_allele.end() ? found->second : -1;
                    };
                    int s0 = site_allele_of_slot(0);
                    int s1 = site_allele_of_slot(1);
                    if (a >= 0 && s0 >= 0 && s1 >= 0) {
                        if (s0 == a && s1 == b2) {
                            keep_phase_set = true;
                        } else if (s0 == b2 && s1 == a) {
                            keep_phase_set = true;
                            slot_order[0] = 1;
                            slot_order[1] = 0;
                        }
                    }
                }
            }
        }
        {
            string gt;
            if (nested_strand >= 0) {
                gt = nested_strand_genotype(block_gt[0], nested_strand);
            } else {
                for (size_t s = 0; s < block_gt.size(); ++s) {
                    const int value = block_gt[slot_order[s]];
                    gt += value == MISSING_ALLELE_MARKER ? string(".") : std::to_string(value);
                    if (s + 1 != block_gt.size()) {
                        gt += keep_phase_set ? '|' : '/';
                    }
                }
            }
            b_var.samples[sample_name]["GT"] = {gt};
        }
        if (!keep_phase_set) {
            // PS labels a phase block, so it means nothing on an unphased genotype.
            b_var.samples[sample_name].erase("PS");
            b_var.format.erase(std::remove(b_var.format.begin(), b_var.format.end(), "PS"),
                               b_var.format.end());
        }

        // AD and GL are taken from the site, so every block of a snarl reports the same evidence;
        // INFO/SB lets a consumer avoid counting it more than once.
        auto& fmt = b_var.samples[sample_name];
        auto site_fmt = site.samples.find(sample_name);
        if (site_fmt != site.samples.end()) {
            auto site_ad = site_fmt->second.find("AD");
            if (site_ad != site_fmt->second.end()) {
                vector<string> ad;
                for (size_t a = 0; a < alleles.size(); ++a) {
                    size_t si = (size_t)site_of_block[a];
                    ad.push_back(si < site_ad->second.size() ? site_ad->second[si] : "0");
                }
                fmt["AD"] = ad;
            }
            auto site_gl = site_fmt->second.find("GL");
            bool has_marker = std::any_of(block_gt.begin(), block_gt.end(),
                                          [](int a) { return a < 0; });
            if (site_gl != site_fmt->second.end() && !has_marker) {
                const size_t k = alleles.size();
                const size_t ks = site.alleles.size();
                vector<string> gl;
                bool ok = true;
                if (block_gt.size() == 1) {
                    gl.resize(k);
                    for (size_t a = 0; a < k && ok; ++a) {
                        size_t si = (size_t)site_of_block[a];
                        ok = si < site_gl->second.size();
                        if (ok) {
                            gl[a] = site_gl->second[si];
                        }
                    }
                } else if (block_gt.size() == 2) {
                    gl.resize(k * (k + 1) / 2);
                    for (size_t j = 0; j < k && ok; ++j) {
                        for (size_t i = 0; i <= j && ok; ++i) {
                            size_t si = (size_t)site_of_block[i];
                            size_t sj = (size_t)site_of_block[j];
                            size_t src = gl_genotype_index(std::min(si, sj), std::max(si, sj), ks,
                                                           gl_layout);
                            size_t dst = gl_genotype_index(i, j, k, gl_layout);
                            ok = src < site_gl->second.size() && dst < gl.size();
                            if (ok) {
                                gl[dst] = site_gl->second[src];
                            }
                        }
                    }
                } else {
                    ok = false;
                }
                if (ok) {
                    fmt["GL"] = gl;
                } else {
                    fmt.erase("GL");
                    b_var.format.erase(std::remove(b_var.format.begin(), b_var.format.end(), "GL"),
                                       b_var.format.end());
                }
            }
        }

        b_var.updateAlleleIndexes();
        vg::flatten_common_allele_ends(b_var, true, 0);
        vg::flatten_common_allele_ends(b_var, false, 0);
        built.push_back(std::move(b_var));
        built_one_to_one.push_back(set<int>(site_of_block.begin(), site_of_block.end()).size()
                                   == site_of_block.size());
    }

    if (built.empty()) {
        ++counters.refuse[(size_t)AtomizeRefusal::ReferenceBasesOnly];
        return -1;
    }
    // One block replaces the site record only where the site record says more than the block. A
    // block spells a strand's own visits inside its difference blocks and the reference's over the
    // steps the strand matches, so where a strand crosses a matched child chain, the block spells the
    // reference's route through it, and the chain's own record reports the strand's route. The site
    // record spells each strand's site allele in full, so it says more exactly where some strand's
    // allele takes a route through a matched chain that spells other bases, and there it would repeat
    // the chain's record.
    if (built.size() == 1) {
        // A chain is genotyped from each allele's first crossing of it, so where the reference or a
        // strand crosses a chain more than once, the chain's record need not report the crossing the
        // block leaves to it, and the site record stands.
        auto crosses_a_chain_twice = [](const SymbolicAllele& sym) {
            set<pair<nid_t, nid_t>> chains;
            for (const SymbolicStep& step : sym) {
                if (step.is_chain() && !chains.insert(make_pair(step.id, step.end_id)).second) {
                    return true;
                }
            }
            return false;
        };
        bool chain_crossed_twice = crosses_a_chain_twice(sref);
        bool site_says_more = false;
        for (size_t s = 0; s < genotype.size(); ++s) {
            if (haps[s].trav < 0 || haps[s].trav == ref_trav_idx) {
                continue;
            }
            chain_crossed_twice = chain_crossed_twice || crosses_a_chain_twice(haps[s].sym);
            const Traversal& t = called_traversals[haps[s].trav];
            string as_blocks;
            size_t next = 0;
            for (const DiffBlock& b : haps[s].blocks) {
                as_blocks += seq_of(ref_trav, visit_of_step(ref_ranges, next, ref_trav),
                                    visit_of_step(ref_ranges, (size_t)b.ref_begin, ref_trav));
                as_blocks += seq_of(t, visit_of_step(haps[s].ranges, (size_t)b.alt_begin, t),
                                    visit_of_step(haps[s].ranges, (size_t)b.alt_end, t));
                next = (size_t)b.ref_end;
            }
            as_blocks += seq_of(ref_trav, visit_of_step(ref_ranges, next, ref_trav),
                                (int)ref_trav.size());
            // The site record spells the strand's site allele, which is the reference for a route that
            // differs from it only inside child chains.
            auto allele = trav_to_allele.find(haps[s].trav);
            const string site_allele =
                allele != trav_to_allele.end() && allele->second == 0
                    ? seq_of(ref_trav, visit_of_step(ref_ranges, 0, ref_trav), (int)ref_trav.size())
                    : seq_of(t, visit_of_step(haps[s].ranges, 0, t), (int)t.size());
            site_says_more = site_says_more || as_blocks != site_allele;
        }
        if (!site_says_more) {
            ++counters.refuse[(size_t)AtomizeRefusal::SameAsSiteRecord];
            return -1;
        }
        if (chain_crossed_twice) {
            ++counters.refuse[(size_t)AtomizeRefusal::ChainCrossedTwice];
            return -1;
        }
        // The block takes its genotype, likelihoods and phase from the site, through the site allele
        // each of its alleles stands for, so it cannot be written where two of its alleles stand for
        // one. Two routes can spell one site allele and still differ here, where a difference outside
        // a child chain is cancelled by one inside it.
        if (!built_one_to_one[0]) {
            ++counters.refuse[(size_t)AtomizeRefusal::StrandsDisagree];
            return -1;
        }
    }

    int added = 0;
    size_t block_index = 0;
    for (vcflib::Variant& b_var : built) {
        // The number of records the site writes, which leaves out a cluster that every strand
        // spells as the reference.
        b_var.info["SB"][1] = std::to_string(built.size());
        // VCF allows an ID on one record only, so each block record takes the site's ID with its
        // index appended. A snarl's name has no '_', so the site's ID can be read back from it
        // (see block_site_name).
        b_var.id += "_" + b_var.info["SB"][0];
        if (add_line(b_var, block_index)) {
            ++added;
        }
        ++block_index;
    }
    ++counters.split_sites;
    counters.split_lines += (size_t)added;
    return added;
}

}

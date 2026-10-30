#include <algorithm>
#include <array>
#include <limits>
#include <optional>
#include <stdexcept>

#include "block_records.hpp"
#include "../src/utility.hpp"

namespace vg::multipass {

// The names of the AtomizeRefusal reasons, in order. The initializer sets the size, so that the
// check below fails when a name is missing as well as when one is extra.
static constexpr auto g_atomize_refuse_name = std::to_array<const char*>({
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
});
static_assert(g_atomize_refuse_name.size() == (size_t)AtomizeRefusal::Count,
              "each AtomizeRefusal must have a name in g_atomize_refuse_name, "
              "and each name must have a reason");

namespace {

/// Return whether `step` is the symbol of the chain with boundary nodes `chain_start` and
/// `chain_end`.
bool is_chain_step(const SymbolicStep& step, nid_t chain_start, nid_t chain_end) {
    return step.is_chain() && step.id == chain_start && step.end_id == chain_end;
}

/// Return whether step `j` of an ALT projection lies inside one of `blocks`.
bool in_a_block(const vector<DiffBlock>& blocks, size_t j) {
    return std::ranges::any_of(blocks, [&](const DiffBlock& b) {
        return (size_t)b.alt_begin <= j && j < (size_t)b.alt_end;
    });
}

/// Return whether a projection crosses one chain more than once.
bool crosses_a_chain_twice(const SymbolicAllele& sym) {
    set<pair<nid_t, nid_t>> chains;
    for (const SymbolicStep& step : sym) {
        if (step.is_chain() && !chains.insert(make_pair(step.id, step.end_id)).second) {
            return true;
        }
    }
    return false;
}

/// Return the index of the handle of walk `t` at which step boundary `step` of its projection
/// falls.
///
/// `ranges` gives the range of handles that each step of the projection covers. The ranges cover
/// the walk in order, so boundary k falls at the first handle of step k, and the boundary after
/// the last step falls at the end of the walk.
int visit_of_step(const vector<pair<int, int>>& ranges, size_t step, const Traversal& t) {
    return step < ranges.size() ? ranges[step].first : (int)t.size();
}

/// One slot of the genotype of a site, aligned to the reference traversal by symbolic projection.
struct StrandAlignment {
    /// The traversal the slot carries, or an allele marker.
    int trav = -1;
    /// The projection of the traversal.
    SymbolicAllele sym;
    /// For each step of `sym`, the range of handles of the traversal it covers.
    vector<pair<int, int>> ranges;
    /// The difference blocks between the reference projection and `sym`, in reference order.
    vector<DiffBlock> blocks;
};

/// What the reference and each slot of the genotype spell over one cluster of difference blocks.
struct ClusterSpelling {
    /// The first visit of the reference traversal in the cluster.
    int vb = 0;
    /// The visit of the reference traversal after the cluster.
    int ve = 0;
    /// The sequence of the reference over the cluster.
    string ref_str;
    /// For each slot, the visits that spell its allele over the cluster.
    vector<Traversal> slot_span;
    /// For each slot, the sequence of its allele over the cluster.
    vector<string> slot_str;
    /// For each slot, whether it holds an allele marker, which spells nothing.
    vector<bool> slot_marker;
};

/**
 * The slots of the genotype of a site, each aligned to the reference traversal by symbolic
 * projection.
 *
 * Block records are cut from this alignment. Its positions are boundaries between the steps of
 * the reference projection, which `visit_of_step` turns into visits of the reference traversal.
 */
class SiteAlignment {
public:
    /// Align each slot of `genotype` to the reference traversal `called_traversals[ref_trav_idx]`.
    ///
    /// If the alignment cannot be made, `refusal()` says why. An allele of a slot is spelled by
    /// its step sequences in `sequences`, where it has them.
    ///
    /// @throws std::out_of_range if an entry of `genotype` is neither an allele marker nor the
    /// number of a traversal in `called_traversals`.
    SiteAlignment(const PathPositionHandleGraph& graph, const SiteChildren& children,
                  const vector<Traversal>& called_traversals, const vector<int>& genotype,
                  int ref_trav_idx, const AlleleSequences* sequences)
        : graph(graph), called_traversals(called_traversals), ref_trav_idx(ref_trav_idx),
          sequences(sequences), ref_trav(called_traversals[ref_trav_idx]) {
        sref = symbolic_allele(graph, ref_trav, children, &ref_ranges);
        if (sref.empty() || ref_ranges.size() != sref.size()) {
            failure = AtomizeRefusal::EmptyReferenceProjection;
            return;
        }
        // The reference traversal is consecutive steps of the reference path, so the running sum
        // of node lengths is the offset of each visit from the first base of the site.
        ref_visit_off.assign(ref_trav.size() + 1, 0);
        for (size_t v = 0; v < ref_trav.size(); ++v) {
            ref_visit_off[v + 1] = ref_visit_off[v] + graph.get_length(ref_trav[v]);
        }
        slots.resize(genotype.size());
        for (size_t s = 0; s < genotype.size(); ++s) {
            slots[s].trav = genotype[s];
            if (genotype[s] >= 0 && (size_t)genotype[s] >= called_traversals.size()) {
                throw std::out_of_range("genotype entry " + std::to_string(genotype[s])
                                        + " names no called traversal");
            }
            // A marker has nothing to align, and the reference matches itself everywhere.
            if (genotype[s] < 0 || genotype[s] == ref_trav_idx) {
                continue;
            }
            slots[s].sym = symbolic_allele(graph, called_traversals[genotype[s]], children,
                                           &slots[s].ranges);
            vector<int> alt_before_ref;
            slots[s].blocks = symbolic_diff(sref, slots[s].sym, &alt_before_ref);
            if (alt_before_ref.size() != sref.size() + 1) {
                failure = AtomizeRefusal::StepMapLength;
                return;
            }
        }
    }

    /// Get why the alignment could not be made, if it could not.
    optional<AtomizeRefusal> refusal() const { return failure; }

    /// Get the number of slots of the genotype.
    size_t slot_count() const { return slots.size(); }

    /// Get the traversal slot `s` carries, or its allele marker.
    int traversal_of(size_t s) const { return slots[s].trav; }

    /// Get the offset of visit `v` of the reference traversal from the first base of the site.
    size_t ref_offset(int v) const { return ref_visit_off[v]; }

    /**
     * Return the clusters of the difference blocks of all slots, as ranges of reference steps.
     *
     * Blocks whose ranges overlap or touch go in one cluster, so that a deletion on one strand
     * next to an insertion on the other is one record. As two records, their alleles would
     * disagree about the same reference span. The clusters are in reference order.
     */
    vector<pair<int, int>> difference_clusters() const {
        vector<pair<int, int>> ivs;
        for (const StrandAlignment& slot : slots) {
            for (const DiffBlock& b : slot.blocks) {
                ivs.emplace_back(b.ref_begin, b.ref_end);
            }
        }
        if (ivs.empty()) {
            return ivs;
        }
        std::ranges::sort(ivs);
        vector<pair<int, int>> clusters;
        clusters.push_back(ivs[0]);
        for (size_t k = 1; k < ivs.size(); ++k) {
            if (ivs[k].first <= clusters.back().second) {
                clusters.back().second = std::max(clusters.back().second, ivs[k].second);
            } else {
                clusters.push_back(ivs[k]);
            }
        }
        return clusters;
    }

    /**
     * Return what the reference and each slot spell over the cluster from reference step `rb` to
     * reference step `re`.
     *
     * A slot spells its own visits inside its difference blocks, and the visits of the reference
     * over the steps it matches. A matched chain step is only the same chain, which the slot may
     * cross by another route, and that route is for the record of the chain to report.
     */
    ClusterSpelling spell_cluster(size_t rb, size_t re) const {
        ClusterSpelling spelling;
        spelling.vb = ref_visit(rb);
        spelling.ve = ref_visit(re);
        spelling.ref_str = spell(ref_trav, spelling.vb, spelling.ve);
        spelling.slot_span.resize(slots.size());
        spelling.slot_str.resize(slots.size());
        spelling.slot_marker.assign(slots.size(), false);
        for (size_t s = 0; s < slots.size(); ++s) {
            if (slots[s].trav < 0) {
                spelling.slot_marker[s] = true;
                continue;
            }
            Traversal& span = spelling.slot_span[s];
            if (slots[s].trav == ref_trav_idx) {
                append_visits(span, ref_trav, spelling.vb, spelling.ve);
                spelling.slot_str[s] = spelling.ref_str;
                continue;
            }
            spelling.slot_str[s] = spell_slot_over(s, rb, re, span);
        }
        return spelling;
    }

    /// Return the sequence of visits [vb, ve) of the reference traversal.
    string spell_reference(int vb, int ve) const { return spell(ref_trav, vb, ve); }

    /// Return visits [vb, ve) of the reference traversal.
    Traversal reference_visits(int vb, int ve) const {
        Traversal ref_span;
        for (int v = vb; v < ve && v < (int)ref_trav.size(); ++v) {
            ref_span.push_back(ref_trav[v]);
        }
        return ref_span;
    }

    /// Return whether the reference, or the allele of a slot, crosses one chain more than once.
    bool crosses_a_chain_twice() const {
        if (multipass::crosses_a_chain_twice(sref)) {
            return true;
        }
        return std::ranges::any_of(slots, [&](const StrandAlignment& slot) {
            return slot.trav >= 0 && slot.trav != ref_trav_idx
                   && multipass::crosses_a_chain_twice(slot.sym);
        });
    }

    /**
     * Return whether the site record spells the allele of some slot differently from one block.
     *
     * The site record spells the site allele of each slot in full. A block spells the visits of
     * the reference over the steps the slot matches, including any child chain there, whose own
     * record reports the route of the slot through it. So the site record says more exactly where
     * the route of some slot through a matched chain spells other bases. `trav_to_allele` gives
     * the site allele of each traversal.
     */
    bool site_says_more(const map<int, int>& trav_to_allele) const {
        bool says_more = false;
        for (size_t s = 0; s < slots.size(); ++s) {
            if (slots[s].trav < 0 || slots[s].trav == ref_trav_idx) {
                continue;
            }
            const Traversal& t = called_traversals[slots[s].trav];
            const StepSequences* own = own_sequences(sequences, slots[s].trav);
            string as_blocks;
            size_t next = 0;
            for (const DiffBlock& b : slots[s].blocks) {
                as_blocks += spell(ref_trav, ref_visit(next), ref_visit((size_t)b.ref_begin));
                as_blocks += spell(t, slot_visit(s, (size_t)b.alt_begin),
                                   slot_visit(s, (size_t)b.alt_end), own);
                next = (size_t)b.ref_end;
            }
            as_blocks += spell(ref_trav, ref_visit(next), (int)ref_trav.size());
            // A route that differs from the reference only inside child chains has the reference
            // as its site allele.
            auto allele = trav_to_allele.find(slots[s].trav);
            const string site_allele =
                allele != trav_to_allele.end() && allele->second == 0
                    ? spell(ref_trav, ref_visit(0), (int)ref_trav.size())
                    : spell(t, slot_visit(s, 0), (int)t.size(), own);
            says_more = says_more || as_blocks != site_allele;
        }
        return says_more;
    }

private:
    /// Return the visit of the reference traversal at reference step boundary `step`.
    int ref_visit(size_t step) const { return visit_of_step(ref_ranges, step, ref_trav); }

    /// Return the visit of the traversal of slot `s` at boundary `step` of its projection.
    int slot_visit(size_t s, size_t step) const {
        return visit_of_step(slots[s].ranges, step, called_traversals[slots[s].trav]);
    }

    /// Return the sequence of visits [vb, ve) of walk `t`, spelled by `own` where it is given.
    string spell(const Traversal& t, int vb, int ve, const StepSequences* own = nullptr) const {
        string s;
        // Start at 0 even if `vb` is negative, so that we never read t[-1].
        for (int v = std::max(vb, 0); v < ve && v < (int)t.size(); ++v) {
            s += own != nullptr ? own->at(v) : graph.get_sequence(t[v]);
        }
        return s;
    }

    /// Append visits [from, to) of walk `t` to `span`, and what they spell to `spelled`.
    void append_visits(Traversal& span, const Traversal& t, int from, int to,
                       string* spelled = nullptr, const StepSequences* own = nullptr) const {
        for (int v = std::max(from, 0); v < to && v < (int)t.size(); ++v) {
            span.push_back(t[v]);
            if (spelled != nullptr) {
                *spelled += own != nullptr ? own->at(v) : graph.get_sequence(t[v]);
            }
        }
    }

    /// Append to `span` the visits that spell the ALT of slot `s` over the cluster of reference
    /// steps [rb, re], and return what they spell.
    string spell_slot_over(size_t s, size_t rb, size_t re, Traversal& span) const {
        const Traversal& t = called_traversals[slots[s].trav];
        const StepSequences* own = own_sequences(sequences, slots[s].trav);
        string spelled;
        // Every block that overlaps or touches the cluster lies inside it, since the clusters are
        // the unions of all blocks, so `next` walks the reference steps of the cluster in order.
        size_t next = rb;
        for (const DiffBlock& b : slots[s].blocks) {
            if ((size_t)b.ref_begin <= re && rb <= (size_t)b.ref_end) {
                append_visits(span, ref_trav, ref_visit(next), ref_visit((size_t)b.ref_begin),
                              &spelled);
                append_visits(span, t, slot_visit(s, (size_t)b.alt_begin),
                              slot_visit(s, (size_t)b.alt_end), &spelled, own);
                next = (size_t)b.ref_end;
            }
        }
        append_visits(span, ref_trav, ref_visit(next), ref_visit(re), &spelled);
        return spelled;
    }

    const PathPositionHandleGraph& graph;
    /// The traversals the genotype numbers.
    const vector<Traversal>& called_traversals;
    /// The number of the reference traversal.
    int ref_trav_idx;
    /// The step sequences of the alleles that have them, or null.
    const AlleleSequences* sequences;
    /// The reference traversal.
    const Traversal& ref_trav;
    /// The projection of the reference traversal.
    SymbolicAllele sref;
    /// For each step of `sref`, the range of visits of the reference traversal it covers.
    vector<pair<int, int>> ref_ranges;
    /// For each visit boundary of the reference traversal, its offset from the first base of the
    /// site.
    vector<size_t> ref_visit_off;
    /// The alignment of each slot of the genotype.
    vector<StrandAlignment> slots;
    /// Why the alignment could not be made, if it could not.
    optional<AtomizeRefusal> failure;
};

/// Return whether the alleles of a cluster need an anchor base, which VCF requires of an indel,
/// since it has no empty allele.
bool needs_anchor_base(const ClusterSpelling& spelling) {
    bool needs_anchor = spelling.ref_str.empty();
    for (size_t s = 0; s < spelling.slot_str.size(); ++s) {
        if (!spelling.slot_marker[s] && spelling.slot_str[s].empty()) {
            needs_anchor = true;
        }
    }
    return needs_anchor;
}

/// The alleles of a block record, numbered as VCF numbers them, and its genotype.
struct BlockAlleles {
    /// The sequence of each allele, the reference first.
    vector<string> alleles;
    /// For each allele, the site allele it takes its evidence from: the one the same strand
    /// carries at the site.
    vector<int> site_of_block;
    /// For each slot of the genotype, its allele, or MISSING_ALLELE_MARKER.
    vector<int> block_gt;

    /// Return whether no two alleles stand for one site allele.
    bool one_to_one() const {
        return set<int>(site_of_block.begin(), site_of_block.end()).size() == site_of_block.size();
    }
};

/// Number the alleles of a block record, giving slots that spell the same sequence one allele,
/// as the site record does.
///
/// Two strands that take different routes to the same sequence are then homozygous.
/// `trav_to_allele` gives the site allele of each traversal.
BlockAlleles number_block_alleles(const ClusterSpelling& spelling, const SiteAlignment& alignment,
                                  const map<int, int>& trav_to_allele) {
    BlockAlleles block;
    map<string, int> allele_to_gt;
    allele_to_gt[spelling.ref_str] = 0;
    block.alleles.push_back(spelling.ref_str);
    block.site_of_block.push_back(0);
    block.block_gt.assign(alignment.slot_count(), MISSING_ALLELE_MARKER);
    for (size_t s = 0; s < alignment.slot_count(); ++s) {
        if (spelling.slot_marker[s]) {
            continue;
        }
        if (auto found = allele_to_gt.find(spelling.slot_str[s]); found != allele_to_gt.end()) {
            block.block_gt[s] = found->second;
            continue;
        }
        auto a = (int)block.alleles.size();
        allele_to_gt[spelling.slot_str[s]] = a;
        block.alleles.push_back(spelling.slot_str[s]);
        auto sa = trav_to_allele.find(alignment.traversal_of(s));
        block.site_of_block.push_back(sa != trav_to_allele.end() ? sa->second : 0);
        block.block_gt[s] = a;
    }
    return block;
}

/// Return the first GT entry of `sample_name` in `site`, or null if it has none.
const string* site_gt_text(const vcflib::Variant& site, const string& sample_name) {
    auto site_gt = site.samples.find(sample_name);
    if (site_gt == site.samples.end()) {
        return nullptr;
    }
    auto gt_field = site_gt->second.find("GT");
    if (gt_field == site_gt->second.end() || gt_field->second.empty()) {
        return nullptr;
    }
    return &gt_field->second[0];
}

/// Match the slots of a diploid block to the phased site GT "left|right".
///
/// `s0` and `s1` are the site alleles the strands of slots 0 and 1 carry, or -1. Return false if
/// the slots are in the order of the GT, true if they are swapped, or nothing if they do not match
/// it.
optional<bool> slots_swapped(const string& left, const string& right, int s0, int s1) {
    int a = -1;
    int b = -1;
    try {
        a = std::stoi(left);
        b = std::stoi(right);
    } catch (const std::logic_error&) {
        // std::stoi throws std::invalid_argument or std::out_of_range.
        a = -1;
    }
    if (a < 0 || s0 < 0 || s1 < 0) {
        return std::nullopt;
    }
    if (s0 == a && s1 == b) {
        return false;
    }
    if (s0 == b && s1 == a) {
        return true;
    }
    return std::nullopt;
}

/// How the phase of a site GT carries over to a block record.
struct BlockPhase {
    /// Whether the block GT is phased, so that PS still applies.
    bool keep_phase_set = false;
    /// For a nested ploidy-1 block, which side of "a|." its allele sits on, or -1.
    int nested_strand = -1;
    /// For a diploid block, whether its two slots are in the reverse order of the site GT.
    bool swapped = false;
};

/// Work out how the phase of the site GT `site_gt` carries over to a block with genotype
/// `block_gt`, as `phase_block_genotype` describes.
BlockPhase carry_site_phase(const string* site_gt, const vector<int>& block_gt,
                            const vector<int>& slot_site_allele) {
    BlockPhase phase;
    if (site_gt == nullptr || site_gt->find('/') != string::npos) {
        return phase;
    }
    const bool one_allele = block_gt.size() == 1 && block_gt[0] >= 0;
    const size_t bar = site_gt->find('|');
    if (bar == string::npos) {
        phase.keep_phase_set = one_allele;
        return phase;
    }
    const string left = site_gt->substr(0, bar);
    const string right = site_gt->substr(bar + 1);
    if (one_allele && (left == "." || right == ".")) {
        phase.keep_phase_set = true;
        phase.nested_strand = left == "." ? 1 : 0;
    } else if (block_gt.size() == 2 && block_gt[0] >= 0 && block_gt[1] >= 0) {
        const optional<bool> swapped =
            slots_swapped(left, right, slot_site_allele[0], slot_site_allele[1]);
        phase.keep_phase_set = swapped.has_value();
        phase.swapped = swapped.value_or(false);
    }
    return phase;
}

/// The GT of a block record, and whether the record keeps the phase set of the site.
struct BlockGenotype {
    /// The GT text.
    string gt;
    /// Whether the GT is phased, so that PS still applies.
    bool keep_phase_set = false;
};

/**
 * Return the GT of a block record, keeping the phase of the site GT `site_gt` where it can.
 *
 * The site GT is in site alleles, and the slots of the block are in genotyper order, so each
 * slot is matched to the site allele its strand carries, which `slot_site_allele` gives. Three
 * forms of site GT carry a phase set:
 *
 *   - "a|b", a phased diploid pair, whose order carries over by slot.
 *   - "a|." or ".|a", a nested chain at ploidy 1, on one strand of its parent. A block of it is
 *     part of the same allele on the same strand, so the strand carries over.
 *   - "a", a haploid locus (chrY, or chrX outside the pseudoautosomal regions), with no order,
 *     whose PS still labels its phase set.
 *
 * A GT with "/" is unphased.
 */
BlockGenotype phase_block_genotype(const string* site_gt, const vector<int>& block_gt,
                                   const vector<int>& slot_site_allele) {
    const BlockPhase phase = carry_site_phase(site_gt, block_gt, slot_site_allele);
    BlockGenotype result;
    result.keep_phase_set = phase.keep_phase_set;
    if (phase.nested_strand >= 0) {
        result.gt = nested_strand_genotype(block_gt[0], phase.nested_strand);
        return result;
    }
    for (size_t s = 0; s < block_gt.size(); ++s) {
        // With the slots swapped, the two entries of a diploid GT are written in reverse.
        const size_t slot = phase.swapped ? block_gt.size() - 1 - s : s;
        const int value = block_gt[slot];
        result.gt += value == MISSING_ALLELE_MARKER ? string(".") : std::to_string(value);
        if (s + 1 != block_gt.size()) {
            result.gt += result.keep_phase_set ? '|' : '/';
        }
    }
    return result;
}

/// Return the GL of a block record, taken from the GL of the site, or nothing if an entry it
/// needs is missing.
///
/// `site_allele_count` is the number of alleles of the site, and `gl_layout` the order of both
/// GLs. Only a haploid or diploid genotype has a GL.
optional<vector<string>> block_likelihoods(const vector<string>& site_gl, size_t site_allele_count,
                                           const BlockAlleles& block, GLLayout gl_layout) {
    const size_t k = block.alleles.size();
    vector<string> gl;
    if (block.block_gt.size() == 1) {
        gl.resize(k);
        for (size_t a = 0; a < k; ++a) {
            auto si = (size_t)block.site_of_block[a];
            if (si >= site_gl.size()) {
                return std::nullopt;
            }
            gl[a] = site_gl[si];
        }
        return gl;
    }
    if (block.block_gt.size() != 2) {
        return std::nullopt;
    }
    gl.resize(k * (k + 1) / 2);
    for (size_t j = 0; j < k; ++j) {
        for (size_t i = 0; i <= j; ++i) {
            auto si = (size_t)block.site_of_block[i];
            auto sj = (size_t)block.site_of_block[j];
            size_t src = gl_genotype_index(std::min(si, sj), std::max(si, sj), site_allele_count,
                                           gl_layout);
            size_t dst = gl_genotype_index(i, j, k, gl_layout);
            if (src >= site_gl.size() || dst >= gl.size()) {
                return std::nullopt;
            }
            gl[dst] = site_gl[src];
        }
    }
    return gl;
}

/// Copy AD and GL from the site onto a block record, through the site allele of each block
/// allele.
///
/// Every block of a site so reports the same evidence. INFO/SB lets a consumer count it once.
void copy_site_evidence(vcflib::Variant& b_var, const vcflib::Variant& site,
                        const string& sample_name, const BlockAlleles& block,
                        GLLayout gl_layout) {
    auto& fmt = b_var.samples[sample_name];
    auto site_fmt = site.samples.find(sample_name);
    if (site_fmt == site.samples.end()) {
        return;
    }
    if (auto site_ad = site_fmt->second.find("AD"); site_ad != site_fmt->second.end()) {
        vector<string> ad;
        for (size_t a = 0; a < block.alleles.size(); ++a) {
            auto si = (size_t)block.site_of_block[a];
            ad.push_back(si < site_ad->second.size() ? site_ad->second[si] : "0");
        }
        fmt["AD"] = ad;
    }
    auto site_gl = site_fmt->second.find("GL");
    if (site_gl == site_fmt->second.end()
        || std::ranges::any_of(block.block_gt, [](int a) { return a < 0; })) {
        return;
    }
    optional<vector<string>> gl =
        block_likelihoods(site_gl->second, site.alleles.size(), block, gl_layout);
    if (gl) {
        fmt["GL"] = *gl;
    } else {
        fmt.erase("GL");
        std::erase(b_var.format, "GL");
    }
}

/**
 * A cutter that cuts the record of a site into block records, one for each cluster of difference
 * blocks that some strand spells differently from the reference.
 *
 * Each block record copies every field of the site record that it does not redefine.
 */
class BlockCutter {
public:
    /// Make a cutter for the site record `record`, aligned as `alignment`.
    ///
    /// The block records are for `sample_name`, with nodes named by `translation` if it is not
    /// null. `gl_layout` is the order of GL in the site record.
    BlockCutter(const PathPositionHandleGraph& graph, const SiteAlignment& alignment,
                const SiteRecord& record, const string& sample_name,
                const NodeTranslation* translation, GLLayout gl_layout)
        : graph(graph), alignment(alignment), record(record), sample_name(sample_name),
          translation(translation), gl_layout(gl_layout) {}

    /// Cut each of `clusters`, from `difference_clusters()` of the alignment, into a block record,
    /// and add it to `built`. Return why not, if a cluster cannot be cut.
    ///
    /// A cluster every strand spells as the reference gets no record. For each record, add to
    /// `one_to_one` whether no two of its alleles stand for one site allele.
    optional<AtomizeRefusal> cut(const vector<pair<int, int>>& clusters,
                                 vector<vcflib::Variant>& built, vector<bool>& one_to_one) const {
        for (const auto& [rb, re] : clusters) {
            ClusterSpelling spelling = alignment.spell_cluster((size_t)rb, (size_t)re);
            int64_t pos = record.unflattened_position + (int64_t)alignment.ref_offset(spelling.vb);
            if (needs_anchor_base(spelling)) {
                if (optional<AtomizeRefusal> refused = add_anchor_base(spelling)) {
                    return refused;
                }
                pos -= 1;
            }
            BlockAlleles block = number_block_alleles(spelling, alignment, record.trav_to_allele);
            if (block.alleles.size() < 2) {
                continue;
            }
            built.push_back(block_record(spelling, block, pos, built.size()));
            one_to_one.push_back(block.one_to_one());
        }
        return std::nullopt;
    }

private:
    /// Put the reference base before the cluster at the start of the reference and of each
    /// allele of `spelling`. Return why not, if there is no such base.
    optional<AtomizeRefusal> add_anchor_base(ClusterSpelling& spelling) const {
        if (spelling.vb <= 0) {
            // This also keeps us from reading the visit before the first, which a cluster at
            // the first visit has when the start node appears twice in the reference traversal.
            return AtomizeRefusal::NoAnchorBase;
        }
        string left = alignment.spell_reference(spelling.vb - 1, spelling.vb);
        if (left.empty()) {
            return AtomizeRefusal::AnchorWithoutSequence;
        }
        string base(1, left.back());
        spelling.ref_str = base + spelling.ref_str;
        for (size_t s = 0; s < spelling.slot_str.size(); ++s) {
            if (!spelling.slot_marker[s]) {
                spelling.slot_str[s] = base + spelling.slot_str[s];
            }
        }
        return std::nullopt;
    }

    /// Build block record number `index` of the site, at `pos`, and trim its alleles.
    vcflib::Variant block_record(const ClusterSpelling& spelling, const BlockAlleles& block,
                                 int64_t pos, size_t index) const {
        // Inherit INFO, FORMAT and every field of the site record.
        vcflib::Variant b_var = record.variant;
        b_var.position = pos;
        b_var.ref = block.alleles[0];
        b_var.alt.assign(block.alleles.begin() + 1, block.alleles.end());
        b_var.alleles = block.alleles;
        b_var.info["AT"].clear();
        b_var.info["AT"].resize(block.alleles.size());
        // The second field, the number of records, is known once every cluster is cut.
        b_var.info["SB"] = {std::to_string(index), ""};

        add_block_paths(b_var, spelling, block);

        vector<int> slot_site_allele(block.block_gt.size());
        for (size_t s = 0; s < block.block_gt.size(); ++s) {
            auto found = record.trav_to_allele.find(alignment.traversal_of(s));
            slot_site_allele[s] = found != record.trav_to_allele.end() ? found->second : -1;
        }
        const BlockGenotype genotype = phase_block_genotype(
            site_gt_text(record.variant, sample_name), block.block_gt, slot_site_allele);
        b_var.samples[sample_name]["GT"] = {genotype.gt};
        if (!genotype.keep_phase_set) {
            // PS labels a phase block, so it means nothing on an unphased genotype.
            b_var.samples[sample_name].erase("PS");
            std::erase(b_var.format, "PS");
        }

        copy_site_evidence(b_var, record.variant, sample_name, block, gl_layout);

        b_var.updateAlleleIndexes();
        vg::flatten_common_allele_ends(b_var, true, 0);
        vg::flatten_common_allele_ends(b_var, false, 0);
        return b_var;
    }

    /// Write AT for each allele of a block record, as the visits that spell it over the cluster.
    void add_block_paths(vcflib::Variant& b_var, const ClusterSpelling& spelling,
                         const BlockAlleles& block) const {
        add_allele_path_to_info(b_var, 0,
                                visits_of(graph, alignment.reference_visits(spelling.vb, spelling.ve)),
                                false, translation);
        for (size_t a = 1; a < block.alleles.size(); ++a) {
            // Spell the allele by the first slot that carries it.
            Traversal span;
            for (size_t s = 0; s < block.block_gt.size(); ++s) {
                if (!spelling.slot_marker[s] && block.block_gt[s] == (int)a) {
                    span = spelling.slot_span[s];
                    break;
                }
            }
            add_allele_path_to_info(b_var, (int)a, visits_of(graph, span), false, translation);
        }
    }

    const PathPositionHandleGraph& graph;
    /// The alignment the clusters are cut from.
    const SiteAlignment& alignment;
    /// The finished site record.
    const SiteRecord& record;
    /// The sample the records are written for.
    const string& sample_name;
    /// Names for node IDs, or null to write the IDs.
    const NodeTranslation* translation;
    /// The order of GL in the site record.
    GLLayout gl_layout;
};

/// Return why a site that cuts into one block record keeps its site record instead, if it does.
///
/// `one_to_one` says whether no two alleles of the block stand for one site allele.
optional<AtomizeRefusal> single_block_refusal(const SiteAlignment& alignment,
                                              const SiteRecord& record, bool one_to_one) {
    if (!alignment.site_says_more(record.trav_to_allele)) {
        return AtomizeRefusal::SameAsSiteRecord;
    }
    // A chain is genotyped from the first crossing of it by each allele, so where the reference
    // or a strand crosses a chain more than once, the record of the chain need not report the
    // crossing the block leaves to it.
    if (alignment.crosses_a_chain_twice()) {
        return AtomizeRefusal::ChainCrossedTwice;
    }
    // The block takes its genotype, likelihoods and phase from the site, through the site allele
    // each of its alleles stands for, so it cannot be written where two of its alleles stand for
    // one. Two routes can spell one site allele and still differ here, where a difference outside
    // a child chain is cancelled by one inside it.
    if (!one_to_one) {
        return AtomizeRefusal::StrandsDisagree;
    }
    return std::nullopt;
}

}

void BlockRecordWriter::report() const {
    // Report whenever block emission ran, even if it refused nothing.
    if (counters.sites.load() == 0) {
        return;
    }

    // The first number should be zero when every site is a snarl the manager knows, since a
    // reversed one resolves through its reversed bounds. A chain called as a site is a
    // constructed snarl the manager does not know.
    cerr << "[vg call] atomize: " << counters.site_unresolvable.load()
         << " sites where projection is inert because the snarl does not resolve, "
         << counters.site_reversed.load()
         << " resolved as the reversal flip_snarl produces" << endl;

    if (counters.child_inlined.load() > 0) {
        cerr << "[vg call] atomize: " << counters.child_inlined.load()
             << " child chains left without a line because a block ALT already spells them" << endl;
    }
    report_refusals();
    if (counters.split_sites.load() > 0) {
        cerr << "[vg call] atomize: " << counters.split_sites.load()
             << " sites written as their difference blocks rather than their site record, "
             << counters.split_lines.load()
             << " lines" << endl;
    }
}

void BlockRecordWriter::report_refusals() const {
    size_t total = 0;
    for (const std::atomic<size_t>& refused : counters.refuse) {
        total += refused.load();
    }
    if (total == 0) {
        return;
    }
    cerr << "[vg call] atomize: " << total << " sites declined block emission, so the site"
         << " record stands:";
    // List only the reasons that occurred.
    bool first = true;
    for (size_t i = 0; i < counters.refuse.size(); ++i) {
        size_t n = counters.refuse[i].load();
        if (n > 0) {
            cerr << (first ? " " : "; ") << n << " " << g_atomize_refuse_name[i];
            first = false;
        }
    }
    cerr << endl;
}

BlockRecordWriter::ChainInlineContext BlockRecordWriter::chain_inline_context(
    const HandleGraph& graph, const SiteChildren& children, const vector<Traversal>& travs,
    const vector<int>& genotype, int ref_trav_idx) const {
    ChainInlineContext ctx;
    // With one record per snarl, no chain is inside a block.
    if (!enabled || !nested) {
        return ctx;
    }
    if (ref_trav_idx < 0 || (size_t)ref_trav_idx >= travs.size() || genotype.empty()) {
        return ctx;
    }
    // A snarl whose projection has no symbols would answer that no child is reported.
    if (!children.known) {
        return ctx;
    }
    // A strand that carries the reference allele matches every reference step, including the
    // chain.
    if (std::ranges::find(genotype, ref_trav_idx) != genotype.end()) {
        return ctx;
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
    const nid_t chain_start = graph.get_id(chain.start);
    const nid_t chain_end = graph.get_id(chain.end);

    // A chain the reference does not cross is for the caller to handle.
    if (std::ranges::none_of(ctx.sref, [&](const SymbolicStep& step) {
            return is_chain_step(step, chain_start, chain_end);
        })) {
        return false;
    }

    bool any_crossing = false;
    for (const ChainInlineContext::Alt& alt : ctx.alts) {
        // Look for the crossings of the chain in the projection of this strand itself. A strand
        // that deletes the chain has none, so no block can report a copy of it, even when the
        // step of the reference for the chain lies inside a block.
        for (size_t j = 0; j < alt.salt.size(); ++j) {
            if (!is_chain_step(alt.salt[j], chain_start, chain_end)) {
                continue;
            }
            any_crossing = true;
            if (!in_a_block(alt.blocks, j)) {
                // The strand matches the reference here, so the record of the chain reports it.
                return false;
            }
        }
    }

    if (!any_crossing) {
        // No called strand crosses the chain, so no block ALT spells it.
        return false;
    }

    // Every crossing by every called strand lies inside a difference block, whose ALT spells the
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
    if (!children.known) {
        // The projection would be a bare node list here.
        ++counters.site_unresolvable;
        return;
    }
    if (children.reversed) {
        // Count it here, once per record, rather than where the projection is made, which can
        // happen many times per record.
        ++counters.site_reversed;
    }
}

int BlockRecordWriter::refuse(AtomizeRefusal reason) const {
    ++counters.refuse[(size_t)reason];
    return -1;
}

int BlockRecordWriter::write(const PathPositionHandleGraph& graph, const SiteChildren& children,
                             const vector<Traversal>& called_traversals,
                             const vector<int>& genotype, int ref_trav_idx,
                             const string& sample_name, const NodeTranslation* translation,
                             const SiteRecord& record, GLLayout gl_layout, bool genotype_snarls,
                             const function<bool(vcflib::Variant&, size_t)>& add_line,
                             const AlleleSequences* sequences) const {
    // Block emission being off is not a refusal, so it is not counted.
    if (!enabled || !nested || genotype_snarls) {
        return -1;
    }
    if (genotype.empty()) {
        return refuse(AtomizeRefusal::NoGenotype);
    }
    if (record.alleles_merged) {
        return refuse(AtomizeRefusal::AllelesMerged);
    }
    if (ref_trav_idx < 0 || (size_t)ref_trav_idx >= called_traversals.size()) {
        return refuse(AtomizeRefusal::NoReferenceTraversal);
    }
    if (!children.known) {
        // The projection would see no child chains here.
        return refuse(AtomizeRefusal::Unresolvable);
    }

    const SiteAlignment alignment(graph, children, called_traversals, genotype, ref_trav_idx,
                                  sequences);
    if (optional<AtomizeRefusal> refused = alignment.refusal()) {
        return refuse(*refused);
    }
    const vector<pair<int, int>> clusters = alignment.difference_clusters();
    if (clusters.empty()) {
        return refuse(AtomizeRefusal::NoDifferenceBlocks);
    }

    // Build every record before writing any, so that the decision to split is made on the
    // finished set.
    vector<vcflib::Variant> built;
    vector<bool> built_one_to_one;
    const BlockCutter cutter(graph, alignment, record, sample_name, translation, gl_layout);
    if (optional<AtomizeRefusal> refused = cutter.cut(clusters, built, built_one_to_one)) {
        return refuse(*refused);
    }
    if (built.empty()) {
        return refuse(AtomizeRefusal::ReferenceBasesOnly);
    }
    if (built.size() == 1) {
        if (optional<AtomizeRefusal> refused = single_block_refusal(alignment, record,
                                                                     built_one_to_one[0])) {
            return refuse(*refused);
        }
    }

    int added = 0;
    size_t block_index = 0;
    for (vcflib::Variant& b_var : built) {
        // Count only the records the site writes, not the clusters every strand spells as the
        // reference.
        b_var.info["SB"][1] = std::to_string(built.size());
        // VCF allows an ID on one record only, so each block record takes the ID of the site
        // with its index appended. The name of a snarl has no '_', so the ID of the site can be
        // read back from it (see block_site_name).
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

#ifndef VG_OFF_PANEL_HPP_INCLUDED
#define VG_OFF_PANEL_HPP_INCLUDED

/** \file
 * Off-panel detection, stage A: finding the substitutions and short indels that recur in the reads
 * but that no graph walk spells, while a window of reads is in memory. Each surviving candidate
 * later becomes an edit site (`edit_site_of`), genotyped, phased and written like any other site.
 */

#include <compare>
#include <cstdint>
#include <functional>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include <vg/vg.pb.h>

#include "../src/handle.hpp"

namespace vg {

using namespace std;

/// The shortest indel counted as large: matched loosely, since the reads rarely spell it alike.
constexpr uint32_t LARGE_INDEL_MIN_LENGTH = 50;

/// How the reads are matched to an indel of `LARGE_INDEL_MIN_LENGTH` bases or more. Long reads show
/// one as an indel inside the read, placed and sized a little differently by each read's errors;
/// short reads mostly as clips where the read crosses a deletion's ends. A read's indels of the
/// same kind near the breakpoint are summed, so that one the aligner split counts whole.
struct LargeIndelParams {
    /// The longest indel counted.
    uint32_t max_length = 100000;
    /// How far from the breakpoint, in bases, a read's indel or clip may lie and still match.
    int64_t slop = 10;
    /// How far, as a share of the length, a read's indel may differ in length and still match.
    double length_ratio = 0.2;
    /// A read is REF if it spans a breakpoint by at least this many aligned bases on each side,
    /// with no clip and no matching indel there.
    int64_t flank = 20;
    /// The shortest clip counted, and the bases of it matched to the graph to find the far end
    /// of the deletion it shows.
    uint32_t min_clip = 20;
    /// The longest deletion found from clips alone.
    uint32_t max_clip_deletion = 10000;
    /// A deletion found from clips alone, with no read holding it whole, needs this many
    /// fragments clipped at each of its ends.
    size_t min_fragments_per_end = 2;
};

struct EditCallBar;

/// The floors a candidate must pass to be kept, set low so that the call bar can be fitted from
/// the dump. Where `bar` is set, the candidates it merges are dumped as their merged candidate.
struct EditCountParams {
    /// Reads below this mapping quality are not counted.
    int min_mapq = 5;
    /// Bases below this quality are not counted, as the REF or the ALT. A read with no base
    /// qualities is counted whole.
    int min_base_quality = 20;
    /// The fewest fragments carrying the ALT.
    size_t min_fragments = 3;
    /// The smallest share of the fragments covering the base that carry the ALT.
    double min_fraction = 0.15;
    /// The longest indel counted, unless `large_indels`.
    uint32_t max_indel_length = 49;
    /// Whether to count indels of `LARGE_INDEL_MIN_LENGTH` bases or more, matched as `large` says.
    bool large_indels = false;
    LargeIndelParams large;
    /// Where not null, the candidates on one node that it admits and whose changed bases touch are
    /// merged into one, of kind 'M'.
    const EditCallBar* bar = nullptr;
};

/// A change the reads make to a walk of the graph: a substitution, the node's base `ref` read as
/// `alt` at one base of one node, or an insertion or deletion of any length. The offset and the
/// bases are along the node's forward strand. Fragments are counted once however many of their
/// reads cover the change, by read name, which mates share; a fragment whose reads disagree about
/// it is counted in `discordant_fragments` and nowhere else.
///
/// An indel is placed as far left along the walk as it stays the same change, and `node` and
/// `offset` give its anchor, the base before it, so that an indel the aligner put anywhere in its
/// repeat tract has one key. A read covering the anchor and the first base after the tract counts
/// as the ALT if its one indel there is this one, as the REF if it has none there, and as other if
/// it has another. A large indel (`is_large`) is matched loosely instead (`observe_large_indel`),
/// and its ALT is one read's: for an insertion, the medoid of the reads that carry it.
///
/// A merged candidate (kind 'M') replaces the bases `sequence` from `offset` on along `walk` with
/// each of its ALTs, which are what the reads spell there; see `count_edits`.
struct EditCandidate {
    nid_t node = 0;
    uint32_t offset = 0;
    char ref = 'N';
    /// For a substitution, the base read; 0 for an indel.
    char alt = 'N';
    /// 'S' for a substitution, 'D' for a deletion, 'I' for an insertion, 'M' for a merged change.
    char kind = 'S';
    /// For an indel, the bases it deletes or inserts, along the forward strand; for a merged
    /// change, the bases it replaces.
    string sequence;
    /// For a merged change, its ALTs, along the forward strand.
    vector<string> alts;
    /// For an indel, the nodes from the anchor's to that of the first base after the tract, each
    /// read forward: the walk its edit site takes.
    vector<nid_t> walk;
    /// For an indel, the bases from the anchor to the first base after the tract, the tract's
    /// length plus 1; for a merged change, the bases it replaces.
    uint32_t span = 0;
    /// For an indel, the length of the run of one base its tract starts with.
    uint32_t homopolymer = 0;
    uint32_t alt_fragments = 0;
    uint32_t ref_fragments = 0;
    /// Fragments with a third base.
    uint32_t other_fragments = 0;
    uint32_t discordant_fragments = 0;
    /// The mean mapping quality of the reads carrying the ALT.
    float mean_alt_mapq = 0.0f;

    /// The share of the fragments covering the base, discordant ones aside, that carry the ALT.
    double fraction() const;

    bool is_indel() const {
        return kind == 'D' || kind == 'I';
    }

    /// Whether it is an indel of `LARGE_INDEL_MIN_LENGTH` bases or more.
    bool is_large() const {
        return is_indel() && sequence.size() >= LARGE_INDEL_MIN_LENGTH;
    }

    /// Order candidates by the change they make: by node, offset, kind, ALT base and sequence.
    std::strong_ordering operator<=>(const EditCandidate& other) const {
        return std::tie(node, offset, kind, alt, sequence)
               <=> std::tie(other.node, other.offset, other.kind, other.alt, other.sequence);
    }

    /// Whether two candidates make the same change, whatever their counts.
    bool operator==(const EditCandidate& other) const {
        return (*this <=> other) == 0;
    }
};

/// The call bar: which candidates become edit sites, and which of those are staged once genotyped.
///
/// It applies after the floors of `EditCountParams`. The defaults were fitted on chr20 at about
/// 30x, for short reads and ONT reads alike, from a dump taken with the floors as the bar.
struct EditCallBar {
    /// The fewest fragments carrying the ALT.
    size_t min_alt_fragments = 7;
    /// The smallest ALT fraction.
    double min_fraction = 0.25;
    /// The lowest GQ at which a site called with the ALT is staged.
    double min_gq = 0.0;
    /// The highest depth ratio (DR) at which a site called with the ALT is staged, or 0 for no
    /// limit.
    ///
    /// DR is the number of reads at the site over the number its called genotype predicts. A
    /// collapsed repeat piles the reads of several copies onto one walk, and the edits that tell
    /// the copies apart then look like variants.
    double max_depth_ratio = 3.0;
    /// The fewest fragments carrying the ALT of an indel.
    size_t indel_min_alt_fragments = 7;
    /// The smallest ALT fraction of an indel.
    double indel_min_fraction = 0.25;
    /// No indel is called whose tract starts with a run of this many of one base or more; 0 for
    /// no limit, and 1 calls no indel.
    uint32_t indel_max_homopolymer = 0;
    /// The fewest fragments carrying the ALT of an indel of `LARGE_INDEL_MIN_LENGTH` bases or more.
    /// No homopolymer limit applies to such an indel.
    size_t large_min_alt_fragments = 7;
    /// The smallest ALT fraction of an indel of `LARGE_INDEL_MIN_LENGTH` bases or more.
    double large_min_fraction = 0.25;

    /// Whether `candidate` clears the bar for its kind.
    ///
    /// A merged change always does, since each of its ALTs cleared it already.
    bool admits(const EditCandidate& candidate) const {
        if (candidate.kind == 'M') {
            return true;
        }
        if (candidate.is_large()) {
            return candidate.alt_fragments >= large_min_alt_fragments
                   && candidate.fraction() >= large_min_fraction;
        }
        if (candidate.is_indel()) {
            return candidate.alt_fragments >= indel_min_alt_fragments
                   && candidate.fraction() >= indel_min_fraction
                   && (indel_max_homopolymer == 0
                       || candidate.homopolymer < indel_max_homopolymer);
        }
        return candidate.alt_fragments >= min_alt_fragments
               && candidate.fraction() >= min_fraction;
    }
};

/// What became of one candidate that passed the call bar: where it was placed, and its direct-pass
/// call.
struct EditOutcome {
    nid_t node = 0;
    uint32_t offset = 0;
    char alt = 'N';
    char kind = 'S';
    string sequence;
    /// "top" or "nested" where the site was genotyped, or why it was not.
    string placed;
    string site;
    vector<int> genotype;
    double gq = 0.0;
    /// The direct call's depth ratio, or -1 where it has none.
    double depth_ratio = -1.0;
    /// Whether it was staged, to be phased and written like any other site.
    bool staged = false;
};

/// A ledger of edit outcomes, which each OpenMP thread files into its own list, so that filing
/// needs no lock.
class EditOutcomes {
public:
    /// Make a ledger for thread numbers below `threads`.
    explicit EditOutcomes(size_t threads) : by_thread(threads) {}
    /// File an outcome from the calling OpenMP thread.
    void add(EditOutcome&& outcome);
    /// Every outcome, by node, offset and ALT.
    vector<EditOutcome> sorted() const;
private:
    vector<vector<EditOutcome>> by_thread;
};

/// Count the substitutions and indels that the reads `for_each_read` visits make, keyed on the
/// nodes `owns` accepts, and return those that pass `params`, sorted. A change that another walk
/// of the graph spells, as the other branch of a bubble does, is a panel allele the mapper did not
/// take, and is left out. An indel is keyed on its anchor once left-aligned, and every read that
/// covers the anchor is fetched with the anchor's node, so each window counts the indels anchored
/// on its own nodes in full, wherever its reads first placed them.
///
/// With `params.bar`, the candidates on one node that it admits, large indels aside, whose changed
/// bases touch or overlap are merged: a substitution changes its base, an indel its repeat tract.
/// The ALTs are the spellings the reads on the longest of their walks give between the bases on
/// either side, each admitted as a candidate of its own would be, and changing the length only
/// with a walk base on each side. They are not merged where they are all indels at one anchor, or
/// where one is in no ALT.
vector<EditCandidate> count_edits(
    const HandleGraph& graph,
    const function<void(const function<void(const Alignment&)>&)>& for_each_read,
    const function<bool(nid_t)>& owns, const EditCountParams& params);

/// An insertion or deletion along a walk, in walk coordinates, the bases counted from the walk's
/// first: `deleted` bases removed from `at` on, or `inserted` placed before the base at `at`.
struct WalkIndel {
    int64_t at = 0;
    uint32_t deleted = 0;
    string inserted;

    bool operator==(const WalkIndel& other) const = default;
    /// Its length.
    uint32_t length() const {
        return deleted > 0 ? deleted : (uint32_t)inserted.size();
    }
};

/// A read's alignment along a walk: the walk bases it covers, its indels in walk order, and where
/// its aligned bases lie in the read.
struct WalkAlignment {
    /// The walk bases [begin, end) from the read's first aligned base to its last.
    int64_t begin = 0;
    int64_t end = 0;
    vector<WalkIndel> indels;
    /// A run of aligned bases: its first walk base, its first read base, and its length.
    struct Block {
        int64_t walk;
        size_t read;
        size_t length;
    };
    vector<Block> blocks;
    /// The read bases inserted before the first aligned base and after the last: the read's
    /// clips, where these mappings begin and end the read.
    size_t clip_begin = 0;
    size_t clip_end = 0;

    /// The read offset of the base aligned to walk base `x`, or -1 where none is.
    int64_t read_offset_at(int64_t x) const;
};

/// The alignment along a walk of the consecutive mappings `mappings`, the i-th on a node read the
/// way the walk reads it and starting at walk base `node_start[i]`, the first mapping's first base
/// at read offset `read_offset`. Deletions that meet are stitched into one, across mappings too.
/// Insertions before the first aligned base or after the last, which are clips, are left out. The
/// indels are where the aligner put them; see `left_align`.
WalkAlignment walk_alignment(const vector<const Mapping*>& mappings,
                             const vector<int64_t>& node_start, size_t read_offset,
                             string_view sequence);

/// Move `indel` as far left as it stays the same change, keeping it at or after `floor`. `base(x)`
/// is the walk's base at x.
void left_align(WalkIndel& indel, int64_t floor, const function<char(int64_t)>& base);

/// `left_align` each of `aln`'s indels, keeping a base before each that the read covers and that
/// no earlier indel deletes, so that each has an anchor.
void left_align(WalkAlignment& aln, const function<char(int64_t)>& base);

/// The first walk base after the repeat tract of the left-aligned `indel`: the base after it once
/// moved as far right as it stays the same change, no further than `limit`, the end of the bases
/// known.
int64_t tract_end(const WalkIndel& indel, int64_t limit, const function<char(int64_t)>& base);

/// What a read says about an indel.
enum class IndelObservation { NONE, REF, ALT, OTHER };

/// What the left-aligned alignment `aln` says about the left-aligned `indel`, whose tract ends at
/// walk base `after`: NONE if the read does not cover the anchor (`indel.at - 1`) through `after`,
/// ALT if its one indel touching that stretch is `indel`, REF if it has none, and OTHER if it has
/// another. `touching`, if given, receives the indels touching the stretch.
IndelObservation observe_indel(const WalkAlignment& aln, const WalkIndel& indel, int64_t after,
                               vector<const WalkIndel*>* touching = nullptr);

/// What the left-aligned alignment `aln` says about the left-aligned indel `indel` of
/// `LARGE_INDEL_MIN_LENGTH` bases or more, whose tract ends at walk base `after`, matched as
/// `params` says. ALT if the read's indels of its kind within `slop` of the tract sum to its length
/// within `length_ratio`, or, for a deletion, if the read is clipped by `min_clip` bases or more
/// at either of its ends; REF if the read spans either end by `flank` bases each side with neither;
/// OTHER if it spans one with something else; and NONE if it spans neither. The read reaches from
/// walk base `reach_begin` to `reach_end`, which may lie past the bases `aln` covers where the
/// read carries on beyond the walk. `by_clip`, if given, is set to 1 for an ALT that was a clip at
/// the deletion's left end, 2 for one at its right end, and 0 otherwise.
IndelObservation observe_large_indel(const WalkAlignment& aln, const WalkIndel& indel,
                                     int64_t after, const LargeIndelParams& params,
                                     int64_t reach_begin, int64_t reach_end,
                                     int* by_clip = nullptr);

/// Visit the candidates in `sorted` that lie on the nodes `nodes`. `sorted` is sorted by node,
/// and `nodes` is sorted and free of duplicates.
void for_each_candidate_on(const vector<EditCandidate>& sorted, const vector<nid_t>& nodes,
                           const function<void(const EditCandidate&)>& iteratee);

}

#endif

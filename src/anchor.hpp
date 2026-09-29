#ifndef VG_ANCHOR_HPP_INCLUDED
#define VG_ANCHOR_HPP_INCLUDED

/** \file anchor.hpp
 *
 * Assembly anchors (--anchors-out): for each genotyped site, which reads support which of the
 * sample's strands, for pangenome-guided assembly.
 *
 * An anchor is a pin: a point between two adjacent bases of the graph, with no sequence of its
 * own, together with the reads that cross it and where each crosses it. Each site has two pins,
 * both reading in the site's direction:
 *
 *   - the S pin, just after the last base of the site's start boundary node;
 *   - the E pin, just before the first base of the site's end boundary node.
 *
 * A pin is at one end of a node, not at the node, so two snarls that share a boundary node in a
 * chain pin at opposite ends of it and never coincide, however short the node is.
 *
 * In a read, a pin is at a junction between two consecutive nodes of its walk. A junction can
 * hold only the S pin of the snarl starting at its left node and the E pin of the snarl ending at
 * its right node, and if it holds both, they belong to the same snarl. That happens for a read
 * that crosses straight from a snarl's start boundary to its end boundary, carrying an allele that
 * deletes the site's interior. Such a read is kept at the S pin and dropped at the E pin; see
 * `AnchorCounters::coincident`.
 */

#include <atomic>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "handle.hpp"
#include "site_read_source.hpp"
#include "vg/vg.pb.h"

namespace vg {

using namespace std;

/**
 * Where one read crosses one pin.
 *
 * `offset` is the 0-based index, in the read as sequenced, of the last base before the pin,
 * reading in the site's direction. So for `strand == 0` the pin lies just after `read[offset]`,
 * and for `strand == 1` just before it.
 *
 * vg represents a reverse-strand alignment by reverse node visits rather than by storing the read
 * reverse-complemented, so the read as sequenced is the frame of the read file and of every other
 * offset in vg. A consumer that wants the oriented read uses `len - 1 - offset`.
 */
struct AnchorPlacement {
    /// Negative when the read does not cross this pin.
    int64_t offset = -1;
    /// 0 when the read crosses the pin in the site's direction, 1 when against it.
    uint8_t strand = 0;
    bool placed() const { return offset >= 0; }
};

/// Why a read could not be pinned, and whether the invariant held. Shared across threads.
struct AnchorCounters {
    /// Pins resolved and checked against the graph's own base.
    atomic<size_t> verified{0};
    /// Pins whose read base did not match the graph's base, which indicates a coordinate or
    /// strand error; see `resolve_anchor_pin`. Should be zero.
    atomic<size_t> verify_failed{0};
    /// The read's alignment never visits the pin's node.
    atomic<size_t> no_visit{0};
    /// It visits it more than once, so which visit the pin belongs to is ambiguous.
    atomic<size_t> repeat_visit{0};
    /// The node base the pin is defined against has no read base aligned to it, because a deletion
    /// covers it or the alignment stops before the node's end. The pin is not moved to an earlier
    /// base, which on a 1 bp boundary node would be the neighbouring snarl's pin.
    atomic<size_t> unaligned_base{0};
    /// E pin only: the read has no base just upstream on its own walk, because it begins here or
    /// deletes the node upstream. The pin is not moved to a base beyond that node, which would be
    /// the neighbouring snarl's pin.
    atomic<size_t> no_neighbour{0};
    /// A read whose two pins at one snarl resolved to the same position; kept at S, dropped at E.
    atomic<size_t> coincident{0};
    /// Reads excluded because their best-fitting allele was not a called one.
    atomic<size_t> off_call{0};
    /// Snarls whose two boundary nodes are the same node, where only the S pin is emitted.
    atomic<size_t> degenerate_site{0};
    /// End anchors suppressed because every read on them is also on their slot's start anchor.
    atomic<size_t> single_pin{0};
    /// Reads sharing a name with another read at the same site, as paired mates do. A consumer that
    /// keys on the name would merge them.
    atomic<size_t> shared_name{0};

    /// Homozygous sites split into two slots by the reads' strand log-odds (--anchors-hom-split),
    /// and those left as one slot because the reads did not divide.
    atomic<size_t> hom_split{0};
    atomic<size_t> hom_unsplit{0};
    /// Reads at a split homozygous site whose strand log-odds are zero. Both slots carry the same
    /// allele, so nothing places such a read in either.
    atomic<size_t> hom_split_no_opinion{0};
    /// Of those, the ones placed by a coin flip derived from the read, since the two slots spell
    /// the same allele. The rest, reads seen in more than one phase set, are dropped.
    atomic<size_t> hom_split_coin{0};
    /// Read placements at a heterozygous site whose slot weights --anchors-phase-hets changed by the
    /// read's strand log-odds. Reads with no strand log-odds, or seen in more than one phase set,
    /// are not counted, since the option does not affect them.
    atomic<size_t> het_phase_tilted{0};
    /// Under --anchors-strict-hets, placements that the sign of the strand log-odds moved off the
    /// slot the allele match would have chosen.
    atomic<size_t> het_strict_moved{0};

    atomic<size_t> phase_checked{0};
    atomic<size_t> phase_agree{0};
    atomic<size_t> phase_confident{0};
    atomic<size_t> phase_confident_agree{0};
    atomic<size_t> phase_no_opinion{0};

    void report(ostream& out) const;
};

/**
 * Resolve where one read crosses one pin, and check the result against the graph.
 *
 * `site_backward` is the orientation in which the site's alleles visit `node_id`, from the
 * snarl's own start or end visit, since the pin reads in the site's direction. `exit_pin` selects
 * the S pin, at the junction leaving the node, whose reference base is the node's last base in
 * the site's direction, rather than the E pin, at the junction entering it, whose reference base is
 * the node's first.
 *
 * The check is made here, while the read sequence is available: the read base the pin is resolved
 * against, complemented if the read visits the node in reverse, must equal the graph's base. A
 * mismatch is counted in `counters.verify_failed`.
 */
AnchorPlacement resolve_anchor_pin(const SiteRead& read, const HandleGraph& graph,
                                   nid_t node_id, bool site_backward, bool exit_pin,
                                   AnchorCounters& counters);


/// One read's contribution to a site's anchors.
struct AnchorRead {
    string name;
    /// The MAPQ-derived mismapping probability, clamped as the genotype model clamps it.
    float mismap = 0.0f;
    AnchorPlacement start_pin;
    AnchorPlacement end_pin;
};

/**
 * What a site's anchors need, kept from the sweep until the record is rendered.
 *
 * No alignment is kept: each read's pins are resolved to `(strand, offset)` while its alignment
 * is available. Held on the site's `ReadLikelihoodCallInfo`.
 */
struct AnchorSiteEvidence {
    vector<AnchorRead> reads;
    /// rel(r, a), row major, reads x alleles, row-normalised into [0,1].
    vector<float> rel;
    size_t n_alleles = 0;
    /// The alleles' full lengths and the site's mean read length, for the allele-length weights of
    /// each read's confidence. Full lengths rather than the unique lengths the genotype model
    /// uses, since unique lengths depend on the pair of alleles, which is not known until the
    /// genotype is.
    vector<uint32_t> allele_length;
    float mean_read_length = 0.0f;
    bool length_weighted = true;
    nid_t start_node = 0;
    nid_t end_node = 0;

    float rel_at(size_t read, size_t allele) const {
        return rel[read * n_alleles + allele];
    }
    /// Bytes held, for the progress line.
    size_t bytes() const;
};

/// Thresholds and switches, from the --anchors-* options.
struct AnchorParams {
    bool enabled = false;
    /// Write anchors only at heterozygous sites (--anchors-het-only). Homozygous and haploid sites
    /// carry no haplotype information, but they connect the reads that cross them.
    bool het_only = false;

    /// Write a slot's end anchor only where it holds at least this many reads that the slot's start
    /// anchor does not (--anchors-end-new); 0 writes it always.
    ///
    /// Both pins of a site divide the reads the same way, so an end anchor adds only another
    /// point, and if every read at it is also at the start anchor, it adds no link between reads.
    /// Which reads reach each pin depends on depth and read length, so the anchors written depend
    /// on them too.
    size_t end_pin_min_new = 0;
    /// Write anchors only at sites with no child chains (--anchors-leaf-only).
    bool leaf_only = false;
    size_t min_reads = 2;
    double min_gqn = 0.0;
    double min_read_score = 0.0;
    /// Keep reads whose best-fitting allele, over all candidate alleles, was not called
    /// (--anchors-keep-off-call). Off by default, so that a read that fits neither called allele is
    /// not placed on one.
    bool keep_off_call = false;

    /// Where this run's counters live. Not owned; the likelihood calculator counts into the same
    /// ones. Null means do not count.
    AnchorCounters* counters = nullptr;

    /// Divide a homozygous site's reads between two slots by the sign of their strand log-odds
    /// (--anchors-hom-split).
    ///
    /// A homozygous site's alleles say nothing about which strand a read came from, so its reads
    /// can be divided only by the heterozygous sites they also cross. Off by default, since a
    /// wrongly divided site is worse for an assembler than an undivided one.
    bool hom_split = false;

    /// Minimum size of a read's tempered strand log-odds for it to count as confidently placed
    /// when deciding whether a homozygous site may be split (--split-min-q). 0.5 is a strand
    /// probability of about 62%.
    ///
    /// This decides whether the site is split, not which reads are kept: every read at a split
    /// site is placed, and its confidence is written in its row so that a consumer can filter.
    double phase_min = 0.5;

    /// Minimum number of confidently placed reads on each strand before a homozygous site may be
    /// split (--split-min-side). A site whose reads all point to one strand has not been divided.
    ///
    /// A lower value splits more sites, which joins longer runs of anchors that each name a
    /// strand, but a run that joins across a phasing error puts the sequence after it on the wrong
    /// strand. The count is of reads with a resolved pin, a superset of the reads written to the
    /// file, so it cannot be recomputed exactly from the file. It is a raw count, not scaled to
    /// depth.
    size_t phase_min_side = 10;

    /// At a heterozygous site, choose a read's slot using its strand log-odds as well as its
    /// allele match (on unless --no-anchors-phase-hets).
    ///
    /// The strand log-odds, computed leaving this site out, change the read's slot weights as
    /// `phase_aware_correction` does in re-genotyping. Without them each site places a read from
    /// its own alleles alone, so neighbouring sites can disagree about a read's strand. With them,
    /// the anchors agree with the phasing where the strand is confident, so they are no longer an
    /// independent check on it.
    bool phase_hets = true;

    /// At a heterozygous site, choose a read's slot by the sign of its strand log-odds alone,
    /// ignoring the allele match (--anchors-strict-hets), as a comparison for `phase_hets`. A read
    /// with no strand log-odds keeps the slot of its allele match, so both rules place the same
    /// reads.
    bool strict_hets = false;
};

/**
 * Collects anchors while records are rendered, and writes the TSV.
 *
 * Rows are ordered by node ID, each anchor row followed by its read rows. Read rows do not repeat
 * the anchor's key, so a read row means nothing without the anchor row above it; the header says
 * so.
 */
class AnchorWriter {
public:
    struct ReadRow {
        string name;
        int64_t offset = -1;
        float score = 0.0f;
        uint8_t strand = 0;
    };
    struct Anchor {
        nid_t node = 0;
        string snarl;
        /// Which strand of the settled genotype this anchor's reads are on, and which candidate
        /// traversal that strand carries.
        ///
        /// `slot` indexes the settled phased pair, so slot i is field i of the record's `GT` at
        /// the same site ID: slot 0 the left allele, slot 1 the right. A homozygote has one slot
        /// holding every read. A nested chain at ploidy 1 also has one slot, which is the parent's
        /// strand that carries the chain, the one `nested_strand` names; the VCF writes it as `a|.`
        /// or `.|a`, so its slot can be 1.
        ///
        /// `allele` indexes the site's candidate traversals, not the VCF's ALT alleles, which are
        /// chosen later.
        int slot = 0;
        int allele = -1;
        double gqn = -1.0;
        double explained = 1.0;
        /// The site's reliability: the mean confidence of the reads written for it, each read
        /// counted once, as `PhaseSite::reliability` is. Written although it could be derived
        /// from the read rows, which give confidences to one decimal only. -1 where the site wrote
        /// no read.
        double reliability = -1.0;
        vector<ReadRow> reads;
    };

    explicit AnchorWriter(size_t threads);

    /// Add an anchor to the queue of the calling OpenMP thread, found by `omp_get_thread_num()`.
    void add(Anchor&& anchor);

    size_t anchor_count() const;
    size_t read_row_count() const;

    /// Sort by (node, snarl, slot) and write. Returns false if the file could not be opened.
    bool write(const string& path, const string& graph_name, const string& sample,
               const string& reads_source, double mismap_min, const AnchorParams& params);

private:
    vector<vector<Anchor>> queues;
};

/// The allele-length weights of the two slots, for a settled slot-to-allele mapping. Shared with
/// read phasing, so that the anchors and the phasing split a read between the strands the same
/// way.
vector<double> site_slot_weights(const vector<uint32_t>& allele_length, size_t n_alleles,
                                 float mean_read_length, bool length_weighted,
                                 const vector<int>& slot_allele);

/// Turn one site's evidence and its settled genotype into anchors.
///
/// Each read goes to the slot of the called allele with the larger
///
///     x_i = (1 - e_r) * v_i * rel(r, a_i)
///
/// where v_i are the allele-length weights, and its confidence is
/// -10 log10(1 - max_i x_i / (sum_i x_i + e_r)). At a heterozygous site the confidence can be at
/// most the heterozygous score ceiling, which --mismap-min sets.
///
/// `genotype` is phase-ordered. `haploid_slot` is the strand, 0 or 1, that a one-allele
/// `genotype` sits on, and is ignored otherwise; the caller supplies it from the phasing.
void build_site_anchors(const AnchorSiteEvidence& evidence, const vector<int>& genotype,
                        const string& snarl_id, double gqn, double explained, int haploid_slot,
                        const AnchorParams& params, AnchorCounters& counters,
                        vector<AnchorWriter::Anchor>& out,
                        const vector<double>* read_strand = nullptr);

}

#endif

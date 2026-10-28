#ifndef VG_ANCHOR_HPP_INCLUDED
#define VG_ANCHOR_HPP_INCLUDED

/** \file anchor.hpp
 *
 * Assembly anchors: for each genotyped site, which reads support which of the sample's strands,
 * for pangenome-guided assembly.
 *
 * A *pin* is a point between two adjacent bases of the graph. Each site has two, both read in the
 * site's direction: the *start pin*, just after its start boundary node, and the *end pin*, just
 * before its end boundary node. A read crosses a pin where its walk passes between those two
 * bases, and the anchors record where, as an offset into the read.
 *
 * A site's reads are divided among its *slots*, one for each distinct allele of its chosen
 * genotype, in phase order. An *anchor* is one pin together with the reads of one slot that cross
 * it.
 */

#include <array>
#include <atomic>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "../src/handle.hpp"
#include "site_read_source.hpp"
#include "vg/vg.pb.h"

namespace vg {

namespace multipass {
class ReadStrandTable;
}

using namespace std;

/**
 * Where one read crosses one pin.
 *
 * `offset` is the 0-based index, in the read as sequenced, of the last read base before the pin,
 * reading in the site's direction. So for `direction == 0` the pin lies just after
 * `read[offset]`, and for `direction == 1` just before it.
 *
 * vg stores a reverse-strand alignment as reverse node visits, with the read as sequenced. So
 * offsets here are in the frame of the read file, as are other read offsets in vg. A consumer that
 * wants the offset in the read oriented along the site uses `len - 1 - offset`.
 */
struct AnchorPlacement {
    /// The read base before the pin, or a negative number when the read does not cross the pin.
    int64_t offset = -1;
    /// 0 when the read crosses the pin in the site's direction, 1 when against it.
    ///
    /// The anchor file writes it in its `strand` column.
    uint8_t direction = 0;
    /// A hash of where the alignment starts.
    ///
    /// It tells apart the alignments that share one read name, such as paired-end mates.
    uint32_t alignment = 0;
    /// Whether the read crosses the pin.
    bool placed() const { return offset >= 0; }
};

/// A tally of the pins resolved, by outcome. Shared across threads.
struct AnchorPinCounters {
    /// Pins resolved and checked against the graph's own base.
    atomic<size_t> verified{0};
    /// Pins whose read base did not match the graph's base. Should be zero.
    ///
    /// A mismatch means a coordinate or strand error in `resolve_anchor_pin`.
    atomic<size_t> verify_failed{0};
    /// Reads whose alignment never visits the pin's node.
    atomic<size_t> no_visit{0};
    /// Reads that visit the pin's node more than once, so that the pin's visit is ambiguous.
    atomic<size_t> repeat_visit{0};
    /// Reads with no base aligned to the node base the pin is defined against.
    ///
    /// A deletion covers that base, or the alignment stops before it.
    atomic<size_t> unaligned_base{0};
    /// Reads at an end pin with no base just before the pin's node on their own walk.
    ///
    /// The read begins at the pin's node, or deletes the node before it.
    atomic<size_t> no_neighbour{0};

    /// Write the tally to a log.
    void report(ostream& out) const;
};

/**
 * A tally of how often a read's strand log-odds put it on the slot its allele chose.
 *
 * At a heterozygous site the alleles show which strand each read is on, so this measures how well
 * the strand log-odds would divide a homozygous site under `AnchorParams::hom_split`. Shared across
 * threads.
 */
struct AnchorPhaseCheckCounters {
    /// Reads with strand log-odds, checked.
    atomic<size_t> checked{0};
    /// Checked reads whose strand log-odds favour their slot.
    atomic<size_t> agree{0};
    /// Checked reads whose strand log-odds are at least `AnchorParams::phase_min` in size.
    atomic<size_t> confident{0};
    /// Confident reads whose strand log-odds favour their slot.
    atomic<size_t> confident_agree{0};
    /// Reads with no strand log-odds, or a strand from another phase set, which were not checked.
    atomic<size_t> no_opinion{0};

    /**
     * Count one read placed in `slot` by its allele, whose strand log-odds are `lo`.
     *
     * @param phase_min The size of strand log-odds that counts as confident.
     */
    void count(double lo, int slot, double phase_min);

    /// Write the tally to a log, if any read was checked.
    void report(ostream& out) const;
};

/// A tally of what happened to reads and sites while anchors were made. Shared across threads.
struct AnchorCounters {
    /// How the pins were resolved.
    AnchorPinCounters pins;

    /// Reads whose two pins at one site resolved to the same read position. Such a read is kept at
    /// the start pin only.
    ///
    /// A read that crosses straight from the start boundary node to the end boundary node, carrying
    /// an allele that deletes the interior, has both pins at one position.
    atomic<size_t> coincident{0};
    /// Reads left out because their best-fitting allele was not a called one.
    atomic<size_t> off_call{0};
    /// Sites whose two pins are one point. They get a start pin only.
    ///
    /// A snarl whose two boundary nodes are the same node is one.
    atomic<size_t> degenerate_site{0};
    /// End anchors left out because too few of their reads were missing from their slot's start
    /// anchor.
    atomic<size_t> single_pin{0};
    /// Reads sharing a name with another read at the same site, as paired mates do.
    ///
    /// The file names reads through one table, so such reads share one index there.
    atomic<size_t> shared_name{0};

    /// Homozygous sites split into two slots by the strand log-odds of their reads.
    atomic<size_t> hom_split{0};
    /// Homozygous sites considered for a split and left as one slot, because too few reads were
    /// confidently placed on one of the strands.
    atomic<size_t> hom_unsplit{0};
    /// Reads at a split homozygous site whose strand log-odds name no strand of its phase set.
    ///
    /// Those whose strand is from another phase set, or from more than one, are dropped.
    atomic<size_t> hom_split_no_opinion{0};
    /// Reads at a split homozygous site that have no strand log-odds at all.
    ///
    /// Both slots spell the same allele, so such a read is placed by a coin flip derived from its
    /// name, which puts it in the same slot at every site.
    atomic<size_t> hom_split_coin{0};
    /// Read placements at a heterozygous site whose slot choice used the strand log-odds of the read.
    ///
    /// Reads with no strand log-odds, or with a strand from another phase set, are left out.
    atomic<size_t> het_phase_tilted{0};
    /// Read placements that `AnchorParams::strict_hets` moved off the slot their allele match
    /// chose.
    atomic<size_t> het_strict_moved{0};

    /// How often the strand log-odds agree with the alleles at heterozygous sites.
    AnchorPhaseCheckCounters phase_check;

    /// Write the tally to a log.
    void report(ostream& out) const;
};

/**
 * Resolve where one read crosses one pin, and check the result against the graph.
 *
 * The check is made here because the read sequence is available only here. The read base the pin
 * is resolved against, complemented if the read visits the node in reverse, must equal the graph's
 * base. A mismatch is counted in `counters.pins.verify_failed`, and the read is not placed.
 *
 * @param site_backward The orientation in which the site's alleles visit `node_id`. The pin reads
 * in the site's direction.
 * @param exit_pin True for the start pin, where the walk leaves the node. Its reference base is the
 * last base of the node in the site's direction. False for the end pin, where the walk enters the
 * node. Its reference base is the first base of the node.
 * @param pin For a pin inside the node, the number of node bases before the pin, along the node's
 * forward strand, or -1 for a pin at a node end. An exit pin inside the node is resolved against
 * the node base just before it in the site's direction. An entry pin inside the node is resolved
 * against the node base just after it, and its offset is the read base before that one, which can
 * be an inserted base but not a soft-clipped one.
 */
AnchorPlacement resolve_anchor_pin(const SiteRead& read, const HandleGraph& graph,
                                   nid_t node_id, bool site_backward, bool exit_pin,
                                   AnchorCounters& counters, int64_t pin = -1);


/**
 * A table of read names, which gives each distinct name a 32-bit index.
 *
 * A site's anchor evidence and the anchor rows refer to a read by its index rather than keeping a
 * copy of its name, since on a whole genome a read is placed at hundreds of sites. Names can be
 * added from several threads at once, and looked up without a lock. A name, once added, never
 * moves.
 */
class ReadNameTable {
public:
    /// Return the index of a name, adding the name if it is new.
    uint32_t intern(string_view name);

    /// Return the name at an index this table gave out.
    string_view name(uint32_t index) const;

private:
    /// The number of shards, among which names are spread by hash.
    ///
    /// Each shard has its own lock, so that threads adding names rarely wait for each other. An
    /// index is the position of the name in its shard times SHARDS, plus the shard.
    static constexpr size_t SHARDS = 64;
    /// The number of names in one chunk of a shard.
    static constexpr size_t CHUNK = 1 << 16;
    /// The number of chunks a shard can hold.
    static constexpr size_t MAX_CHUNKS = 1 << 12;
    /**
     * One lock's share of the names.
     *
     * A shard stores its names in chunks that are never moved or freed until the shard is
     * destroyed, found through a fixed array of pointers. So a lookup needs no lock while another
     * thread adds a name.
     */
    struct Shard {
        /// Held while a name is added.
        std::mutex mutex;
        /// The index of each name in the shard, keyed by views of the names in `chunks`.
        unordered_map<string_view, uint32_t> index;
        /// The chunks, each an array of CHUNK names, or null where not yet allocated. Owned.
        std::array<std::atomic<string*>, MAX_CHUNKS> chunks{};
        /// The number of names in the shard.
        size_t count = 0;
        ~Shard();
    };
    std::array<Shard, SHARDS> shards;
};

/// Return the read names of the run.
ReadNameTable& read_names();

/// One read's contribution to a site's anchors.
struct AnchorRead {
    /// The read, as its index in read_names().
    uint32_t read = 0;
    /// The MAPQ-derived mismapping probability, clamped between the `min_mismap_prob` and
    /// `max_mismap_prob` of `AlleleLikelihoodParams`.
    float mismap = 0.0f;
    /// Where the read crosses the start pin of the site.
    AnchorPlacement start_pin;
    /// Where the read crosses the end pin of the site.
    AnchorPlacement end_pin;
};

/**
 * What a site's anchors need to know about its reads, kept until the record is rendered.
 *
 * Each read's pins are resolved while its alignment is available, so no alignment is kept. Held on
 * the site's `ReadLikelihoodCallInfo`.
 */
struct AnchorSiteEvidence {
    /// The site's reads.
    vector<AnchorRead> reads;
    /// rel(r, a), row major, reads x alleles, each row normalised to sum to 1.
    vector<float> rel;
    /// The number of candidate alleles, which is the row length of `rel`.
    size_t n_alleles = 0;
    /// The full length of each candidate allele, for the allele-length weights of each read's
    /// confidence.
    ///
    /// These are full lengths rather than the unique lengths the genotype model uses, since unique
    /// lengths depend on the pair of alleles called. They are kept before the genotype is known,
    /// and the linkage model can still change the pair.
    vector<uint32_t> allele_length;
    /// The mean read length R that the site's likelihood matrix used for its length weights.
    ///
    /// See AlleleReadLikelihoods::set_length_weights.
    float mean_read_length = 0.0f;
    /// Whether the allele-length weights depend on the allele lengths, or are flat.
    bool length_weighted = true;

    /// The node the start pin lies on.
    nid_t start_node = 0;
    /// The node the end pin lies on.
    nid_t end_node = 0;
    /// Where the start pin lies in its node: the number of node bases before it, along the node's
    /// forward strand.
    ///
    /// For a snarl the start pin follows its node in the site's direction.
    uint32_t start_pin = 0;
    /// Where the end pin lies in its node: the number of node bases before it, along the node's
    /// forward strand.
    ///
    /// For a snarl the end pin precedes its node in the site's direction.
    uint32_t end_pin = 0;
    /// Whether the site lies inside nodes, as an off-panel site does. Its pins then bracket the
    /// bases its alleles change.
    bool inside_node = false;

    /// Return rel(read, allele).
    float rel_at(size_t read, size_t allele) const {
        return rel[read * n_alleles + allele];
    }
    /// Return the bytes held, for the progress line.
    size_t bytes() const;
};

/// Thresholds and switches for the assembly anchors.
struct AnchorParams {
    /// Whether anchors are made at all.
    bool enabled = false;
    /// Write anchors only at diploid heterozygous sites.
    ///
    /// Homozygous sites and haploid ones, nested haploid chains included, are left out, though
    /// they connect the reads that cross them.
    bool het_only = false;

    /// Write a slot's end anchor only where it holds at least this many reads that the slot's start
    /// anchor does not. 0 writes it always.
    ///
    /// Both pins of a site divide the reads the same way, so an end anchor adds only another
    /// point. If every read at it is also at the start anchor, it adds no link between reads.
    /// Which reads reach each pin depends on depth and read length, so the anchors written depend
    /// on them too.
    size_t end_pin_min_new = 0;
    /// Write anchors only at sites with no child chains.
    bool leaf_only = false;
    /// The fewest reads an anchor must hold to be written.
    size_t min_reads = 2;
    /// The smallest GQN a site must have for its anchors to be written, or 0 for no limit.
    double min_gqn = 0.0;
    /// The smallest confidence a read must have to be placed in an anchor.
    double min_read_score = 0.0;
    /// Whether to keep reads whose best-fitting allele, over all candidate alleles, was not
    /// called.
    ///
    /// Off by default, so that a read that fits neither called allele is not placed on one.
    bool keep_off_call = false;

    /// Where this run's counters live, or null to not count.
    ///
    /// Not owned. The likelihood calculator counts into the same ones.
    AnchorCounters* counters = nullptr;

    /// Whether to divide a homozygous site's reads between two slots by the sign of their strand
    /// log-odds.
    ///
    /// A homozygous site's alleles say nothing about which strand a read came from, so its reads
    /// can be divided only by the heterozygous sites they also cross. Off by default, since a
    /// wrongly divided site is worse for an assembler than an undivided one.
    bool hom_split = false;

    /// Minimum size of a read's tempered strand log-odds for it to count as confidently placed,
    /// when deciding whether a homozygous site may be split.
    ///
    /// 0.5 is a strand probability of about 62%. This decides whether the site is split, not
    /// which reads are kept. At a split site every read is placed, except one whose strand is from
    /// another phase set, and its confidence is written in its row so that a consumer can filter.
    double phase_min = 0.5;

    /// Minimum number of confidently placed reads on each strand before a homozygous site may be
    /// split.
    ///
    /// A lower value splits more sites, which joins longer runs of anchors that each name a
    /// strand. But a run that joins across a phasing error puts the sequence after it on the wrong
    /// strand. The count is of reads with a resolved pin, a superset of the reads written to the
    /// file, so it cannot be recomputed exactly from the file. It is a raw count, not scaled to
    /// depth.
    size_t phase_min_side = 10;

    /// Whether to choose a read's slot at a heterozygous site using its strand log-odds as well as
    /// its allele match. On by default.
    ///
    /// The strand log-odds, computed leaving this site out, multiply the weight of the slot they
    /// disfavour by exp(-|log-odds|), as `phase_aware_correction` does in re-genotyping. Without
    /// them each site places a read from its own alleles alone, so neighbouring sites can disagree
    /// about a read's strand. With them, the anchors agree with the phasing where the strand is
    /// confident, so they are no longer an independent check on it.
    bool phase_hets = true;

    /// Whether to choose a read's slot at a heterozygous site by the sign of its strand log-odds
    /// alone, ignoring the allele match.
    ///
    /// This is a comparison for `phase_hets`. A read with no strand log-odds keeps the slot of its
    /// allele match, so both rules place the same reads.
    bool strict_hets = false;
};

/**
 * A collector of anchors, which writes them to the anchor TSV once all are in.
 *
 * Rows are ordered by node ID, each anchor row followed by its read rows. Read rows do not repeat
 * the key of their anchor, so a read row means nothing without the anchor row above it. The file
 * header says so.
 */
class AnchorWriter {
public:
    /// One read placed in an anchor, as the file writes it.
    struct ReadRow {
        /// The read, as its index in read_names().
        uint32_t read = 0;
        /// Which alignment of the read it is; see `AnchorPlacement::alignment`.
        uint32_t alignment = 0;
        /// The read base before the pin; see `AnchorPlacement::offset`.
        int64_t offset = -1;
        /// The read's confidence in its slot, phred-scaled.
        float score = 0.0f;
        /// The direction the read crosses the pin in; see `AnchorPlacement::direction`.
        uint8_t direction = 0;
    };
    /// One pin of one site, with the reads of one slot that cross it.
    struct Anchor {
        /// The node the pin lies on.
        nid_t node = 0;
        /// Where the pin lies in `node`: the number of node bases before it, along the node's
        /// forward strand.
        uint32_t pin = 0;
        /// The site's ID in the VCF.
        string snarl;
        /// Which slot of the chosen phased pair this anchor holds.
        ///
        /// Slot i is field i of the record's `GT` at the same site ID, so slot 0 is the left allele
        /// and slot 1 the right. A homozygote has one slot, 0, holding every read, unless
        /// `AnchorParams::hom_split` divides it. A nested chain at ploidy 1 also has one slot, the
        /// parent's strand that carries it, which the VCF writes as `a|.` or `.|a`, so its slot can
        /// be 1.
        int slot = 0;
        /// The candidate traversal the slot carries, as an index into the site's candidate
        /// traversals. This is not a VCF allele number.
        int allele = -1;
        /// The site's GQN, or NaN where it has none.
        double gqn = -1.0;
        /// The share of the site's reads that its genotype explains.
        double explained = 1.0;
        /// The site's reliability: the mean confidence of the reads written for it, each read
        /// counted once, or -1 where the site wrote no read.
        ///
        /// The read rows give confidences to one decimal only, so it is written here as well.
        double reliability = -1.0;
        /// Whether the anchor belongs to an off-panel site, which lies inside nodes.
        bool edit = false;
        /// The reads placed in the anchor.
        vector<ReadRow> reads;
    };
    /// What `write` wrote and left out.
    struct Tally {
        /// The anchors written.
        size_t anchors_written = 0;
        /// The read rows written.
        size_t placements_written = 0;
        /// The read rows of off-panel anchors left out, because another anchor holds their point.
        size_t placements_left_out = 0;
        /// The off-panel anchors left out, because too few of their reads were left.
        size_t anchors_left_out = 0;
    };

    /// Make a writer that collects anchors from up to `threads` OpenMP threads at once.
    explicit AnchorWriter(size_t threads);

    /// Add an anchor to the queue of the calling OpenMP thread, found by `omp_get_thread_num()`.
    void add(Anchor&& anchor);

    /// Return the number of anchors collected and not yet written.
    size_t anchor_count() const;

    /**
     * Sort the collected anchors by (node, snarl, slot) and write them to a file.
     *
     * No point of a read alignment is written in two anchors. An off-panel anchor leaves out a
     * read whose point a snarl's anchor holds, or an earlier off-panel anchor in the file. It is
     * then left out itself if it holds fewer than `params.min_reads` reads.
     *
     * @returns false if the file could not be opened.
     */
    bool write(const string& path, const string& graph_name, const string& sample,
               const string& reads_source, double mismap_min, const AnchorParams& params);

    /// Return what the last `write` wrote and left out.
    const Tally& tally() const { return written; }

private:
    /// Remove and return every collected anchor, in the order the file writes them.
    vector<Anchor> take_sorted();

    /// The anchors collected, one queue for each OpenMP thread.
    vector<vector<Anchor>> queues;
    /// What the last `write` wrote and left out.
    Tally written;
};

/**
 * A record of the phased call at a site and of the confidence that weights its anchors.
 */
struct AnchorCall {
    /// The called alleles, in phase order.
    vector<int> genotype;
    /// The strand, 0 or 1, that a one-allele `genotype` sits on, from the phasing.
    ///
    /// It is ignored for other genotypes.
    int haploid_slot = 0;
    /// The value for the gqn column, or NaN where the site had no gap to normalise.
    double gqn = 0.0;
    /// The share of the site's reads the genotype explains, for the explained column.
    double explained = 1.0;
};

/**
 * Turn one site's evidence and its chosen genotype into anchors.
 *
 * Each read goes to the slot of the called allele with the larger
 *
 *     x_i = (1 - e_r) * v_i * rel(r, a_i)
 *
 * where v_i are the allele-length weights. Its confidence is
 * -10 log10(1 - max_i x_i / (sum_i x_i + e_r)). The e_r in the denominator bounds the
 * confidence. For two alleles of equal length the bound is the heterozygous score ceiling,
 * which the floor on e_r sets.
 *
 * `read_strand`, indexed like `evidence.reads`, holds each read's tempered strand log-odds with
 * this site left out. A positive value favours slot 0, 0 means the read has none, and NaN means
 * its strand is from another phase set, or from more than one (see `read_strand_usable`). It
 * changes the slot choice in three cases:
 *
 * - At a heterozygous site under `params.phase_hets` (the default), the disfavoured slot's v_i
 *   is multiplied by exp(-|log-odds|), as `phase_aware_correction` weights it.
 * - Under `params.strict_hets`, the sign alone chooses a heterozygous site's slot.
 * - Under `params.hom_split`, a diploid homozygous site with enough confidently placed reads on
 *   each strand gets two slots, and the sign divides its reads between them. A read with no
 *   log-odds goes by a hash of its name, and a read whose log-odds are NaN is dropped.
 *
 * A read's confidence is computed without the strand log-odds in every case.
 *
 * @pre `params.counters` is not null.
 *
 * @param out Where the anchors are appended.
 */
void build_site_anchors(const AnchorSiteEvidence& evidence, const AnchorCall& call,
                        string_view snarl_id, const AnchorParams& params,
                        vector<AnchorWriter::Anchor>& out,
                        const vector<double>* read_strand = nullptr);

/**
 * A collector that turns each chosen site into anchors as the sites are rendered, and writes the
 * anchor file at the end of the run.
 *
 * It collects nothing until `configure` gives it a path. It is not copyable, since its parameters
 * can point at its own counters.
 */
class AnchorCollector {
public:
    AnchorCollector() = default;
    AnchorCollector(const AnchorCollector&) = delete;
    AnchorCollector& operator=(const AnchorCollector&) = delete;

    /**
     * Set up to write anchors to a file. Call before any thread collects.
     *
     * Where `anchor_params` has no counters, as in a unit test, this collector counts into its
     * own.
     *
     * @param anchor_params Which sites and reads qualify.
     * @param graph Names the graph in the file header.
     * @param reads Names the reads in the file header.
     * @param mismap_floor The floor on the mismapping probability, which bounds a read's anchor
     * confidence, for the file header.
     */
    void configure(string_view anchor_path, const AnchorParams& anchor_params, string_view graph,
                   string_view reads, double mismap_floor);

    /// Return whether anchors are being collected.
    bool is_enabled() const { return !path.empty(); }

    /// Return whether `collect` needs each site's leaf status, which has a cost to find.
    bool wants_leaf_test() const { return !path.empty() && params.leaf_only; }

    /**
     * Turn one chosen site into anchors.
     *
     * Safe to call from many threads, but not from inside a nested parallel region.
     *
     * @param evidence The reads' evidence the genotyper kept for the site.
     * @param call The phased call at the site. A NaN `gqn` is written as `.`.
     * @param site_name Names the site in the file.
     * @param is_leaf Whether the site has no child sites. Used only when `wants_leaf_test`.
     * @param strands Gives each read's strand log-odds from the sites other than `record_key`.
     */
    void collect(const AnchorSiteEvidence& evidence, const AnchorCall& call,
                 const string& site_name, bool is_leaf, const multipass::ReadStrandTable& strands,
                 size_t record_key);

    /// Write the anchor file, naming `sample_name`, and report the counters. Do nothing unless
    /// anchors are being collected.
    void write(const string& sample_name);

private:
    /// The anchor file, or empty when anchors are not collected.
    string path;
    /// Which sites and reads qualify. Its counters are never null after `configure`.
    AnchorParams params;
    /// The counters, when `configure` was given none.
    AnchorCounters owned_counters;
    /// Names the graph in the file header.
    string graph_name;
    /// Names the reads in the file header.
    string reads_source;
    /// Collects the anchors and writes the file. Null until `configure` is given a path.
    unique_ptr<AnchorWriter> writer;
    /// The floor on the mismapping probability, for the file header.
    double mismap_min = 0.0;
};

}

#endif

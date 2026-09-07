#include "anchor.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <unordered_map>

#include <omp.h>

#include "path.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Counters
////////////////////////////////////////////////////////////////////////////////

AnchorCounters& anchor_counters() {
    static AnchorCounters counters;
    return counters;
}

void AnchorCounters::report(ostream& out) const {
    out << "[vg call] anchors: " << verified.load() << " pins verified against the graph, "
        << verify_failed.load() << " failed" << endl;
    out << "[vg call] anchor pins refused: " << no_visit.load()
        << " read does not visit the boundary node, " << repeat_visit.load()
        << " visits it more than once, " << unaligned_base.load()
        << " no read base aligned to the pin's node base, " << no_neighbour.load()
        << " nothing upstream of an entry pin" << endl;
    out << "[vg call] anchor reads: " << coincident.load()
        << " with both pins at one position (kept at S, dropped at E), " << off_call.load()
        << " fitting no called allele";
    if (shared_name.load() > 0) {
        out << ", " << shared_name.load()
            << " sharing a read name with another read at the same site (paired mates do; the name "
               "then does not identify an alignment)";
    }
    if (degenerate_site.load() > 0) {
        out << ", " << degenerate_site.load() << " sites whose two boundaries are one node";
    }
    if (single_pin.load() > 0) {
        out << ", " << single_pin.load()
            << " sites given only their start pin because the end pin held no reads of its own";
    }
    out << endl;
}

////////////////////////////////////////////////////////////////////////////////
// Pin resolution
////////////////////////////////////////////////////////////////////////////////

/// The first and last read base of one mapping that consumes a node base, or -1 for a mapping that
/// consumes none -- a node deleted outright in this read.
///
/// Per MAPPING rather than as one merged run over the read, because the entry pin has to step
/// exactly one node along the read's walk, and a merged run cannot tell "the previous node's last
/// base" from "the last base before a stretch of deleted nodes".
struct MappingExtent {
    int64_t first_consuming = -1;
    int64_t last_consuming = -1;
};

static MappingExtent extent_of(const Mapping& m, size_t read_start) {
    MappingExtent out;
    size_t read_pos = read_start;
    for (int64_t j = 0; j < m.edit_size(); ++j) {
        const Edit& e = m.edit(j);
        if (e.from_length() > 0 && e.to_length() > 0 && e.from_length() == e.to_length()) {
            if (out.first_consuming < 0) {
                out.first_consuming = (int64_t)read_pos;
            }
            out.last_consuming = (int64_t)(read_pos + (size_t)e.to_length() - 1);
        }
        read_pos += (size_t)e.to_length();
    }
    return out;
}

AnchorPlacement resolve_anchor_pin(const SiteRead& read, const HandleGraph& graph,
                                   nid_t node_id, bool site_backward, bool exit_pin,
                                   AnchorCounters& counters) {
    AnchorPlacement out;
    const Alignment& aln = *read.aln;
    const Path& path = aln.path();

    // Find the pin node's mapping and where in the read it starts. The pin node is one of the
    // site's boundaries, so every visit to it is listed in the site's index when there is one --
    // which matters, since the count is what refuses an ambiguous repeat visit.
    int64_t hit = -1;
    size_t hit_read_start = 0;
    size_t visits = 0;
    // Read start of the mapping before the pin's, needed only by the entry pin. Tracked here rather
    // than recovered later because the unindexed walk passes it once and cannot go back.
    size_t before_read_start = 0;
    if (read.indexed()) {
        for (size_t k = 0; k < read.mapping_count; ++k) {
            int64_t i = (int64_t)read.mappings[k];
            if (path.mapping(i).position().node_id() != node_id) {
                continue;
            }
            ++visits;
            hit = i;
            hit_read_start = read.read_offsets[i];
        }
        if (hit > 0) {
            before_read_start = read.read_offsets[hit - 1];
        }
    } else {
        size_t read_pos = 0;
        for (int64_t i = 0; i < path.mapping_size(); ++i) {
            const Mapping& m = path.mapping(i);
            if (m.position().node_id() == node_id) {
                ++visits;
                hit = i;
                hit_read_start = read_pos;
            }
            read_pos += (size_t)mapping_to_length(m);
        }
        if (hit > 0) {
            before_read_start = hit_read_start
                                - (size_t)mapping_to_length(path.mapping(hit - 1));
        }
    }

    if (visits == 0) {
        // Expected rather than exceptional: a read that starts inside the site never visits the
        // start boundary, and one that ends inside it never visits the end boundary.
        ++counters.no_visit;
        return out;
    }
    if (visits > 1) {
        // Which visit the pin belongs to is genuinely ambiguous, so it is refused rather than
        // guessed at.
        ++counters.repeat_visit;
        return out;
    }

    const Mapping& m = path.mapping(hit);
    const bool is_rev = m.position().is_reverse();
    handle_t handle = graph.get_handle(node_id, false);
    const size_t node_len = graph.get_length(handle);
    if (node_len == 0) {
        ++counters.unaligned_base;
        return out;
    }

    // The node base the pin is defined against, in node-forward coordinates.
    //   exit  (S) pin: the node's LAST base in the site's direction; the pin follows it.
    //   entry (E) pin: the node's FIRST base in the site's direction; the pin precedes it.
    const size_t q_fwd = exit_pin ? (site_backward ? 0 : node_len - 1)
                                  : (site_backward ? node_len - 1 : 0);
    // The same base in the orientation the READ visits the node in, which is what the mapping's
    // offsets are measured in.
    const size_t q_vis = is_rev ? (node_len - 1 - q_fwd) : q_fwd;

    // Walk the mapping's edits for the read base aligned to it.
    size_t node_pos = (size_t)m.position().offset();
    size_t rp = hit_read_start;
    int64_t base_index = -1;
    bool exact_edit = false;
    char edit_base = 0;
    for (int64_t j = 0; j < m.edit_size(); ++j) {
        const Edit& e = m.edit(j);
        const size_t from = (size_t)e.from_length();
        const size_t to = (size_t)e.to_length();
        if (from > 0 && q_vis >= node_pos && q_vis < node_pos + from) {
            if (from == to) {
                base_index = (int64_t)(rp + (q_vis - node_pos));
                if (e.sequence().empty()) {
                    exact_edit = true;
                } else if (e.sequence().size() == to) {
                    edit_base = e.sequence()[q_vis - node_pos];
                }
            }
            // Otherwise the node base is deleted (to == 0) or sits inside a replacement of unequal
            // length, and no read base corresponds to it.
            break;
        }
        node_pos += from;
        rp += to;
    }
    if (base_index < 0 || (size_t)base_index >= aln.sequence().size()) {
        ++counters.unaligned_base;
        return out;
    }

    // The invariant, checked here because the read is gone by the time anchors are written.
    //
    // On a match edit the read base must be the graph's base, complemented when the read visits the
    // node in reverse -- one comparison that catches every off-by-one AND every strand inversion. On
    // a mismatch edit the base legitimately differs from the graph, so the edit's own recorded
    // sequence is the reference instead, which still catches the offset arithmetic.
    const char read_base = (char)toupper(aln.sequence()[base_index]);
    if (exact_edit) {
        const char node_base = graph.get_base(handle, q_fwd);
        const char want = (char)toupper(is_rev ? reverse_complement(node_base) : node_base);
        if (read_base != want) {
            ++counters.verify_failed;
            return out;
        }
        ++counters.verified;
    } else if (edit_base != 0) {
        if (read_base != (char)toupper(edit_base)) {
            ++counters.verify_failed;
            return out;
        }
        ++counters.verified;
    }

    // Does the read cross the pin in the site's direction? `is_rev` is relative to the node, the
    // site's own orientation is `site_backward`, and the strand is the comparison of the two.
    out.strand = (is_rev == site_backward) ? 0 : 1;

    if (exit_pin) {
        // The reference base is immediately upstream of the pin, which is what `offset` reports.
        out.offset = base_index;
        return out;
    }

    // Entry pin: the reference base is DOWNSTREAM of the pin, so the offset is the read base
    // immediately upstream -- which is the adjacent base on the previous node of the read's WALK,
    // not simply the previous read base that happens to consume a node.
    //
    // Stepping by read index instead is wrong, and wrong in a way that produces a plausible file:
    // a node deleted outright in this read contributes no read bases at all, so an index walk
    // crosses it invisibly and lands on the base before it -- which is the neighbouring snarl's own
    // pin position. That is the collision this design exists to rule out, and it is exactly the
    // hazard the 1 bp boundary node was flagged for, arriving through the deletion rather than
    // through the node's length. Found by the guarantee's own check on real data: 11 collisions in
    // 4.07 M placements, all of this shape.
    //
    // So: one step along the read's walk. Insertions and soft clips inside the neighbouring mapping
    // are skipped, because they consume no node base; a neighbouring mapping that consumes no node
    // base AT ALL means the graph step upstream is deleted here, and the read is refused instead.
    const int64_t neighbour = (out.strand == 0) ? hit - 1 : hit + 1;
    if (neighbour < 0 || neighbour >= (int64_t)path.mapping_size()) {
        // The read begins (or ends) here, so it does not cross the pin.
        ++counters.no_neighbour;
        return out;
    }
    const size_t neighbour_read_start =
        (out.strand == 0) ? before_read_start
                          : hit_read_start + (size_t)mapping_to_length(path.mapping(hit));
    const MappingExtent adjacent = extent_of(path.mapping(neighbour), neighbour_read_start);
    const int64_t candidate =
        (out.strand == 0) ? adjacent.last_consuming : adjacent.first_consuming;
    if (candidate < 0 || candidate >= (int64_t)aln.sequence().size()) {
        ++counters.no_neighbour;
        return out;
    }
    out.offset = candidate;
    return out;
}

////////////////////////////////////////////////////////////////////////////////
// Retention
////////////////////////////////////////////////////////////////////////////////

size_t AnchorSiteEvidence::bytes() const {
    size_t total = sizeof(AnchorSiteEvidence);
    total += reads.capacity() * sizeof(AnchorRead);
    for (const AnchorRead& r : reads) {
        total += r.name.capacity();
    }
    total += rel.capacity() * sizeof(float);
    total += allele_length.capacity() * sizeof(uint32_t);
    return total;
}

////////////////////////////////////////////////////////////////////////////////
// Building anchors from a settled genotype
////////////////////////////////////////////////////////////////////////////////

void build_site_anchors(const AnchorSiteEvidence& evidence, const vector<int>& genotype,
                        const string& snarl_id, double gqn, double explained,
                        const AnchorParams& params, AnchorCounters& counters,
                        vector<AnchorWriter::Anchor>& out) {

    if (genotype.empty() || evidence.reads.empty() || evidence.n_alleles == 0) {
        return;
    }
    for (int a : genotype) {
        if (a < 0 || (size_t)a >= evidence.n_alleles) {
            // A star or missing allele: there is no traversal to assign reads to.
            return;
        }
    }
    bool hom = true;
    for (size_t i = 1; i < genotype.size(); ++i) {
        if (genotype[i] != genotype[0]) {
            hom = false;
        }
    }
    if (hom && params.het_only) {
        return;
    }
    // A negative GQN means the site offered no gap to normalise, which is not the same as zero and
    // must not be filtered as though it were.
    if (params.min_gqn > 0.0 && (gqn < 0.0 || gqn < params.min_gqn)) {
        return;
    }

    // One slot per haplotype of the genotype, collapsed for a homozygote: the pair degenerates to a
    // single anchor holding every read, which carries connectivity but no haplotype information.
    vector<int> slot_allele;
    if (hom) {
        slot_allele.push_back(genotype[0]);
    } else {
        slot_allele.assign(genotype.begin(), genotype.end());
    }
    const size_t n_slots = slot_allele.size();

    // Expected share of the site's reads per haplotype, from the alleles' spelled lengths: the
    // number of read start positions yielding a read overlapping the site from an allele of length L
    // is L + R - 1. Flat when the lengths are unavailable or --flat-mixture is in force.
    vector<double> weight(n_slots, 1.0 / (double)n_slots);
    if (evidence.length_weighted && evidence.mean_read_length > 0.0
        && evidence.allele_length.size() == evidence.n_alleles && n_slots > 1) {
        double total = 0.0;
        vector<double> raw(n_slots, 0.0);
        for (size_t i = 0; i < n_slots; ++i) {
            raw[i] = (double)evidence.allele_length[slot_allele[i]]
                     + (double)evidence.mean_read_length - 1.0;
            if (raw[i] < 1.0) {
                raw[i] = 1.0;
            }
            total += raw[i];
        }
        if (total > 0.0) {
            for (size_t i = 0; i < n_slots; ++i) {
                weight[i] = raw[i] / total;
            }
        }
    }

    // Two anchors per slot: one at each pin. The end pin is skipped where the site's two boundaries
    // are the same node, since the anchors would then be indistinguishable.
    bool degenerate = evidence.start_node == evidence.end_node;
    if (degenerate) {
        ++counters.degenerate_site;
    }
    vector<AnchorWriter::Anchor> start_anchors(n_slots), end_anchors(n_slots);
    for (size_t i = 0; i < n_slots; ++i) {
        start_anchors[i].node = evidence.start_node;
        start_anchors[i].snarl = snarl_id;
        start_anchors[i].slot = (int)i;
        start_anchors[i].gqn = gqn;
        start_anchors[i].explained = explained;
        end_anchors[i] = start_anchors[i];
        end_anchors[i].node = evidence.end_node;
    }

    // Reads sharing a name at this site. Cheap, and the only place the collision is visible: by the
    // time a file exists the two are indistinguishable.
    {
        unordered_map<string, size_t> name_count;
        name_count.reserve(evidence.reads.size() * 2);
        for (const AnchorRead& read : evidence.reads) {
            ++name_count[read.name];
        }
        for (const auto& entry : name_count) {
            if (entry.second > 1) {
                counters.shared_name += entry.second;
            }
        }
    }

    // Reads each slot's end anchor holds that its own start anchor does not. If there are none, the
    // two are joined by all the same reads and the end anchor can offer no linkage the start anchor
    // does not -- see AnchorParams::end_pin_min_new.
    //
    // Per SLOT, not per site. A site-level count is not the same question: anchors are emitted and
    // min-reads-filtered per slot, so a site can hold a read that reaches only the end boundary while
    // the slot carrying it falls below the threshold and the *other* slot's end anchor -- entirely
    // redundant -- survives on the strength of it. Measured on chr20, a site-level count left 4,522
    // of 41,763 surviving end pins carrying nothing of their own.
    vector<size_t> new_at_end(n_slots, 0);

    for (size_t r = 0; r < evidence.reads.size(); ++r) {
        const AnchorRead& read = evidence.reads[r];
        if (!read.start_pin.placed() && !read.end_pin.placed()) {
            continue;
        }
        const double mismap = (double)read.mismap;

        // Responsibilities over the called alleles, plus the mismapping outcome. The mismap term is
        // what bounds the score: at a 0.02 floor a perfectly discriminating read on a balanced het
        // reaches about 14 phred, not 60.
        double best_resp = -1.0;
        size_t best_slot = 0;
        double total = mismap;
        for (size_t i = 0; i < n_slots; ++i) {
            double resp = (1.0 - mismap) * weight[i]
                          * (double)evidence.rel_at(r, (size_t)slot_allele[i]);
            total += resp;
            if (resp > best_resp) {
                best_resp = resp;
                best_slot = i;
            }
        }
        if (!(best_resp > 0.0) || !(total > 0.0)) {
            // The read fits no called allele at all.
            ++counters.off_call;
            continue;
        }

        // Does the read's best-fitting allele over EVERY scored traversal lie inside the call? A
        // read that prefers an allele the genotype does not carry should not be asserted onto one.
        double row_best = 0.0;
        for (size_t a = 0; a < evidence.n_alleles; ++a) {
            row_best = max(row_best, (double)evidence.rel_at(r, a));
        }
        bool called_is_best = false;
        for (size_t i = 0; i < n_slots; ++i) {
            if ((double)evidence.rel_at(r, (size_t)slot_allele[i]) >= row_best) {
                called_is_best = true;
            }
        }
        if (!called_is_best) {
            ++counters.off_call;
            if (!params.keep_off_call) {
                continue;
            }
        }

        const double share = best_resp / total;
        double score = 99.0;
        if (share < 1.0) {
            score = -10.0 * log10(1.0 - share);
        }
        if (score < params.min_read_score) {
            continue;
        }

        AnchorWriter::ReadRow row;
        row.name = read.name;
        row.score = (float)score;

        if (read.end_pin.placed() && !read.start_pin.placed()) {
            ++new_at_end[best_slot];
        }

        bool wrote_start = false;
        if (read.start_pin.placed()) {
            row.offset = read.start_pin.offset;
            row.strand = read.start_pin.strand;
            start_anchors[best_slot].reads.push_back(row);
            wrote_start = true;
        }
        if (!degenerate && read.end_pin.placed()) {
            // A read whose walk goes straight from the start boundary to the end boundary carries an
            // allele deleting the whole site, so it has no bases inside it and both pins land on one
            // position. It is kept at S and dropped here: a pinned position must not appear in more
            // than one anchor, and for such a read the two pins say the same thing anyway.
            if (wrote_start && read.end_pin.offset == read.start_pin.offset
                && read.end_pin.strand == read.start_pin.strand) {
                ++counters.coincident;
            } else {
                row.offset = read.end_pin.offset;
                row.strand = read.end_pin.strand;
                end_anchors[best_slot].reads.push_back(row);
            }
        }
    }

    for (size_t i = 0; i < n_slots; ++i) {
        if (start_anchors[i].reads.size() >= params.min_reads) {
            out.push_back(std::move(start_anchors[i]));
        }
        // Drop this slot's end anchor where it brings too few reads of its own. Decided after the
        // loop because it is a property of what the site actually placed, not of the graph.
        bool drop_end = degenerate;
        if (!drop_end && params.end_pin_min_new > 0
            && new_at_end[i] < params.end_pin_min_new) {
            drop_end = true;
            ++counters.single_pin;
        }
        if (!drop_end && end_anchors[i].reads.size() >= params.min_reads) {
            out.push_back(std::move(end_anchors[i]));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Writer
////////////////////////////////////////////////////////////////////////////////

AnchorWriter::AnchorWriter(size_t threads) : queues(max((size_t)1, threads)) {
}

void AnchorWriter::add(Anchor&& anchor) {
    size_t thread = (size_t)max(0, omp_get_thread_num());
    queues[thread % queues.size()].push_back(std::move(anchor));
}

size_t AnchorWriter::anchor_count() const {
    size_t total = 0;
    for (const auto& q : queues) {
        total += q.size();
    }
    return total;
}

size_t AnchorWriter::read_row_count() const {
    size_t total = 0;
    for (const auto& q : queues) {
        for (const Anchor& a : q) {
            total += a.reads.size();
        }
    }
    return total;
}

bool AnchorWriter::write(const string& path, const string& graph_name, const string& sample,
                         const string& reads_source, double mismap_min,
                         const AnchorParams& params) {
    vector<Anchor> all;
    all.reserve(anchor_count());
    for (auto& q : queues) {
        std::move(q.begin(), q.end(), back_inserter(all));
        q.clear();
    }
    // Ordered by node ID, then by snarl and slot: a node ID does not identify an anchor, because a
    // boundary node is shared by two snarls at 23% of them.
    sort(all.begin(), all.end(), [](const Anchor& a, const Anchor& b) {
        if (a.node != b.node) {
            return a.node < b.node;
        }
        if (a.snarl != b.snarl) {
            return a.snarl < b.snarl;
        }
        return a.slot < b.slot;
    });
    // Within an anchor, by read name, so the file is byte-reproducible however the sweep was
    // scheduled across threads.
    for (Anchor& a : all) {
        sort(a.reads.begin(), a.reads.end(), [](const ReadRow& x, const ReadRow& y) {
            if (x.name != y.name) {
                return x.name < y.name;
            }
            return x.offset < y.offset;
        });
    }

    // Intern the read names.
    //
    // Names are 70.8% of an un-interned chr20 file and each read is placed 3.86 times on average, so
    // this is the file's dominant cost and it is paid over and over. Long reads, the target here,
    // repeat far more than that -- one 15 kb read crosses tens of sites -- so the saving grows with
    // exactly the input this is for.
    //
    // Ids are assigned in sorted-name order rather than first-seen order, so the file is identical
    // however the sweep was scheduled across threads, and the table can be binary-searched.
    map<string, size_t> name_id;
    for (const Anchor& a : all) {
        for (const ReadRow& r : a.reads) {
            name_id.emplace(r.name, 0);
        }
    }
    {
        size_t next = 0;
        for (auto& entry : name_id) {
            entry.second = next++;
        }
    }

    ofstream out(path);
    if (!out) {
        cerr << "error [vg call]: could not open " << path << " for the anchor output" << endl;
        return false;
    }
    out << "#anchors-version\t2\n";
    out << "#graph\t" << graph_name << "\n";
    out << "#sample\t" << sample << "\n";
    out << "#reads\t" << reads_source << "\n";
    out << "#mismap-min\t" << mismap_min << "\n";
    out << "#sites\t" << (params.het_only ? "het" : "het+hom")
        << (params.leaf_only ? ",leaf" : "") << "\n";
    // What was filtered out, so a consumer seeing no anchor at a site it expected one at can tell a
    // threshold from an absence. Written whether or not the defaults were changed.
    out << "#filters\tmin-reads=" << params.min_reads << " min-gqn=" << params.min_gqn
        << " min-read-score=" << params.min_read_score
        << " off-call=" << (params.keep_off_call ? "kept" : "dropped")
        << " end-pin-min-new=" << params.end_pin_min_new << "\n";
    out << "#note\tan anchor is a zero-length pin at the junction between a snarl boundary node "
           "and the site interior. There is no anchor sequence and no length.\n";
    out << "#note\tR rows name their read by an integer id into the #read table above them, which "
           "is sorted by name. The table is the file's only copy of each name.\n";
    out << "#note\toffset is the 0-based index, in the read AS SEQUENCED, of the last base before "
           "the pin read in the site's direction. strand 0: the pin follows that base. strand 1: "
           "it precedes it. For the oriented read, use (read length - 1 - offset).\n";
    out << "#note\tR rows belong to the A row above them and do not repeat its key. The file is "
           "already ordered by node; do not sort it.\n";
    out << "#note\tno (read, strand, offset) appears in more than one anchor -- for distinct "
           "alignments. A read NAME shared by two alignments (paired-end mates share one) can "
           "repeat a position; the shared_name counter reports how many such reads there were.\n";
    out << "#note\tgqn and explained are the site's, so they repeat across its slots. Anything else "
           "about the site is in the VCF, joinable on the snarl column, which is its ID.\n";
    out << "#reads-interned\t" << name_id.size() << "\n";
    out << "#H\tA\tnode\tsnarl\tslot\tgqn\texplained\n";
    out << "#H\tR\tread_id\tstrand\toffset\tscore\n";
    for (const auto& entry : name_id) {
        out << "#read\t" << entry.second << "\t" << entry.first << "\n";
    }
    out << std::fixed;
    for (const Anchor& a : all) {
        out << "A\t" << a.node << "\t" << a.snarl << "\t" << a.slot << "\t";
        if (a.gqn < 0.0) {
            out << ".";
        } else {
            out << std::setprecision(3) << a.gqn;
        }
        out << "\t" << std::setprecision(3) << a.explained << "\n";
        for (const ReadRow& r : a.reads) {
            out << "R\t" << name_id[r.name] << "\t" << (int)r.strand << "\t" << r.offset << "\t"
                << std::setprecision(1) << r.score << "\n";
        }
    }
    return (bool)out;
}

}

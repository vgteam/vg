#include "haplotype_assignment.hpp"

#include <algorithm>
#include <optional>
#include <set>

#include <gbwtgraph/utils.h>
#include <vg/vg.pb.h>

#include "handle.hpp"

namespace vg {

using namespace std;

HaplotypeAssigner::HaplotypeAssigner(const gbwtgraph::GBZ& gbz,
                                     const SnarlDistanceIndex& dist_index)
    : gbz_(&gbz), dist_index_(&dist_index) {
}

vector<net_handle_t> HaplotypeAssigner::enclosing_snarls(nid_t node_id) const {
    vector<net_handle_t> snarls;
    // A root snarl is a root handle, so the walk stops before it.
    net_handle_t net = dist_index_->get_parent(dist_index_->get_node_net_handle(node_id));
    while (!dist_index_->is_root(net)) {
        if (dist_index_->is_snarl(net)) {
            snarls.push_back(dist_index_->canonical(net));
        }
        net = dist_index_->get_parent(net);
    }
    return snarls;
}

bool HaplotypeAssigner::in_snarl(nid_t node_id, const net_handle_t& snarl) const {
    for (const net_handle_t& s : enclosing_snarls(node_id)) {
        if (s == snarl) {
            return true;
        }
    }
    return false;
}

HaplotypeAssignmentResult HaplotypeAssigner::assign(const Alignment& aln,
                                                    bool extend_to_snarls) const {
    HaplotypeAssignmentResult result;

    // Collect the path's nodes as GBWT nodes, with their lengths.
    const auto& path = aln.path();
    if (path.mapping_size() == 0) {
        return result;
    }
    vector<gbwt::node_type> path_nodes;
    vector<size_t> node_bp;
    path_nodes.reserve(path.mapping_size());
    node_bp.reserve(path.mapping_size());
    for (int i = 0; i < path.mapping_size(); ++i) {
        const auto& pos = path.mapping(i).position();
        path_nodes.push_back(gbwt::Node::encode(pos.node_id(), pos.is_reverse()));
        node_bp.push_back(gbz_->graph.get_length(
            gbz_->graph.get_handle(pos.node_id(), pos.is_reverse())));
    }
    const size_t n = path_nodes.size();

    // Split the path into segments. From each start we extend a GBWT search
    // forward until no haplotype goes on to the next node. That node starts
    // the next segment.
    size_t i = 0;
    while (i < n) {
        gbwt::BidirectionalState state = gbz_->index.bdFind(path_nodes[i]);
        if (state.empty()) {
            // No haplotype visits this node in this orientation. Skip it, so it
            // is in no segment.
            ++i;
            continue;
        }
        size_t j = i + 1;
        while (j < n) {
            gbwt::BidirectionalState next =
                gbz_->index.bdExtendForward(state, path_nodes[j]);
            if (next.empty()) {
                break;
            }
            state = next;
            ++j;
        }

        // When one segment is the whole path, walk out from each end while
        // the haplotypes all take the same node. Each step drops the
        // haplotypes that stop there. The walk goes through child snarls and
        // ends at the first node outside the smallest snarl that contains the
        // whole path. That node is one of the snarl's boundary nodes.
        const bool whole_path = (i == 0 && j == n);
        optional<net_handle_t> path_snarl;
        if (extend_to_snarls && whole_path) {
            // The snarls that contain every node, innermost first. Each node's
            // enclosing snarls run from innermost to outermost, so the common
            // ones are a tail of the first node's list.
            vector<net_handle_t> common = enclosing_snarls(gbwt::Node::id(path_nodes.front()));
            for (size_t k = 1; k < n && !common.empty(); ++k) {
                vector<net_handle_t> mine = enclosing_snarls(gbwt::Node::id(path_nodes[k]));
                auto inner = std::find_first_of(common.begin(), common.end(),
                                                mine.begin(), mine.end());
                common.erase(common.begin(), inner);
            }
            // A path in no snarl has no boundary nodes to find.
            if (!common.empty()) {
                path_snarl = common.front();
            }
        }
        if (path_snarl) {
            const net_handle_t& start_snarl = *path_snarl;

            // Walk left.
            for (size_t step = 0; step < max_extension_nodes; ++step) {
                vector<gbwt::BidirectionalState> nexts;
                gbz_->graph.follow_paths(state, /*backward=*/true,
                    [&](const gbwt::BidirectionalState& next) -> bool {
                        if (!next.empty()) {
                            nexts.push_back(next);
                        }
                        return true;
                    });
                if (nexts.empty() || nexts.size() > 1) {
                    break;
                }
                // The backward search holds the new node reversed.
                gbwt::node_type pred = gbwt::Node::reverse(nexts[0].backward.node);
                nid_t pred_id = gbwt::Node::id(pred);
                if (!in_snarl(pred_id, start_snarl)) {
                    result.left_boundary_node = pred_id;
                    break;
                }
                state = nexts[0];
            }

            // Walk right.
            for (size_t step = 0; step < max_extension_nodes; ++step) {
                vector<gbwt::BidirectionalState> nexts;
                gbz_->graph.follow_paths(state, /*backward=*/false,
                    [&](const gbwt::BidirectionalState& next) -> bool {
                        if (!next.empty()) {
                            nexts.push_back(next);
                        }
                        return true;
                    });
                if (nexts.empty() || nexts.size() > 1) {
                    break;
                }
                gbwt::node_type succ = nexts[0].forward.node;
                nid_t succ_id = gbwt::Node::id(succ);
                if (!in_snarl(succ_id, start_snarl)) {
                    result.right_boundary_node = succ_id;
                    break;
                }
                state = nexts[0];
            }
        }

        // Record the segment and the sequences that carry it.
        HaplotypeAssignmentResult::Segment seg;
        seg.start_mapping = i;
        seg.end_mapping   = j;
        for (size_t k = i; k < j; ++k) {
            seg.covered_bp += node_bp[k];
        }
        seg.candidate_count = state.size();
        // gbwt::size_type is not size_t on macOS, so std::min needs the type named.
        const size_t cap = std::min<size_t>(state.size(), max_haplotypes_to_report);
        const gbwt::SearchState& fwd = state.forward;
        seg.seq_ids.reserve(cap);
        for (size_t t = 0; t < cap; ++t) {
            seg.seq_ids.push_back(gbz_->index.locate(fwd.node, fwd.range.first + t));
        }
        result.segments.push_back(std::move(seg));

        i = j;
    }

    return result;
}

const HaplotypeAssignmentResult::Segment* longest_segment(const HaplotypeAssignmentResult& result) {
    // max_element returns the first of equal segments.
    auto longest = std::max_element(result.segments.begin(), result.segments.end(),
        [](const HaplotypeAssignmentResult::Segment& a, const HaplotypeAssignmentResult::Segment& b) {
            return a.covered_bp < b.covered_bp;
        });
    return longest == result.segments.end() ? nullptr : &*longest;
}

string haplotype_name(const gbwt::Metadata& metadata, gbwt::size_type path_id) {
    gbwt::PathName path = metadata.path(path_id);
    string sample = (metadata.samples() > 0 && path.sample < metadata.samples())
                    ? metadata.sample(path.sample) : std::to_string(path.sample);
    string contig = (metadata.contigs() > 0 && path.contig < metadata.contigs())
                    ? metadata.contig(path.contig) : std::to_string(path.contig);
    // Only a reference path's name has no fragment, so name it as one.
    return PathMetadata::create_path_name(PathSense::REFERENCE, sample, contig, path.phase,
                                          PathMetadata::NO_PHASE_BLOCK, PathMetadata::NO_SUBRANGE);
}

vector<string> haplotype_names(const gbwtgraph::GBZ& gbz,
                               const vector<gbwt::size_type>& seq_ids) {
    std::set<string> uniq;
    const gbwt::GBWT& index = gbz.index;
    const bool have_md = index.hasMetadata();
    const gbwt::Metadata& md = index.metadata;
    for (gbwt::size_type s : seq_ids) {
        gbwt::size_type path_id = s / 2;
        if (!have_md || path_id >= md.paths()) {
            uniq.insert(std::to_string(s));
            continue;
        }
        uniq.insert(haplotype_name(md, path_id));
    }
    return vector<string>(uniq.begin(), uniq.end());
}

pair<string, bool> choose_representative(const gbwtgraph::GBZ& gbz,
                                         const vector<gbwt::size_type>& seq_ids) {
    const gbwt::GBWT& index = gbz.index;
    if (seq_ids.empty() || !index.hasMetadata()) {
        return {string(), false};
    }
    const gbwt::Metadata& md = index.metadata;
    const gbwtgraph::sample_name_set& refs = gbz.get_reference_samples();

    bool any = false, ref_found = false;
    gbwt::size_type best_any = 0, best_ref = 0;
    for (gbwt::size_type s : seq_ids) {
        gbwt::size_type path_id = s / 2;
        if (path_id >= md.paths()) continue;
        gbwt::PathName pn = md.path(path_id);
        string sample = (md.samples() > 0 && pn.sample < md.samples())
                        ? md.sample(pn.sample) : string();
        bool is_ref = !sample.empty() &&
                      (refs.count(sample) > 0 ||
                       sample == gbwtgraph::GENERIC_PATH_SAMPLE_NAME);
        if (!any || s < best_any) { best_any = s; any = true; }
        if (is_ref && (!ref_found || s < best_ref)) { best_ref = s; ref_found = true; }
    }
    if (!any) {
        return {string(), false};
    }
    gbwt::size_type chosen = ref_found ? best_ref : best_any;
    return {haplotype_name(md, chosen / 2), ref_found};
}

} // namespace vg

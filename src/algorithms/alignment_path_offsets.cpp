#include "alignment_path_offsets.hpp"


namespace vg {
namespace algorithms {

unordered_map<path_handle_t, vector<pair<size_t, bool> > >
alignment_path_offsets(const PathPositionHandleGraph& graph,
                       const Alignment& aln,
                       bool just_min,
                       bool nearby,
                       size_t search_limit) {
    if (nearby && search_limit == 0) {
        // Fill in the search limit
        search_limit = aln.sequence().size();
    }
    unordered_map<path_handle_t, vector<pair<size_t, bool> > > offsets;
    for (auto& mapping : aln.path().mapping()) {
        // How many bases does this Mapping cover over?
        size_t mapping_width = mapping_from_length(mapping);
        if (mapping_width == 0 && !nearby) {
            // Just skip over this mapping; it touches no bases.
            continue;
        }
        // We may have to consider both the starts and ends of mappings
        vector<bool> end = {false};
        if (just_min && !nearby) {
            // We want the min actually touched position along each path. It
            // could come from the Mapping start or the Mapping end.
            end.push_back(true);
        }
        // Find the position of this end of this mapping
        pos_t mapping_pos = make_pos_t(mapping.position());
        // Find the positions for this end of this Mapping
        auto pos_offs = algorithms::nearest_offsets_in_paths(&graph, mapping_pos, nearby ? search_limit : -1);
        for (auto look_at_end : end) {
            // For the start and the end of the Mapping, as needed
            for (auto& p : pos_offs) {
                // For each path, splice the list of path positions for this Mapping
                // onto the end of the list of positions we found in that path
                auto& v = offsets[p.first];
                for (pair<size_t, bool>& y : p.second) {
                    v.emplace_back(y.second ? y.first - mapping_width : y.first,
                                   y.second);
                }
            }
        }
    }
    if (!nearby && offsets.empty()) {
        // find the nearest if we couldn't find any before
        return alignment_path_offsets(graph, aln, just_min, true, search_limit);
    }
    if (just_min) {
        // We need the minimum position for each path
        for (auto& p : offsets) {
            auto& v = p.second;
            auto m = *min_element(v.begin(), v.end(),
                                  [](const pair<size_t, bool>& a,
                                     const pair<size_t, bool>& b)
                                  { return a.first < b.first; });
            v.clear();
            v.push_back(m);
        }
    }
    return offsets;
}

void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, Alignment& aln, size_t search_limit) {
    if (!aln.refpos_size()) {
        unordered_map<path_handle_t, vector<pair<size_t, bool> > > positions = alignment_path_offsets(graph, aln, true, false, search_limit);
        for (auto& path : positions) {
            Position* refpos = aln.add_refpos();
            refpos->set_name(graph.get_path_name(path.first));
            refpos->set_offset(path.second.front().first);
            refpos->set_is_reverse(path.second.front().second);
        }
    }
}

void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, vector<Alignment>& alns, size_t search_limit) {
    for (auto& aln : alns) annotate_with_initial_path_positions(graph, aln, search_limit);
}

}
}

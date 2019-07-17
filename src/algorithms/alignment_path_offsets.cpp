#include "alignment_path_offsets.hpp"


namespace vg {
namespace algorithms {

// gives just the first positions we encounter
unordered_map<path_handle_t, vector<pair<size_t, bool> > > alignment_path_offsets(const PathPositionHandleGraph& graph, const Alignment& aln) {
    unordered_map<path_handle_t, vector<pair<size_t, bool> > > positions;
    for (auto& mapping : aln.path().mapping()) {
        if (mapping.has_position()) {
            graph.get_handle(mapping.position().node_id());
            pos_t pos = make_pos_t(mapping.position());
            auto offsets = algorithms::nearest_offsets_in_paths(&graph, pos, -1);
            for (auto& path : offsets) {
                auto& path_positions = positions[path.first];
                for (auto& pos : path.second) {
                    path_positions.push_back(pos);
                }
            }
            break;
        }
    }
    return positions;
}

void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, Alignment& aln) {
    unordered_map<path_handle_t, vector<pair<size_t, bool> > > positions = alignment_path_offsets(graph, aln);
    for (auto& path : positions) {
        Position* refpos = aln.add_refpos();
        refpos->set_name(graph.get_path_name(path.first));
        refpos->set_offset(path.second.front().first);
        refpos->set_is_reverse(path.second.front().second);
    }
}

void annotate_with_initial_path_positions(const PathPositionHandleGraph& graph, vector<Alignment>& alns) {
    for (auto& aln : alns) annotate_with_initial_path_positions(graph, aln);
}

}
}

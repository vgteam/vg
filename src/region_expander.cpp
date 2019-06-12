#include "region_expander.hpp"

namespace vg {

    RegionExpander::RegionExpander(XG* xg_index, const SnarlManager* snarl_manager) :
        xg_index(xg_index), snarl_manager(snarl_manager)
    {
        // Nothing to do
    }

    map<pair<id_t, bool>, pair<uint64_t,uint64_t >> RegionExpander::expanded_subgraph(const GFFRecord& gff_record) {
        
        map<pair<id_t, bool>, pair<uint64_t, uint64_t>> return_val;

        vector<pair<id_t, bool>> interval_subpath;
        
        assert(gff_record.start != -1 && gff_record.end != -1 && gff_record.start <= gff_record.end);
        
        if (!xg_index->has_path(gff_record.sequence_id)) {
            cerr << "error [RegionExpander] cannot expand genomic interval, graph does not contain path with name: " << gff_record.sequence_id << endl;
            exit(1);
        }
        const XGPath& path = xg_index->get_path(gff_record.sequence_id);
        
        // walk along the path for the interval and add the corresponding nodes to the subgraph
        
        size_t offset = path.offset_at_position(gff_record.start);
        id_t node_id = path.node(offset);
        bool is_rev = path.directions[offset];
        size_t node_length = xg_index->node_length(node_id);
        
        size_t at_pos = path.positions[offset];
        
        interval_subpath.emplace_back(node_id, is_rev);
        return_val[make_pair(node_id, is_rev)] = pair<uint64_t, uint64_t>(gff_record.start - at_pos,
                                                                          node_length);
        at_pos += node_length;
        while (at_pos <= gff_record.end) {
            offset++;
            node_id = path.node(offset);
            is_rev = path.directions[offset];
            node_length = xg_index->node_length(node_id);
            
            interval_subpath.emplace_back(node_id, is_rev);
            return_val[make_pair(node_id, is_rev)] = pair<uint64_t, uint64_t>(0, node_length);
            
            at_pos += node_length;
        }
        
        return_val[make_pair(node_id, is_rev)].second = gff_record.end - (at_pos - node_length) + 1;
        
        // walk along the path and identify snarls that have both ends on the path
        
        unordered_set<const Snarl*> entered_snarls;
        unordered_set<const Snarl*> completed_snarls;
        
        for (size_t i = 0; i + 1 < interval_subpath.size(); i++) {
            const Snarl* snarl_out = snarl_manager->into_which_snarl(interval_subpath[i].first,
                                                                     !interval_subpath[i].second);
            
            if (snarl_out) {
                if (entered_snarls.count(snarl_out)) {
                    completed_snarls.insert(snarl_out);
                }
            }
            
            const Snarl* snarl_in = snarl_manager->into_which_snarl(interval_subpath[i].first,
                                                                    interval_subpath[i].second);
            
            if (snarl_in) {
                entered_snarls.insert(snarl_in);
            }
        }
        
        // traverse the subgraph in each snarl and add it to the annotation
        for (const Snarl* snarl : completed_snarls) {
            // orient the snarl to match the orientation of the annotation
            handle_t oriented_start, oriented_end;
            if (return_val.count(pair<id_t, bool>(snarl->start().node_id(),
                                                  snarl->start().backward()))) {
                
                oriented_start = xg_index->get_handle(snarl->start().node_id(),
                                                      snarl->start().backward());
                
                oriented_end = xg_index->get_handle(snarl->end().node_id(),
                                                    snarl->end().backward());
                
            }
            else {
                
                oriented_start = xg_index->get_handle(snarl->end().node_id(),
                                                      !snarl->end().backward());
                
                oriented_end = xg_index->get_handle(snarl->start().node_id(),
                                                    !snarl->start().backward());
            }
            
            // mark all the exits and the entry point as untraversable
            unordered_set<handle_t> stacked{oriented_start, oriented_end, xg_index->flip(oriented_start)};
            vector<handle_t> stack{oriented_start};
            
            // make an index for jumping over the inside of child snarls
            unordered_map<handle_t, handle_t> child_snarl_skips;
            for (const Snarl* child : snarl_manager->children_of(snarl)) {
                handle_t start = xg_index->get_handle(child->start().node_id(),
                                                      child->start().backward());
                handle_t end = xg_index->get_handle(child->end().node_id(),
                                                    child->end().backward());
                
                child_snarl_skips[start] = end;
                child_snarl_skips[xg_index->flip(end)] = xg_index->flip(start);
            }
            
            // traverse the subgraph and add it to the return value
            while (!stack.empty()) {
                handle_t handle = stack.back();
                stack.pop_back();
                
                pair<id_t, bool> trav = make_pair(xg_index->get_id(handle),
                                                  xg_index->get_is_reverse(handle));
                
                if (!return_val.count(trav)) {
                    return_val[trav] = pair<uint64_t, uint64_t>(0, xg_index->get_length(handle));
                }
                
                if (child_snarl_skips.count(handle)) {
                    // skip over the internals of the child snarl we're pointing into
                    handle_t next = child_snarl_skips[handle];
                    if (!stacked.count(next)) {
                        stack.push_back(next);
                        stacked.insert(next);
                    }
                }
                else {
                    // traverse edges
                    xg_index->follow_edges(handle, false, [&](const handle_t& next) {
                        if (!stacked.count(next)) {
                            stack.push_back(next);
                            stacked.insert(next);
                        }
                    });
                }
            }
        }
        
        if (gff_record.strand_is_rev) {
            // we did all our queries with respect to the forward strand, flip it back to the reverse
            map<pair<id_t, bool>, pair<uint64_t, uint64_t>> reversed_map;
            for (const auto& record : return_val) {
                uint64_t node_length = xg_index->node_length(record.first.first);
                reversed_map[make_pair(record.first.first, !record.first.second)] = make_pair(node_length - record.second.second,
                                                                                              node_length - record.second.first);
            }
            return_val = move(reversed_map);
        }
        
        return return_val;
    }
}

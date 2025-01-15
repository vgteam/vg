#include "masker.hpp"

#include <algorithm>
#include <unordered_set>

#include "integrated_snarl_finder.hpp"

//#define debug

namespace vg {

Masker::Masker(gbwtgraph::GBWTGraph& _graph, const SnarlManager* _snarl_manager) : graph(&_graph), snarl_manager(_snarl_manager) {
    if (!snarl_manager) {
        init_snarls();
    }
}

Masker::Masker(MutablePathDeletableHandleGraph& _graph, const SnarlManager* _snarl_manager) : graph(&_graph), snarl_manager(_snarl_manager) {
    if (!snarl_manager) {
        init_snarls();
    }
}

void Masker::init_snarls() {
    
    // construct our own snarl finder
    IntegratedSnarlFinder snarl_finder(*graph);
    own_snarl_manager.reset(new SnarlManager(std::move(snarl_finder.find_snarls_parallel())));
    snarl_manager = own_snarl_manager.get();
}

void Masker::mask_sequences(const std::vector<std::tuple<std::string, size_t, size_t>>& regions) const {
    
#ifdef debug
    std::cerr << "input regions:" << std::endl;
    for (const auto& region : regions) {
        std::cerr << '\t' << std::get<0>(region) << '\t' << std::get<1>(region) << '\t' << std::get<2>(region) << std::endl;
    }
#endif
    
    auto sorted_regions = regions;
    std::sort(sorted_regions.begin(), sorted_regions.end());
    
    // consolidate the regions in place
    size_t i = 0;
    for (size_t j = 1; j < sorted_regions.size(); ++j) {
        auto& into = sorted_regions[i];
        auto& here = sorted_regions[j];
        if (std::get<0>(here) != std::get<0>(into) || std::get<1>(here) > std::get<2>(into)) {
            // they do not overlap
            ++i;
            if (i != j) {
                // move merged the next interval into the prefix of the set
                sorted_regions[i] = std::move(here);
            }
        }
        else {
            // they overlap
            std::get<2>(into) = std::max(std::get<2>(into), std::get<2>(here));
        }
    }
    // trim off remnants of the intervals that were merged away
    if (i + 1 < sorted_regions.size()) {
        sorted_regions.resize(i + 1);
    }
    
#ifdef debug
    std::cerr << "sorted and consolidated regions:" << std::endl;
    for (const auto& region : sorted_regions) {
        std::cerr << '\t' << std::get<0>(region) << '\t' << std::get<1>(region) << '\t' << std::get<2>(region) << std::endl;
    }
#endif
    
    // collect the path names intervals associated with each path handle
    std::vector<std::tuple<std::string, size_t, path_handle_t>> path_intervals;
    graph->for_each_path_matching({PathSense::REFERENCE, PathSense::GENERIC}, {}, {}, [&](const path_handle_t& path_handle) {
        
        std::string annotated_path_name = graph->get_path_name(path_handle);
        PathSense sense;
        std::string sample;
        std::string locus;
        size_t haplotype;
        size_t phase_block;
        subrange_t subrange;
        graph->parse_path_name(annotated_path_name, sense, sample, locus, haplotype, phase_block, subrange);
        
        auto non_range_path_name = graph->create_path_name(sense, sample, locus, haplotype, phase_block,
                                                           PathMetadata::NO_SUBRANGE);
        
        path_intervals.emplace_back(non_range_path_name,
                                    subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first,
                                    path_handle);
    });
    std::sort(path_intervals.begin(), path_intervals.end());
    
#ifdef debug
    std::cerr << "path intervals:" << std::endl;
    for (const auto& interval : path_intervals) {
        std::cerr << '\t' << std::get<0>(interval) << '\t' << std::get<1>(interval) << std::endl;
    }
#endif
    
    // records of (first step, offset on first, last step, offset on last)
    std::vector<std::tuple<step_handle_t, size_t, step_handle_t, size_t>> step_intervals;
    
    // iterate along both lists to find the bounding steps where they overlap
    size_t p_idx = 0, r_idx = 0;
    while (r_idx < sorted_regions.size() && p_idx < path_intervals.size()) {
        const auto& region_here = sorted_regions[r_idx];
        const auto& path_interval_here = path_intervals[p_idx];
        // check if they are on different paths
        if (std::get<0>(region_here) < std::get<0>(path_interval_here)) {
            ++r_idx;
        }
        else if (std::get<0>(region_here) > std::get<0>(path_interval_here)) {
            ++p_idx;
        }
        else {
            size_t p_end = std::numeric_limits<size_t>::max();
            if (p_idx + 1 < path_intervals.size() &&
                std::get<0>(path_intervals[p_idx + 1]) == std::get<0>(path_interval_here)) {
                p_end = std::get<1>(path_intervals[p_idx + 1]);
            }
            // check if they are on disjoint intervals
            if (std::get<2>(region_here) < std::get<1>(path_interval_here)) {
                ++r_idx;
            }
            else if (p_end < std::get<1>(region_here)) {
                ++p_idx;
            }
            else {
#ifdef debug
                std::cerr << "possible overlap at path interval " << p_idx << " (" << graph->get_path_name(std::get<2>(path_interval_here)) << ") and region " << r_idx << std::endl ;
#endif
                // we may have reached an overlap, scan the path
                size_t offset = std::get<1>(path_interval_here);
                bool inside_region = false;
                auto step = graph->path_begin(std::get<2>(path_interval_here));
                auto prev_step = step; // just a placeholder for now
                auto end = graph->path_end(std::get<2>(path_interval_here));
                while (step != end && r_idx < sorted_regions.size() && std::get<1>(sorted_regions[r_idx]) < p_end) {

                    size_t end_offset = offset + graph->get_length(graph->get_handle_of_step(step));
                                        
                    if (!inside_region &&
                        ((offset <= std::get<1>(sorted_regions[r_idx]) && std::get<1>(sorted_regions[r_idx]) < end_offset) ||
                         (std::get<1>(sorted_regions[r_idx]) <= offset && offset < std::get<2>(sorted_regions[r_idx])))) {
                        // we're overlapping a new region for this path interval
                        
                        size_t node_offset = std::get<1>(sorted_regions[r_idx]) >= offset ? std::get<1>(sorted_regions[r_idx]) - offset : 0;
                        
                        // the latter two are placeholders for now
                        step_intervals.emplace_back(step, node_offset, step, -1);
                        
                        inside_region = true;
                    }
                    
                    if (inside_region && std::get<2>(sorted_regions[r_idx]) <= end_offset) {
                        // we're exiting the region we've been traversing
                        std::get<2>(step_intervals.back()) = step;
                        std::get<3>(step_intervals.back()) = std::get<2>(sorted_regions[r_idx]) - offset;
                        
                        inside_region = false;
                        ++r_idx;
                    }
                    else {
                        prev_step = step;
                        step = graph->get_next_step(step);
                        offset = end_offset;
                    }
                }
                
                // we never left and closed out this region (it must end after this path fragment)
                if (inside_region) {
                    std::get<2>(step_intervals.back()) = prev_step;
                    std::get<3>(step_intervals.back()) = graph->get_length(graph->get_handle_of_step(prev_step));
                }
                
                ++p_idx;
            }
        }
    }
    
#ifdef debug
    std::cerr << "step intervals:" << std::endl;
    for (const auto& interval : step_intervals) {
        std::cerr << '\t' << graph->get_id(graph->get_handle_of_step(std::get<0>(interval))) << '\t' << std::get<1>(interval) << '\t' << graph->get_id(graph->get_handle_of_step(std::get<2>(interval))) << '\t' << std::get<3>(interval) << std::endl;
    }
#endif
    
    MutablePathDeletableHandleGraph* mutable_graph = dynamic_cast<MutablePathDeletableHandleGraph*>(graph);
    gbwtgraph::GBWTGraph* gbwt_graph = dynamic_cast<gbwtgraph::GBWTGraph*>(graph);
    if (mutable_graph) {
        // lambda to change sequence of a mutable handle graph
        std::function<void(handle_t, const std::string&)> apply = [&](handle_t handle, const std::string& sequence) {
            mutable_graph->change_sequence(handle, sequence);
        };
        apply_mask_sequences(step_intervals, apply);
    }
    else if (gbwt_graph) {
        // lambda to change sequence of a gbwt graph
        std::function<void(handle_t, const std::string&)> apply = [&](handle_t handle, const std::string& sequence) {
            
            auto rev_sequence = reverse_complement(sequence);
            
            auto fwd_handle = gbwt_graph->forward(handle);
            auto rev_handle = gbwt_graph->flip(fwd_handle);
            // GBWT graph stores sequence in both orientations, so we have to modify both
            for (bool fwd : {true, false}) {
                auto local_handle = fwd ? fwd_handle : rev_handle;
                const auto& local_seq = fwd ? sequence : rev_sequence;
                // FIXME: have to replicate GBWTGraph::node_offset's function because it's private
                size_t node_offset = gbwt_graph->handle_to_node(local_handle) - gbwt_graph->index->firstNode();
                size_t seq_offset = gbwt_graph->sequences.index[node_offset];
                for (size_t i = 0; i < sequence.size(); ++i) {
                    gbwt_graph->sequences.strings[seq_offset + i] = local_seq[i];
                }
            }
            
            
        };
        
        apply_mask_sequences(step_intervals, apply);
    }
    else {
        std::cerr << "error: could not identify input graph to Masker as a GBZ/GBWTGraph or a DeletableHandleGraph" << std::endl;
        exit(1);
    }
}

void Masker::apply_mask_sequences(const std::vector<std::tuple<step_handle_t, size_t, step_handle_t, size_t>>& mask_intervals,
                                  const std::function<void(handle_t, const std::string&)>& apply_mask) const {
    
    // lambda to mask out part of a node, with an interval given in the handle's orientation
    auto partial_mask = [&](handle_t handle, size_t begin, size_t end) {
        auto sequence = graph->get_sequence(graph->forward(handle));
        size_t i, n;
        if (graph->get_is_reverse(handle)) {
            i = sequence.size() - end;
            n = sequence.size() - begin;
        }
        else {
            i = begin;
            n = end;
        }
        for (; i < n; ++i) {
            sequence[i] = 'N';
        }
        apply_mask(handle, sequence);
    };
    
    for (const auto& mask_interval : mask_intervals) {
        
        if (std::get<0>(mask_interval) == std::get<2>(mask_interval)) {
            // the entire interval is contained within one node
            
            partial_mask(graph->get_handle_of_step(std::get<0>(mask_interval)),
                         std::get<1>(mask_interval), std::get<3>(mask_interval));
        }
        else {
            
            // returns true if the handle points into a non-trivial snarl
            auto into_nontrivial_snarl = [&](handle_t handle) -> bool {
                if (graph->get_degree(handle, false) < 2) {
                    // any snarl that this points into will be trivial
                    return false;
                }
                return snarl_manager->into_which_snarl(graph->get_id(handle), graph->get_is_reverse(handle));
            };
            
            // the depth (relative to where we enter) in the snarl tree
            int64_t level = 0;
            // the prefix of recorded snarls that cannot be contained in a future snarl
            size_t uncontained_prefix = 0;
            // snarls for which both ends lie on this path interval
            std::vector<std::tuple<int64_t, handle_t, handle_t>> snarls_discovered;
            // stack of the snarl entrypoints we've passed
            std::vector<handle_t> snarl_source_stack;
            
            // handle snarls that we may be exiting
            auto process_incoming = [&](handle_t handle) {
                if (into_nontrivial_snarl(graph->flip(handle))) {
                    --level;
                    if (!snarl_source_stack.empty()) {
                        // clear out any snarls contained in this one
                        while (snarls_discovered.size() > uncontained_prefix && std::get<0>(snarls_discovered.back()) > level) {
                            snarls_discovered.pop_back();
                        }
                        // record the snarl
                        snarls_discovered.emplace_back(level, snarl_source_stack.back(), handle);
                        snarl_source_stack.pop_back();
#ifdef debug
                        std::cerr << "found traversed snarl between " << graph->get_id(std::get<1>(snarls_discovered.back())) << " and " << graph->get_id(std::get<2>(snarls_discovered.back())) << '\n';
#endif
                        
                        if (snarl_source_stack.empty()) {
                            // nothing can bracket this snarl on both sides anymore
                            uncontained_prefix = snarls_discovered.size();
                        }
                    }
                }
            };
            
            // handle snarls we may be entering
            auto process_outgoing = [&](handle_t handle) {
                if (into_nontrivial_snarl(handle)) {
                    ++level;
                    snarl_source_stack.push_back(handle);
                }
            };
            
            std::unordered_set<handle_t> to_mask;
            
            // walk the path interval
            
#ifdef debug
            std::cerr << "starting interval traversal on node " << graph->get_id(graph->get_handle_of_step(std::get<0>(mask_interval))) << '\n';
#endif
            
            // first partial node
            partial_mask(graph->get_handle_of_step(std::get<0>(mask_interval)),
                         std::get<1>(mask_interval),
                         graph->get_length(graph->get_handle_of_step(std::get<0>(mask_interval))));
            process_outgoing(graph->get_handle_of_step(std::get<0>(mask_interval)));
            
            // full nodes
            for (auto step = graph->get_next_step(std::get<0>(mask_interval));
                 step != std::get<2>(mask_interval); step = graph->get_next_step(step)) {
                
                handle_t handle = graph->get_handle_of_step(step);
                to_mask.insert(graph->forward(handle));
                
#ifdef debug
                std::cerr << "at node " << graph->get_id(handle) << '\n';
                std::cerr << "snarl stack state:\n";
                for (auto handle : snarl_source_stack) {
                    std::cerr << '\t' << graph->get_id(handle) << '\n';
                }
#endif
                
                process_incoming(handle);
                process_outgoing(handle);
            }
            
#ifdef debug
            std::cerr << "ending interval traversal on node " << graph->get_id(graph->get_handle_of_step(std::get<2>(mask_interval))) << '\n';
#endif
            
            // last partial node
            process_incoming(graph->get_handle_of_step(std::get<2>(mask_interval)));
            partial_mask(graph->get_handle_of_step(std::get<2>(mask_interval)),
                         0, std::get<3>(mask_interval));
            
            // traverse the snarls
            for (const auto& snarl : snarls_discovered) {
                
                // DFS starting from one end
                std::unordered_set<handle_t> seen{graph->forward(std::get<1>(snarl)), graph->forward(std::get<2>(snarl))};
                std::vector<handle_t> stack;
                auto lambda = [&](const handle_t& next) {
                    auto fwd = graph->forward(next);
                    bool is_new = seen.insert(fwd).second;
                    if (is_new) {
                        stack.push_back(fwd);
                        to_mask.insert(fwd);
#ifdef debug
                        std::cerr << "inserting mask node from snarl traversal: " << graph->get_id(fwd) << '\n';
#endif
                    }
                };
                graph->follow_edges(std::get<1>(snarl), false, lambda);
                while (!stack.empty()) {
                    handle_t here = stack.back();
                    stack.pop_back();
                    for (bool to_left : {true, false}) {
                        graph->follow_edges(here, to_left, lambda);
                    }
                }
            }
            
#ifdef debug
            std::cerr << "final masking nodes:\n";
            for (auto handle : to_mask) {
                std::cerr << '\t' << graph->get_id(handle) << '\n';
            }
#endif
            
            for (auto handle : to_mask) {
                apply_mask(handle, std::string(graph->get_length(handle), 'N'));
            }
        }
    }
}

}

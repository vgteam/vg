//
//  multipath_alignment.cpp
//  
//
//  Created by Jordan Eizenga on 10/10/16.
//
//

#include "multipath_alignment.hpp"

using namespace std;
namespace vg {
    
    void identify_start_subpaths(MultipathAlignment& multipath_aln) {
        
        // remove start nodes if there are any (avoids doubling them if this function is used liberally)
        multipath_aln.clear_start();
        
        // label nodes with incoming edges
        vector<bool> has_incoming_edge = vector<bool>(multipath_aln.subpath_size(), false);
        for (int i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            for (int j = 0; j < subpath.next_size(); j++) {
                has_incoming_edge[subpath.next(j)] = true;
            }
        }
        
        // construct list of starts
        for (int i = 0; i < multipath_aln.subpath_size(); i++) {
            if (!has_incoming_edge[i]) {
                multipath_aln.add_start(i);
            }
        }
    }
    
    // TODO: should I switch to a custom linked list implementation so that I can keep static locations
    // in the node-to-location map instead of iterators when I swap segments of haplotypes?
    
    int32_t optimal_score_on_path(const MultipathAlignment& multipath_aln, const VG& graph, const list<NodeTraversal>& path,
                                  const unordered_map<int64_t, vector<list<NodeTraversal>::iterator>& node_locations) {
        
        int32_t optimal_score = 0;
        
        // find the places in the path where the alignment might start
        set<list<NodeTraversal>::iterator> candidate_start_positions;
        for (int i = 0; i < multipath_aln.start_size(); i++) {
            // a starting subpath in the multipath alignment
            Subpath* start_subpath = multipath_aln.subpath(multipath_aln.start(i));
            const Position& start_pos = start_subpath.path().mapping(0).position();
            
            // add each location the start nodes occur in the path to the candidate starts
            if (node_locations.count(start_pos.node_id())) {
                const vector<list<NodeTraversal>::iterator>& subpath_matches = node_locations[start_pos.node_id()]
                candidate_start_positions.insert(subpath_matches.begin(), subpath_matches.end());
            }
        }
        
        // check alignments starting at each node in the path that has a source subpath starting on it
        for (list<NodeTraversal>::iterator& path_start_node : candidate_start_positions) {
            
            // match up forward and backward traversal on the path to forward and backward traversal through
            // the multipath alignment
            auto move_forward = [path](list<NodeTraversal>::iterator& path_node) { path_node++; };
            auto move_backward = [path](list<NodeTraversal>::iterator& path_node) {
                if (path_node == path.begin()) {
                    path_node = path.end();
                }
                else {
                    path--;
                }
            };
            if (start_pos.is_reverse() != (*path_node).backward) {
                swap(move_forward, move_backward);
            }
            
            // initialize 2 dynamic programming structures:
            // pointer to place in path corresponding to the beginning of a subpath
            vector<list<NodeTraversal>::iterator> subpath_nodes = vector<list<NodeTraversal>::iterator>(multipath_aln.subpath_size(),
                                                                                                        path.end());
            // score of the best preceding path before this subpath
            vector<int32_t> subpath_prefix_score = vector<int32_t>(multipath_aln.subpath_size(), 0);
            
            // set DP base case with the subpaths that start at this path node
            for (int i = 0; i < multipath_aln.start_size(); i++) {
                Subpath* start_subpath = multipath_aln.subpath(multipath_aln.start(i));
                const Position& start_pos = first_path_position(start_subpath.path());
                
                if (start_pos.node_id() == (*path_start_node).node->id()) {
                    subpath_nodes[multipath_aln.start(i)] = path_start_node;
                }
            }
            
            for (int i = 0; i < multipath_aln.subpath_size(); i++) {
                list<NodeTraversal>::iterator subpath_node = subpath_nodes[i];
                
                // this subpath may be unreachable from subpaths consistent with the path
                if (subpath_node == path.end()) {
                    continue;
                }
                
                const Subpath& subpath = multipath_aln.subpath(i);
                
                // iterate through mappings in this subpath (assumes one mapping per node)
                bool subpath_follows_path = true;
                for (int j = 0; j < subpath.path().mapping_size(); j++, move_forward(subpath_node)) {
                    // check if mapping corresponds to the next node in the path
                    if (subpath.path().mapping(j).position().node_id() != (*subpath_node).node->id()) {
                        subpath_follows_path = false;
                        break;
                    }
                }
                
                // update subsequent subpaths or record completed alignment if this subpath didn't fall off the path
                if (subpath_follows_path) {
                    int32_t extended_prefix_score = subpath_prefix_score[j] + subpath.score();
                    if (subpath.next_size() == 0) {
                        // reached a sink subpath (thereby completing an alignment), check for optimality
                        if (extended_prefix_score > optimal_score) {
                            optimal_score = extended_prefix_score;
                        }
                    }
                    else {
                        // edge case: check if subpath_node was improperly incremented from a mapping that ended in the
                        // middle of a node
                        Position end_pos = last_path_position(subpath.path());
                        if (end_pos.offset() != graph.get_node(end_pos.node_id())->sequence().length()) {
                            move_backward(subpath_node);
                        }
                        
                        for (int j = 0; j < subpath.next_size(); j++) {
                            if (subpath_prefix_score[subpath.next(j)] < extended_prefix_score) {
                                subpath_prefix_score[subpath.next(j)] = extended_prefix_score;
                                subpath_nodes[subpath.next(j)] = subpath_node;
                            }
                        }
                    }
                }
            }
        }
        
        return optimal_score;
    }
    
    void optimal_alignment(const MultipathAlignment& multipath_aln, Alignment& aln_out) {
        
        // transfer read information over to alignment
        aln_out.set_sequence(multipath_aln.sequence());
        aln_out.set_quality(multipath_aln.set_quality());
        aln_out.set_name(multipath_aln.name());
        aln_out.set_sample_name(multipath_aln.sample_name());
        aln_out.set_read_group(multipath_aln.read_group());
        
        // initialize DP structures
        
        // score of the optimal alignment ending in this subpath
        vector<int> prefix_score = vector<int>(multipath_aln.subpath_size());
        // previous subpath for traceback (we refer to subpaths by their index)
        vector<int> prev_subpath = vector<int>(multipath_aln.subpath_size());
        
        // identify source nodes if necessary
        if (multipath_aln.start_size() == 0) {
            identify_start_subpaths(multipath_aln);
        }
        
        // add sentinel for alignment start as base case
        for (int i = 0; i < multipath_aln.start_size(); i++) {
            prev_subpath[multipath_aln.start(i)] = -1;
        }
        
        int opt_score = 0;
        int opt_subpath = 0; we refer to subpaths by their index
        
        for (int i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath subpath& = multipath_aln.subpath(i);
            int32_t extended_score = prefix_score[i] + subpath.score();
            // are there outgoing subpaths?
            if (subpath.next_size() > 0) {
                for (int j = 0; j < subpath.next_size(); j++) {
                    int next = subpath.next(j);
                    // can we improve prefix score on following subpath through this one?
                    if (extended_score > prefix_score[next]) {
                        prev_subpath[next] = i;
                        prefix_score[next] = extended_score;
                    }
                }
            }
            else {
                // check if optimal alignment ends here
                if (extended_score > opt_score) {
                    opt_score = extended_score;
                    final_subpath = i;
                }
            }
        }
        
        // traceback the optimal subpaths until hitting sentinel (-1)
        vector<Path*> optimal_path_chunks;
        while (opt_subpath >= 0) {
            optimal_path_chunks.push_back(multipath_aln.subpath(opt_subpath).mutable_path());
            opt_subpath = prev_subpath[opt_subpath];
        }
        
        // merge the subpaths into one optimal path in the Alignment object
        aln_out.set_path(*(optimal_path_chunks[0]));
        Path* opt_path = aln_out.mutable_path();
        for (int i = 1; i < optimal_path_chunks.size(); i++) {
            Mapping* curr_end_mapping = opt_path->mutable_mapping(opt_path->mapping_size() - 1);
            
            // get the first mapping of the next path
            const Path& next_path = *(optimal_path_chunks[i]);
            const Mapping& next_start_mapping = next_path.mapping(0);
            
            // merge mappings if they occur on the same node and same strand
            if (curr_end_mapping->position().node_id() == next_start_mapping->position().node_id()
                && curr_end_mapping->position().is_reverse() == next_start_mapping->position().is_reverse()) {
                *curr_end_mapping = concat_mappings(*curr_end_mapping, *next_start_mapping);
            }
            else {
                opt_path.add_mapping(next_start_mapping);
            }
            
            // append the rest of the mappings
            for (int j = 1; j < next_path.mapping_size(); j++) {
                opt_path->add_mapping(next_path.mapping(j));
            }
        }
    }
    
    // stores the reverse complement of a Subpath in another Subpath
    //
    // note: this is not included in the header because reversing a subpath without going through
    // the multipath alignment can break invariants related to the edge lists
    //
    //  Args:
    //    subpath           subpath to reverse complement
    //    node_length       a function that returns the length of a node sequence from its node ID
    //    rev_comp_out      empty subpath to store reverse complement in (data will be overwritten
    //                      if not empty)
    //
    inline void rev_comp_subpath(const Subpath& subpath, const function<int64_t(int64_t)>& node_length,
                          Subpath& rev_comp_out) {
        
        rev_comp_out.set_path(reverse_complement_path(subpath.path(), node_length))
        rev_comp_out.set_score(subpath.score());
        // leave reversing the edges to the multipath alignment
    }
    
    void rev_comp_multipath_alignment(const MultipathAlignment& multipath_aln, const function<int64_t(int64_t)>& node_length,
                                      MultipathAlignment& rev_comp_out) {
        
        // reverse complement sequence
        rev_comp_out.set_sequence(reverse_complement(multipath_aln.sequence()));
        // reverse base qualities
        rev_comp_out.set_quality(string(multipath_aln.quality().rbegin(), multipath_aln.quality().rend()));
        
        vector<vector<int>> reverse_edge_lists = vector<vector<size_t>>(multipath_aln.subpath_size());
        vector<int> reverse_starts;
        
        // remove subpaths to avoid duplicating
        rev_comp_out.clear_subpath();
        
        // add subpaths in reverse order to maintain topological ordering
        for (int i = multipath_aln.subpath_size() - 1; i >= 0; i--) {
            const Subpath& subpath = multipath_aln.subpath(i);
            Subpath* rev_comp_subpath = rev_comp_out.add_subpath();
            rev_comp_subpath(subpath, *rev_comp_subpath);
            
            if (subpath.next_size() > 0) {
                // collect edges by their target (for reversing)
                for (int j = 0; j < subpath.next_size(); j++) {
                    reverse_edge_lists[subpath.next(j)].push_back(i);
                }
            }
            else {
                // sink subpaths become sources in reverse
                reverse_starts.push_back(i);
            }
        }
        
        // add reversed edges
        for (int i = 0; i < multipath_aln.subpath_size(); i++) {
            Subpath* rev_comp_subpath = rev_comp_out.mutable_subpath(i);
            vector<int>& reverse_edge_list = reverse_edge_lists[i]
            for (int j = 0; j < reverse_edge_list.size(); j++) {
                rev_comp_subpath->add_next(reverse_edge_list[j]);
            }
        }
        
        // remove start nodes that are invalid in reverse
        rev_comp_out.clear_start();
        
        // assume that if the original multipath alignment had its starts labeled they want them
        // labeled in the reverse complement too
        if (multipath_aln.start_size() > 0) {
            for (int i = 0; i < reverse_starts.size(); i++) {
                rev_comp_out.add_start(reverse_starts[i]);
            }
        }
    }
}










//
//  multipath_alignment.cpp
//

#include "multipath_alignment.hpp"

using namespace std;
namespace vg {
    
    void identify_start_subpaths(MultipathAlignment& multipath_aln) {
        
        // remove start nodes if there are any (avoids doubling them if this function is used liberally)
        multipath_aln.clear_start();
        
        // label nodes with incoming edges
        vector<bool> has_incoming_edge(multipath_aln.subpath_size(), false);
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            for (size_t j = 0; j < subpath.next_size(); j++) {
                has_incoming_edge[subpath.next(j)] = true;
            }
        }
        
        // construct list of starts
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            if (!has_incoming_edge[i]) {
                multipath_aln.add_start(i);
            }
        }
    }
    
    int32_t optimal_alignment_internal(const MultipathAlignment& multipath_aln, Alignment* aln_out) {
        // initialize DP structures
                
        // score of the optimal alignment ending in this subpath
        vector<int32_t> prefix_score(multipath_aln.subpath_size(), numeric_limits<int32_t>::min());
        // previous subpath for traceback (we refer to subpaths by their index)
        vector<int64_t> prev_subpath(multipath_aln.subpath_size());
        
        // add sentinel for alignment start as base case
        for (size_t i = 0; i < multipath_aln.start_size(); i++) {
            prev_subpath[multipath_aln.start(i)] = -1;
            prefix_score[multipath_aln.start(i)] = 0;
        }
        
        int32_t opt_score = 0;
        int64_t opt_subpath = -1; // we refer to subpaths by their index
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            int32_t extended_score = prefix_score[i] + subpath.score();
            // are there outgoing subpaths?
            if (subpath.next_size() > 0) {
                for (size_t j = 0; j < subpath.next_size(); j++) {
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
                    opt_subpath = i;
                }
            }
        }
        
        // are we constructing the alignment, or just getting the score?
        if (aln_out) {
            
            // traceback the optimal subpaths until hitting sentinel (-1)
            vector<const Path*> optimal_path_chunks;
            while (opt_subpath >= 0) {
                optimal_path_chunks.push_back(&(multipath_aln.subpath(opt_subpath).path()));
                opt_subpath = prev_subpath[opt_subpath];
            }
            
            // merge the subpaths into one optimal path in the Alignment object
            auto iter = optimal_path_chunks.rbegin();
            if (iter != optimal_path_chunks.rend()) {
                Path* opt_path = aln_out->mutable_path();
                *opt_path = *(*iter);
                iter++;
                int next_rank = opt_path->mapping_size() + 1;
                for (; iter != optimal_path_chunks.rend(); iter++) {
                    Mapping* curr_end_mapping = opt_path->mutable_mapping(opt_path->mapping_size() - 1);
                    
                    // get the first mapping of the next path
                    const Path& next_path = *(*iter);
                    const Mapping& next_start_mapping = next_path.mapping(0);
                    
                    size_t mapping_start_idx = 0;
                    // merge mappings if they occur on the same node and same strand
                    if (curr_end_mapping->position().node_id() == next_start_mapping.position().node_id()
                        && curr_end_mapping->position().is_reverse() == next_start_mapping.position().is_reverse()) {
                        
                        Edit* last_edit = curr_end_mapping->mutable_edit(curr_end_mapping->edit_size() - 1);
                        const Edit& first_edit = next_start_mapping.edit(0);
                        
                        // merge the first edit if it is the same type
                        size_t edit_start_idx = 0;
                        if ((last_edit->from_length() > 0) == (first_edit.from_length() > 0)
                            && (last_edit->to_length() > 0) == (first_edit.to_length() > 0)
                            && (last_edit->sequence().empty()) == (first_edit.sequence().empty())) {
                            
                            last_edit->set_from_length(last_edit->from_length() + first_edit.from_length());
                            last_edit->set_to_length(last_edit->to_length() + first_edit.to_length());
                            last_edit->set_sequence(last_edit->sequence() + first_edit.sequence());
                            
                            edit_start_idx++;
                            
                        }
                        
                        // append the rest of the edits
                        for (size_t j = edit_start_idx; j < next_start_mapping.edit_size(); j++) {
                            *curr_end_mapping->add_edit() = next_start_mapping.edit(j);
                        }
                        
                        mapping_start_idx++;
                    }
                    
                    // append the rest of the mappings
                    for (size_t j = mapping_start_idx; j < next_path.mapping_size(); j++) {
                        Mapping* next_mapping = opt_path->add_mapping();
                        *next_mapping = next_path.mapping(j);
                        next_mapping->set_rank(next_rank);
                        next_rank++;
                    }
                }
            }
        }
        
        return opt_score;
    }
    
    void optimal_alignment(const MultipathAlignment& multipath_aln, Alignment& aln_out) {
        
        // transfer read information over to alignment
        transfer_read_metadata(multipath_aln, aln_out);
        aln_out.set_mapping_quality(multipath_aln.mapping_quality());
        
        // do dynamic programming and traceback the optimal alignment
        int32_t score = optimal_alignment_internal(multipath_aln, &aln_out);
        
        aln_out.set_score(score);
    }
    
    
    int32_t optimal_alignment_score(const MultipathAlignment& multipath_aln){
        // do dynamic programming without traceback
        return optimal_alignment_internal(multipath_aln, nullptr);
    }
    
    /// Stores the reverse complement of a Subpath in another Subpath
    ///
    /// note: this is not included in the header because reversing a subpath without going through
    /// the multipath alignment can break invariants related to the edge lists
    ///
    ///  Args:
    ///    subpath           subpath to reverse complement
    ///    node_length       a function that returns the length of a node sequence from its node ID
    ///    rev_comp_out      empty subpath to store reverse complement in (data will be overwritten
    ///                      if not empty)
    ///
    inline void rev_comp_subpath(const Subpath& subpath, const function<int64_t(int64_t)>& node_length,
                          Subpath& rev_comp_out) {
        
        *(rev_comp_out.mutable_path()) = reverse_complement_path(subpath.path(), node_length);
        rev_comp_out.set_score(subpath.score());
        // leave reversing the edges to the multipath alignment
    }
    
    void rev_comp_multipath_alignment(const MultipathAlignment& multipath_aln, const function<int64_t(int64_t)>& node_length,
                                      MultipathAlignment& rev_comp_out) {
        
        // reverse complement sequence
        rev_comp_out.set_sequence(reverse_complement(multipath_aln.sequence()));
        // reverse base qualities
        rev_comp_out.set_quality(string(multipath_aln.quality().rbegin(), multipath_aln.quality().rend()));
        
        vector< vector<size_t> > reverse_edge_lists(multipath_aln.subpath_size());
        vector<size_t> reverse_starts;
        
        // remove subpaths to avoid duplicating
        rev_comp_out.clear_subpath();
        
        // add subpaths in reverse order to maintain topological ordering
        for (int64_t i = multipath_aln.subpath_size() - 1; i >= 0; i--) {
            const Subpath& subpath = multipath_aln.subpath(i);
            Subpath* rc_subpath = rev_comp_out.add_subpath();
            rev_comp_subpath(subpath, node_length, *rc_subpath);
            
            if (subpath.next_size() > 0) {
                // collect edges by their target (for reversing)
                for (size_t j = 0; j < subpath.next_size(); j++) {
                    reverse_edge_lists[subpath.next(j)].push_back(i);
                }
            }
            else {
                // sink subpaths become sources in reverse
                reverse_starts.push_back(i);
            }
        }
        
        // add reversed edges
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            Subpath* rc_subpath = rev_comp_out.mutable_subpath(i);
            vector<size_t>& reverse_edge_list = reverse_edge_lists[multipath_aln.subpath_size() - i - 1];
            for (size_t j = 0; j < reverse_edge_list.size(); j++) {
                rc_subpath->add_next(multipath_aln.subpath_size() - reverse_edge_list[j] - 1);
            }
        }
        
        // remove start nodes that are invalid in reverse
        rev_comp_out.clear_start();
        
        // assume that if the original multipath alignment had its starts labeled they want them
        // labeled in the reverse complement too
        if (multipath_aln.start_size() > 0) {
            for (size_t i = 0; i < reverse_starts.size(); i++) {
                rev_comp_out.add_start(multipath_aln.subpath_size() - reverse_starts[i] - 1);
            }
        }
    }
    
    void rev_comp_multipath_alignment_in_place(MultipathAlignment* multipath_aln,
                                               const function<int64_t(int64_t)>& node_length) {
        
        // reverse complement sequence
        reverse_complement_in_place(*multipath_aln->mutable_sequence());
        // reverse base qualities
        string* quality = multipath_aln->mutable_quality();
        std::reverse(quality->begin(), quality->end());
        
        // current reverse edges
        vector< vector<int64_t> > reverse_edge_lists(multipath_aln->subpath_size());
        // current sink nodes (will be starts)
        vector<int64_t> reverse_starts;
        
        int64_t subpath_swap_size = multipath_aln->subpath_size() / 2;
        int64_t last = multipath_aln->subpath_size() - 1;
        for (int64_t i = 0, j = last; i < subpath_swap_size; i++, j--) {
            Subpath* subpath_1 = multipath_aln->mutable_subpath(i);
            Subpath* subpath_2 = multipath_aln->mutable_subpath(j);
            
            // add reverse edges for first subpath
            if (subpath_1->next_size() > 0) {
                for (int64_t k = 0; k < subpath_1->next_size(); k++) {
                    reverse_edge_lists[subpath_1->next(k)].push_back(i);
                }
            }
            else {
                reverse_starts.push_back(i);
            }
            
            // add reverse edges for second subpath
            if (subpath_2->next_size() > 0) {
                for (int64_t k = 0; k < subpath_2->next_size(); k++) {
                    reverse_edge_lists[subpath_2->next(k)].push_back(j);
                }
            }
            else {
                reverse_starts.push_back(j);
            }
            
            // clear current edges
            subpath_1->clear_next();
            subpath_2->clear_next();
            
            // reverse complement the paths
            reverse_complement_path_in_place(subpath_1->mutable_path(), node_length);
            reverse_complement_path_in_place(subpath_2->mutable_path(), node_length);
            
            // swap their positions (to maintain topological ordering)
            std::swap(*subpath_1, *subpath_2);
        }
        
        // repeat process for the middle subpath if there is an odd number
        if (multipath_aln->subpath_size() % 2) {
            Subpath* subpath = multipath_aln->mutable_subpath(subpath_swap_size);
            if (subpath->next_size() > 0) {
                for (int64_t k = 0; k < subpath->next_size(); k++) {
                    reverse_edge_lists[subpath->next(k)].push_back(subpath_swap_size);
                }
            }
            else {
                reverse_starts.push_back(subpath_swap_size);
            }
            
            subpath->clear_next();
            reverse_complement_path_in_place(subpath->mutable_path(), node_length);
        }
        
        // add reversed edges
        for (int64_t i = 0, j = last; i < multipath_aln->subpath_size(); i++, j--) {
            vector<int64_t> edges = reverse_edge_lists[j];
            Subpath* subpath = multipath_aln->mutable_subpath(i);
            for (int64_t k : edges) {
                subpath->add_next(last - k);
            }
        }
        
        // if we had starts labeled before, label them again
        if (multipath_aln->start_size() > 0) {
            multipath_aln->clear_start();
            for (int64_t i : reverse_starts) {
                multipath_aln->add_start(last - i);
            }
        }
    }
    
    void to_multipath_alignment(const Alignment& aln, MultipathAlignment& multipath_aln_out) {
        
        // clear repeated fields
        multipath_aln_out.clear_subpath();
        multipath_aln_out.clear_start();
        
        // transfer read and alignment metadata
        transfer_read_metadata(aln, multipath_aln_out);
        multipath_aln_out.set_mapping_quality(aln.mapping_quality());
        
        // transfer alignment and score
        if (aln.has_path() || aln.score()) {
            Subpath* subpath = multipath_aln_out.add_subpath();
            subpath->set_score(aln.score());
            *(subpath->mutable_path()) = aln.path();
        }
        
    }
    
    void transfer_read_metadata(const Alignment& from, MultipathAlignment& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        to.set_read_group(from.read_group());
        to.set_name(from.name());
        to.set_sample_name(from.sample_name());
        
        // no difference in these fields for MultipathAlignments
        if (from.has_fragment_prev()) {
            to.set_paired_read_name(from.fragment_prev().name());
        }
        else {
            to.set_paired_read_name(from.fragment_next().name());
        }
    }
    
    void transfer_read_metadata(const MultipathAlignment& from, Alignment& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        to.set_read_group(from.read_group());
        to.set_name(from.name());
        to.set_sample_name(from.sample_name());
        
        // note: not transferring paired_read_name because it is unclear whether
        // it should go into fragment_prev or fragment_next
    }
}










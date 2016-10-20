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
    
    void optimal_alignment(const MultipathAlignment& multipath_aln, Alignment& aln_out) {
        
        // transfer read information over to alignment
        aln_out.set_sequence(multipath_aln.sequence());
        aln_out.set_quality(multipath_aln.quality());
        aln_out.set_name(multipath_aln.name());
        aln_out.set_sample_name(multipath_aln.sample_name());
        aln_out.set_read_group(multipath_aln.read_group());
        aln_out.set_mapping_quality(multipath_aln.mapping_quality());
        
        // initialize DP structures
        
        // score of the optimal alignment ending in this subpath
        vector<int> prefix_score = vector<int>(multipath_aln.subpath_size());
        // previous subpath for traceback (we refer to subpaths by their index)
        vector<int> prev_subpath = vector<int>(multipath_aln.subpath_size());
        
        // add sentinel for alignment start as base case
        for (int i = 0; i < multipath_aln.start_size(); i++) {
            prev_subpath[multipath_aln.start(i)] = -1;
        }
        
        int opt_score = 0;
        int opt_subpath = 0; // we refer to subpaths by their index
        
        for (int i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
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
                    opt_subpath = i;
                }
            }
        }
        
        // traceback the optimal subpaths until hitting sentinel (-1)
        vector<const Path*> optimal_path_chunks;
        while (opt_subpath >= 0) {
            optimal_path_chunks.push_back(&(multipath_aln.subpath(opt_subpath).path()));
            opt_subpath = prev_subpath[opt_subpath];
        }
        
        // merge the subpaths into one optimal path in the Alignment object
        auto iter = optimal_path_chunks.rbegin();
        Path* opt_path = aln_out.mutable_path();
        *opt_path = *(*iter);
        iter++;
        for (int i = 0; iter != optimal_path_chunks.rend(); iter++, i++) {
            Mapping* curr_end_mapping = opt_path->mutable_mapping(opt_path->mapping_size() - 1);
            
            // get the first mapping of the next path
            const Path& next_path = *(*iter);
            const Mapping& next_start_mapping = next_path.mapping(0);
            
            // merge mappings if they occur on the same node and same strand
            if (curr_end_mapping->position().node_id() == next_start_mapping.position().node_id()
                && curr_end_mapping->position().is_reverse() == next_start_mapping.position().is_reverse()) {
                *curr_end_mapping = concat_mappings(*curr_end_mapping, next_start_mapping);
            }
            else {
                *(opt_path->add_mapping()) = next_start_mapping;
            }
            
            // append the rest of the mappings
            for (int j = 1; j < next_path.mapping_size(); j++) {
                *(opt_path->add_mapping()) = next_path.mapping(j);
            }
        }
        
        aln_out.set_score(opt_score);
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
        
        vector< vector<int> > reverse_edge_lists = vector< vector<int> >(multipath_aln.subpath_size());
        vector<int> reverse_starts;
        
        // remove subpaths to avoid duplicating
        rev_comp_out.clear_subpath();
        
        // add subpaths in reverse order to maintain topological ordering
        for (int i = multipath_aln.subpath_size() - 1; i >= 0; i--) {
            const Subpath& subpath = multipath_aln.subpath(i);
            Subpath* rc_subpath = rev_comp_out.add_subpath();
            rev_comp_subpath(subpath, node_length, *rc_subpath);
            
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
            Subpath* rc_subpath = rev_comp_out.mutable_subpath(i);
            vector<int>& reverse_edge_list = reverse_edge_lists[i];
            for (int j = 0; j < reverse_edge_list.size(); j++) {
                rc_subpath->add_next(reverse_edge_list[j]);
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
    
    void to_multipath_alignment(const Alignment& aln, MultipathAlignment& multipath_aln_out) {
        
        // clear repeated fields
        multipath_aln_out.clear_subpath();
        multipath_aln_out.clear_start();
        
        // transfer read and alignment metadata
        multipath_aln_out.set_sequence(aln.sequence());
        multipath_aln_out.set_quality(aln.quality());
        multipath_aln_out.set_read_group(aln.read_group());
        multipath_aln_out.set_name(aln.name());
        multipath_aln_out.set_sample_name(aln.sample_name());
        multipath_aln_out.set_mapping_quality(aln.mapping_quality());
        
        // transfer alignment and score
        Subpath* subpath = multipath_aln_out.add_subpath();
        subpath->set_score(aln.score());
        *(subpath->mutable_path()) = aln.path();
        
    }
}










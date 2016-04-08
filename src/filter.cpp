#include "filter.hpp"

using namespace std;

namespace vg{

    Filter::Filter(){

    }

    Filter::~Filter(){

    }

    void Filter::set_min_depth(int depth){
        min_depth = depth;
    }
    void Filter::set_min_qual(int qual){
        min_qual = qual;
    }

    void Filter::set_min_percent_identity(double pct_id){
        min_percent_identity = pct_id;
    }

    void Filter::set_avg_qual(double avg_qual){
        min_avg_qual = avg_qual;
    }

    void Filter::set_filter_matches(bool fm){
        filter_matches = fm;
    }
    void Filter::set_remove_failing_alignments(bool fm){
        remove_failing_alignments = fm;
    }

    Alignment Filter::path_divergence_filter(Alignment& aln){
        
    }

    void Filter::set_softclip_filter(bool dsc){
        do_softclip = dsc;        
    }

    void Filter::set_splitread_filter(bool dsr){
        do_splitread = dsr;
    }

    /**
     * CLI: vg filter -d 10 -q 40 -r -R
     * -r: track depth of both novel variants and those in the graph.
     * -R: remove edits that fail the filter (otherwise toss the whole alignment)
     */
    Alignment Filter::depth_filter(Alignment& aln){
        Path path = aln.path();
        //TODO handle reversing mappings

        for (int i = 0; i < path.mapping_size(); i++){
            Mapping mapping = path.mapping(i);
            Position start_pos = mapping.position();
            int64_t start_node = start_pos.node_id();
            int64_t start_offset = start_pos.offset();
            int64_t curr_offset_in_graph = 0;
            int64_t curr_offset_in_alignment = 0;
            stringstream pst;

            pst << start_node << "_" << curr_offset_in_graph;
            string p_hash = pst.str();
            for (int j = 0; j < mapping.edit_size(); j++){
                Edit ee = mapping.edit(j);
                if (ee.from_length() == ee.to_length() && ee.sequence() == ""){
                    if (!filter_matches){
                        continue;
                    }
                }
                stringstream est;
                est <<  ee.from_length() << "_" << ee.to_length() << "_" + ee.sequence();
                string e_hash = est.str();
                pos_to_edit_to_depth[p_hash][e_hash] += 1;
                /**
                 * If an edit fails the filter, either return a new empty alignment
                 * OR
                 * return a new alignment identical to the old one EXCEPT where
                 * the offending edit has been replaced by a match to the reference.
                 */
                if (pos_to_edit_to_depth[p_hash][e_hash] < min_depth){
                    if (remove_failing_alignments){
                        return Alignment();
                    }

                    else {
                        Alignment edited_aln = Alignment(aln);
                        edited_aln.mutable_path()->mutable_mapping(i)->mutable_edit(j)->set_sequence("");
                        edited_aln.mutable_path()->mutable_mapping(i)->mutable_edit(j)->set_from_length(ee.from_length());
                        edited_aln.mutable_path()->mutable_mapping(i)->mutable_edit(j)->set_to_length(ee.from_length());
                        return edited_aln;
                    }
                }
            }
            return aln;
        }


    }
    Alignment Filter::qual_filter(Alignment& aln){

    }
    Alignment Filter::coverage_filter(Alignment& aln){

    }
    Alignment Filter::avg_qual_filter(Alignment& aln){
        double total_qual = 0.0;

        cerr << "UNTESTED" << endl;
        exit(1);

        std::function<double(int64_t, int64_t)> calc_avg_qual = [](int64_t total_qual, int64_t length){
            return ((double) total_qual / (double) length);
        };
        
    

        Path path = aln.path();
        //TODO: handle reversing alignments
        
        for (int i = 0; i < aln.quality().size(); i++){
            cerr << aln.quality()[i] << endl;
        }
        if (calc_avg_qual(total_qual, aln.sequence().size()) < min_avg_qual){
            return Alignment();
        }

        return aln;

    }


    Alignment Filter::soft_clip_filter(Alignment& aln){

    }
    /**
     * Split reads map to two separate paths OR vastly separated non-consecutive
     * nodes in a single path.
     *
     * They're super important for detecting structural variants, so we may want to
     * filter them out or collect only split reads.
     */
    Alignment Filter::split_read_filter(Alignment& aln){

        //TODO binary search for breakpoint in read.

    }
    /**
     * Filter reads that are less than <PCTID> reference.
     * I.E. if a read matches the reference 80% alogn 80% of its
     * length, and your cutoff is 90% PCTID, throw it out.
     */
    Alignment Filter::percent_identity_filter(Alignment& aln){
        double read_pctid = 0.0;
        //read pct_id = len(matching sequence / len(total sequence)

        int64_t aln_total_len = aln.sequence().size();
        int64_t aln_match_len = 0;

        std::function<double(int64_t, int64_t)> calc_pct_id = [](int64_t rp, int64_t ttlp){
            return ((double) rp / (double) ttlp);
        };



        Path path = aln.path();
        //TODO handle reversing mappings

        for (int i = 0; i < path.mapping_size(); i++){
            Mapping mapping = path.mapping(i);

            for (int j = 0; j < mapping.edit_size(); j++){
                Edit ee = mapping.edit(j);
                if (ee.from_length() == ee.to_length() && ee.sequence() == ""){
                    aln_match_len += ee.to_length();        
                }

            }    
        }
        if (calc_pct_id(aln_match_len, aln_total_len) < min_percent_identity){
            return Alignment();
        }
        
        return aln;
        

    }
}

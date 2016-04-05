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

    void Filter::set_min_pct_identity(double pct_id){
        min_pct_identity = pct_id;
    }

    void Filter::set_avg_qual(double avg_qual){
        min_avg_qual = avg_qual;
    }
    /**
     * CLI: vg filter -d 10 -q 40 -r -R
     * -r: track depth of both novel variants and those in the graph.
     * -R: remove edits that fail the filter (otherwise toss the whole alignment)
     */
    Alignment Filter::depth_filter(Alignment& aln){
        Path path = aln.path();
        bool depth_filter_matches = false;
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
                    if (depth_filter_matches){

                    }
                    else {
                        continue;
                    }
                }
                stringstream est;
                est <<  ee.from_length() << "_" << ee.to_length() << "_" + ee.sequence();
                string e_hash = est.str();
                pos_to_edit_to_depth[p_hash][e_hash] += 1;
                if (pos_to_edit_to_depth[p_hash][e_hash] < min_depth){
                    if (remove_failing_edits){
                        return Alignment();
                    }

                    else {
                        Alignment edited_aln = Alignment();
                        //Alignment.set_sequence(aln.sequence());
                        //Edit new_ed = Edit();
                        //Alignment.path().mapping(i).edit(j) = new_ed;
                        return edited_aln;
                    }
                }
            }
            return aln;
        }


    }
    Alignment& Filter::qual_filter(Alignment& aln){

    }
    Alignment& Filter::coverage_filter(Alignment& aln){

    }
    Alignment& Filter::avg_qual_filter(Alignment& aln){

    }
    Alignment& Filter::soft_clip_filter(Alignment& aln){

    }
    Alignment& Filter::percent_identity_filter(Alignment& aln){
        double pctid = 0.0;
        int64_t graph_len = my_vg->total_length_of_nodes();

    }
}

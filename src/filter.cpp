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

    void Filter::set_avg(bool is_avg){
        use_avg = is_avg;
    }

    void Filter::set_filter_matches(bool fm){
        filter_matches = fm;
    }
    void Filter::set_remove_failing_edits(bool fm){
        remove_failing_edits = fm;
    }

    void Filter::set_soft_clip_limit(int max_clip){
        soft_clip_limit = max_clip;
    }

    void Filter::set_split_read_limit(int sr){
        split_read_limit = sr;
    }

    void Filter::set_window_length(int wind_len){
        window_length = wind_len;
    }

    void Filter::set_my_vg(vg::VG* vg){
        my_vg = vg;
    }

    void Filter::set_my_xg_idx(xg::XG* idx){
        my_xg_index = idx;
    }

    void Filter::set_inverse(bool do_inv){
        inverse = do_inv;
    }

    Alignment Filter::depth_filter(Alignment& aln){
        if (use_avg && window_length != 0){

        }
        else if (use_avg != 0){

        }
        else{

        }

        Path path = aln.path();
        //TODO handle reversing mappings
        vector<int>* qual_window;
        if (window_length > 0){
            qual_window = new vector<int>();
        }

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
#pragma omp critical(write)
                pos_to_edit_to_depth[p_hash][e_hash] += 1;
                /**
                 * If an edit fails the filter, either return a new empty alignment
                 * OR
                 * return a new alignment identical to the old one EXCEPT where
                 * the offending edit has been replaced by a match to the reference.
                 */
                if (pos_to_edit_to_depth[p_hash][e_hash] < min_depth){
                    if (!remove_failing_edits){
                        return inverse ? aln : Alignment();
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
            return inverse ? Alignment() : aln;
        }


    }

    Alignment Filter::path_length_filter(Alignment& aln){
        for (int i = 0; i < aln.fragment_size(); i++){
            Path t = aln.fragment(i);
            if (t.length() > this->max_path_length){
                return inverse ? Alignment() : aln;
            }
            else{
                return inverse ? aln : Alignment();
            }
        }

        return inverse ? aln : Alignment();

    }

    std::pair<Alignment, Alignment> Filter::path_length_filter(Alignment& aln_first, Alignment& aln_second){
        Alignment x = path_length_filter(aln_first);
        Alignment y = path_length_filter(aln_second);
        if (x.name().empty() || y.name().empty()){
            return inverse ? make_pair(x, y) : make_pair(Alignment(), Alignment());
        }
        else{
            return  inverse ? make_pair(Alignment(), Alignment()) : make_pair(x, y);
        }
    }




    pair<Alignment, Alignment> Filter::depth_filter(Alignment& aln_first, Alignment& aln_second){
        aln_first = depth_filter(aln_first);
        aln_second = depth_filter(aln_second);
        if (!(aln_first.name() == "") && !(aln_first.name() == "")){
            return inverse ? make_pair(aln_first, aln_second) : make_pair(Alignment(), Alignment());
        }
        else{
            return inverse ? make_pair(Alignment(), Alignment()) : make_pair(aln_first, aln_second);
        }
    }
    pair<Alignment, Alignment> qual_filter(Alignment& aln_first, Alignment& aln_second){


        if (aln_first.name() == "" || aln_first.name() == ""){
            return make_pair(Alignment(), Alignment());
        }
        else{
            return make_pair(aln_first, aln_second);
        }

    }
    pair<Alignment, Alignment> percent_identity_filter(Alignment& aln_first, Alignment& aln_second){


        if (aln_first.name() == "" || aln_first.name() == ""){
            return make_pair(Alignment(), Alignment());
        }
        else{
            return make_pair(aln_first, aln_second);
        }

    }
    pair<Alignment, Alignment> Filter::soft_clip_filter(Alignment& aln_first, Alignment& aln_second){
        Alignment a_check = soft_clip_filter(aln_first);
        Alignment b_check = soft_clip_filter(aln_second);

        if (a_check.name() == "" || b_check.name() == ""){
            return inverse ? make_pair(Alignment(), Alignment()) : make_pair(aln_first, aln_second) ;
        }
        else{
            return inverse ? make_pair(aln_first, aln_second) : make_pair(Alignment(), Alignment()) ;
        }

    }
    pair<Alignment, Alignment> split_read_filter(Alignment& aln_first, Alignment& aln_second){

    }
    pair<Alignment, Alignment> path_divergence_filter(Alignment& aln_first, Alignment& aln_second){

    }
    pair<Alignment, Alignment> reversing_filter(Alignment& aln, Alignment& aln_second){

    }


    /* PE functions using fragment_prev and fragment_next */
    Alignment Filter::one_end_anchored_filter(Alignment& aln){
        if (aln.fragment_prev().name() != ""){
            if (aln.path().name() == "" || aln.fragment_prev().path().name() == ""){
                inverse ? Alignment() : aln;
            }
            else{
                inverse ? aln : Alignment();
            }
        }
        else{
            return inverse ? aln : Alignment();
        }
    }

    Alignment Filter::interchromosomal_filter(Alignment& aln){
        bool fails = aln.path().name() != aln.fragment_prev().path().name();
        if (fails){
            return inverse ? Alignment() : aln;
        }
        else{
            return inverse ? aln : Alignment();
        }
    }

    Alignment Filter::insert_size_filter(Alignment& aln){
        // Get mapping distance from XG index
        // between last node in aln and first node in aln.fragment_prev

    }

    Alignment Filter::orientation_filter(Alignment& aln){
        bool f_rev = false;
        bool s_rev = false;
        Path f_path = aln.path();
        Path s_path = aln.fragment_prev().path();
        for (int i = 0; i < f_path.mapping_size(); i++){
            if (f_path.mapping(i).position().is_reverse()){
                f_rev = true;
            }
        }

        for (int j = 0; j < s_path.mapping_size(); j++){
            if (s_path.mapping(j).position().is_reverse()){
                s_rev = true;
            }
        }

        if (f_rev & s_rev){
            return inverse ? Alignment() : aln;
        }
        else{
            return inverse ? aln : Alignment();
        }
    }



    /*PE Functions*/
    pair<Alignment, Alignment> Filter::one_end_anchored_filter(Alignment& aln_first, Alignment& aln_second){
        if (aln_first.mapping_quality() == 0 | aln_second.mapping_quality() == 0){
            return inverse ? std::make_pair(Alignment(), Alignment()) : std::make_pair(aln_first, aln_second);
        }
        else{
            return inverse ? std::make_pair(aln_first, aln_second) : std::make_pair(Alignment(), Alignment());
        }
    }

    pair<Alignment, Alignment> Filter::interchromosomal_filter(Alignment& aln_first, Alignment& aln_second){
        if (aln_first.path().name() != aln_second.path().name()){
            return std::make_pair(aln_first, aln_second);
        }
        else{
            return std::make_pair(Alignment(), Alignment());
        }
    }

    pair<Alignment, Alignment> Filter::insert_size_filter(Alignment& aln_first, Alignment& aln_second){
        // TODO: gret positions from aln_first and aln_second
        int distance = my_xg_index->approx_path_distance(aln_first.path().name(), 1, 1);
        if (distance > my_max_distance){
            return std::make_pair(aln_first, aln_second);
        }
        else{
            return std::make_pair(Alignment(), Alignment());
        }
    }

    pair<Alignment, Alignment> Filter::orientation_filter(Alignment& aln_first, Alignment& aln_second){

        bool f_rev = false;
        bool s_rev = false;
        Path f_path = aln_first.path();
        Path s_path = aln_second.path();
        for (int i = 0; i < f_path.mapping_size(); i++){
            if (f_path.mapping(i).position().is_reverse()){
                f_rev = true;
            }
        }

        for (int j = 0; j < s_path.mapping_size(); j++){
            if (s_path.mapping(j).position().is_reverse()){
                s_rev = true;
            }
        }



        if (!s_rev != !f_rev){
            return inverse ? std::make_pair(aln_first, aln_second) : std::make_pair(Alignment(), Alignment());
        }
        else{
            return inverse ? std::make_pair(Alignment(), Alignment()) : std::make_pair(aln_first, aln_second);
        }

    }

    pair<Alignment, Alignment> Filter::deletion_filter(Alignment& aln_first, Alignment& aln_second){
        // path_length, split read
        // REFINE USING SOFT CLIPS
        pair<Alignment, Alignment> ret_alns = path_length_filter(aln_first, aln_second);
        //pair<Alignment, Alignment> x_alns = split_read_filter(ret_alns.first, ret_alns.second);
        
        bool found = false;
        if (inverse | (ret_alns.first.name() != "" || ret_alns.second.name() != "")){
            found = true;
        }

        
        if (found){
            return ret_alns;
        }
        else {
            return make_pair(Alignment(), Alignment());
        }

    }

    /**
    * Filters insertion-characterizing reads based upon 
    * discordant pairs (fragment length too short)
    * one-end anchored / softclipped portions
    * and read depth, one day
    */
    pair<Locus, Locus> Filter::insertion_filter(Alignment& aln_first, Alignment& aln_second){
        
    }

    pair<Locus, Locus> Filter::duplication_filter(Alignment& aln_first, Alignment& aln_second){

    }

    pair<Locus, Locus> Filter::inversion_filter(Alignment& aln_first, Alignment& aln_second){

    }
    pair<Locus, Locus> Filter::breakend_filter(Alignment& aln_first, Alignment& aln_second){

    }



    /**
     *
     * Looks for alignments that transition from one path to another
     * over their length. This may occur for one of several reasons:
     * 1. The read covers a translocation
     * 2. The read looks a lot like two different (but highly-similar paths)
     * 3. The read is shattered (e.g. as in chromothripsis)
     *
     * Default behavior: if the Alignment is path divergent, return an empty Alignment, else return aln
     * Inverse behavior: if the Alignment is path divergent, return aln, else return an empty Alignment
     */
    Alignment Filter::path_divergence_filter(Alignment& aln){
        Path path = aln.path();
        for (int i = 1; i < path.mapping_size(); i++){
            Mapping mapping = path.mapping(i);
            Position pos = mapping.position();
            id_t current_node = pos.node_id();
            id_t prev_node = path.mapping(i - 1).position().node_id();
            bool paths_match = false;
            vector<size_t> paths_of_prev = my_xg_index->paths_of_node(prev_node);
            for (int i = 0; i < paths_of_prev.size(); i++){
                string p_name = my_xg_index->path_name(paths_of_prev[i]);
                if (my_xg_index->path_contains_node(p_name, current_node)){
                    paths_match = true;
                }
            }
            if (!paths_match){
                return inverse ? aln : Alignment();
            }

        }
        return inverse ? Alignment() : aln;
    }


    /**
     * Looks for alignments that change direction over their length.
     * This may happen because of:
     * 1. Mapping artifacts
     * 2. Cycles
     * 3. Highly repetitive regions
     * 4. Inversions (if you're lucky enough)
     *
     * Default behavior: if the Alignment reverses, return an empty Alignment.
     * inverse behavior: if the Alignment reverses, return the Alignment.
     */
    Alignment Filter::reversing_filter(Alignment& aln){

        Path path = aln.path();
        bool prev = false;

        for (int i = 1; i < path.mapping_size(); i++){
            Mapping mapping = path.mapping(i);
            Position pos = mapping.position();
            bool prev = path.mapping(i - 1).position().is_reverse();
            if (prev != pos.is_reverse()){
                return inverse ? aln : Alignment();
            }

        }
        return inverse ? Alignment() : aln;

    }

    /**
     * Looks for Alignments that have large overhangs at the end of them.
     *
     * Default behavior: if an alignment has a right- or left- clip that is longer
     * than the maximum allowed, return an empty alignment.
     *
     * Inverse Behavior: if the alignment has a clip that is larger than the 
     * maximum allowed at either end, return the alignment.
     */


    /**
     * CLI: vg filter -d 10 -q 40 -r -R
     * -r: track depth of both novel variants and those in the graph.
     * -R: remove edits that fail the filter (otherwise toss the whole alignment)
     */
    Alignment Filter::qual_filter(Alignment& aln){
        string quals = aln.quality();
        int offset = qual_offset;
        for (int i = 0; i < quals.size(); i++){
            if (((int) quals[i] - offset) < min_qual){
                if (!remove_failing_edits){
                    return inverse ? aln : Alignment();
                }
                else{
                    //TODO: Should we really edit out poor-quality bases??
                    // It's complex and awkward.
                    return inverse ? aln : Alignment();
                }
            }
        }

        return inverse ? Alignment() : aln;
    }


    Alignment Filter::coverage_filter(Alignment& aln){

    }
    Alignment Filter::avg_qual_filter(Alignment& aln){
        double total_qual = 0.0;
        // If the parameter for window size is zero, set the local equivalent
        // to the entire size of the qual scores.
        int win_len = window_length; //>= 1 ? window_length : aln.quality().size();

        cerr << "UNTESTED" << endl;
        exit(1);

        std::function<double(int64_t, int64_t)> calc_avg_qual = [](int64_t total_qual, int64_t length){
            return ((double) total_qual / (double) length);
        };

        /**
         * Helper function.
         * Sums qual scores from start : start + window_length
         * if start + window_len > qualstr.size(), return a negative float
         * otherwise return the average quality within the window
         *
         */
        std::function<double(string, int, int)> window_qual = [&](string qualstr, int start, int window_length){
            if (start + window_length > qualstr.size()){
                return -1.0;
            }

            total_qual = 0.0;
            for (int i = start; i < start + window_length; i++){
                total_qual += (int) qualstr[i];
            }

            return calc_avg_qual(total_qual, window_length);
        };



        //TODO: handle reversing alignments
        string quals = aln.quality();
        for (int i = 0; i < quals.size() - win_len; i++){
            double w_qual = window_qual(quals, i, win_len);
            if ( w_qual < 0.0){
                break;
            }
            if (w_qual < min_avg_qual){
                return inverse ? aln : Alignment();
            }
        }

        return inverse ? Alignment() : aln;

    }


    Alignment Filter::soft_clip_filter(Alignment& aln){
        //Find overhangs - portions of the read that
        // are inserted at the ends.
        if (aln.path().mapping_size() > 0){
            Path path = aln.path();
            Edit left_edit = path.mapping(0).edit(0);
            Edit right_edit = path.mapping(path.mapping_size() - 1).edit(path.mapping(path.mapping_size() - 1).edit_size() - 1);
            int left_overhang = left_edit.to_length() - left_edit.from_length();
            int right_overhang = right_edit.to_length() - right_edit.from_length();
            if (left_overhang > soft_clip_limit || right_overhang > soft_clip_limit){
                return inverse ? Alignment() : aln;
            }
            else{
                return inverse ?  aln : Alignment();
            }
        }
        else{
            if (aln.sequence().length() > soft_clip_limit){
                return inverse ? Alignment() : aln;
            }
            cerr << "WARNING: SHORT ALIGNMENT: " << aln.sequence().size() << "bp" << endl
                << "WITH NO MAPPINGS TO REFERENCE" << endl
                << "CONSIDER REMOVING IT FROM ANALYSIS" << endl;
            return inverse ? Alignment() : aln;
        }

    }
    /**
     * Split reads map to two separate paths in the graph OR vastly separated non-consecutive
     * nodes in a single path.
     *
     * They're super important for detecting structural variants, so we may want to
     * filter them out or collect only split reads.
     */
    Alignment Filter::split_read_filter(Alignment& aln){

        //TODO binary search for breakpoint in read would be awesome.
        Path path = aln.path();
        //check if nodes are on same path(s)

        int top_side = path.mapping_size() - 1;
        int bottom_side = 0;

        Mapping bottom_mapping;
        Mapping top_mapping;

        string main_path = "";
        while (top_side > bottom_side){
            //main_path = path_of_node(path.mapping(bottom_side);
            //
            //Check if paths are different
            //if (divergent(node1, node2){
            //    return inverse ? aln : Alignment();
            //}
            top_mapping = path.mapping(top_side);
            bottom_mapping = path.mapping(bottom_side);
            Position top_pos = top_mapping.position();
            Position bot_pos = bottom_mapping.position();
            id_t top_id = top_pos.node_id();
            id_t bottom_id = bot_pos.node_id();

            // TODO USE THE XG
            if (abs(top_id - bottom_id) > 10){
                return inverse ? aln : Alignment();
            }

            // Check if two mappings are far apart 
            //
            // Check if a single mapping has a huge indel




            top_side--;
            bottom_side++;
        }

        return inverse ? Alignment() : aln;

    }


    /**
     * Filter reads that are less than <PCTID> reference.
     * I.E. if a read matches the reference along 80% of its
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
            return inverse ? aln : Alignment();
        }

        return inverse ? Alignment() : aln;


    }
}

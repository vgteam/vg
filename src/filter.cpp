#include "filter.hpp"

using namespace std;

namespace vg{

    Filter::Filter(){

    }

    Filter::~Filter(){

    }


    // to expand to multiple paths, we'll need to maintain a map of maps
    // pathname -> map<node_id, int> >
    int64_t Filter::distance_between_positions(Position first, Position second){
        int64_t f_node_pos = node_to_position[first.node_id()];
        int64_t r_node_pos = node_to_position[second.node_id()];
        int64_t fp = f_node_pos + first.offset();
        int64_t rp = r_node_pos + second.offset();
        return abs(rp - fp);
    }

    pair<Alignment, int> Filter::refactor_split_alignment(Alignment& a){
        // Get the clipped portion
        string clipseq = get_clipped_seq(a);
        // Remap it
        // TODO switch out for Banded Aligner to speed things up
        vector<Alignment> remaps = remap(clipseq);
        if (remaps.size() < 1){
            return make_pair(a, 0);
        }
        Alignment ret = a;
        // Tack the new split's path onto the anchored bit, depending on which side it came from.
        bool l_clipped = is_left_clipped(a);
        Path anchor = trim_hanging_ends(a.path());
        Path remap_path = remaps[0].path();
        Path* new_path = ret.mutable_path();
        if (l_clipped){
            // Prepend our remapped alignment portion
            new_path->CopyFrom( append_path(remap_path, anchor));
        }
        else if (!l_clipped){
            // Append our remapped alignment portion
            new_path->CopyFrom( append_path(anchor, remap_path));
        }

        // Determine if the new alignment's path represents an inversion, deletion, or insertion
        // We do this by checking the distance and relative orientation between Anchor
        // and remap_path.
        // 0: Unset, 1: INS, 2: DEL, 3: INV, 4: DUP

        int sv_indicator = 0;
        Alignment bogus;
        Path* b_path = bogus.mutable_path();
        b_path->CopyFrom(anchor);
        Alignment bogus_right;
        Path* br_path = bogus_right.mutable_path();
        br_path->CopyFrom(remap_path);
        if ( pair_orientation_filter(bogus, bogus_right)){
            sv_indicator = 3;
        }
        if (!node_to_position.empty()){
            Position other = l_clipped ? last_path_position(remap_path) : first_path_position(remap_path);
            int dist = distance_between_positions(get_clipped_position(a), other);
            int internal_homology_limit = 25;
            if (dist > internal_homology_limit){
                sv_indicator = 2;
            }
            else{
                sv_indicator = 1;
            }
        }

        return make_pair(ret, sv_indicator);

    }

    bool Filter::is_left_clipped(Alignment& a){
        Path path = a.path();
        Edit left_edit = path.mapping(0).edit(0);
        Edit right_edit = path.mapping(path.mapping_size() - 1).edit(path.mapping(path.mapping_size() - 1).edit_size() - 1);
        int left_overhang = left_edit.to_length() - left_edit.from_length();
        int right_overhang = right_edit.to_length() - right_edit.from_length();
        if (left_overhang > soft_clip_limit){
            return true;
        }
        else if (right_overhang > soft_clip_limit){
            return false;
        }
        else{
            cerr << "ALignment is not softclipped" << endl;
            exit(999);
        }
    }

    string Filter::get_clipped_seq(Alignment& a){
        if (a.path().mapping_size() > 0){
            Path path = a.path();
            Edit left_edit = path.mapping(0).edit(0);
            Edit right_edit = path.mapping(path.mapping_size() - 1).edit(path.mapping(path.mapping_size() - 1).edit_size() - 1);
            int left_overhang = left_edit.to_length() - left_edit.from_length();
            int right_overhang = right_edit.to_length() - right_edit.from_length();
            if (left_overhang > soft_clip_limit){
                return left_edit.sequence();
            }
            else if (right_overhang > soft_clip_limit){
                return right_edit.sequence();
            }
            else{
                cerr << "WARNING: BOTH ENDS CLIPPED" << endl
                << "IGNORING READ";
                return "";
            }

        }
        else{
            return "";
        }
    }

    Alignment Filter::remove_clipped_portion(Alignment& a){
        Alignment ret = a;
        ret.clear_path();
        Path clipped = trim_hanging_ends(a.path());
        ret.mutable_path()->CopyFrom(clipped);

        return ret;
    }


    Position Filter::get_clipped_position(Alignment& a){
        Position rPos;
        if (a.path().mapping_size() > 0){
            Path path = a.path();
            Edit left_edit = path.mapping(0).edit(0);
            Edit right_edit = path.mapping(path.mapping_size() - 1).edit(path.mapping(path.mapping_size() - 1).edit_size() - 1);
            int left_overhang = left_edit.to_length() - left_edit.from_length();
            int right_overhang = right_edit.to_length() - right_edit.from_length();
            
            if (left_overhang > soft_clip_limit){
                /// Get the position of the first match
                bool tripped = false;
                int i = 0; 
                for (i = 0; i < path.mapping_size(); i++){
                    Mapping m = path.mapping(i);
                    Position pp = m.position();
                    for (int j = 0; j < m.edit_size(); j++){
                        Edit e = m.edit(j);
                        if (e.to_length() == e.from_length() && e.sequence().empty()){
                            tripped = true;
                            rPos = pp;
                            break;
                        }
                    }
                    if (tripped){
                        break;
                    }
                }
                if (i >= path.mapping_size()){
                    cerr << "ERROR: ALIGNMENT MATCHES ALONG WHOLE LENGTH" << a.name() << endl;
                    exit(394);
                }

            }
            else if (right_overhang > soft_clip_limit){
                bool tripped = false;
                /// Get the position of the first match
                int i; 
                for (i = path.mapping_size() - 1; i >= 0; i--){
                    Mapping m = path.mapping(i);
                    Position pp = m.position();
                    for (int j = m.edit_size() - 1; j >= 0; j--){
                        Edit e = m.edit(j);
                        if (e.to_length() == e.from_length() && e.sequence().empty()){
                            tripped = true;
                            rPos = pp;
                            break;
                        }
                    }
                    if (tripped){
                        break;
                    }
                }

            }
            else{
                cerr << "WARNING: BOTH ENDS CLIPPED" << endl
                << "IGNORING READ";
            }

        }
        return rPos;
    }

    int64_t Filter::get_clipped_ref_position(Alignment& a){
        if (a.path().mapping_size() > 0){
            Path path = a.path();
            Edit left_edit = path.mapping(0).edit(0);
            Edit right_edit = path.mapping(path.mapping_size() - 1).edit(path.mapping(path.mapping_size() - 1).edit_size() - 1);
            int left_overhang = left_edit.to_length() - left_edit.from_length();
            int right_overhang = right_edit.to_length() - right_edit.from_length();
            
            if (left_overhang > soft_clip_limit){
                /// Get the position of the first match
                bool tripped = false;
                int i = 0; 
                for (i = 0; i < path.mapping_size(); i++){
                    Mapping m = path.mapping(i);
                    for (int j = 0; j < m.edit_size(); j++){
                        Edit e = m.edit(j);
                        if (e.to_length() == e.from_length() && e.sequence().empty()){
                            tripped = true;
                        }
                    }
                    if (tripped){
                        break;
                    }
                }
                if (i >= path.mapping_size()){
                    cerr << "ERROR: ALIGNMENT MATCHES ALONG WHOLE LENGTH" << a.name() << endl;
                    exit(394);
                }
                return node_to_position[path.mapping(i).position().node_id()] + 
                (path.mapping(i).position().is_reverse() ? (-1 * path.mapping(i).position().offset()) : path.mapping(i).position().offset());
            }
            else if (right_overhang > soft_clip_limit){
                bool tripped = false;
                /// Get the position of the first match
                int i; 
                for (i = path.mapping_size() - 1; i >= 0; i--){
                    Mapping m = path.mapping(i);
                    for (int j = m.edit_size() - 1; j >= 0; j--){
                        Edit e = m.edit(j);
                        if (e.to_length() == e.from_length() && e.sequence().empty()){
                            tripped = true;
                        }
                    }
                    if (tripped){
                        break;
                    }
                }
                return node_to_position[path.mapping(i).position().node_id()] + 
                (path.mapping(i).position().is_reverse() ? (-1 * path.mapping(i).position().offset()) : path.mapping(i).position().offset());

            }
            else{
                cerr << "WARNING: BOTH ENDS CLIPPED" << endl
                << "IGNORING READ";
                return 0;
            }

        }
        return 0;
    }

    void Filter::fill_node_to_position(string pathname){
        if (my_vg == NULL){
            cerr << "VG must be provided to use node_to_position" << endl;
            exit(1);
        }
        if (my_vg->paths._paths.count(pathname)){
            auto& maps = my_vg->paths._paths[pathname];
            int64_t dist_from_start = 0;
            for (auto m : maps){
                node_to_position[m.node_id()] = dist_from_start;
                dist_from_start += my_vg->get_node(m.node_id())->sequence().size();
            }
        }
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

    void Filter::set_my_path_position_graph(PathPositionHandleGraph* graph){
        my_path_position_graph = graph;
    }

    void Filter::set_inverse(bool do_inv){
        inverse = do_inv;
    }

    void Filter::init_mapper(){
        if (my_path_position_graph == NULL || gcsa_ind == NULL || lcp_ind == NULL){
            cerr << "Must provide an xg and gcsa to inititiate mapper for split read mapping."
            << endl;
            exit(1);
        }
        my_mapper = new Mapper(my_path_position_graph, gcsa_ind, lcp_ind);
    }

    bool Filter::perfect_filter(Alignment& aln){
        for (int i = 0; i < aln.path().mapping_size(); i++){
            Mapping m = aln.path().mapping(i);
            for (int j = 0; j < m.edit_size(); j++){
                Edit e = m.edit(j);
                if (e.to_length() != e.from_length() || !e.sequence().empty()){
                    return false;
                }
            }
        }
        return true;
    }

    bool Filter::mark_smallVariant_alignments(Alignment& a, Alignment& b){
        return (perfect_filter(a) || perfect_filter(b) || anchored_filter(a) || anchored_filter(b) );
    }

    bool Filter::mark_sv_alignments(Alignment& a, Alignment& b){
        bool ret = false;
        ret = (ret | soft_clip_filter(a));
        ret = (ret | soft_clip_filter(b));
        ret = (ret | interchromosomal_filter(a, b));
        ret = (ret | one_end_anchored_filter(a, b));
        ret = (ret | (unmapped_filter(a) & unmapped_filter(b)));
        ret = (ret | insert_size_filter(a, b));
        ret = (ret | pair_orientation_filter(a, b));
        return ret;
    }
    bool Filter::anchored_filter(Alignment& aln){
        int min_match_len = 10;
        if (aln.path().mapping_size() > 2){
            return false;
        }
        Mapping map_end = aln.path().mapping( aln.path().mapping_size() - 1);
        Edit front = aln.path().mapping(0).edit(0);
        Edit back = map_end.edit( map_end.edit_size() - 1);

        bool valid = true;
        if (front.to_length() != front.from_length() ||
             front.to_length() < min_match_len ||
             !front.sequence().empty()){
                valid = false;
        }
        if (back.to_length() != back.from_length() ||
             back.to_length() < min_match_len ||
             !back.sequence().empty()){
                valid = false;
        }
        return valid;
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

        return inverse ? Alignment() : aln;

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

    bool Filter::unmapped_filter(Alignment& aln){
        if (aln.score() == 0 || aln.path().mapping_size() == 0){
            aln.set_read_mapped(false);
            return true;
        }
        aln.set_read_mapped(true);
        return false;
    
    }

    /*PE Functions*/
    bool Filter::one_end_anchored_filter(Alignment& aln_first, Alignment& aln_second){
        bool f = unmapped_filter(aln_first);
        bool s = unmapped_filter(aln_second);
        if ( (f && !s)){ 
            aln_first.set_read_mapped(false);
            aln_first.set_mate_unmapped(false);
            aln_second.set_read_mapped(true);
            aln_second.set_mate_unmapped(true);
            return true;
        }
        else if((s && !f)){
            aln_first.set_read_mapped(true);
            aln_first.set_mate_unmapped(true);
            aln_second.set_read_mapped(false);
            aln_second.set_mate_unmapped(false);
            return true;
        }
        else{
            return false;
        }
    }

    bool Filter::interchromosomal_filter(Alignment& aln_first, Alignment& aln_second){
        if (aln_first.path().name() != aln_second.path().name() && !aln_first.path().name().empty() && !aln_second.path().name().empty()){
            return true;
        }
        else{
            return false;
        }
    }

    bool Filter::insert_size_filter(Alignment& aln_first, Alignment& aln_second){

        double zed;
        bool check_first = true;
        bool check_second = true;
        if (aln_first.fragment_size() < 1){
            check_first = false;
        }
        if (aln_second.fragment_size() < 1){
            check_second = false;
        }

        if (check_first && !check_second){
            zed = ((double) aln_first.fragment(0).length() - (double) insert_mean) / (double) insert_sd;
        }
        else if (check_second){
            zed = ((double) aln_second.fragment(0).length() - (double) insert_mean) / (double) insert_sd;
        }else{
            return false;
        }

        if (zed >= 1.95 || zed <= -1.95){
            return true;
        }
        else{
            return false;
        }
    }

    bool Filter::pair_orientation_filter(Alignment& aln_first, Alignment& aln_second){

        // TODO need to check the innie/outie case
        // --->    <--- normal
        // and
        // <---    ---> not so normal
        // plus the reversing edge
        // -->  -->
        // <--  <--
        bool f_rev = false;
        bool s_rev = false;
        
        if (! (aln_first.mapping_quality() > 0 && aln_second.mapping_quality() > 0)){
            return false;
        }

        Path f_path = aln_first.path();
        Path s_path = aln_second.path();
        // for (int i = 0; i < aln_first.fragment_size(); i++){
        //     cerr << aln_first.name() << " " << aln_first.fragment(i).length() << endl;
        // }
        // for (int i = 0; i < aln_second.fragment_size(); i++){
        //     cerr << aln_second.name() << " " << aln_second.fragment(i).length() << endl;
        // }
        bool flipped = false;
        if (aln_first.fragment_size() > 0){
            flipped = aln_first.fragment(0).length() < 0;
        }

        for (int i = 0; i < f_path.mapping_size(); i++){
            if (f_path.mapping(i).position().is_reverse()){
                f_rev = true;
                aln_first.set_read_on_reverse_strand(true);
                aln_second.set_mate_on_reverse_strand(true);
                aln_first.set_read_mapped(true);
                aln_second.set_mate_unmapped(false);
            }
        }

        for (int j = 0; j < s_path.mapping_size(); j++){
            if (s_path.mapping(j).position().is_reverse()){
                s_rev = true;
                aln_second.set_read_mapped(true);
                aln_first.set_mate_unmapped(false);
                aln_second.set_read_on_reverse_strand(true);
                aln_first.set_mate_on_reverse_strand(true);
            }
        }

        if (!f_rev){
            aln_second.set_mate_on_reverse_strand(false);
            aln_first.set_read_on_reverse_strand(false);
        }
        if (!s_rev){
            aln_first.set_mate_on_reverse_strand(false);
            aln_second.set_read_on_reverse_strand(false);
        }
        if (f_rev == s_rev){
            return true;
        }
        else if ( ((f_rev != s_rev) && flipped)){
            return true;
        }
        else{
            return false;
        }

    }

    pair<Alignment, Alignment> Filter::deletion_filter(Alignment& aln_first, Alignment& aln_second){
        // path_length, split read
        // REFINE USING SOFT CLIPS
        //pair<Alignment, Alignment> ret_alns = path_length_filter(aln_first, aln_second);
        //pair<Alignment, Alignment> x_alns = split_read_filter(ret_alns.first, ret_alns.second);
        
        bool found = false;
        //if (inverse | (ret_alns.first.name() != "" || ret_alns.second.name() != "")){
        //    found = true;
        //}

        
        if (found){
            //return ret_alns;
            throw runtime_error("Unimplemented!");
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
    pair<Locus, Locus> Filter::insertion_filter(Alignment& aln_first, Alignment& aln_second) {
        throw runtime_error("Unimplemented!");
    }

    /**
    * Find reads that support duplications
    *
    */
    pair<Locus, Locus> Filter::duplication_filter(Alignment& aln_first, Alignment& aln_second) {
        throw runtime_error("Unimplemented!");
    }

    /**
    * Find reads that support inversions
    * split reads
    * discordant insert size
    * bad orientation
    * instead of ---->    <-----
    * we'll see  <----    <----- or ---->    ----->
    */
    bool Filter::inversion_filter(Alignment& aln_first, Alignment& aln_second){
        if (pair_orientation_filter(aln_first, aln_second)){
            return true;
        }
        // TODO: Full implementation
        return false;
    }

    /**
    * split reads or discordant reads/insert size may indicate a breakend but not a clean SV type
    * we'd like to report all possible breakends, even if that don't match an SV type very well.
    */
    pair<Locus, Locus> Filter::breakend_filter(Alignment& aln_first, Alignment& aln_second){
        throw runtime_error("Unimplemented!");
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
        const Path& path = aln.path();
        if (path.mapping_size() > 0) {
            
            handle_t first_handle = my_path_position_graph->get_handle(path.mapping(0).position().node_id());
            unordered_set<path_handle_t> prev_node_paths;
            my_path_position_graph->for_each_step_on_handle(first_handle, [&](const step_handle_t& step) {
                prev_node_paths.insert(my_path_position_graph->get_path_handle_of_step(step));
            });
            
            for (int i = 1; i < path.mapping_size(); i++){
                
                handle_t handle = my_path_position_graph->get_handle(path.mapping(i).position().node_id());
                unordered_set<path_handle_t> curr_node_paths;
                my_path_position_graph->for_each_step_on_handle(handle, [&](const step_handle_t& step) {
                    curr_node_paths.insert(my_path_position_graph->get_path_handle_of_step(step));
                });
                
                bool paths_match = false;
                
                for (const path_handle_t& path : curr_node_paths) {
                    if (prev_node_paths.count(path)) {
                        paths_match = true;
                        break;
                    }
                }
                
                if (!paths_match){
                    return inverse ? aln : Alignment();
                }
                
                prev_node_paths = std::move(curr_node_paths);
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
                return aln;
            }

        }
        return Alignment();

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
        throw runtime_error("Unimplemented!");
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


    bool Filter::soft_clip_filter(Alignment& aln){
        //Find overhangs - portions of the read that
        // are inserted at the ends.
        if (aln.path().mapping_size() > 0){
            Path path = aln.path();
            Edit left_edit = path.mapping(0).edit(0);
            Edit right_edit = path.mapping(path.mapping_size() - 1).edit(path.mapping(path.mapping_size() - 1).edit_size() - 1);
            int left_overhang = left_edit.to_length() - left_edit.from_length();
            int right_overhang = right_edit.to_length() - right_edit.from_length();
            if (left_overhang > soft_clip_limit || right_overhang > soft_clip_limit){
                aln.set_soft_clipped(true);
                return true;
            }
            else{
                aln.set_soft_clipped(false);
                return false;
            }
        }
        else{
            if (aln.sequence().length() > soft_clip_limit){
                aln.set_soft_clipped(false);
                return false;
            }
            cerr << "WARNING: SHORT ALIGNMENT: " << aln.sequence().size() << "bp" << endl
                << "WITH NO MAPPINGS TO REFERENCE" << endl
                << "CONSIDER REMOVING IT FROM ANALYSIS" << endl;
            return false;
        }

    }

    

    vector<Alignment> Filter::remap(Alignment& aln){
        if (this->my_path_position_graph == NULL || this->gcsa_ind == NULL || this->my_mapper == NULL){
            cerr << "An XG and GCSA are required for remapping." << endl;
            exit(1337);
        }

        vector<Alignment> match = this->my_mapper->align_multi(aln);
        return match;
    }

    vector<Alignment> Filter::remap(string seq){
        if (this->my_path_position_graph == NULL || this->gcsa_ind == NULL){
            cerr << "An XG and GCSA are required for remapping." << endl;
            exit(1337);
        }

        Alignment x;
        x.set_sequence(seq);

        vector<Alignment> ret;
        ret = this->my_mapper->align_multi(x);

        

        return ret;
    }
    /**
     * Split reads map to two separate paths in the graph OR vastly separated non-consecutive
     * nodes in a single path.
     *
     * They're super important for detecting structural variants, so we may want to
     * filter them out or collect only split reads.
     */
    bool Filter::split_read_filter(Alignment& aln){

        if (this->my_path_position_graph == NULL || this->gcsa_ind == NULL){
            cerr << "An XG and GCSA are required for split read processing." << endl;
            exit(1337);
        }
        bool flagged = false;
        // Check softclips

        if (soft_clip_filter(aln)){
            flagged = true;
            string clipseq = get_clipped_seq(aln);
            Alignment clipmatch = my_mapper->align(clipseq);
            cerr << clipmatch.path().mapping_size() << endl;
            if (clipmatch.path().mapping_size() < 1 || clipmatch.mapping_quality() < 5){
                flagged = false;
            }
            else{
                clipmatch.set_name(aln.name() + "_split" );
            }
        }

        if (flagged){
            return true;
        }

        return false;


    
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

#include "filter.hpp"

using namespace std;

namespace vg{

    Filter::Filter(){

    }

    Filter::~Filter(){

    }
    
    bool Filter::mark_smallVariant_alignments(Alignment& a, Alignment& b){
        return (is_perfect(a) || is_perfect(b) || anchored_filter(a) || anchored_filter(b) );
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

    bool Filter::unmapped_filter(Alignment& aln){
        if (!is_mapped(aln)){
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
}

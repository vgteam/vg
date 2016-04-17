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
        if (prev | pos.is_reverse()){
                return (inverse ? aln : Alignment());
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
  void Filter::set_softclip_filter(int max_clip){
    max_softclip = max_clip;
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
        #pragma omp critical(write)
        pos_to_edit_to_depth[p_hash][e_hash] += 1;
        /**
        * If an edit fails the filter, either return a new empty alignment
        * OR
        * return a new alignment identical to the old one EXCEPT where
        * the offending edit has been replaced by a match to the reference.
        */
        if (pos_to_edit_to_depth[p_hash][e_hash] < min_depth){
          if (remove_failing_alignments){
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


  Alignment Filter::qual_filter(Alignment& aln){
    string quals = aln.quality();
    int offset = 0;
    for (int i = 0; i < quals.size(); i++){
      if (((int) quals[i] - offset) < min_qual){
        if (remove_failing_alignments){
          return Alignment();
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
    int win_len = window_len; // >= 1 ? window_len : aln.quality().size();

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
        Edit right_edit = path.mapping(path.mapping_size() - 1).edit(path.mapping(0).edit_size() - 1);
        int left_overhang = left_edit.to_length() - left_edit.from_length();
        int right_overhang = right_edit.to_length() - right_edit.from_length();
        if (left_overhang > max_softclip || right_overhang > max_softclip){
            return Alignment();
        }
        return inverse ? Alignment() : aln;
    }
    else{
        if (aln.sequence().length() > max_softclip){
            return Alignment();
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
    
    int top_side = path.mapping_size();
    int bottom_side = 0;
    bool diverges = false;

    string main_path = "";
    while (top_side > bottom_side){
       //main_path = path_of_node(path.mapping(bottom_side);
       //
        //Check if paths are different
        //if (divergent(node1, node2){
        //    diverges = true;
        //}
        top_side--;
        bottom_side++;
    }

    //TODO check this; I'm not sure this is a fully safe conditional:
    // reg: diverges: return ALN()
    // reg: !diverges: return aln
    // inv: diverges: return return aln;
    // inv: !diverges: return ALN()
    if ((diverges)){
        return inverse ? aln : Alignment();
    }
    else{
        return inverse ? Alignment() : aln;
    }

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

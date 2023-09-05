#ifndef VG_FILTER
#define VG_FILTER

#include <vector>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <string>
#include <math.h>
#include "vg.hpp"
#include "mapper.hpp"
#include <vg/vg.pb.h>

/** \file
 * Provides a way to filter Edits contained
 * within Alignments. This can be used to clean out
 * sequencing errors and to find high-quality candidates
 * for variant calling.
 *
 */
namespace vg{

struct BREAKPOINT{
    string name;
    Position position;
    vector<BREAKPOINT> mates;
    
    string contig;
    int64_t start = -1;
    int64_t upper_bound = 100;
    int64_t lower_bound = 100;
    
    // Does the breakpoint point this way --->>
    // or this way <<---
    bool isForward;
    // 0: Unset, 1: INS, 2: DEL, 3: INV, 4: DUP
    int SV_TYPE = 0;
    //
    
    int normal_supports = 0;
    int tumor_supports = 0;
    
    int fragl_supports = 0;
    int split_supports = 0;
    int other_supports = 0;
    
    inline int total_supports(){
        return fragl_supports + split_supports + other_supports;
    }
    inline bool overlap(BREAKPOINT p, int dist){
        
        if (start > -1 ){
            if ( abs(start - p.start) < dist){
                return true;
            }
        }
        else{
            if (position.node_id() == p.position.node_id() && abs(position.offset() - p.position.offset()) < dist){
                return true;
            }
        }
        
        return false;
    }
    inline string to_string(){
        stringstream x;
        x << "Pos: " << start << " u: " << upper_bound << " l: " << lower_bound << " s: " << total_supports();
        return x.str();
    }
    
};

class Filter{
    public:
        Filter();
        ~Filter();
        // map<string, interval> sv-to-interval
         /* Filter functions.
         * Take an Alignment and walk it.
         * Perform the desired summary statistics.
         * Output the Alignment if it passes OR
         * A new, modified Alignment if we allow it OR
         * an empty Alignment if the alignment fails and we don't allow
         * modified alignments.
         */
        bool perfect_filter(Alignment& aln);
        bool anchored_filter(Alignment& aln);
        bool mark_sv_alignments(Alignment& a, Alignment& b);
        bool mark_smallVariant_alignments(Alignment& a, Alignment& b);
        Alignment depth_filter(Alignment& aln);
        Alignment qual_filter(Alignment& aln);
        Alignment coverage_filter(Alignment& aln);
        Alignment avg_qual_filter(Alignment& aln);
        Alignment percent_identity_filter(Alignment& aln);
       
        bool soft_clip_filter(Alignment& aln);
        bool unmapped_filter(Alignment& aln);
        bool split_read_filter(Alignment& aln);
        Alignment path_divergence_filter(Alignment& aln);
        Alignment reversing_filter(Alignment& aln);

        vector<Alignment> remap(Alignment& aln);
        vector<Alignment> remap(string seq);

        bool is_left_clipped(Alignment& a);

        // Returns a new alignment which doesn't have a split path
        // and an integer indicating the SV type it might support
        // 0: Unset, 1: INS, 2: DEL, 3: INV, 4: DUP
        pair<Alignment, int> refactor_split_alignment(Alignment& a);

        // Returns whether a split read supports a DEL, INS, or INV
        // int split_read_type(Alignment& a);

        Alignment path_length_filter(Alignment& aln);

        /*PE Functions*/
        bool one_end_anchored_filter(Alignment& aln_first, Alignment& aln_second);
       bool interchromosomal_filter(Alignment& aln_first, Alignment& aln_second);
        
        // // TODO should give this one an insert size arg
        bool insert_size_filter(Alignment& aln_first, Alignment& aln_second);
        bool pair_orientation_filter(Alignment& aln_first, Alignment& aln_second);

        // pair<Alignment, Alignment> path_length_filter(Alignment& aln_first, Alignment& aln_second);

        // SV filters
        // Take in paired GAM and return Locus records
        pair<Alignment, Alignment> deletion_filter(Alignment& aln_first, Alignment& aln_second);
        pair<Locus, Locus> insertion_filter(Alignment& aln_first, Alignment& aln_second);
        pair<Locus, Locus> duplication_filter(Alignment& aln_first, Alignment& aln_second);
        bool inversion_filter(Alignment& aln_first, Alignment& aln_second);
        pair<Locus, Locus> breakend_filter(Alignment& aln_first, Alignment& aln_second);

        void set_min_depth(int depth);
        //void set_min_kmer_depth(int d);

        void set_min_qual(int qual);
        void set_min_percent_identity(double pct_id);
        void set_avg(bool do_avg);
        void set_filter_matches(bool fm);
        void set_remove_failing_edits(bool fm);
        void set_soft_clip_limit(int max_clip);
        void set_split_read_limit(int split_limit);
        void set_window_length(int window_length);
        void set_my_vg(vg::VG* vg);
        void set_my_path_position_graph(PathPositionHandleGraph* graph);
        void set_inverse(bool do_inv);

        void init_mapper();

        vg::VG* my_vg = NULL;
        PathPositionHandleGraph* my_path_position_graph = NULL;
        gcsa::GCSA* gcsa_ind;
        gcsa::LCPArray * lcp_ind;
        Mapper* my_mapper;
        
        map<int64_t, int64_t> node_to_position;
        void fill_node_to_position(string pathname);
        int64_t distance_between_positions(Position first, Position second);
        string get_clipped_seq(Alignment& a);
        int64_t get_clipped_ref_position(Alignment& a);
        Position get_clipped_position(Alignment& a);
        Alignment remove_clipped_portion(Alignment& a);
        //Position: NodeID + offset
        // different edits may be present at each position.
        // is there some way to just hash the mappings?
        unordered_map<string, unordered_map<string, int> > pos_to_edit_to_depth;
        unordered_map<int, int> pos_to_qual;

        // Map position/interval to locus
        // map Locus to SV evidence, so that we can augment
        // it with each additional read
    public:
        // we really need a reservoir sampling method /
        // some way to effectively calculate less-biased moving averages.
        bool inverse = false;
        bool do_remap = false;
        bool remove_failing_edits = false;
        bool filter_matches = false;
        bool use_avg = false;;
        int min_depth = 0;
        int min_qual = 0;
        int min_cov = 0;
        int window_length = 0;
        int qual_offset = 0;
        int soft_clip_limit = -1;
        int split_read_limit = -1;
        double min_percent_identity = 0.0;
        double min_avg_qual = 0.0;

        int max_path_length = 0;

        int my_max_distance = 1000;

        float insert_mean = 1000;
        float insert_sd = 100;


        };
}

#endif

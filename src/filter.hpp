#ifndef VG_FILTER
#define VG_FILTER

#include <vector>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <string>
#include "vg.hpp"
#include "xg.hpp"
#include "vg.pb.h"

/** \file
 * Provides a way to filter Edits contained
 * within Alignments. This can be used to clean out
 * sequencing errors and to find high-quality candidates
 * for variant calling.
 *
 */
namespace vg{

struct SV_EVIDENCE{
    int SR = 0;
    int PE = 0;
    bool PRECISE = false;

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
        bool simple_filter(Alignment& aln);
        Alignment depth_filter(Alignment& aln);
        Alignment qual_filter(Alignment& aln);
        Alignment coverage_filter(Alignment& aln);
        Alignment avg_qual_filter(Alignment& aln);
        Alignment percent_identity_filter(Alignment& aln);
        Alignment soft_clip_filter(Alignment& aln);
        Alignment split_read_filter(Alignment& aln);
        Alignment path_divergence_filter(Alignment& aln);
        Alignment reversing_filter(Alignment& aln);

        Alignment path_length_filter(Alignment& aln);

        /*PE Functions*/
        pair<Alignment, Alignment> one_end_anchored_filter(Alignment& aln_first, Alignment& aln_second);
        pair<Alignment, Alignment> interchromosomal_filter(Alignment& aln_first, Alignment& aln_second);
        
        // // TODO should give this one an insert size arg
        pair<Alignment, Alignment> insert_size_filter(Alignment& aln_first, Alignment& aln_second);
        pair<Alignment, Alignment> pair_orientation_filter(Alignment& aln_first, Alignment& aln_second);

        // pair<Alignment, Alignment> path_length_filter(Alignment& aln_first, Alignment& aln_second);

        // SV filters
        // Take in paired GAM and return Locus records
        pair<Alignment, Alignment> deletion_filter(Alignment& aln_first, Alignment& aln_second);
        pair<Locus, Locus> insertion_filter(Alignment& aln_first, Alignment& aln_second);
        pair<Locus, Locus> duplication_filter(Alignment& aln_first, Alignment& aln_second);
        pair<Locus, Locus> inversion_filter(Alignment& aln_first, Alignment& aln_second);
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
        void set_my_xg_idx(xg::XG* xg_idx);
        void set_inverse(bool do_inv);

        vg::VG* my_vg;
        xg::XG* my_xg_index = NULL;
        gcsa::GCSA* gcsa_ind;
        gcsa::LCPArray * lcp_ind;
 
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

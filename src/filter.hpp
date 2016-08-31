#ifndef VG_FILTER_HPP
#define VG_FILTER_HPP

#include <vector>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <string>
#include "vg.hpp"
#include "xg.hpp"
#include "vg.pb.h"

/**
 * Provides a way to filter Edits contained
 * within Alignments. This can be used to clean out
 * sequencing errors and to find high-quality candidates
 * for variant calling.
 *
 */
namespace vg{
class Filter{
    public:
        Filter();
        ~Filter();
         /* Filter functions.
         * Take an Alignment and walk it.
         * Perform the desired summary statistics.
         * Output the Alignment if it passes OR
         * A new, modified Alignment if we allow it OR
         * an empty Alignment if the alignment fails and we don't allow
         * modified alignments.
         */
        Alignment depth_filter(Alignment& aln);
        Alignment qual_filter(Alignment& aln);
        Alignment coverage_filter(Alignment& aln);
        Alignment avg_qual_filter(Alignment& aln);
        Alignment percent_identity_filter(Alignment& aln);
        Alignment soft_clip_filter(Alignment& aln);
        Alignment split_read_filter(Alignment& aln);
        Alignment path_divergence_filter(Alignment& aln);
        Alignment reversing_filter(Alignment& aln);
        Alignment kmer_filter(Alignment& aln);
        void set_min_depth(int depth);
        //void set_min_kmer_depth(int d);
        void set_min_qual(int qual);
        void set_min_percent_identity(double pct_id);
        void set_avg_qual(double avg_qual);
        void set_filter_matches(bool fm);
        void set_remove_failing_edits(bool fm);
        void set_soft_clip_limit(int max_clip);
        void set_split_read_limit(int split_limit);
        void set_reversing(bool do_reversing_filter);
        void set_path_divergence(bool do_path_divergence);
        void set_window_length(int window_length);
        void set_my_vg(vg::VG* vg);
        void set_my_xg_idx(xg::XG* xg_idx);
        void set_inverse(bool do_inv);

        int get_min_depth();
        int get_min_qual();
        int get_window_length();
        int get_soft_clip_limit();
        int get_split_read_limit();
        double get_min_percent_identity();
        double get_min_avg_qual();
        bool get_inverse();
        bool get_filter_matches();
        bool get_remove_failing_edits();
        bool get_do_path_divergence();
        bool get_do_reversing();


    private:
        vg::VG* my_vg;
        xg::XG* my_xg_idx;
        //Position: NodeID + offset
        // different edits may be present at each position.
        // is there some way to just hash the mappings?
        unordered_map<string, unordered_map<string, int> > pos_to_edit_to_depth;
        unordered_map<int, int> pos_to_qual;
        bool inverse = false;
        bool remove_failing_edits = false;
        bool filter_matches = false;
        bool do_path_divergence;
        bool do_reversing;
        int min_depth = 0;
        int min_qual = 0;
        int min_cov = 0;
        int window_length = 0;
        int qual_offset = 0;
        int soft_clip_limit = -1;
        int split_read_limit = -1;
        double min_percent_identity = 0.0;
        double min_avg_qual = 0.0;
        };
}
#endif

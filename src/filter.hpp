#include <vector>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <string>
#include "vg.hpp"
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
         * a NULLPTR if the alignment fails and we don't allow
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
        void set_min_depth(int depth);
        void set_min_qual(int qual);
        void set_min_percent_identity(double pct_id);
        void set_avg_qual(double avg_qual);
        void set_filter_matches(bool fm);
        void set_remove_failing_alignments(bool fm);
        void set_softclip_filter(int max_clip);
        void set_splitread_filter(bool do_splitread);
        void set_reversing_filter(bool do_reversing_filter);

    private:
        vg::VG* my_vg;
        //Position: NodeID + offset
        // different edits may be present at each position.
        // is there some way to just hash the mappings?
        unordered_map<string, unordered_map<string, int> > pos_to_edit_to_depth;
        unordered_map<int, int> pos_to_qual;
        bool inverse = false;
        bool do_splitread = false;
        bool remove_failing_alignments = true;
        bool filter_matches = true;
        int min_depth = 0;
        int min_qual = 0;
        int min_cov = 0;
        int window_len = 0;
        int max_softclip = -1;
        double min_percent_identity = 0.0;
        double min_avg_qual = 0.0;
        };
}

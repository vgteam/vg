#include <vector>
#include <iostream>
#include <unordered_map>
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
        /**
         * Filter functions.
         * Take an Alignment and walk it.
         * Perform the desired summary statistics.
         * Output the Alignment if it passes OR
         * A new, modified Alignment if we allow it OR
         * a NULLPTR if the alignment fails and we don't allow
         * modified alignments.
         */
        Alignment depth_filter(Alignment& aln);
        Alignment& qual_filter(Alignment& aln);
        Alignment& coverage_filter(Alignment& aln);
        Alignment& avg_qual_filter(Alignment& aln);
        Alignment& percent_identity_filter(Alignment& aln);
        Alignment& soft_clip_filter(Alignment& aln);
        void set_min_depth(int depth);
        void set_min_qual(int qual);
        void set_min_pct_identity(double pct_id);
        void set_avg_qual(double avg_qual);

    private:
        unordered_map<int, int> pos_to_depth;
        unordered_map<int, int> pos_to_qual;
        int min_depth;
        int min_qual;
        int min_cov;
        double avg_qual;
        double min_pct_identity;
        bool remove_failing_edits = false;


};
}

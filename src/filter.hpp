#include <vector>
#include <iostream>
#include <unordered_map>

/**
 * Provides a way to filter Edits contained
 * within Alignments. This can be used to clean out
 * sequencing errors and to find high-quality candidates
 * for variant calling.
 *
 */
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
        Alignment& depth_filter(Alignment& aln, int depth);
        Alignment& qual_filter(Alignment& aln, int qual);
        Alignment& coverage_filter(Alignment& aln, int depth);
        Alignment& cov_qual_filter(Alignment& aln, int avg_qual);
        Alignment& soft_clip_filter(Alignment& aln);
    private:
        unordered_map<int, int> pos_to_depth;
        unordered_map<int, int> pos_to_qual;

}

#include "filter.hpp"
#include "index.hpp"
#include "vg.pb.h"

/**
 * Overview:
 *      Use the GAM/GAM index and a filter to locate Alignments
 *      which may indicate the presence of
 *      structural variants at a given site.
 *
 *      Signatures include:
 *          Deletions/Insertions: Stacked soft clips (tips)
 *          Inversions: mismatched P/E reads(  <-- && -->  rather than the expected ( -->   <-- )
 *          Duplications: Read depth signals
 *          Translocations: Distant read pairs
 */

using namespace std;
namespace vg{
class SRPE{
    public:
        /* Filter an entire GAM, using the index to look for supporting evidence when a signature is seen
         * in a single Alignment.*/
        void apply(Alignment& aln, vector<Alignment>& suspect);
        void apply(Alignment& aln, vector<Position>& suspect);
        bool in_suspects(Position p);
        string bucket(Position p);
        void pre_depth(Alignment& aln);
        string pos_to_string(Position& p);
        void call_depth_variants(void);
        void call(void);
        
    private:
        Filter my_filter;
        int avg_depth;
        // int depth_window;
        // Consider variants within my_pos_window_size of each other on either
        // side (WITHIN THE SAME NODE) the same variant.

        int my_bucket_size = 50;
        map<string, int> pos_to_depth;
        // Use a map and some modulo division to bucket variants.

};
}

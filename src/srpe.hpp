#include "filter.hpp"
#include "index.hpp"
#include "vg.pb.h"
#include "vg.cpp"


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
//    vcfparse::variant locus_to_sv_vcf(Locus ll);
    string locus_to_sv_vcf(Locus ll);
    pair<Locus, Locus> pe_aln_to_locus(Alignment& aln_one, Alignment& aln_two);
    
    vector<Locus> refine(vector<Locus> loci, vg::VG* graph, index gam_ind);

    private:
        Filter ff;

};
}

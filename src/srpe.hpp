#ifndef VG_SRPE
#define VG_SRPE

#include "filter.hpp"
#include "index.hpp"
#include "vg.pb.h"
#include "vg.hpp"
#include "gcsa.h"
#include "alignment.hpp"
#include "mapper.hpp"
using namespace std;
namespace vg{

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


    class SRPE{

    string locus_to_sv_vcf(Locus ll);

   void remap(vg::VG* graph, Index gam_index, vector<pair<Alignment, Alignment> >& remapped);

    void filter(vector<Alignment>& in_alns, vector<Alignment>& out_alns);
    public:
    vg::Filter ff;
        vg::VG* vg;


};
}
#endif

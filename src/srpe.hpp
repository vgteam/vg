#include "filter.hpp"
#include "index.hpp"
#include "vg.pb.h"
#include "vg.hpp"
#include "gcsa.h"
#include "alignment.hpp"
#include "mapper.hpp"

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
        void remap(vg::VG* graph, xg::XG* xg_index, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_index,  string gam_file, Index gam_index, vector<Alignment>& remapped);
        void transform(vector<Alignment>& alignments, vector<Locus>& sigs, vector<string>& svtypes);
        void call(vector<Locus>& sigs, vector<string>& vars);
        void filter(vector<Alignment>& in_alns, vector<Alignment>& out_alns);

        Filter ff;
        vg::VG* vg;


};
}

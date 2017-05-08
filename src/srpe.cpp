#include "srpe.hpp"

using namespace std;
namespace vg{

    double SRPE::discordance_score(vector<Alignment> alns, VG* subgraph){
    // Sum up the mapping scores
    // subtract the soft clips
    // the discordant insert size reads
    // and a flat penalty for discordant orientation

    }

    void SRPE::aln_to_bseq(Alignment& a, bseq1_t* read){
        read->seq = (char*) a.sequence().c_str();
        read->qual = (char*) a.quality().c_str();
        read->l_seq = a.sequence().length();
    }

    /*
        typedef struct {
	    int32_t len;      // length of sequence
	    int32_t nsr;      // number of supporting reads
	    char *seq;        // unitig sequence
	    char *cov;        // cov[i]-33 gives per-base coverage at i
	    int n_ovlp[2];    // number of 5'-end [0] and 3'-end [1] overlaps
	    fml_ovlp_t *ovlp; // overlaps, of size n_ovlp[0]+n_ovlp[1]
        } fml_utg_t;
    */

    /**
    * function assemble
    * inputs: a vector of Alignments to be assembled (based on their sequences)
    * outputs: 
    */
    void SRPE::assemble(vector<Alignment> alns, vector<fml_utg_t>& unitigs){
        int n_seqs, n_utgs;
        n_seqs = alns.size();
        bseq1_t* mr_bseqs = new bseq1_t [alns.size()];
        for (int i = 0; i < n_seqs; ++i){
            aln_to_bseq( alns[i], mr_bseqs + i );
        }
        fml_utg_t *utgs;
        fml_opt_t opt;
        fml_opt_init(&opt);
        utgs = fml_assemble(&opt, n_seqs, mr_bseqs, &n_utgs);
        for (int i = 0; i < n_utgs; ++i){
            unitigs.push_back( *(utgs + i) );
        }

        fml_utg_destroy(n_utgs, utgs);
    }

    void SRPE::assemble(string refpath, int64_t start_pos, int64_t end_pos, vector<fml_utg_t>& unitigs){
        // Get all alignments on <refpath> from <startpos> to <endpos>
    }
    void SRPE::assemble(int64_t node_id, int64_t pos, int window_size){

    }

    void SRPE::breakpoints_to_intervals(vector<BREAKPOINT> bps, vector<INS_INTERVAL>& ret, vector<INS_INTERVAL> existing){

    }

    void SRPE::intervals_to_variants(vector<INS_INTERVAL> intervals, vector<vcflib::Variant>& vars){

    }

    // void SRPE::remap(vg::VG* graph, Index gam_index, vector<pair<Alignment, Alignment> >& remapped){

    // }
    // void SRPE::filter(vector<Alignment>& in_alns, vector<Alignment>& out_alns){
    
    // }

    // check if our insert size for a given set of alignments falls within our expected range.
    /** bool SRPE::reasonable_insert_size(pair<Alignment, Alignment> mates){

    }

    void SRPE::normalize_pairs(vector<pair<Alignment&, Alignment&> > alns, vector<Flags> offenses, vector<pair<Alignment, Alignment> > normals){
        
        for (auto a : alns){
            
        }

    }
    void SRPE::normalize_singles(vector<Alignment&> alns, vector<Alignment&> normals){

    }

    **/

}


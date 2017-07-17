#include "srpe.hpp"

using namespace std;
namespace vg{

    double SRPE::discordance_score(vector<Alignment> alns, VG* subgraph){
    // Sum up the mapping scores
    // subtract the soft clips
    // the discordant insert size reads
    // and a flat penalty for discordant orientation

    }

    void SRPE::call_svs_paired_end(vg::VG* graph, ifstream& gamstream, vector<BREAKPOINT>& bps, string refpath){

    }

    void SRPE::call_svs_split_read(vg::VG* graph, ifstream& gamstream, vector<BREAKPOINT>& bps, string refpath){
        // We're going to do a bunch of split-read mappings now,
        // then decide if our orientations support an inversion, an insertion,
        // or a deletion.
    }


    void SRPE::call_svs(string graphfile, string gamfile, string refpath){
        vg::VG* graph;
        if (!graphfile.empty()){
            ifstream in(graphfile);
            graph = new VG(in, false);
        }
        ifstream gamstream;
        gamstream.open(gamfile);
        // Set up path index
        ff.set_my_vg(graph);
        ff.soft_clip_limit = 20;
        ff.fill_node_to_position(refpath);
    
        std::function<vector<BREAKPOINT> (vector<BREAKPOINT>)> merge_breakpoints = [](vector<BREAKPOINT> bps){
        vector<BREAKPOINT> ret;
        BREAKPOINT sent;
        sent.start = -100;
        ret.push_back(sent);
        for (int i = 0; i < bps.size(); i++){
            BREAKPOINT a = bps[i];
            bool merged = false;
            for (int j = 0; j < ret.size(); j++){
                if (ret[j].overlap(a, 20)){
                    ret[j].other_supports += 1;
                    merged = true;
                }
            }
            if (!merged){
                ret.push_back(a);
            }
        }
        return ret;
    };

    vector<BREAKPOINT> pe_bps;
    vector<BREAKPOINT> sr_bps;

    call_svs_paired_end(graph, gamstream, pe_bps, refpath);
    call_svs_split_read(graph, gamstream, sr_bps, refpath);
    vector<BREAKPOINT> pe_merged = merge_breakpoints(pe_bps);
    vector<BREAKPOINT> sr_merged = merge_breakpoints(sr_bps);
    vector<BREAKPOINT> merged;
    merged.insert(merged.begin(), pe_merged.begin(), pe_merged.end());
    merged.insert(merged.begin(), sr_merged.begin(), sr_merged.end());
    merged = merge_breakpoints(merged);

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



}


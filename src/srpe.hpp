#ifndef VG_SRPE
#define VG_SRPE
#include <string>
#include <cstdint>
#include <Variant.h>
#include "filter.hpp"
#include "index.hpp"
#include "IntervalTree.h"
#include "vg.pb.h"
#include "fml.h"
#include "vg.hpp"
#include "gcsa.h"
#include "alignment.hpp"
#include "genotypekit.hpp"
using namespace std;
namespace vg{


    struct BREAKPOINT{
        string name;
        Position position;
        vector<BREAKPOINT> mates;

        string contig;
        int64_t start = -1;
        int64_t upper_bound = 100;
        int64_t lower_bound = 100;

        bool isForward;

        int fragl_supports = 0;
        int split_supports = 0;
        int other_supports = 0;

        inline int total_supports(){
            return fragl_supports + split_supports + other_supports;
        }
        inline bool overlap(BREAKPOINT p, int dist){

            if (start > -1 ){
                if ( abs(start - p.start) < dist){
                    return true;
                }
            }
            else{
                if (position.node_id() == p.position.node_id() && abs(position.offset() - p.position.offset()) < dist){
                    return true;
                }
            }
            
            return false;
        }
        inline string to_string(){
            stringstream x;
            x << "Pos: " << start << " u: " << upper_bound << " l: " << lower_bound << " s: " << total_supports();
            return x.str();
        }

    };



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
          
class DepthMap {
    /**
    *  Map <node_id : depth>
    *  or
    *  Map <node_id : offset : depth>
    */
public:
  int8_t* depths;
  uint64_t size;
  inline DepthMap(int64_t sz) { depths = new int8_t[sz]; };
  inline DepthMap() { depths = new int8_t[1000]; };
  inline int8_t get_depth(int64_t node_id) { return depths[node_id]; };
  inline void set_depth(int64_t node_id, int8_t d) { depths[node_id] = d; };
  inline void fill(vector<Alignment> alns){
    for (auto a : alns){
        for (int i = 0; i < a.path().mapping_size(); ++i){
            Mapping m = a.path().mapping(i);
            for (int j = 0; j < m.edit_size(); j++){
                if (m.edit(j).to_length() == m.edit(j).from_length() &&
                     m.edit(j).sequence().empty()){
                        depths[m.position().node_id()] += 1;
                        break;
                     }
            }
        }
    }
  };
  inline void fill(vector<Path> paths){
    #pragma omp parallel for
    for (int p_ind = 0; p_ind < paths.size(); ++p_ind){
        Path p = paths[p_ind];
        for (int i=0; i<p.mapping_size(); ++i){
            Mapping m = p.mapping(i);
            bool match = false;
            for (int j = 0; j < m.edit_size(); ++j){
                if (m.edit(j).from_length() == m.edit(j).to_length() &&
                        m.edit(j).sequence().empty()){
                            match = true;
                        }
            }
            if (match){
                #pragma omp atomic
                depths[m.position().node_id()] += 1;
            }
            
        }
    }
  };
  inline void fill(Path path){
        for (int i=0; i<path.mapping_size(); ++i){
            Mapping m = path.mapping(i);
            for (int j = 0; j < m.edit_size(); ++j){
                if (m.edit(j).from_length() == m.edit(j).to_length() &&
                        m.edit(j).sequence().empty()){
                            depths[m.position().node_id()] == 1;
                            break;
                        }
            }
            
        }
  };
  inline void fill(Mapping m){
    #pragma omp atomic
    depths[m.position().node_id()]++;
  };
  inline void fill(int64_t node_id) {
    if (node_id < size) {
        #pragma omp atomic
        depths[node_id] += 1;
    } else {
        cerr << "WARNING: INVALID NODE" << node_id << endl;
    }
  }
};


    class SRPE{

        
           

        public:
            vector<string> ref_names;

            vector<pair<int, int> > intervals;

            // Calculate a proxy for discordance between a set of Alginments
            // and a subgraph (e.g. one that's been modified with a candidate variant)
            // Useful for deciding which variant is closest to what's represented in reads
            double discordance_score(vector<Alignment> alns, VG* subgraph);

            // Convert Alignments to the read-like objects Fermi-lite uses in assembly
            void aln_to_bseq(Alignment& a, bseq1_t* read);

            // Assemble a set of alignments into a set of unitigs
            // UNITIGS ARE GRAPH ELEMENTS - you could make them subgraphs.
            // Alignments need not map to the graph (e.g. they could be unmapped reads)
            void assemble(vector<Alignment> alns, vector<fml_utg_t>& unitigs);

            // Assemble a set of Alignments that map along <refpath> between <startpos> and <endpos>,
            // which are reference-relative coordinates (a.k.a your standard, linear ref coordinates)
            void assemble(string refpath, int64_t start_pos, int64_t end_pos, vector<fml_utg_t>& unitigs);

            // Assemble all reads that overlap a given position (within window_size bp)
            void assemble(int64_t node_id, int64_t offset, int window_size);

            // Are multiple references present in the same subgraph?
            bool overlapping_refs = false;

            // Maps from node-id to read depth
            DepthMap depth;

            // Every SRPE gets its own filter
            vg::Filter ff;

            // Every SRPE also gets its own name->alignment map
            // and a name->mate map
            map<string, Alignment> name_to_aln;
            map<string, string> aln_to_mate;

            // A graph (or subgraph) for the region this SRPE is handling.
            vg::VG* graph;
            // xg::XG* xindex;
            // gcsa::GCSA* gindex;
            // gcsa::LCPArray * lcp_ind;

            // Cap the total coverage at a given position
            int max_reads = 125;



};
}
#endif

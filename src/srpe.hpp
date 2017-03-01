#ifndef VG_SRPE
#define VG_SRPE
#include <string>
#include <cstdint>
#include <Variant.h>
#include "filter.hpp"
#include "index.hpp"
#include "vg.pb.h"
#include "vg.hpp"
#include "gcsa.h"
#include "alignment.hpp"
#include "genotyper.hpp"
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

//  class BPDepthMap{
//     public:
//         // A whole bunch of pointers
//         // node_id -> bp -> int4_t
//         uint8_t** bp_depths;
//         uint64_t size;
//         inline BPDepthMap(int64_t sz) {bp_depths = new uint8_t**[sz]};

//  };

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
  inline void add_node(int64_t node_id) {
    if (node_id < size) {

    } else {

    }
  }
};

/**
 class SnarlIndex{
    map<Snarl, vector<string> > snarl_alleles;
    map<Snarl, vector<uint8_t> > allele_to_depth 
 };
 **/


    class SRPE{

            // Looks like locus_to_vcf, but handles the SVTYPE /SVLEN info fields
            // vcflib::Variant locus_to_sv_vcf(Locus ll);
            // vcflib::Variant locus_to_vcf(Locus ll);

            // void remap(vg::VG* graph, Index gam_index, vector<pair<Alignment, Alignment> >& remapped);

            // void filter(vector<Alignment>& in_alns, vector<Alignment>& out_alns);

            // void normalize_pairs(vector<pair<Alignment&, Alignment&> > alns);
            // void normalize_singles(vector<Alignment&> alns);
            // // Performs graph augmentation in a way that won't shatter the graph into single
            // // base pairs. Integrates one path per variant; does not integrate a path for
            // // perfect matches to the reference path.
            // void single_smart_augment(vector<Alignment&> alns);
            // void pair_smart_augment(vector<pair<Alignment&, Alignment&> > alns);
            // void make_depth_map(Alignment& aln);
            // // Write a signed 16-bit int for every node id in the graph.
            // // Every 16 bits is a node id. -1 signifies there is no node POSITION
            // // in the graph.
            // void export_depth_map(string outname);
            // void import_depth_map(string mapname);
            // void export_reference_names(vector<string> names);
            /**
    struct LongVariantAffinity{
        bool consistent_enough = false;
        double affinity = 0.0;
        bool is_reverse = false;
        double score = 0.0;
        double likelihood_ln = 0.0;
        LongVariantAffinity() = default;
        LongVariantAffinity(double affinity, bool is_reverse) : consistent_enough(affinity >= 0.95),
            affinity(affinity), is_reverse(is_reverse), score(affinity);
    };
            */

            // Remove ultrabubbles/superbubbles/snarls that represent simple edges
            // between adjacent nodes in the graph.
            /*void remove_node_junctions(vector<Snarls> snarls);
            //void SRPE::remove_node_junctions(vector<Snarls> snarls){
                vector<Snarls> ret;
                for (auto x : snarls){
                    if (outdegree(start) == 1 && indegree(end) == 1 && end - start == 1{
                        continue;
                    }
                    else{
                        ret.push_back(x);
                    }
                }
            }
            */
        public:
            vector<string> ref_names;
            //DepthMap depth;
            vg::Filter ff;
            Genotyper gg;
            
            // A graph (or subgraph) for the region this SRPE is handling.
            vg::VG* vg;


};
}
#endif

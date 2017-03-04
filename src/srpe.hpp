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
#include "genotypekit.hpp"
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
           
    
    struct LongVariantAffinity{
        bool consistent_enough = true;
        double affinity = 0.0;
        bool is_reverse = false;
        double score = 0.0;
        double likelihood_ln = 0.0;
        LongVariantAffinity() = default;
    };


struct Flags{
    /**
    * 1 = mapped
    * 2 = mate_unmapped
    * 3 = read is soft-clipped
    * 4 = read is split
    * 5 = discordant orientation
    * 6 = discordant insert length
    * 7 = 
    * 8 = 
    */
    bool* flags;
};

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

            // Looks like locus_to_vcf, but handles the SVTYPE /SVLEN info fields
            // vcflib::Variant locus_to_sv_vcf(Locus ll);
            // vcflib::Variant locus_to_vcf(Locus ll);

            // void restrict_to_interval(string pathname, uint64_t start, uint64_t end, vg::VG* graph);

            // clear any region restrictions, delete all the pointers, clear depth maps and settings.
            // void reset();

            // Remake the depth map for a graph / subgraph.
            // void regenerate_depth_map(Alignment& a);
            // void regenerate_depth_map(pair<Alignment, Alignment> alns);

            // void remap(vg::VG* graph, Index gam_index, vector<pair<Alignment, Alignment> >& remapped);

            // Use Filter to locate read-pair signatures
            // void detect_paired_signatures(vector<Alignment>& in_alns, vector<Alignment>& out_alns);

            // Use Filter to detect paired-read signatures
            // void detect_split_read_signatures(Alignment& aln);

            // std::pair<int,int> calculate_insert_size_bounds(vector<pair<Alignment, Alignment> > alns);
            // std::pair<double, double> calculate_insert_size_sigma(vector<pair<Alignment, Alignment> alns);

            // check if our insert size for a given set of alignments falls within our expected range.
            // bool reasonable_insert_size(pair<Alignment, Alignment> mates);

            // void normalize_pairs(vector<pair<Alignment&, Alignment&> > alns, vector<pair<Alignment, Alignment> > normals);
            // void normalize_singles(vector<Alignment&> alns, vector<Alignment&> normals);

            // void load_read(vector<Alignment> a);
            // void flush_reads(vector<Alignment> a);

            
            // // Performs graph augmentation in a way that won't shatter the graph into single
            // // base pairs. Integrates one path per variant; does not integrate a path for
            // // perfect matches to the reference path. Does not integrate unanchored paths.
            // void single_smart_augment(vector<Alignment&> alns);
            // void pair_smart_augment(vector<pair<Alignment&, Alignment&> > alns);


            // // Write a signed 16-bit int for every node id in the graph.
            // // Every 16 bits is a node id. -1 signifies there is no node POSITION
            // // in the graph.
            // void export_depth_map(string outname);
            // void import_depth_map(string mapname);
            // void export_reference_names(vector<string> names);
           

        public:
            vector<string> ref_names;
            // Are multiple references present in the same subgraph?
            bool overlapping_refs = false;
            DepthMap depth;
            vg::Filter ff;

            // Mapper mapper;

            // A graph (or subgraph) for the region this SRPE is handling.
            vg::VG* vg;

            // Cap the total coverage at a given position
            // int max_reads = 125;


};
}
#endif

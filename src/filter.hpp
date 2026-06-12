#ifndef VG_FILTER
#define VG_FILTER

#include <vector>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <string>
#include <math.h>
#include "vg.hpp"
#include "mapper.hpp"
#include <vg/vg.pb.h>

/** \file
 * Provides a way to filter Edits contained
 * within Alignments. This can be used to clean out
 * sequencing errors and to find high-quality candidates
 * for variant calling.
 *
 */
namespace vg{

struct BREAKPOINT{
    string name;
    Position position;
    vector<BREAKPOINT> mates;
    
    string contig;
    int64_t start = -1;
    int64_t upper_bound = 100;
    int64_t lower_bound = 100;
    
    // Does the breakpoint point this way --->>
    // or this way <<---
    bool isForward;
    // 0: Unset, 1: INS, 2: DEL, 3: INV, 4: DUP
    int SV_TYPE = 0;
    //
    
    int normal_supports = 0;
    int tumor_supports = 0;
    
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

class Filter{
    public:
        Filter();
        ~Filter();
         /* Filter functions.
         * Take an Alignment and walk it.
         * Most of these were used by the now-removed "vg sift"
         * What remains is used by "vg genotype"
         */
        bool anchored_filter(Alignment& aln);
        bool mark_sv_alignments(Alignment& a, Alignment& b);
        bool mark_smallVariant_alignments(Alignment& a, Alignment& b);
       
        bool soft_clip_filter(Alignment& aln);
        bool unmapped_filter(Alignment& aln);
        bool split_read_filter(Alignment& aln);

        /*PE Functions*/
        bool one_end_anchored_filter(Alignment& aln_first, Alignment& aln_second);
        bool interchromosomal_filter(Alignment& aln_first, Alignment& aln_second);
        
        // // TODO should give this one an insert size arg
        bool insert_size_filter(Alignment& aln_first, Alignment& aln_second);
        bool pair_orientation_filter(Alignment& aln_first, Alignment& aln_second);

    public:
        int soft_clip_limit = -1;
        float insert_mean = 1000;
        float insert_sd = 100;
        };
}

#endif

#include "ggsv.hpp"

GGSV::GGSV(){

}

GGSV::GGSV(string indexfile){
    
}

GGSV::~GGSV(){

}

void GGSV::set_index(string indexfile){

}


vector<StructuralVariant> GGSV::alignment_to_known_sv(Alignment aln){
    vector<StructuralVariant> ret;

}
    
vector<StructuralVariant> GGSV::alignment_to_putative_sv(Alignment aln){
    vector<StructuralVariant> ret;

}


/**
 * Take in an alignment
 * scan it for evidence of deletions, 
 * including:
 *  - softclips
 *  - 
 */
vector<Deletion> GGSV::alignment_to_deletion(Alignment& aln){
    vector<Deletion> ret;
    
}

//CALL INDELS
void GGSV::call_split_read(string gamfile){
    
    // Open gam file
    

    //For each alignment in gam file:
    

        // if there's a split read (as defined by the filter),
        // then add a split read to a map<WigglePosition, SV>

        // If a given WigglePosition contains a sufficient number of SV calls,
        // consider it a legitimate SV (a deletion or an insertion

}

//CALL INVERSIONS
void GGSV::call_inversion(string gamfile){

    /**
     * Iterate through Alignments in the gam file
     * and place putative breakpoints in a global map
     * at their WigglePosition if the alignments support
     * such a breakpoint
     *
     * Once enough have been accumulated, call the SV at the
     * given position.
     */

}

//CALL INVERSIONS
void GGSV::call_duplications(string gamfile){

}


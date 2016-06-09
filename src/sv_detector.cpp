#include "sv_detector.hpp"

using namespace vg;
using namespace std;
SV_DETECTOR::SV_DETECTOR(){

}

SV_DETECTOR::SV_DETECTOR(VG* g){
    my_graph = g;
}


SV_DETECTOR::~SV_DETECTOR(){

}



vector<StructuralVariant> SV_DETECTOR::alignment_to_putative_sv(Alignment aln){
    vector<StructuralVariant> ret;

}
/**
void SV_DETECTOR::cache_paths_of_graph(){
    regex is_alt("_alt_.+_[0-9]+");
    for (auto p : my_graph->paths._paths){
        if (regex_match(p.first, is_alt)){
            name_to_paths[p.first] = p.second;
            for (auto m : p.second){
                vg::id_t n_id = m.position().node_id();
                node_id_to_path[n_id] = p.first;
            }
        }
    }
}
*/

vector<StructuralVariant> SV_DETECTOR::gam_to_known_sv(string gamfile){
    vector<StructuralVariant> ret;
    
    vector<string> matches;
    std::function<void(Alignment&)> detect = [&, &matches](Alignment& aln){
        Path path = aln.path();
        for (int i = 0; i < path.mapping_size(); i++){
            Position x = path.mapping(i).position();
            vg::id_t node_id = x.node_id();
            map<vg::id_t, string>::iterator t_path = node_id_to_path.find(node_id);
            if (t_path != node_id_to_path.end()){
                // We have a node_id which supports a known SV,
                // But we haven't determined yet if it's ref or alt
                // or put it in a count...
            }
        }
    };

    if (gamfile == "-"){
        stream::for_each_parallel(cin, detect);
    }   
    else{
        ifstream in;
        if (in.good()){
            stream::for_each_parallel(in, detect);
        }
    }
}

vector<string> SV_DETECTOR::alignment_to_known_sv(Alignment aln){

}


/**
 * Take in an alignment
 * scan it for evidence of deletions, 
 * including:
 *  - softclips
 *  - 

 vector<Deletion> SV_DETECTOR::alignment_to_deletion(Alignment& aln){
 vector<Deletion> ret;

 }

//CALL INDELS
void SV_DETECTOR::call_split_read(string gamfile){

// Open gam file


//For each alignment in gam file:


// if there's a split read (as defined by the filter),
// then add a split read to a map<WigglePosition, SV>

// If a given WigglePosition contains a sufficient number of SV calls,
// consider it a legitimate SV (a deletion or an insertion

}

//CALL INVERSIONS
void SV_DETECTOR::call_inversion(string gamfile){

/**
 * Iterate through Alignments in the gam file
 * and place putative breakpoints in a global map
 * at their WigglePosition if the alignments support
 * such a breakpoint
 *
 * Once enough have been accumulated, call the SV at the
 * given position.
 *

 }

//CALL INVERSIONS
void SV_DETECTOR::call_duplications(string gamfile){

}
*/

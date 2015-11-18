#include "deconstructor.hpp"

using namespace std;

namespace vg {
    Deconstructor::Deconstructor() : index(nullptr), reference(NULL), xg(NULL){
    }

    Deconstructor::~Deconstructor(){
        clear();
    }

    void Deconstructor::clear(){

    }

    void Deconstructor::set_reference(string ref_file){
      if (!ref_file.empty()){
        FastaReference fr;
        fr.open(ref_file);
        reference = &fr;
      }
      else{
        cerr << "Error: fasta file not provided" << endl;
        exit(1);
      }
    }

    /**
     * Set the deconstructor's index, be it a rocksdb-backed disk index
     * or an XG + gcsa index
     */
    void Deconstructor::set_index(Index* i){
      index = i;
    }

    void Deconstructor::set_xg(xg::XG* x){
      xg = x;
    }

    void Deconstructor::enumerate_paths_in_index(){
        // TODO Load index

        // TODO Build a vector of reference paths
        // from a pathfile or fasta file

        // TODO Lookup paths in index and return the name if present


    }

    /**
     * This one is still a bit unsettled in implementation, but the gist
     * is to:
     * 1. Open the index or xg
     * 2. One path represents the reference and the other represents a variant
     *    - BFS along the reference path
     *    - at each x-furcation, process the different variant path(s) into a variant
     *    record.
     */
    Path Deconstructor::relative_mapping(Path& p1, Path& p2){
        VG graph;
        Path ret;

    }

    Path Deconstructor::relative_mapping(Path& p1, string p2){
        VG graph;
        Path ret;
    }

    vcflib::Variant Deconstructor::mapping_to_variant(Path variant, Path ref){

    }

    vcflib::VariantCallFile Deconstructor::write_variants(string filename, vector<vcflib::Variant> variants){

    }

}

#include "deconstructor.hpp"

using namespace std;

namespace vg {
    Deconstructor::Deconstructor() : index_file(""), reference(""), xg_file(""), vgraph(nullptr) {
    }

    Deconstructor::~Deconstructor(){
        clear();
    }

    void Deconstructor::clear(){

    }

    void Deconstructor::set_graph(VG* v){
      vgraph = v;
    }

    void Deconstructor::set_reference(string ref_file){
      cerr << "Setting reference to " << ref_file << endl;
      // if (!ref_file.empty()){
      //   FastaReference fr;
      //   fr.open(ref_file);
      //   for (auto seq : fr.index->sequenceNames){
      //     cout << seq << endl;
      //   }
      //   reference = &fr;
      // }
      // else{
      //   cerr << "Error: fasta file not provided" << endl;
      //   exit(1);
      // }
      reference = ref_file;
    }

    /**
     * Set the deconstructor's index, be it a rocksdb-backed disk index
     * or an XG + gcsa index
     */
    void Deconstructor::set_index(string i){
      cerr << "Setting index to " << i << "." << endl;
      index_file = i;
    }

    void Deconstructor::set_xg(string x){
      cerr << "Setting XG index to " << x << "." << endl;
      xg_file = x;
    }

    void Deconstructor::enumerate_path_names_in_index(){
        FastaReference fr;
        fr.open(reference);
        for (auto& seq : fr.index->sequenceNames){
          cout << seq << endl;
        }
        // TODO Load index
        if (this->xg_file != ""){

        }
        else if (this->index_file != ""){
          Index ind;
          ind.open_read_only(index_file);
          auto path_ids = ind.paths_by_id();
          for (auto& id : path_ids){
            cerr << id << endl;
          }
        }
        else{
          cerr << "An XG- or rocksdb-index must be provided for deconstruction." << endl;
          exit(1);
        }
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
        Path ret;
        int steps_back = 0;
        int max_steps = 1000;
        VG g = *vgraph;

    }

    Path Deconstructor::relative_mapping(Path& p1, string path_name2){
        if ((*vgraph).paths.has_path(path_name2)){
            cerr << "Path: " << path_name2 << " in graph " << endl;
        }
       list<Mapping> p2_mapping = (*vgraph).paths.get_path(path_name2);
       Path p2;
        return relative_mapping(p1, p2);
    }

    vcflib::Variant Deconstructor::path_to_variant(Path variant, Path ref){

    }

    vcflib::Variant Deconstructor::pathname_to_variant(string variant, Path ref){

    }

    vcflib::VariantCallFile Deconstructor::write_variants(string filename, vector<vcflib::Variant> variants){

    }

}

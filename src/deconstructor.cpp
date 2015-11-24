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

    /**
    * A drop-in replacement for index.surject_alignment that uses an XG index.
    */
    // bool Deconstructor::surject_alignment(Alignment& source,
    //                         set<string> path_names,
    //                         Alignment& surjection,
    //                         string& path_name,
    //                         int64_t& path_pos,
    //                         int window){
    //   VG graph;
    //
    //   if (!source.has_path() || source.path().mapping_size() == 0) {
    //       return false;
    //   }
    //   int64_t from_id = source.path().mapping(0).position().node_id() - window;
    //   int64_t to_id = source.path().mapping(source.path().mapping_size()-1).position().node_id() + window;
    //   //TODO Not implemented get_range(max((int64_t)0, from_id), to_id, graph);
    //   graph.remove_orphan_edges();
    //   set<string> kept_paths;
    //   graph.keep_paths(path_names, kept_paths);
    //   surjection = source;
    //   surjection.clear_path();
    //   graph.align(surjection);
    //
    //   if (surjection.path().mapping_size() > 0 && kept_paths.size() == 1) {
    //     // determine the paths of the node we mapped into
    //     //  ... get the id of the first node, get the pahs of it
    //     assert(kept_paths.size() == 1);
    //     path_name = *kept_paths.begin();
    //
    //     int64_t path_id = get_path_id(path_name);
    //     int64_t hit_id = surjection.path().mapping(0).position().node_id();
    //     Node hit_node;
    //     get_node(hit_id, hit_node);
    //     bool hit_backward = surjection.path().mapping(0).position().is_reverse();
    //     int64_t pos = surjection.path().mapping(0).position().offset();
    //     // we pick up positional information using the index
    //     int64_t prev_pos=0, next_pos=0;
    //     bool prev_orientation, next_orientation;
    //     list<pair<int64_t, bool>> path_prev, path_next;
    //     //TODO NOT IMPLEMENTED get_node_path_relative_position(hit_id, hit_backward, path_id, path_prev, prev_pos, prev_orientation,
    //     //                                path_next, next_pos, next_orientation);
    //     // as a surjection this node must be in the path
    //     // the nearest place we map should be the head of this node
    //     assert(path_prev.size()==1 && hit_id == path_prev.front().first && hit_backward == path_prev.front().second);
    //     // we then add the position of this node against the path to our offset in the node (from the left side)
    //     // to get the final chrom+position for the surjection
    //     path_pos = prev_pos + (hit_backward ? hit_node.sequence().size() - pos - 1 : pos);
    //     // we need the cigar, but this comes from a function on the alignment itself
    //     return true;
    // } else {
    //     return false;
    // }
    //
    // }

    // void Deconstructor::get_range_in_xg(int64_t from_id, int64_t to_id, VG& graph){
    //
    // }
    //
    // void Deconstructor::get_node_in_xg(int64_t node_id, Node& node){
    //
    // }

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
        vector<string> paths_in_ref;
        map<string, int64_t> intersection_ref_and_index;
        for (auto& seq : fr.index->sequenceNames){
          //cout << seq << endl;
          paths_in_ref.push_back(seq);
        }
        // TODO Same logic for index file applies to xg, roughly speaking
        if (this->xg_file != ""){

        }
        else if (this->index_file != ""){
          Index ind;
          ind.open_read_only(index_file);
          map<string, int64_t> path_ids = ind.paths_by_id();
          map<string, int64_t> ::iterator it;
          for (it = path_ids.begin(); it != path_ids.end(); it++){
            //cerr << it->first << endl;
            intersection_ref_and_index[it->first] = it->second;
          }
        }
        else{
          cerr << "An XG- or rocksdb-index must be provided for deconstruction." << endl;
          exit(1);
        }

        ref_paths = paths_in_ref;
        inter_ref_and_index = intersection_ref_and_index;


    }

    vector<vcflib::Variant> Deconstructor::get_variants(string region_name = ""){
      vector<vcflib::Variant> variants;

      vector<string> paths_to_project;
      if (region_name != ""){
        paths_to_project = vector<string>();
        //TODO check if region in reference paths TODO
        paths_to_project.push_back(region_name);
      }
      else{
        map<string, int64_t>::iterator it;
        for (it = inter_ref_and_index.begin(); it != inter_ref_and_index.end(); it++){
          paths_to_project.push_back(it->first);
        }
      }

      //#pragma omp parallel for
      int i;
      vcflib::Variant v;
      for (i = 0; i < paths_to_project.size(); i++){
        // v = pathname_to_variants(paths_to_project[i]);
        cerr << "Extracting variant paths for " << paths_to_project[i] << endl;
        //for ()
        //#pragma omp critical
        variants.push_back(v);
      }


      return variants;

    }

    /**
     * This one is still a bit unsettled in implementation, but the gist
     * is to:
     * 1. Open the index or xg
     * 2. One path represents the reference and the other represents a variant
     *    - BFS along the reference path
     *    - at each x-furcation, process the different variant path(s) into a variant
     *    record (essentially an alignment against the reference)
     */
    Path Deconstructor::relative_mapping(Path& p1, Path& p2){
        Path ret;
        Index vindex;
        vindex.open_read_only(index_file);
        //(*index).get_nearest_node_next_path_member();
        //(*index).get_nearest_node_prev_path_member();
        //(*index).path_relative_mapping()

        pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> node_and_path_taken_prev;
        pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> node_and_path_taken_next;
        // TODO Get previous node in path1
        // could also grab the Position from the first/last mappings
        const Position p = path_start(p1);
        pos_t p_pos_t = make_pos_t(p);
        cerr << get_id(p_pos_t) << endl;
        int64_t p_node_id = get_id(p_pos_t);
        // Position q = path_start(p1);
        // int64_t q_node_id = get_id(q);
        // int64_t future_path_pos;
        // bool future_rel_orientation;
        // node_and_path_taken_prev = vindex.get_nearest_node_prev_path_member(p_node_id, is_rev(p), vindex.get_path_id(p2.name()),
        //                             future_path_pos, future_rel_orientation, 4);
        // node_and_path_taken_next = vindex.get_nearest_node_next_path_member(q_node_id, is_rev(q), vindex.get_path_id(p2.name()),
        //                             future_path_pos, future_rel_orientation, 4);

        return ret;
    }

    Path Deconstructor::relative_mapping(Path& p1, string path_name2){
        if ((*vgraph).paths.has_path(path_name2)){
            cerr << "Path: " << path_name2 << " in graph " << endl;
        }
       list<Mapping> p2_mapping = (*vgraph).paths.get_path(path_name2);
       Path p2;
        return relative_mapping(p1, p2);
    }

    vcflib::Variant Deconstructor::path_to_variants(Path variant, Path ref){

    }

    vcflib::Variant Deconstructor::pathname_to_variants(string variant, Path ref){

    }

    void Deconstructor::write_variants(string filename, vector<vcflib::Variant> variants){
        cerr << "writing variants." << endl;
        vcflib::VariantCallFile v;

    }

}

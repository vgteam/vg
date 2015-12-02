#include "deconstructor.hpp"

using namespace std;

namespace vg {
    Deconstructor::Deconstructor() : index_file(""), reference(""), xg_file(""), vgraph(nullptr) {
    }

    Deconstructor::~Deconstructor(){
        clear();
    }

    void Deconstructor::clear(){
      vgraph = nullptr;
      index_file = "";
      reference = "";
      xg_file = "";

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

    void Deconstructor::enumerate_graph_paths(){

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

    // TODO
    vector<vcflib::Variant> Deconstructor::get_variants(string region_file){
      vector<vcflib::Variant> vars;
      cerr << "Not implemented" << endl;
      exit(1);
      return vars;
    }

    /**
    * If a region name is given, extract variants for that path in the reference.
    * Otherwise, try to get variants for all paths. If a pathname is not in the
    * reference, skip it.
    * First, locate list<Mapping> for nodes that connect to reference paths but
    * that lie off of them.
    * Next, transform this information into a vcf record, append it to a vector,
    * and return.
    *
    */
    vector<vcflib::Variant> Deconstructor::get_variants(string region_name, int start, int end){
      vector<vcflib::Variant> variants;
      //This function must be called, as it takes the intersection of paths in
      // the reference and index, which is used a few lines later.
      enumerate_path_names_in_index();
      vector<string> paths_to_project;
      if (region_name != ""){
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
      int j;
      vcflib::Variant v;
      for (i = 0; i < paths_to_project.size(); i++){
        // Get a mapping for each node in the graph that lies attached, but not on,
        // a reference path.
        string p = paths_to_project.at(i);
        list<Mapping> m = get_mappings_off_reference(p);
        cerr << "Retrieved mappings of non-reference nodes. Converting to VCF..." << endl;
        cerr << "There are " << m.size() << " mappings to convert." << endl;
        list<Mapping>::iterator it;
        for (it = m.begin(); it != m.end(); it++){
          v = mapping_to_variant(*it);
          variants.push_back(v);
        }
      }
      return variants;

    }

    list<Mapping> Deconstructor::get_mappings_off_reference(string pathname){
      Index vindex;
      //vindex.open_read_only(index_file);
      //map<string, int64_t> path_ids = vindex.paths_by_id();
      //int64_t path_id = path_ids.at(pathname);
      Path ref = (*vgraph).paths.path(pathname);
      //vindex.close();
      list<Mapping> m = get_mappings_off_reference(ref);
      return m;
    }

    /**
     * Returns a path, which is a list of mappings.
     * A mapping is a list of edits that transform a position
     * within a path to another path at that same position.
     *
     * So, for each position in the returns list of mappings, it's
     * possible to reconstruct the variant at that position by transforming
     * the edits.
     */
     list<Mapping> Deconstructor::get_mappings_off_reference(Path& ref){
        Index ix;
        ix.open_read_only(index_file);
        std::list<Mapping> mapping_list = std::list<Mapping>();
        (*vgraph).for_each_node([this, &mapping_list, ref, &ix](Node* n) mutable {
            //cerr << n->id() << endl;
            if ((*vgraph).paths.has_node_mapping(n->id())){
              //cerr << "In mapping" << endl;
            }
            else{
              //pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> (*index)
              //Index ix;
              //ix.open_read_only(index_file);
              map<string, int64_t> paths = ix.paths_by_id();
              bool backward = false;

              // Alright, here's the meat of it.
              // Get the previous node in the path and the next node in the path
              // Then get a relative mapping using these.
              // Then devise a way to translate this to a vcf record.
              pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> prev_node_in_named_path;
              pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> next_node_in_named_path;
              int64_t path_pos;
              bool rel;
              Mapping m; //mapping = {edits} + position, where position P is
              //

              auto ref_id = paths.at(ref.name());
              //auto alt_id = paths.at(alt.name());

              prev_node_in_named_path = ix.get_nearest_node_prev_path_member((int64_t) n->id(), backward,
                                                ref_id, path_pos, rel, 4);
              next_node_in_named_path = ix.get_nearest_node_next_path_member((int64_t) n->id(), backward,
                                                ref_id, path_pos,rel, 4);
              m = ix.path_relative_mapping((int64_t) n->id(), backward, ref_id,
                                        prev_node_in_named_path.first, prev_node_in_named_path.second.first, prev_node_in_named_path.second.second,
                                        next_node_in_named_path.first, next_node_in_named_path.second.first, next_node_in_named_path.second.second);
              mapping_list.push_back(m);
            }
    });

    return mapping_list;

  }



    /**
    * Transforms a Mapping (a list of edits on a path) to a vcf entry.
    * The reference path that the variant originates from can be grabbed
    * from the Mapping's position.
    */
    vcflib::Variant Deconstructor::mapping_to_variant(Mapping m){
      vcflib::Variant v;
      int64_t n_id = (int64_t) m.position().node_id();
      Node* n = (*vgraph).get_node(n_id);
      const string x = mapping_sequence(m, *n);
      cerr << x << " " << m.edit(0).sequence() << endl;


      return v;
    }



    // vcflib::Variant Deconstructor::pathname_to_variants(string ref_path){
    //   Path alt;
    //
    //   //list<Mapping> mapping_list = relative_mapping()
    // }

    void Deconstructor::write_variants(string filename, vector<vcflib::Variant> variants){
        cerr << "writing variants." << endl;
        vcflib::VariantCallFile v;

    }

}

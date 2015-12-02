#include "deconstructor.hpp"

using namespace std;

namespace vg {
    Deconstructor::Deconstructor() : index_file(""), reference(""), xg_file(""), vgraph(nullptr) {
    }

    Deconstructor::~Deconstructor(){
        clear();
    }


    /**
     * Nulls out the important class members.
     * Probably extraneous for destruction but might be useful
     * if we ever wanted to make this class a singleton.
     */
    void Deconstructor::clear(){
      vgraph = nullptr;
      index_file = "";
      reference = "";
      xg_file = "";

    }

    /**
     * Stores a pointer to the graph so that we can access things
     * like paths_by_id. 
     */
    void Deconstructor::set_graph(VG* v){
      vgraph = v;
    }

    /**
     * Stores the name of the reference fasta file.
     * The reference file itself is opened with FastaHack.
     */
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


    /**
     * Sets the name of the xg index file. Right now this isn't used TODO
     * but in the future it will behave much like Index currently does.
     */
    void Deconstructor::set_xg(string x){
      cerr << "Setting XG index to " << x << "." << endl;
      xg_file = x;
    }

    /**
     * This function takes in the reference file and the index.
     * It then compares them to see which paths are in both.
     * This way, we never try grabbing a path that's in the
     * reference but not the index.
     */
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
    /**
     *
     */
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
      // Create a list of 
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
        std::list<Mapping> mapping_list;
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
              bool backward = true;

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
      //cerr << (*n).sequence();
      cerr << (int64_t) m.position().node_id() << " " << endl;
      cerr << m.edit(0).sequence() << " " << m.edit(0).from_length() << " " << m.edit(0).to_length() << endl;
      //cerr << x << " " << m.position().node_id() << endl;
      //int64_t prenode = (((*vgraph).edges_on_start[ (int64_t) m.position().node_id()])[0]).first;


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

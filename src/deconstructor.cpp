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
    // Just argument handling here - either take the single region given
    // or the entire intersection of the reference and index
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

    int i;
    int j;
    vcflib::Variant v;
    for (i = 0; i < paths_to_project.size(); i++){
      // Get a mapping for each node in the graph that lies attached, but not on,
      // a reference path.
      string p = paths_to_project.at(i);
      //list<Mapping> m = get_mappings_off_reference(p);
      vector<int64_t> variant_node_ids = get_variant_node_ids(p);
      int j;
      //#pragma omp parallel for
      for (j = 0; j < variant_node_ids.size(); j++){
        cerr << variant_node_ids[j] << endl;
        //Mapping m = node_id_to_mapping(variant_node_ids[j]);
        //v = mapping_to_simple_variant(m, variant_node_ids[j]);
        //#pragma omp critical
        //variants.push_back(v);
      }
      cerr << "Retrieved mappings of non-reference nodes. Converting to VCF..." << endl;
      //cerr << "There are " << m_list.size() << " mappings to convert." << endl;
      //list<Mapping>::iterator it;
      //for (it = m_list.begin(); it != m_list.end(); it++){
      //v = mapping_to_variant(*it);
      //  variants.push_back(v);
      //}
    }
    return variants;

  }

  vector<int64_t> Deconstructor::get_variant_node_ids(string pathname){
    Path path = (*vgraph).paths.path(pathname);
    vector<int64_t> ret;

    Index ix;
    ix.open_read_only(index_file);

    int i;
    list<Mapping> m_list = (*vgraph).paths._paths[pathname];
    list<Mapping>::iterator it;
    for (it = m_list.begin(); it != m_list.end(); it++){
      //get the node id of the mapping we've retrieved.
      int64_t n_id = (*it).position().node_id();
      bool orientation = false;
      // get the node's neighbors and check if they are path members
      vector<pair<int64_t, bool>> destinations;
      ix.get_nodes_next(n_id, orientation, destinations);
      int j;
      for (j = 0; j < destinations.size(); j++){
        int64_t neighbor_id = destinations[j].first;
        if (! (*vgraph).paths.has_node_mapping(neighbor_id)){
          ret.push_back(neighbor_id);
        }
        else{
          continue;
        }
      }
    }
    return ret;
  }

  void Deconstructor::print_using_edges(string pathname){
    Path path = (*vgraph).paths.path(pathname);
    Index ix;
    ix.open_read_only(index_file);
    bool backward = false;
    int64_t path_id = ix.paths_by_id().at(pathname);

    //int i;
    int j;
    list<Mapping> m_list = (*vgraph).paths._paths[pathname];
    list<Mapping>::iterator it;
    for(it = m_list.begin(); it != m_list.end(); it++){
      int64_t current_node_id = (*it).position().node_id();
      vector<Edge> edges_ahead;
      if(backward) {
        // "next" = right = start
        ix.get_edges_on_start(current_node_id, edges_ahead);
      } else {
        // "next" = right = end
        ix.get_edges_on_end(current_node_id, edges_ahead);
      }
      //cerr << current_node_id << endl;

      // base case: no variation, reference has simply been split up.
      // We ignore self-cycling as a possibility.
      if (edges_ahead.size() == 1){
         continue;
      }
      // Otherwise, we have found an n-furcation of the graph,
      // indicating that variation has been inserted.
      else{
        // We need to cross the variation simultaneously along all paths
        // to reach the next common node.
        //
        // TODO this should really use recursive backtracking but I'm not feeling that quite yet.

        //update edges.
        vector<pair<list<pair<int64_t, bool>>, pair<int64_t, bool>>> paths_to_common_node_on_path (edges_ahead.size());
        stack<int64_t> farthest;
        // int max_iters = 8;
        // int steps = 0;
        //while (steps < max_iters){
          // for (int edge_ind = 0; edge_ind < edges_ahead.size(); edge_ind++){
          //   Edge e = edges_ahead[edge_ind];
          //   int64_t next_n_id = e.to();
          //   Node n;
          //   vector<Edge> n_edges_ahead;
          //   ix.get_node(next_n_id, n);
          //
          //   int64_t path_pos; Mapping m; bool n_backward;
          //   if (ix.get_node_path(next_n_id, path_id, path_pos, n_backward, m) <= 0){
          //     if(!n_backward) {
          //       // "next" = right = start
          //       ix.get_edges_on_start(next_n_id, n_edges_ahead);
          //     } else {
          //       // "next" = right = end
          //       ix.get_edges_on_end(next_n_id, n_edges_ahead);
          //     }
          //   }
          //   if (n_edges_ahead.size() > 1){
          //     farthest.push(next_n_id);
          //     cerr << next_n_id << endl;
          //   }
          // }
          //steps++;
        //}
        list<Mapping> alts;
        Mapping ref;
        for (int edge_ind = 0; edge_ind < edges_ahead.size(); edge_ind++){
            // base case: from_node for each path is the same
            // and to_node for each path is the same.
            Edge e = edges_ahead[edge_ind];
            int64_t next_n_id = e.to();
            Node n;
            ix.get_node(next_n_id, n);

            // Node is not part of graph yet, so it must represent a variant
            if (!(*vgraph).paths.has_node_mapping(next_n_id)){
              pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> next_on_path;
              pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> prev_on_path;
              int64_t path_pos;
              bool rel;
              next_on_path = ix.get_nearest_node_next_path_member(next_n_id, backward,
                              path_id, path_pos, rel, 4);
              prev_on_path = ix.get_nearest_node_prev_path_member(next_n_id, backward,
                              path_id, path_pos, rel, 4);
              Mapping m = ix.path_relative_mapping(next_n_id, backward, path_id,
              prev_on_path.first, prev_on_path.second.first, prev_on_path.second.second,
              next_on_path.first, next_on_path.second.first, next_on_path.second.second);
              //cerr << "Alt node: " << m.position().node_id() << " " << m.edit(0).sequence() << endl;
              alts.push_back(m);
            }
            else {
              vector<Edge> n_e_in;
              vector<Edge> n_e_out;
              ix.get_edges_on_start(next_n_id, n_e_in);
              ix.get_edges_on_end(next_n_id, n_e_out);
              //if ((n_e_in.size() == 1 && n_e_out.size() == 1)){
                // Reference equivalent of a snp
                pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> next_on_path;
                pair<list<pair<int64_t, bool>>, pair<int64_t, bool>> prev_on_path;
                int64_t path_pos;
                bool rel;
                next_on_path = ix.get_nearest_node_next_path_member(next_n_id, backward,
                                path_id, path_pos, rel, 4);
                prev_on_path = ix.get_nearest_node_prev_path_member(next_n_id, backward,
                                path_id, path_pos, rel, 4);
                Mapping m = ix.path_relative_mapping(next_n_id, backward, path_id,
                prev_on_path.first, prev_on_path.second.first, prev_on_path.second.second,
                next_on_path.first, next_on_path.second.first, next_on_path.second.second);
                ref = m;
                //TODO currently captures the correct reference node but yields an incorrect
                //sequence.
                Node ref_seq_n; ix.get_node((int64_t) m.position().node_id(), ref_seq_n);
                //cerr << "Ref node: " << m.position().node_id() << " " << ref_seq_n.sequence() << endl;
              // }
              // else if (n_e_in.size() > 1 && n_e_out.size() == 1){
              //   // This represents a deletion TODO I think?
              //   cerr << "Deletion case " << next_n_id << endl;
              // }
              // else{
              //   cerr << "Highly complex case. " << next_n_id << endl;
              // }
            }
        }
        if (alts.size() > 0){
          Node ii;
          ix.get_node((int64_t) ref.position().node_id(), ii);
          cerr << "CHROM: " << pathname << " Pos: " << ref.position().node_id() <<
          " REF: " << ii.sequence() << " Alt: " << alts.front().edit(0).sequence() << endl;
        }

    }
    }
  }

  // stack<int64_t> Deconstructor::find_next_path_branch(int64_t c_id, Index ix, string pathname, int steps=4){
  //   int64_t path_id = ix.paths_by_id().at(pathname);
  //   if(backward) {
  //     // "next" = right = start
  //     ix.get_edges_on_end(c_id, edges_ahead);
  //   } else {
  //     // "next" = right = end
  //     ix.get_edges_on_start(c_id, edges_ahead);
  //   }
  //
  // }



  Mapping Deconstructor::node_id_to_mapping(int64_t alt_id){
    Mapping ret;
    return ret;

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
        //cerr << "prev node: " << prev_node_in_named_path.second.first << endl;
        //cerr << "current node: " << n->id() << " and its sequence: " << n->sequence() << endl;
        //cerr << "next node: " << next_node_in_named_path.second.first << endl;
        m = ix.path_relative_mapping((int64_t) n->id(), backward, ref_id,
        prev_node_in_named_path.first, prev_node_in_named_path.second.first, prev_node_in_named_path.second.second,
        next_node_in_named_path.first, next_node_in_named_path.second.first, next_node_in_named_path.second.second);
        //cerr << "Mapping node pos: " << m.position().node_id() << " and edit size: " << m.edit().size() << " " << m.edit(0).sequence() << endl;
        //cerr << "Edit lens; from: " << m.edit(0).from_length() << " and to: " << m.edit(0).to_length() << endl;
        mapping_list.push_back(m);
      }
    });

    int i;
    for (int i = 0; i < mapping_list.size(); i++){
      cerr << "Maps on maps." << endl;
    }

    return mapping_list;

  }



  /**
  * Transforms a Mapping (a list of edits on a path) to a vcf entry.
  * Currently, the mapping contains the reference allele.
  * The reference path that the variant originates from can be grabbed
  * from the Mapping's position.
  */
  vcflib::Variant Deconstructor::mapping_to_simple_variant(Mapping m, int64_t alt_id){
    vcflib::Variant v;
    int64_t n_id = (int64_t) m.position().node_id();
    Node* n = (*vgraph).get_node(n_id);
    const string x = mapping_sequence(m, *n);
    // Need: ref, alt
    // Quality, Genotype, etc
    //cerr << (*n).sequence();
    //cerr << (int64_t) m.position().node_id() << " " << endl;
    //cerr << m.edit(0).sequence() << " " << m.edit(0).from_length() << " " << m.edit(0).to_length() << endl;
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

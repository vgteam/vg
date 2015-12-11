#include "deconstructor.hpp"

using namespace std;

namespace vg {
    Deconstructor::Deconstructor() : ref_index(""), reference(""), xg_file(""), vgraph(nullptr) {
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
        ref_index = "";
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
        ref_index = i;
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
        else if (this->ref_index != ""){
            Index ind;
            ind.open_read_only(ref_index);
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



void get_variants_using_edges_from_file(string pathfile){

}

void Deconstructor::b_call(string pathname){
  vector<string> paths_to_project;
  if (pathname != ""){
      //TODO check if region in reference paths TODO
      paths_to_project.push_back(pathname);
  }
  else{
      map<string, int64_t>::iterator it;
      for (it = inter_ref_and_index.begin(); it != inter_ref_and_index.end(); it++){
          paths_to_project.push_back(it->first);
      }}
      for (auto pathname : paths_to_project){
          Path path = (*vgraph).paths.path(pathname);
          Index ix;
          ix.open_read_only(ref_index);
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

              // base case: no variation, reference has simply been split up.
              // We ignore self-cycling as a possibility.
              if (edges_ahead.size() <= 1){
                  continue;
              }
              // Otherwise, we have found an n-furcation of the graph,
              // indicating that variation has been inserted.
              else{
                Node current;
                ix.get_node(current_node_id, current);
                //cerr << "busted here " << endl;
                Node n = get_anchor_node(current, path_id);
                cerr << "Anchor: " << current_node_id << " Bow: " << n.id() << endl;
                  /**
                   * This should catch multi-allelic snps as-is,
                   * but it seems to be a bit off for INDELS.
                   */

                  // if (alts.size() > 0){
                  //     Node ii;
                  //     ix.get_node((int64_t) ref.position().node_id(), ii);
                  //     //TODO Need to cut sequence using node offset.
                  //     cerr << "CHROM: " << pathname << " Pos: " << ref.position().node_id() <<
                  //         " REF: " << ii.sequence() << " Alt: " << alts.front().edit(0).sequence() << ", " << alts.front().position().offset() <<
                  //         ", from_len: " << alts.front().edit(0).from_length() << ", to_len: " << alts.front().edit(0).to_length() << endl;
                  // }

              }
          }
      }


}

int Deconstructor::inDegree(Node n){
    return (*vgraph).edges_on_start[n.id()].size();
    // Index ix;
    // ix.open_read_only( ref_index);
    // vector<Edge> edges_ahead;
    // ix.get_edges_on_start(n.id(), edges_ahead);
    // return edges_ahead.size();
}

bool Deconstructor::on_ref(Node n, int64_t path){
    return (*vgraph).paths.has_node_mapping(n.id());
    // Index ix;
    // ix.open_read_only(ref_index);
    // Mapping m;
    // int64_t pp;
    // bool b;
    // return (ix.get_node_path(n.id(), path, pp, b, m) > 0);
}

bool Deconstructor::on_ref(int64_t node_id, int64_t path){
    return (*vgraph).paths.has_node_mapping(node_id);
    // Index ix;
    // ix.open_read_only(ref_index);
    // Mapping m;
    // int64_t pp;
    // bool b;
    // return (ix.get_node_path(n.id(), path, pp, b, m) > 0);
}

bool Deconstructor::beenVisited(Node n, map<int64_t, int> node_to_level){
  return (node_to_level.find(n.id()) != node_to_level.end());
}


/**
* Uses a BFS to find "Anchor" nodes, nodes on the reference which
* bring branches (variation) back to the reference path.
*
*/
Node Deconstructor::get_alleles(Node n, int64_t path_id, map<int64_t, int>& node_to_level){
    Node current;
    queue<Node> nq;
    nq.push(n);
    int level = 0;
    //cerr << "here " << endl;

    while (!nq.empty()){
        current = nq.front(); nq.pop();
        vector<pair<int64_t, bool>> edges;
        edges = (*vgraph).edges_on_end[(int64_t) current.id()];
        for (int i = 0; i < edges.size(); i++){
            Node n;
            n = *((*vgraph).get_node(edges[i].first));
            if ((inDegree(n) >= 2) && on_ref(n, path_id) && beenVisited(n, node_to_level)){
                //cerr << "Found anchor " << n.id() << endl;
                return n;
            }
            if (!beenVisited(n, node_to_level)){
                node_to_level[n.id()] = level;
                nq.push(n);
            }

            if (level > 4){
              break;
            }
            //cerr << n.id() <<  "level" << level << endl;
        }
        level += 1;
    }
}

Node Deconstructor::get_anchor_node(Node current, int64_t path){
    /** Start at the node which begins the xfurcation
      * Traverse along all branches
      *
      */
      map<int64_t, int> node_to_level;
      Node n = get_alleles(current, path, node_to_level);
      return n;
}

//TODO
map<int64_t, long> Deconstructor::cache_path_positions(string path){
    map <int64_t, long> node_to_position;
    long position = 0;
    queue<Node> nq;
    Path p = (*vgraph).paths.path(path);
    int64_t n_id = (path_start(p)).node_id();
    Node* n_p = (*vgraph).get_node(n_id);
    Node n = *n_p;
    map<int64_t, int> node_to_level;
    int level = 0;

    // Start at first node.
    nq.push(n);
    // BFS forward.
    while (!nq.empty()){
      // For each node, add the length of the nodes previously found.
      n = nq.front(); nq.pop();
      vector<pair<int64_t, bool>> edges;
      edges = (*vgraph).edges_on_end[(int64_t) n.id()];
      long t_len = (long) n.sequence().size();
      node_to_position[n.id()] = position;
      for (int i = 0; i < edges.size(); i++){
        n = *((*vgraph).get_node(edges[i].first));
        if (!beenVisited(n, node_to_level)){
          node_to_level[n.id()] = level;
          nq.push(n);
        }
    }
    position += t_len;
    ++level;
  }
  return node_to_position;

}

//TODO
/**
* Takes a node_to_level map (i.e. a distance from the start node)
* and the anchor node at the end of the bifurcation and returns
* a string representing the type of mutation at a site and
* a map from "ref" to a list of nodes and "alt" to a list of nodes.
*/
pair< string, map<string, list<int64_t>>> Deconstructor::parse_node_level_to_mutation(map<int64_t, int> node_to_level, Node anchor){
  map<string, list<int64_t>> ret;
  ret["ref"] = list<int64_t>();
  ret["alt"] = list<int64_t>();
  bool snp = false;
  bool ins = false;
  bool del = false;
  bool complex = false;

  map<int64_t, int>::iterator it;
  for(it = node_to_level.begin(); it != node_to_level.end(); it++){
    // Ref part of SNP
    if (it->first != anchor.id() &&
     it->second < node_to_level[anchor.id() &&
     on_ref(it->first, 100L)]){
       ret["ref"].push_back(it->first);
       snp = true;
    }
    // Alt portion of snp
    else if (it->first != anchor.id() &&
     it->second < node_to_level[anchor.id()] &&
     !on_ref(it->first, 100L)){
       ret["alt"].push_back(it->first);
       snp = true;
       complex = (ins || del) ? true : false;
    }
    // Deletion, del portion
    else if (it->first != anchor.id() &&
      it->second == node_to_level[anchor.id()] &&
      on_ref(it->first, 100L)){
        ret["ref"].push_back(it->first);
        del = true;
        complex = (ins || snp) ? true : false;
      }
      //Insertion, inserted portion
    else if (it->first != anchor.id() &&
      it->second == node_to_level[anchor.id()] &&
      !on_ref(it->first, 100L)){
        ret["alt"].push_back(it->first);
        ins = true;
        complex = (snp || del) ? true : false;
      }
    else if (it->first == anchor.id()){
      continue;
    }
    else{
      complex = true;
    }
  }

  if (complex){
    return pair< string, map<string, list<int64_t>>>("complex", ret);
  }
  else if (snp){
    return pair< string, map<string, list<int64_t>>> ("snp", ret);
  }
  else if (del){
    return pair< string, map<string, list<int64_t>>> ("del", ret);
  }
  else if (ins){
    return pair< string, map<string, list<int64_t>>>("ins", ret);
  }
  else{
    cerr << "Unknown variation found." << endl;
    return pair< string, map<string, list<int64_t>>>("unknown", ret);
  }
  // SNP case, three or more nodes, same level,
  // neither is anchor (which is up a level), one Ref


  // Del case, two nodes or more, all on the same level as anchor
  // all nodes are reference.

  // Ins case, two nodes or more, all on the same level as anchor,
  // at least one alternate sequence.

}


/**
* Assumes a is before b in the path.
* Takes two nodes and returns a pair with the reference node
* and a Mapping for the
* mutation that occured between them.
*/
map<Node, list<Mapping>> Deconstructor::map_between_nodes(Node a, Node b){
  Node* ret;
  list<Mapping> m_list;
  Mapping mapping; // Pos Edit(s) Rank
  //Position p; // NodeID OffSet IsReverse
  //Edit e; // ToLen FromLen Seq
  bool isBackward = false;
  map<Node, list<Mapping>> ref_to_mutations;
  // *matches* from_length == to_length
  // *snps* from_length == to_length; sequence = alt
  // *deletions* from_length > to_length; sequence may be unset or empty
  // *insertions* from_length < to_length; sequence contains relative insertion
  // *skip* from_length == 0, to_length > 0; implies "soft clip" or sequence skip
  //
  // BFS between the nodes.
  // Keep track of the sequences along the edges.
  // two nodes at the same level, immediately before the anchor: SNP
  // One node, ref: Deletion
  // One node, alt: Insertion
  map<int64_t, int> node_to_level;
  int level = 0;
  queue<Node> nq;
  nq.push(a);

  while (!nq.empty()){
    Node current;
    current = nq.front(); nq.pop();
    vector<pair<int64_t, bool>> edges;
    edges = (*vgraph).edges_on_end[(int64_t) current.id()];
    for (int i = 0; i < edges.size(); i++){
      current = *((*vgraph).get_node(edges[i].first));
      if (!beenVisited(current, node_to_level)){
        node_to_level[current.id()] = level;
        nq.push(current);
      }
      if (current.id() == b.id() && beenVisited(b, node_to_level)){
        map<int64_t, int>::iterator it;
        for (it = node_to_level.begin(); it != node_to_level.end(); it++){
          cerr << it->first << " lev: " << it->second << endl;
        }
        //pair< string, map<string, list<int64_t>>> mut_type_to_nodes;
        //mut_type_to_nodes = parse_node_level_to_mutation(node_to_level, current);
        // //mapping.mutable_position()->set_node_id(node_id);
        // //mapping.mutable_position()->set_is_reverse(isBackward);
        // Edit* e = mapping.add_edit();
        // // SNP
        // if (mut_type_to_nodes.first == "snp"){
        //   ret = (*vgraph).get_node(mut_type_to_nodes.second["ref"].front());
        //   e->set_from_length(0);
        //   e->set_to_length(0);
        //   Node* alt = (*vgraph).get_node(mut_type_to_nodes.second["alt"].front());
        //   e->set_sequence(alt->sequence());
        // }

        // // Insertion
        // else if (mut_type_to_nodes.first == "ins"){
        //   ret = (*vgraph).get_node(mut_type_to_nodes.second["ref"].first());
        //   e->set_from_length(0);
        //   e->set_to_length(1);
        //   e->set_sequence((*vgraph).get_node(mut_type_to_nodes.second["alt"].first()).sequence());
        // }
        //
        // // Deletion
        // else if (mut_type_to_nodes.first == "del"){
        //   ret = (*vgraph).get_node(mut_type_to_nodes.second["ref"].first());
        //   e->set_from_length(1);
        //   e->set_to_length(0);
        //   e->set_sequence((*vgraph).get_node(mut_type_to_nodes.second["alt"].first()).sequence());
        // }
        // else if (mut_type_to_nodes.first == "complex"){
        //   ret = (*vgraph).get_node(mut_type_to_nodes.second["ref"].first());
        //   e->set_from_length(0);
        //   e->set_to_length();
        //   e->set_sequence((*vgraph).get_node(mut_type_to_nodes.second["alt"].first()).sequence());
        // }
        // TODO translocations

        // Unknown???? TODO
        // else{
        //   cerr << "You've called a really complex variant." <<
        //   endl << "We're not sure what to tell you. We'll return " <<
        //   "an empty mapping." << endl;
        // }

        return ref_to_mutations;
      }
    }
    level += 1;
  }

}

void Deconstructor::indel_caller(string pathname){
  enumerate_path_names_in_index();

  //Just argument handling here - either take the single region given
  // or the entire intersection of the reference and index
  vector<string> paths_to_project;
  if (pathname != ""){
    //TODO check if region in reference paths TODO
    paths_to_project.push_back(pathname);
  }
  else{
    map<string, int64_t>::iterator it;
    for (it = inter_ref_and_index.begin(); it != inter_ref_and_index.end(); it++){
      paths_to_project.push_back(it->first);
    }
  }
  for (auto pathname : paths_to_project){
    map<int64_t, long> node_id_to_pos = cache_path_positions(pathname);
    Path path = (*vgraph).paths.path(pathname);
    bool backward = false;
    list<Mapping> m_list = (*vgraph).paths._paths[pathname];
    list<Mapping>::iterator it;
    for(it = m_list.begin(); it != m_list.end(); it++){
        int64_t current_node_id = (*it).position().node_id();
        vector<pair<int64_t, bool>> edges_ahead;
        if(backward) {
            // "next" = right = start
            edges_ahead = (*vgraph).edges_on_start[current_node_id];
        } else {
            // "next" = right = end
            edges_ahead = (*vgraph).edges_on_end[current_node_id];
        }

        // base case: no variation, reference has simply been split up.
        // We ignore self-cycling as a possibility.
        if (edges_ahead.size() < 2){
            continue;
        }
        else{
          //cerr << "nfurc" << endl;
          Node n = (*(*vgraph).get_node(current_node_id));
          //cerr << "Retrieved node: " << n.id() << endl;
          Node anchor = get_anchor_node(n, 100L);
          cerr << pathname << " " << node_id_to_pos[n.id()] << " " <<
           n.sequence() << " " << anchor.sequence() << endl;
          map_between_nodes(n, anchor);
        }
      }

  }
}

void Deconstructor::get_variants_using_edges(string pathname){
    enumerate_path_names_in_index();
    //Just argument handling here - either take the single region given
    // or the entire intersection of the reference and index
    vector<string> paths_to_project;
    Index ix;
    if (pathname != ""){
        //TODO check if region in reference paths TODO
        paths_to_project.push_back(pathname);
    }
    else{
        map<string, int64_t>::iterator it;
        for (it = inter_ref_and_index.begin(); it != inter_ref_and_index.end(); it++){
            paths_to_project.push_back(it->first);
        }
      }
        for (auto pathname : paths_to_project){
            Path path = (*vgraph).paths.path(pathname);
            ix.open_read_only(ref_index);
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

                // base case: no variation, reference has simply been split up.
                // We ignore self-cycling as a possibility.
                if (edges_ahead.size() == 1){
                    continue;
                }
                // Otherwise, we have found an n-furcation of the graph,
                // indicating that variation has been inserted.
                else{
                    // TODO this should really use recursive backtracking or something clever

                    /**
                     * This should catch multi-allelic snps as-is,
                     * but it seems to be a bit off for INDELS.
                     */
                    list<Mapping> alts;
                    Mapping ref;
                    for (int edge_ind = 0; edge_ind < edges_ahead.size(); edge_ind++){
                        // base case: from_node for each path is the same
                        // and to_node for each path is the same.
                        Edge e = edges_ahead[edge_ind];
                        int64_t next_n_id = e.to();
                        Node n;
                        ix.get_node(next_n_id, n);

                        // Case: SNPs and Insertions.
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
                        // Case: reference sequence for SNPs, deletions and insertions.
                        else {
                            vector<Edge> n_e_in;
                            vector<Edge> n_e_out;
                            ix.get_edges_on_start(next_n_id, n_e_in);
                            ix.get_edges_on_end(next_n_id, n_e_out);
                            //if ((n_e_in.size() == 1)){
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
                            //TODO Not happy with indels
                            Node ref_seq_n; ix.get_node((int64_t) m.position().node_id(), ref_seq_n);
                            //cerr << "Ref node: " << m.position().node_id() << " " << ref_seq_n.sequence() << endl;
                            //}
                            // else if (n_e_in.size() > 1){
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
                        //TODO Need to cut sequence using node offset.
                        // Pos: get_node_path_relative_position(int64_t node_id, bool backward, int64_t path_id,
                        // list<pair<int64_t, bool>>& path_prev, int64_t& prev_pos, bool& prev_orientation,
                        // list<pair<int64_t, bool>>& path_next, int64_t& next_pos, bool& next_orientation)
                        list<pair<int64_t, bool>> pp; int64_t ppos; bool p_or;
                        list<pair<int64_t, bool>> nn; int64_t npos; bool n_or;
                        // ref.position().node_id()
                        string pos = "POS NOT IMPLEMENTED";
                        cerr << "CHROM: " << pathname << " Pos: " <<  pos <<
                            " REF: " << ii.sequence() << " Alt: " << alts.front().edit(0).sequence() << ", " << alts.front().position().offset() <<
                            ", from_len: " << alts.front().edit(0).from_length() << ", to_len: " << alts.front().edit(0).to_length() << endl;
                    }

                }
            }
        }

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


        return v;
    }

    void Deconstructor::write_variants(string filename, vector<vcflib::Variant> variants){
        cerr << "writing variants." << endl;
        vcflib::VariantCallFile v;

    }

}

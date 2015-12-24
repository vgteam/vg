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



map<string, vector<vcflib::Variant>> Deconstructor::get_variants_on_paths_from_file(string pathfile){
  map<string, vector<vcflib::Variant>> pathname_to_variants;
  ifstream pfi(pathfile);
  string line;
  while(getline(pfi, line)){
    //TODO probably need to strip the input lines.
    map<string, vector<vcflib::Variant>> t = indel_caller(line);
    map<string, vector<vcflib::Variant>>::iterator it;
    for (it = t.begin(); it != t.end(); it++){
      pathname_to_variants[it->first] = it->second;
    }
  }

  return pathname_to_variants;
}

map<string, vector<vcflib::Variant>> Deconstructor::get_variants_on_path(string pathname){
  map<string, vector<vcflib::Variant>> pathname_to_variants;
  //TODO check if path on reference.
  if (pathname == ""){
    FastaReference fr;
    fr.open(reference);
    vector<string> paths_in_ref;
    for (auto& seq : fr.index->sequenceNames){
        //cout << seq << endl;
        if ((*vgraph).paths.has_path(seq)){
          map<string, vector<vcflib::Variant>> t = indel_caller(seq);
          map<string, vector<vcflib::Variant>>::iterator it;
          for (it = t.begin(); it != t.end(); it++){
            pathname_to_variants[it->first] = it->second;
          }
        }
    }
  }
  else{
    if ((*vgraph).paths.has_path(pathname)){
      pathname_to_variants = indel_caller(pathname);
    }
  }
  return pathname_to_variants;
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

//TODO straight up wrong.
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
    long t_len = 0;
    // BFS forward.
    while (!nq.empty()){
      // For each node, add the length of the nodes previously found.
      n = nq.front(); nq.pop();
      vector<pair<int64_t, bool>> edges;
      edges = (*vgraph).edges_on_end[(int64_t) n.id()];
      //long t_len = 0;
      //TODO t_len is added at the wrong time as the moment
      t_len = (long) n.sequence().size();
      position += t_len;
      for (int i = 0; i < edges.size(); i++){
        n = *((*vgraph).get_node(edges[i].first));
        if (!beenVisited(n, node_to_level)){
          node_to_position[n.id()] = position;
          node_to_level[n.id()] = level;
          nq.push(n);
        }
        //t_len = 0;
    }

    ++level;
  }
  node_to_pos = node_to_position;
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
     it->second < node_to_level[anchor.id()] &&
     on_ref(it->first, 100L)){
       ret["ref"].emplace_back(it->first);
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
map<int64_t, list<Mapping>> Deconstructor::map_between_nodes(Node a, Node b){
  int64_t ret;
  list<Mapping> m_list;
  //Mapping mapping; // Pos Edit(s) Rank
  //Position p; // NodeID OffSet IsReverse
  //Edit e; // ToLen FromLen Seq
  bool isBackward = false;
  map<int64_t, list<Mapping>> ref_to_mutations;
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
        bool snp = false;
        bool del = false;
        bool ins = false;
        bool complex = false;
        for (it = node_to_level.begin(); it != node_to_level.end(); it++){
          //cerr << it->first << " lev: " << it->second << endl;

          if (it->first != b.id() &&
           it->second < node_to_level[b.id()] &&
           on_ref(it->first, 100L)){
             //ret["ref"].emplace_back(it->first);
             //cerr << it->first << " is the node id";

             ret = (int64_t) (*((*vgraph).get_node(it->first))).id();
             //cerr << endl << "And the ret is: " << ret.id() << endl;
             snp = true;

          }
          else if (it->first != b.id() &&
           it->second < node_to_level[b.id()] &&
           !on_ref(it->first, 100L)){
             Node n = *((*vgraph).get_node(it->first));
             Mapping m;
             m.mutable_position()->set_node_id(it->first);
             m.mutable_position()->set_is_reverse(isBackward);
             Edit* e = m.add_edit();
             int32_t to_length = n.sequence().size();
             int32_t from_length = n.sequence().size();
             e->set_to_length(to_length);
             e->set_from_length(from_length);
             e->set_sequence(n.sequence());
             m_list.push_back(m);
             //ret["alt"].push_back(it->first);
             snp = true;
             complex = (ins || del) ? true : false;
             //ref_to_mutations[n] =
          }
          //Deletion
          else if (it->first != b.id() &&
            it->second == node_to_level[b.id()] &&
            on_ref(it->first, 100L)){
              //ret["ref"].push_back(it->first);
              Node n = *((*vgraph).get_node(it->first));
              Mapping m;
              m.mutable_position()->set_node_id(it->first);
              m.mutable_position()->set_is_reverse(isBackward);
              Edit* e = m.add_edit();
              int32_t from_length = n.sequence().size();
              int32_t to_length = 0;
              e->set_from_length(from_length);
              //Node t = *((*vgraph).get_node(ret));
              ret = a.id();
              e->set_to_length(to_length);
              e->set_sequence(n.sequence());
              m_list.push_back(m);
              del = true;
              complex = (ins || snp) ? true : false;
            }
            //Insertion, inserted portion
          else if (it->first != b.id() &&
            it->second == node_to_level[b.id()] &&
            !on_ref(it->first, 100L)){
              //TODO hack: set ref node to left anchor
              // what should this realy be?
              ret = a.id();
              Node n = *((*vgraph).get_node(it->first));
              Mapping m;
              m.mutable_position()->set_node_id(it->first);
              m.mutable_position()->set_is_reverse(isBackward);
              Edit* e = m.add_edit();
              int32_t to_length = n.sequence().size();
              int32_t from_length = 0;
              e->set_to_length(to_length);
              e->set_from_length(from_length);
              e->set_sequence(n.sequence());
              m_list.push_back(m);
              ins = true;
              complex = (snp || del) ? true : false;
            }
          else if (it->first == b.id()){
            continue;
          }

        }


        // TODO translocations

        // Unknown???? TODO
        // else{
        //   cerr << "You've called a really complex variant." <<
        //   endl << "We're not sure what to tell you. We'll return " <<
        //   "an empty mapping." << endl;
        // }
        // if (ret == NULL){
        //   cerr << "Caught error during deconstruction. No reference node found." <<
        //   endl << "Please file a bug report at github/edawson/vg." << endl;
        // }
        //cerr << ret.id() << endl;
        ref_to_mutations[ret] = m_list;

        return ref_to_mutations;
      }
    }
    level += 1;
  }

}

map<string, vector<vcflib::Variant>> Deconstructor::indel_caller(string pathname){
  //enumerate_path_names_in_index();

  //Just argument handling here - either take the single region given
  // or the entire intersection of the reference and index
  vector<string> paths_to_project;
  // if (pathname != ""){
  //   //TODO check if region in reference paths TODO
  //   paths_to_project.push_back(pathname);
  // }
  // else{
  //   map<string, int64_t>::iterator it;
  //   for (it = inter_ref_and_index.begin(); it != inter_ref_and_index.end(); it++){
  //     paths_to_project.push_back(it->first);
  //   }
  // }
  paths_to_project.push_back(pathname);

  map<string, vector<vcflib::Variant>> pathname_to_variants;

  for (auto pathname : paths_to_project){
    vector<vcflib::Variant> variants;
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
          // cerr << pathname << " " << node_id_to_pos[n.id()] << " " <<
          //  n.sequence() << " " << anchor.sequence() << endl;
          map<int64_t, list<Mapping>> node_to_mappings = map_between_nodes(n, anchor);

          map<int64_t, list<Mapping>>::iterator nm;
          //cerr << "Map size" << node_to_mappings.size();

          for (nm = node_to_mappings.begin(); nm != node_to_mappings.end(); nm++){
            mapping_to_simple_variant(pathname, nm->first, nm->second, variants);
          }
          pathname_to_variants[pathname] = variants;
          // for (int k = 0; k < variants.size(); k++){
          //   cerr << variants[k] << endl;
          // }
        }
      }

  }
  return pathname_to_variants;
}



    /**
     * Transforms a Mapping (a list of edits on a path) to a vcf entry.
     * Currently, the mapping contains the reference allele.
     * The reference path that the variant originates from can be grabbed
     * from the Mapping's position.
     */
    void Deconstructor::mapping_to_simple_variant(string pathname,
                            int64_t ref_id, list<Mapping> mappings,
                            vector<vcflib::Variant>& variants){
        //vector<vcflib::Variant> variants;
        Node ref = *((*vgraph).get_node(ref_id));
        map<string, vector<vcflib::VariantAllele> > variantAlleles;
        long pos;
        vector<string> alleles;
        // map<int64_t, long> node_to_pos = Deconstructor::cache_path_positions
        vcflib::Variant v;
        v.alt = vector<string>();
        v.sequenceName = pathname;
        v.position = node_to_pos[ref_id];
        // REF
        vcflib::VariantAllele ref_var(ref.sequence(), ref.sequence(), 1);
        string seq;
        string flank;
        flank = ref.sequence().substr(ref.sequence().length() - 1, ref.sequence().length());
        //v.ref = flank;

        list<Mapping>::iterator it;
        if (mappings.size() == 0){
          cerr << "Mapping size 0" << endl;

        }
        //cerr << "Mapping variant " << endl;
        for (it = mappings.begin(); it != mappings.end(); it++){
          Mapping m = *it;
          // Check the to/from length to determine what kind of variant this is.
          // TODO to handle complex variants, we'll need to use multiple edits.
          // SNPs
          //cerr << m.edit(0).to_length();
          if (m.edit(0).to_length() == m.edit(0).from_length()){
            seq = m.edit(0).sequence();
            v.alt.push_back(seq);
            v.ref = flank;
            node_to_pos[ref.id()];
            vcflib::VariantAllele alt_var(ref.sequence(), seq, node_to_pos[ref.id()]);
            //cerr << alt_var << endl;
          }
          //Ins
          else if (m.edit(0).to_length() > m.edit(0).from_length()){
            seq = m.edit(0).sequence();
            // grab a flanking base from the tail of the ref.
            flank = ref.sequence().substr(ref.sequence().length() - 1, ref.sequence().length());
            v.ref = flank;
            pos = node_to_pos[ref.id()] + ref.sequence().length() - 1;
            v.alt.push_back(flank + seq);
            vcflib::VariantAllele alt_var(flank, flank + seq, pos);
            //cerr << alt_var << endl;
          }
          //Del
          else if (m.edit(0).to_length() < m.edit(0).from_length()){
            seq = m.edit(0).sequence();
            flank = ref.sequence().substr(ref.sequence().length() - 1, ref.sequence().length());
            pos = node_to_pos[ref.id()] + ref.sequence().length() - 1;
            v.ref = flank + seq;
            v.alt.push_back(seq);
            vcflib::VariantAllele alt_var(flank + seq, seq, node_to_pos[ref.id()]);
            //cerr << alt_var << endl;
          }
          //complex TODO unimplemented
          else{
            cerr << "Complex variation found but handling not yet implemented." << endl;
          }

        }
        //v.alt = alts;
        //v.alleles = vector<string>(alts);
        //vector<string>::iterator vit;
        //v.alleles.insert(vit, ref.sequence());
        variants.push_back(v);
        //cout << v << endl;

      //  return variants;
    }

    void Deconstructor::write_variants(string filename, map<string, vector<vcflib::Variant>> pathname_to_variants){
        Header h;
        h.set_date();
        h.set_version("VCF4.1");
        h.set_reference(reference);
        h.set_contig("NA");
        h.set_source("vg_deconstruct");
        //TODO use a var call file instead of just an ostream...
        vcflib::VariantCallFile v;
        cout << h << endl;
        map<string, vector<vcflib::Variant>>::iterator it;
        for (it = pathname_to_variants.begin(); it != pathname_to_variants.end(); it++){
          for (auto v : it->second){
            cout << v << endl;
          }
        }

    }

}

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
     *
     vector<vcflib::Variant> Deconstructor::get_variants(string region_file){
     vector<vcflib::Variant> vars;
     cerr << "Not implemented" << endl;
     exit(1);
     return vars;
     }
     */

    /**
     * If a region name is given, extract variants for that path in the reference.
     * Otherwise, try to get variants for all paths. If a pathname is not in the
     * reference, skip it.
     * First, locate list<Mapping> for nodes that connect to reference paths but
     * that lie off of them.
     * Next, transform this information into a vcf record, append it to a vector,
     * and return.
     *
     *
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
*/

void get_variants_using_edges_from_file(string pathfile){

}

void Deconstructor::get_variants_using_edges(string pathname){
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
        }}
        for (auto pathname : paths_to_project){
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

                        // Case: SNPs and Deletions.
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
                        // Case: reference half of snps, deletions.
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
                        cerr << "CHROM: " << pathname << " Pos: " << ref.position().node_id() <<
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

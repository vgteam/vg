#include "homogenizer.hpp"

using namespace std;
using namespace vg;

void Homogenizer::homogenize(vg::VG* o_graph, xg::XG* xindex, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_index, vg::Index reads_index){
    /**
     * Pattern for SV homogenization
     * 1. Locate SV-indicating reads with Sift. Save them in a gam file
     * 2. index that gam with index -d dbname -N
     * 3. Send those reads here - for each possible position
     *      Find reads supporting that position
     *      Generate candidate edges and nodes
     *      Remap reads locally (w/in some subgraph containing the SV)
     *      Score it somehow
     *      Check the reads again for SV signatures.
     */
    
    }


void Homogenizer::homogenize(vg::VG* o_graph, xg::XG* xindex, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_index, Paths cached_paths, int kmer_size){

    bool in_mem_path_only = true;

    vector<id_t> tips = find_non_ref_tips(o_graph);

    /* TODO filter by whether a read is on the ref path
    // 2. Cache the reference path(s)
    // (Not sure how to grab these, so for now just grab the first path in the graph)
    map<string, list<Mapping> > cached_paths;
    set<string> kept_paths;
    */
    string ref_path = o_graph->paths.all_path_names()[0];

    //o_graph->unchop();
    /*
    Paths p = o_graph->paths;
    //
    // 3. Remove all paths in the graph, except the reference
    (o_graph->paths).keep_paths(kept_paths);
    */

    /* Generate edges/nodes to add to graph */
    //vector<MaximalExactMatch> find_smems(const string& seq);
    Mapper* mapper;
    mapper= new Mapper(xindex, gcsa_index, lcp_index);

    map<vg::id_t, vector<MaximalExactMatch> > matches;
    map<vg::id_t, string> ref_node_to_clip;
    vector<string> seqs;
    for (int i = 0; i < tips.size(); i++){
        Node* n = o_graph->get_node(tips[i]);
        if (n->sequence().length() < 4 || (o_graph->paths).has_node_mapping(n->id())){
            continue;
        }
        vector<MaximalExactMatch> m = mapper->find_mems(n->sequence().begin(),
                                                        n->sequence().end(),
                                                        200, 0);
        // Why >1? Because we need to match the node AND somewhere else in the graph.
        if (m.size() > 1){
            cerr << "POTENTIAL NEW EDGE" << endl;
            matches[tips[i]] = m;
            // map<id_t, map<string, set<Mapping*>>> node_mapping;
            // Get paths of tip
            set<string> paths_of_tip = cached_paths.of_node(tips[i]);
            // Find the closest reference node to the tip
            vg::id_t ref_node = -1;

            for (auto m : paths_of_tip){
                cerr << m << endl;
                Path path_of_tip = cached_paths.path(m);
                if (m == ref_path){
                    continue;     
                }
                else{
                    bool on_ref = false;
                    for (int j = 0; j < path_of_tip.mapping_size(); j++){
                        Mapping nearby_mapping = path_of_tip.mapping(j);
                        Position pos = nearby_mapping.position();
                        vg::id_t n_id = pos.node_id();
                        if (o_graph->paths.has_node_mapping(n_id)){
                            ref_node = n_id;
                            on_ref = true;
                        }
                        else if (on_ref == false && n_id != tips[i]){
                            ref_node = n_id;
                        }
                        else{
                            continue;
                        }

                    }
                }
            }
            if (ref_node != -1){
                ref_node_to_clip[ref_node] = n->sequence();
            }
        }
    }

    cerr << ref_node_to_clip.size() << " candidate edges generated. Modifying graph" << endl;

    //need to remove the tips sequences first.
    //cut_tips(tips, o_graph);

    vector<Path> new_p_vec;
    for (auto x : ref_node_to_clip){
        Alignment clip_aln;
        cerr << "Length of softclip: " << x.second.size() << endl;
        clip_aln = mapper->align(x.second);
        if (clip_aln.score() < 30){
            continue;
        }
        cerr << clip_aln.DebugString();
        Path new_aln_p = clip_aln.path();
        //new_p_vec.clear();
        new_p_vec.push_back(new_aln_p);
        //for (int i = 0; i < new_aln_p.mapping_size(); i++){
        //vector<Translation> tras = o_graph->edit(new_p_vec);
        //translator.load(tras);
        //o_graph->paths.rebuild_mapping_aux();
        //    Edge * e = o_graph->create_edge(x.first, new_aln_p.mapping(i).position().node_id(), false, false);
         //   o_graph->add_edge(*e);
         //   cerr << "Edge made from " << x.first << " to " << new_aln_p.mapping(i).position().node_id() << endl;
        //}
        
        /** Reindex graph and reset mapper **/
        //delete xindex;
        //xindex = new xg::XG(o_graph->graph);
        //delete gcsa_index;
        //delete lcp_index;
        //o_graph->build_gcsa_lcp(gcsa_index, lcp_index, kmer_size, in_mem_path_only, false, 2);
        //delete mapper;
        //mapper = new Mapper(xindex, gcsa_index, lcp_index);


    }

    //vector<vg::id_t> after_tips = find_tips(o_graph);
    //cut_tips(after_tips, o_graph);



    //o_graph->unchop();

    // Remap the paths (reads) that pass through our current set of tips, and
    // see if the overall score of the graph improves.
    //
    //map<id_t, map<string, set<Mapping*>>> node_mapping;
    //

}

void Homogenizer::cut_tips(vector<vg::id_t> tip_ids, vg::VG* graph){
    for (auto i : tip_ids){
        graph->destroy_node(i);
    }
    graph->remove_orphan_edges();
}

void Homogenizer::cut_tips(vg::VG* graph){
    vector<vg::id_t> tips = find_tips(graph);
    cut_tips(tips, graph);
}

void Homogenizer::cut_nonref_tips(vg::VG* graph){

}

vector<vg::id_t> Homogenizer::find_non_ref_tips(vg::VG* graph){
    vector<vg::id_t> ret;
    std::function<void(Node*)> is_tip = [graph, &ret](Node* n){
    if ((graph->start_degree(n) == 0 | graph->end_degree(n) == 0) &&
            !(graph->paths.has_node_mapping(n))){
            #pragma omp critical
            ret.push_back(n->id());
        }
    };
    graph->for_each_node_parallel(is_tip);
    return ret;

}

vector<vg::id_t> Homogenizer::find_tips(vg::VG* graph){
    vector<vg::id_t> ret;
    std::function<void(Node*)> is_tip = [graph, &ret](Node* n){
        if ((graph->start_degree(n) == 0 | graph->end_degree(n) == 0)){
            #pragma omp critical
            ret.push_back(n->id());
        }
    };
    graph->for_each_node_parallel(is_tip);
    return ret;
}

#include "homogenizer.hpp"

using namespace std;
using namespace vg;

void Homogenizer::homogenize(vg::VG* o_graph, xg::XG* xindex, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_index){
    vector<id_t> tips = find_tips(o_graph);

    // TODO filter by whether a read is on the ref path
    // i.e.:
    // 1. Copy the graph
    //    VG* o_graph = new VG(*graph);
    //    cerr << "GRAPH COPIED" << endl;
    // 2. Cache the reference path(s)
    // (Not sure how to grab these, so for now just grab the first path in the graph)
    map<string, list<Mapping> > cached_paths;
    set<string> kept_paths;
    for (auto x : o_graph->paths._paths){
        cached_paths[x.first] = x.second;
        kept_paths.insert(x.first);
        break;
    }
    //
    Paths p = o_graph->paths;
    //
    // 3. Remove all paths in the graph, except the reference
    (o_graph->paths).keep_paths(kept_paths);
    //


    // Remove tiny softclips
    /*    int sc_threshold = 10;
          vector<id_t> remove_me;
          vector<id_t> keeps;
          for (auto t : tips){
          if ((o_graph->get_node(t))->sequence().length() < sc_threshold){
          remove_me.push_back(t);
          }
          else if (!o_graph->paths.has_node_mapping(t)){
          remove_me.push_back(t);
          }
          else{
    //if (o_graph->paths.has_node_mapping(t) && (o_graph->get_node(t))->sequence().length() < sc_threshold){
    keeps.push_back(t);
    }
    //}   
    }
    cut_tips(remove_me, o_graph);
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
        if (n->sequence().length() < 4){
            continue;
        }
        vector<MaximalExactMatch> m = mapper->find_smems(n->sequence());
        if (m.size() >= 1){
            matches[tips[i]] = m;
            // map<id_t, map<string, set<Mapping*>>> node_mapping;
            // Get paths of tip
            map<string, set<Mapping*>> paths_of_tip = (o_graph->paths).get_node_mapping(tips[i]);
            // Find the closest reference node to the tip
            vg::id_t ref_node = -1;
            for (auto m : paths_of_tip){
                for (auto n : m.second){
                    vg::id_t n_id = n->position().node_id();
                    if ((o_graph->paths).has_node_mapping(n_id)){
                        ref_node = n_id;
                    }
                }
            }
            if (ref_node != -1){
                ref_node_to_clip[ref_node] = (o_graph->get_node(tips[i]))->sequence();
            }
        }
    }

    // TODO: need to remove the tips sequences first.
    for (auto x : ref_node_to_clip){
        Alignment clip_aln;
        clip_aln = mapper->align(x.second);
        Edge * e = o_graph->create_edge(x.first, clip_aln.path().mapping(0).position().node_id(), false, false);
        o_graph->add_edge(*e);
    }

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

vector<vg::id_t> Homogenizer::find_tips(vg::VG* graph){
    vector<vg::id_t> ret;
    std::function<void(Node*)> is_tip = [graph, &ret](Node* n){
        if (graph->start_degree(n) == 0 | graph->end_degree(n) == 0){
#pragma omp critical
            ret.push_back(n->id());
        }
    };
    graph->for_each_node_parallel(is_tip);
    return ret;
}

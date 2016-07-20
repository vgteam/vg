#include "homogenizer.hpp"

using namespace std;
using namespace vg;

void Homogenizer::homogenize(vg::VG* graph, xg::XG* xindex, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_index){
    vector<id_t> tips = find_tips(graph);

    // TODO filter by whether a read is on the ref path
    // i.e.:
    // 1. Copy the graph
    // VG* o_graph = new VG(*graph)
    // 2. Cache the reference path(s)
    // (Not sure how to grab these, so for now just grab the first path in the graph)
    // map<string, list<Mapping> > cached_paths;
    // set<string> kept_paths;
    // for (auto x : o_graph.paths()){
    //  cached_paths[x.first] = x.second;
    //  kept_paths.insert(x.first);
    //  break;
    // }
    //
    // Paths p = graph->_paths;
    //
    // 3. Remove all paths in the graph
    // o_graph.keep_paths(kept_paths);
    //


    // Remove tiny softclips
    int sc_threshold = 10;
    vector<id_t> remove_me;
    vector<id_t> keeps;
    for (auto t : tips){
        if ((graph->get_node(t))->sequence().length() < sc_threshold){
            remove_me.push_back(t);
        }
        else if (!graph->paths.has_node_mapping(t)){
            remove_me.push_back(t);
        }
        else{
            keeps.push_back(t);
        }   
    }
    cut_tips(remove_me, graph);

    /* Generate edges/nodes to add to graph */
    //vector<MaximalExactMatch> find_smems(const string& seq);
    Mapper* mapper;
    mapper= new Mapper(xindex, gcsa_index, lcp_index);

    map<vg::id_t, vector<MaximalExactMatch> > matches;
    for (int i = 0; i < keeps.size(); i++){
        Node* n = graph->get_node(keeps[i]);
        vector<MaximalExactMatch> m = mapper->find_smems(n->sequence());
        for (auto y : m){
            cerr << y.sequence() << endl;
        }
        matches[keeps[i]] = m;
    }


    // TODO: need to remove the tips sequences first.
    for (auto x : matches){
        Alignment clip_aln;
        clip_aln = mapper->align((graph->get_node(x.first))->sequence());
        Edge * e = graph->create_edge(x.first, clip_aln.path().mapping(0).position().node_id(), false, false);
        graph->add_edge(*e);
    }

    // Remap the paths (reads) that pass through our current set of tips, and
    // see if the overall score of the graph improves.
    //
    //map<id_t, map<string, set<Mapping*>>> node_mapping;
    //

    //delete o_graph TODO


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

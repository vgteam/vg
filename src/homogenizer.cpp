#include "homogenizer.hpp"

using namespace vg;

void Homogenizer::homogenize(vg::VG* graph, xg::XG* xindex, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_index){
    vector<id_t> tips = find_tips(graph);

    // TODO filter by whether a read is on the ref path

    // Remove tiny softclips
    int sc_threshold = 10;
    vector<id_t> remove_me;
    for (auto t : tips){
        if ((graph->get_node(t))->sequence().length() < sc_threshold){
            remove_me.push_back(t);
        }
    }
    cut_tips(remove_me, graph);

    /* Generate edges/nodes to add to graph */
    //vector<MaximalExactMatch> find_smems(const string& seq);
    //Mapper* mapper;
    //mapper= new Mapper(xindex, gcsa_index, lcp_index);

    //mapper->find_smems(s);



}

void cut_tips(vector<vg::id_t> tip_ids, vg::VG* graph){
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

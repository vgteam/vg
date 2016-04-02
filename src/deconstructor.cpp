#include "deconstructor.hpp"

using namespace std;


namespace vg {
	Deconstructor::Deconstructor(){

	}

	Deconstructor::Deconstructor(VG* v) {
		my_vg = v;;
  	
		init();
	}

	Deconstructor::~Deconstructor(){
	}

	void Deconstructor::init(){
	}

/**
*   BFS through a superbubble and fill out the corresponding SuperBubble
*   struct.
*/
SuperBubble Deconstructor::report_superbubble(int64_t start, int64_t end){
	
	SuperBubble sb;
	return sb;
}

bool Deconstructor::is_nested(SuperBubble sb){
    return false;
}

/**
* Uses a BFS between nodes in the graph
* labeled as the endpoints of superbubbles
* to enumerate the nodes between them.
*TODO: the dagify transform records the node translation
*/
vector<SuperBubble> Deconstructor::get_all_superbubbles(){
    
    pair<id_t, id_t> endpoints;
    
    
    vector<pair<id_t, id_t> > supbubs = my_vg->get_superbubbles();
    SuperBubble x;
    vector<SuperBubble> ret(supbubs.size(), x);

    /**
    * For each superbubble, BFS through it and record
    * possible paths.
    */
    
    for (int i = 0; i < supbubs.size(); i++){
        vector<id_t> alleles;
        queue<id_t> nq;
        endpoints = supbubs[i];
        nq.push(endpoints.first);
        cerr << endpoints.first << "\t" << endpoints.second << endl;
        while(!nq.empty()){
            id_t current = nq.front(); nq.pop();
            vector<pair<id_t, bool>> edges_on_end = my_vg->edges_end(current);
            for (int j = 0; j < edges_on_end.size(); j++){
                nq.push(edges_on_end[j].first);    
            }
             id_t next_id = (edges_on_end[j].first);
             if (next_id == endpoints.second){
                    ret[i].nodes = alleles;
                    ret[i].start_node = endpoints.first;
                    ret[i].end_node = endpoints.second;
               }
            
        }
        
    }


	return ret;
}


vector<int64_t> Deconstructor::nt_to_ids(deque<NodeTraversal>& nt){
	vector<int64_t> ret = vector<int64_t>(nt.size(), 0);
	int64_t count = 0;
	for (auto n: nt){
		ret[(n.node->id() - 1)] = count;
		count++;
	}

	return ret;
}

}

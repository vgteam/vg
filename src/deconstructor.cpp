#include "deconstructor.hpp"

using namespace std;


namespace vg {
    Deconstructor::Deconstructor(){

    }

    Deconstructor::Deconstructor(VG* v) {
        my_vg = v;
        

        init();
    }

    Deconstructor::~Deconstructor(){
    }

    void Deconstructor::init(){
    }

    void Deconstructor::unroll_my_vg(int steps){
       *my_vg = my_vg->unfold(steps, my_unroll_translation);
        if (my_translation.size() > 0){
            my_vg->overlay_node_translations(my_unroll_translation, my_translation);
        }

    }

    void Deconstructor::dagify_my_vg(int steps){
        *my_vg = my_vg->dagify(steps, my_dagify_translation, 0, my_max_length);
        if (my_translation.size() > 0){
            my_vg->overlay_node_translations(my_dagify_translation, my_translation);
        }
    }
    
    /**
     * For each superbubble in the graph:
     *  If a superbubble is nested and simple (contains no superbubbles),
     *  transform it into a node.
     *  Record the translation from new node in the graph -> old superbubble
     *  map<id_t, SuperBubble>
     *
     *  At each step, find the new superbubbles of the graph and continue with this process.
     *
     *
     */
    vg::VG* Deconstructor::compact(int compact_steps){
        
        
    }

    SuperBubble Deconstructor::translate_id(id_t transformed){
        return id_to_bub[transformed];

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

    void Deconstructor::sb2vcf(string outfile){
        Header h;
        h.set_date();
        h.set_source("VG");
        h.set_reference("");
        h.set_version("VCF4.2");

        cout << h << endl;
        for (auto s : my_superbubbles){
          map<int, vector<id_t> >::iterator it;
          for (it = s.level_to_nodes.begin(); it != s.level_to_nodes.end(); ++it){
                vcflib::Variant v;
                v.sequenceName = "x";
                v.ref = "A";
                v.position = 100;
                v.alt = vector<string>();
                v.alt.push_back("T");
                cout << v << endl;
            }
            
        }

        
    }

    /**
     * Uses a BFS between nodes in the graph
     * labeled as the endpoints of superbubbles
     * to enumerate the nodes between them.
     *TODO: the dagify transform records the node translation
     */
    vector<SuperBubble> Deconstructor::get_all_superbubbles(){

        pair<id_t, id_t> endpoints;
        //map<id_t, pair<id_t, bool> > node_translation;
        //uint32_t dag_len = 1;
        //my_vg->dagify(dag_len, node_translation);


        vector<pair<id_t, id_t> > supbubs = my_vg->get_superbubbles();
        SuperBubble x;
        vector<SuperBubble> ret(supbubs.size(), x);
        map<int, vector<id_t> > level_to_node;
        int level = 0;

        /**
         * For each superbubble, BFS through it and record
         * possible paths.
         *
         */
        map<int, vector<id_t> > level_to_nodes;
        for (int i = 0; i < supbubs.size(); i++){
            vector<id_t> alleles;
            queue<id_t> nq;
            endpoints = supbubs[i];
            nq.push(endpoints.first);
            while(!nq.empty()){
                id_t current = nq.front(); nq.pop();
                vector<pair<id_t, bool>> edges_on_end = my_vg->edges_end(current);
                for (int j = 0; j < edges_on_end.size(); j++){
                    nq.push(edges_on_end[j].first);
                    id_t next_id = (edges_on_end[j].first);
                    if (next_id == endpoints.second){
                        ret[i].level_to_nodes = level_to_nodes;
                        ret[i].start_node = endpoints.first;
                        ret[i].end_node = endpoints.second;
                        break;
                    }
                    alleles.push_back(next_id);
                }
                level++;
            }

        }
        
        //Use the DFS interface instead
        
        // vector<
        // function<void(Node*)> is_end = [](Node* node){
            
        // };
        // function<void(Node*)> is_start = [](Node* node){
            
        // };
        // my_vg->dfs(is_start, is_end);


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

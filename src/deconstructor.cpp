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
    
    void Deconstructor::set_xg(xg::XG* xindex){
        my_xg = xindex;
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
     * detect if there are superbubbles contained within the current superbubble
     * (defined by Start and End)
     *
     * This is easiest done using a simple linear search between the nodes
     * in topologically order.
     *
     */
    bool Deconstructor::contains_nested(pair<int64_t, int64_t> start_and_end){
        /*for (auto i: my_sbs){
            if (i.first > start_and_end.first & i.second < start_and_end.second ){
                return true;
            }
        }*/
        return false;
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
        map<id_t, SuperBubble> node_to_sb;
        for (int i = 0; i < compact_steps; i++){
            vector<pair<int64_t, int64_t> > supbubs = my_vg->get_superbubbles();
            if (supbubs.size() == 0){
                break;
            }
            for (pair<int64_t, int64_t> s : supbubs){
                if (!contains_nested(s)){
                    //Generate a SuperBubble object,
                    //which contains the startnode, endnode, and a trace through the interior.
                    SuperBubble bub = report_superbubble(s.first, s.second);

                    // Swap out the superbubble for its start node.

                    // Record the translation.


                    // Link up the right edges


                }

            }

        }

        return my_vg;


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
        sb.start_node = start;
        sb.end_node = end;

        //BFS through interior and fill out level_to_nodes

        queue<id_t> nq;
        nq.push(start);
        int level = 0;
        map<int, vector<id_t> > level_to_nodes;
        while(!nq.empty()){
            id_t current = nq.front(); nq.pop();
            level_to_nodes[level].push_back(current);
            vector<pair<id_t, bool>> edges_on_end = my_vg->edges_end(current);
            for (int j = 0; j < edges_on_end.size(); j++){
                nq.push(edges_on_end[j].first);
                id_t next_id = (edges_on_end[j].first);
                if (next_id == end){
                    sb.level_to_nodes = level_to_nodes;
                    break;
                }
                
            }
            level++;
        }
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

        // for each superbubble:
        // Fill out a vcflib Variant
        // Check if it is masked by an input vcf
        // if not, print it to stdout
        
        
        
        map<id_t, vcflib::Variant> node_to_var;
        vcflib::VariantCallFile mask;
        if (!mask_file.empty()){
            //node_to_var = my_vg->get_node_to_variant(mask);
        }
        for (auto s : my_sbs){
            vector<vcflib::Variant> varis;
            vcflib::Variant var;
            /*set ref = str(ref)
             * First fill out alleles[0] with ref
             * then alts as 1....N
             *
             * push the alternate alleles into alts[] too
             *
             * Position: your guess is as good as mine
             *
             * id: need annotations
             *
             *
             */

            // Make subgraphs out of the superbubble:
            // Operating on a pair<id_t, id_t>, vector<id_t>
            // then enumerate k_paths through the SuperBubbles
            set<Node*> nodes;
            set<Edge*> edges;
            //nodes.insert(my_vg->get_node(s.first.first));
            //vector<Edge*> edges_entr = my_vg->edges_from(my_vg->get_node(s.first.first));
            //edges.insert(edges_entr.begin(), edges_entr.end());

            for (int i = 0; i < s.second.size(); i++){
                id_t n_id = s.second[i];
                //cerr << n_id << endl;
                Node* n_node = my_vg->get_node(n_id);
                vector<Edge*> e_end = my_vg->edges_from(n_node);
                nodes.insert(n_node);
                if (i < s.second.size() - 1){
                    edges.insert(e_end.begin(), e_end.end());
                }
            }
            //nodes.insert(my_vg->get_node(s.first.second));

            vg::VG t_graph = vg::VG(nodes, edges);

            vector<Path> paths;
            map<id_t, Path> node_id_to_path;

            std::function<void(NodeTraversal)> no_op = [](NodeTraversal n){};
            std::function<void(size_t, Path&)> extract_path = [&paths](size_t x_size, Path& path){
                paths.push_back(path);
            };
            
            t_graph.for_each_kpath(10000, false, 100, no_op, no_op, extract_path);

            std::function<std::vector<Path>(vector<Path>)> uniquify = [](vector<Path> v){
                map<string, Path> unqs;
                vector<Path> ret;
                for (auto x: v){
                    unqs[to_string(x)] = x;
                }

                for (auto y : unqs){
                    ret.push_back(y.second);
                }
                return ret;
            };

            paths = uniquify(paths);

            /*
             * This means we now have vectors for the superbubble
             * that have the paths through the nodes within it (including end nodes)
             * however, these paths are repeated several times.
             * We should find a way to prevent them being inserted once for each node.
             *
             * Next on the agenda: use the get_path_dist thing from vg call / vg stats
             * to get the distance to the head node.
             * Might need an XG index for this.
             *
             * Also need a way to deal with GAMs for this i.e. a way to 
             * count the number of times we see something come up in the gam
             */

            for (auto x : paths){
                //cerr << "Path: ";
                for (int m_i = 1; m_i < x.mapping_size() -1 ; m_i++){
                    Mapping m = x.mapping(m_i);
                    id_t pos_id = m.position().node_id();
                    //cerr << pos_id << " ";
                    map<string, set<Mapping*> > path_to_mappings =  my_vg->paths.get_node_mapping(pos_id);
                    if (path_to_mappings.size() > 0){
                        //TODO check positions; if not equal, create a fresh variant.
                        for (auto y : path_to_mappings){
                            var.sequenceName = y.first;
                            //cerr << y.first << "; ";
                        //int position = my_xg.approx_path_distance(var.sequenceName, 1, pos_id);
                        }
                        var.alleles.push_back((my_vg->get_node(pos_id))->sequence());
                        var.ref = (my_vg->get_node(pos_id))->sequence();
                        var.position = my_xg->approx_path_distance(var.sequenceName, 1, pos_id);
                    }
                    else{
                        
                        var.alt.push_back((my_vg->get_node(pos_id))->sequence());
                    }
                    

                }

                //cerr << endl;
            }

            var.id = ".";
            

            var.position = 0;
            // Get the reference path
             

            // Get the list of alts


            // find the position along the reference path

            //var.sequenceName = "seq";
            /*map<int, vector<id_t> >::iterator it;
            for (it = s.level_to_nodes.begin(); it != s.level_to_nodes.end(); it++){
            //cout << s.start_node << "_" << s.end_node << "\t";
                vector<id_t> middle_nodes = it->second;
                for (int ind = 0; ind < middle_nodes.size(); ind++){
                    Node* n = my_vg->get_node(middle_nodes[ind]);
                    if (my_vg->paths.has_node_mapping(n)){
                        var.ref = n->sequence();
                        var.alleles.push_back(n->sequence());
                    }
                    else{
                        var.alt.push_back(n->sequence());
                    }
                    //cout << middle_nodes[ind] << "\t";
                }
            }*/

            cout << var << endl;


        }
        //cerr << endl;


    }

    /**
     * Uses a BFS between nodes in the graph
     * labeled as the endpoints of superbubbles
     * to enumerate the nodes between them.
     *TODO: the dagify transform records the node translation

     * IDEALLY: return the topological order, the starts/ends of superbubbles,
     * and an index from node -> location in topo order. This makes
     * checking if things are nested trivial.
     */
    map<pair<id_t, id_t>, vector<id_t> >  Deconstructor::get_all_superbubbles(){

        my_sbs = my_vg->superbubbles();
        return my_sbs;

        /*
        pair<id_t, id_t> endpoints;
        //map<id_t, pair<id_t, bool> > node_translation;
        //uint32_t dag_len = 1;
        //my_vg->dagify(dag_len, node_translation);

        vector<SuperBubble> ret;
        unordered_map<id_t, SuperBubble> entrance_to_SB;
        unordered_map<id_t, SuperBubble> exit_to_SB;

        vector<pair<id_t, id_t> > supbubs = my_vg->get_superbubbles();
        for (auto pp : supbubs){
            SuperBubble bub;
            bub.start_node = pp.first;
            bub.end_node = pp.second;

            entrance_to_SB[bub.start_node] = bub;
            exit_to_SB[bub.end_node] = bub;

            //ret.push_back(bub);
        }
        map<int, vector<id_t> > level_to_node;
        int level = 0;
        vector<id_t> rto;


        function<bool(id_t)> is_exit_node = [&exit_to_SB](id_t id){
            try{
                exit_to_SB.at(id);
                return true;
            }
            catch (const std::out_of_range& oor){
                return false;
            }
        };

        function<bool(id_t)> is_entrance_node = [&entrance_to_SB](id_t id){
            try{
                entrance_to_SB.at(id);
                return true;
            }
            catch (const std::out_of_range& oor){
                return false;
            }
        };

        function<void(Node*)> on_node_end = [&rto](Node* node){
            //cerr << "is_end: " << node->id() << endl;
            rto.push_back(node->id());
         };
         function<void(Node*)> on_node_start = [&ret](Node* node){
            //cerr << "is_start: " << node->id() << endl;
            
         };
         my_vg->dfs(on_node_start, on_node_end);
        
         /**
          * Go through nodes in topo order
          * If they are entrance/exits, switch to concat mode
          * concat middle nodes to superbubble
          * if node is exit, switch off concat mode
          * Need to maintain a pointer to SBs in map
          * Then, do a linear pass over map, create the set of pointers
          * and return
          *
          *for (int i = 0; i < rto.size(); i++){
            if (is_entrance_node(rto[i])){
                ret.push_back(current);
                in_bub = false;
            }
            else if (in_bub){
                current.level_to_nodes[0].push_back(rto[i]);
            }

            if (is_exit_node(rto[i])){
                in_bub = true;
                current = exit_to_SB[rto[i]];
            }

            
        }
          
        bool in_bub = false;
        SuperBubble current;
        for (int i = rto.size() - 1; i > 0; i--){
            if (is_exit_node(rto[i])){
                ret.push_back(current);
                in_bub = false;
            }
            else if (in_bub){
                current.level_to_nodes[0].push_back(rto[i]);
            }

            if (is_entrance_node(rto[i])){
                in_bub = true;
                current = entrance_to_SB[rto[i]];
            }

            
        }

        reverse_topo_order = rto;
        

        return ret;*/
    }

/*    map<id_t, int> cache_path_distances(){

    } */


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

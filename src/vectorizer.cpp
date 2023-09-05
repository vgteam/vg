#include "vectorizer.hpp"

using namespace std;
using namespace vg;
using namespace sdsl;
Vectorizer::Vectorizer(const PathPositionHandleGraph* x) : my_xg(x){
    size_t rank = 1;
    my_xg->for_each_handle([&](handle_t handle) {
            id_to_rank[my_xg->get_id(handle)] = rank++;
        });
}

Vectorizer::~Vectorizer(){
}

string Vectorizer::output_wabbit_map(){
    unordered_map<string, int>::iterator wab_it;
    stringstream sout;
    for (wab_it = wabbit_map.begin(); wab_it != wabbit_map.end(); wab_it++){
        sout << wab_it->second << "\t" << wab_it->first << "\n";
    }
    return sout.str();
}

void Vectorizer::emit(ostream &out, bool r_format=false, bool annotate=false){
    /**TODO print header*/
    //size_t ent_size = my_xg.node_count + my_xg.edge_count;
    // if (annotate){
    //   out << "Alignments" << "\t";
    //   for (int i = 0; i < my_vectors[0].size(); i++){
    //     if (my_xg.entity_is_node(i)){
    //         out << my_xg.rank_to_id(my_xg.entity_rank_as_node_rank(i));
    //     }
    //     else{
    //       out << "edge";
    //     }
    //       if (i < my_vectors[0].size() - 1){
    //         out << "\t";
    //       }
    //   }
    // }

    if (annotate){
        r_format = true;
        assert(my_names.size() == my_vectors.size());
    }

    int count = 0;
    for (auto v : my_vectors){
        if (annotate){
            out << my_names[count] << "\t";
        }
        if (r_format){
            out << format(v) << endl;
        }
        else{
            out << v << endl;
        }

        count += 1;
    }
}

void Vectorizer::add_bv(bit_vector v){
    my_vectors.push_back(v);
}

void Vectorizer::add_name(string n){
    my_names.push_back(n);
}

vector<int> Vectorizer::alignment_to_a_hot(Alignment a){
    int64_t entity_size = my_xg->get_node_count();
    vector<int> ret(entity_size, 0);
    Path path = a.path();
    for (int i = 0; i < path.mapping_size(); i++){
        Mapping mapping = path.mapping(i);
        if(! mapping.has_position()){
            continue;
        }
        Position pos = mapping.position();
        int64_t node_id = pos.node_id();
        if (!node_id) continue;
        int64_t key = id_to_rank[node_id];
        vector<step_handle_t> node_steps;
        my_xg->for_each_step_on_handle(my_xg->get_handle(node_id), [&](step_handle_t step) {
                node_steps.push_back(step);
            });
        if (node_steps.size() > 0){
            ret[key - 1] = 2;
        }
        else{
            ret[key - 1] = 1;
        }
    }
    return ret;
}

vector<double> Vectorizer::alignment_to_identity_hot(Alignment a){
    int64_t entity_size = my_xg->get_node_count();
    vector<double> ret(entity_size, 0.0);
    Path path = a.path();
    for (int i = 0; i < path.mapping_size(); i ++){
        Mapping mapping = path.mapping(i);
        if(! mapping.has_position()){
            continue;
        }
        Position pos = mapping.position();
        int64_t node_id = pos.node_id();
        if (!node_id) continue;
        int64_t key = id_to_rank[node_id];

        //Calculate % identity by walking the edits and counting matches.
        double pct_id = 0.0;
        double match_len = 0.0;
        double total_len = 0.0;

        for (int j = 0; j < mapping.edit_size(); j++){
            Edit e = mapping.edit(j);
            total_len += e.from_length();
            if (e.from_length() == e.to_length() && e.sequence() == ""){
                match_len += (double) e.to_length();
            }
            else if (e.from_length() == e.to_length() && e.sequence() != ""){
                // TODO if we map but don't match exactly, add half the average length to match_length
                //match_len += (double) (0.5 * ((double) e.to_length()));
            }
            else{
                
            }
            
        }
        pct_id = (match_len == 0.0 && total_len == 0.0) ? 0.0 : (match_len / total_len);
        ret[key - 1] = pct_id;
    }
    return ret;
}

bit_vector Vectorizer::alignment_to_onehot(Alignment a){
    int64_t entity_size = my_xg->get_node_count();
    bit_vector ret(entity_size, 0);
    Path path = a.path();
    for (int i = 0; i < path.mapping_size(); i++){
        Mapping mapping = path.mapping(i);
        if(! mapping.has_position()){
            continue;
        }
        Position pos = mapping.position();
        int64_t node_id = pos.node_id();
        if (!node_id) continue;
        int64_t key = id_to_rank[node_id];
        //Find entity rank of edge
        ret[key - 1] = 1;
    }
    return ret;
}

vector<double> Vectorizer::alignment_to_custom_score(Alignment a, std::function<double(Alignment)> lambda ){
    vector<double> ret;
    

    return ret;
}

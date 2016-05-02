#include "vectorizer.hpp"

using namespace std;
using namespace vg;
using namespace sdsl;
using namespace stream;
Vectorizer::Vectorizer(xg::XG* x) : my_xg(x){

}

Vectorizer::~Vectorizer(){
    delete my_xg;
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
    int64_t entity_size = my_xg->node_count + my_xg->edge_count;
    vector<int> ret(entity_size, 0);
    Path path = a.path();
    for (int i = 0; i < path.mapping_size(); i++){
        Mapping mapping = path.mapping(i);
        Position pos = mapping.position();
        int64_t node_id = pos.node_id();
        int64_t key = my_xg->node_rank_as_entity(node_id);
        // Okay, solved the previous out of range errors:
        // We have to use an entity-space that is |nodes + edges + 1|
        // as nodes are indexed from 1, not from 0.
        //TODO: this means we may one day have to do the same bump up
        // by one for edges, as I assume they are also indexed starting at 1.
        //cerr << key << " - " << entity_size << endl;

        //Find edge by current / previous node ID
        // we can check the orientation, though it shouldn't **really** matter
        // whether we catch them in the forward or reverse direction.
        if (i > 0){
            Mapping prev_mapping = path.mapping(i - 1);
            Position prev_pos = prev_mapping.position();
            int64_t prev_node_id = prev_pos.node_id();
            if (my_xg->has_edge(prev_node_id, false, node_id, false)){
                int64_t edge_key = my_xg->edge_rank_as_entity(prev_node_id, false, node_id, false);
                vector<size_t> edge_paths = my_xg->paths_of_entity(edge_key);
                if (edge_paths.size() > 0){
                    ret[edge_key - 1] = 1;
                }
                else{
                    ret[edge_key - 1] = 2;
                }
            }
        }
        //Check if the node of interest is on a path
        vector<size_t> node_paths = my_xg->paths_of_node(node_id);
        if (node_paths.size() > 0){
            ret[key - 1] = 2;
        }
        else{
            ret[key - 1] = 1;
        }

    }

    return ret;

}

vector<double> Vectorizer::alignment_to_identity_hot(Alignment a){
    int64_t entity_size = my_xg->node_count + my_xg->edge_count;
    vector<double> ret(entity_size, 0.0);

    Path path = a.path();
    for (int i = 0; i < path.mapping_size(); i ++){
        Mapping mapping = path.mapping(i);
        Position pos = mapping.position();
        int64_t node_id = pos.node_id();
        int64_t key = my_xg->node_rank_as_entity(node_id);

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
        pct_id = match_len / total_len;
        ret[key - 1] = pct_id;

        if (i > 0){
            Mapping prev_mapping = path.mapping(i - 1);
            Position prev_pos = prev_mapping.position();
            int64_t prev_node_id = prev_pos.node_id();
            if (my_xg->has_edge(prev_node_id, false, node_id, false)){
                int64_t edge_key = my_xg->edge_rank_as_entity(prev_node_id, false, node_id, false);
                ret[edge_key - 1] = 1.0;
            }
        }
    }
    return ret;
}

bit_vector Vectorizer::alignment_to_onehot(Alignment a){
    // Make a vector as large as the | |nodes| + |edges| | space
    // TODO handle edges
    int64_t entity_size = my_xg->node_count + my_xg->edge_count;
    bit_vector ret(entity_size, 0);
    Path path = a.path();
    for (int i = 0; i < path.mapping_size(); i++){
        Mapping mapping = path.mapping(i);
        Position pos = mapping.position();
        int64_t node_id = pos.node_id();
        int64_t key = my_xg->node_rank_as_entity(node_id);
        // Okay, solved the previous out of range errors:
        // We have to use an entity-space that is |nodes + edges + 1|
        // as nodes are indexed from 1, not from 0.
        //TODO: this means we may one day have to do the same bump up
        // by one for edges, as I assume they are also indexed starting at 1.
        //cerr << key << " - " << entity_size << endl;

        //Find edge by current / previous node ID
        // we can check the orientation, though it shouldn't **really** matter
        // whether we catch them in the forward or reverse direction.
        if (i > 0){
            Mapping prev_mapping = path.mapping(i - 1);
            Position prev_pos = prev_mapping.position();
            int64_t prev_node_id = prev_pos.node_id();
            if (my_xg->has_edge(prev_node_id, false, node_id, false)){
                int64_t edge_key = my_xg->edge_rank_as_entity(prev_node_id, false, node_id, false);
                ret[edge_key - 1] = 1;
            }
        }
        //Find entity rank of edge

        ret[key - 1] = 1;
    }

    return ret;
}

vector<double> Vectorizer::alignment_to_custom_score(Alignment a, std::function<double(Alignment)> lambda ){
    vector<double> ret;
    

    return ret;
}

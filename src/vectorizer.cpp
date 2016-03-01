#include "vectorizer.hpp"

using namespace std;
using namespace vg;
using namespace sdsl;
using namespace stream;
Vectorizer::Vectorizer(xg::XG x){
  my_xg = x;
}

Vectorizer::~Vectorizer(){

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

string Vectorizer::format(bit_vector v){
  stringstream sout;
  for (int i = 0; i < v.size(); i++){
    sout << v[i];
    if (i < v.size() - 1){
      sout << "\t";
    }
  }
  return sout.str();
}


bit_vector Vectorizer::alignment_to_onehot(Alignment a){
  // Make a vector as large as the | |nodes| + |edges| | space
  // TODO handle edges
  size_t entity_size = my_xg.node_count + my_xg.edge_count;
  bit_vector ret(entity_size, 0);
  Path path = a.path();
  for (int i = 0; i < path.mapping_size(); i++){
    Mapping mapping = path.mapping(i);
    Position pos = mapping.position();
    int64_t node_id = pos.node_id();
    size_t key = my_xg.node_rank_as_entity(node_id);
    ret[key] = 1;
  }
  //aln_to_onehot[a] = ret;
  return ret;
}

vector<double> Vectorizer::alignment_to_custom_score(Alignment a, std::function<double(Alignment)> lambda ){
  vector<double> ret;


  return ret;
}

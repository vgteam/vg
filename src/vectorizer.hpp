#ifndef VECVEC
#define VECVEC
#include <iostream>
#include <sstream>
#include "sdsl/bit_vectors.hpp"
#include <vector>
#include <unordered_map>
#include "vg.hpp"
#include "xg.hpp"
#include <vg/vg.pb.h>

/**
* This class provides a way to transform
* an alignment (or set of alignments) into
* various vector formats useful in downstream
* analysis. For example, the most basic is a
* one-hot vector displaying coverage at each node/edge
* entity in the graph.
*/
using namespace std;
using namespace sdsl;
using namespace vg;
class Vectorizer{
  public:
    Vectorizer(const PathPositionHandleGraph* x);
    ~Vectorizer();
    void add_bv(bit_vector v);
    void add_name(string n);
    void emit(ostream& out, bool r_format, bool annotate);
    bit_vector alignment_to_onehot(Alignment a);
    vector<int> alignment_to_a_hot(Alignment a);
    vector<double> alignment_to_custom_score(Alignment a, std::function<double(Alignment)> lambda);
    vector<double> alignment_to_identity_hot(Alignment a);
    string output_wabbit_map();
    template<typename T> string format(T v){
        stringstream sout;
        for (int i = 0; i < v.size(); i++){
            sout << v[i];
            if (i < v.size() - 1){
                sout << "\t";
            }   
        }
        return sout.str();
    }
    template<typename T> string wabbitize(string name, T v){
        stringstream sout;
        if (!(wabbit_map.count(name) > 0)){
            wabbit_map[name] = wabbit_map.size();
        }
        sout << wabbit_map[name] << " " << "1.0" << " " << "'" << name
            << " " << "|" << " " << "vectorspace" << " ";
        for (int i = 0; i < v.size(); i++){
            sout << i << ":" << v[i];
            if (i < v.size() - 1){
                sout << " ";
            }
        }
        return sout.str();
    }
  private:
    const PathPositionHandleGraph* my_xg;
    unordered_map<vg::id_t, size_t> id_to_rank;
    //We use vectors for both names and bit vectors because we want to allow the use of duplicate
    // names. This allows things like generating simulated data with true cluster as the name.
    vector<bit_vector> my_vectors;
    vector<string> my_names;
    bool output_tabbed = false;
    bool output_names = false;
    //bool output_wabbit = false;
    unordered_map<string, int> wabbit_map;

};

#endif

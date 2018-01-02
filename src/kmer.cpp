#include "kmer.hpp"

namespace vg {

using namespace std;

void for_each_kmer(const HandleGraph& graph, size_t k, const function<void(const string&, const vector<handle_t>&, const size_t)>& lambda) {
    // for each position on the forward and reverse of the graph
    // todo -- no parallel interface in handlegraph
//#pragma omp parallel for schedule(dynamic, 1)
    graph.for_each_handle([&](const handle_t& handle) {
            // we walk out the kmer length k in a depth first search over the graph
            
            // building a tree
            
            // for each position in the node
            // passing each kmer we encounter and its path (including offset in first node) through the graph to the callback
            // 
        });
}

}

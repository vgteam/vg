#include "mapper.h"

namespace vg {

Mapper::~Mapper(void) {
    // noop
}

Alignment Mapper::align(string& sequence, int kmer_size) {
    Alignment alignment;
    alignment.set_sequence(sequence);
    return align(alignment, kmer_size);
}

Alignment& Mapper::align(Alignment& alignment, int kmer_size) {

    if (index == NULL) {
        cerr << "error:[vg::Mapper] no index loaded, cannot map alignment!" << endl;
        exit(1);
    }

    // establish kmers

    const string& sequence = alignment.sequence();
    vector<string> kmers;
    if (!sequence.empty()) {
        for (int i = 0; i < sequence.size()-kmer_size; ++i) {
            kmers.push_back(sequence.substr(i,kmer_size));
        }
    }

    VG graph; // to return

    for (vector<string>::iterator k = kmers.begin(); k != kmers.end(); ++k) {
        VG g;
        index->get_kmer_subgraph(*k, g);
        graph.extend(g);
    }

    int max_iter = 10;
    int iter = 0;
    int context_step = 1;
    int max_subgraph_size = 0;

    while (max_subgraph_size < sequence.size() && iter < max_iter) {

        index->expand_context(graph, context_step);
        list<VG> subgraphs;
        graph.disjoint_subgraphs(subgraphs);

        map<int, vector<VG*> > subgraphs_by_size;
        // these are topologically-sorted
        for (list<VG>::iterator s = subgraphs.begin(); s != subgraphs.end(); ++s) {
            VG& subgraph = *s;
            int64_t length = subgraph.total_length_of_nodes();
            subgraphs_by_size[length].push_back(&*s);
        }
        max_subgraph_size = subgraphs_by_size.begin()->first;
        ++iter;
    }

    graph.join_heads();

    return graph.align(alignment);

}

}

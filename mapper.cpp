#include "mapper.hpp"

namespace vg {

Mapper::Mapper(Index* idex)
    : index(idex)
    , best_n_graphs(0) {
    kmer_sizes = index->stored_kmer_sizes();
    if (kmer_sizes.empty()) {
        cerr << "error:[vg::Mapper] the index (" 
             << index->name << ") does not include kmers" << endl;
        exit(1);
    }
}

Mapper::~Mapper(void) {
    // noop
}

Alignment Mapper::align(string& sequence, int stride) {
    Alignment alignment;
    alignment.set_sequence(sequence);
    return align(alignment, stride);
}

Alignment& Mapper::align(Alignment& alignment, int stride) {

    if (index == NULL) {
        cerr << "error:[vg::Mapper] no index loaded, cannot map alignment!" << endl;
        exit(1);
    }

    // establish kmers

    const string& sequence = alignment.sequence();
    auto kmers = kmers_of(sequence, stride);

    VG* graph = new VG;

    for (auto& k : kmers) {
        VG g;
        index->get_kmer_subgraph(k, g);
        graph->extend(g);
    }

    int max_iter = 10;
    int iter = 0;
    int context_step = 1;
    int max_subgraph_size = 0;

    do {

        list<VG> subgraphs;
        graph->disjoint_subgraphs(subgraphs);

        // turn me into a lambda
        map<int, vector<VG*> > subgraphs_by_size;
        // these are topologically-sorted
        for (list<VG>::iterator s = subgraphs.begin(); s != subgraphs.end(); ++s) {
            VG& subgraph = *s;
            int64_t length = subgraph.total_length_of_nodes();
            subgraphs_by_size[length].push_back(&*s);
        }
        max_subgraph_size = subgraphs_by_size.begin()->first;

        // pick only the best to work with
        delete graph; graph = new VG;
        auto it = subgraphs_by_size.begin();
        for (int i = 0; (best_n_graphs == 0 || i < best_n_graphs) && it != subgraphs_by_size.end(); ++i, ++it) {
            for (auto g : it->second) {
                graph->extend(*g);
            }
        }

        index->expand_context(*graph, context_step);

        ++iter;
    } while (max_subgraph_size < sequence.size() && iter < max_iter);

    graph->join_heads();
    graph->sort();

    /*
    ofstream f("vg_align.vg");
    graph->serialize_to_ostream(f);
    f.close();
    */

    return graph->align(alignment);

}

set<string> Mapper::kmers_of(const string& seq, const int stride) {
    set<string> kmers;
    if (!seq.empty()) {
        for (int kmer_size : kmer_sizes) {
            for (int i = 0; i < seq.size()-kmer_size; i+=stride) {
                kmers.insert(seq.substr(i,kmer_size));
            }
        }
    }
    return kmers;
}

}

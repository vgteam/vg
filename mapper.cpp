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
    map<int64_t, int> kmer_count;

    for (auto& k : kmers) {
        VG g;
        index->get_kmer_subgraph(k, g);
        g.for_each_node([&kmer_count](Node* node) {
                kmer_count[node->id()]++;
            });
        index->get_connected_nodes(g);
        graph->extend(g);
    }

    /*
    ofstream f("vg_align.vg");
    graph->serialize_to_ostream(f);
    f.close();
    */

    int max_iter = 10;
    int iter = 0;
    int context_step = 1;
    int max_subgraph_size = 0;

    auto update_graph = [this, &max_subgraph_size, &graph, &kmer_count]() {

        list<VG> subgraphs;
        graph->disjoint_subgraphs(subgraphs);

        // turn me into a lambda
        map<int, set<VG*> > subgraphs_by_size;
        map<VG*, int> subgraph_kmer_count;
        // TODO
        map<double, set<VG*> > subgraphs_by_kmer_density;
        // these are topologically-sorted
        for (list<VG>::iterator s = subgraphs.begin(); s != subgraphs.end(); ++s) {
            VG& subgraph = *s;
            int64_t length = subgraph.total_length_of_nodes();
            subgraph.for_each_node([&s, &subgraph_kmer_count, &kmer_count](Node* node) {
                subgraph_kmer_count[&*s] += kmer_count[node->id()];
            });
            subgraphs_by_size[length].insert(&*s);
            subgraphs_by_kmer_density[(double)length/(double)subgraph_kmer_count[&*s]].insert(&*s);
        }
        max_subgraph_size = subgraphs_by_size.begin()->first;

        // pick only the best to work with

        set<VG*> passing_density;
        auto it = subgraphs_by_kmer_density.begin();
        for (int i = 0; (best_n_graphs == 0 || i < best_n_graphs)
                 && it != subgraphs_by_kmer_density.end(); ++i, ++it) {
            for (auto g : it->second) {
                passing_density.insert(g);
            }
        }
        delete graph; graph = new VG;
        for (auto g : passing_density) {
            graph->extend(*g);
        }

    };

    update_graph();

    while (max_subgraph_size < sequence.size() && iter < max_iter) {
        index->expand_context(*graph, context_step);
        update_graph();
        ++iter;
    }

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

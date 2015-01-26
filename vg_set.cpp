#include "vg_set.hpp"

namespace vg {
// sets of VGs on disk

void VGset::transform(std::function<void(VG*)> lambda) {
    for (auto& name : filenames) {
        // load
        VG* g = NULL;
        if (name == "-") {
            g = new VG(std::cin, show_progress);
        } else {
            ifstream in(name.c_str());
            g = new VG(in, show_progress);
            in.close();
        }
        g->name = name;
        // apply
        lambda(g);
        // write to the same file
        ofstream out(name.c_str());
        g->serialize_to_ostream(out);
        out.close();
        delete g;
    }
}

void VGset::for_each(std::function<void(VG*)> lambda) {
    for (auto& name : filenames) {
        // load
        VG* g = NULL;
        if (name == "-") {
            g = new VG(std::cin, show_progress);
        } else {
            ifstream in(name.c_str());
            g = new VG(in, show_progress);
            in.close();
        }
        g->name = name;
        // apply
        lambda(g);
        delete g;
    }
}

int64_t VGset::merge_id_space(void) {
    int64_t max_id = 0;
    auto lambda = [&max_id](VG* g) {
        if (max_id > 0) g->increment_node_ids(max_id);
        max_id = g->max_node_id();
    };
    transform(lambda);
    return max_id;
}

void VGset::store_in_index(Index& index) {
    for_each([&index, this](VG* g) {
        g->show_progress = show_progress;
        index.load_graph(*g);
    });
}

// stores kmers of size kmer_size with stride over paths in graphs in the index
void VGset::index_kmers(Index& index, int kmer_size, int edge_max, int stride) {

    for_each([&index, kmer_size, edge_max, stride, this](VG* g) {

        auto keep_kmer = [&index, this](string& kmer, Node* n, int p) {
            if (allATGC(kmer)) {
                index.put_kmer(kmer, n->id(), p);
            }
        };

        g->create_progress("indexing kmers of " + g->name, g->size());
        g->for_each_kmer_parallel(kmer_size, edge_max, keep_kmer, stride);
        g->destroy_progress();
        index.remember_kmer_size(kmer_size);

    });

}

void VGset::for_each_kmer_parallel(function<void(string&, Node*, int)>& lambda,
                                   int kmer_size, int edge_max, int stride) {
    for_each([&lambda, kmer_size, edge_max, stride, this](VG* g) {
        g->show_progress = show_progress;
        g->progress_message = "processing kmers of " + g->name;
        g->for_each_kmer_parallel(kmer_size, edge_max, lambda, stride);
    });
}


}

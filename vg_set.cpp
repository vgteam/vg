#include "vg_set.hpp"

namespace vg {
// sets of VGs on disk

void VGset::transform(std::function<void(VG*)> lambda) {
    for (auto& name : filenames) {
        // load
        ifstream in(name.c_str());
        VG* g = new VG(in);
        in.close();
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
        ifstream in(name.c_str());
        VG* g = new VG(in);
        in.close();
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
    for_each([&index](VG* g) { index.load_graph(*g); });
}

// stores kmers of size kmer_size with stride over paths in graphs in the index
void VGset::index_kmers(vector<string>& filenames, Index& index, int kmer_size, int stride) {
    auto lambda = [&index, kmer_size](VG* g) {
        string_hash_map<string, hash_map<Node*, int> > kmer_map;
        g->kmers_of(kmer_map, kmer_size);
        index.store_kmers(kmer_map);
    };
    for_each(lambda);
}


}

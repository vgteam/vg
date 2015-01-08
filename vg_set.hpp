#ifndef VG_SET_H
#define VG_SET_H

#include <set>
#include "vg.hpp"
#include "index.hpp"
#include "hash_map.hpp"

namespace vg {

// for dealing with collections of VGs on disk
class VGset {
public:

    set<string> filenames;

    VGset() { };
    VGset(set<string>& files) : filenames(files) { }

    void transform(std::function<void(VG*)> lambda);
    void for_each(std::function<void(VG*)> lambda);

    // merges the id space of a set of graphs on-disk
    // necessary when storing many graphs in the same index
    int64_t merge_id_space(void);

    // stores the nodes in the VGs identified by the filenames into the index
    void store_in_index(Index& index);

    // stores kmers of size kmer_size with stride over paths in graphs in the index
    void index_kmers(vector<string>& filenames, Index& index, int kmer_size, int stride);
    
};

}


#endif

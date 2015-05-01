#ifndef VG_SET_H
#define VG_SET_H

#include <set>
#include <stdlib.h>
#include "vg.hpp"
#include "index.hpp"
#include "hash_map.hpp"

namespace vg {

// for dealing with collections of VGs on disk
class VGset {
public:

    vector<string> filenames;

    VGset()
        : show_progress(false)
        { };

    VGset(vector<string>& files)
        : filenames(files)
        , show_progress(false)
        { };

    void transform(std::function<void(VG*)> lambda);
    void for_each(std::function<void(VG*)> lambda);

    // merges the id space of a set of graphs on-disk
    // necessary when storing many graphs in the same index
    int64_t merge_id_space(void);

    // stores the nodes in the VGs identified by the filenames into the index
    void store_in_index(Index& index);
    void store_paths_in_index(Index& index);

    // stores kmers of size kmer_size with stride over paths in graphs in the index
    void index_kmers(Index& index, int kmer_size, int edge_max, int stride = 1);
    void for_each_kmer_parallel(function<void(string&, Node*, int, list<Node*>&, VG&)>& lambda,
                                int kmer_size, int edge_max, int stride, bool allow_dups);
    void write_gcsa_out(ostream& out, int kmer_size, int edge_max, int stride, bool allow_dups = true);

    bool show_progress;
    
};

}


#endif

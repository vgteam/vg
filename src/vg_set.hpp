#ifndef VG_SET_HPP_INCLUDED
#define VG_SET_HPP_INCLUDED

#include <set>
#include <regex>
#include <stdlib.h>
#include <gcsa/gcsa.h>
#include "vg.hpp"
#include "index.hpp"
#include "xg.hpp"
#include "kmer.hpp"


namespace vg {

// for dealing with collections of VGs on disk
class VGset {
public:

    vector<string> filenames;

    VGset() { };

    VGset(vector<string>& files)
        : filenames(files)
        { };

    void transform(std::function<void(VG*)> lambda);
    void for_each(std::function<void(VG*)> lambda);
    void for_each_graph_chunk(std::function<void(Graph&)> lamda);

    /// Stream through the files and determine the max node id
    id_t max_node_id(void);
    
    /// merges the id space of a set of graphs on-disk
    /// necessary when storing many graphs in the same index
    int64_t merge_id_space(void);

    /// Transforms to a succinct, queryable representation
    void to_xg(xg::XG& index, bool store_threads = false);
    /// As above, except paths with names matching the given regex are removed.
    /// They are returned separately by inserting them into the provided map if not null.
    void to_xg(xg::XG& index, bool store_threads, const regex& paths_to_take, map<string, Path>* removed_paths = nullptr);

    // stores the nodes in the VGs identified by the filenames into the index
    void store_in_index(Index& index);
    void store_paths_in_index(Index& index);

    // stores kmers of size kmer_size with stride over paths in graphs in the index
    void index_kmers(Index& index, int kmer_size, bool path_only, int edge_max, int stride = 1, 
                     bool allow_negatives = false);
    void for_each_kmer_parallel(int kmer_size, const function<void(const kmer_t&)>& lambda);
    
    /**
     * Write out kmer lines to GCSA2.
     * size_limit is the maximum space usage for the kmer files in bytes. When the
     * function returns, size_limit is the total size of the kmer files in bytes.
     */
    void write_gcsa_kmers_ascii(ostream& out, int kmer_size,
                                int64_t head_id=0, int64_t tail_id=0);
    void write_gcsa_kmers_binary(ostream& out, int kmer_size, size_t& size_limit,
                                 int64_t head_id=0, int64_t tail_id=0);
    vector<string> write_gcsa_kmers_binary(int kmer_size, size_t& size_limit,
                                           int64_t head_id=0, int64_t tail_id=0);

    // Should we show our progress running through each graph?             
    bool show_progress = false;

};

}


#endif

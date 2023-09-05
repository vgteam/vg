#ifndef VG_SET_HPP_INCLUDED
#define VG_SET_HPP_INCLUDED

#include <set>
#include <regex>
#include <functional>
#include <stdlib.h>
#include <gcsa/gcsa.h>
#include "handle.hpp"
#include "xg.hpp"
#include "vg.hpp"
#include "kmer.hpp"


namespace vg {

// for dealing with collections of HandleGraphs on disk
class VGset {
public:

    vector<string> filenames;

    VGset() { };

    VGset(vector<string>& files)
        : filenames(files)
        { };

    void transform(std::function<void(MutableHandleGraph*)> lambda);
    void for_each(std::function<void(HandleGraph*)> lambda);

    /// Stream through the files and determine the max node id
    id_t max_node_id(void);
    
    /// merges the id space of a set of graphs on-disk
    /// necessary when storing many graphs in the same index
    int64_t merge_id_space(void);

    /// Transforms to a succinct, queryable representation
    void to_xg(xg::XG& index);

    /// As above, except paths with names matching the given predicate are removed.
    /// They are returned separately by inserting them into the provided map if not null.
    void to_xg(xg::XG& index, const function<bool(const string&)>& paths_to_take,
               map<string, Path>* removed_paths = nullptr);

    /// Iterate over all kmers in the graph.
    void for_each_kmer_parallel(size_t kmer_size, const function<void(const kmer_t&)>& lambda);

    /**
     * Write out kmer lines to GCSA2.
     * size_limit is the maximum space usage for the kmer files in bytes. When the
     * function returns, size_limit is the total size of the kmer files in bytes. If the size
     * limit is exceeded, throws SizeLimitExceededException, and (if relevant) deletes
     * temporary files
     */
    void write_gcsa_kmers_ascii(ostream& out, int kmer_size,
                                nid_t head_id=0, nid_t tail_id=0);
    void write_gcsa_kmers_binary(ostream& out, int kmer_size, size_t& size_limit,
                                 nid_t head_id=0, nid_t tail_id=0);
    vector<string> write_gcsa_kmers_binary(int kmer_size, size_t& size_limit,
                                           nid_t head_id=0, nid_t tail_id=0);

    // Should we show our progress running through each graph?             
    bool show_progress = false;

    // Use progress bars if show_progress is true?
    bool progress_bars = true;
};

}


#endif

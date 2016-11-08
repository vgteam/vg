#ifndef VG_SET_H
#define VG_SET_H

#include <set>
#include <regex>
#include <stdlib.h>
#include "gcsa.h"
#include "vg.hpp"
#include "index.hpp"
#include "xg.hpp"


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

    // merges the id space of a set of graphs on-disk
    // necessary when storing many graphs in the same index
    int64_t merge_id_space(void);

    // Transforms to a succinct, queryable representation
    xg::XG to_xg(bool store_threads = false);
    // As above, except paths with names matching the given regex are removed
    // and returned separately by inserting them into the provided map.
    xg::XG to_xg(bool store_threads, const regex& paths_to_take, map<string, Path>& removed_paths);

    // stores the nodes in the VGs identified by the filenames into the index
    void store_in_index(Index& index);
    void store_paths_in_index(Index& index);

    // stores kmers of size kmer_size with stride over paths in graphs in the index
    void index_kmers(Index& index, int kmer_size, bool path_only, int edge_max, int stride = 1, 
                     bool allow_negatives = false);
    void for_each_kmer_parallel(
        const function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)>& lambda,
        int kmer_size, bool path_only, int edge_max, int stride, 
        bool allow_dups, bool allow_negatives = false);
    
    // Write out kmer lines to GCSA2
    void write_gcsa_out(ostream& out, int kmer_size,
                        bool path_only, bool forward_only,
                        int64_t head_id=0, int64_t tail_id=0);

    void write_gcsa_kmers_binary(ostream& out,
                                 int kmer_size,
                                 bool path_only, bool forward_only,
                                 int64_t head_id=0, int64_t tail_id=0);

    // gets all the kmers in GCSA's internal format.
    void get_gcsa_kmers(int kmer_size,
                        bool path_only, bool forward_only,
                        const function<void(vector<gcsa::KMer>&, bool)>& handle_kmers,
                        int64_t head_id=0, int64_t tail_id=0);

    vector<string> write_gcsa_kmers_binary(int kmer_size,
                                           bool path_only, bool forward_only,
                                           int64_t head_id=0, int64_t tail_id=0);
                              
    // Should we show our progress running through each graph?             
    bool show_progress = false;

private:

    void for_each_gcsa_kmer_position_parallel(int kmer_size,
                                              bool path_only, bool forward_only,
                                              int64_t& head_id, int64_t& tail_id,
                                              function<void(KmerPosition&)> lambda);
    
};

}


#endif

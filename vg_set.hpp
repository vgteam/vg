#ifndef VG_SET_H
#define VG_SET_H

#include <set>
#include <stdlib.h>
#include "gcsa.h"
#include "vg.hpp"
#include "index.hpp"
#include "hash_map.hpp"
#include "xg.hpp"


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

    // saves as a succinct, queryable representation
    void to_xg(const string& xg_db_name);

    // stores the nodes in the VGs identified by the filenames into the index
    void store_in_index(Index& index);
    void store_paths_in_index(Index& index);

    // stores kmers of size kmer_size with stride over paths in graphs in the index
    void index_kmers(Index& index, int kmer_size, int edge_max, int stride = 1, 
                     bool allow_negatives = false);
    void for_each_kmer_parallel(function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)>& lambda,
                                int kmer_size, int edge_max, int stride, 
                                bool allow_dups, bool allow_negatives = false);
    
    // Write out kmer lines to GCSA2
    void write_gcsa_out(ostream& out, int kmer_size, int edge_max, int stride,
                        int64_t head_id=0, int64_t tail_id=0);
    
    // gets all the kmers in GCSA's internal format.
    void get_gcsa_kmers(int kmer_size, int edge_max, int stride,
                        vector<gcsa::KMer>& kmers_out,
                        int64_t head_id=0, int64_t tail_id=0);

    bool show_progress;
    
private:
    
    // We create a struct that represents each kmer record we want to send to gcsa2
    struct KmerPosition {
        string kmer;
        string pos;
        set<char> prev_chars;
        set<char> next_chars;
        set<string> next_positions;
    };
    
    // We can loop over these in order to implement the other gcsa-related
    // functions above. GCSA kmers are the kmers in the graph with each node
    // existing in both its forward and reverse-complement orientation. Node IDs
    // in the GCSA graph are 2 * original node ID, +1 if the GCSA node
    // represents the reverse complement, and +0 if it does not. Non-reversing
    // edges link the forward copy of the from node to the forward copy of the
    // to node, and similarly for the reverse complement copies, while reversing
    // edges link the forward copy of the from node to the *reverse complement*
    // copy of the to node, and visa versa. This allows us to index both the
    // forward and reverse strands of every node, and to deal with GCSA's lack
    // of support for reversing edges, with the same trick. Note that
    // start_tail_id, if zero, will be replaced with the ID actually used for the
    // start/end node before lambda is ever called.
    void for_each_gcsa_kmer_position_parallel(int kmer_size, int edge_max, int stride,
                                              int64_t& head_id, int64_t& tail_id,
                                              function<void(KmerPosition&)> lambda);
    
};

}


#endif

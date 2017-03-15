#ifndef VG_SIMULATOR_H
#define VG_SIMULATOR_H

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "vg.hpp"
#include "xg.hpp"
#include "alignment.hpp"
#include "path.hpp"
#include "position.hpp"
#include "lru_cache.h"
#include "json2pb.h"

namespace vg {

using namespace std;

/**
 * Generate Alignments (with or without mutations, and in pairs or alone) from
 * an XG index.
 */
class Sampler {

public:

    xg::XG* xgidx;
    // We need this so we don't re-load the node for every character we visit in
    // it.
    LRUCache<id_t, Node> node_cache;
    mt19937 rng;
    int64_t nonce;
    // If set, only sample positions/start reads on the forward strands of their
    // nodes.
    bool forward_only;
    // A flag that we set if we don't want to generate sequences with Ns (on by dfault)
    bool no_Ns;
    Sampler(xg::XG* x,
            int seed = 0,
            bool forward_only = false,
            bool allow_Ns = false)
        : xgidx(x),
          node_cache(100),
          forward_only(forward_only),
          no_Ns(!allow_Ns),
          nonce(0) {
        if (!seed) {
            seed = time(NULL);
        }
        rng.seed(seed);
    }

    pos_t position(void);
    string sequence(size_t length);
    Alignment alignment(size_t length);
    Alignment alignment_with_error(size_t length,
                                   double base_error,
                                   double indel_error);
    vector<Alignment> alignment_pair(size_t read_length,
                                     size_t fragment_length,
                                     double fragment_std_dev,
                                     double base_error,
                                     double indel_error);
    size_t node_length(id_t id);
    char pos_char(pos_t pos);
    map<pos_t, char> next_pos_chars(pos_t pos);

    Alignment mutate(const Alignment& aln,
                     double base_error,
                     double indel_error);

    /**
     * Mutate the given edit, producing a vector of edits that should replace
     * it. Position is the position of the start of the edit, and is updated to
     * point to the next base after the mutated edit.
     */
    vector<Edit> mutate_edit(const Edit& edit,
                             const pos_t& position,
                             double base_error,
                             double indel_error,
                             const string& bases,
                             uniform_real_distribution<double>& rprob,
                             uniform_int_distribution<int>& rbase);

    string alignment_seq(const Alignment& aln);
    
    /// Return true if the alignment is semantically valid against the XG index
    /// we wrap, and false otherwise. Checks from_lengths on mappings to make
    /// sure all node bases are accounted for. Won't accept alignments with
    /// internal jumps between graph locations or regions; all skipped bases
    /// need to be accounted for by deletions.
    bool is_valid(const Alignment& aln);

};

}

#endif

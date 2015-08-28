#ifndef MAPPER_H
#define MAPPER_H

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "vg.hpp"
#include "xg.hpp"
#include "index.hpp"
#include "gcsa.h"
#include "alignment.hpp"
#include "path.hpp"
#include "json2pb.h"
#include "entropy.hpp"

namespace vg {

using namespace std;

class Mapper {


private:

    // Private constructor to delegate everything to. It might have all these
    // indexing structures null, for example if being called from the default
    // constructor.
    Mapper(Index* idex, gcsa::GCSA* g, xg::XG* xidex);

public:
    // Make a Mapper that pulls from a RocksDB index and optionally a GCSA2 kmer index.
    Mapper(Index* idex, gcsa::GCSA* g = nullptr);
    // Make a Mapper that pulls from an XG succinct graph and a GCSA2 kmer index.
    Mapper(xg::XG* xidex, gcsa::GCSA* g);
    Mapper(void);
    ~Mapper(void);
    Index* index;
    gcsa::GCSA* gcsa;
    xg::XG* xindex;

    // Align the given string and return an Alignment.
    Alignment align(string& seq, int kmer_size = 0, int stride = 0, int band_width = 1000);
    // Align the given read and return an aligned copy. Does not modify the input Alignment.
    Alignment align(Alignment& read, int kmer_size = 0, int stride = 0, int band_width = 1000);
    // Align the given read with multi-mapping. Returns the alignments in score
    // order, up to max_multimaps. Does not update the alignment passed in.
    // If the sequence is longer than the band_width, will only produce a single best banded alignment.
    // All alignments but the first are marked as secondary.
    vector<Alignment> align_multi(Alignment& aln, int kmer_size = 0, int stride = 0, int band_width = 1000);

    // Align read2 to the subgraph near the alignment of read1.
    // TODO: support banded alignment and intelligently use orientation heuristics
    void align_mate_in_window(Alignment& read1, Alignment& read2, int pair_window);

    // Return the one best banded alignment.
    Alignment align_banded(Alignment& read, int kmer_size = 0, int stride = 0, int band_width = 1000);

    // paired-end based
    
    // Both vectors of alignments will be sorted in order of increasing score.
    // All alignments but the first in each vector are marked as secondary.
    // Alignments at corresponding positions in the two vectors may or may not
    // be corresponding paired alignments. If a read does not map, its vector
    // will be empty.
    pair<vector<Alignment>, vector<Alignment>> 
        align_paired_multi(Alignment& read1,
                           Alignment& read2,
                           int kmer_size = 0,
                           int stride = 0,
                           int band_width = 1000,
                           int pair_window = 64);
    
    // Paired-end alignment ignoring multi-mapping. Returns either the two
    // highest-scoring reads if no rescue was required, or the highest-scoring
    // read and its corresponding rescue result if rescue was used.
    pair<Alignment, Alignment> align_paired(Alignment& read1,
                                            Alignment& read2,
                                            int kmer_size = 0,
                                            int stride = 0,
                                            int band_width = 1000,
                                            int pair_window = 64);

    // base algorithm for above Update the passed-in Alignment with a highest-
    // score alignment, and return all good alignments sorted by score up to
    // max_multimaps. If the read does not map, the returned vector will be
    // empty. No alignments will be marked as secondary; the caller must do that
    // if they plan to produce GAM output.
    vector<Alignment> align_threaded(Alignment& read,
                                     int& hit_count,
                                     int kmer_size = 0,
                                     int stride = 0,
                                     int attempt = 0);

    // not used
    Alignment& align_simple(Alignment& alignment, int kmer_size = 0, int stride = 0);

    set<int> kmer_sizes;
    bool debug;
    int best_clusters;
    int cluster_min;
    int hit_max;
    int hit_size_threshold;
    int kmer_min;
    int kmer_threshold;
    int max_thread_gap;
    int kmer_sensitivity_step;
    int thread_extension;
    int thread_extension_max;
    int context_depth;
    int max_multimaps;
    int max_attempts;
    int softclip_threshold;
    float target_score_per_bp;
    bool prefer_forward;
    bool greedy_accept;
    float min_kmer_entropy;

};

// utility
int softclip_start(Alignment& alignment);
int softclip_end(Alignment& alignment);
const vector<string> balanced_kmers(const string& seq, int kmer_size, int stride);


}

#endif

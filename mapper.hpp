#ifndef MAPPER_H
#define MAPPER_H

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "vg.hpp"
#include "index.hpp"
#include "alignment.hpp"
#include "path.hpp"
#include "pb2json.h"

namespace vg {

using namespace std;

class Mapper {

public:

    Mapper(Index* idex);
    Mapper(void) : index(NULL), best_clusters(0) { }
    ~Mapper(void);
    Index* index;

    Alignment align(string& seq, int kmer_size = 0, int stride = 0, int band_width = 1000);
    Alignment align(Alignment& read, int kmer_size = 0, int stride = 0, int band_width = 1000);

    void align_mate_in_window(Alignment& read1, Alignment& read2, int pair_window);

    Alignment align_banded(Alignment& read, int kmer_size = 0, int stride = 0, int band_width = 1000);

    // paired-end based
    pair<Alignment, Alignment> align_paired(Alignment& read1,
                                            Alignment& read2,
                                            int kmer_size = 0,
                                            int stride = 0,
                                            int band_width = 1000,
                                            int pair_window = 64);

    // base algorithm for above
    Alignment& align_threaded(Alignment& read,
                              int& hit_count,
                              int kmer_size = 0,
                              int stride = 0,
                              int attempt = 0);

    // not used
    Alignment& align_simple(Alignment& alignment, int kmer_size = 0, int stride = 0);

    set<int> kmer_sizes;
    bool debug;
    int best_clusters;
    int hit_max;
    int hit_size_threshold;
    int kmer_min;
    int kmer_threshold;
    int kmer_sensitivity_step;
    int thread_extension;
    int thread_extension_max;
    int max_attempts;
    int softclip_threshold;
    float target_score_per_bp;
    bool prefer_forward;

};

// utility
int softclip_start(Alignment& alignment);
int softclip_end(Alignment& alignment);
const vector<string> balanced_kmers(const string& seq, int kmer_size, int stride);


}

#endif

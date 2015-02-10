#ifndef MAPPER_H
#define MAPPER_H

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "vg.hpp"
#include "index.hpp"
#include "pb2json.h"

namespace vg {

using namespace std;

class Mapper {

public:

    Mapper(Index* idex);
    Mapper(void) : index(NULL), best_clusters(0) { }
    ~Mapper(void);
    Index* index;

    Alignment align(string& seq, int kmer_size = 0, int stride = 0);
    Alignment& align_threaded(Alignment& read, int kmer_size = 0, int stride = 0,
                              int attempt = 0, int hit_count = 0);
    Alignment& align_simple(Alignment& alignment, int kmer_size = 0, int stride = 0);

    set<int> kmer_sizes;

    int best_clusters;
    int hit_max;
    int hit_size_threshold;
    int kmer_min;
    int thread_extension;
    int thread_extension_max;

};

// utility
int softclip_start(Alignment& alignment);
int softclip_end(Alignment& alignment);
const vector<string> balanced_kmers(const string& seq, int kmer_size, int stride);


}

#endif

#ifndef MAPPER_H
#define MAPPER_H

#include <iostream>
#include <map>
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

    Alignment align(string& seq, int stride);
    Alignment& align(Alignment& read, int stride);

    set<int> kmer_sizes;
    vector<string> kmers_of(const string& seq, int stride = 1);

    int best_clusters;

};

}

#endif

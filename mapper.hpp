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
    Mapper(void) {
        index = NULL;
    }
    ~Mapper(void);
    Index* index;

    Alignment align(string& seq);
    Alignment& align(Alignment& read);

    set<int> kmer_sizes;
    set<string> kmers_of(const string& seq);

};

}

#endif

#ifndef MAPPER_H
#define MAPPER_H

#include <iostream>
#include <map>
#include "vg.h"
#include "index.h"
#include "pb2json.h"
#include "leveldb/db.h"

namespace vg {

using namespace std;

class Mapper {

public:

    Mapper(Index* idex) {
        index = idex;
    }
    Mapper(void) {
        index = NULL;
    }
    ~Mapper(void);
    Index* index;

    Alignment align(string& seq, int kmer_size);
    Alignment& align(Alignment& read, int kmer_size);

};

}

#endif

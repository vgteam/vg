#ifndef PILEUP_H
#define PILEUP_H

#include <iostream>
#include <algorithm>
#include <functional>
#include <set>
#include <list>
#include <unordered_set>
#include "vg.pb.h"
#include "vg.hpp"
#include "hash_map.hpp"

namespace vg {

using namespace std;

// This is a collection of protobuf Pileup records that are indexed
// on their position. Pileups can be merged and streamed, and computed
// from Alignments.  The pileup records themselves are essentially
// protobuf versions of lines in Samtools pileup format.  
class Pileups {
public:
    
    Pileups() {}

    // from protobuf stream
    Pileups(istream& in, bool showp = false);
    
    // copy constructor
    Pileups(const Pileups& other) {
        if (this != &other) {
            for (auto& p : other._pileups) {
                insert(new Pileup(*p));
            }
        }
    }

    // move constructor
    Pileups(Pileups&& other) noexcept {
        _pileups = other._pileups;
        other._pileups.clear();
    }

    // copy assignment operator
    Pileups& operator=(const Pileups& other) {
        Pileups tmp(other);
        *this = move(tmp);
        return *this;
    }

    // move assignment operator
    Pileups& operator=(Pileups&& other) noexcept {
        swap(_pileups, other._pileups);
        other._pileups.clear();
        return *this;
    }

    // delete contents of table
    ~Pileups() {
        clear();
    }
    void clear();

    // hash a Pileup* by its position (using pair hash from hash_mapp.hpp)
    // don't store NULLs!!
    typedef pair<google::protobuf::int64, google::protobuf::int64> PBI64Pair;
    struct PileupHashFn {
        size_t operator()(const Pileup* p) const {
            return std::hash<PBI64Pair>()(PBI64Pair(p->position().node_id(),
                                                    p->position().offset()));
        }
    };
    typedef unordered_set<Pileup*, PileupHashFn> PileupHash;
    
    // This maps from Position to Pileup.
    PileupHash _pileups;

    // write to JSON
    void to_json(ostream& out);
    // read from protobuf
    void load(istream& in);
    // write to protobuf
    void write(ostream& out, uint64_t buffer_size = 1000);

    // apply function to each pileup in table
    void for_each(function<void(Pileup&)>& lambda);

    // lookup pileup. by using a set, we require having a Pileup object
    // in hand to search (map would let us search by Position which is nicer).
    // memory saved should be worth it tho..
    Pileup* get(const Pileup& pileup) {
        Pileup* t = const_cast<Pileup*>(&pileup);
        PileupHash::const_iterator p = _pileups.find(t);
        return p != _pileups.end() ? *p : NULL;
    }

    // get a pileup.  if it's null, create a new one and insert it.
    Pileup* get_create(const Pileup& pileup) {
        Pileup* p = get(pileup);
        if (p == NULL) {
            p = new Pileup(pileup);
            insert(p);
        }
        return p;
    }   

    // insert a pileup into the table. it will be deleted by ~Pileups()!!!
    // return true if new pileup inserted, false if merged into existing one
    bool insert(Pileup* pileup);

    // create / update all pileups from a single alignment
    void compute_from_alignment(VG& graph, Alignment& alignment);

    // create / update all pileups from an edit (called by above).
    // query stores the current position (and nothing else).  
   void compute_from_edit(Pileup& query, int64_t& read_offset,
                          const Node& node, const Alignment& alignment,
                          const Mapping& mapping, const Edit& edit);
            
    // move all entries in other object into this one.
    // if two positions collide, they are merged.
    // other will be left empty. this is returned
    Pileups& merge(Pileups& other);
};

// merge p2 into p1 and return 1
Pileup& merge_pileup(Pileup& p1, const Pileup& p2);

}

#endif

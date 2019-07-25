#ifndef VG_MEM_HPP_INCLUDED
#define VG_MEM_HPP_INCLUDED

#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>
#include "position.hpp"
#include "utility.hpp"

#include <functional>
#include <string>
#include <vector>
#include <map>


/**
 * \file mem.hpp
 *
 * Maximal Exact Matches (MEMs) header
 */

namespace vg {
    
using namespace std;

class MaximalExactMatch {

public:

    //const string* source;
    string::const_iterator begin;
    string::const_iterator end;
    gcsa::range_type range;
    size_t match_count;
    size_t queried_count;
    int fragment;
    bool primary; // if not a sub-MEM
    std::vector<gcsa::node_type> nodes;
    unordered_map<path_handle_t, vector<pair<size_t, bool> > > positions;
    
    MaximalExactMatch(string::const_iterator b,
                      string::const_iterator e,
                      gcsa::range_type r,
                      size_t m = 0)
        : begin(b), end(e), range(r), match_count(m) { }

    // construct the sequence of the MEM; useful in debugging
    string sequence(void) const;
    // get the length of the MEM
    int length(void) const;
    // tells if the MEM contains an N
    size_t count_Ns(void) const;
    // filter out every ceil(nodes.size()/n) hits
    size_t filter_hits_to(int limit);

    friend bool operator==(const MaximalExactMatch& m1, const MaximalExactMatch& m2);
    friend bool operator<(const MaximalExactMatch& m1, const MaximalExactMatch& m2);
    friend ostream& operator<<(ostream& out, const MaximalExactMatch& m);

    MaximalExactMatch(void) = default;                                      // Copy constructor
    MaximalExactMatch(const MaximalExactMatch&) = default;               // Copy constructor
    MaximalExactMatch(MaximalExactMatch&&) = default;                    // Move constructor
    MaximalExactMatch& operator=(const MaximalExactMatch&) & = default;  // Copy assignment operator
    MaximalExactMatch& operator=(MaximalExactMatch&&) & = default;       // Move assignment operator
    //virtual ~MaximalExactMatch() { }                     // Destructor
};

const string mems_to_json(const vector<MaximalExactMatch>& mems);

// helpers for computing the number of bases in the query covered by a cluster
vector<string::const_iterator> cluster_cover(const vector<MaximalExactMatch>& cluster);
int cluster_coverage(const vector<MaximalExactMatch>& cluster);
// helper to tell if mems are ovelapping
bool mems_overlap(const MaximalExactMatch& mem1,
                  const MaximalExactMatch& mem2);
// distance of overlap, or 0 if there is no overlap
int mems_overlap_length(const MaximalExactMatch& mem1,
                        const MaximalExactMatch& mem2);
// helper to tell if clusters have any overlap
bool clusters_overlap_in_read(const vector<MaximalExactMatch>& cluster1,
                              const vector<MaximalExactMatch>& cluster2);
bool clusters_overlap_in_graph(const vector<MaximalExactMatch>& cluster1,
                               const vector<MaximalExactMatch>& cluster2);
int clusters_overlap_length(const vector<MaximalExactMatch>& cluster1,
                            const vector<MaximalExactMatch>& cluster2);
vector<pos_t> cluster_nodes(const vector<MaximalExactMatch>& cluster);
vector<MaximalExactMatch> translate_mems(const vector<MaximalExactMatch>& mems,
                                         const unordered_map<id_t, pair<id_t, bool> >& trans);

// returns the min distance and the relative orientation implied by the position annotations on the mems
pair<int64_t, int64_t> mem_min_oriented_distances(const MaximalExactMatch& m1, const MaximalExactMatch& m2);


}

#endif

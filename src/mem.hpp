#ifndef VG_MEM_HPP_INCLUDED
#define VG_MEM_HPP_INCLUDED

#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"
#include "position.hpp"

#include <functional>
#include <string>
#include <vector>
#include <map>


/**
 * \file mem.hpp
 *
 * Maximal Exact Matches (MEMs) and the chaining and clustering tools to work
 * with them.
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
    int fragment;
    bool primary; // if not a sub-MEM
    std::vector<gcsa::node_type> nodes;
    map<string, vector<size_t> > positions;
    
    MaximalExactMatch(string::const_iterator b,
                      string::const_iterator e,
                      gcsa::range_type r,
                      size_t m = 0)
        : begin(b), end(e), range(r), match_count(m) { }

    // construct the sequence of the MEM; useful in debugging
    string sequence(void) const;
    // get the length of the MEM
    int length(void) const;
    // uses a callback to fill out the MEM positions
    void fill_positions(const function<map<string, vector<size_t>>(gcsa::node_type)>& node_positions_in_paths);
    // tells if the MEM contains an N
    size_t count_Ns(void) const;

    friend bool operator==(const MaximalExactMatch& m1, const MaximalExactMatch& m2);
    friend bool operator<(const MaximalExactMatch& m1, const MaximalExactMatch& m2);
    friend ostream& operator<<(ostream& out, const MaximalExactMatch& m);

    MaximalExactMatch(void) = default;                                      // Copy constructor
    MaximalExactMatch(const MaximalExactMatch&) = default;               // Copy constructor
    MaximalExactMatch(MaximalExactMatch&&) = default;                    // Move constructor
    MaximalExactMatch& operator=(const MaximalExactMatch&) & = default;  // MaximalExactMatchopy assignment operator
    MaximalExactMatch& operator=(MaximalExactMatch&&) & = default;       // Move assignment operator
    //virtual ~MaximalExactMatch() { }                     // Destructor
};

const string mems_to_json(const vector<MaximalExactMatch>& mems);

// helper for computing the number of bases in the query covered by a cluster
int cluster_coverage(const vector<MaximalExactMatch>& cluster);
// helper to tell if mems are ovelapping
bool mems_overlap(const MaximalExactMatch& mem1,
                  const MaximalExactMatch& mem2);
// distance of overlap, or 0 if there is no overlap
int mems_overlap_length(const MaximalExactMatch& mem1,
                        const MaximalExactMatch& mem2);
// helper to tell if clusters have any overlap
bool clusters_overlap(const vector<MaximalExactMatch>& cluster1,
                      const vector<MaximalExactMatch>& cluster2);

class MEMChainModelVertex {
public:
    MaximalExactMatch mem;
    vector<pair<MEMChainModelVertex*, double> > next_cost; // for forward
    vector<pair<MEMChainModelVertex*, double> > prev_cost; // for backward
    double weight;
    double score;
    int approx_position;
    MEMChainModelVertex* prev;
    MEMChainModelVertex(void) = default;                                      // Copy constructor
    MEMChainModelVertex(const MEMChainModelVertex&) = default;               // Copy constructor
    MEMChainModelVertex(MEMChainModelVertex&&) = default;                    // Move constructor
    MEMChainModelVertex& operator=(const MEMChainModelVertex&) & = default;  // MEMChainModelVertexopy assignment operator
    MEMChainModelVertex& operator=(MEMChainModelVertex&&) & = default;       // Move assignment operator
    virtual ~MEMChainModelVertex() { }                     // Destructor
};

class MEMChainModel {
public:
    vector<MEMChainModelVertex> model;
    map<int, vector<vector<MEMChainModelVertex>::iterator> > approx_positions;
    set<vector<MEMChainModelVertex>::iterator> redundant_vertexes;
    MEMChainModel(
        const vector<size_t>& aln_lengths,
        const vector<vector<MaximalExactMatch> >& matches,
        const function<int(pos_t)>& approx_position,
        const function<double(const MaximalExactMatch&, const MaximalExactMatch&)>& transition_weight,
        int band_width = 10,
        int position_depth = 1,
        int max_connections = 10);
    void score(const set<MEMChainModelVertex*>& exclude);
    MEMChainModelVertex* max_vertex(void);
    vector<vector<MaximalExactMatch> > traceback(int alt_alns, bool paired, bool debug);
    void display(ostream& out);
    void clear_scores(void);
};

}

#endif

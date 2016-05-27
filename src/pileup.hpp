#ifndef PILEUP_H
#define PILEUP_H

#include <iostream>
#include <algorithm>
#include <functional>
#include "vg.pb.h"
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"

namespace vg {

using namespace std;

// This is a collection of protobuf NodePileup records that are indexed
// on their position, as well as EdgePileup records.
// Pileups can be merged and streamed, and computed
// from Alignments.  The pileup records themselves are essentially
// protobuf versions of lines in Samtools pileup format, with deletions
// represented using a graph-based notation. 
class Pileups {
public:
    
    Pileups(VG* graph, int min_quality = 0, int max_mismatches = 1, int window_size = 0,
            int max_depth = 1000) :
        _graph(graph),
        _min_quality(min_quality),
        _max_mismatches(max_mismatches),
        _window_size(window_size),
        _max_depth(max_depth),
        _min_quality_count(0),
        _max_mismatch_count(0),
        _bases_count(0)
{}
    
    // copy constructor
    Pileups(const Pileups& other) {
        if (this != &other) {
            _graph = other._graph;
            for (auto& p : other._node_pileups) {
                insert_node_pileup(new NodePileup(*p.second));
            }
            _min_quality = other._min_quality;
            _max_mismatches = other._max_mismatches;
            _window_size = other._window_size;
            _max_depth = other._max_depth;
            _min_quality_count = other._min_quality_count;
            _max_mismatch_count = other._max_mismatch_count;
            _bases_count = other._bases_count;
        }
    }

    // move constructor
    Pileups(Pileups&& other) noexcept {
        _graph = other._graph;
        _node_pileups = other._node_pileups;
        other._node_pileups.clear();
        _min_quality = other._min_quality;
        _max_mismatches = other._max_mismatches;
        _window_size = other._window_size;
        _max_depth = other._max_depth;
        _min_quality_count = other._min_quality_count;
        _max_mismatch_count = other._max_mismatch_count;
        _bases_count = other._bases_count;
    }

    // copy assignment operator
    Pileups& operator=(const Pileups& other) {
        Pileups tmp(other);
        *this = move(tmp);
        return *this;
    }

    // move assignment operator
    Pileups& operator=(Pileups&& other) noexcept {
        _graph = other._graph;
        swap(_node_pileups, other._node_pileups);
        other._node_pileups.clear();
        _min_quality = other._min_quality;
        _max_mismatches = other._max_mismatches;
        _window_size = other._window_size;
        _max_depth = other._max_depth;
        _min_quality_count = other._min_quality_count;
        _max_mismatch_count = other._max_mismatch_count;
        _bases_count = other._bases_count;
        return *this;
    }

    // delete contents of table
    ~Pileups() {
        clear();
    }
    void clear();

    typedef hash_map<int64_t, NodePileup*> NodePileupHash;
    typedef pair_hash_map<pair<NodeSide, NodeSide>, EdgePileup*> EdgePileupHash;

    VG* _graph;
    
    // This maps from Position to Pileup.
    NodePileupHash _node_pileups;
    EdgePileupHash _edge_pileups;

    // Ignore bases with quality less than this
    int _min_quality;
    // max mismatches within window_size
    int _max_mismatches;
    // number of bases to scan in each direction for mismatches
    int _window_size;
    // prevent giant protobufs
    int _max_depth;
    // Keep count of bases filtered by quality
    mutable uint64_t _min_quality_count;
    // keep count of bases filtered by mismatches
    mutable uint64_t _max_mismatch_count;
    // overall count for perspective on above
    mutable uint64_t _bases_count;

    // write to JSON
    void to_json(ostream& out);
    // read from protobuf
    void load(istream& in);
    // write to protobuf
    void write(ostream& out, uint64_t buffer_size = 5);

    // apply function to each pileup in table
    void for_each_node_pileup(const function<void(NodePileup&)>& lambda);

    // search hash table for node id
    NodePileup* get_node_pileup(int64_t node_id) {
        auto p = _node_pileups.find(node_id);
        return p != _node_pileups.end() ? p->second : NULL;
    }
        
    // get a pileup.  if it's null, create a new one and insert it.
    NodePileup* get_create_node_pileup(const Node* node) {
      NodePileup* p = get_node_pileup(node->id());
        if (p == NULL) {
            p = new NodePileup();
            p->set_node_id(node->id());
            for (int i = 0; i < node->sequence().length(); ++i) {
                BasePileup* b = p->add_base_pileup();
                b->set_num_bases(0);
                b->set_ref_base((int)node->sequence()[i]);
            }
            _node_pileups[node->id()] = p;
        }
        return p;
    }

    void for_each_edge_pileup(const function<void(EdgePileup&)>& lambda);

    // search hash table for edge id
    EdgePileup* get_edge_pileup(pair<NodeSide, NodeSide> sides);
            
    // get a pileup.  if it's null, create a new one and insert it.
    EdgePileup* get_create_edge_pileup(pair<NodeSide, NodeSide> sides);
    
    void extend(Pileup& pileup);

    // insert a pileup into the table. it will be deleted by ~Pileups()!!!
    // return true if new pileup inserted, false if merged into existing one
    bool insert_node_pileup(NodePileup* pileup);
    bool insert_edge_pileup(EdgePileup* edge_pileup);
    
    // create / update all pileups from a single alignment
    void compute_from_alignment(Alignment& alignment);

    // create / update all pileups from an edit (called by above).
    // query stores the current position (and nothing else).  
    void compute_from_edit(NodePileup& pileup, int64_t& node_offset, int64_t& read_offset,
                           const Node& node, const Alignment& alignment,
                           const Mapping& mapping, const Edit& edit,
                           const Edit* next_edit,
                           const vector<int>& mismatch_counts,
                           pair<const Mapping*, int64_t>& last_match,
                           pair<const Mapping*, int64_t>& last_del,
                           pair<const Mapping*, int64_t>& open_del);

    // do one pass to count all mismatches in read, so we can do
    // mismatch filter efficiently in 2nd path.
    // mismatches[i] stores number of mismatches in range (0, i)
    static void count_mismatches(VG& graph, const Path& path, vector<int>& mismatches,
                                 bool skipIndels = false);

    // check base quality as well as miss match filter
    bool pass_filter(const Alignment& alignment, int64_t read_offset,
                     int64_t length,
                     const vector<int>& mismatches) const;
            
    // move all entries in other object into this one.
    // if two positions collide, they are merged.
    // other will be left empty. this is returned
    Pileups& merge(Pileups& other);

    // merge p2 into p1 and return 1. p2 is left an empty husk
    BasePileup& merge_base_pileups(BasePileup& p1, BasePileup& p2);

    // merge p2 into p1 and return 1. p2 is lef an empty husk
    NodePileup& merge_node_pileups(NodePileup& p1, NodePileup& p2);
    
    // merge p2 into p1 and return 1. p2 is lef an empty husk
    EdgePileup& merge_edge_pileups(EdgePileup& p1, EdgePileup& p2);

    // get ith BasePileup record
    static BasePileup* get_base_pileup(NodePileup& np, int64_t offset) {
        assert(offset < np.base_pileup_size());
        return np.mutable_base_pileup(offset);
    }
    static const BasePileup* get_base_pileup(const NodePileup& np, int64_t offset) {
        assert(offset < np.base_pileup_size());
        return &np.base_pileup(offset);
    }

    // get ith BasePileup record, create if doesn't exist
    static BasePileup* get_create_base_pileup(NodePileup& np, int64_t offset) {
        for (int64_t i = np.base_pileup_size(); i <= offset; ++i) {
            np.add_base_pileup();
        }
        return get_base_pileup(np, offset);
    }

    // the bases string in BasePileup doesn't allow random access.  This function
    // will parse out all the offsets of snps, insertions, and deletions
    // into one array, each offset is a pair of indexes in the bases and qualities arrays
    static void parse_base_offsets(const BasePileup& bp,
                                   vector<pair<int64_t, int64_t> >& offsets);

    // transform case of every character in string
    static void casify(string& seq, bool is_reverse);

    // make the sam pileup style token
    static void make_match(string& seq, int64_t from_length, bool is_reverse);
    static void make_insert(string& seq, bool is_reverse);
    static void make_delete(string& seq, bool is_reverse,
                            const pair<const Mapping*, int64_t>& last_match,
                            const Mapping& mapping, int64_t node_offset);
    static void make_delete(string& seq, bool is_reverse,
                            int64_t from_id, int64_t from_offset, bool from_start,
                            int64_t to_id, int64_t to_offset, bool to_end);
        
    static void parse_insert(const string& tok, int64_t& len, string& seq, bool& is_reverse);
    static void parse_delete(const string& tok, bool& is_reverse,
                             int64_t& from_id, int64_t& from_offset, bool& from_start,
                             int64_t& to_id, int64_t& to_offset, bool& to_end);

    static bool base_equal(char c1, char c2, bool is_reverse);
    
    // get a pileup value on forward strand
    static char extract_match(const BasePileup& bp, int64_t offset);

    // get arbitrary value from offset on forward strand
    static string extract(const BasePileup& bp, int64_t offset);
};



}

#endif

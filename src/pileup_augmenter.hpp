// Augments a graph using a pileup (made with vg pileup)

#ifndef VG_CALLER_HPP_INCLUDED
#define VG_CALLER_HPP_INCLUDED

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include <vg/vg.pb.h>
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "pileup.hpp"
#include "path_index.hpp"
#include "genotypekit.hpp"
#include "option.hpp"

namespace vg {

using namespace std;

// container for storing pairs of support for calls (value for each strand)
struct StrandSupport {
    int fs; // forward support
    int rs; // reverse support
    double qual; // phred score (derived from sum log p-err over all observations)
    StrandSupport(int f = 0, int r = 0, double q = 0) :
        fs(f), rs(r), qual(q) {}
    bool operator<(const StrandSupport& other) const {
        if ((fs + rs) == (other.fs + other.rs)) {
            // more strand bias taken as less support
            return abs(fs - rs) > abs(other.fs - rs);
        } 
        return fs + rs < other.fs + other.rs;
    }
    bool operator>=(const StrandSupport& other) const {
        return !(*this < other);
    }
    bool operator==(const StrandSupport& other) const {
        return fs == other.fs && rs == other.rs && qual == other.qual;
    }
    // min out at 0
    StrandSupport operator-(const StrandSupport& other) const {
        return StrandSupport(max(0, fs - other.fs), max(0, rs - other.rs),
                             max(0., qual - other.qual));
    }
    StrandSupport& operator+=(const StrandSupport& other) {
        fs += other.fs;
        rs += other.rs;
        qual += other.qual;
        return *this;
    }
    int total() { return fs + rs; }
};

inline StrandSupport minSup(vector<StrandSupport>& s) {
    if (s.empty()) {
        return StrandSupport();
    }
    return *min_element(s.begin(), s.end());
}
inline StrandSupport maxSup(vector<StrandSupport>& s) {
    if (s.empty()) {
        return StrandSupport();
    }
    return *max_element(s.begin(), s.end());
}
inline StrandSupport avgSup(vector<StrandSupport>& s) {
    StrandSupport ret;
    if (!s.empty()) {
        for (auto sup : s) {
            ret.fs += sup.fs;
            ret.rs += sup.rs;
            ret.qual += sup.qual;
        }
        ret.fs /= s.size();
        ret.rs /= s.size();
        ret.qual /= s.size();
    }
    return ret;
}
inline StrandSupport totalSup(vector<StrandSupport>& s) {
    StrandSupport ret;
    if (!s.empty()) {
        for (auto sup : s) {
            ret.fs += sup.fs;
            ret.rs += sup.rs;
            ret.qual += sup.qual;
        }
    }
    return ret;
}

inline ostream& operator<<(ostream& os, const StrandSupport& sup) {
    return os << sup.fs << ", " << sup.rs << ", " << sup.qual;
}

// We need to break apart nodes but remember where they came from to update edges.
// Wrap all this up in this class.  For a position in the input graph, we can have
// up to three nodes in the augmented graph (Ref, Alt1, Alt2), so we map to node
// triplets (Entry struct below).  Note we never *call* all three nodes due to
// diploid assumption, but the augmented graph stores everything. 
struct NodeDivider {
    // up to three fragments per position in augmented graph (basically a Node 3-tuple,
    // avoiding aweful C++ tuple syntax)
    enum EntryCat {Ref = 0, Alt1, Alt2, Last};
    struct Entry { Entry(Node* r = 0, vector<StrandSupport> sup_r = vector<StrandSupport>(),
                         Node* a1 = 0, vector<StrandSupport> sup_a1 = vector<StrandSupport>(),
                         Node* a2 = 0, vector<StrandSupport> sup_a2 = vector<StrandSupport>()) : ref(r), alt1(a1), alt2(a2),
                                                               sup_ref(sup_r), sup_alt1(sup_a1), sup_alt2(sup_a2){}
        Node* ref; Node* alt1; Node* alt2;
        vector<StrandSupport> sup_ref;
        vector<StrandSupport> sup_alt1;
        vector<StrandSupport> sup_alt2;
        Node*& operator[](int i) {
            assert(i >= 0 && i <= 2);
            return i == EntryCat::Ref ? ref : (i == EntryCat::Alt1 ? alt1 : alt2);
        }
        vector<StrandSupport>& sup(int i) {
            assert(i >= 0 && i <= 2);
            return i == EntryCat::Ref ? sup_ref : (i == EntryCat::Alt1 ? sup_alt1 : sup_alt2);
        }
    };
    // offset in original graph node -> up to 3 nodes in call graph
    typedef map<int, Entry> NodeMap;
    // Node id in original graph to map above
    typedef hash_map<int64_t, NodeMap> NodeHash;
    NodeHash index;
    int64_t* _max_id;
    // map given node to offset i of node with id in original graph
    // this function can never handle overlaps (and should only be called before break_end)
    void add_fragment(const Node* orig_node, int offset, Node* subnode, EntryCat cat, vector<StrandSupport> sup);
    // break node if necessary so that we can attach edge at specified side
    // this function wil return NULL if there's no node covering the given location
    Entry break_end(const Node* orig_node, VG* graph, int offset, bool left_side);
    // assuming input node is fully covered, list of nodes that correspond to it in call graph
    // if node not in structure at all, just return input (assumption uncalled nodes kept as is)
    list<Mapping> map_node(int64_t node_id, int64_t start_offset, int64_t length, bool reverse);
    // erase everything (but don't free any Node pointers, they belong to the graph)
    void clear();
};
ostream& operator<<(ostream& os, const NodeDivider::NodeMap& nm);
ostream& operator<<(ostream& os, NodeDivider::Entry entry);

/**
 * Super simple graph augmentor/caller.
 * Idea: Idependently process Pileup records, using simple model to make calls that
 *       take into account read errors with diploid assumption.  Edges and node positions
 *       are called independently for now.  
 * Outputs either a sample graph (only called nodes and edges) or augmented graph
 * (include uncalled nodes and edges too).
 */
class PileupAugmenter {
public:

    // log of zero
    static const double Log_zero;
    // use this score when pileup is missing quality
    static const char Default_default_quality;
    // don't augment graph without minimum support
    static const int Default_min_aug_support;
    
    PileupAugmenter(VG* graph,
           int default_quality = Default_default_quality,
           int min_aug_support = Default_min_aug_support);

   ~PileupAugmenter();
    void clear();

    // input graph
    VG* _graph;
    // Output augmented graph with annotations
    SupportAugmentedGraph _augmented_graph;

    // buffer for base calls for each position in the node
    // . = reference
    // - = missing
    typedef pair<string, string> Genotype;
    vector<Genotype> _node_calls;
    vector<pair<StrandSupport, StrandSupport> > _node_supports;
    // separate structure for isnertion calls since they
    // don't really have reference coordinates (instead happen just to
    // right of offset).  
    vector<Genotype> _insert_calls;
    vector<pair<StrandSupport, StrandSupport> > _insert_supports;
    // buffer for current node;
    const Node* _node;
    // max id in call_graph
    int64_t _max_id;
    // link called nodes back to the original graph. needed
    // to figure out where edges go
    NodeDivider _node_divider;
    unordered_set<int64_t> _visited_nodes;
    unordered_map<pair<NodeSide, NodeSide>, StrandSupport> _called_edges; // map to support
    // deletes can don't necessarily need to be in incident to node ends
    // so we throw in an offset into the mix. 
    typedef pair<NodeSide, int> NodeOffSide;
    // map a call category to an edge
    typedef unordered_map<pair<NodeOffSide, NodeOffSide>, char> EdgeHash;
    EdgeHash _augmented_edges;
    // keep track of inserted nodes for tsv output
    struct InsertionRecord {
        Node* node;
        StrandSupport sup;
        int64_t orig_id;
        int orig_offset;
    };
    typedef unordered_map<int64_t, InsertionRecord> InsertionHash;
    InsertionHash _inserted_nodes;
    // hack for better estimating support for edges that go around
    // insertions (between the adjacent ref nodes)
    typedef unordered_map<pair<NodeOffSide, NodeOffSide>, StrandSupport> EdgeSupHash;
    EdgeSupHash _insertion_supports;

    // need to keep track of support for augmented deletions
    // todo: generalize augmented edge support
    EdgeSupHash _deletion_supports;

    // maximum number of nodes to call before writing out output stream
    int _buffer_size;
    // if we don't have a mapping quality for a read position, use this
    char _default_quality;
    // minimum support to augment graph
    int _min_aug_support;

    // write the augmented graph
    void write_augmented_graph(ostream& out, bool json);

    // call every position in the node pileup
    void call_node_pileup(const NodePileup& pileup);

    // call an edge.  remembering it in a table for the whole graph
    void call_edge_pileup(const EdgePileup& pileup);

    // fill in edges in the augmented graph (those that are incident to 2 call
    // nodes) and add uncalled nodes (optionally)
    void update_augmented_graph();

    // map a path (can have edits, ie from Alignment) from base graph to augmented graph
    // aug_path parameter is empty path that will be written to
    void map_path(const Path& base_path, list<mapping_t>& aug_path, bool expect_edits);

    // Apply edits from base_mapping to corresponding augmented mappings that share same
    // from interval, but don't yet have edits (called by map_path);
    void apply_mapping_edits(const Mapping& base_mapping, list<Mapping>& aug_mappings);

    // TODO:
    // method to normalize mapped paths back onto the augmented graph.  ie check each
    // non-match edit to see if it can be turned into a match on the augmented graph.

    // map paths from input graph into called (augmented) graph
    void map_paths();
    // make sure mapped paths generate same strings as input paths
    void verify_path(const Path& in_path, const list<mapping_t>& call_path);
    
    // call position at given base
    // if insertion flag set to true, call insertion between base and next base
    void call_base_pileup(const NodePileup& np, int64_t offset, bool insertions);
    
    // Find the top-two bases in a pileup, along with their counts
    // Last param toggles whether we consider only inserts or everything else
    // (do not compare all at once since inserts do not have reference coordinates)
    void compute_top_frequencies(const BasePileup& bp,
                                 const vector<pair<int64_t, int64_t> >& base_offsets,
                                 string& top_base, int& top_count, int& top_rev_count,
                                 string& second_base, int& second_count, int& second_rev_count,
                                 int& total_count, bool inserts);

    // Sum up the qualities of a given symbol in a pileup
    double total_base_quality(const BasePileup& pb,
                              const vector<pair<int64_t, int64_t> >& base_offsets,
                              const string& val);

    // write graph structure corresponding to all the calls for the current
    // node.  
    void create_node_calls(const NodePileup& np);

    void create_augmented_edge(Node* node1, int from_offset, bool left_side1, bool aug1,
                               Node* node2, int to_offset, bool left_side2, bool aug2, char cat,
                               StrandSupport support);

    // Annotate nodes and edges in the augmented graph with call info.
    void annotate_augmented_node(Node* node, char call, StrandSupport support, int64_t orig_id, int orig_offset);
    void annotate_augmented_edge(Edge* edge, char call, StrandSupport support);
    void annotate_augmented_nodes();

    // Add nodes that are passed through as-is (ie not augmented at all) to the translation table
    void annotate_non_augmented_nodes();

    // log function that tries to avoid 0s
    static double safe_log(double v) {
        return v == 0. ? Log_zero : ::log10(v);
    }

    // call missing
    static bool missing_call(const Genotype& g) {
        return g.first == "-" &&  g.second == "-";
    }

    // call is reference
    static bool ref_call(const Genotype& g) {
        return g.first == "." && (g.second == "." || g.second == "-");
    }

    // classify call as 0: missing 1: reference 2: snp
    // (holdover from before indels)
    static int call_cat(const Genotype&g) {
        if (missing_call(g)) {
            return 0;
        } else if (ref_call(g)) {
            return 1;
        }
        return 2;
    }
};

ostream& operator<<(ostream& os, const PileupAugmenter::NodeOffSide& no);


}

#endif

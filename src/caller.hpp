#ifndef CALLER_H
#define CALLER_H

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <tuple>
#include "vg.pb.h"
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "pileup.hpp"

namespace vg {

using namespace std;

// container for storing pairs of support for calls (value for each strand)
struct StrandSupport {
    int fs; // forward support
    int rs; // reverse support
    int os; // support for other stuff (ie errors)
    double likelihood; // log likelihood from caller (0 if not available)
    StrandSupport(int f = 0, int r = 0, int o = 0, double ll = -1e100) :
        fs(f), rs(r), os(o), likelihood(ll) {}
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
        return fs == other.fs && rs == other.rs && os == other.os && likelihood == other.likelihood;
    }
    // min out at 0
    StrandSupport operator-(const StrandSupport& other) const {
        return StrandSupport(max(0, fs - other.fs), max(0, rs - other.rs),
                             max(0, os - other.os), likelihood);
    }
    StrandSupport& operator+=(const StrandSupport& other) {
        fs += other.fs;
        rs += other.rs;
        os += other.os;
        likelihood = max(likelihood, other.likelihood);
        return *this;
    }
    int depth() { return fs + rs + os; }
    int total() { return fs + rs; }
};

inline ostream& operator<<(ostream& os, const StrandSupport& sup) {
    return os << sup.fs << ", " << sup.rs << ", " << sup.os << ", " << sup.likelihood;
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
    struct Entry { Entry(Node* r = 0, StrandSupport sup_r = StrandSupport(),
                         Node* a1 = 0, StrandSupport sup_a1 = StrandSupport(),
                         Node* a2 = 0, StrandSupport sup_a2 = StrandSupport()) : ref(r), alt1(a1), alt2(a2),
                                                               sup_ref(sup_r), sup_alt1(sup_a1), sup_alt2(sup_a2){}
        Node* ref; Node* alt1; Node* alt2;
        StrandSupport sup_ref; StrandSupport sup_alt1; StrandSupport sup_alt2;
        Node*& operator[](int i) {
            assert(i >= 0 && i <= 2);
            return i == EntryCat::Ref ? ref : (i == EntryCat::Alt1 ? alt1 : alt2);
        }
        StrandSupport& sup(int i) {
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
    void add_fragment(const Node* orig_node, int offset, Node* subnode, EntryCat cat, StrandSupport sup);
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

// Super simple variant caller, for now written to get bakeoff evaluation bootstrapped.
// Idea: Idependently process Pileup records, using simple model to make calls that
//       take into account read errors with diploid assumption.  Edges and node positions
//       are called independently for now.  
// Outputs either a sample graph (only called nodes and edges) or augmented graph
// (include uncalled nodes and edges too).
// The augmented graph (leave_uncalled), with optional text annotation output (text_calls)
// is required to convert calls to VCF. 
class Caller {
public:

    // log of zero
    static const double Log_zero;
    // heterzygous prior (from r MAQ paper)
    static const double Default_het_prior;
    // minimum size of pileup to call a snp
    static const int Default_min_depth;
    // maximum size of pileup to call a snp
    static const int Default_max_depth;
    // minimum number of reads that support snp required to call it
    static const int Default_min_support;
    // same as min_support, but as fraction of total depth
    static const double Default_min_frac;
    // minimum likelihood to call a snp
    static const double Default_min_log_likelihood;
    // use this score when pileup is missing quality
    static const char Default_default_quality;
    // use to balance alignments to forward and reverse strand
    static const double Default_max_strand_bias;
    
    Caller(VG* graph,
           double het_prior = Default_het_prior,
           int min_depth = Default_min_depth,
           int max_depth = Default_max_depth,
           int min_support = Default_min_support,
           double min_frac = Default_min_frac,
           double min_log_likelihood = Default_min_log_likelihood, 
           bool leave_uncalled = false,
           int default_quality = Default_default_quality,
           double max_strand_bias = Default_max_strand_bias,
           ostream* text_calls = NULL,
           bool bridge_alts = false);
    ~Caller();
    void clear();

    // input graph
    VG* _graph;
    // output called graph
    VG _call_graph;
    // optional text file of calls
    ostream* _text_calls;

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

    // used to favour homozygous genotype (r from MAQ paper)
    double _het_log_prior;
    double _hom_log_prior;
    // maximum number of nodes to call before writing out output stream
    int _buffer_size;
    // minimum depth of pileup to call variants on
    int _min_depth;
    // maximum depth of pileup to call variants on
    int _max_depth;
    // min reads supporting snp to call it
    int _min_support;
    // minimum fraction of bases in pileup that nucleotide must have to be snp
    double _min_frac;
    // minimum log likelihood for a snp to be called
    double _min_log_likelihood;
    // don't erase positions that aren't aligned to. so ouput will be original
    // graph plus snps. 
    bool _leave_uncalled;
    // if we don't have a mapping quality for a read position, use this
    char _default_quality;
    // min deviation from .5 in proportion of negative strand reads
    double _max_strand_bias;
    // the base-by-base calling is very limited, and adjacent
    // variants are not properly phased according to the reads.
    // so we choose either to add all edges between neighboring
    // positions (true) or none except via reference (false)
    // (default to latter as most haplotypes rarely contain
    // pairs of consecutive alts). 
    bool _bridge_alts;

    // write the call graph
    void write_call_graph(ostream& out, bool json);

    // call every position in the node pileup
    void call_node_pileup(const NodePileup& pileup);

    // call an edge.  remembering it in a table for the whole graph
    void call_edge_pileup(const EdgePileup& pileup);

    // fill in edges in the call graph (those that are incident to 2 call nodes)
    // and add uncalled nodes (optionally)
    void update_call_graph();

    // map paths from input graph into called graph
    void map_paths();
    // make sure mapped paths generate same strings as input paths
    void verify_path(const Path& in_path, const list<Mapping>& call_path);
    
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
    
    // compute a likelihood from the pileup qualities
    // "first" and "second" are used to virtually split the pileup across two nodes:
    // all bases == "first" are kept
    // all bases == "second" are ignored
    // all otherse are squarerooted (to split their probabilities evenly between the two virtual pileups)
    // returns pair of (likelihood, effective depth), where the effective depth is the number
    // of pileup entries that were considered in computing the likelihood
    pair<double, int> base_log_likelihood(const BasePileup& pb,
                                             const vector<pair<int64_t, int64_t> >& base_offsets,
                                             const string& val, const string& first, const string& second);

    // write graph structure corresponding to all the calls for the current
    // node.  
    void create_node_calls(const NodePileup& np);

    void create_augmented_edge(Node* node1, int from_offset, bool left_side1, bool aug1,
                               Node* node2, int to_offset, bool left_side2, bool aug2, char cat,
                               StrandSupport support);

    // write calling info to tsv to help with VCF conversion
    void write_node_tsv(Node* node, char call, StrandSupport support, int64_t orig_id, int orig_offset);
    void write_edge_tsv(Edge* edge, char call, StrandSupport support);
    void write_nd_tsv();

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

ostream& operator<<(ostream& os, const Caller::NodeOffSide& no);
}

namespace glenn2vcf {
// old glenn2vcf interface
// todo: integrate more smoothly into caller class
int call2vcf(
    // Augmented graph
    vg::VG& vg,
    // "glennfile" as string (relic from old pipeline)
    const string& glennfile,
    // Option variables
    // What's the name of the reference path in the graph?
    string refPathName,
    // What name should we give the contig in the VCF file?
    string contigName,
    // What name should we use for the sample in the VCF file?
    string sampleName,
    // How far should we offset positions of variants?
    int64_t variantOffset,
    // How many nodes should we be willing to look at on our path back to the
    // primary path? Keep in mind we need to look at all valid paths (and all
    // combinations thereof) until we find a valid pair.
    int64_t maxDepth,
    // What should the total sequence length reported in the VCF header be?
    int64_t lengthOverride,
    // Should we load a pileup and print out pileup info as comments after
    // variants?
    string pileupFilename,
    // What fraction of average coverage should be the minimum to call a variant (or a single copy)?
    // Default to 0 because vg call is still applying depth thresholding
    double minFractionForCall,
    // What fraction of the reads supporting an alt are we willing to discount?
    // At 2, if twice the reads support one allele as the other, we'll call
    // homozygous instead of heterozygous. At infinity, every call will be
    // heterozygous if even one read supports each allele.
    double maxHetBias,
    // Like above, but applied to ref / alt ratio (instead of alt / ref)
    double maxRefBias,
    // What's the minimum integer number of reads that must support a call? We
    // don't necessarily want to call a SNP as het because we have a single
    // supporting read, even if there are only 10 reads on the site.
    size_t minTotalSupportForCall,
    // Bin size used for counting coverage along the reference path.  The
    // bin coverage is used for computing the probability of an allele
    // of a certain depth
    size_t refBinSize,
    // On some graphs, we can't get the coverage because it's split over
    // parallel paths.  Allow overriding here
    size_t expCoverage,
    // Should we drop variants that would overlap old ones? TODO: we really need
    // a proper system for accounting for usage of graph material.
    bool suppress_overlaps);
}

#endif

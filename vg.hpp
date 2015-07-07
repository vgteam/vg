#ifndef VG_H
#define VG_H

#include <vector>
#include <set>
#include <string>
#include <deque>
#include <list>
#include <omp.h>
#include <unistd.h>
#include <limits.h>
#include <algorithm>
#include <random>

#include "gssw.h"
#include "gssw_aligner.hpp"
#include "region.hpp"
#include "path.hpp"
#include "utility.hpp"
#include "json.hpp"

#include "vg.pb.h"
#include "hash_map.hpp"

#include "progress_bar.hpp"
#include "lru_cache.h"

#include "Variant.h"
#include "Fasta.h"

#include "swap_remove.hpp"

// uncomment to enable verbose debugging to stderr
//#define debug

namespace vg {


class VG {

public:

    // protobuf-based representation
    // NB: we can't subclass this safely, so it's best as a member
    Graph graph;

    // manages paths of the graph
    // initialized by setting paths._paths = graph.paths
    Paths paths;

    // name
    string name;

    // current id
    int64_t current_id;
    // todo
    //int64_t min_id;
    //int64_t max_id;

    // nodes by id
    hash_map<int64_t, Node*> node_by_id;

    // edges by nodes they connect
    pair_hash_map<pair<int64_t, int64_t>, Edge*> edge_by_id;

    // nodes by position in nodes repeated field
    // this is critical to allow fast deletion of nodes
    hash_map<Node*, int> node_index;

    // edges by position in edges repeated field
    // same as for nodes, this allows fast deletion
    hash_map<Edge*, int> edge_index;

    // edges indexed by nodes they connect
    hash_map<int64_t, vector<int64_t> > edges_from_to;
    hash_map<int64_t, vector<int64_t> > edges_to_from;

    // set the edge indexes through this function
    void set_edge(int64_t from, int64_t to, Edge*);
    void print_edges(void);

    // convenience accessors
    vector<int64_t>& edges_from(Node* node);
    vector<int64_t>& edges_from(int64_t id);
    vector<int64_t>& edges_to(Node* node);
    vector<int64_t>& edges_to(int64_t id);
    void remove_edge_fti(int64_t from, int64_t to);
    void remove_edge_tfi(int64_t from, int64_t to);

    // properties of the graph
    size_t size(void); // number of nodes
    size_t length(void);

    // clear everything
    //void clear(void);

    // constructors

    // default
    VG(void);

    // construct from protobufs
    VG(istream& in, bool showp = false);

    // construct from sets of nodes and edges (e.g. subgraph of another graph)
    VG(set<Node*>& nodes, set<Edge*>& edges);

    // construct from VCF
    VG(vcflib::VariantCallFile& variantCallFile,
       FastaReference& reference,
       string& target,
       int vars_per_region,
       int max_node_size = 0,
       bool showprog = false);
    void from_alleles(const map<long, set<vcflib::VariantAllele> >& altp,
                      string& seq,
                      string& chrom);
    void vcf_records_to_alleles(vector<vcflib::Variant>& records,
                                map<long, set<vcflib::VariantAllele> >& altp,
                                int start_pos,
                                int stop_pos,
                                int max_node_size = 0);
    void slice_alleles(map<long, set<vcflib::VariantAllele> >& altp,
                       int start_pos,
                       int stop_pos,
                       int max_node_size);
    void dice_nodes(int max_node_size);

    void from_gfa(istream& in, bool showp = false);


    // default constructor, destructor
    ~VG(void);

    // copy constructor
    VG(const VG& other) {
        init();
        if (this != &other) {
            // cleanup
            clear_indexes();
            // assign
            graph = other.graph;
            paths = other.paths;
            // re-index
            rebuild_indexes();
        }
    }

    // move constructor
    VG(VG&& other) noexcept {
        init();
        graph = other.graph;
        paths = other.paths;
        other.graph.Clear();
        rebuild_indexes();
        // should copy over indexes
    }

    // copy assignment operator
    VG& operator=(const VG& other) {
        VG tmp(other);
        *this = std::move(tmp);
        return *this;
    }

    // move assignment operator
    VG& operator=(VG&& other) noexcept {
        std::swap(graph, other.graph);
        rebuild_indexes();
        return *this;
    }

    // todo
    // set vg up to not build indexes
    // providing very light runtime
    // this would disable a ton of functions
    // and change the deserialization semantics
    //
    //void no_indexes(void);
    //void yes_indexes(void);

    void build_indexes(void);
    void index_paths(void);
    void clear_indexes(void);
    void clear_indexes_no_resize(void);
    void resize_indexes(void);
    void rebuild_indexes(void);

    // literally merge protobufs
    void merge(Graph& g);
    void merge(VG& g);

    // merge protobufs after removing overlaps
    // good when there aren't many overlaps
    void merge_union(VG& g);
    // helper to merge_union
    void remove_duplicated_in(VG& g);

    // limit the local complexity of the graph, connecting pruned components to a head and tail node
    // depending on the direction which we come into the node when the edge_max is passed
    void prune_complex_paths(int length, int edge_max, Node* head_node, Node* tail_node);
    void prune_short_subgraphs(size_t min_size);

    // write to a stream in chunked graphs
    void serialize_to_ostream(ostream& out, int64_t chunk_size = 1000);
    void serialize_to_file(const string& file_name, int64_t chunk_size = 1000);

    // can we handle this with merge?
    //void concatenate(VG& g);

    int64_t max_node_id(void);
    int64_t min_node_id(void);
    void compact_ids(void);
    void increment_node_ids(int64_t increment);
    void decrement_node_ids(int64_t decrement);
    void swap_node_id(int64_t node_id, int64_t new_id);
    void swap_node_id(Node* node, int64_t new_id);

    // iteratively add when nodes and edges are novel
    // good when there are very many overlaps
    void extend(VG& g);
    void extend(Graph& graph);

    // modify ids of the second graph to ensure we don't have conflicts
    // then attach tails of this graph to the heads of the other, and extend(g)
    void append(VG& g);

    // don't append or join the nodes in the graphs
    // just ensure that ids are unique, then apply extend
    void combine(VG& g);

    // edit the graph to include the path
    void include(const Path& path);
    // or a set of mappings against one node
    void edit_node(int64_t node_id, const vector<Mapping>& mappings);
    // for each node, modify it with the associated mappings
    void edit(const map<int64_t, vector<Mapping> >& mappings);

    void add_node(Node& node);
    void add_nodes(vector<Node>& nodes);
    void add_edge(Edge& edge);
    void add_edges(vector<Edge>& edges);
    void add_nodes(set<Node*>& nodes);
    void add_edges(set<Edge*>& edges);

    int64_t node_count(void);
    int64_t edge_count(void);
    int64_t total_length_of_nodes(void);
    int in_degree(Node* node);
    int out_degree(Node* node);
    void edges_of_node(Node* node, vector<Edge*>& edges);
    void edges_of_nodes(set<Node*>& nodes, set<Edge*>& edges);

    // use the VG class to generate ids
    Node* create_node(string seq);
    Node* get_node(int64_t id);
    void node_context(Node* node, VG& g);
    void destroy_node(Node* node);
    void destroy_node(int64_t id);
    bool has_node(int64_t id);
    bool has_node(Node* node);
    bool has_node(Node& node);
    void for_each_node(function<void(Node*)> lambda);
    void for_each_node_parallel(function<void(Node*)> lambda);

    // is the graph empty?
    bool empty(void);

    // remove nodes with no sequence
    // these are created in some cases during the process of graph construction
    void remove_null_nodes(void);
    // remove a node but connect all of its predecessor and successor nodes with new edges 
    void remove_node_forwarding_edges(Node* node);
    // remove null nodes but connect predecessors and successors, preserving structure
    void remove_null_nodes_forwarding_edges(void);

    // remove edges for which one of the nodes is not present
    void remove_orphan_edges(void);
    void keep_paths(set<string>& path_names, set<string>& kept_names);
    void keep_path(string& path_name);

    // path stats
    int path_edge_count(list<Node*>& path, int32_t offset, int path_length);
    int path_end_node_offset(list<Node*>& path, int32_t offset, int path_length);

    // edges
    Edge* create_edge(Node* from, Node* to);
    Edge* create_edge(int64_t from, int64_t to);
    Edge* get_edge(int64_t from, int64_t to);
    void destroy_edge(Edge* edge);
    void destroy_edge(int64_t from, int64_t to);
    bool has_edge(int64_t from, int64_t to);
    bool has_edge(Edge* edge);
    bool has_edge(Edge& edge);
    void for_each_edge(function<void(Edge*)> lambda);
    void for_each_edge_parallel(function<void(Edge*)> lambda);

    // connect node -> nodes
    void connect_node_to_nodes(Node* node, vector<Node*>& nodes);
    // connect nodes -> node
    void connect_nodes_to_node(vector<Node*>& nodes, Node* node);

    // utilities
    void divide_node(Node* node, int pos, Node*& left, Node*& right);
    void divide_path(map<long, int64_t>& path, long pos, Node*& left, Node*& right);
    //void node_replace_prev(Node* node, Node* before, Node* after);
    //void node_replace_next(Node* node, Node* before, Node* after);

    void to_dot(ostream& out, vector<Alignment> alignments = {});
    void to_gfa(ostream& out);
    bool is_valid(void);

    // topologically orders nodes
    void sort(void);
    // helper function, not really meant for external use
    void topological_sort(deque<Node*>& l);
    void swap_nodes(Node* a, Node* b);

    Alignment& align(Alignment& alignment);
    Alignment align(string& sequence);
    //Alignment& align(Alignment& alignment);
    void destroy_alignable_graph(void);

    GSSWAligner* gssw_aligner;

    // returns all node-crossing paths with up to length across node boundaries
    void for_each_kpath(int k, int edge_max,
                        function<void(Node*)> handle_prev_maxed,
                        function<void(Node*)> handle_next_maxed,
                        function<void(Node*,list<Node*>&)> lambda);
    void for_each_kpath_parallel(int k, int edge_max,
                                 function<void(Node*)> handle_prev_maxed,
                                 function<void(Node*)> handle_next_maxed,
                                 function<void(Node*,list<Node*>&)> lambda);
    void for_each_kpath(int k, int edge_max,
                        function<void(Node*)> handle_prev_maxed,
                        function<void(Node*)> handle_next_maxed,
                        function<void(Node*,Path&)> lambda);
    void for_each_kpath_parallel(int k, int edge_max,
                                 function<void(Node*)> handle_prev_maxed,
                                 function<void(Node*)> handle_next_maxed,
                                 function<void(Node*,Path&)> lambda);
    void for_each_kpath_of_node(Node* node, int k, int edge_max,
                                function<void(Node*)> handle_prev_maxed,
                                function<void(Node*)> handle_next_maxed,
                                function<void(Node*,list<Node*>&)> lambda);
    void for_each_kpath_of_node(Node* n, int k, int edge_max,
                                function<void(Node*)> handle_prev_maxed,
                                function<void(Node*)> handle_next_maxed,
                                function<void(Node*,Path&)> lambda);

    void kpaths(set<list<Node*> >& paths, int length, int edge_max,
                function<void(Node*)> prev_maxed, function<void(Node*)> next_maxed);
    void kpaths(vector<Path>& paths, int length, int edge_max,
                function<void(Node*)> prev_maxed, function<void(Node*)> next_maxed);

    void kpaths_of_node(Node* node, set<list<Node*> >& paths,
                        int length, int edge_max,
                        function<void(Node*)> prev_maxed, function<void(Node*)> next_maxed);
    void kpaths_of_node(Node* node, vector<Path>& paths,
                        int length, int edge_max,
                        function<void(Node*)> prev_maxed, function<void(Node*)> next_maxed);
    void kpaths_of_node(int64_t node_id, vector<Path>& paths, int length, int edge_max,
                        function<void(Node*)> prev_maxed, function<void(Node*)> next_maxed);
    void prev_kpaths_from_node(Node* node, int length, int edge_max,
                               list<Node*> postfix, set<list<Node*> >& paths,
                               function<void(Node*)>& maxed_nodes);
    void next_kpaths_from_node(Node* node, int length, int edge_max,
                               list<Node*> prefix, set<list<Node*> >& paths,
                               function<void(Node*)>& maxed_nodes);

    void paths_between(Node* from, Node* to, vector<Path>& paths);
    void paths_between(int64_t from, int64_t to, vector<Path>& paths);
    void likelihoods(vector<Alignment>& alignments, vector<Path>& paths, vector<long double>& likelihoods);

    string path_sequence(const Path& path);

    // traversal
    void nodes_prev(Node* n, vector<Node*>& nodes);
    void nodes_next(Node* n, vector<Node*>& nodes);
    int node_count_prev(Node* n);
    int node_count_next(Node* n);

    // paths
    Path create_path(const list<Node*>& nodes);
    Path create_path(const vector<Node*>& nodes);
    string path_string(const list<Node*>& nodes);
    string path_string(Path& path);
    void expand_path(const list<Node*>& path, vector<Node*>& expanded);
    void node_starts_in_path(const list<Node*>& path,
                             map<Node*, int>& node_start);

    // kmers
    void for_each_kmer_parallel(int kmer_size,
                                int edge_max,
                                function<void(string&, Node*, int, list<Node*>&, VG&)> lambda,
                                int stride = 1,
                                bool allow_dups = false,
                                bool allow_negatives = false);
    void for_each_kmer(int kmer_size,
                       int edge_max,
                       function<void(string&, Node*, int, list<Node*>&, VG&)> lambda,
                       int stride = 1,
                       bool allow_dups = false,
                       bool allow_negatives = false);

    // for gcsa2
    // For the given kmer of the given length starting at the given offset into
    // the given Node along the given path, fill in prev_chars with the
    // characters that preceed it, next_chars with the characters that follow
    // it, and next_positions with the (node ID, offset) pairs of the places you
    // can go next (from the right end of the kmer). Refuses to follow more than
    // edge_max edges
    void kmer_context(string& kmer,
                      int kmer_size,
                      int edge_max,
                      list<Node*>& path,
                      Node* node,
                      int32_t offset,
                      set<char>& prev_chars,
                      set<char>& next_chars,
                      set<pair<int64_t, int32_t> >& next_positions);
    // for pruning graph prior to indexing with gcsa2
    // takes all nodes that would introduce paths of > edge_max edge crossings, removes them, and links their neighbors to
    // head_node or tail_node depending on which direction the path extension was stopped
    void prune_complex(int path_length, int edge_max, Node* head_node, Node* tail_node);

private:
    void _for_each_kmer(int kmer_size,
                        int edge_max,
                        function<void(string&, Node*, int, list<Node*>&, VG&)> lambda,
                        bool parallel,
                        int stride,
                        bool allow_dups,
                        bool allow_negatives);


public:

    // reads
    string random_read(int length, mt19937& rng, int64_t min_id, int64_t max_id, bool either_strand);

    // subgraphs
    void disjoint_subgraphs(list<VG>& subgraphs);
    void head_nodes(vector<Node*>& nodes);
    vector<Node*> head_nodes(void);
    void tail_nodes(vector<Node*>& nodes);
    vector<Node*> tail_nodes(void);
    void collect_subgraph(Node* node, set<Node*>& subgraph);

    // join head nodes of graph to common null node
    Node* join_heads(void);
    // or heads and tails to common
    void join_heads(Node* node);
    void join_tails(Node* node);

    // add singular head and tail null nodes to graph
    void wrap_with_null_nodes(void);
    // to prepare for indexing with GCSA; length is at least kmer length to be used
    void add_start_and_end_markers(int length, char start_char, char end_char,
                                   Node*& head_node, Node*& tail_node);

    bool show_progress;
    string progress_message;
    long progress_count;
    long last_progress;
    ProgressBar* progress;
    void create_progress(const string& message, long count);
    void create_progress(long count);
    void update_progress(long i);
    void destroy_progress(void);

    // for managing parallel construction
    struct Plan {
        VG* graph;
        map<long, set<vcflib::VariantAllele> >* alleles;
        string seq;
        string name;
        Plan(VG* g,
             map<long, set<vcflib::VariantAllele> >* a,
             string s,
             string n)
            : graph(g)
            , alleles(a)
            , seq(s)
            , name(n) { };
        ~Plan(void) { delete alleles; }
    };


private:

    void init(void); // setup, ensures that gssw == NULL on startup
    // placeholder for empty
    vector<int64_t> empty_ids;

};

// utility functions

bool allATGC(string& s);
void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar);
string cigar_string(vector<pair<int, char> >& cigar);
string mapping_string(const string& source, const Mapping& mapping);
void divide_invariant_mapping(Mapping& orig, Mapping& left, Mapping& right, int offset, Node* nl, Node* nr);

} // end namespace vg

#endif

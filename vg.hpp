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

#include "gssw.h"
#include "gssw_aligner.hpp"
#include "region.hpp"

#include "vg.pb.h"
#include "stream.hpp"
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

    // name
    string name;

    // current id
    int64_t current_id;

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

    // constructors

    // default
    VG(void);

    // construct from protobufs
    VG(istream& in);

    // construct from sets of nodes and edges (e.g. subgraph of another graph)
    VG(set<Node*>& nodes, set<Edge*>& edges);

    // construct from VCF
    VG(vcf::VariantCallFile& variantCallFile,
       FastaReference& reference,
       string& target,
       int vars_per_region,
       int max_node_size = 0,
       bool showprog = false);
    void from_alleles(const map<long, set<vcf::VariantAllele> >& altp,
                      string& seq,
                      string& chrom);
    void vcf_records_to_alleles(vector<vcf::Variant>& records,
                                map<long, set<vcf::VariantAllele> >& altp,
                                int start_pos,
                                int stop_pos,
                                int max_node_size = 0);


    // default constructor, destructor
    ~VG(void);
    VG& operator=(const VG& other) {
        if (this != &other) {
            // cleanup
            clear_indexes();
            // assign
            graph = other.graph;
            // re-index
            build_indexes();
        }
        return *this;
    }

    void build_indexes(void);
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

    // write to a stream in chunked graphs
    void serialize_to_ostream(ostream& out, int64_t chunk_size = 1000);

    // can we handle this with merge?
    //void concatenate(VG& g);

    int64_t max_node_id(void);
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

    // utilities
    void divide_node(Node* node, int pos, Node*& left, Node*& right);
    void divide_path(map<long, Node*>& path, long pos, Node*& left, Node*& right);
    //void node_replace_prev(Node* node, Node* before, Node* after);
    //void node_replace_next(Node* node, Node* before, Node* after);

    void to_dot(ostream& out);
    void to_gfa(ostream& out);
    bool is_valid(void);

    void topologically_sort_graph(void);
    void topological_sort(deque<Node*>& l);
    void swap_nodes(Node* a, Node* b);

    Alignment& align(Alignment& alignment);
    Alignment align(string& sequence);
    //Alignment& align(Alignment& alignment);
    void destroy_alignable_graph(void);

    GSSWAligner* gssw_aligner;

    // returns all node-crossing paths with up to length across node boundaries
    void for_each_kpath(int k, function<void(list<Node*>&)> lambda);
    void for_each_kpath_parallel(int k, function<void(list<Node*>&)> lambda);
    void for_each_kpath(int k, function<void(Path&)> lambda);
    void for_each_kpath_parallel(int k, function<void(Path&)> lambda);

    void for_each_kpath_of_node(Node* n, int k, function<void(list<Node*>&)> lambda);
    void for_each_kpath_of_node(Node* n, int k, function<void(Path&)> lambda);

    void kpaths(vector<Path>& paths, int length);
    void kpaths(set<list<Node*> >& paths, int length);

    void kpaths_of_node(Node* node, set<list<Node*> >& paths, int length);
    void kpaths_of_node(Node* node, vector<Path>& paths, int length);
    void kpaths_of_node(int64_t node_id, vector<Path>& paths, int length);
    void prev_kpaths_from_node(Node* node, int length, list<Node*> postfix, set<list<Node*> >& paths);
    void next_kpaths_from_node(Node* node, int length, list<Node*> prefix, set<list<Node*> >& paths);

    void paths_between(Node* from, Node* to, vector<Path>& paths);
    void paths_between(int64_t from, int64_t to, vector<Path>& paths);
    void likelihoods(vector<Alignment>& alignments, vector<Path>& paths, vector<long double>& likelihoods);

    string path_sequence(Path& path);

    // traversal
    void nodes_prev(Node* n, vector<Node*>& nodes);
    void nodes_next(Node* n, vector<Node*>& nodes);

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
                                function<void(string&, Node*, int)> lambda,
                                int stride = 1);
    void for_each_kmer(int kmer_size,
                       function<void(string&, Node*, int)> lambda,
                       int stride = 1);
    
private:
    void _for_each_kmer(int kmer_size,
                        function<void(string&, Node*, int)> lambda,
                        bool parallel,
                        int stride);

public:

    // subgraphs
    void disjoint_subgraphs(list<VG>& subgraphs);
    void head_nodes(vector<Node*>& nodes);
    vector<Node*> head_nodes(void);
    void tail_nodes(vector<Node*>& nodes);
    vector<Node*> tail_nodes(void);
    void collect_subgraph(Node* node, set<Node*>& subgraph);

    // join head nodes of graph to common null node
    Node* join_heads(void);

    // add singular head and tail null nodes to graph
    void wrap_with_null_nodes(void);

    bool show_progress;
    string progress_message;
    long progress_count;
    long last_progress;
    ProgressBar* progress;
    void create_progress(const string& message, long count);
    void create_progress(long count);
    void update_progress(long i);
    void destroy_progress(void);

private:

    void init(void); // setup, ensures that gssw == NULL on startup
    // placeholder for empty
    vector<int64_t> empty_ids;

};

// utility functions

bool allATGC(string& s);

} // end namespace vg

#endif

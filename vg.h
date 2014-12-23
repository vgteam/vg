#ifndef VG_H
#define VG_H

#include <vector>
#include <set>
#include <string>
#include <deque>
#include <list>
#include <omp.h>
#include <unistd.h>

#include "gssw.h"
#include "gssw_aligner.h"
#include "region.h"

#include "vg.pb.h"
#include "stream.h"

extern "C" {
#include "progressbar.h"
}

#include "Variant.h"
#include "Fasta.h"


namespace vg {

class VG {

public:

    // protobuf-based representation
    // NB: we can't subclass this safely, so it's best as a member
    Graph graph;

    // current id
    int64_t current_id;

    // nodes by id
    map<int64_t, Node*> node_by_id;

    // nodes by position in nodes repeated field
    // this is critical to allow fast deletion of nodes
    map<Node*, int> node_index;

    // edges indexed by nodes they connect
    map<int64_t, map<int64_t, Edge*> > edge_from_to;
    map<int64_t, map<int64_t, Edge*> > edge_to_from;

    // convenience accessors
    map<int64_t, Edge*>& edges_from(Node* node);
    map<int64_t, Edge*>& edges_from(int64_t id);
    map<int64_t, Edge*>& edges_to(Node* node);
    map<int64_t, Edge*>& edges_to(int64_t id);

    // edges by position in edges repeated field
    // same as for nodes, this allows fast deletion
    map<Edge*, int> edge_index;

    // constructors

    // default
    VG(void);

    // construct from protobufs
    VG(istream& in);

    // construct from sets of nodes and edges (e.g. subgraph of another graph)
    VG(set<Node*>& nodes, set<Edge*>& edges);

    // construct from VCF records
    VG(vcf::VariantCallFile& variantCallFile, FastaReference& reference, string& target, int vars_per_region);
    VG(vector<vcf::Variant>& records, string seq, string chrom, int offset);
    void from_vcf_records(vector<vcf::Variant>* records, string seq, string chrom, int offset);

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
    void rebuild_indexes(void);

    // literally merge protobufs
    void merge(Graph& g);

    // merge protobufs after removing overlaps
    // good when there aren't many overlaps
    void merge(VG& g);
    // helper to merge
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
    void bounded_paths(vector<Path>& paths, int length);
    void bounded_paths(set<list<Node*> >& paths, int length);
    void bounded_paths(Node* node, set<list<Node*> >& paths, int length);
    void bounded_paths(Node* node, vector<Path>& paths, int length);
    void bounded_paths(int64_t node_id, vector<Path>& paths, int length);
    void bounded_prev_paths_from_node(Node* node, int length, list<Node*> postfix, set<list<Node*> >& paths);
    void bounded_next_paths_from_node(Node* node, int length, list<Node*> prefix, set<list<Node*> >& paths);

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
    void kmers_of(map<string, map<Node*, int> >& kmer_map, int kmer_size);

    // subgraphs
    void disjoint_subgraphs(list<VG>& subgraphs);
    void head_nodes(vector<Node*>& nodes);
    void tail_nodes(vector<Node*>& nodes);
    void collect_subgraph(Node* node, set<Node*>& subgraph);

    // join head nodes of graph to common null node
    Node* join_heads(void);

    bool show_progress;

private:

    void init(void); // setup, ensures that gssw == NULL on startup

};


} // end namespace vg

#endif

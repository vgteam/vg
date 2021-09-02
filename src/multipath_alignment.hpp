/// \file multipath_alignment.hpp
///
/// utility functions for the multipath_alignment_t object
///

#ifndef multipath_alignment_hpp
#define multipath_alignment_hpp

#include <stdio.h>
#include <vector>
#include <list>
#include <unordered_map>
#include <algorithm>

#include <vg/vg.pb.h>
#include "path.hpp"
#include "alignment.hpp"
#include "utility.hpp"
#include "handle.hpp"
#include "annotation.hpp"

// Declare the haplo::ScoreProvider we use for haplotype-aware traceback generation.
namespace haplo {
    class ScoreProvider;
}

namespace vg {

    class connection_t {
    public:
        connection_t() = default;
        connection_t(const connection_t&) = default;
        connection_t(connection_t&&) = default;
        ~connection_t() = default;
        connection_t& operator=(const connection_t&) = default;
        connection_t& operator=(connection_t&&) = default;
        inline int32_t next() const;
        inline void set_next(int32_t n);
        inline int32_t score() const;
        inline void set_score(int32_t s);
    private:
        uint32_t _next;
        int32_t _score;
    };

    /*
     * STL implementations of the protobuf object for use in in-memory operations
     */
    class subpath_t {
    public:
        subpath_t() = default;
        subpath_t(const subpath_t&) = default;
        subpath_t(subpath_t&&) = default;
        ~subpath_t() = default;
        subpath_t& operator=(const subpath_t&) = default;
        subpath_t& operator=(subpath_t&&) = default;
        inline const path_t& path() const;
        inline path_t* mutable_path();
        inline bool has_path() const;
        inline const vector<uint32_t>& next() const;
        inline uint32_t next(size_t i) const;
        inline vector<uint32_t>* mutable_next();
        inline void set_next(size_t i, uint32_t n);
        inline void add_next(uint32_t n);
        inline void clear_next();
        inline size_t next_size() const;
        inline bool has_next() const;
        inline int32_t score() const;
        inline void set_score(int32_t s);
        inline const vector<connection_t>& connection() const;
        inline const connection_t& connection(size_t i) const;
        inline vector<connection_t>* mutable_connection();
        inline connection_t* mutable_connection(size_t i);
        inline void set_connection(size_t i, const connection_t& c);
        inline connection_t* add_connection();
        inline void clear_connection();
        inline size_t connection_size() const;
        inline bool has_connection() const;
    private:
        path_t _path;
        vector<uint32_t> _next;
        int32_t _score;
        vector<connection_t> _connection;
    };

    // TODO: the metadata could be removed and only added to the protobuf at serialization time
    class multipath_alignment_t {
    public:
        multipath_alignment_t();
        multipath_alignment_t(const multipath_alignment_t& other);
        multipath_alignment_t(multipath_alignment_t&& other);
        ~multipath_alignment_t();
        multipath_alignment_t& operator=(const multipath_alignment_t& other);
        multipath_alignment_t& operator=(multipath_alignment_t&& other);
        inline const string& sequence() const;
        inline string* mutable_sequence();
        inline void set_sequence(const string& s);
        inline const string& quality() const;
        inline string* mutable_quality();
        inline void set_quality(const string& q);
        inline const vector<subpath_t>& subpath() const;
        inline const subpath_t& subpath(size_t i) const;
        inline vector<subpath_t>* mutable_subpath();
        inline subpath_t* mutable_subpath(size_t i);
        inline subpath_t* add_subpath();
        inline void clear_subpath();
        inline size_t subpath_size() const;
        inline int32_t mapping_quality() const;
        inline void set_mapping_quality(int32_t q);
        inline const vector<uint32_t>& start() const;
        inline uint32_t start(size_t i) const;
        inline vector<uint32_t>* mutable_start();
        inline void set_start(size_t i, uint32_t s);
        inline void add_start(uint32_t s);
        inline void clear_start();
        inline size_t start_size() const;
        inline bool has_start() const;
        
        // annotation interface
        // TODO: add List and Struct from https://github.com/protocolbuffers/protobuf/blob/master/src/google/protobuf/struct.proto
        enum anno_type_t {Null = 0, Double = 2, Bool = 3, String = 4};
        void set_annotation(const string& annotation_name);
        void set_annotation(const string& annotation_name, double value);
        void set_annotation(const string& annotation_name, bool value);
        void set_annotation(const string& annotation_name, const string& value);
        void clear_annotation(const string& annotation_name);
        pair<anno_type_t, const void*> get_annotation(const string& annotation_name) const;
        void for_each_annotation(function<void(const string&, anno_type_t, const void*)> lambda) const;
    private:
        string _sequence;
        string _quality;
        vector<subpath_t> _subpath;
        int32_t _mapping_quality;
        vector<uint32_t> _start;
        map<string, pair<anno_type_t, void*>> _annotation;
    };

    string debug_string(const connection_t& connection);
    string debug_string(const subpath_t& subpath);
    string debug_string(const multipath_alignment_t& multipath_aln);
    
    /// Put subpaths in topological order (assumed to be true for other algorithms)
    void topologically_order_subpaths(multipath_alignment_t& multipath_aln);
    
    /// Finds the start subpaths (i.e. the source nodes of the multipath DAG) and stores
    /// them in the 'start' field of the multipath_alignment_t
    void identify_start_subpaths(multipath_alignment_t& multipath_aln);
    
    /// Stores the highest scoring alignment contained in the multipath_alignment_t in an Alignment
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal path through
    ///    aln_out           empty alignment to store optimal alignment in (data will be
    ///                      overwritten if not empty)
    ///    subpath_global    if true, only allows alignments that source subpath to sink subpath
    ///                      in the multipath DAG, else allows any start and end subpath
    ///
    void optimal_alignment(const multipath_alignment_t& multipath_aln, Alignment& aln_out,
                           bool subpath_global = false);
    
    /// Returns the score of the highest scoring alignment contained in the multipath_alignment_t
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal score in
    ///    subpath_global    if true, only allows alignments that source subpath to sink subpath
    ///                      in the multipath DAG, else allows any start and end subpath
    ///
    int32_t optimal_alignment_score(const multipath_alignment_t& multipath_aln,
                                    bool subpath_global = false);

    /// Returns the score of the lowest-scoring source-to-sink alignment in the multipath_alignment_t.
    /// Assumes that subpaths are topologically ordered and starts have been identified.
    int32_t worst_alignment_score(const multipath_alignment_t& multipath_aln);
    
    /// Returns the top k highest-scoring alignments contained in the multipath_alignment_t.
    /// Note that some or all of these may be duplicate Alignments, which were spelled out
    /// by tracebacks through different sequences of subpaths that shared alignment material.
    ///
    /// If the best alignment is no alignment (i.e. the read is unmapped), returns an empty vector.
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal paths through
    ///    count             maximum number of top alignments to return
    ///
    vector<Alignment> optimal_alignments(const multipath_alignment_t& multipath_aln, size_t count);
    
    /// Finds k or fewer top-scoring alignments using only distinct subpaths.
    /// Asymmetrical: the optimal alignment for each end subpath is found, greedily, subject to the constraint,
    /// but the other subpaths are first-come first-serve. Also, distinct subpaths may not guarantee distinct
    /// actual alignments, so alignments may need deduplication.
    ///
    /// If the best alignment is no alignment (i.e. the read is unmapped), returns an empty vector.
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal paths through
    ///    count             maximum number of top alignments to return
    ///
    vector<Alignment> optimal_alignments_with_disjoint_subpaths(const multipath_alignment_t& multipath_aln, size_t count);
    

    /// Finds all alignments consistent with haplotypes available by incremental search with the given haplotype
    /// score provider. Pads to a certain count with haplotype-inconsistent alignments that are population-scorable
    /// (i.e. use only edges used by some haplotype in the index), and then with unscorable alignments if scorable
    /// ones are unavailable. This may result in an empty vector.
    ///
    /// Output Alignments may not be unique. The input multipath_alignment_t may have exponentially many ways to
    /// spell the same Alignment, and we will look at all of them. We also may have duplicates of the optimal
    /// alignment if we are asked to produce it unconsitionally.
    ///
    /// Note: Assumes that each subpath's Path object uses one Mapping per node and that
    /// start subpaths have been identified
    ///
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to find optimal paths through
    ///    score_provider    a haplo::ScoreProvider that supports incremental search over its haplotype database (such as a GBWTScoreProvider)
    ///    soft_count        maximum number of haplotype-inconsistent alignments to pad to
    ///    hard_count        maximum number of alignments, including haplotype-consistent (0 if no limit)
    ///    optimal_first     always compute and return first the optimal alignment, even if not haplotype-consistent
    ///
    vector<Alignment> haplotype_consistent_alignments(const multipath_alignment_t& multipath_aln, const haplo::ScoreProvider& score_provider,
        size_t soft_count, size_t hard_count, bool optimal_first = false);
    
    /// The indexes on the read sequence of the portion of the read that is aligned outside of soft clips
    pair<int64_t, int64_t> aligned_interval(const multipath_alignment_t& multipath_aln);

    /// Stores the reverse complement of a multipath_alignment_t in another multipath_alignment_t
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to reverse complement
    ///    node_length       a function that returns the length of a node sequence from its node ID
    ///    rev_comp_out      empty multipath alignment to store reverse complement in (some data may
    ///                      be overwritten if not empty)
    ///
    void rev_comp_multipath_alignment(const multipath_alignment_t& multipath_aln,
                                      const function<int64_t(int64_t)>& node_length,
                                      multipath_alignment_t& rev_comp_out);
    
    /// Stores the reverse complement of a multipath_alignment_t in another multipath_alignment_t
    ///
    ///  Args:
    ///    multipath_aln     multipath alignment to reverse complement in place
    ///    node_length       a function that returns the length of a node sequence from its node ID
    ///
    void rev_comp_multipath_alignment_in_place(multipath_alignment_t* multipath_aln,
                                               const function<int64_t(int64_t)>& node_length);

    /// Replaces all U's in the sequence and the aligned Paths with T's
    void convert_Us_to_Ts(multipath_alignment_t& multipath_aln);

    /// Replaces all T's in the sequence and the aligned Paths with U's
    void convert_Ts_to_Us(multipath_alignment_t& multipath_aln);
    
    /// Convert an STL-based multipath_alignment_t to a protobuf MultipathAlignment
    void to_proto_multipath_alignment(const multipath_alignment_t& multipath_aln,
                                      MultipathAlignment& proto_multipath_aln_out);

    /// Convert a protobuf MultipathAlignment to an STL-based multipath_alignment_t
    void from_proto_multipath_alignment(const MultipathAlignment& proto_multipath_aln,
                                        multipath_alignment_t& multipath_aln_out);

    /// Converts a Alignment into a multipath_alignment_t  with one subpath and stores it in an object
    ///
    ///  Args:
    ///    aln               alignment to convert
    ///    multipath_aln     empty multipath alignment to store converted alignment in (data may be
    ///                      be overwritten if not empty)
    ///
    void to_multipath_alignment(const Alignment& aln, multipath_alignment_t& multipath_aln_out);
    
    // TODO: these metadata functions should also transfer annotations

    /// Copies metadata from an Alignment object and transfers it to a multipath_alignment_t
    ///
    ///  Args:
    ///    from    copy metadata from this
    ///    to      into this
    ///
    void transfer_read_metadata(const Alignment& from, multipath_alignment_t& to);
    
    /// Copies metadata from an multipath_alignment_t object and transfers it to a Alignment
    ///
    ///  Args:
    ///    from    copy metadata from this
    ///    to      into this
    ///
    void transfer_read_metadata(const multipath_alignment_t& from, Alignment& to);
    
    /// Copies metadata from an multipath_alignment_t object and transfers it to another multipath_alignment_t
    ///
    ///  Args:
    ///    from    copy metadata from this
    ///    to      into this
    ///
    void transfer_read_metadata(const multipath_alignment_t& from, multipath_alignment_t& to);
    
    /// Copies metadata from an Alignment object and transfers it to another Alignment
    ///
    ///  Args:
    ///    from    copy metadata from this
    ///    to      into this
    ///
    void transfer_read_metadata(const Alignment& from, Alignment& to);

    void transfer_read_metadata(const MultipathAlignment& from, multipath_alignment_t& to);

    void transfer_read_metadata(const multipath_alignment_t& from, MultipathAlignment& to);

    void transfer_proto_metadata(const Alignment& from, MultipathAlignment& to);

    void transfer_proto_metadata(const MultipathAlignment& from, Alignment& to);

    /// Merges non-branching paths in a multipath alignment in place
    /// Does not assume topological order among subpaths
    void merge_non_branching_subpaths(multipath_alignment_t& multipath_aln,
                                      const unordered_set<size_t>* prohibited_merges = nullptr);

    /// Removes all edit, mappings, and subpaths that have no aligned bases, and introduces transitive edges
    /// to preserve connectivity through any completely removed subpaths
    void remove_empty_alignment_sections(multipath_alignment_t& multipath_aln);
    
    /// Returns a vector whose elements are vectors with the indexes of the Subpaths in
    /// each connected component. An unmapped multipath_alignment_t with no subpaths produces an empty vector.
    vector<vector<int64_t>> connected_components(const multipath_alignment_t& multipath_aln);
    
    /// Extract the multipath_alignment_t consisting of the Subpaths with the given indexes
    /// into a new multipath_alignment_t object
    void extract_sub_multipath_alignment(const multipath_alignment_t& multipath_aln,
                                         const vector<int64_t>& subpath_indexes,
                                         multipath_alignment_t& sub_multipath_aln);

    /// Add the subpaths of one multipath alignment onto another
    void append_multipath_alignment(multipath_alignment_t& multipath_aln,
                                    const multipath_alignment_t& to_append);

    /// Returns true if any subpath has a connection adjacency
    bool contains_connection(const multipath_alignment_t& multipath_aln);

    /// Returns all of the positions where a given sequence index occurs at a given graph
    /// graph position (if any), where positions are represented as tuples of
    /// (subpath index, mapping index, edit index, index within edit)
    vector<tuple<int64_t, int64_t, int64_t, int64_t>>
    search_multipath_alignment(const multipath_alignment_t& multipath_aln,
                               const pos_t& graph_pos, int64_t seq_pos);

    /// Returns a pair of (mapping, edit, base) and possibly multiple (subpath, mapping, edit, base),of the furthest position
    /// that can be traced through the multipath alignment along the pathstarting the indicated position in the multipath
    /// alignment. The path can be traced rightward starting at the beginning, or leftward starting.
    /// Search is limited to not passing a given mapping on the path.
    pair<tuple<int64_t, int64_t, int64_t>, vector<tuple<int64_t, int64_t, int64_t, int64_t>>>
    trace_path(const multipath_alignment_t& multipath_aln, const Path& path,
               int64_t subpath_idx, int64_t mapping_idx, int64_t edit_idx, int64_t base_idx, bool search_left,
               int64_t search_limit);

    /// Returns true if the multipath alignment contains a match of a given length starting at the graph and
    /// read position
    bool contains_match(const multipath_alignment_t& multipath_aln, const pos_t& pos,
                        int64_t read_pos, int64_t match_length);

    /// Convert a surjected multipath alignment into a CIGAR sequence against a path. Splicing will be allowed
    /// at connections and at any silent deletions of path sequence. Surjected multipath alignment graph must
    /// consist of a single non-branching path
    vector<pair<int, char>> cigar_against_path(const multipath_alignment_t& multipath_aln, const string& path_name, bool rev,
                                               int64_t path_pos, const PathPositionHandleGraph& graph,
                                               int64_t min_splice_length = numeric_limits<int64_t>::max());

    /// Debugging function to check that multipath alignment meets the formalism's basic
    /// invariants. Returns true if multipath alignment is valid, else false. Does not
    /// validate alignment score.
    bool validate_multipath_alignment(const multipath_alignment_t& multipath_aln, const HandleGraph& handle_graph);
    
    /// Send a formatted string representation of the multipath_alignment_t into the ostream
    void view_multipath_alignment(ostream& out, const multipath_alignment_t& multipath_aln, const HandleGraph& handle_graph);
    
    /// Converts a multipath_alignment_t to a GraphViz Dot representation, output to the given ostream.
    void view_multipath_alignment_as_dot(ostream& out, const multipath_alignment_t& multipath_aln, bool show_graph = false);
    
    // TODO: function for adding a graph augmentation to an existing multipath alignment

    /*
     * Implementations of inline methods
     */

    /*
     * connection_t
     */
    inline int32_t connection_t::next() const {
        return _next;
    }
    inline void connection_t::set_next(int32_t n) {
        _next = n;
    }
    inline int32_t connection_t::score() const {
        return _score;
    }
    inline void connection_t::set_score(int32_t s) {
        _score = s;
    }

    /*
     * subpath_t
     */
    inline const path_t& subpath_t::path() const {
        return _path;
    }
    inline path_t* subpath_t::mutable_path() {
        return &_path;
    }
    inline bool subpath_t::has_path() const {
        return _path.mapping_size();
    }
    inline const vector<uint32_t>& subpath_t::next() const {
        return _next;
    }
    inline uint32_t subpath_t::next(size_t i) const {
        return _next[i];
    }
    inline vector<uint32_t>* subpath_t::mutable_next() {
        return &_next;
    }
    inline void subpath_t::set_next(size_t i, uint32_t n) {
        _next[i] = n;
    }
    inline void subpath_t::add_next(uint32_t n) {
        _next.emplace_back(n);
    }
    inline void subpath_t::clear_next() {
        _next.clear();
    }
    inline size_t subpath_t::next_size() const {
        return _next.size();
    }
    inline bool subpath_t::has_next() const {
        return !_next.empty();
    }
    inline int32_t subpath_t::score() const {
        return _score;
    }
    inline void subpath_t::set_score(int32_t s) {
        _score = s;
    }
    inline const vector<connection_t>& subpath_t::connection() const {
        return _connection;
    }
    inline const connection_t& subpath_t::connection(size_t i) const {
        return _connection[i];
    }
    inline vector<connection_t>* subpath_t::mutable_connection() {
        return &_connection;
    }
    inline connection_t* subpath_t::mutable_connection(size_t i) {
        return &_connection[i];
    }
    inline void subpath_t::set_connection(size_t i, const connection_t& c) {
        _connection[i] = c;
    }
    inline connection_t* subpath_t::add_connection() {
        _connection.emplace_back();
        return &_connection.back();
    }
    inline void subpath_t::clear_connection() {
        _connection.clear();
    }
    inline size_t subpath_t::connection_size() const {
        return _connection.size();
    }
    inline bool subpath_t::has_connection() const {
        return !_connection.empty();
    }

    /*
     * multipath_alignment_t
     */
    inline const string& multipath_alignment_t::sequence() const {
        return _sequence;
    }
    inline string* multipath_alignment_t::mutable_sequence() {
        return &_sequence;
    }
    inline void multipath_alignment_t::set_sequence(const string& s) {
        _sequence = s;
    }
    inline const string& multipath_alignment_t::quality() const {
        return _quality;
    }
    inline string* multipath_alignment_t::mutable_quality() {
        return &_quality;
    }
    inline void multipath_alignment_t::set_quality(const string& q) {
        _quality = q;
    }
    inline const vector<subpath_t>& multipath_alignment_t::subpath() const {
        return _subpath;
    }
    inline const subpath_t& multipath_alignment_t::subpath(size_t i) const {
        return _subpath[i];
    }
    inline vector<subpath_t>* multipath_alignment_t::mutable_subpath() {
        return &_subpath;
    }
    inline subpath_t* multipath_alignment_t::mutable_subpath(size_t i) {
        return &_subpath[i];
    }
    inline subpath_t* multipath_alignment_t::add_subpath() {
        _subpath.emplace_back();
        return &_subpath.back();
    }
    inline void multipath_alignment_t::clear_subpath() {
        _subpath.clear();
    }
    inline size_t multipath_alignment_t::subpath_size() const {
        return _subpath.size();
    }
    inline int32_t multipath_alignment_t::mapping_quality() const {
        return _mapping_quality;
    }
    inline void multipath_alignment_t::set_mapping_quality(int32_t q) {
        _mapping_quality = q;
    }
    inline const vector<uint32_t>& multipath_alignment_t::start() const {
        return _start;
    }
    inline uint32_t multipath_alignment_t::start(size_t i) const {
        return _start[i];
    }
    inline vector<uint32_t>* multipath_alignment_t::mutable_start() {
        return &_start;
    }
    inline void multipath_alignment_t::set_start(size_t i, uint32_t s) {
        _start[i] = s;
    }
    inline void multipath_alignment_t::add_start(uint32_t s) {
        _start.emplace_back(s);
    }
    inline void multipath_alignment_t::clear_start() {
        _start.clear();
    }
    inline size_t multipath_alignment_t::start_size() const {
        return _start.size();
    }
    inline bool multipath_alignment_t::has_start() const {
        return !_start.empty();
    }
}


#endif /* multipath_alignment_hpp */





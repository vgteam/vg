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

// Declare the haplo::ScoreProvider we use for haplotype-aware traceback generation.
namespace haplo {
    class ScoreProvider;
}

namespace vg {

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
        inline int32_t score() const;
        inline void set_score(int32_t s);
    private:
        path_t _path;
        vector<uint32_t> _next;
        int32_t _score;
    };

    // TODO: the metadata could be removed and only added to the protobuf at serialization time
    class multipath_alignment_t {
    public:
        multipath_alignment_t() = default;
        multipath_alignment_t(const multipath_alignment_t&) = default;
        multipath_alignment_t(multipath_alignment_t&&) = default;
        ~multipath_alignment_t() = default;
        multipath_alignment_t& operator=(const multipath_alignment_t&) = default;
        multipath_alignment_t& operator=(multipath_alignment_t&&) = default;
        inline const string& sequence() const;
        inline string* mutable_sequence();
        inline void set_sequence(const string& s);
        inline const string& quality() const;
        inline string* mutable_quality();
        inline void set_quality(const string& q);
        inline const string& name() const;
        inline string* mutable_name();
        inline void set_name(const string& n);
        inline const string& sample_name() const;
        inline string* mutable_sample_name();
        inline void set_sample_name(const string& n);
        inline const string& read_group() const;
        inline string* mutable_read_group();
        inline void set_read_group(const string& g);
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
        inline const string& paired_read_name() const;
        inline string* mutable_paired_read_name();
        inline void set_paired_read_name(const string& n);
        // TODO: how to handle annotations?
        
    private:
        string _sequence;
        string _quality;
        string _name;
        string _sample_name;
        string _read_group;
        vector<subpath_t> _subpath;
        int32_t _mapping_quality;
        vector<uint32_t> _start;
        string _paired_read_name;
        map<string, string> _annotation;
    };

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

    /// Merges non-branching paths in a multipath alignment in place
    void merge_non_branching_subpaths(multipath_alignment_t& multipath_aln);

    /// Removes subpaths that have no aligned bases and adds in any implied edges crossing through them
    void remove_empty_subpaths(multipath_alignment_t& multipath_aln);
    
    /// Returns a vector whose elements are vectors with the indexes of the Subpaths in
    /// each connected component. An unmapped multipath_alignment_t with no subpaths produces an empty vector.
    vector<vector<int64_t>> connected_components(const multipath_alignment_t& multipath_aln);
    
    /// Extract the multipath_alignment_t consisting of the Subpaths with the given indexes
    /// into a new multipath_alignment_t object
    void extract_sub_multipath_alignment(const multipath_alignment_t& multipath_aln,
                                         const vector<int64_t>& subpath_indexes,
                                         multipath_alignment_t& sub_multipath_aln);
    
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
    inline int32_t subpath_t::score() const {
        return _score;
    }
    inline void subpath_t::set_score(int32_t s) {
        _score = s;
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
    inline const string& multipath_alignment_t::name() const {
        return _name;
    }
    inline string* multipath_alignment_t::mutable_name() {
        return &_name;
    }
    inline void multipath_alignment_t::set_name(const string& n) {
        _name = n;
    }
    inline const string& multipath_alignment_t::sample_name() const {
        return _sample_name;
    }
    inline string* multipath_alignment_t::mutable_sample_name() {
        return &_sample_name;
    }
    inline void multipath_alignment_t::set_sample_name(const string& n) {
        _sample_name = n;
    }
    inline const string& multipath_alignment_t::read_group() const {
        return _read_group;
    }
    inline string* multipath_alignment_t::mutable_read_group() {
        return &_read_group;
    }
    inline void multipath_alignment_t::set_read_group(const string& g) {
        _read_group = g;
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
    inline const string& multipath_alignment_t::paired_read_name() const {
        return _paired_read_name;
    }
    inline string* multipath_alignment_t::mutable_paired_read_name() {
        return &_paired_read_name;
    }
    inline void multipath_alignment_t::set_paired_read_name(const string& n) {
        _paired_read_name = n;
    }
}


#endif /* multipath_alignment_hpp */





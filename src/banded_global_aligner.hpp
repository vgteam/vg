//
//  banded_global_aligner.hpp
//  
//  Contains objects to support Aligner in performing banded global alignment against
//  a graph.
//

#ifndef banded_global_aligner_hpp
#define banded_global_aligner_hpp

#include <stdio.h>
#include <ctype.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <list>
#include <exception>

#include "handle.hpp"


using namespace std;

namespace vg {
    
    /**
     * This gets thrown when the aligner can't find any valid alignment in
     * the band that was requested.
     */
    class NoAlignmentInBandException : public exception {
        virtual const char* what() const noexcept;
        static const string message;
        int get_count();
    };
    
    /**
     * The outward-facing interface for banded global graph alignment. It computes optimal alignment
     * of a DNA sequence to a DAG with POA. The alignment will start at any source node in the graph and
     * end at any sink node. It is also restricted to falling within a certain diagonal band from the
     * start node. Any signed integer type can be used for the dynamic programming matrices, but there
     * are no checks for overflow.
     *
     * THIS IS A COMPONENT OF THE ALIGNER CLASS.
     *
     * Use Aligner::align_global_banded() instead.
     *
     */
    template <class IntType>
    class BandedGlobalAligner {
    public:
        /// Initializes banded alignment
        ///
        /// Args:
        ///  alignment                   empty alignment with a sequence (and possibly base qualities)
        ///  g                           graph to align to; no dangling edges are permitted
        ///  band_padding                width to expand band by
        ///  permissive_banding          expand band, not necessarily symmetrically, to allow all node paths
        ///  adjust_for_base_quality     perform base quality adjusted alignment (see QualAdjAligner)
        ///
        BandedGlobalAligner(Alignment& alignment, const HandleGraph& g,
                            int64_t band_padding, bool permissive_banding = false,
                            bool adjust_for_base_quality = false);
        
        
        /// Initializes banded multi-alignment, which computes the top scoring alternate alignments in addition
        /// to the optimal alignment
        ///
        /// Args:
        ///  alignment                   empty alignment with a sequence (and possibly base qualities)
        ///  g                           graph to align to; no dangling edges are permitted
        ///  alt_alignments              an empty vector to store alternate alignments in, the first element
        ///                              will be a copy of the primary alignment
        ///  max_multi_alns              the maximum number of alternate alignments (including the primary)
        ///  band_padding                width to expand band by
        ///  permissive_banding          expand band, not necessarily symmetrically, to allow all node paths
        ///  adjust_for_base_quality     perform base quality adjusted alignment (see QualAdjAligner)
        BandedGlobalAligner(Alignment& alignment, const HandleGraph& g,
                            vector<Alignment>& alt_alignments, int64_t max_multi_alns,
                            int64_t band_padding, bool permissive_banding = false,
                            bool adjust_for_base_quality = false);
        
        ~BandedGlobalAligner();
        
        /// Adds path and score to the alignment object given in the constructor. If a multi-alignment
        /// vector was also supplied, fills the vector with the top scoring alignments as well.
        ///
        /// Note: the score arguments are not <IntType> so that they can interact easily with the Aligner
        ///
        /// Args:
        ///  score_mat   matrix of match/mismatch scores from Aligner (if performing base quality
        ///              adjusted alignment, use QualAdjAligner's adjusted score matrix)
        ///  nt_table    table of indices by DNA char from Aligner
        ///  gap_open    gap open penalty from Algner (if performing base quality adjusted alignment,
        ///              use QualAdjAligner's scaled penalty)
        ///  gap_extend  gap extension penalty from Algner (if performing base quality adjusted alignment,
        ///              use QualAdjAligner's scaled penalty)
        void align(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend);
        
        
    private:
        
        class BAMatrix;
        class BABuilder;
        class AltTracebackStack;
        
        /// Matrices used in Smith-Waterman-Gotoh alignment algorithm
        enum matrix_t {Match, InsertCol, InsertRow};
        
        const HandleGraph& graph;
        /// The primary alignment
        Alignment& alignment;
        /// Vector for alternate alignments, or null if not making any
        vector<Alignment>* alt_alignments;
        /// Number of alignments to compute, including the optimal alignment
        int64_t max_multi_alns;
        /// Use base quality adjusted scoring for alignments?
        bool adjust_for_base_quality;
        
        /// Dynamic programming matrices for each node
        vector<BAMatrix*> banded_matrices;
        
        /// Map from node IDs to the index used in internal vectors
        unordered_map<int64_t, int64_t> node_id_to_idx;
        /// A topological ordering of the nodes
        vector<handle_t> topological_order;
        /// Source nodes in the graph
        vector<handle_t> source_nodes;
        /// Sink nodes in the graph
        vector<handle_t> sink_nodes;
        
        /// Internal constructor that the public constructors funnel into
        BandedGlobalAligner(Alignment& alignment, const HandleGraph& g,
                            vector<Alignment>* alt_alignments, int64_t max_multi_alns,
                            int64_t band_padding, bool permissive_banding = false,
                            bool adjust_for_base_quality = false);
        
        /// Traceback through dynamic programming matrices to compute alignment
        void traceback(int8_t* score_mat, int8_t* nt_table, int8_t gap_open, int8_t gap_extend, IntType min_inf);
        
        /// Constructor helper function: compute the longest and shortest path to a sink for each node
        void path_lengths_to_sinks(vector<int64_t>& shortest_path_to_sink, vector<int64_t>& longest_path_to_sink);
        /// Constructor helper function: compute which diagonals the bands cover on each node's matrix
        void find_banded_paths(bool permissive_banding, int64_t band_padding, vector<bool>& node_masked,
                               vector<pair<int64_t, int64_t>>& band_ends);
        /// Constructor helper function: compute the shortest path from a source to each node
        void shortest_seq_paths(vector<int64_t>& seq_lens_out);
    };

    /**
     * Represents the band from the DP matrix for one node in the graph
     *
     */
    template <class IntType>
    class BandedGlobalAligner<IntType>::BAMatrix {
        
    public:
        BAMatrix(Alignment& alignment, handle_t node, int64_t top_diag, int64_t bottom_diag,
                 const vector<BAMatrix*>& seeds, int64_t cumulative_seq_len);
        ~BAMatrix();
        
        /// Use DP to fill the band with alignment scores
        void fill_matrix(const HandleGraph& graph, int8_t* score_mat, int8_t* nt_table, int8_t gap_open,
                         int8_t gap_extend, bool qual_adjusted, IntType min_inf);
        
        void init_traceback_indexes(const HandleGraph& graph, int64_t& i, int64_t& j);
        
        void traceback(const HandleGraph& graph, BABuilder& builder, AltTracebackStack& traceback_stack,
                       int64_t& i, int64_t& j, matrix_t& mat, bool& in_lead_gap,
                       const int8_t* score_mat, const int8_t* nt_table, const int8_t gap_open, const int8_t gap_extend,
                       const bool qual_adjusted, IntType const min_inf);
        
        void traceback_over_edge(const HandleGraph& graph, BABuilder& builder, AltTracebackStack& traceback_stack,
                                 int64_t& i, int64_t& j, matrix_t& mat, bool& in_lead_gap, int64_t& node_id,
                                 const int8_t* score_mat, const int8_t* nt_table, const int8_t gap_open,
                                 const int8_t gap_extend, const bool qual_adjusted, IntType const min_inf);
        
        /// Debugging function
        void print_full_matrices(const HandleGraph& graph);
        /// Debugging function
        void print_rectangularized_bands(const HandleGraph& graph);
        
    private:

        /// The diagonals in the DP matrix that the band passes through with the bottom index inclusive
        int64_t top_diag;
        int64_t bottom_diag;
        
        handle_t node;
        
        Alignment& alignment;
        
        /// Length of shortest sequence leading to matrix from a source node
        int64_t cumulative_seq_len;
        
        /// Matrices for nodes with edges into this node
        vector<BAMatrix*> seeds;
        
        /// DP matrix
        IntType* match;
        /// DP matrix
        IntType* insert_col;
        /// DP matrix
        IntType* insert_row;
        
        /// Debugging function
        void print_matrix(const HandleGraph& graph, matrix_t which_mat);
        /// Debugging function
        void print_band(const HandleGraph& graph, matrix_t which_mat);
        
        friend class BABuilder;
        friend class AltTracebackStack; // not a fan of this one, but constructor ugly without it
        friend class BandedGlobalAligner; // also not a fan of this one but i have to refactor some
                                          // debug statements without it
    };
    
    /**
     * Specialized linked list stack that finds and keeps track of the top-scoring non-optimal alignments.
     * Sub-optimal alignments are represented by the locations where their traceback deviates from the optimal
     * traceback, with other steps of the traceback implicitly following the optimal path. New tracebacks can
     * be discovered while performing any of the tracebacks and added to the stack.
     */
    template <class IntType>
    class BandedGlobalAligner<IntType>::AltTracebackStack {
    public:
        AltTracebackStack(const HandleGraph& graph, int64_t max_multi_alns, int32_t empty_score,
                          unordered_set<BAMatrix*>& source_node_matrices,
                          unordered_set<BAMatrix*>& sink_node_matrices,
                          int8_t gap_open,
                          int8_t gap_extend,
                          IntType min_inf);
        ~AltTracebackStack();
        
        /// Get the start position of the current alignment and advance deflection pointer to the first deflection
        inline void get_alignment_start(int64_t& node_id, matrix_t& matrix);
        
        
        /// Are there any more alternate alignments in the stack?
        inline bool has_next();
        
        /// Advance to the next alternate trac
        inline void next_traceback_alignment();
        
        /// Is the next highest scoring alignment an empty path with no DP matrices?
        inline bool next_is_empty();
        
        /// Get the next empty path alignment
        inline void next_empty_alignment(Alignment& alignment);
        
        /// Check if a deflection from the current traceback is better than any alignments currently in the
        /// stack and if so insert it in the correct position
        inline void propose_deflection(const IntType score, const int64_t from_node_id, const int64_t row_idx,
                                       const int64_t col_idx, const int64_t to_node_id, const matrix_t to_matrix);
        
        /// Score of the current traceback
        inline IntType current_traceback_score();
        
        /// The prefix of the current traceback that consists of only empty nodes
        inline const list<int64_t>& current_empty_prefix();
        
        /// Are these the coordinates of the next deflection?
        inline bool at_next_deflection(int64_t node_id, int64_t row_idx, int64_t col_idx);
        
        /// Get the matrix to deflect to and advance to the next deflection
        inline BandedGlobalAligner::matrix_t deflect_to_matrix();
        /// Get the matrix and the node id to deflect to (for use at a node boundary)
        inline BandedGlobalAligner::matrix_t deflect_to_matrix(int64_t& to_node_id);
        
    private:
        class Deflection;
        
        /// Maximum number of alternate alignments to keep track of (including the optimal alignment)
        int64_t max_multi_alns;
        /// List of tuples that contain the scores of alternate alignments, the places where their
        /// traceback deviates from the optimum (Deflections), and all empty nodes that occur at the suffix
        /// of the traceback (if any). Maintains an invariant where stack is sorted in descending score order.
        list<tuple<vector<Deflection>, IntType, list<int64_t>>> alt_tracebacks;
        /// All of the paths through the graph that take only empty nodes
        list<list<int64_t>> empty_full_paths;
        int32_t empty_score;
        const HandleGraph& graph;
        
        /// Pointer to the traceback directions for the alignment we are currently tracing back
        typename list<tuple<vector<Deflection>, IntType, list<int64_t>>>::iterator curr_traceback;
        
        /// Pointer to the next deviation from the optimal traceback we will take
        typename vector<Deflection>::iterator curr_deflxn;
        
        /// Internal method for propose_deflection
        inline void insert_traceback(const vector<Deflection>& traceback_prefix, const IntType score,
                                     const int64_t from_node_id, const int64_t row_idx,
                                     const int64_t col_idx, const int64_t to_node_id, const matrix_t to_matrix,
                                     const list<int64_t>& empty_node_prefix);
    };
    
    /**
     * Represents a deviation from the optimal traceback that a sub-optimal traceback takes.
     */
    template <class IntType>
    class BandedGlobalAligner<IntType>::AltTracebackStack::Deflection {
    public:
        Deflection(const int64_t from_node_id, const int64_t row_idx, const int64_t col_idx,
                   const int64_t to_node_id, const matrix_t to_matrix);
        ~Deflection();
        
        /// Node ID where deflection occurs
        const int64_t from_node_id;
        /// Coordinate from rectangularized band (not original matrix) where deflection occurs
        const int64_t row_idx;
        /// Coordinate from rectangularized band (not original matrix) where deflection occurs
        const int64_t col_idx;
        /// Node to deflect to
        const int64_t to_node_id;
        /// Dynamic programming matrix deflect to
        const matrix_t to_matrix;
    };
    
    /**
     * Translates a traceback path into a Path object and stores it in an Alignment object
     */
    template <class IntType>
    class BandedGlobalAligner<IntType>::BABuilder {
    public:
        BABuilder(Alignment& alignment);
        ~BABuilder();
        
        /// Add next step in traceback
        void update_state(const HandleGraph& graph, matrix_t matrix, const handle_t& node, int64_t read_idx,
                          int64_t node_idx, bool empty_node_seq = false);
        /// Call after concluding traceback to finish adding edits to alignment
        void finalize_alignment(const list<int64_t>& empty_prefix);
        
    private:
        Alignment& alignment;
        
        list<Mapping> node_mappings;
        list<Edit> mapping_edits;
        
        matrix_t matrix_state;
        bool matching = false;
        id_t current_node_id = 0;
        string current_node_sequence = "";
        int64_t edit_length = 0;
        int64_t edit_read_end_idx = 0;
        
        void finish_current_edit();
        void finish_current_node();
    };

}

#endif /* banded_global_aligner_hpp */


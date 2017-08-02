#ifndef VG_GSSW_ALIGNER_HPP_INCLUDED
#define VG_GSSW_ALIGNER_HPP_INCLUDED

#include <algorithm>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include "gssw.h"
#include "vg.pb.h"
#include "vg.hpp"
#include "Variant.h"
#include "Fasta.h"
#include "path.hpp"
#include "utility.hpp"
#include "banded_global_aligner.hpp"

namespace vg {

    static const int8_t default_match = 1;
    static const int8_t default_mismatch = 4;
    static const int8_t default_gap_open = 6;
    static const int8_t default_gap_extension = 1;
    static const int8_t default_full_length_bonus = 0;
    static const int8_t default_max_scaled_score = 32;
    static const uint8_t default_max_qual_score = 255;
    static const double default_gc_content = 0.5;

    /**
     * The interface that any Aligner should implement, with some default implementations.
     */
    class BaseAligner {
    protected:
        BaseAligner() = default;
        ~BaseAligner();
        
        // for construction
        // needed when constructing an alignable graph from the nodes
        gssw_graph* create_gssw_graph(Graph& g, bool add_pinning_node, gssw_node** gssw_pinned_node_out);
        void topological_sort(list<gssw_node*>& sorted_nodes);
        void visit_node(gssw_node* node,
                        list<gssw_node*>& sorted_nodes,
                        set<gssw_node*>& unmarked_nodes,
                        set<gssw_node*>& temporary_marks);
        
        // create a reversed graph for left-pinned alignment
        void reverse_graph(Graph& g, Graph& reversed_graph_out);
        // reverse all node sequences (other aspects of graph object not unreversed)
        void unreverse_graph(Graph& graph);
        // convert graph mapping back into unreversed node positions
        void unreverse_graph_mapping(gssw_graph_mapping* gm);
        
        // alignment functions
        void gssw_mapping_to_alignment(gssw_graph* graph,
                                       gssw_graph_mapping* gm,
                                       Alignment& alignment,
                                       bool pinned,
                                       bool pin_left,
                                       bool print_score_matrices = false);
        string graph_cigar(gssw_graph_mapping* gm);
        
        double maximum_mapping_quality_exact(vector<double>& scaled_scores, size_t* max_idx_out);
        double maximum_mapping_quality_approx(vector<double>& scaled_scores, size_t* max_idx_out);
        double estimate_next_best_score(int length, double min_diffs);
        
        // must be called before querying mapping_quality
        void init_mapping_quality(double gc_content);
        
        // TODO: this algorithm has numerical problems, just removing it for now
        //vector<double> all_mapping_qualities_exact(vector<double> scaled_scores);
        
        
    public:

        double estimate_max_possible_mapping_quality(int length, double min_diffs, double next_min_diffs);
        
        /// Store optimal local alignment against a graph in the Alignment object.
        /// Gives the full length bonus separately on each end of the alignment.
        /// Assumes that graph is topologically sorted by node index.
        virtual void align(Alignment& alignment, Graph& g, bool print_score_matrices = false) = 0;
        
        // store optimal alignment against a graph in the Alignment object with one end of the sequence
        // guaranteed to align to a source/sink node
        //
        // pinning left means that that the alignment starts with the first base of the read sequence and
        // the first base of a source node sequence, pinning right means that the alignment starts with
        // the final base of the read sequence and the final base of a sink node sequence
        //
        // Gives the full length bonus only on the non-pinned end of the alignment.
        //
        // assumes that graph is topologically sorted by node index
        virtual void align_pinned(Alignment& alignment, Graph& g, bool pin_left) = 0;
        
        // store the top scoring pinned alignments in the vector in descending score order up to a maximum
        // number of alignments (including the optimal one). if there are fewer than the maximum number in
        // the return value, then it includes all alignments with a positive score. the optimal alignment
        // will be stored in both the vector and in the main alignment object
        //
        // assumes that graph is topologically sorted by node index
        virtual void align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, Graph& g,
                                        bool pin_left, int32_t max_alt_alns) = 0;
        
        // store optimal global alignment against a graph within a specified band in the Alignment object
        // permissive banding auto detects the width of band needed so that paths can travel
        // through every node in the graph
        virtual void align_global_banded(Alignment& alignment, Graph& g,
                                         int32_t band_padding = 0, bool permissive_banding = true) = 0;
        
        // store top scoring global alignments in the vector in descending score order up to a maximum number
        // of alternate alignments (including the optimal alignment). if there are fewer than the maximum
        // number of alignments in the return value, then the vector contains all possible alignments. the
        // optimal alignment will be stored in both the vector and the original alignment object
        virtual void align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments,
                                               Graph& g, int32_t max_alt_alns, int32_t band_padding = 0,
                                               bool permissive_banding = true) = 0;
                        
        /// Compute the score of an exact match in the given alignment, from the
        /// given offset, of the given length.
        virtual int32_t score_exact_match(const Alignment& aln, size_t read_offset, size_t length) = 0;
        
        // stores -10 * log_10(P_err) in alignment mapping_quality field where P_err is the
        // probability that the alignment is not the correct one (assuming that one of the alignments
        // in the vector is correct). alignments must have been created with this Aligner for quality
        // score to be valid
        void compute_mapping_quality(vector<Alignment>& alignments,
                                     int max_mapping_quality,
                                     bool fast_approximation,
                                     double cluster_mq,
                                     bool use_cluster_mq,
                                     int overlap_count,
                                     double mq_estimate,
                                     double identity_weight);
        // same function for paired reads, mapping qualities are stored in both alignments in the pair
        void compute_paired_mapping_quality(pair<vector<Alignment>, vector<Alignment>>& alignment_pairs,
                                            const vector<double>& frag_weights,
                                            int max_mapping_quality1,
                                            int max_mapping_quality2,
                                            bool fast_approximation,
                                            double cluster_mq,
                                            bool use_cluster_mq,
                                            int overlap_count1,
                                            int overlap_count2,
                                            double mq_estimate1,
                                            double mq_estimate2,
                                            double identity_weight);
        
        // Convert a score to an unnormalized log likelihood for the sequence.
        // Requires log_base to have been set.
        double score_to_unnormalized_likelihood_ln(double score);
        
        /// The longest gap detectable from a read position without soft-clipping
        size_t longest_detectable_gap(const Alignment& alignment, const string::const_iterator& read_pos) const;
        
        /// The longest gap detectable from any read position without soft-clipping
        size_t longest_detectable_gap(const Alignment& alignment) const;
        
        /// Use the score values in the aligner to score the given alignment,
        /// scoring gaps caused by jumping between between nodes using a custom
        /// gap length estimation function (which takes the from position, the
        /// to position, and a search limit in bp that happens to be the read
        /// length).
        ///
        /// May include full length bonus or not. TODO: bool flags are bad.
        virtual int32_t score_alignment(const Alignment& aln,
            const function<size_t(pos_t, pos_t, size_t)>& estimate_distance,
            bool strip_bonuses = false);
            
        /// Without necessarily rescoring the entire alignment, return the score
        /// of the given alignment with bonuses removed. Assumes that bonuses
        /// are actually included in the score.
        /// Needs to know if the alignment was pinned-end or not, and, if so, which end was pinned.
        virtual int32_t remove_bonuses(const Alignment& aln, bool pinned = false, bool pin_left = false);
        
        // members
        int8_t* nt_table = nullptr;
        int8_t* score_matrix = nullptr;
        int8_t match;
        int8_t mismatch;
        int8_t gap_open;
        int8_t gap_extension;
        int8_t full_length_bonus;
        
        // log of the base of the logarithm underlying the log-odds interpretation of the scores
        double log_base = 0.0;
        
    };
    
    /**
     * An ordinary aligner.
     */
    class Aligner : public BaseAligner {
    private:
        
        // internal function interacting with gssw for pinned and local alignment
        void align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, Graph& g,
                            bool pinned, bool pin_left, int32_t max_alt_alns,
                            bool print_score_matrices = false);
    public:
        
        Aligner(int8_t _match = default_match,
                int8_t _mismatch = default_mismatch,
                int8_t _gap_open = default_gap_open,
                int8_t _gap_extension = default_gap_extension,
                int8_t _full_length_bonus = default_full_length_bonus,
                double _gc_content = default_gc_content);
        ~Aligner(void) = default;
        
        /// Store optimal local alignment against a graph in the Alignment object.
        /// Gives the full length bonus separately on each end of the alignment.
        /// Assumes that graph is topologically sorted by node index.
        void align(Alignment& alignment, Graph& g, bool print_score_matrices = false);
        
        // store optimal alignment against a graph in the Alignment object with one end of the sequence
        // guaranteed to align to a source/sink node
        //
        // pinning left means that that the alignment starts with the first base of the read sequence and
        // the first base of a source node sequence, pinning right means that the alignment starts with
        // the final base of the read sequence and the final base of a sink node sequence
        //
        // Gives the full length bonus only on the non-pinned end of the alignment.
        //
        // assumes that graph is topologically sorted by node index
        void align_pinned(Alignment& alignment, Graph& g, bool pin_left);
                
        // store the top scoring pinned alignments in the vector in descending score order up to a maximum
        // number of alignments (including the optimal one). if there are fewer than the maximum number in
        // the return value, then it includes all alignments with a positive score. the optimal alignment
        // will be stored in both the vector and in the main alignment object
        //
        // assumes that graph is topologically sorted by node index
        void align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, Graph& g,
                                bool pin_left, int32_t max_alt_alns);
        
        // store optimal global alignment against a graph within a specified band in the Alignment object
        // permissive banding auto detects the width of band needed so that paths can travel
        // through every node in the graph
        void align_global_banded(Alignment& alignment, Graph& g,
                                 int32_t band_padding = 0, bool permissive_banding = true);
        
        // store top scoring global alignments in the vector in descending score order up to a maximum number
        // of alternate alignments (including the optimal alignment). if there are fewer than the maximum
        // number of alignments in the return value, then the vector contains all possible alignments. the
        // optimal alignment will be stored in both the vector and the original alignment object
        void align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, Graph& g,
                                       int32_t max_alt_alns, int32_t band_padding = 0, bool permissive_banding = true);
        
        
        int32_t score_exact_match(const Alignment& aln, size_t read_offset, size_t length);
        int32_t score_exact_match(const string& sequence) const;
        int32_t score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end) const;

    };

    /**
     * An aligner that uses read base qualities to adjust its scores and alignments.
     */
    class QualAdjAligner : public BaseAligner {
    public:
        
        QualAdjAligner(int8_t _match = default_match,
                       int8_t _mismatch = default_mismatch,
                       int8_t _gap_open = default_gap_open,
                       int8_t _gap_extension = default_gap_extension,
                       int8_t _full_length_bonus = default_full_length_bonus,
                       int8_t _max_scaled_score = default_max_scaled_score,
                       uint8_t _max_qual_score = default_max_qual_score,
                       double gc_content = default_gc_content);

        ~QualAdjAligner(void) = default;

        // base quality adjusted counterparts to functions of same name from Aligner
        void align(Alignment& alignment, Graph& g, bool print_score_matrices = false);
        void align_global_banded(Alignment& alignment, Graph& g,
                                 int32_t band_padding = 0, bool permissive_banding = true);
        void align_pinned(Alignment& alignment, Graph& g, bool pin_left);
        void align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, Graph& g,
                                       int32_t max_alt_alns, int32_t band_padding = 0, bool permissive_banding = true);
        void align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, Graph& g,
                                bool pin_left, int32_t max_alt_alns);
        
        
        void init_mapping_quality(double gc_content);
        
        int32_t score_exact_match(const Alignment& aln, size_t read_offset, size_t length);
        int32_t score_exact_match(const string& sequence, const string& base_quality) const;
        int32_t score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end,
                                  string::const_iterator base_qual_begin) const;
        
        uint8_t max_qual_score;
        int8_t scale_factor;
        
    private:

        void align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, Graph& g,
                            bool pinned, bool pin_left, int32_t max_alt_alns,
                            bool print_score_matrices = false);
        

    };
} // end namespace vg

#endif

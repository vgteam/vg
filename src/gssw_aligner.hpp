#ifndef GSSW_ALIGNER_H
#define GSSW_ALIGNER_H

#include <vector>
#include <set>
#include <string>
#include "gssw.h"
#include "vg.pb.h"
#include "Variant.h"
#include "Fasta.h"
#include "path.hpp"

static const int8_t default_match = 1;
static const int8_t default_mismatch = 4;
static const int8_t default_gap_open = 6;
static const int8_t default_gap_extension = 1;
static const int8_t default_max_scaled_score = 32;
static const uint8_t default_max_qual_score = 255;
static const double default_gc_content = 0.5;

namespace vg {

    class Aligner {
    protected:
        gssw_graph* create_gssw_graph(Graph& g);
        
    public:
        
        Aligner(int32_t _match = default_match,
                int32_t _mismatch = default_mismatch,
                int32_t _gap_open = default_gap_open,
                int32_t _gap_extension = default_gap_extension);
        ~Aligner(void);

        // for construction
        // needed when constructing an alignable graph from the nodes
        void topological_sort(list<gssw_node*>& sorted_nodes);
        void visit_node(gssw_node* node,
                        list<gssw_node*>& sorted_nodes,
                        set<gssw_node*>& unmarked_nodes,
                        set<gssw_node*>& temporary_marks);

        // alignment functions
        void align(Alignment& alignment, Graph& g, bool print_score_matrices = false);
        void gssw_mapping_to_alignment(gssw_graph* graph,
                                       gssw_graph_mapping* gm,
                                       Alignment& alignment,
                                       bool print_score_matrices = false);
        string graph_cigar(gssw_graph_mapping* gm);
        
        // must be called before querying mapping_quality
        void init_mapping_quality(double gc_content);
        // returns -log_10(P_err) where P_err is the probability that the maximum scoring alignment
        // in this vector is not the correct one (assuming that one of the alignments is correct)
        // alignments must have been created with this Aligner for quality score to be valid
        uint8_t mapping_quality(vector<Alignment>& alignments);
        // fast approximation of mapping_quality
        uint8_t mapping_quality_approx(vector<Alignment>& alignments);

        // members
        map<int64_t, gssw_node*> nodes;
        int8_t* nt_table;
        int8_t* score_matrix;
        int32_t match;
        int32_t mismatch;
        int32_t gap_open;
        int32_t gap_extension;
        
        // log of the base of the logarithm underlying the log-odds interpretation of the scores
        double log_base;

    };

    class QualAdjAligner : public Aligner {
    public:
        
        QualAdjAligner(int8_t _match = default_match,
                       int8_t _mismatch = default_mismatch,
                       int8_t _gap_open = default_gap_open,
                       int8_t _gap_extension = default_gap_extension,
                       int8_t _max_scaled_score = default_max_scaled_score,
                       uint8_t _max_qual_score = default_max_qual_score,
                       double gc_content = default_gc_content);

        ~QualAdjAligner(void);

        void align(Alignment& alignment, Graph& g, bool print_score_matrices = false);
        void init_mapping_quality(double gc_content);

        uint8_t max_qual_score;
        int8_t scaled_gap_open;
        int8_t scaled_gap_extension;
        int8_t* adjusted_score_matrix;
        
    private:
        void init_quality_adjusted_scores(int8_t _max_scaled_score,
                                          uint8_t _max_qual_score,
                                          double gc_content);

    };
} // end namespace vg

#endif

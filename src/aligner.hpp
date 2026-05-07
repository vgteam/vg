#ifndef VG_ALIGNER_HPP_INCLUDED
#define VG_ALIGNER_HPP_INCLUDED

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vg/vg.pb.h>

#include "gssw.h"
#include "Variant.h"
#include "Fasta.h"
#include "handle.hpp"
#include "path.hpp"
#include "dozeu_interface.hpp"
#include "deletion_aligner.hpp"
#include "alignment_scorer.hpp"
#include "mapping_quality_calculator.hpp"

// #define BENCH
// #include "bench.h"

namespace vg {

    /**
     * The abstract interface that any Aligner should implement.
     */
    class BaseAligner {
    public:

        /// Store optimal local alignment against a graph in the Alignment object.
        /// Gives the full length bonus separately on each end of the alignment.
        virtual void align(Alignment& alignment, const HandleGraph& g, bool traceback_aln) const = 0;
    };

    /**
     * The basic GSSW-based core aligner implementation, which can then be
     * quality-adjusted or not. Owns its scorer and MAPQ calculator; all
     * scoring and MAPQ work goes through those, not through the aligner.
     *
     * You aren't meant to instantiate this directly; you should have an
     * Aligner or a QualAdjAligner as the concrete class instead.
     */
    class GSSWAligner : public BaseAligner {
    
    protected:
        /// Build a GSSWAligner around a MatrixAlignmentScorer.
        ///
        /// This is the protected constructor menat to be called by subclasses.
        GSSWAligner(std::unique_ptr<MatrixAlignmentScorer> owned_scorer, double gc_content);
        ~GSSWAligner();

    public:
        // Move-only: the owned scorer/MQ calc are unique_ptrs.
        GSSWAligner(GSSWAligner&&) = default;
        GSSWAligner& operator=(GSSWAligner&&) = default;
        GSSWAligner(const GSSWAligner&) = delete;
        GSSWAligner& operator=(const GSSWAligner&) = delete;
        
    protected: 

        // for construction
        // needed when constructing an alignable graph from the nodes
        gssw_graph* create_gssw_graph(const HandleGraph& g) const;

        // identify the IDs of nodes that should be used as pinning points in GSSW for pinned
        // alignment ((i.e. non-empty nodes as close as possible to sinks))
        unordered_set<id_t> identify_pinning_points(const HandleGraph& graph) const;

        // convert graph mapping back into unreversed node positions
        void unreverse_graph_mapping(gssw_graph_mapping* gm) const;
        // convert from graph sequences back into unrereversed form
        void unreverse_graph(gssw_graph* graph) const;

        // alignment functions
        void gssw_mapping_to_alignment(gssw_graph* graph,
                                       gssw_graph_mapping* gm,
                                       Alignment& alignment,
                                       bool pinned,
                                       bool pin_left) const;
        string graph_cigar(gssw_graph_mapping* gm) const;

    public:
        /// store optimal alignment against a graph in the Alignment object with one end of the sequence
        /// guaranteed to align to a source/sink node. if xdrop is selected, use the xdrop heuristic, which
        /// does not guarantee an optimal alignment.
        ///
        /// pinning left means that that the alignment starts with the first base of the read sequence and
        /// the first base of a source node sequence, pinning right means that the alignment starts with
        /// the final base of the read sequence and the final base of a sink node sequence
        ///
        /// Gives the full length bonus only on the non-pinned end of the alignment.
        virtual void align_pinned(Alignment& alignment, const HandleGraph& g, bool pin_left, bool xdrop = false,
                                  uint16_t xdrop_max_gap_length = default_xdrop_max_gap_length) const = 0;

        /// store the top scoring pinned alignments in the vector in descending score order up to a maximum
        /// number of alignments (including the optimal one). if there are fewer than the maximum number in
        /// the return value, then it includes all alignments with a positive score. the optimal alignment
        /// will be stored in both the vector and in the main alignment object
        virtual void align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                        bool pin_left, int32_t max_alt_alns) const = 0;

        /// Store optimal global alignment against a graph within a specified band in the Alignment object.
        /// Permissive banding auto detects the width of band needed so that paths can travel
        /// through every node in the graph.
        ///
        /// Throws BandMatricesTooBigException if the max_cells limit on DP matric size is hit.
        virtual void align_global_banded(Alignment& alignment, const HandleGraph& g,
                                         int32_t band_padding = 0, bool permissive_banding = true,
                                         uint64_t max_cells = std::numeric_limits<uint64_t>::max()) const = 0;

        /// Store top scoring global alignments in the vector in descending score order up to a maximum number
        /// of alternate alignments (including the optimal alignment). If there are fewer than the maximum
        /// number of alignments in the return value, then the vector contains all possible alignments. The
        /// optimal alignment will be stored in both the vector and the original alignment object.
        ///
        /// When multiple alignments have the same score, they are ordered deterministically but
        /// arbitrarily.
        ///
        /// Throws BandMatricesTooBigException if the max_cells limit on DP matric size is hit.
        virtual void align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments,
                                               const HandleGraph& g, int32_t max_alt_alns, int32_t band_padding = 0,
                                               bool permissive_banding = true,
                                               uint64_t max_cells = std::numeric_limits<uint64_t>::max()) const = 0;
        /// xdrop aligner
        virtual void align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<MaximalExactMatch>& mems,
                                 bool reverse_complemented, uint16_t max_gap_length = default_xdrop_max_gap_length) const = 0;

        /// xdrop aligner, but with a precomputed topological order on the graph, which need not include
        /// all of the graph's handles and which may contain both orientations of a handle
        virtual void align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<handle_t>& order,
                                 const vector<MaximalExactMatch>& mems, bool reverse_complemented,
                                 uint16_t max_gap_length = default_xdrop_max_gap_length) const = 0;
    
   
        /// Alignment scorer that represents the scoring scheme used for
        /// computing alignments.
        std::unique_ptr<MatrixAlignmentScorer> scorer;

        /// Mapping quality calculator for computing mapping qualities under
        /// the scorign scheme we use. Not used internally, but exposed so that
        /// aligner clinets have easy access to a MAPQ calculator for free and
        /// don't need to build their own.
        std::unique_ptr<MappingQualityCalculator> mapq_calc;

    protected:
        /// Widget for producing alignments that we know in advance will be
        /// complete deletions.
        DeletionAligner deletion_aligner;

    };

    /**
     * An ordinary aligner.
     */
    class Aligner : public GSSWAligner {

    public:

        Aligner(const int8_t* _score_matrix = default_score_matrix,
                int8_t _gap_open = default_gap_open,
                int8_t _gap_extension = default_gap_extension,
                int8_t _full_length_bonus = default_full_length_bonus,
                double _gc_content = default_gc_content);
        // Move-only (owned scorer/MQ calc are unique_ptrs).
        Aligner(Aligner&&) = default;
        Aligner& operator=(Aligner&&) = default;

        /// Store optimal local alignment against a graph in the Alignment object.
        /// Gives the full length bonus separately on each end of the alignment.
        void align(Alignment& alignment, const HandleGraph& g, bool traceback_aln) const;

        /// Align against a subgraph induced by a subset of nodes. The topological
        /// order of the handles in the subgraph must be provided.
        /// Store optimal local alignment in the Alignment object.
        /// Gives the full length bonus separately on each end of the alignment.
        void align(Alignment& alignment, const HandleGraph& g,
                   const std::vector<handle_t>& topological_order) const;

        void align_pinned(Alignment& alignment, const HandleGraph& g, bool pin_left, bool xdrop = false,
                          uint16_t xdrop_max_gap_length = default_xdrop_max_gap_length) const;

        void align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                bool pin_left, int32_t max_alt_alns) const;

        void align_global_banded(Alignment& alignment, const HandleGraph& g,
                                 int32_t band_padding = 0, bool permissive_banding = true,
                                 uint64_t max_cells = std::numeric_limits<uint64_t>::max()) const;

        void align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                       int32_t max_alt_alns, int32_t band_padding = 0, bool permissive_banding = true,
                                       uint64_t max_cells = std::numeric_limits<uint64_t>::max()) const;

        void align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<MaximalExactMatch>& mems,
                         bool reverse_complemented, uint16_t max_gap_length = default_xdrop_max_gap_length) const;

        void align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<handle_t>& order,
                         const vector<MaximalExactMatch>& mems, bool reverse_complemented,
                         uint16_t max_gap_length = default_xdrop_max_gap_length) const;

    private:

        // internal function interacting with gssw for pinned and local alignment
        void align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, const HandleGraph& g,
                            bool pinned, bool pin_left, int32_t max_alt_alns,
                            bool traceback_aln) const;

        // members
        vector<XdropAligner> xdrops;
    };

    /**
     * An aligner that uses read base qualities to adjust its scores and alignments.
     */
    class QualAdjAligner : public GSSWAligner {
    public:

        QualAdjAligner(const int8_t* _score_matrix = default_score_matrix,
                       int8_t _gap_open = default_gap_open,
                       int8_t _gap_extension = default_gap_extension,
                       int8_t _full_length_bonus = default_full_length_bonus,
                       double _gc_content = default_gc_content);
        // Move-only.
        QualAdjAligner(QualAdjAligner&&) = default;
        QualAdjAligner& operator=(QualAdjAligner&&) = default;

        void align(Alignment& alignment, const HandleGraph& g, bool traceback_aln) const;
        void align_global_banded(Alignment& alignment, const HandleGraph& g,
                                 int32_t band_padding = 0, bool permissive_banding = true,
                                 uint64_t max_cells = std::numeric_limits<uint64_t>::max()) const;
        void align_pinned(Alignment& alignment, const HandleGraph& g, bool pin_left, bool xdrop = false,
                          uint16_t xdrop_max_gap_length = default_xdrop_max_gap_length) const;
        void align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                       int32_t max_alt_alns, int32_t band_padding = 0, bool permissive_banding = true,
                                       uint64_t max_cells = std::numeric_limits<uint64_t>::max()) const;
        void align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                bool pin_left, int32_t max_alt_alns) const;

        void align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<MaximalExactMatch>& mems,
                         bool reverse_complemented, uint16_t max_gap_length = default_xdrop_max_gap_length) const;
        void align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<handle_t>& order,
                         const vector<MaximalExactMatch>& mems, bool reverse_complemented,
                         uint16_t max_gap_length = default_xdrop_max_gap_length) const;

    protected:

        // internal function interacting with gssw for pinned and local alignment
        void align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, const HandleGraph& g,
                            bool pinned, bool pin_left, int32_t max_alt_alns,
                            bool traceback_aln) const;

        /// Convenience accessor: scorer downcast to its concrete type so we
        /// can reach the qual-adjusted full-length bonus table without going
        /// through a virtual.
        const QualAdjAlignmentScorer* qa_scorer() const {
            return static_cast<const QualAdjAlignmentScorer*>(scorer.get());
        }

        // members
        vector<QualAdjXdropAligner> xdrops;
    };


    /**
     * Holds a set of alignment scores, and has methods to produce aligners of various types on demand, using those scores.
     * Provides a get_aligner() method to get ahold of a useful, possibly quality-adjusted Aligner.
     * Base functionality that is shared between alignment and surjections
     */
    class AlignerClient {
    protected:

        /// Create an AlignerClient, which creates the default aligner instances,
        /// which can depend on a GC content estimate.
        AlignerClient(double gc_content_estimate = vg::default_gc_content);

        /// Get the appropriate aligner to use, based on
        /// adjust_alignments_for_base_quality. By setting have_qualities to false,
        /// you can force the non-quality-adjusted aligner, for reads that lack
        /// quality scores.
        const GSSWAligner* get_aligner(bool have_qualities = true) const;

        // Sometimes you really do need the two kinds of aligners, to pass to code
        // that expects one or the other.
        const QualAdjAligner* get_qual_adj_aligner() const;
        const Aligner* get_regular_aligner() const;

    public:

        /// Set all the aligner scoring parameters and create the stored aligner instances.
        virtual void set_alignment_scores(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);

        /// Set the algner scoring parameters and create the stored aligner instances. The
        /// stream should contain a 4 x 4 whitespace-separated substitution matrix (in the
        /// order ACGT)
        virtual void set_alignment_scores(std::istream& matrix_stream, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);

        /// Set the algner scoring parameters and create the stored aligner instances. The
        /// score matrix should by a 4 x 4 array in the order (ACGT).
        /// Other overloads of set_alignment_scores all call this one.
        /// Note that an override of this method can't be called from the
        /// constructor, so when overriding it, make sure to also do your extra
        /// work in the constructor.
        virtual void set_alignment_scores(const int8_t* score_matrix, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus);

        /// Allocates an array to hold a 4x4 substitution matrix and returns it
        static int8_t* parse_matrix(std::istream& matrix_stream);

        bool adjust_alignments_for_base_quality = false; // use base quality adjusted alignments

    private:

        // GSSW aligners
        unique_ptr<QualAdjAligner> qual_adj_aligner;
        unique_ptr<Aligner> regular_aligner;

    protected:
        // GC content estimate that we need for building the aligners.
        double gc_content_estimate;
    };

} // end namespace vg

#endif

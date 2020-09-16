
#ifndef VG_TRANSCRIPTOME_HPP_INCLUDED
#define VG_TRANSCRIPTOME_HPP_INCLUDED

#include <algorithm>
#include <mutex>

#include <google/protobuf/util/message_differencer.h>
#include <gbwt/dynamic_gbwt.h>
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include "bdsg/overlays/path_position_overlays.hpp"

#include "../vg.hpp"
#include "../types.hpp"

namespace vg {

using namespace std;


typedef vector<gbwt::node_type> exon_nodes_t;
typedef vector<gbwt::size_type> thread_ids_t;


/**
 * Data structure that defines a transcript annotation.
 */
struct Exon {

    /// Exon coordinates (start and end) on the chromosome/contig.
    pair<int32_t, int32_t> coordinates;

    /// Exon border node offsets (last position in upstream intron and
    /// first position in downstream intron) on a variation graph. 
    pair<uint32_t, uint32_t> border_offsets;

    /// Exon border reference path steps (last position in upstream intron and
    /// first position in downstream intron) on a variation graph. 
    pair<step_handle_t, step_handle_t> border_steps;
};


/**
 * Data structure that defines a transcript annotation.
 */
struct Transcript {

    /// Transcript name.
    const string name;

    /// Is transcript in reverse direction (strand == '-').
    const bool is_reverse;

    /// Name of chromosome/contig where transcript exist.
    const string chrom;

    /// Length of chromosome/contig where transcript exist.
    const uint32_t chrom_length;

    /// Transcript exons.
    vector<Exon> exons;

    Transcript(const string & name_in, const bool is_reverse_in, const string & chrom_in, const uint32_t & chrom_length_in) : name(name_in), is_reverse(is_reverse_in), chrom(chrom_in), chrom_length(chrom_length_in) {

        assert(chrom_length > 0);
    }
};


/**
 * Data structure that defines a base transcript path.
 */ 
struct TranscriptPath {

    /// Transcript path name.
    string name;

    /// Transcript origin name.
    const string transcript_origin;

    /// Reference origin name.
    string reference_origin;

    /// Haplotype origin ids.
    vector<uint32_t> haplotype_origin_ids;

    /// Embedded path origin names.
    string path_origin_names;
    
    TranscriptPath(const string & transcript_origin_in) : transcript_origin(transcript_origin_in) {

        name = "";
        reference_origin = "";
        path_origin_names = "";
    }
};

/**
 * Data structure that defines an edited transcript path.
 */ 
struct EditedTranscriptPath : public TranscriptPath {

    /// Transcript path.
    Path path;

    EditedTranscriptPath(const string & transcript_origin_in) : TranscriptPath(transcript_origin_in) {}
};

/**
 * Data structure that defines a completed transcript path.
 */ 
struct CompletedTranscriptPath : public TranscriptPath {

    /// Transcript path.
    vector<handle_t> path;

    CompletedTranscriptPath(const string & transcript_origin_in) : TranscriptPath(transcript_origin_in) {}
};



/**
 * Class that defines a transcriptome represented by a set of transcript paths.
 */
class Transcriptome {

    public:

        Transcriptome(const string &, const bool);   

        /// Number of threads used for transcript path construction. 
        int32_t num_threads = 1;

        /// Feature type to parse in the gtf/gff file. Parse all types if empty. 
        string feature_type; 

        /// Attribute tag used to parse the transcript id/name in the gtf/gff file. 
        string transcript_tag;

        /// Use all paths embedded in the graph for transcript path construction. 
        bool use_all_paths = false;

        /// Use reference paths embedded in the graph for transcript path construction. 
        bool use_reference_paths = false;

        /// Collapse identical transcript paths.
        bool collapse_transcript_paths = true;

        /// Add splice-junstions from a intron BED file. 
        /// Returns number of parsed introns. 
        int32_t add_intron_splice_junctions(istream & intron_stream, unique_ptr<gbwt::GBWT> & haplotype_index);

        /// Add splice-junstions from a transcript gtf/gff3 file. 
        /// Returns number of parsed transcripts.  
        int32_t add_transcript_splice_junctions(istream & transcript_stream, unique_ptr<gbwt::GBWT> & haplotype_index);

        /// Constructs transcript paths by projecting transcripts from a gtf/gff file onto 
        /// embedded paths in a variation graph and/or haplotypes in a GBWT index. Augments 
        /// graph with transcriptome splice-junctions. Returns number of transcript paths added.   
        int32_t add_transcripts(istream & transcript_stream, const gbwt::GBWT & haplotype_index);

        /// Returns transcript paths.
        const vector<CompletedTranscriptPath> & transcript_paths() const;

        /// Returns number of transcript paths.
        int32_t size() const;

        /// Returns spliced variation graph.
        const MutablePathDeletableHandleGraph & splice_graph() const; 

        /// Returns true if nodes in the spliced variation graph 
        /// have been updated (e.g. split) since parsed.
        bool splice_graph_node_updated() const;

        /// Removes non-transcribed (not in transcript paths) nodes.
        /// Optionally create new reference paths that only include
        /// trancribed nodes and edges.
        void remove_non_transcribed(const bool new_reference_paths);

        /// Topological sort and compact graph.
        void compact_ordered();

        /// Embeds transcript paths in spliced variation graph.  
        /// Returns number of paths embedded.
        int32_t embed_transcript_paths(const bool add_reference_paths, const bool add_non_reference_paths);

        /// Add transcript paths as threads in GBWT index.
        /// Returns number of added threads.
        int32_t add_transcripts_to_gbwt(gbwt::GBWTBuilder * gbwt_builder, const bool output_reference_transcripts, const bool add_bidirectional) const;
        
        /// Writes transcript path sequences to a fasta file.  
        /// Returns number of written sequences.
        int32_t write_sequences(ostream * fasta_ostream, const bool output_reference_transcripts) const;

        /// Writes origin info on transcripts to tsv file.
        /// Returns number of written transcripts.
        int32_t write_info(ostream * tsv_ostream, const gbwt::GBWT & haplotype_index, const bool output_reference_transcripts) const;

        /// Writes spliced variation graph to vg file
        void write_splice_graph(ostream * graph_ostream) const;

    private:

        /// Transcriptome represented by a set of transcript paths. 
        vector<CompletedTranscriptPath> _transcript_paths;
        mutex mutex_transcript_paths;

        /// Spliced variation graph.
        unique_ptr<MutablePathDeletableHandleGraph> _splice_graph;
        mutex mutex_splice_graph;

        /// Have nodes in the spliced variation graph been
        /// updated (e.g. split) since parsed.
        bool _splice_graph_node_updated;

        /// Parse BED file of introns.
        vector<Transcript> parse_introns(istream & intron_stream, const bdsg::PositionOverlay & graph_path_pos_overlay) const;

        /// Parse gtf/gff3 file of transcripts.
        vector<Transcript> parse_transcripts(istream & transcript_stream, const bdsg::PositionOverlay & graph_path_pos_overlay) const;

        /// Returns the mean node length of the graph
        float mean_node_length() const;

        /// Finds the position of each end of a exon on a path in the  
        /// variation graph and adds the exon to a transcript.
        void add_exon(Transcript * transcript, const pair<int32_t, int32_t> & exon_pos, const bdsg::PositionOverlay & graph_path_pos_overlay) const;

        /// Reverses exon order if the transcript is on the reverse strand and the exons 
        /// are ordered in reverse.
        void reorder_exons(Transcript * transcript) const;

        /// Constructs edited transcript paths from a set of 
        /// reference transcripts. 
        list<EditedTranscriptPath> construct_edited_transcript_paths(const vector<Transcript> & transcripts, const bdsg::PositionOverlay & graph_path_pos_overlay) const;

        /// Threaded edited transcript path construction.
        void construct_edited_transcript_paths_callback(list<EditedTranscriptPath> * edited_transcript_paths, mutex * edited_transcript_paths_mutex, const bdsg::PositionOverlay & graph_path_pos_overlay, const int32_t thread_idx, const vector<Transcript> & transcripts) const;

        /// Constructs transcript paths by projecting transcripts onto embedded paths in
        /// a variation graph and/or haplotypes in a GBWT index. 
        void project_and_add_transcripts(const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const bdsg::PositionOverlay & graph_path_pos_overlay, const float mean_node_length);

        /// Threaded transcript projecting.
        void project_and_add_transcripts_callback(const int32_t thread_idx, const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const bdsg::PositionOverlay & graph_path_pos_overlay, const float mean_node_length);

        /// Projects transcripts onto haplotypes in a GBWT index and returns resulting transcript paths.
        list<EditedTranscriptPath> project_transcript_gbwt(const Transcript & cur_transcript, const gbwt::GBWT & haplotype_index, const float mean_node_length) const;

        /// Extracts all unique haplotype paths between two nodes from a GBWT index and returns the 
        /// resulting paths and the corresponding haplotype ids for each path.
        vector<pair<exon_nodes_t, thread_ids_t> > get_exon_haplotypes(const vg::id_t start_node, const vg::id_t end_node, const gbwt::GBWT & haplotype_index, const int32_t expected_length) const;

        /// Projects transcripts onto embedded paths in a variation graph and returns resulting transcript paths.
        list<EditedTranscriptPath> project_transcript_embedded(const Transcript & cur_transcript, const bdsg::PositionOverlay & graph_path_pos_overlay, const bool reference_only) const;

        /// Adds new transcript paths to current set. Has argument to only add unique paths.
        void append_transcript_paths(list<CompletedTranscriptPath> * completed_transcript_path, list<CompletedTranscriptPath> * new_completed_transcript_paths, const bool add_unqiue_paths_only) const;

        /// Constructs completed transcripts paths from 
        /// edited transcript paths. Checks that the
        /// paths contain no edits compared to the graph.
        list<CompletedTranscriptPath> construct_completed_transcript_paths(const list<EditedTranscriptPath> & edited_transcript_paths) const;

        /// Convert a path to a vector of handles. Checks that 
        /// the path is complete (i.e. only consist of whole nodes).
        vector<handle_t> path_to_handles(const Path & path, const EditedTranscriptPath & test) const;

        /// Checks whether transcript path only consist of 
        /// whole nodes (complete). 
        bool has_novel_exon_boundaries(const list<EditedTranscriptPath> & edited_transcript_paths, const bool include_transcript_ends) const;

        /// Augments the variation graph with transcript path exon boundaries and 
        /// splice-junctions. Updates threads in gbwt index to match the augmented graph. 
        void augment_splice_graph(list<EditedTranscriptPath> * edited_transcript_paths, unique_ptr<gbwt::GBWT> & haplotype_index, const bool break_at_transcript_ends);

        /// Update threads in gbwt index using graph translations. 
        void update_haplotype_index(unique_ptr<gbwt::GBWT> & haplotype_index, const vector<Translation> & translation) const;

        /// Adds transcript path splice-junction edges to splice graph
        void add_splice_junction_edges(const list<EditedTranscriptPath> & edited_transcript_paths);
        void add_splice_junction_edges(const vector<CompletedTranscriptPath> & completed_transcript_paths);
};

}

#endif

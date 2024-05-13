
#ifndef VG_TRANSCRIPTOME_HPP_INCLUDED
#define VG_TRANSCRIPTOME_HPP_INCLUDED

#include <algorithm>
#include <mutex>
#include <functional>

#include <google/protobuf/util/message_differencer.h>
#include <gbwt/dynamic_gbwt.h>
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include <bdsg/overlays/path_position_overlays.hpp>
#include <sparsepp/spp.h>

#include "vg.hpp"
#include "types.hpp"
#include "gbwt_helper.hpp"

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
    /// first position in downstream intron) on a graph. 
    pair<uint32_t, uint32_t> border_offsets;

    /// Exon border reference path steps (last position in upstream intron and
    /// first position in downstream intron) on a graph. 
    pair<step_handle_t, step_handle_t> border_steps;
};

/**
 * Data structure that defines a transcript annotation.
 */
struct Transcript {

    /// Transcript name.
    string name;

    /// Is transcript in reverse direction (strand == '-').
    bool is_reverse;

    /// Name of chromosome/contig where transcript exist.
    string chrom;

    /// Length of chromosome/contig where transcript exist.
    uint32_t chrom_length;

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

    /// Transcript names.
    vector<string> transcript_names;

    /// Embedded path names and reference path origin.
    vector<pair<string, bool> > embedded_path_names;

    /// Haplotype gbwt ids and reference path origin.
    vector<pair<gbwt::size_type, bool> > haplotype_gbwt_ids;

    /// Copy id of transcript path 
    uint32_t copy_id;

    /// Is it a reference and/or haplotype-specific transcript path.
    bool is_reference;
    bool is_haplotype;

    TranscriptPath() {}

    TranscriptPath(const string & transcript_name, const string & embedded_path_name, const bool is_reference_in, const bool is_haplotype_in) : is_reference(is_reference_in), is_haplotype(is_haplotype_in) {

        assert(!transcript_name.empty());
        assert(!embedded_path_name.empty());

        assert(is_reference != is_haplotype);

        transcript_names.emplace_back(transcript_name);
        embedded_path_names.emplace_back(embedded_path_name, is_reference);

        copy_id = 1;
    }

    TranscriptPath(const string & transcript_name, const gbwt::size_type & haplotype_gbwt_id, const bool is_reference_in, const bool is_haplotype_in) : is_reference(is_reference_in), is_haplotype(is_haplotype_in) {

        assert(!transcript_name.empty());

        assert(is_reference != is_haplotype);

        transcript_names.emplace_back(transcript_name);
        haplotype_gbwt_ids.emplace_back(haplotype_gbwt_id, is_reference);

        copy_id = 1;
    }

    virtual ~TranscriptPath() {};

    string get_name() const;
};

/**
 * Data structure that defines an edited transcript path.
 */ 
struct EditedTranscriptPath : public TranscriptPath {

    /// Transcript path.
    Path path;

    EditedTranscriptPath(const string & transcript_name, const string & embedded_path_name, const bool is_reference_in, const bool is_haplotype_in) : TranscriptPath(transcript_name, embedded_path_name, is_reference_in, is_haplotype_in) {}
    EditedTranscriptPath(const string & transcript_name, const gbwt::size_type & haplotype_gbwt_id, const bool is_reference_in, const bool is_haplotype_in) : TranscriptPath(transcript_name, haplotype_gbwt_id, is_reference_in, is_haplotype_in) {}

    ~EditedTranscriptPath() {};

    handle_t get_first_node_handle(const HandleGraph & graph) const;

};

/**
 * Data structure that defines a completed transcript path.
 */ 
struct CompletedTranscriptPath : public TranscriptPath {

    /// Transcript path.
    vector<handle_t> path;

    CompletedTranscriptPath(const EditedTranscriptPath & edited_transcript_path);
    CompletedTranscriptPath(const EditedTranscriptPath & edited_transcript_path, const HandleGraph & graph);
    ~CompletedTranscriptPath() {};

    handle_t get_first_node_handle(const HandleGraph & graph) const;
};

struct MappingHash
{
    size_t operator()(const Mapping & mapping) const
    {
        size_t seed = 0;

        spp::hash_combine(seed, mapping.position().node_id());
        spp::hash_combine(seed, mapping.position().offset());
        spp::hash_combine(seed, mapping.position().is_reverse());

        for (auto & edit: mapping.edit()) {

            spp::hash_combine(seed, edit.to_length());
        }

        return seed;
    }
 };

/**
 * Class that defines a transcriptome represented by a set of transcript paths.
 */
class Transcriptome {

    public:
    
        Transcriptome(unique_ptr<MutablePathDeletableHandleGraph>&& graph_in); 

        /// Write progress to stderr.
        bool show_progress = false;

        /// Number of threads used for transcript path construction. 
        int32_t num_threads = 1;

        /// Feature type to parse in the gtf/gff file. Parse all types if empty. 
        string feature_type = "exon";

        /// Attribute tag used to parse the transcript id/name in the gtf/gff file. 
        string transcript_tag = "transcript_id";

        /// Speicifies which paths should be compared when collapsing identical paths. 
        /// Can be no, haplotype or all.
        string path_collapse_type = "haplotype";

        /// Treat a missing path in the transcripts/introns as a data error
        bool error_on_missing_path = true;

        /// Adds splice-junstions from intron BED files to the graph. 
        /// Optionally update haplotype GBWT index with new splice-junctions. 
        /// Returns the number of introns parsed. 
        int32_t add_intron_splice_junctions(vector<istream *> intron_streams, unique_ptr<gbwt::GBWT> & haplotype_index, const bool update_haplotypes);

        /// Adds splice-junstions from transcript gtf/gff3 files to the graph and
        /// creates reference transcript paths. Optionally update haplotype GBWT 
        /// index with new splice-junctions. Returns the number of transcripts parsed. 
        int32_t add_reference_transcripts(vector<istream *> transcript_streams, unique_ptr<gbwt::GBWT> & haplotype_index, const bool use_haplotype_paths, const bool update_haplotypes);

        /// Adds haplotype-specific transcript paths by projecting transcripts in
        /// gtf/gff3 files onto either non-reference embedded paths and/or haplotypes
        /// in a GBWT index. Returns the number of haplotype transcript paths projected.   
        int32_t add_haplotype_transcripts(vector<istream *> transcript_streams, const gbwt::GBWT & haplotype_index, const bool proj_emded_paths);

        /// Returns transcript paths.
        const vector<CompletedTranscriptPath> & transcript_paths() const;

        /// Returns the reference transcript paths.
        vector<CompletedTranscriptPath> reference_transcript_paths() const;

        /// Returns the haplotype transcript paths.
        vector<CompletedTranscriptPath> haplotype_transcript_paths() const;

        /// Returns the graph.
        const MutablePathDeletableHandleGraph & graph() const; 

        /// Removes non-transcribed (not in transcript paths) nodes.
        void remove_non_transcribed_nodes();

        /// Chop nodes so that they are not longer than the supplied 
        /// maximum node length. Returns number of chopped nodes.
        void chop_nodes(const uint32_t max_node_length);

        /// Topological sorts graph and compacts node ids. Only works for 
        /// graphs in the PackedGraph format. Return false if not sorted.
        bool sort_compact_nodes();

        /// Embeds transcriptome transcript paths in the graph.  
        /// Returns the number of paths embedded.
        void embed_transcript_paths(const bool add_reference_transcripts, const bool add_haplotype_transcripts);

        /// Adds transcriptome transcript paths as threads to a GBWT index.
        /// Returns the number of added threads.
        void add_transcripts_to_gbwt(gbwt::GBWTBuilder * gbwt_builder, const bool add_bidirectional, const bool exclude_reference_transcripts) const;

        /// Writes transcriptome transcript path sequences to a fasta file.  
        /// Returns the number of written sequences.
        void write_transcript_sequences(ostream * fasta_ostream, const bool exclude_reference_transcripts) const;

        /// Writes info on transcriptome transcript paths to tsv file.
        /// Returns the number of written transcripts.
        void write_transcript_info(ostream * tsv_ostream, const gbwt::GBWT & haplotype_index, const bool exclude_reference_transcripts) const;

        /// Writes the graph to a file.
        void write_graph(ostream * graph_ostream) const;
    
    private:

        /// Transcript paths representing the transcriptome. 
        vector<CompletedTranscriptPath> _transcript_paths;
        mutex mutex_transcript_paths;

        /// Spliced pangenome graph.
        unique_ptr<MutablePathDeletableHandleGraph> _graph;
        mutex mutex_graph;

        /// Parse BED file of introns.
        void parse_introns(vector<Transcript> * introns, istream * intron_stream, const bdsg::PositionOverlay & graph_path_pos_overlay) const;

        /// Parse gtf/gff3 file of transcripts. Returns the number of non-header lines in the parsed file.
        int32_t parse_transcripts(vector<Transcript> * transcripts, uint32_t * number_of_excluded_transcripts, istream * transcript_stream, const bdsg::PositionOverlay & graph_path_pos_overlay, const gbwt::GBWT & haplotype_index, const bool use_haplotype_paths) const;

        /// Returns gbwt path name without phaseblock and subrange. 
        string get_base_gbwt_path_name(const gbwt::GBWT & haplotype_index, const size_t path_id, const unordered_set<string> & gbwt_reference_samples) const;

        /// Parse gtf/gff3 attribute value.
        string parse_attribute_value(const string & attribute, const string & name) const;

        /// Returns the mean node length of the graph
        float mean_node_length() const;

        /// Adds the exon coordinates to a transcript.
        void add_exon(Transcript * transcript, const pair<int32_t, int32_t> & exon_pos) const;

        /// Adds the exon coordinates to a transcript and finds the 
        /// position of each end of a exon on the contig path in the graph.
        void add_exon(Transcript * transcript, const pair<int32_t, int32_t> & exon_pos, const bdsg::PositionOverlay & graph_path_pos_overlay) const;

        /// Reverses exon order if the transcript is on the reverse strand and the exons 
        /// are ordered in reverse.
        void reorder_exons(Transcript * transcript) const;

        /// Checks whether any adjacent exons overlap.
        bool has_overlapping_exons(const vector<Exon> & exons) const;

        /// Constructs edited reference transcript paths from a set of 
        /// transcripts using embedded graph paths.
        list<EditedTranscriptPath> construct_reference_transcript_paths_embedded(const vector<Transcript> & transcripts, const bdsg::PositionOverlay & graph_path_pos_overlay) const;

        /// Threaded reference transcript path construction using embedded paths.
        void construct_reference_transcript_paths_embedded_callback(list<EditedTranscriptPath> * edited_transcript_paths, spp::sparse_hash_map<handle_t, vector<EditedTranscriptPath *> > * edited_transcript_paths_index, mutex * edited_transcript_paths_mutex, const int32_t thread_idx, const vector<Transcript> & transcripts, const bdsg::PositionOverlay & graph_path_pos_overlay) const;

        /// Projects transcripts onto embedded paths in a graph and returns the resulting transcript paths.
        list<EditedTranscriptPath> project_transcript_embedded(const Transcript & cur_transcript, const bdsg::PositionOverlay & graph_path_pos_overlay, const bool use_reference_paths, const bool use_haplotype_paths) const;

        /// Constructs edited reference transcript paths from a set of 
        /// transcripts using haplotype paths in a GBWT index.
        list<EditedTranscriptPath> construct_reference_transcript_paths_gbwt(const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index) const;

        /// Threaded reference transcript path construction using GBWT haplotype paths.
        void construct_reference_transcript_paths_gbwt_callback(list<EditedTranscriptPath> * edited_transcript_paths, spp::sparse_hash_map<handle_t, vector<EditedTranscriptPath *> > * edited_transcript_paths_index, uint32_t * excluded_transcripts, mutex * edited_transcript_paths_mutex, const int32_t thread_idx, const vector<pair<uint32_t, uint32_t> > & chrom_transcript_sets, const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const spp::sparse_hash_map<string, map<uint32_t, uint32_t> > & haplotype_name_index) const;

        /// Constructs haplotype transcript paths by projecting transcripts onto
        /// embedded paths in a graph and/or haplotypes in a GBWT index. 
        /// Adds haplotype transcript to transcriptome.
        void project_haplotype_transcripts(const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const bdsg::PositionOverlay & graph_path_pos_overlay, const bool proj_emded_paths, const float mean_node_length);

        /// Threaded haplotype transcript projecting.
        void project_haplotype_transcripts_callback(list<CompletedTranscriptPath> * completed_transcript_paths, spp::sparse_hash_map<handle_t, vector<CompletedTranscriptPath *> > * completed_transcript_paths_index,  mutex * completed_transcript_paths_mutex, const int32_t thread_idx, const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const bdsg::PositionOverlay & graph_path_pos_overlay, const bool proj_emded_paths, const float mean_node_length);

        /// Projects transcripts onto haplotypes in a GBWT index and returns the resulting transcript paths.
        list<EditedTranscriptPath> project_transcript_gbwt(const Transcript & cur_transcript, const gbwt::GBWT & haplotype_index,
                                                           const unordered_set<string>& reference_samples, const float mean_node_length) const;

        /// Extracts all unique haplotype paths between two nodes from a GBWT index and returns the 
        /// resulting paths and the corresponding haplotype ids for each path.
        vector<pair<exon_nodes_t, thread_ids_t> > get_exon_haplotypes(const vg::id_t start_node, const vg::id_t end_node, const gbwt::GBWT & haplotype_index,  const unordered_set<string>& reference_samples, const int32_t expected_length) const;

        /// Remove redundant transcript paths and update index.
        template <class T>
        void remove_redundant_transcript_paths(list<T> * new_transcript_paths, spp::sparse_hash_map<handle_t, vector<T*> > * transcript_paths_index) const;

        /// Constructs completed transcripts paths from 
        /// edited transcript paths. Checks that the
        /// paths contain no edits compared to the graph.
        list<CompletedTranscriptPath> construct_completed_transcript_paths(const list<EditedTranscriptPath> & edited_transcript_paths) const;

        /// Adds edited transcript paths to transcriptome
        /// Checks that the paths contain no edits compared to the graph.
        void add_edited_transcript_paths(const list<EditedTranscriptPath> & edited_transcript_paths);

        /// Checks whether transcript path only consist of 
        /// whole nodes (complete). 
        bool has_novel_exon_boundaries(const list<EditedTranscriptPath> & edited_transcript_paths, const bool include_transcript_ends) const;

        /// Augments the graph with transcript path exon boundaries and 
        /// splice-junctions. Updates threads in gbwt index to match the augmented graph. 
        /// Optinally adds transcript paths to the transcriptome.
        void augment_graph(const list<EditedTranscriptPath> & edited_transcript_paths, const bool is_introns, unique_ptr<gbwt::GBWT> & haplotype_index, const bool update_haplotypes, const bool add_reference_transcript_paths);

        /// Update threads in gbwt index using graph translations. 
        void update_haplotype_index(unique_ptr<gbwt::GBWT> & haplotype_index, const spp::sparse_hash_map<handle_t, vector<pair<int32_t, handle_t> > > & update_index) const;

        /// Update/split node handles in transcriptome transcript paths according to index.
        void update_transcript_paths(const spp::sparse_hash_map<handle_t, vector<pair<int32_t, handle_t> > > & update_index);

        /// Adds transcript path splice-junction edges to the graph
        void add_splice_junction_edges(const list<EditedTranscriptPath> & edited_transcript_paths);
        void add_splice_junction_edges(const list<CompletedTranscriptPath> & completed_transcript_paths);
        void add_splice_junction_edges(const vector<CompletedTranscriptPath> & completed_transcript_paths);

        /// Collects all unique nodes in transcriptome transcript paths.
        void collect_transcribed_nodes(spp::sparse_hash_set<nid_t> * transcribed_nodes) const;

        /// Sort transcriptome transcript paths by name and is reference,
        /// and update their copy ids.
        void sort_transcript_paths_update_copy_id();

};

}


#endif

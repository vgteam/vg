
#ifndef VG_TRANSCRIPTOME_HPP_INCLUDED
#define VG_TRANSCRIPTOME_HPP_INCLUDED

/** \file transcriptome.hpp
 * 
 */

#include <algorithm>
#include <mutex>

#include <google/protobuf/util/message_differencer.h>
#include <gbwt/dynamic_gbwt.h>

#include "../vg.hpp"
#include "../path_index.hpp"
#include "../types.hpp"

namespace vg {

using namespace std;

typedef vector<gbwt::node_type> exon_nodes_t;

/**
 * 
 */
class Transcriptome {

    public:

        Transcriptome() {};   

        int32_t num_threads = 1;
        string transcript_tag = "transcript_id";
        bool use_embedded_paths = false;
        bool collapse_transcript_paths = true;
        bool filter_reference_transcript_paths = false;

        struct TranscriptPath {

            Path path;

            int32_t num_total;           
            int32_t num_reference;

            TranscriptPath(const bool is_reference) : num_total(1), num_reference(is_reference) {}
        };

        void add_transcripts(istream & transcript_stream, VG & graph, const gbwt::GBWT & haplotype_index);
        const vector<Path> & transcript_paths() const;
        int32_t size() const;
        void edit_graph(VG * graph, const bool add_paths);

        void construct_gbwt(gbwt::GBWTBuilder * gbwt_builder) const;
        void write_gam_alignments(ostream * gam_ostream) const;
        void write_fasta_sequences(ostream * fasta_ostream, VG & graph) const;
   
    private:

        vector<Path> _transcriptome;
        mutex trancriptome_mutex;

        struct Transcript {

            const string name;
            const bool is_reverse;
            const string chrom;
            
            vector<pair<int32_t, int32_t> > exons;
            vector<pair<Position, Position> > exon_nodes;

            Transcript(const string & name_in, const bool is_reverse_in, const string & chrom_in) : name(name_in), is_reverse(is_reverse_in), chrom(chrom_in) {}
        };

        void add_exon(Transcript * transcript, const pair<int32_t, int32_t> & exon_pos, const PathIndex & chrom_path_index) const;
        void reorder_exons(Transcript * transcript) const;
        void project_transcripts(const vector<Transcript> & transcripts, VG & graph, const gbwt::GBWT & haplotype_index, const float mean_node_length);
        void project_transcripts_callback(const int32_t thread_idx, const vector<Transcript> & transcripts, VG & graph, const gbwt::GBWT & haplotype_index, const float mean_node_length);
        list<TranscriptPath> project_transcript_gbwt(const Transcript & cur_transcript, VG & graph, const gbwt::GBWT & haplotype_index, const float mean_node_length) const;
        vector<pair<exon_nodes_t, vector<gbwt::size_type> > > get_exon_haplotypes(const vg::id_t start_node, const vg::id_t end_node, const gbwt::GBWT & haplotype_index, const int32_t expected_length) const;
        list<TranscriptPath> project_transcript_embedded(const Transcript & cur_transcript, VG & graph) const;
        void collapse_identical_paths(list<TranscriptPath> * cur_transcript_paths) const;
};
}

#endif

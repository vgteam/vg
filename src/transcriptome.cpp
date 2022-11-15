
#include <thread>

#include <gbwtgraph/utils.h>

#include "../io/save_handle_graph.hpp"

#include "transcriptome.hpp"
#include "../augment.hpp"
#include "../utility.hpp"

namespace vg {

using namespace vg::io;

using namespace std;

#define transcriptome_debug

bool operator==(const Exon & lhs, const Exon & rhs) { 

    return (lhs.coordinates == rhs.coordinates && 
            lhs.border_offsets == rhs.border_offsets && 
            lhs.border_steps == rhs.border_steps);
}

bool operator!=(const Exon & lhs, const Exon & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const Exon & lhs, const Exon & rhs) { 

    if (lhs.coordinates.first != rhs.coordinates.first) {

        return (lhs.coordinates.first < rhs.coordinates.first);    
    } 

    if (lhs.coordinates.second != rhs.coordinates.second) {

        return (lhs.coordinates.second < rhs.coordinates.second);    
    } 

    return false;
}

bool operator==(const Transcript & lhs, const Transcript & rhs) { 

    return (lhs.name == rhs.name &&
            lhs.is_reverse == rhs.is_reverse &&
            lhs.chrom == rhs.chrom && 
            lhs.chrom_length == rhs.chrom_length && 
            lhs.exons == rhs.exons);
}

bool operator!=(const Transcript & lhs, const Transcript & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const Transcript & lhs, const Transcript & rhs) { 

    if (lhs.chrom != rhs.chrom) {

        return (lhs.chrom < rhs.chrom);    
    } 

    if (!lhs.exons.empty() && !rhs.exons.empty()) {

        if (lhs.exons.front() != rhs.exons.front()) {

            return (lhs.exons.front() < rhs.exons.front());    
        } 
    
    } else if (!lhs.exons.empty() || !rhs.exons.empty()) {

        return (lhs.exons.size() < rhs.exons.size());    
    }

    if (lhs.is_reverse != rhs.is_reverse) {

        return (lhs.is_reverse < rhs.is_reverse);    
    } 

    return false;
}

bool operator==(const Mapping & lhs, const Mapping & rhs) { 

    return google::protobuf::util::MessageDifferencer::Equals(lhs, rhs);
}

bool operator!=(const Mapping & lhs, const Mapping & rhs) { 

    return !(lhs == rhs);
}

bool operator==(const Path & lhs, const Path & rhs) { 

    return google::protobuf::util::MessageDifferencer::Equals(lhs, rhs);
}

bool operator!=(const Path & lhs, const Path & rhs) { 

    return !(lhs == rhs);
}

bool sort_pair_by_second(const pair<uint32_t, uint32_t> & lhs, const pair<uint32_t, uint32_t> & rhs) {
    
    return (lhs.second < rhs.second);
}

handle_t mapping_to_handle(const Mapping & mapping, const HandleGraph & graph) {

    return (graph.get_handle(mapping.position().node_id(), mapping.position().is_reverse()));
}

string TranscriptPath::get_name() const {

    assert(!transcript_names.empty());
    assert(!transcript_names.front().empty());

    return (transcript_names.front() + "_HST" + to_string(copy_id));
}

handle_t EditedTranscriptPath::get_first_node_handle(const HandleGraph & graph) const {

    assert(path.mapping_size() > 0);

    return graph.get_handle(path.mapping()[0].position().node_id(), path.mapping()[0].position().is_reverse());
}

CompletedTranscriptPath::CompletedTranscriptPath(const EditedTranscriptPath & edited_transcript_path_in) {

    transcript_names = edited_transcript_path_in.transcript_names;
    embedded_path_names = edited_transcript_path_in.embedded_path_names;
    haplotype_gbwt_ids = edited_transcript_path_in.haplotype_gbwt_ids;

    copy_id = edited_transcript_path_in.copy_id;

    is_reference = edited_transcript_path_in.is_reference;
    is_haplotype = edited_transcript_path_in.is_haplotype;
}

CompletedTranscriptPath::CompletedTranscriptPath(const EditedTranscriptPath & edited_transcript_path_in, const HandleGraph & graph) {

    transcript_names = edited_transcript_path_in.transcript_names;
    embedded_path_names = edited_transcript_path_in.embedded_path_names;
    haplotype_gbwt_ids = edited_transcript_path_in.haplotype_gbwt_ids;

    copy_id = edited_transcript_path_in.copy_id;

    is_reference = edited_transcript_path_in.is_reference;
    is_haplotype = edited_transcript_path_in.is_haplotype;

    path.reserve(edited_transcript_path_in.path.mapping_size());

    for (auto mapping: edited_transcript_path_in.path.mapping()) {

        auto handle = mapping_to_handle(mapping, graph);

        // Check that the path only consist of whole nodes (complete).        
        assert(mapping.edit_size() == 1);
        assert(edit_is_match(mapping.edit(0)));
        assert(mapping.position().offset() == 0);
        assert(mapping.edit(0).from_length() == graph.get_length(handle));

        path.emplace_back(handle);
    }
}

handle_t CompletedTranscriptPath::get_first_node_handle(const HandleGraph & graph) const {

    return path.front();
}

Transcriptome::Transcriptome(unique_ptr<MutablePathDeletableHandleGraph>&& graph_in) : _graph(move(graph_in)) {
    
    if (!_graph) {
        cerr << "\tERROR: Could not load graph." << endl;
        exit(1);
    }
}

int32_t Transcriptome::add_intron_splice_junctions(vector<istream *> intron_streams, unique_ptr<gbwt::GBWT> & haplotype_index, const bool update_haplotypes) {

#ifdef transcriptome_debug
    double time_parsing_1 = gcsa::readTimer();
    cerr << "\tDEBUG Parsing start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    // Create path position overlay of graph
    bdsg::PositionOverlay graph_path_pos_overlay(_graph.get());

    vector<Transcript> introns;

    for (auto & intron_stream: intron_streams) {

        // Parse introns in BED format.
        parse_introns(&introns, intron_stream, graph_path_pos_overlay);
    }

    if (introns.empty()) {

        cerr << "\tERROR: No intron parsed" << endl;
        exit(1);        
    }

    sort(introns.begin(), introns.end());

    cerr << "\tParsed " << introns.size() << " introns" << endl;

#ifdef transcriptome_debug
    cerr << "\tDEBUG: " << gcsa::readTimer() - time_parsing_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_project_1 = gcsa::readTimer();
    cerr << "\tDEBUG Construction start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    // Construct edited reference intron paths using embedded graph paths.
    auto edited_transcript_paths = construct_reference_transcript_paths_embedded(introns, graph_path_pos_overlay);

    cerr << "\tConstructed " << edited_transcript_paths.size() << " intron paths" << endl;

#ifdef transcriptome_debug
    cerr << "\tDEBUG: " << gcsa::readTimer() - time_project_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_augment_1 = gcsa::readTimer();
    cerr << "\tDEBUG Updating start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    if (has_novel_exon_boundaries(edited_transcript_paths, false)) {

        // Augment graph with new exon boundaries and splice-junction edges. 
        augment_graph(edited_transcript_paths, false, haplotype_index, update_haplotypes, false);
    
    } else {

        // Augment graph with new splice-junction edges. 
        add_splice_junction_edges(edited_transcript_paths);
    }

    assert(_transcript_paths.empty());

    cerr << "\tUpdated graph with intron paths" << endl;

#ifdef transcriptome_debug
    cerr << "\tDEBUG: " << gcsa::readTimer() - time_augment_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    return introns.size();
}

int32_t Transcriptome::add_reference_transcripts(vector<istream *> transcript_streams, unique_ptr<gbwt::GBWT> & haplotype_index, const bool use_haplotype_paths, const bool update_haplotypes) {

#ifdef transcriptome_debug
    double time_parsing_1 = gcsa::readTimer();
    cerr << "\tDEBUG Parsing start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    bdsg::PositionOverlay graph_path_pos_overlay;

    if (!use_haplotype_paths) {

        // Create path position overlay of graph if using embedded path references.        
        graph_path_pos_overlay = bdsg::PositionOverlay(_graph.get());  
    }

    vector<Transcript> transcripts;

    int64_t lines_parsed = 0;
    for (auto & transcript_stream: transcript_streams) {

        // Parse transcripts in gtf/gff3 format.
        lines_parsed += parse_transcripts(&transcripts, transcript_stream, graph_path_pos_overlay, *haplotype_index, use_haplotype_paths);
    }

    if (transcripts.empty() && lines_parsed != 0) {

        cerr << "\tERROR: No transcripts parsed (remember to set feature type \"-y\" in vg rna or \"-f\" in vg autoindex)" << endl;
        exit(1);        
    }

    sort(transcripts.begin(), transcripts.end());

    cerr << "\tParsed " << transcripts.size() << " transcripts" << endl;

#ifdef transcriptome_debug
    cerr << "\tDEBUG: " << gcsa::readTimer() - time_parsing_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_project_1 = gcsa::readTimer();
    cerr << "\tDEBUG Construction start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    list<EditedTranscriptPath> edited_transcript_paths;

    if (use_haplotype_paths) {

        // Construct edited reference transcript paths using haplotype GBWT paths.
        edited_transcript_paths = construct_reference_transcript_paths_gbwt(transcripts, *haplotype_index);

    } else {

        // Construct edited reference transcript paths using embedded graph paths.
        edited_transcript_paths = construct_reference_transcript_paths_embedded(transcripts, graph_path_pos_overlay);
    } 

    cerr << "\tConstructed " << edited_transcript_paths.size() << " reference transcript paths" << endl;

#ifdef transcriptome_debug
    cerr << "\tDEBUG: " << gcsa::readTimer() - time_project_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_augment_1 = gcsa::readTimer();
    cerr << "\tDEBUG Updating start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    if (has_novel_exon_boundaries(edited_transcript_paths, true)) {

        // Augment graph with new exon boundaries and splice-junction edges. 
        // Adds the edited transcript paths as reference transcript paths.
        augment_graph(edited_transcript_paths, true, haplotype_index, update_haplotypes, true);
    
    } else {

        // Augment graph with new splice-junction edges and add reference transcript paths.
        add_edited_transcript_paths(edited_transcript_paths);
    }

    // Update copy id of transcript paths.
    update_copy_id();

    cerr << "\tUpdated graph with reference transcript paths" << endl;

#ifdef transcriptome_debug
    cerr << "\tDEBUG: " << gcsa::readTimer() - time_augment_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    return transcripts.size();
}

int32_t Transcriptome::add_haplotype_transcripts(vector<istream *> transcript_streams, const gbwt::GBWT & haplotype_index, const bool proj_emded_paths) {

#ifdef transcriptome_debug
    double time_parsing_1 = gcsa::readTimer();
    cerr << "\tDEBUG Parsing start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    // Create path position overlay of splice graph
    bdsg::PositionOverlay graph_path_pos_overlay(_graph.get());

    vector<Transcript> transcripts;

    int64_t lines_parsed = 0;
    for (auto & transcript_stream: transcript_streams) {

        // Parse transcripts in gtf/gff3 format.
        lines_parsed += parse_transcripts(&transcripts, transcript_stream, graph_path_pos_overlay, haplotype_index, false);
    }

    if (transcripts.empty() && lines_parsed != 0) {

        cerr << "\tERROR: No transcripts parsed (remember to set feature type \"-y\" in vg rna or \"-f\" in vg autoindex)" << endl;
        exit(1);        
    }
    
    sort(transcripts.begin(), transcripts.end());

    cerr << "\tParsed " << transcripts.size() << " transcripts" << endl;

#ifdef transcriptome_debug
    cerr << "\tDEBUG: " << gcsa::readTimer() - time_parsing_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_project_1 = gcsa::readTimer();
    cerr << "\tDEBUG Projection start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    // Save number of transcript paths before adding new.
    auto pre_num_transcript_paths = _transcript_paths.size();

    // Project and add transcripts to transcriptome.
    project_haplotype_transcripts(transcripts, haplotype_index, graph_path_pos_overlay, proj_emded_paths, mean_node_length());

    // Augment splice graph with new splice-junction edges.    
    add_splice_junction_edges(_transcript_paths);

    // Update copy id of transcript paths.
    update_copy_id();

    assert(_transcript_paths.size() >= pre_num_transcript_paths);

    cerr << "\tProjected " << _transcript_paths.size() - pre_num_transcript_paths << " haplotype-specific transcript paths" << endl;

#ifdef transcriptome_debug
    cerr << "\tDEBUG: " << gcsa::readTimer() - time_project_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    return (_transcript_paths.size() - pre_num_transcript_paths);
}

void Transcriptome::parse_introns(vector<Transcript> * introns, istream * intron_stream, const bdsg::PositionOverlay & graph_path_pos_overlay) const {

    int32_t line_number = 0;

    string chrom;
    string pos;
    string end;
    string strand;

    while (intron_stream->good()) {

        line_number += 1;
        getline(*intron_stream, chrom, '\t');

        // Skip header.
        if (chrom.empty() || chrom.front() == '#') {

            intron_stream->ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }

        assert(_graph->has_path(chrom) == graph_path_pos_overlay.has_path(chrom));

        if (!_graph->has_path(chrom)) {

            if (error_on_missing_path) {
                cerr << "\tERROR: Chromomsome path \"" << chrom << "\" not found in graph (line " << line_number << ")." << endl;
                exit(1);
            }
            else {
                // seek to the end of the line
                intron_stream->ignore(numeric_limits<streamsize>::max(), '\n');
                continue;
            }
        }

        // Parse start and end intron position and convert end to inclusive.
        assert(getline(*intron_stream, pos, '\t'));
        int32_t spos = stoi(pos);
        
        assert(getline(*intron_stream, end));
        auto end_ss = stringstream(end);
        assert(getline(end_ss, pos, '\t'));
        int32_t epos = stoi(pos) - 1;

        assert(spos <= epos);

        getline(end_ss, strand, '\t');
        getline(end_ss, strand, '\t');

        bool is_reverse = false;

        if (getline(end_ss, strand, '\t')) {

            assert(strand == "+" || strand == "-");
            is_reverse = (strand == "-") ? true : false;
        }

        // Create "intron" transcript.
        introns->emplace_back(Transcript("", is_reverse, chrom, graph_path_pos_overlay.get_path_length(_graph->get_path_handle(chrom))));

        // Add intron boundaries as flanking exons to current "intron" transcript.
        add_exon(&(introns->back()), make_pair(spos - 1, spos - 1), graph_path_pos_overlay);
        add_exon(&(introns->back()), make_pair(epos + 1, epos + 1), graph_path_pos_overlay);
    }
}

int32_t Transcriptome::parse_transcripts(vector<Transcript> * transcripts, istream * transcript_stream, const bdsg::PositionOverlay & graph_path_pos_overlay, const gbwt::GBWT & haplotype_index, const bool use_haplotype_paths) const {

    spp::sparse_hash_map<string, uint32_t> chrom_lengths;

    if (use_haplotype_paths) {

        assert(haplotype_index.bidirectional());
        assert(haplotype_index.hasMetadata());
        
        assert(haplotype_index.metadata.hasPathNames());
        assert(haplotype_index.metadata.hasContigNames());
        
        for (size_t i = 0; i < haplotype_index.sequences(); i++) {

            // Skip reverse threads in bidirectional gbwt index.
            if (i % 2 == 1) {

                continue;
            }

            chrom_lengths.emplace(haplotype_index.metadata.contig(haplotype_index.metadata.path(gbwt::Path::id(i)).contig), numeric_limits<uint32_t>::max());
        }      

    } else {

        assert(_graph->for_each_path_handle([&](const path_handle_t & path_handle) {

            assert(graph_path_pos_overlay.has_path(_graph->get_path_name(path_handle)));
            chrom_lengths.emplace(_graph->get_path_name(path_handle), graph_path_pos_overlay.get_path_length(path_handle));
        }));
    }

    spp::sparse_hash_map<string, uint32_t> transcripts_index;

    int32_t line_number = 0;
    int32_t parsed_lines = 0;

    string chrom;
    string feature;
    string pos;
    string strand;
    string attributes;
    string attribute;

    bool zero_based_exon_number = false;

    while (transcript_stream->good()) {

        line_number += 1;
        getline(*transcript_stream, chrom, '\t');

        // Skip header.
        if (chrom.empty() || chrom.front() == '#') {

            transcript_stream->ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }
        parsed_lines += 1;

        auto chrom_lengths_it = chrom_lengths.find(chrom);

        if (chrom_lengths_it == chrom_lengths.end()) {

            if (error_on_missing_path) {

                cerr << "\tERROR: Chromomsome path \"" << chrom << "\" not found in graph or haplotypes index (line " << line_number << ")." << endl;
                exit(1);
            
            } else {

                // Seek to the end of the line.
                transcript_stream->ignore(numeric_limits<streamsize>::max(), '\n');
                continue;
            }
        }

        transcript_stream->ignore(numeric_limits<streamsize>::max(), '\t');         
        assert(getline(*transcript_stream, feature, '\t'));

        // Select only relevant feature types.
        if (feature != feature_type && !feature_type.empty()) {

            transcript_stream->ignore(numeric_limits<streamsize>::max(), '\n');  
            continue;
        }

        // Parse start and end exon position and convert to 0-base.
        assert(getline(*transcript_stream, pos, '\t'));
        int32_t spos = stoi(pos) - 1;
        assert(getline(*transcript_stream, pos, '\t'));
        int32_t epos = stoi(pos) - 1;

        assert(spos <= epos);

        // Skip score column.
        transcript_stream->ignore(numeric_limits<streamsize>::max(), '\t');  
        
        // Parse strand and set whether it is reverse.
        assert(getline(*transcript_stream, strand, '\t'));
        assert(strand == "+" || strand == "-");
        bool is_reverse = (strand == "-") ? true : false;

        // Skip frame column.
        transcript_stream->ignore(numeric_limits<streamsize>::max(), '\t');  

        string transcript_id = "";
        int32_t exon_number = -1;

        assert(getline(*transcript_stream, attributes, '\n'));
        stringstream attributes_ss(attributes);    

        while (getline(attributes_ss, attribute, ';')) {

            if (attribute.empty()) {

                break;
            }

            // Parse transcript ID.
            if (transcript_id.empty()) {

                transcript_id = parse_attribute_value(attribute, transcript_tag);
            }

            // Parse exon number.
            if (exon_number < 0) {

                auto exon_number_str = parse_attribute_value(attribute, "exon_number");

                if (exon_number_str.empty()) {

                    // If not exon_number attribute try ID.
                    auto exon_id = parse_attribute_value(attribute, "ID");

                    if (count(exon_id.begin(), exon_id.end(), ':') == 2) {

                        auto exon_id_ss = stringstream(exon_id);

                        string element;
                        getline(exon_id_ss, element, ':');   
                        
                        if (element == "exon") {

                            getline(exon_id_ss, element, ':');
                            getline(exon_id_ss, element);   

                            exon_number = stoi(element);
                        }  
                    }   

                } else {

                    exon_number = stoi(exon_number_str);
                }
            }

            if (!transcript_id.empty() && exon_number >= 0) {

                break;
            }
        }

        if (transcript_id.empty()) {

            cerr << "\tERROR: Tag \"" << transcript_tag << "\" not found in attributes (line " << line_number << ")." << endl;
            exit(1);
        }

        auto transcripts_index_it = transcripts_index.emplace(transcript_id, transcripts->size());

        // Is this a new transcript.
        if (transcripts_index_it.second) {

            transcripts->emplace_back(Transcript(transcript_id, is_reverse, chrom, chrom_lengths_it->second));
        }

        Transcript * transcript = &(transcripts->at(transcripts_index_it.first->second));

        assert(transcript->name == transcript_id);
        assert(transcript->is_reverse == is_reverse);
        assert(transcript->chrom == chrom);
        assert(transcript->chrom_length == chrom_lengths_it->second);

        if (use_haplotype_paths) {
            
            // Add exon to current transcript.
            add_exon(transcript, make_pair(spos, epos));

        } else {

            // Add exon to current transcript.
            add_exon(transcript, make_pair(spos, epos), graph_path_pos_overlay);
        }

        // Check if exons are in correct order in file. 
        if (exon_number >= 0) {

            // If first transcript and exon, set whether exon numbering is zero-based. 
            if (transcripts_index.size() == 1 && transcript->exons.size() == 1) {

                zero_based_exon_number = (exon_number == 0) ? true : false;
            }

            if (transcript->exons.size() - static_cast<uint32_t>(zero_based_exon_number) != exon_number) {

                cerr << "\tERROR: Exons are not in correct order according to attributes (line " << line_number << ")." << endl;
                exit(1); 
            } 
        }
    }

    for (auto & transcript_idx: transcripts_index) {

        // Reorder reversed order exons.
        reorder_exons(&(transcripts->at(transcript_idx.second)));
    }
    
    return parsed_lines;
}

string Transcriptome::parse_attribute_value(const string & attribute, const string & name) const {

    string value = "";

    const uint32_t attribute_start_pos = (attribute.front() == ' ');

    if (attribute.substr(attribute_start_pos, name.size()) == name) {

        // Is gff3 format.
        if (attribute.substr(name.size(), 1) == "=") {

            assert(attribute_start_pos == 0);
            value = attribute.substr(name.size() + 1);

        } else {

            // Is value in quotes (""). 
            if (attribute.substr(attribute_start_pos + name.size() + 1, 1) == "\"") {

                value = attribute.substr(attribute_start_pos + name.size() + 2);

                assert(value.back() == '\"');
                value.pop_back();

            } else {

                value = attribute.substr(attribute_start_pos + name.size() + 1);
            }
        }
    }

    return value;
}

float Transcriptome::mean_node_length() const {

    return static_cast<float>(_graph->get_total_length()) / _graph->get_node_count();
}

void Transcriptome::add_exon(Transcript * transcript, const pair<int32_t, int32_t> & exon_pos) const {

    assert(exon_pos.first >= 0);
    assert(exon_pos.second < transcript->chrom_length);

    transcript->exons.emplace_back(Exon());
    transcript->exons.back().coordinates = exon_pos;
}

void Transcriptome::add_exon(Transcript * transcript, const pair<int32_t, int32_t> & exon_pos, const bdsg::PositionOverlay & graph_path_pos_overlay) const {

    add_exon(transcript, exon_pos);

    // Exon border positions (last position in upstream intron and 
    // first position in downstream intron). The positions are not 
    // offset if it is the first or last on the path.
    const pair<uint32_t, uint32_t> exon_border_pos = make_pair(max(0, exon_pos.first - 1), min(static_cast<int32_t>(transcript->chrom_length) - 1, exon_pos.second + 1));

    auto path_handle = _graph->get_path_handle(transcript->chrom);

    // Find path positions of exon node borders (start - 1 and end + 1).
    auto chrom_path_start_step = graph_path_pos_overlay.get_step_at_position(path_handle, exon_border_pos.first);
    auto chrom_path_end_step = graph_path_pos_overlay.get_step_at_position(path_handle, exon_border_pos.second);

    assert(chrom_path_start_step != graph_path_pos_overlay.path_end(path_handle));
    assert(chrom_path_end_step != graph_path_pos_overlay.path_end(path_handle));

    // Find the start position of the exon border nodes.
    auto chrom_path_start_node_pos = graph_path_pos_overlay.get_position_of_step(chrom_path_start_step);
    auto chrom_path_end_node_pos = graph_path_pos_overlay.get_position_of_step(chrom_path_end_step);

    assert(chrom_path_start_node_pos <= exon_border_pos.first);
    assert(chrom_path_end_node_pos <= exon_border_pos.second);

    // Add node offsets of exon border boundaries. 
    transcript->exons.back().border_offsets = make_pair(exon_border_pos.first - chrom_path_start_node_pos, exon_border_pos.second  - chrom_path_end_node_pos);

    // Add path steps of exon border boundaries. 
    transcript->exons.back().border_steps = make_pair(chrom_path_start_step, chrom_path_end_step);
}

void Transcriptome::reorder_exons(Transcript * transcript) const {

    if (transcript->is_reverse) {

        // Is exons in reverse order.
        bool is_reverse_order = true;
        for (size_t i = 1; i < transcript->exons.size(); i++) {

            if (transcript->exons.at(i).coordinates.second > transcript->exons.at(i - 1).coordinates.first) { 

                is_reverse_order = false; 
            }
        }

        // Reverse if exons are in reverse order.
        if (is_reverse_order) { 

            reverse(transcript->exons.begin(), transcript->exons.end()); 
        }
    }
}

list<EditedTranscriptPath> Transcriptome::construct_reference_transcript_paths_embedded(const vector<Transcript> & transcripts, const bdsg::PositionOverlay & graph_path_pos_overlay) const {

    list<EditedTranscriptPath> edited_transcript_paths;
    spp::sparse_hash_map<handle_t, vector<EditedTranscriptPath *> > edited_transcript_paths_index;

    mutex edited_transcript_paths_mutex;

    vector<thread> construction_threads;
    construction_threads.reserve(num_threads);

    // Spawn construction threads.
    for (size_t thread_idx = 0; thread_idx < num_threads; thread_idx++) {

        construction_threads.push_back(thread(&Transcriptome::construct_reference_transcript_paths_embedded_callback, this, &edited_transcript_paths, &edited_transcript_paths_index, &edited_transcript_paths_mutex, thread_idx, ref(transcripts), ref(graph_path_pos_overlay)));
    }

    // Join construction threads.   
    for (auto & thread: construction_threads) {
        
        thread.join();
    }

    return edited_transcript_paths;
}

void Transcriptome::construct_reference_transcript_paths_embedded_callback(list<EditedTranscriptPath> * edited_transcript_paths, spp::sparse_hash_map<handle_t, vector<EditedTranscriptPath *> > * edited_transcript_paths_index, mutex * edited_transcript_paths_mutex, const int32_t thread_idx, const vector<Transcript> & transcripts, const bdsg::PositionOverlay & graph_path_pos_overlay) const {

    list<EditedTranscriptPath> thread_edited_transcript_paths;

    int32_t transcripts_idx = thread_idx;

    while (transcripts_idx < transcripts.size()) {

        // Get next transcript belonging to current thread.
        const Transcript & transcript = transcripts.at(transcripts_idx);

        // Construct edited transcript paths.
        auto new_edited_transcript_paths = project_transcript_embedded(transcript, graph_path_pos_overlay, true, false);

        if (!new_edited_transcript_paths.empty()) {

            assert(new_edited_transcript_paths.size() == 1);
            thread_edited_transcript_paths.emplace_back(new_edited_transcript_paths.front());
        }

        transcripts_idx += num_threads;
    }

    edited_transcript_paths_mutex->lock();
    remove_redundant_transcript_paths<EditedTranscriptPath>(&thread_edited_transcript_paths, edited_transcript_paths_index);
    edited_transcript_paths->splice(edited_transcript_paths->end(), thread_edited_transcript_paths);
    edited_transcript_paths_mutex->unlock();
}

list<EditedTranscriptPath> Transcriptome::project_transcript_embedded(const Transcript & cur_transcript, const bdsg::PositionOverlay & graph_path_pos_overlay, const bool use_reference_paths, const bool use_haplotype_paths) const {

    assert(use_reference_paths != use_haplotype_paths);

    vector<multimap<path_handle_t, step_handle_t> > exon_start_node_path_steps;
    vector<multimap<path_handle_t, step_handle_t> > exon_end_node_path_steps;

    exon_start_node_path_steps.reserve(cur_transcript.exons.size());
    exon_end_node_path_steps.reserve(cur_transcript.exons.size());

    // Get embedded path ids and node mappings for all exon border nodes in transcript.
    for (auto & exon: cur_transcript.exons) {

        auto exon_path_handle = _graph->get_path_handle_of_step(exon.border_steps.first);
        assert(exon_path_handle == _graph->get_path_handle_of_step(exon.border_steps.second));

        assert(cur_transcript.chrom == _graph->get_path_name(exon_path_handle));

        exon_start_node_path_steps.emplace_back(multimap<path_handle_t, step_handle_t>());
        exon_start_node_path_steps.back().emplace(exon_path_handle, exon.border_steps.first);

        exon_end_node_path_steps.emplace_back(multimap<path_handle_t, step_handle_t>());
        exon_end_node_path_steps.back().emplace(exon_path_handle, exon.border_steps.second);

        if (use_haplotype_paths) {

            auto start_border_is_reverse = _graph->get_is_reverse(_graph->get_handle_of_step(exon.border_steps.first));

            _graph->for_each_step_on_handle(_graph->get_handle_of_step(exon.border_steps.first), [&](const step_handle_t & step) {

                // Do not allow multiple lengths due to cycles for reference exon.
                auto step_path_handle = _graph->get_path_handle_of_step(step);
                if (step_path_handle != exon_path_handle && _graph->get_is_reverse(_graph->get_handle_of_step(step)) == start_border_is_reverse) {

                    exon_start_node_path_steps.back().emplace(step_path_handle, step);
                }
            });

            auto end_border_is_reverse = _graph->get_is_reverse(_graph->get_handle_of_step(exon.border_steps.second));

            _graph->for_each_step_on_handle(_graph->get_handle_of_step(exon.border_steps.second), [&](const step_handle_t & step) {

                // Do not allow multiple lengths due to cycles for reference exon.
                auto step_path_handle = _graph->get_path_handle_of_step(step);
                if (step_path_handle != exon_path_handle && _graph->get_is_reverse(_graph->get_handle_of_step(step)) == end_border_is_reverse) {

                    exon_end_node_path_steps.back().emplace(step_path_handle, step);
                }
            });
        }
    }

    list<EditedTranscriptPath> edited_transcript_paths;

    // Loop over all paths that contain the transcript start node.
    for (auto & path_steps_start: exon_start_node_path_steps.front()) {

        // Skip path if transcript end node is not in the current path.
        if (exon_end_node_path_steps.back().find(path_steps_start.first) == exon_end_node_path_steps.back().end()) {

            continue;
        }

        const auto path_origin_name = _graph->get_path_name(path_steps_start.first);

        // Skip alternative allele paths (_alt).
        if (Paths::is_alt(path_origin_name)) {

            continue;
        }

        // Do not construct transcript paths originating from a reference chromosome/contig.
        if (path_origin_name == cur_transcript.chrom && !use_reference_paths) {

            continue;
        }

        // Do not construct transcript paths originating from a haplotype.
        if (path_origin_name != cur_transcript.chrom && !use_haplotype_paths) {

            continue;
        }

        list<EditedTranscriptPath> cur_edited_transcript_paths;

        // Construct transcript path and set transcript origin name.
        cur_edited_transcript_paths.emplace_back(cur_transcript.name, path_origin_name, use_reference_paths, use_haplotype_paths);

        bool is_partial = false;

        for (size_t exon_idx = 0; exon_idx < exon_start_node_path_steps.size(); ++exon_idx) {

            if (is_partial) { break; }

            // Transcripts with cycles at both exon boundaries are currently 
            // not supported.
            // TODO: Add support for this.
            if (exon_start_node_path_steps.at(exon_idx).count(path_steps_start.first) > 1 && exon_end_node_path_steps.at(exon_idx).count(path_steps_start.first) > 1) {

                assert(use_haplotype_paths);

                is_partial = true;
                break;
            }

            auto haplotype_path_start_it_range = exon_start_node_path_steps.at(exon_idx).equal_range(path_steps_start.first);
            auto haplotype_path_end_it_range = exon_end_node_path_steps.at(exon_idx).equal_range(path_steps_start.first);

            // Transcript paths are partial if either the start or end exon path 
            // step is empty. Partial transcripts are currently not supported.
            // TODO: Add support for partial transcript paths.
            if (haplotype_path_start_it_range.first == haplotype_path_start_it_range.second || haplotype_path_end_it_range.first == haplotype_path_end_it_range.second) {

                is_partial = true;
                break;
            }

            haplotype_path_start_it_range.second--; 
            haplotype_path_end_it_range.second--; 

            auto cur_edited_transcript_paths_base_it = cur_edited_transcript_paths.begin();
            auto cur_edited_transcript_paths_base_eit = cur_edited_transcript_paths.end();

            assert(cur_edited_transcript_paths_base_it != cur_edited_transcript_paths_base_eit);
            cur_edited_transcript_paths_base_eit--;

            while (true) {

                auto border_offsets = cur_transcript.exons.at(exon_idx).border_offsets;

                // Get path step at exon start if exon start node is in the current path.
                auto haplotype_path_start_step = haplotype_path_start_it_range.first->second;

                // Get path mapping at exon end if exon end node is in the current path.
                auto haplotype_path_end_step = haplotype_path_end_it_range.first->second;

                // Exons with different border node orientations are currently not supported
                // for haplotype path projection.
                // TODO: Add support for this special case.
                if (use_haplotype_paths && (_graph->get_is_reverse(_graph->get_handle_of_step(haplotype_path_start_step)) != _graph->get_is_reverse(_graph->get_handle_of_step(haplotype_path_end_step)))) {

                    is_partial = true;
                    break;
                }

                // Swap start and end steps if in reverse order on path
                if (graph_path_pos_overlay.get_position_of_step(haplotype_path_start_step) > graph_path_pos_overlay.get_position_of_step(haplotype_path_end_step)) {

                    assert(border_offsets.first + 1 == _graph->get_length(_graph->get_handle_of_step(haplotype_path_start_step)));
                    assert(border_offsets.second == 0);

                    swap(haplotype_path_start_step, haplotype_path_end_step);
                    border_offsets.first = _graph->get_length(_graph->get_handle_of_step(haplotype_path_start_step)) - 1;
                }

                Path exon_path;
                bool is_first_step = true;

                while (true) {

                    auto node_length = _graph->get_length(_graph->get_handle_of_step(haplotype_path_start_step));
                    int32_t offset = 0;

                    // Adjust start position from exon border (last position in upstream intron)
                    // to first position in exon. Do not adjust if first position in path.
                    if ((cur_transcript.exons.at(exon_idx).coordinates.first > 0) && is_first_step) {

                        if (border_offsets.first + 1 == node_length) {

                            assert(haplotype_path_start_step != haplotype_path_end_step);
                            haplotype_path_start_step = _graph->get_next_step(haplotype_path_start_step);

                            is_first_step = false;
                            continue;
                        
                        } else {

                            offset = border_offsets.first + 1;
                        }
                    }

                    int32_t edit_length = node_length - offset;

                    // Adjust end position from exon border (first position in downstream intron)
                    // to last position in exon. Do not adjust if last position in path.
                    if ((cur_transcript.exons.at(exon_idx).coordinates.second < cur_transcript.chrom_length - 1) && (haplotype_path_start_step == haplotype_path_end_step)) {

                        if (border_offsets.second == 0) {

                            break;

                        } else {

                            edit_length = border_offsets.second - offset;
                        }
                    }

                    assert(0 <= offset && offset < node_length);
                    assert(0 < edit_length && edit_length <= node_length);

                    // Add new mapping in forward direction. Later the whole path will
                    // be reverse complemented if transcript is on the '-' strand.
                    auto new_mapping = exon_path.add_mapping();
                    new_mapping->set_rank(exon_path.mapping_size());

                    new_mapping->mutable_position()->set_node_id(_graph->get_id(_graph->get_handle_of_step(haplotype_path_start_step)));
                    new_mapping->mutable_position()->set_offset(offset);
                    new_mapping->mutable_position()->set_is_reverse(_graph->get_is_reverse(_graph->get_handle_of_step(haplotype_path_start_step)));

                    // Add new edit representing a complete match.
                    auto new_edit = new_mapping->add_edit();
                    new_edit->set_from_length(edit_length);
                    new_edit->set_to_length(edit_length);
                    
                    if (haplotype_path_start_step == haplotype_path_end_step) { break; }

                    haplotype_path_start_step = _graph->get_next_step(haplotype_path_start_step);
                    is_first_step = false;
                }

                if (haplotype_path_start_it_range.first == haplotype_path_start_it_range.second && haplotype_path_end_it_range.first == haplotype_path_end_it_range.second) {

                    auto exon_cur_edited_transcript_paths_base_it = cur_edited_transcript_paths_base_it;

                    while (true) {

                        exon_cur_edited_transcript_paths_base_it->path = concat_paths(exon_cur_edited_transcript_paths_base_it->path, exon_path);

                        if (exon_cur_edited_transcript_paths_base_it == cur_edited_transcript_paths_base_eit) { break; }
                        ++exon_cur_edited_transcript_paths_base_it;
                    }

                    break;
                
                } else {

                    auto exon_cur_edited_transcript_paths_base_it = cur_edited_transcript_paths_base_it;

                    while (true) {

                        // If not last boundary combination copy current base transcipt path.
                        cur_edited_transcript_paths.emplace_back(*exon_cur_edited_transcript_paths_base_it);
                        cur_edited_transcript_paths.back().path = concat_paths(cur_edited_transcript_paths.back().path, exon_path);

                        if (exon_cur_edited_transcript_paths_base_it == cur_edited_transcript_paths_base_eit) { break; }
                        ++exon_cur_edited_transcript_paths_base_it;
                    }
                }

                assert(haplotype_path_start_it_range.first == haplotype_path_start_it_range.second || haplotype_path_end_it_range.first == haplotype_path_end_it_range.second);

                if (haplotype_path_start_it_range.first != haplotype_path_start_it_range.second) {

                    ++haplotype_path_start_it_range.first;
                
                } else {

                    assert(haplotype_path_end_it_range.first != haplotype_path_end_it_range.second);
                    ++haplotype_path_end_it_range.first;
                }
            }
        }

        if (!is_partial) {

            auto cur_edited_transcript_paths_it = cur_edited_transcript_paths.begin();

            while (cur_edited_transcript_paths_it != cur_edited_transcript_paths.end()) {

                if (cur_edited_transcript_paths_it->path.mapping_size() == 0) {
                
                    // Delete empty paths.
                    cur_edited_transcript_paths_it = cur_edited_transcript_paths.erase(cur_edited_transcript_paths_it);

                } else {

                    // Reverse complement transcript paths that are on the '-' strand.
                    if (cur_transcript.is_reverse) {

                        reverse_complement_path_in_place(&(cur_edited_transcript_paths_it->path), [&](vg::id_t node_id) {return _graph->get_length(_graph->get_handle(node_id, false));});
                    } 
                }

                ++cur_edited_transcript_paths_it;
            }

            edited_transcript_paths.splice(edited_transcript_paths.end(), cur_edited_transcript_paths); 
        } 
    } 

    return edited_transcript_paths;
}

list<EditedTranscriptPath> Transcriptome::construct_reference_transcript_paths_gbwt(const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index) const {

    vector<pair<uint32_t, uint32_t> > chrom_transcript_sets;
    string cur_chrom = "";

    // Create sets of transcripts on same chromosome/contig.
    for (size_t i = 0; i < transcripts.size(); ++i) {

        if (cur_chrom != transcripts.at(i).chrom) {

            if (!chrom_transcript_sets.empty()) {

                // Set size of previous set.
                chrom_transcript_sets.back().second = i - chrom_transcript_sets.back().first;
            }

            chrom_transcript_sets.emplace_back(i, 0);
            cur_chrom = transcripts.at(i).chrom;
        }
    }

    // Set size of last set.
    chrom_transcript_sets.back().second = transcripts.size() - chrom_transcript_sets.back().first;    
    sort(chrom_transcript_sets.rbegin(), chrom_transcript_sets.rend(), sort_pair_by_second);

    assert(haplotype_index.bidirectional());
    assert(haplotype_index.hasMetadata());
    
    assert(haplotype_index.metadata.hasPathNames());
    assert(haplotype_index.metadata.hasContigNames());

    spp::sparse_hash_map<string, map<uint32_t, uint32_t> > haplotype_name_index;

    // Create index mapping haplotype contig names to GBWT sequence ids and offsets.
    for (size_t i = 0; i < haplotype_index.sequences(); i++) {

        // Skip reverse threads in bidirectional gbwt index.
        if (i % 2 == 1) {

            continue;
        }

        const auto & path_metadata = haplotype_index.metadata.path(gbwt::Path::id(i));

        auto haplotype_name_index_it = haplotype_name_index.emplace(haplotype_index.metadata.contig(path_metadata.contig), map<uint32_t, uint32_t>());
        assert(haplotype_name_index_it.first->second.emplace(path_metadata.count, i).second);
    }

    list<EditedTranscriptPath> edited_transcript_paths;
    spp::sparse_hash_map<handle_t, vector<EditedTranscriptPath *> > edited_transcript_paths_index;

    mutex edited_transcript_paths_mutex;

    vector<thread> construction_threads;
    construction_threads.reserve(num_threads);

    // Spawn construction threads.
    for (size_t thread_idx = 0; thread_idx < num_threads; thread_idx++) {

        construction_threads.push_back(thread(&Transcriptome::construct_reference_transcript_paths_gbwt_callback, this, &edited_transcript_paths, &edited_transcript_paths_index, &edited_transcript_paths_mutex, thread_idx, ref(chrom_transcript_sets), ref(transcripts), ref(haplotype_index), ref(haplotype_name_index)));
    }

    // Join construction threads.   
    for (auto & thread: construction_threads) {
        
        thread.join();
    }

    return edited_transcript_paths;
}

void Transcriptome::construct_reference_transcript_paths_gbwt_callback(list<EditedTranscriptPath> * edited_transcript_paths, spp::sparse_hash_map<handle_t, vector<EditedTranscriptPath *> > * edited_transcript_paths_index, mutex * edited_transcript_paths_mutex, const int32_t thread_idx, const vector<pair<uint32_t, uint32_t> > & chrom_transcript_sets, const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const spp::sparse_hash_map<string, map<uint32_t, uint32_t> > & haplotype_name_index) const {

    int32_t chrom_transcript_sets_idx = thread_idx;

    while (chrom_transcript_sets_idx < chrom_transcript_sets.size()) {

        list<EditedTranscriptPath> thread_edited_transcript_paths;

        const pair<uint32_t, uint32_t> & transcript_set = chrom_transcript_sets.at(chrom_transcript_sets_idx);

        assert(transcript_set.second > 0);
        uint32_t transcript_idx = transcript_set.first;

        list<pair<EditedTranscriptPath, tuple<uint32_t, uint32_t, bool> > > incomplete_transcript_paths;

        auto haplotype_name_index_it = haplotype_name_index.find(transcripts.at(transcript_idx).chrom);
        assert(haplotype_name_index_it != haplotype_name_index.end());

        for (auto & haplotype_idx: haplotype_name_index_it->second) {

            auto incomplete_transcript_paths_it = incomplete_transcript_paths.begin();

            while (incomplete_transcript_paths_it != incomplete_transcript_paths.end()) {

                // Delete transcripts with exon overlapping haplotype break.
                if (get<2>(incomplete_transcript_paths_it->second)) {
                    
                    cerr << "\tWARNING: Skipping transcript " << transcripts.at(get<0>(incomplete_transcript_paths_it->second)).name << " on " << transcripts.at(get<0>(incomplete_transcript_paths_it->second)).chrom << " since one of its exons overlap a haplotype break." << endl;
                    incomplete_transcript_paths_it = incomplete_transcript_paths.erase(incomplete_transcript_paths_it);
                
                } else {

                    ++incomplete_transcript_paths_it;
                }
            }

            auto node_start_pos = haplotype_idx.first;
            const gbwt::vector_type & gbwt_haplotype = haplotype_index.extract(haplotype_idx.second);

            for (auto & gbwt_node: gbwt_haplotype) {

                auto node_handle = gbwt_to_handle(*_graph, gbwt_node);
                auto node_length = _graph->get_length(node_handle);

                while (transcript_idx < transcript_set.first + transcript_set.second) {

                    const Transcript & cur_transcript = transcripts.at(transcript_idx);

                    // Create new transcript path for transcript with first
                    // exon starting in current node.
                    if (cur_transcript.exons.front().coordinates.first >= node_start_pos && cur_transcript.exons.front().coordinates.first < node_start_pos + node_length) {

                        incomplete_transcript_paths.emplace_back(EditedTranscriptPath(cur_transcript.name, gbwt::Path::id(haplotype_idx.second), true, false), make_tuple(transcript_idx, 0, false));
                    } 

                    if (node_start_pos + node_length <= cur_transcript.exons.front().coordinates.first) {

                        break;
                    }

                    ++transcript_idx;
                }

                incomplete_transcript_paths_it = incomplete_transcript_paths.begin();

                while (incomplete_transcript_paths_it != incomplete_transcript_paths.end()) {

                    const Transcript & cur_transcript = transcripts.at(get<0>(incomplete_transcript_paths_it->second));

                    while (get<1>(incomplete_transcript_paths_it->second) < cur_transcript.exons.size()) {

                        const pair<int32_t, int32_t> & exon_coords = cur_transcript.exons.at(get<1>(incomplete_transcript_paths_it->second)).coordinates;    

                        // Exon is downstream current node.
                        if (node_start_pos + node_length <= exon_coords.first) {

                            break;
                        }

                        int32_t offset = 0;
                        int32_t edit_length = 0;

                        // Exon is starting in current node.
                        if (!get<2>(incomplete_transcript_paths_it->second) && node_start_pos <= exon_coords.first) {

                            offset = exon_coords.first - node_start_pos;
                            get<2>(incomplete_transcript_paths_it->second) = true;
                        }

                        if (get<2>(incomplete_transcript_paths_it->second)) {

                            edit_length = node_length - offset;

                            // Exon is ending in current node.
                            if (exon_coords.second < node_start_pos + node_length) {

                                edit_length = exon_coords.second - node_start_pos - offset + 1;                                

                                get<1>(incomplete_transcript_paths_it->second)++;
                                get<2>(incomplete_transcript_paths_it->second) = false;
                            }

                            assert(0 <= offset && offset < node_length);
                            assert(edit_length > 0 && edit_length <= node_length);

                            // Add new mapping in forward direction. Later the whole path will
                            // be reverse complemented if transcript is on the '-' strand.
                            auto new_mapping = incomplete_transcript_paths_it->first.path.add_mapping();
                            new_mapping->set_rank(incomplete_transcript_paths_it->first.path.mapping_size());

                            new_mapping->mutable_position()->set_node_id(_graph->get_id(node_handle));
                            new_mapping->mutable_position()->set_offset(offset);
                            new_mapping->mutable_position()->set_is_reverse(_graph->get_is_reverse(node_handle));

                            // Add new edit representing a complete match.
                            auto new_edit = new_mapping->add_edit();
                            new_edit->set_from_length(edit_length);
                            new_edit->set_to_length(edit_length);

                            if (node_start_pos + node_length <= exon_coords.second) {

                                break;
                            }

                        } else {

                            break;
                        }
                    }

                    if (get<1>(incomplete_transcript_paths_it->second) == cur_transcript.exons.size()) {

                        // Reverse complement transcript paths that are on the '-' strand.
                        if (cur_transcript.is_reverse) {

                            reverse_complement_path_in_place(&(incomplete_transcript_paths_it->first.path), [&](vg::id_t node_id) {return _graph->get_length(_graph->get_handle(node_id, false));});
                        } 

                        assert(incomplete_transcript_paths_it->first.path.mapping_size() > 0);
                        thread_edited_transcript_paths.emplace_back(move(incomplete_transcript_paths_it->first));

                        incomplete_transcript_paths_it = incomplete_transcript_paths.erase(incomplete_transcript_paths_it);

                    } else {

                        ++incomplete_transcript_paths_it;
                    }
                }

                if (transcript_idx == transcript_set.first + transcript_set.second && incomplete_transcript_paths.empty()) {

                    break;
                }

                node_start_pos += node_length;
            }

            if (transcript_idx == transcript_set.first + transcript_set.second && incomplete_transcript_paths.empty()) {

                break;
            }
        }
    
        for (auto & transcript_path: incomplete_transcript_paths) {

            cerr << "\tWARNING: Skipping transcript " << transcripts.at(get<0>(transcript_path.second)).name << " on " << transcripts.at(get<0>(transcript_path.second)).chrom << " since one of its exons overlap a haplotype break." << endl;
        }

        edited_transcript_paths_mutex->lock();
        remove_redundant_transcript_paths<EditedTranscriptPath>(&thread_edited_transcript_paths, edited_transcript_paths_index);
        edited_transcript_paths->splice(edited_transcript_paths->end(), thread_edited_transcript_paths);
        edited_transcript_paths_mutex->unlock();

        chrom_transcript_sets_idx += num_threads;
    }
}

void Transcriptome::project_haplotype_transcripts(const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const bdsg::PositionOverlay & graph_path_pos_overlay, const bool proj_emded_paths, const float mean_node_length) {

    list<CompletedTranscriptPath> completed_transcript_paths;
    spp::sparse_hash_map<handle_t, vector<CompletedTranscriptPath *> > completed_transcript_paths_index;

    for (auto & transcript_path: _transcript_paths) {

        auto completed_transcript_paths_index_it = completed_transcript_paths_index.emplace(transcript_path.get_first_node_handle(*_graph), vector<CompletedTranscriptPath *>());
        completed_transcript_paths_index_it.first->second.emplace_back(&transcript_path);
    }

    mutex completed_transcript_paths_mutex;

    vector<thread> projection_threads;
    projection_threads.reserve(num_threads);

    // Spawn projection threads.
    for (size_t thread_idx = 0; thread_idx < num_threads; thread_idx++) {

        projection_threads.push_back(thread(&Transcriptome::project_haplotype_transcripts_callback, this, &completed_transcript_paths, &completed_transcript_paths_index, &completed_transcript_paths_mutex, thread_idx, ref(transcripts), ref(haplotype_index), ref(graph_path_pos_overlay), proj_emded_paths, mean_node_length));
    }

    // Join projection threads.   
    for (auto & thread: projection_threads) {
        
        thread.join();
    }

    _transcript_paths.reserve(_transcript_paths.size() + completed_transcript_paths.size());

    for (auto & transcript_path: completed_transcript_paths) {

        _transcript_paths.emplace_back(move(transcript_path));
    }
}

void Transcriptome::project_haplotype_transcripts_callback(list<CompletedTranscriptPath> * completed_transcript_paths, spp::sparse_hash_map<handle_t, vector<CompletedTranscriptPath *> > * completed_transcript_paths_index,  mutex * completed_transcript_paths_mutex, const int32_t thread_idx, const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const bdsg::PositionOverlay & graph_path_pos_overlay, const bool proj_emded_paths, const float mean_node_length) {

    list<CompletedTranscriptPath> thread_completed_transcript_paths;

    int32_t transcripts_idx = thread_idx;

    while (transcripts_idx < transcripts.size()) {

        // Get next transcript belonging to current thread.
        const Transcript & transcript = transcripts.at(transcripts_idx);

        list<CompletedTranscriptPath> completed_transcript_paths;

        if (!haplotype_index.empty()) { 

            // Project transcript onto haplotypes in GBWT index.
            thread_completed_transcript_paths.splice(thread_completed_transcript_paths.end(), construct_completed_transcript_paths(project_transcript_gbwt(transcript, haplotype_index, mean_node_length)));
        }

        if (proj_emded_paths) { 

            // Project transcript onto embedded paths.
            thread_completed_transcript_paths.splice(thread_completed_transcript_paths.end(), construct_completed_transcript_paths(project_transcript_embedded(transcript, graph_path_pos_overlay, false, true)));
        }

        transcripts_idx += num_threads;
    }

    // Add haplotype transcript paths to transcriptome.
    completed_transcript_paths_mutex->lock();
    remove_redundant_transcript_paths<CompletedTranscriptPath>(&thread_completed_transcript_paths, completed_transcript_paths_index);
    completed_transcript_paths->splice(completed_transcript_paths->end(), thread_completed_transcript_paths);
    completed_transcript_paths_mutex->unlock();
}

list<EditedTranscriptPath> Transcriptome::project_transcript_gbwt(const Transcript & cur_transcript, const gbwt::GBWT & haplotype_index, const float mean_node_length) const {

    assert(haplotype_index.bidirectional());

    list<EditedTranscriptPath> edited_transcript_paths;

    vector<pair<vg::id_t, vg::id_t> > exon_node_ids;
    exon_node_ids.reserve(cur_transcript.exons.size());

    vector<pair<vector<exon_nodes_t>, thread_ids_t> > haplotypes;
    multimap<int32_t, pair<int32_t, int32_t> > haplotype_id_index;

    for (size_t exon_idx = 0; exon_idx < cur_transcript.exons.size(); ++exon_idx) {

        const Exon & cur_exon = cur_transcript.exons.at(exon_idx);

        // Add node exon boundary ids
        exon_node_ids.emplace_back(_graph->get_id(_graph->get_handle_of_step(cur_exon.border_steps.first)), _graph->get_id(_graph->get_handle_of_step(cur_exon.border_steps.second)));

        // Calculate expected number of nodes between exon start and end.
        const int32_t expected_length = ceil((cur_exon.coordinates.second - cur_exon.coordinates.first + 1) / mean_node_length);

        // Get all haplotypes in GBWT index between exon start and end border nodes (last position in upstream intron and
        // first position in downstream intron).
        auto exon_haplotypes = get_exon_haplotypes(exon_node_ids.back().first, exon_node_ids.back().second, haplotype_index, expected_length);

        if (haplotypes.empty()) {

            for (auto & exon_haplotype: exon_haplotypes) {

                haplotypes.emplace_back(vector<exon_nodes_t>(1, exon_haplotype.first), exon_haplotype.second);
                haplotypes.back().first.reserve(cur_transcript.exons.size());

                for (auto & haplotype_id: exon_haplotype.second) {

                    haplotype_id_index.emplace(haplotype_id, make_pair(haplotypes.size() - 1, exon_idx + 1));
                }
            }
            
        } else {

            for (auto & exon_haplotype: exon_haplotypes) {

                assert(!exon_haplotype.first.empty());
                spp::sparse_hash_map<int32_t, uint32_t> extended_haplotypes;

                for (auto & haplotype_id: exon_haplotype.second) {

                    auto haplotype_id_index_it_range = haplotype_id_index.equal_range(haplotype_id);
                    auto haplotype_id_index_it = haplotype_id_index_it_range.first;

                    while (haplotype_id_index_it != haplotype_id_index_it_range.second) {

                        if (exon_idx != haplotype_id_index_it->second.second) {

                            assert(haplotype_id_index_it->second.second < exon_idx);

                            haplotype_id_index_it = haplotype_id_index.erase(haplotype_id_index_it);
                            continue;
                        }

                        haplotype_id_index_it->second.second++;
                        pair<vector<exon_nodes_t>, thread_ids_t> * cur_haplotype = &haplotypes.at(haplotype_id_index_it->second.first);

                        if (extended_haplotypes.find(haplotype_id_index_it->second.first) != extended_haplotypes.end()) {

                            assert(cur_haplotype->first.size() == exon_idx + 1);
                            haplotypes.at(extended_haplotypes.at(haplotype_id_index_it->second.first)).second.emplace_back(haplotype_id);
                            haplotype_id_index_it->second.first = extended_haplotypes.at(haplotype_id_index_it->second.first);                        
                        
                        } else if (cur_haplotype->first.size() == exon_idx) {

                            cur_haplotype->first.emplace_back(exon_haplotype.first);
                            cur_haplotype->second = {haplotype_id};
                            assert(extended_haplotypes.emplace(haplotype_id_index_it->second.first, haplotype_id_index_it->second.first).second);
                        
                        } else if (cur_haplotype->first.size() == exon_idx + 1) {

                            haplotypes.emplace_back(vector<exon_nodes_t>(cur_haplotype->first.begin(), cur_haplotype->first.end() - 1), thread_ids_t(1, haplotype_id));
                            haplotypes.back().first.emplace_back(exon_haplotype.first);

                            assert(extended_haplotypes.emplace(haplotype_id_index_it->second.first, haplotypes.size() - 1).second);
                            haplotype_id_index_it->second.first = haplotypes.size() - 1;                
                        
                        } else {

                            haplotype_id_index_it = haplotype_id_index.erase(haplotype_id_index_it);
                            continue;
                        } 

                        ++haplotype_id_index_it;
                    }
                }
            }
        }  
    }

    for (auto & haplotype: haplotypes) {

        // Skip partial transcript paths.
        // TODO: Add support for partial transcript paths.
        if (haplotype.first.size() != cur_transcript.exons.size()) {

            continue;
        }

        auto haplotype_thread_ids_it = haplotype.second.begin();
        assert(haplotype_thread_ids_it != haplotype.second.end());

        // Construct transcript path and set transcript origin name.
        edited_transcript_paths.emplace_back(cur_transcript.name, gbwt::Path::id(*haplotype_thread_ids_it), false, true);
        edited_transcript_paths.back().haplotype_gbwt_ids.reserve(haplotype.second.size());

        // Add haplotype names as origins.
        while (haplotype_thread_ids_it != haplotype.second.end()) {

            // Convert bidirectional path id before finding name. 
            edited_transcript_paths.back().haplotype_gbwt_ids.emplace_back(gbwt::Path::id(*haplotype_thread_ids_it));
            ++haplotype_thread_ids_it;
        }

        for (size_t exon_idx = 0; exon_idx < cur_transcript.exons.size(); ++exon_idx) {

            const Exon & cur_exon = cur_transcript.exons.at(exon_idx);

            assert(gbwt::Node::id(haplotype.first.at(exon_idx).front()) == exon_node_ids.at(exon_idx).first);
            assert(gbwt::Node::id(haplotype.first.at(exon_idx).back()) == exon_node_ids.at(exon_idx).second);

            for (size_t exon_node_idx = 0; exon_node_idx < haplotype.first.at(exon_idx).size(); ++exon_node_idx) {

                assert(haplotype.first.at(exon_idx).at(exon_node_idx) != gbwt::ENDMARKER);

                auto node_id = gbwt::Node::id(haplotype.first.at(exon_idx).at(exon_node_idx));
                auto node_length = _graph->get_length(_graph->get_handle(node_id, false));

                int32_t offset = 0;

                // Adjust start position from exon border (last position in upstream intron)
                // to first position in exon. Do not adjust if first position in path.
                if ((cur_exon.coordinates.first > 0) && (exon_node_idx == 0)) {

                    if (cur_exon.border_offsets.first + 1 == node_length) {

                        assert(haplotype.first.at(exon_idx).size() > 1);
                        assert(node_id != exon_node_ids.at(exon_idx).second);

                        continue;
                    
                    } else {

                        offset = cur_exon.border_offsets.first + 1;
                    }
                }

                int32_t edit_length = node_length - offset;

                // Adjust end position from exon border (first position in downstream intron)
                // to last position in exon. Do not adjust if last position in path.
                if ((cur_exon.coordinates.second < cur_transcript.chrom_length - 1) && (exon_node_idx == haplotype.first.at(exon_idx).size() - 1)) {

                    if (cur_exon.border_offsets.second == 0) {

                        break;

                    } else {

                        edit_length = cur_exon.border_offsets.second - offset;
                    }
                }

                assert(0 <= offset && offset < node_length);
                assert(0 < edit_length && edit_length <= node_length);

                // Add new mapping in forward direction. Later the whole path will
                // be reverse complemented if transcript is on the '-' strand.
                auto new_mapping = edited_transcript_paths.back().path.add_mapping();
                new_mapping->set_rank(edited_transcript_paths.back().path.mapping_size());

                new_mapping->mutable_position()->set_node_id(node_id);
                new_mapping->mutable_position()->set_offset(offset);
                new_mapping->mutable_position()->set_is_reverse(false);

                // Add new edit representing a complete match.
                auto new_edit = new_mapping->add_edit();
                new_edit->set_from_length(edit_length);
                new_edit->set_to_length(edit_length);
            }
        }

        if (edited_transcript_paths.back().path.mapping_size() == 0) {

            // Delete empty paths.
            edited_transcript_paths.pop_back();

        } else {

            if (cur_transcript.is_reverse) {

                // Reverse complement transcript paths that are on the '-' strand.
                reverse_complement_path_in_place(&(edited_transcript_paths.back().path), [&](vg::id_t node_id) {return _graph->get_length(_graph->get_handle(node_id, false));});
            }

            // Copy paths if collapse of identical transcript paths is not wanted.
            if (path_collapse_type == "no" && edited_transcript_paths.back().haplotype_gbwt_ids.size() > 1) {

                auto all_haplotype_gbwt_ids = edited_transcript_paths.back().haplotype_gbwt_ids;
                edited_transcript_paths.back().haplotype_gbwt_ids = {edited_transcript_paths.back().haplotype_gbwt_ids.front()};

                // Create identical copies of all haplotype origins.
                for (size_t i = 1; i < all_haplotype_gbwt_ids.size(); ++i) {

                    edited_transcript_paths.emplace_back(edited_transcript_paths.back());
                    edited_transcript_paths.back().haplotype_gbwt_ids.front() = {all_haplotype_gbwt_ids.at(i)};
                }
            }
        }
    }

    return edited_transcript_paths; 
}

vector<pair<exon_nodes_t, thread_ids_t> > Transcriptome::get_exon_haplotypes(const vg::id_t start_node, const vg::id_t end_node, const gbwt::GBWT & haplotype_index, const int32_t expected_length) const {

    assert(expected_length > 0);

    // Calculate the expected upperbound of the length between the two 
    // nodes (number of nodes). 
    const int32_t expected_length_upperbound = 1.1 * expected_length;

    // Calcuate frequency for how often a check on whether an extension 
    // should be terminated is performed. 
    const int32_t termination_frequency = ceil(0.1 * expected_length);

    // Get ids for haplotypes that contain the end node.
    spp::sparse_hash_set<int32_t> end_haplotype_ids;
    for (auto & haplotype_id: haplotype_index.locate(haplotype_index.find(gbwt::Node::encode(end_node, false)))) {

        end_haplotype_ids.emplace(haplotype_id);
    }

    vector<pair<exon_nodes_t, thread_ids_t> > exon_haplotypes;

    // Initialise haplotype extension queue on the start node.
    std::queue<pair<exon_nodes_t, gbwt::SearchState> > exon_haplotype_queue;
    exon_haplotype_queue.push(make_pair(exon_nodes_t(1, gbwt::Node::encode(start_node, false)), haplotype_index.find(gbwt::Node::encode(start_node, false))));
    exon_haplotype_queue.front().first.reserve(expected_length_upperbound);

    // Empty queue if no haplotypes containing the start node exist. 
    if (exon_haplotype_queue.front().second.empty()) { exon_haplotype_queue.pop(); }

    // Perform depth-first haplotype extension.
    while (!exon_haplotype_queue.empty()) {

        pair<exon_nodes_t, gbwt::SearchState> & cur_exon_haplotype = exon_haplotype_queue.front();

        // Stop current extension if end node is reached.
        if (gbwt::Node::id(cur_exon_haplotype.first.back()) == end_node) {

            exon_haplotypes.emplace_back(cur_exon_haplotype.first, haplotype_index.locate(cur_exon_haplotype.second));
            assert(exon_haplotypes.back().second.size() <= cur_exon_haplotype.second.size());

            if (exon_haplotypes.back().second.size() == cur_exon_haplotype.second.size()) {

                exon_haplotype_queue.pop();
                continue;   
            }          
        }

        // Check whether any haplotypes in the current extension contains the
        // end node. If not, end current extension. This check is only performed
        // after the upperbound on the expected number of nodes is reached. 
        if (cur_exon_haplotype.first.size() >= expected_length_upperbound && (cur_exon_haplotype.first.size() % termination_frequency) == 0) {

            bool has_relevant_haplotype = false;

            for (auto & haplotype_id: haplotype_index.locate(cur_exon_haplotype.second)) {

                if (end_haplotype_ids.find(haplotype_id) != end_haplotype_ids.end()) {

                    has_relevant_haplotype = true;
                    break;
                }
            }

            if (!has_relevant_haplotype) {

                exon_haplotype_queue.pop();
                continue;                
            }
        }
        
        auto out_edges = haplotype_index.edges(cur_exon_haplotype.first.back());

        // End current extension if no outgoing edges exist.
        if (out_edges.empty()) {

            exon_haplotype_queue.pop();
            continue;
        }

        auto out_edges_it = out_edges.begin(); 
        ++out_edges_it;

        while (out_edges_it != out_edges.end()) {

            // Do not extend haplotypes that end within the exon.
            if (out_edges_it->first != gbwt::ENDMARKER) {

                auto extended_search = haplotype_index.extend(cur_exon_haplotype.second, out_edges_it->first);

                // Add new extension to queue if not empty (haplotypes found).
                if (!extended_search.empty()) { 

                    exon_haplotype_queue.push(make_pair(cur_exon_haplotype.first, extended_search));
                    exon_haplotype_queue.back().first.emplace_back(out_edges_it->first);
                }
            }

            ++out_edges_it;
        }

        // Do not extend haplotypes that end within the exon.
        if (out_edges.begin()->first != gbwt::ENDMARKER) {

            cur_exon_haplotype.first.emplace_back(out_edges.begin()->first);
            cur_exon_haplotype.second = haplotype_index.extend(cur_exon_haplotype.second, out_edges.begin()->first);        

            // End current extension if empty (no haplotypes found). 
            if (cur_exon_haplotype.second.empty()) { exon_haplotype_queue.pop(); }
    
        } else {

            exon_haplotype_queue.pop();
        }
    }

    return exon_haplotypes;
}

template <class T>
void Transcriptome::remove_redundant_transcript_paths(list<T> * new_transcript_paths, spp::sparse_hash_map<handle_t, vector<T*> > * transcript_paths_index) const {

    auto new_transcript_paths_it = new_transcript_paths->begin();

    while (new_transcript_paths_it != new_transcript_paths->end()) {

        assert(!new_transcript_paths_it->transcript_names.empty());

        bool unique_transcript_path = true;

        auto transcript_paths_index_it = transcript_paths_index->emplace(new_transcript_paths_it->get_first_node_handle(*_graph), vector<T*>());

        // Add unique transcript paths only.
        if (!transcript_paths_index_it.second && path_collapse_type != "no") {

            assert(!transcript_paths_index_it.first->second.empty());

            for (auto & transcript_path: transcript_paths_index_it.first->second) {

                assert(!transcript_path->transcript_names.empty());

                if (path_collapse_type == "all" || transcript_path->transcript_names.front() == new_transcript_paths_it->transcript_names.front()) {

                    // Check if two paths are identical.
                    if (transcript_path->path == new_transcript_paths_it->path) {

                        if (path_collapse_type == "all") {

                            for (auto & new_transcript_name: new_transcript_paths_it->transcript_names) {

                                if (find(transcript_path->transcript_names.begin(), transcript_path->transcript_names.end(), new_transcript_name) != transcript_path->transcript_names.end()) {

                                    transcript_path->transcript_names.emplace_back(new_transcript_name);
                                }
                            }

                        } else {

                            assert(path_collapse_type == "haplotype");
                            assert(transcript_path->transcript_names.size() == 1);
                            assert(new_transcript_paths_it->transcript_names.size() == 1);
                        }

                        // Merge embedded path origin names.
                        transcript_path->embedded_path_names.insert(transcript_path->embedded_path_names.end(), new_transcript_paths_it->embedded_path_names.begin(), new_transcript_paths_it->embedded_path_names.end());

                        // Merge haplotype gbwt origin ids.
                        transcript_path->haplotype_gbwt_ids.insert(transcript_path->haplotype_gbwt_ids.end(), new_transcript_paths_it->haplotype_gbwt_ids.begin(), new_transcript_paths_it->haplotype_gbwt_ids.end());

                        transcript_path->is_reference = (transcript_path->is_reference || new_transcript_paths_it->is_reference);
                        transcript_path->is_haplotype = (transcript_path->is_haplotype || new_transcript_paths_it->is_haplotype);
                        
                        // Delete non-unique transcript path.
                        new_transcript_paths_it = new_transcript_paths->erase(new_transcript_paths_it);

                        unique_transcript_path = false;
                        break;
                    }
                }
            }
        }

        if (unique_transcript_path) {

            transcript_paths_index_it.first->second.emplace_back(&(*new_transcript_paths_it));
            ++new_transcript_paths_it;
        }        
    } 
}

list<CompletedTranscriptPath> Transcriptome::construct_completed_transcript_paths(const list<EditedTranscriptPath> & edited_transcript_paths) const {

    list<CompletedTranscriptPath> completed_transcript_paths;

    for (auto & transcript_path: edited_transcript_paths) {

        completed_transcript_paths.emplace_back(transcript_path, *_graph);
    }

    return completed_transcript_paths;     
}

void Transcriptome::add_edited_transcript_paths(const list<EditedTranscriptPath> & edited_transcript_paths) {

    add_splice_junction_edges(edited_transcript_paths);

    for (auto & transcript_path: edited_transcript_paths) {

        _transcript_paths.emplace_back(transcript_path, *_graph);
    }
}

bool Transcriptome::has_novel_exon_boundaries(const list<EditedTranscriptPath> & edited_transcript_paths, const bool include_transcript_ends) const {

    for (auto & transcript_path: edited_transcript_paths) {

        for (size_t i = 0; i < transcript_path.path.mapping_size(); i++) {

            auto cur_mapping = transcript_path.path.mapping(i);
            auto cur_handle = mapping_to_handle(cur_mapping, *_graph);

            assert(cur_mapping.edit_size() == 1);
            assert(edit_is_match(cur_mapping.edit(0)));

            // Do not check if left boundary of start exon is novel.
            if (!include_transcript_ends && i == 0) {

                if (cur_mapping.position().offset() + cur_mapping.edit(0).from_length() != _graph->get_length(cur_handle)) {

                    return true;
                }

            // Do not check if right boundary of end exon is novel.
            } else if (!include_transcript_ends && i == transcript_path.path.mapping_size() - 1) {

                if (cur_mapping.position().offset() > 0) {

                    return true;
                }

            // Check if both boundaries are novel.
            } else if (cur_mapping.position().offset() > 0 || cur_mapping.edit(0).from_length() != _graph->get_length(cur_handle)) {

                return true;
            }
        }
    }

    return false;
}

void Transcriptome::augment_graph(const list<EditedTranscriptPath> & edited_transcript_paths, const bool break_at_transcript_ends, unique_ptr<gbwt::GBWT> & haplotype_index, const bool update_haplotypes, const bool add_reference_transcript_paths) {

#ifdef transcriptome_debug
    double time_convert_1 = gcsa::readTimer();
    cerr << "\t\tDEBUG Creation start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    // Create set of exon boundary paths to augment graph with.
    vector<Path> exon_boundary_paths;
    spp::sparse_hash_set<Mapping, MappingHash> exon_boundary_mapping_index;

    for (auto & transcript_path: edited_transcript_paths) {

        for (size_t j = 0; j < transcript_path.path.mapping_size(); ++j) {

            const Mapping & mapping = transcript_path.path.mapping(j);

            const auto mapping_length = mapping_to_length(mapping);
            assert(mapping_length == mapping_from_length(mapping));

            // Add exon boundary path.
            if (mapping.position().offset() > 0 || mapping.position().offset() + mapping_length < _graph->get_length(_graph->get_handle(mapping.position().node_id(), false))) {

                exon_boundary_paths.emplace_back(Path());
                *(exon_boundary_paths.back().add_mapping()) = mapping; 
                exon_boundary_paths.back().mutable_mapping(0)->set_rank(1);

                // Remove if already added.
                if (!exon_boundary_mapping_index.emplace(exon_boundary_paths.back().mapping(0)).second) {
                
                    exon_boundary_paths.pop_back();
                }
            }  
        }
    }   

#ifdef transcriptome_debug
    cerr << "\t\tDEBUG Created " << exon_boundary_paths.size() << " exon boundary paths: " << gcsa::readTimer() - time_convert_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_augment_1 = gcsa::readTimer();
    cerr << "\t\tDEBUG Augmention start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    vector<Translation> translations;

    // Augment graph with edited paths. 
    augment(static_cast<MutablePathMutableHandleGraph *>(_graph.get()), exon_boundary_paths, "GAM", &translations, "", false, break_at_transcript_ends);

#ifdef transcriptome_debug
    cerr << "\t\tDEBUG Augmented graph with " << translations.size() << " translations: " << gcsa::readTimer() - time_augment_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_index_1 = gcsa::readTimer();
    cerr << "\t\tDEBUG Indexing start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    spp::sparse_hash_map<handle_t, vector<pair<int32_t, handle_t> > > translation_index;

    #pragma omp parallel num_threads(num_threads)
    {
        spp::sparse_hash_map<handle_t, vector<pair<int32_t, handle_t> > > thread_translation_index;

        // Create translation index 
        #pragma omp for schedule(static)
        for (size_t i = 0; i < translations.size(); ++i) {

            const Translation & translation = translations.at(i);

            assert(translation.from().mapping_size() == 1);
            assert(translation.to().mapping_size() == 1);

            auto & from_mapping = translation.from().mapping(0);
            auto & to_mapping = translation.to().mapping(0);

            assert(to_mapping.position().offset() == 0);
            assert(from_mapping.position().is_reverse() == to_mapping.position().is_reverse());

            // Only store changes
            if (from_mapping != to_mapping) {

                auto thread_translation_index_it = thread_translation_index.emplace(mapping_to_handle(from_mapping, *_graph), vector<pair<int32_t, handle_t> >());
                thread_translation_index_it.first->second.emplace_back(from_mapping.position().offset(), mapping_to_handle(to_mapping, *_graph));
            }
        }

        #pragma omp critical
        {
            for (auto & translation: thread_translation_index) {

                auto translation_index_it = translation_index.emplace(translation.first, translation.second);

                if (!translation_index_it.second) {

                    translation_index_it.first->second.insert(translation_index_it.first->second.end(), translation.second.begin(), translation.second.end());
                }
            }
        }
    }

    // Sort translation index by offset
    for (auto & translation: translation_index) {

        sort(translation.second.begin(), translation.second.end());
    }

#ifdef transcriptome_debug
    cerr << "\t\tDEBUG Indexed " << translation_index.size() << " translated nodes: " << gcsa::readTimer() - time_index_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    if (!haplotype_index->empty() && update_haplotypes) {

#ifdef transcriptome_debug
    double time_update_1 = gcsa::readTimer();
    cerr << "\t\tDEBUG Updating (GBWT) start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

        // Update haplotypes in gbwt index to match new augmented graph.
        update_haplotype_index(haplotype_index, translation_index);

#ifdef transcriptome_debug
    cerr << "\t\tDEBUG Updated " << haplotype_index->sequences() / 2 << " haplotype paths: " << gcsa::readTimer() - time_update_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    } 

    if (!_transcript_paths.empty()) {

#ifdef transcriptome_debug
    double time_update_2 = gcsa::readTimer();
    cerr << "\t\tDEBUG Updating (transcriptome) start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

        // Update transcript paths in transcriptome to new augmented graph.
        update_transcript_paths(translation_index);

#ifdef transcriptome_debug
    cerr << "\t\tDEBUG Updated " << _transcript_paths.size() << " transcript paths: " << gcsa::readTimer() - time_update_2 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    } 

#ifdef transcriptome_debug
    double time_update_3 = gcsa::readTimer();
    cerr << "\t\tDEBUG Updating (paths) start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    list<CompletedTranscriptPath> updated_transcript_paths;

    // Update paths to match new augmented graph and add them
    // as reference transcript paths.
    for (auto & transcript_path: edited_transcript_paths) {

        updated_transcript_paths.emplace_back(transcript_path);

        for (auto mapping: transcript_path.path.mapping()) {

            auto mapping_handle = mapping_to_handle(mapping, *_graph);
            auto mapping_offset = mapping.position().offset();
            auto mapping_length = mapping_to_length(mapping);

            assert(mapping_length > 0);
            assert(mapping_length == mapping_from_length(mapping));

            auto translation_index_it = translation_index.find(mapping_handle);

            if (translation_index_it != translation_index.end()) {

                // First node id is the same (new node offset is larger than 0). 
                if (mapping_offset == 0 & translation_index_it->second.front().first > 0) {

                    updated_transcript_paths.back().path.emplace_back(mapping_handle);
                }

                // Add new nodes.
                for (auto & new_node: translation_index_it->second) {

                    if (new_node.first >= mapping_offset && new_node.first < mapping_offset + mapping_length) {

                        updated_transcript_paths.back().path.emplace_back(new_node.second);
                    }
                }

            } else {

                updated_transcript_paths.back().path.emplace_back(mapping_handle);
            }
        }
    }

    add_splice_junction_edges(updated_transcript_paths);

    if (add_reference_transcript_paths) {

        if (!_transcript_paths.empty()) {

            spp::sparse_hash_map<handle_t, vector<CompletedTranscriptPath *> > transcript_paths_index;

            for (auto & transcript_path: _transcript_paths) {

                auto transcript_paths_index_it = transcript_paths_index.emplace(transcript_path.get_first_node_handle(*_graph), vector<CompletedTranscriptPath *>());
                transcript_paths_index_it.first->second.emplace_back(&transcript_path);
            }

            remove_redundant_transcript_paths<CompletedTranscriptPath>(&updated_transcript_paths, &transcript_paths_index);
        }

        _transcript_paths.reserve(_transcript_paths.size() + updated_transcript_paths.size());

        for (auto & transcript_path: updated_transcript_paths) {

            _transcript_paths.emplace_back(move(transcript_path));
        }
    }

#ifdef transcriptome_debug
    cerr << "\t\tDEBUG Updated " << updated_transcript_paths.size() << " transcript paths: " << gcsa::readTimer() - time_update_3 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif  
}

void Transcriptome::update_haplotype_index(unique_ptr<gbwt::GBWT> & haplotype_index, const spp::sparse_hash_map<handle_t, vector<pair<int32_t, handle_t> > > & update_index) const {

    assert(haplotype_index->bidirectional());

    // Silence GBWT index construction. 
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT); 
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(_graph->max_node_id(), true)));

    // Transfer metadata
    gbwt_builder.index.addMetadata();
    gbwt_builder.index.metadata = haplotype_index->metadata;

    for (size_t i = 0; i < haplotype_index->sequences(); i++) {

        // Only update forward threads in bidirectional gbwt index.
        if (i % 2 == 1) {

            continue;
        }

        auto cur_gbwt_thread = haplotype_index->extract(i);

        gbwt::vector_type new_gbwt_threads;
        new_gbwt_threads.reserve(cur_gbwt_thread.size());

        for (auto & node: cur_gbwt_thread) {

            auto handle = gbwt_to_handle(*_graph, node);
            assert(_graph->get_is_reverse(handle) == 0);

            auto update_index_it = update_index.find(handle);

            if (update_index_it != update_index.end()) {

                // First node id is the same (new node offset is 
                // larger than 0). 
                if (update_index_it->second.front().first > 0) {

                    new_gbwt_threads.emplace_back(node);
                }

                // Add new nodes.
                for (auto & new_node: update_index_it->second) {

                    assert(_graph->get_is_reverse(new_node.second) == 0);
                    new_gbwt_threads.emplace_back(handle_to_gbwt(*_graph, new_node.second));
                }

            } else {

                new_gbwt_threads.emplace_back(node);
            }
        }

        // Insert thread bidirectionally.
        gbwt_builder.insert(new_gbwt_threads, true);
    }

    // Finish contruction and recode index.
    gbwt_builder.finish();
    haplotype_index.reset(new gbwt::GBWT(gbwt_builder.index));   
}   

void Transcriptome::update_transcript_paths(const spp::sparse_hash_map<handle_t, vector<pair<int32_t, handle_t> > > & update_index) {

    #pragma omp parallel num_threads(num_threads)
    {
        // Update transcript paths 
        #pragma omp for schedule(static)
        for (size_t i = 0; i < _transcript_paths.size(); ++i) {

            vector<handle_t> new_transcript_path;
            new_transcript_path.reserve(_transcript_paths.at(i).path.size());

            for (auto & handle: _transcript_paths.at(i).path) {

                auto update_index_it = update_index.find(handle); 

                if (update_index_it != update_index.end()) {

                    // First handle is the same (new node offset is 
                    // larger than 0). 
                    if (update_index_it->second.front().first > 0) {

                        new_transcript_path.emplace_back(handle);
                    }

                    // Add new handles.
                    for (auto & new_handle: update_index_it->second) {
                    
                        new_transcript_path.emplace_back(new_handle.second);
                    }
                
                } else { 

                    update_index_it = update_index.find(_graph->flip(handle));

                    if (update_index_it != update_index.end()) {

                        // First handle is the same (new node offset is 
                        // larger than 0). 
                        if (update_index_it->second.front().first > 0) {

                            new_transcript_path.emplace_back(handle);
                        }

                        for (auto update_handle_rit = update_index_it->second.rbegin(); update_handle_rit != update_index_it->second.rend(); ++update_handle_rit) {
                        
                            new_transcript_path.emplace_back(_graph->flip(update_handle_rit->second));
                        }

                    } else {

                        new_transcript_path.emplace_back(handle);
                    }
                }
            }

            _transcript_paths.at(i).path = move(new_transcript_path);
        }
    }
}

void Transcriptome::add_splice_junction_edges(const list<EditedTranscriptPath> & edited_transcript_paths) {

    for (auto & transcript_path: edited_transcript_paths) {

        for (size_t i = 1; i < transcript_path.path.mapping_size(); i++) {

            auto & prev_mapping = transcript_path.path.mapping(i - 1);
            auto & cur_mapping = transcript_path.path.mapping(i);

            auto prev_handle = mapping_to_handle(prev_mapping, *_graph);
            auto cur_handle = mapping_to_handle(cur_mapping, *_graph);
            
            // Ensure the edge exists.
            _graph->create_edge(prev_handle, cur_handle);
        }
    }
}

void Transcriptome::add_splice_junction_edges(const list<CompletedTranscriptPath> & completed_transcript_paths) {

    for (auto & transcript_path: completed_transcript_paths) {

        for (size_t i = 1; i < transcript_path.path.size(); i++) {
            
            // Ensure the edge exists.
            _graph->create_edge(transcript_path.path.at(i - 1), transcript_path.path.at(i));
        }
    }
}

  
void Transcriptome::add_splice_junction_edges(const vector<CompletedTranscriptPath> & completed_transcript_paths) {

    for (auto & transcript_path: completed_transcript_paths) {

        for (size_t i = 1; i < transcript_path.path.size(); i++) {
            
            // Ensure the edge exists.
            _graph->create_edge(transcript_path.path.at(i - 1), transcript_path.path.at(i));
        }
    }
}

void Transcriptome::update_copy_id() {

    spp::sparse_hash_map<string, uint32_t> name_copy_id_index;

    for (auto & transcript_path: _transcript_paths) {

        sort(transcript_path.transcript_names.begin(), transcript_path.transcript_names.end());

        assert(!transcript_path.transcript_names.empty());
        assert(!transcript_path.transcript_names.front().empty());

        auto name_copy_id_index_it = name_copy_id_index.emplace(transcript_path.transcript_names.front(), 0);
        name_copy_id_index_it.first->second++;

        transcript_path.copy_id = name_copy_id_index_it.first->second;
    }
}

const vector<CompletedTranscriptPath> & Transcriptome::transcript_paths() const {

    return _transcript_paths;
}

vector<CompletedTranscriptPath> Transcriptome::reference_transcript_paths() const {

    vector<CompletedTranscriptPath> reference_transcript_paths;

    for (auto & transcript_path: _transcript_paths) {

        assert(transcript_path.is_reference || transcript_path.is_haplotype);

        if (transcript_path.is_reference) {

            reference_transcript_paths.emplace_back(transcript_path);
        }
    }

    return reference_transcript_paths;
}

vector<CompletedTranscriptPath> Transcriptome::haplotype_transcript_paths() const {

    vector<CompletedTranscriptPath> haplotype_transcript_paths;

    for (auto & transcript_path: _transcript_paths) {

        assert(transcript_path.is_reference || transcript_path.is_haplotype);

        if (transcript_path.is_haplotype) {

            haplotype_transcript_paths.emplace_back(transcript_path);
        }
    }

    return haplotype_transcript_paths;
}

const MutablePathDeletableHandleGraph & Transcriptome::graph() const {

    return *_graph;
}

void Transcriptome::collect_transcribed_nodes(spp::sparse_hash_set<nid_t> * transcribed_nodes) const {

    for (auto & transcript_path: _transcript_paths) {

        assert(transcript_path.path.size() > 0);
        for (auto & handle: transcript_path.path) {

            transcribed_nodes->emplace(_graph->get_id(handle));
        }    
    } 
}

void Transcriptome::remove_non_transcribed_nodes() {

    vector<path_handle_t> path_handles;
    path_handles.reserve(_graph->get_path_count());

    assert(_graph->for_each_path_handle([&](const path_handle_t & path_handle) {

        path_handles.emplace_back(path_handle);
    }));    

    // Remove all paths.
    for (auto & path_handle: path_handles) {

        _graph->destroy_path(path_handle);
    }

    assert(_graph->get_path_count() == 0);

    // Find all nodes that are in a transcript path.
    spp::sparse_hash_set<nid_t> transcribed_nodes;

    collect_transcribed_nodes(&transcribed_nodes);

    vector<handle_t> non_transcribed_handles;
    non_transcribed_handles.reserve(_graph->get_node_count() - transcribed_nodes.size());

    // Collect all nodes that are not in a transcript path.
    assert(_graph->for_each_handle([&](const handle_t & handle) {
        
        if (transcribed_nodes.count(_graph->get_id(handle)) == 0) {

            non_transcribed_handles.emplace_back(handle); 
        } 
    }));

    for (auto & handle: non_transcribed_handles) {

        // Delete node and in/out edges. 
        _graph->destroy_handle(handle);
    }

    assert(_graph->get_node_count() == transcribed_nodes.size());
}

uint32_t Transcriptome::chop_nodes(const uint32_t max_node_length) {

    spp::sparse_hash_map<handle_t, vector<pair<int32_t, handle_t> > > split_index;

    assert(_graph->for_each_handle([&](const handle_t & handle) {

        const uint32_t handle_length = _graph->get_length(handle);

        if (handle_length > max_node_length) {

            vector<size_t> offsets;
            offsets.reserve(ceil(handle_length / static_cast<float>(max_node_length)));

            uint32_t offset = max_node_length;

            while (offset < handle_length) {
            
                offsets.emplace_back(offset);
                offset += max_node_length;
            }

            auto split_index_it = split_index.emplace(handle, vector<pair<int32_t, handle_t> >());
            assert(split_index_it.second);

            for (auto & div_handle: _graph->divide_handle(handle, offsets)) {

                split_index_it.first->second.emplace_back(0, div_handle);
            }
        }
    }));

    update_transcript_paths(split_index);

    return split_index.size();
}

bool Transcriptome::sort_compact_nodes() {

    if (dynamic_cast<bdsg::PackedGraph*>(_graph.get()) == nullptr) {

        return false;
    }

#ifdef transcriptome_debug
    double time_sort_1 = gcsa::readTimer();
    cerr << "\tDEBUG Sorting start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    auto new_order = handlealgs::topological_order(_graph.get());
    assert(new_order.size() == _graph->get_node_count());

#ifdef transcriptome_debug
    cerr << "\tDEBUG Sorted " << new_order.size() << " nodes: " << gcsa::readTimer() - time_sort_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_index_1 = gcsa::readTimer();
    cerr << "\tDEBUG Indexing start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    spp::sparse_hash_map<handle_t, vector<pair<int32_t, handle_t> > > update_index; 
    uint32_t order_idx = 1;

    for (auto handle: new_order) {
        
        assert(update_index.emplace(handle, vector<pair<int32_t, handle_t> >({make_pair(0, _graph->get_handle(order_idx, _graph->get_is_reverse(handle)))})).second);
        ++order_idx;
    }

#ifdef transcriptome_debug
    cerr << "\tDEBUG Indexed " << update_index.size() << " node updates: " << gcsa::readTimer() - time_index_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_update_1 = gcsa::readTimer();
    cerr << "\tDEBUG Updating (graph) start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    _graph->apply_ordering(new_order, true);

#ifdef transcriptome_debug
    cerr << "\tDEBUG Updated graph: " << gcsa::readTimer() - time_update_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_update_2 = gcsa::readTimer();
    cerr << "\tDEBUG Updating (paths) start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    update_transcript_paths(update_index);

#ifdef transcriptome_debug
    cerr << "\tDEBUG Updated " << _transcript_paths.size() << " transcript paths: " << gcsa::readTimer() - time_update_2 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    return true;
}

int32_t Transcriptome::embed_transcript_paths(const bool add_reference_transcripts, const bool add_haplotype_transcripts) {

    assert(add_reference_transcripts || add_haplotype_transcripts);

    int32_t num_embedded_paths = 0;

    // Add transcript paths to graph
    for (auto & transcript_path: _transcript_paths) {

        assert(transcript_path.is_reference || transcript_path.is_haplotype);

        // Exclude reference or haplotype transcripts
        if (!((transcript_path.is_reference && add_reference_transcripts) || (transcript_path.is_haplotype && add_haplotype_transcripts)))  {

            continue;
        }

        ++num_embedded_paths;

        assert(!_graph->has_path(transcript_path.get_name()));

        auto path_handle = _graph->create_path_handle(transcript_path.get_name());

        for (auto & handle: transcript_path.path) {

            _graph->append_step(path_handle, handle);
        }
    }

    return num_embedded_paths;
}

int32_t Transcriptome::add_transcripts_to_gbwt(gbwt::GBWTBuilder * gbwt_builder, const bool add_bidirectional, const bool output_reference_transcripts) const {

    int32_t num_added_threads = 0;

    vector<string> sample_names;
    sample_names.reserve(_transcript_paths.size());

    if (!gbwt_builder->index.hasMetadata()) {

        gbwt_builder->index.addMetadata();
    }

    // Get current number of haplotypes in GBWT index.
    auto pre_num_haplotypes = gbwt_builder->index.metadata.haplotypes();

    for (auto & transcript_path: _transcript_paths) {

        assert(transcript_path.is_reference || transcript_path.is_haplotype);

        // Exclude reference transcripts
        if (transcript_path.is_reference && !output_reference_transcripts) {

            continue;
        }

        ++num_added_threads;

        // Convert transcript path to GBWT thread.
        gbwt::vector_type gbwt_thread(transcript_path.path.size());
        for (size_t i = 0; i < transcript_path.path.size(); i++) {

            gbwt_thread[i] = handle_to_gbwt(*_graph, transcript_path.path.at(i));
        }

        // Insert transcript path as thread into GBWT index.
        gbwt_builder->insert(gbwt_thread, add_bidirectional);

        // Insert transcript path name into GBWT index.
        gbwt_builder->index.metadata.addPath(pre_num_haplotypes + sample_names.size(), 0, 0, 0);

        sample_names.emplace_back(transcript_path.get_name());
    }

    sample_names.shrink_to_fit();

    // Set number number of haplotypes and transcript path name in metadata.
    gbwt_builder->index.metadata.setHaplotypes(pre_num_haplotypes + sample_names.size());
    gbwt_builder->index.metadata.addSamples(sample_names);

    return num_added_threads;
}

int32_t Transcriptome::write_transcript_sequences(ostream * fasta_ostream, const bool output_reference_transcripts) const {

    int32_t num_written_sequences = 0;

    for (auto & transcript_path: _transcript_paths) {

        assert(transcript_path.is_reference || transcript_path.is_haplotype);

        // Exclude reference transcripts
        if (transcript_path.is_reference && !output_reference_transcripts) {

            continue;
        }

        ++num_written_sequences;

        // Construct transcript path sequence.
        string transcript_path_sequence = "";
        for (auto & handle: transcript_path.path) {

            transcript_path_sequence += _graph->get_sequence(handle);
        }

        // Write transcript path name and sequence.
        write_fasta_sequence(transcript_path.get_name(), transcript_path_sequence, *fasta_ostream);
    }

    return num_written_sequences;
}

int32_t Transcriptome::write_transcript_info(ostream * tsv_ostream, const gbwt::GBWT & haplotype_index, const bool output_reference_transcripts) const {

    *tsv_ostream << "Name\tLength\tTranscripts\tHaplotypes" << endl; 
    
    // Pre-parse GBWT metadata
    auto gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(haplotype_index);

    int32_t num_written_info = 0;

    for (auto & transcript_path: _transcript_paths) {

        assert(transcript_path.is_reference || transcript_path.is_haplotype);

        // Exclude reference transcripts
        if (transcript_path.is_reference && !output_reference_transcripts) {

            continue;
        }

        ++num_written_info;

        // Get transcript path length.
        int32_t transcript_path_length = 0;

        for (auto & handle: transcript_path.path) {

            transcript_path_length += _graph->get_length(handle);
        }

        *tsv_ostream << transcript_path.get_name();
        *tsv_ostream << "\t" << transcript_path_length;
        *tsv_ostream << "\t";

        assert(!transcript_path.transcript_names.empty());

        bool is_first = true;

        for (auto & name: transcript_path.transcript_names) {

            if (!is_first) {

                *tsv_ostream << ",";
            } 

            is_first = false;
            *tsv_ostream << name;
        }

        *tsv_ostream << "\t";

        assert(!transcript_path.embedded_path_names.empty() || !transcript_path.haplotype_gbwt_ids.empty());
        
        is_first = true;

        for (auto & name: transcript_path.embedded_path_names) {

            if (!is_first) {

                *tsv_ostream << ",";
            } 

            is_first = false;
            *tsv_ostream << name;
        }

        is_first = true;

        for (auto & id: transcript_path.haplotype_gbwt_ids) {

            if (!is_first) {

                *tsv_ostream << ",";
            } 

            is_first = false;             
            PathSense sense = gbwtgraph::get_path_sense(haplotype_index, id, gbwt_reference_samples);
            *tsv_ostream << gbwtgraph::compose_path_name(haplotype_index, id, sense);
        }

        *tsv_ostream << endl;
    }

    return num_written_info;
}

void Transcriptome::write_graph(ostream * graph_ostream) const {

    vg::io::save_handle_graph(_graph.get(), *graph_ostream);
}

}



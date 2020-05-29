
#include <thread>

#include "../algorithms/topological_sort.hpp"
#include "../algorithms/apply_bulk_modifications.hpp"
#include "../io/save_handle_graph.hpp"

#include "transcriptome.hpp"
#include "../gbwt_helper.hpp"
#include "../augment.hpp"
#include "../utility.hpp"

namespace vg {

using namespace std;

//#define transcriptome_debug


Transcriptome::Transcriptome(const string & graph_filename, const bool show_progress) {

    // Load variation graph.
    get_input_file(graph_filename, [&](istream& in) {
        _splice_graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
    });

    _splice_graph_node_updated = false;

    if (!_splice_graph) {
        cerr << "[transcriptome] ERROR: Could not load graph." << endl;
        exit(1);
    }
}

int32_t Transcriptome::add_intron_splice_junctions(istream & intron_stream, unique_ptr<gbwt::GBWT> & haplotype_index) {

#ifdef transcriptome_debug
    double time_parsing_1 = gcsa::readTimer();
    cerr << "DEBUG parsing start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    // Parse introns in BED format.
    auto introns = parse_introns(intron_stream);

#ifdef transcriptome_debug
    cerr << "DEBUG parsing end: " << gcsa::readTimer() - time_parsing_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_project_1 = gcsa::readTimer();
    cerr << "DEBUG construct start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    // Construct edited transcript paths.
    auto edited_transcript_paths = construct_edited_transcript_paths(introns);

#ifdef transcriptome_debug
    cerr << "DEBUG construct end: " << gcsa::readTimer() - time_project_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_augment_1 = gcsa::readTimer();
    cerr << "DEBUG augment start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    if (has_novel_exon_boundaries(edited_transcript_paths, false)) {

        // Augment splice graph with new exon boundaries 
        // and splice-junction edges.
        augment_splice_graph(&edited_transcript_paths, haplotype_index, false);
    
    } else {

        // Augment splice graph with new splice-junction edges.
        add_splice_junction_edges(edited_transcript_paths);
    }

#ifdef transcriptome_debug
    cerr << "DEBUG augment end: " << gcsa::readTimer() - time_augment_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    return introns.size();
}

int32_t Transcriptome::add_transcript_splice_junctions(istream & transcript_stream, unique_ptr<gbwt::GBWT> & haplotype_index) {

#ifdef transcriptome_debug
    double time_parsing_1 = gcsa::readTimer();
    cerr << "DEBUG parsing start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    // Parse transcripts in gtf/gff3 format.
    auto transcripts = parse_transcripts(transcript_stream);

#ifdef transcriptome_debug
    cerr << "DEBUG parsing end: " << gcsa::readTimer() - time_parsing_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_project_1 = gcsa::readTimer();
    cerr << "DEBUG construct start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    // Construct edited transcript paths.
    auto edited_transcript_paths = construct_edited_transcript_paths(transcripts);

#ifdef transcriptome_debug
    cerr << "DEBUG construct end: " << gcsa::readTimer() - time_project_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_augment_1 = gcsa::readTimer();
    cerr << "DEBUG augment start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    if (has_novel_exon_boundaries(edited_transcript_paths, true)) {

        // Augment splice graph with new exon boundaries 
        // and splice-junction edges.
        augment_splice_graph(&edited_transcript_paths, haplotype_index, true);
    
    } else {

        // Augment splice graph with new splice-junction edges.
        add_splice_junction_edges(edited_transcript_paths);
    }

#ifdef transcriptome_debug
    cerr << "DEBUG augment end: " << gcsa::readTimer() - time_augment_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    return transcripts.size();
}

int32_t Transcriptome::add_transcripts(istream & transcript_stream, const gbwt::GBWT & haplotype_index) {

#ifdef transcriptome_debug
    double time_parsing_1 = gcsa::readTimer();
    cerr << "DEBUG parsing start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    // Parse transcripts in gtf/gff3 format.
    auto transcripts = parse_transcripts(transcript_stream);

#ifdef transcriptome_debug
    cerr << "DEBUG parsing end: " << gcsa::readTimer() - time_parsing_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_project_1 = gcsa::readTimer();
    cerr << "DEBUG project and add start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    // Save number of transcript paths before adding new.
    auto pre_num_transcript_paths = _transcript_paths.size();

    // Project and add transcripts to transcriptome.
    project_and_add_transcripts(transcripts, haplotype_index, mean_node_length());

    // Augment splice graph with new splice-junction edges.    
    add_splice_junction_edges(_transcript_paths);

#ifdef transcriptome_debug
    cerr << "DEBUG project and add end: " << gcsa::readTimer() - time_project_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    assert(_transcript_paths.size() >= pre_num_transcript_paths);
    return (_transcript_paths.size() - pre_num_transcript_paths);
}

vector<Transcript> Transcriptome::parse_introns(istream & intron_stream) const {

    vector<Transcript> introns;

    pair<string, PathIndex *> chrom_path_index("", nullptr);

    int32_t line_number = 0;

    string chrom;
    string pos;
    string end;
    string strand;

    while (intron_stream.good()) {

        line_number += 1;
        getline(intron_stream, chrom, '\t');

        // Skip header.
        if (chrom.empty() || chrom.front() == '#') {

            intron_stream.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }

        if (!_splice_graph->has_path(chrom)) {
        
            cerr << "[transcriptome] ERROR: Chromomsome path \"" << chrom << "\" not found in graph (line " << line_number << ")." << endl;
            exit(1);

        } else if (chrom_path_index.first != chrom) {

            delete chrom_path_index.second;
            chrom_path_index.first = chrom;

            // Construct path index for chromosome/contig.
            chrom_path_index.second = new PathIndex(*_splice_graph, chrom_path_index.first);
        }

        assert(chrom_path_index.second);

        // Parse start and end intron position and convert end to inclusive.
        assert(getline(intron_stream, pos, '\t'));
        int32_t spos = stoi(pos);
        
        assert(getline(intron_stream, end));
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
        introns.emplace_back(Transcript("", is_reverse, chrom));

        // Add intron boundaries as flanking exons to current "intron" transcript.
        add_exon(&(introns.back()), make_pair(spos - 1, spos - 1), *chrom_path_index.second);
        add_exon(&(introns.back()), make_pair(epos + 1, epos + 1), *chrom_path_index.second);
    }

    if (introns.empty()) {

        cerr << "[transcriptome] ERROR: No intron parsed" << endl;
        exit(1);        
    }

    delete chrom_path_index.second;

    return introns;
}

vector<Transcript> Transcriptome::parse_transcripts(istream & transcript_stream) const {

    vector<Transcript> transcripts;

    pair<string, PathIndex *> chrom_path_index("", nullptr);

    int32_t line_number = 0;

    string chrom;
    string feature;
    string pos;
    string strand;
    string attributes;

    smatch regex_id_match;

    // Regex used to extract transcript name/id from gtf file.
    regex regex_id_exp_gtf(transcript_tag + "\\s{1}\"?([^\"]*)\"?;?");

    // Regex used to extract transcript name/id from gff file.
    regex regex_id_exp_gff(transcript_tag + "={1}([^;]*);?");

    while (transcript_stream.good()) {

        line_number += 1;
        getline(transcript_stream, chrom, '\t');

        // Skip header.
        if (chrom.empty() || chrom.front() == '#') {

            transcript_stream.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }

        if (!_splice_graph->has_path(chrom)) {
        
            cerr << "[transcriptome] ERROR: Chromomsome path \"" << chrom << "\" not found in graph (line " << line_number << ")." << endl;
            exit(1);

        } else if (chrom_path_index.first != chrom) {

            delete chrom_path_index.second;
            chrom_path_index.first = chrom;

            // Construct path index for chromosome/contig.
            chrom_path_index.second = new PathIndex(*_splice_graph, chrom_path_index.first);
        }

        assert(chrom_path_index.second);

        transcript_stream.ignore(numeric_limits<streamsize>::max(), '\t');         
        assert(getline(transcript_stream, feature, '\t'));

        // Select only relevant feature types.
        if (feature != feature_type && !feature_type.empty()) {

            transcript_stream.ignore(numeric_limits<streamsize>::max(), '\n');  
            continue;
        }

        // Parse start and end exon position and convert to 0-base.
        assert(getline(transcript_stream, pos, '\t'));
        int32_t spos = stoi(pos) - 1;
        assert(getline(transcript_stream, pos, '\t'));
        int32_t epos = stoi(pos) - 1;

        assert(spos <= epos);

        // Skip score column.
        transcript_stream.ignore(numeric_limits<streamsize>::max(), '\t');  
        
        // Parse strand and set whether it is reverse.
        assert(getline(transcript_stream, strand, '\t'));
        assert(strand == "+" || strand == "-");
        bool is_reverse = (strand == "-") ? true : false;

        // Skip frame column.
        transcript_stream.ignore(numeric_limits<streamsize>::max(), '\t');  

        assert(getline(transcript_stream, attributes, '\n'));

        string transcript_id = "";

        // Get transcript name/id from gtf attribute column using regex.
        if (std::regex_search(attributes, regex_id_match, regex_id_exp_gtf)) {

            assert(regex_id_match.size() == 2);
            transcript_id = regex_id_match[1];
        }

        // Get transcript name/id from gff attribute column using regex.
        if (std::regex_search(attributes, regex_id_match, regex_id_exp_gff)) {

            assert(regex_id_match.size() == 2);
            transcript_id = regex_id_match[1];
        }

        if (transcript_id.empty()) {

            cerr << "[transcriptome] ERROR: Tag \"" << transcript_tag << "\" not found in attributes \"" << attributes << "\" (line " << line_number << ")." << endl;
            exit(1);
        }

        // Is this a new transcript.
        if (transcripts.empty()) {

            transcripts.emplace_back(Transcript(transcript_id, is_reverse, chrom));
        
        // Is this a new transcript.
        } else if (transcripts.back().name != transcript_id) {

            // Reorder reversed order exons.
            reorder_exons(&transcripts.back());
            transcripts.emplace_back(Transcript(transcript_id, is_reverse, chrom));
        }

        assert(transcripts.back().is_reverse == is_reverse);

        // Add exon to current transcript.
        add_exon(&(transcripts.back()), make_pair(spos, epos), *chrom_path_index.second);
    }

    if (transcripts.empty()) {

        cerr << "[transcriptome] ERROR: No transcripts parsed (remember to set feature type \"-y\")" << endl;
        exit(1);        
    }

    delete chrom_path_index.second;

    reorder_exons(&transcripts.back());

    return transcripts;
}

float Transcriptome::mean_node_length() const {

    return static_cast<float>(_splice_graph->get_total_length()) / _splice_graph->get_node_count();
}

void Transcriptome::add_exon(Transcript * transcript, const pair<int32_t, int32_t> & exon_pos, const PathIndex & chrom_path_index) const {

    transcript->exons.emplace_back(exon_pos);

    // Exon border positions (last position in upstream intron and 
    // first position in downstream intron).
    const pair<int32_t, int32_t> exon_border_pos = make_pair(exon_pos.first - 1, exon_pos.second + 1);

    // Find path positions (node start position and id) of exon node 
    // borders (start - 1 and end + 1) using path index.
    auto chrom_path_index_start_it = chrom_path_index.find_position(exon_border_pos.first);
    auto chrom_path_index_end_it = chrom_path_index.find_position(exon_border_pos.second);

    assert(chrom_path_index_start_it != chrom_path_index.end());
    assert(chrom_path_index_end_it != chrom_path_index.end());

    assert(chrom_path_index_start_it->first <= exon_border_pos.first);
    assert(chrom_path_index_end_it->first <= exon_border_pos.second);

    transcript->exon_border_nodes.emplace_back(Position(), Position());

    // Set node id of exon border boundaries.
    transcript->exon_border_nodes.back().first.set_node_id(chrom_path_index_start_it->second.node);
    transcript->exon_border_nodes.back().second.set_node_id(chrom_path_index_end_it->second.node);

    // Set node offset of exon border boundaries. 
    transcript->exon_border_nodes.back().first.set_offset(exon_border_pos.first - chrom_path_index_start_it->first);
    transcript->exon_border_nodes.back().second.set_offset(exon_border_pos.second  - chrom_path_index_end_it->first);

    // Set whether exon node border boundaries are reverse.
    transcript->exon_border_nodes.back().first.set_is_reverse(transcript->is_reverse);
    transcript->exon_border_nodes.back().second.set_is_reverse(transcript->is_reverse);
}

void Transcriptome::reorder_exons(Transcript * transcript) const {

    if (transcript->is_reverse) {

        // Is exons in reverse order.
        bool is_reverse_order = true;
        for (size_t i = 1; i < transcript->exons.size(); i++) {

            if (transcript->exons[i].second >= transcript->exons[i-1].first) { 

                is_reverse_order = false; 
            }
        }

        // Reverse if exons are in reverse order.
        if (is_reverse_order) { 

            reverse(transcript->exons.begin(), transcript->exons.end()); 
            reverse(transcript->exon_border_nodes.begin(), transcript->exon_border_nodes.end()); 
        }
    }
}

list<EditedTranscriptPath> Transcriptome::construct_edited_transcript_paths(const vector<Transcript> & transcripts) const {

    list<EditedTranscriptPath> edited_transcript_paths;
    mutex edited_transcript_paths_mutex;

    vector<thread> construction_threads;
    construction_threads.reserve(num_threads);

    // Spawn construction threads.
    for (size_t thread_idx = 0; thread_idx < num_threads; thread_idx++) {

        construction_threads.push_back(thread(&Transcriptome::construct_edited_transcript_paths_callback, this, &edited_transcript_paths, &edited_transcript_paths_mutex, thread_idx, ref(transcripts)));
    }

    // Join construction threads.   
    for (auto & thread: construction_threads) {
        
        thread.join();
    }

    return edited_transcript_paths;
}

void Transcriptome::construct_edited_transcript_paths_callback(list<EditedTranscriptPath> * edited_transcript_paths, mutex * edited_transcript_paths_mutex, const int32_t thread_idx, const vector<Transcript> & transcripts) const {

    list<EditedTranscriptPath> thread_edited_transcript_paths;

    int32_t transcripts_idx = thread_idx;

    while (transcripts_idx < transcripts.size()) {

        // Get next transcript belonging to current thread.
        const Transcript & transcript = transcripts.at(transcripts_idx);

        // Construct edited transcript paths.
        auto new_edited_transcript_paths = project_transcript_embedded(transcript, true);
        
        assert(new_edited_transcript_paths.size() == 1);
        assert(!new_edited_transcript_paths.front().reference_origin.empty());

        thread_edited_transcript_paths.splice(thread_edited_transcript_paths.end(), new_edited_transcript_paths);
        transcripts_idx += num_threads;
    }

    edited_transcript_paths_mutex->lock();
    edited_transcript_paths->splice(edited_transcript_paths->end(), thread_edited_transcript_paths);
    edited_transcript_paths_mutex->unlock();
}

void Transcriptome::project_and_add_transcripts(const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const float mean_node_length) {

    vector<thread> projection_threads;
    projection_threads.reserve(num_threads);

    // Spawn projection threads.
    for (size_t thread_idx = 0; thread_idx < num_threads; thread_idx++) {

        projection_threads.push_back(thread(&Transcriptome::project_and_add_transcripts_callback, this, thread_idx, ref(transcripts), ref(haplotype_index), mean_node_length));
    }

    // Join projection threads.   
    for (auto & thread: projection_threads) {
        
        thread.join();
    }
}

void Transcriptome::project_and_add_transcripts_callback(const int32_t thread_idx, const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const float mean_node_length) {

    list<CompletedTranscriptPath> thread_completed_transcript_paths;

    int32_t transcripts_idx = thread_idx;

    while (transcripts_idx < transcripts.size()) {

        // Get next transcript belonging to current thread.
        const Transcript & transcript = transcripts.at(transcripts_idx);

        list<CompletedTranscriptPath> completed_transcript_paths;

        if (!haplotype_index.empty()) { 

            // Project transcript onto haplotypes in GBWT index.
            completed_transcript_paths = construct_completed_transcript_paths(project_transcript_gbwt(transcript, haplotype_index, mean_node_length)); 
        }

        if (use_all_paths || use_reference_paths) { 

            // Project transcript onto embedded paths.
            auto new_completed_transcript_paths = construct_completed_transcript_paths(project_transcript_embedded(transcript, false));

            // Combine and optionally remove redundant paths.
            append_transcript_paths(&completed_transcript_paths, &new_completed_transcript_paths, collapse_transcript_paths);
        }

        auto thread_completed_transcript_paths_it = completed_transcript_paths.begin();
        int32_t transcript_path_idx = 1;

        while (thread_completed_transcript_paths_it != completed_transcript_paths.end()) {

            // Set transcript path name. The name contains the original transcript name/id 
            // and a unique index for each non-reference transcript copy.
            if (thread_completed_transcript_paths_it->reference_origin.empty()) {

                thread_completed_transcript_paths_it->name = thread_completed_transcript_paths_it->transcript_origin + "_" + to_string(transcript_path_idx);
                ++transcript_path_idx;

            } else {

                thread_completed_transcript_paths_it->name = thread_completed_transcript_paths_it->transcript_origin;                
            }

            ++thread_completed_transcript_paths_it;
        }

        thread_completed_transcript_paths.splice(thread_completed_transcript_paths.end(), completed_transcript_paths);
        transcripts_idx += num_threads;
    }

    // Add transcript paths to transcriptome.
    mutex_transcript_paths.lock();

    _transcript_paths.reserve(_transcript_paths.size() + thread_completed_transcript_paths.size());

    for (auto & transcript_path: thread_completed_transcript_paths) {

        _transcript_paths.emplace_back(move(transcript_path));
    }

    mutex_transcript_paths.unlock();
}

list<EditedTranscriptPath> Transcriptome::project_transcript_gbwt(const Transcript & cur_transcript, const gbwt::GBWT & haplotype_index, const float mean_node_length) const {

    assert(haplotype_index.bidirectional());

    list<EditedTranscriptPath> edited_transcript_paths;

    vector<pair<vector<exon_nodes_t>, thread_ids_t> > haplotypes;
    unordered_map<int32_t, pair<int32_t, int32_t> > haplotype_id_index;

    for (size_t exon_idx = 0; exon_idx < cur_transcript.exons.size(); ++exon_idx) {

        // Calculate expected number of nodes between exon start and end.
        const int32_t expected_length = ceil((cur_transcript.exons.at(exon_idx).second - cur_transcript.exons.at(exon_idx).first + 1) / mean_node_length);

        // Get all haplotypes in GBWT index between exon start and end border nodes (last position in upstream intron and
        // first position in downstream intron).
        auto exon_haplotypes = get_exon_haplotypes(cur_transcript.exon_border_nodes.at(exon_idx).first.node_id(), cur_transcript.exon_border_nodes.at(exon_idx).second.node_id(), haplotype_index, expected_length);

        if (haplotypes.empty()) {

            for (auto & exon_haplotype: exon_haplotypes) {

                haplotypes.emplace_back(vector<exon_nodes_t>(1, exon_haplotype.first), exon_haplotype.second);
                haplotypes.back().first.reserve(cur_transcript.exons.size());

                for (auto & haplotype_id: exon_haplotype.second) {

                    assert(haplotype_id_index.emplace(haplotype_id, make_pair(haplotypes.size() - 1, exon_idx + 1)).second);
                }
            }
            
        } else {

            for (auto & exon_haplotype: exon_haplotypes) {

                assert(!exon_haplotype.first.empty());
                unordered_map<int32_t, uint32_t> extended_haplotypes;

                for (auto & haplotype_id: exon_haplotype.second) {

                    auto haplotype_id_index_it = haplotype_id_index.find(haplotype_id);

                    if (haplotype_id_index_it == haplotype_id_index.end()) {

                        continue;         
                    }

                    if (exon_idx != haplotype_id_index_it->second.second) {

                        assert(haplotype_id_index_it->second.second < exon_idx);
                        haplotype_id_index.erase(haplotype_id_index_it);
                        continue;
                    }

                    haplotype_id_index_it->second.second++;
                    pair<vector<exon_nodes_t>, thread_ids_t> * cur_haplotype = &haplotypes.at(haplotype_id_index_it->second.first);

                    if (extended_haplotypes.find(haplotype_id_index_it->second.first) != extended_haplotypes.end()) {

                        assert(cur_haplotype->first.size() == exon_idx + 1);
                        haplotypes.at(extended_haplotypes.at(haplotype_id_index_it->second.first)).second.emplace_back(haplotype_id);
                        haplotype_id_index_it->second.first = extended_haplotypes.at(haplotype_id_index_it->second.first);                       
                        continue;
                    }

                    if (cur_haplotype->first.size() == exon_idx) {

                        cur_haplotype->first.emplace_back(exon_haplotype.first);
                        cur_haplotype->second = {haplotype_id};
                        assert(extended_haplotypes.emplace(haplotype_id_index_it->second.first, haplotype_id_index_it->second.first).second);
                    
                    } else if (cur_haplotype->first.size() == exon_idx + 1) {

                        haplotypes.emplace_back(vector<exon_nodes_t>(cur_haplotype->first.begin(), cur_haplotype->first.end() - 1), thread_ids_t(1, haplotype_id));
                        haplotypes.back().first.emplace_back(exon_haplotype.first);

                        assert(extended_haplotypes.emplace(haplotype_id_index_it->second.first, haplotypes.size() - 1).second);
                        haplotype_id_index_it->second.first = haplotypes.size() - 1;                
                    
                    } else {

                        haplotype_id_index.erase(haplotype_id_index_it);
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

        // Construct transcript path and set transcript origin name.
        edited_transcript_paths.emplace_back(cur_transcript.name);
        edited_transcript_paths.back().haplotype_origin_ids.reserve(haplotype.second.size());

        // Add haplotype names as origins.
        for (auto & thread_id: haplotype.second) {

            // Convert bidirectional path id before finding name. 
            edited_transcript_paths.back().haplotype_origin_ids.emplace_back(gbwt::Path::id(thread_id));
        }

        for (size_t exon_idx = 0; exon_idx < cur_transcript.exons.size(); ++exon_idx) {

            auto exon_border_start_node = cur_transcript.exon_border_nodes.at(exon_idx).first;
            auto exon_border_end_nodes = cur_transcript.exon_border_nodes.at(exon_idx).second;

            assert(gbwt::Node::id(haplotype.first.at(exon_idx).front()) == exon_border_start_node.node_id());
            assert(gbwt::Node::id(haplotype.first.at(exon_idx).back()) == exon_border_end_nodes.node_id());

            for (auto & exon_node: haplotype.first.at(exon_idx)) {

                assert(exon_node != gbwt::ENDMARKER);

                auto node_id = gbwt::Node::id(exon_node);
                auto node_length = _splice_graph->get_length(_splice_graph->get_handle(node_id, false));

                int32_t offset = 0;

                if (node_id == exon_border_start_node.node_id()) {

                    if (exon_border_start_node.offset() + 1 == node_length) {

                        assert(haplotype.first.at(exon_idx).size() > 1);
                        assert(node_id != exon_border_end_nodes.node_id());

                        continue;
                    
                    } else {

                        offset = exon_border_start_node.offset() + 1;
                    }
                }

                int32_t edit_length = node_length - offset;

                // Adjust end position from exon border (first position in downstream intron)
                // to last position in exon.
                if (node_id == exon_border_end_nodes.node_id()) {

                    if (exon_border_end_nodes.offset() == 0) {

                        break;

                    } else {

                        edit_length = exon_border_end_nodes.offset() - offset;
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
                reverse_complement_path_in_place(&(edited_transcript_paths.back().path), [&](size_t node_id) {return _splice_graph->get_length(_splice_graph->get_handle(node_id, false));});
            }

            // Copy paths if collapse of identical transcript paths is not wanted.
            if (!collapse_transcript_paths && edited_transcript_paths.back().haplotype_origin_ids.size() > 1) {

                auto all_haplotype_origin_ids = edited_transcript_paths.back().haplotype_origin_ids;
                edited_transcript_paths.back().haplotype_origin_ids = {edited_transcript_paths.back().haplotype_origin_ids.front()};

                // Create identical copies of all haplotype origins.
                for (size_t i = 1; i < all_haplotype_origin_ids.size(); ++i) {

                    edited_transcript_paths.emplace_back(edited_transcript_paths.back());
                    edited_transcript_paths.back().haplotype_origin_ids.front() = {all_haplotype_origin_ids.at(i)};
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
    unordered_set<int32_t> end_haplotype_ids;
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
            exon_haplotype_queue.pop();
            continue;            
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

list<EditedTranscriptPath> Transcriptome::project_transcript_embedded(const Transcript & cur_transcript, const bool reference_only) const {

    vector<unordered_map<path_handle_t, step_handle_t> > exon_start_node_path_steps;
    vector<unordered_map<path_handle_t, step_handle_t> > exon_end_node_path_steps;

    exon_start_node_path_steps.reserve(cur_transcript.exon_border_nodes.size());
    exon_end_node_path_steps.reserve(cur_transcript.exon_border_nodes.size());

    // Get embedded path ids and node mappings for all exon border nodes in transcript.
    for (auto & exon_node: cur_transcript.exon_border_nodes) {

        exon_start_node_path_steps.emplace_back(unordered_map<path_handle_t, step_handle_t>());
        _splice_graph->for_each_step_on_handle(_splice_graph->get_handle(exon_node.first.node_id(), false), [&](const step_handle_t & step) {
            assert(exon_start_node_path_steps.back().emplace(_splice_graph->get_path_handle_of_step(step), step).second);
        });

        exon_end_node_path_steps.emplace_back(unordered_map<path_handle_t, step_handle_t>());
        _splice_graph->for_each_step_on_handle(_splice_graph->get_handle(exon_node.second.node_id(), false), [&](const step_handle_t & step) {
            assert(exon_end_node_path_steps.back().emplace(_splice_graph->get_path_handle_of_step(step), step).second);
        });
    }

    list<EditedTranscriptPath> edited_transcript_paths;

    // Loop over all paths that contain the transcript start node.
    for (auto & path_steps_start: exon_start_node_path_steps.front()) {

        // Skip path if transcript end node is not in the current path.
        if (exon_end_node_path_steps.back().find(path_steps_start.first) == exon_end_node_path_steps.back().end()) {

            continue;
        }

        const auto path_origin_name = _splice_graph->get_path_name(path_steps_start.first);

        // Skip alternative allele paths (_alt).
        if (Paths::is_alt(path_origin_name)) {

            continue;
        }

        if (path_origin_name != cur_transcript.chrom) {

            // Construct only transcript paths originating from a reference chromosome/contig.        
            if (reference_only || !use_all_paths) {

                continue;
            }
        }

        // Construct transcript path and set transcript origin name.
        edited_transcript_paths.emplace_back(cur_transcript.name);

        // Does transcript path originate from a reference chromosome/contig.
        if (path_origin_name == cur_transcript.chrom) {

            edited_transcript_paths.back().reference_origin = path_origin_name;

        } else {

            if (edited_transcript_paths.back().path_origin_names.empty()) {

                edited_transcript_paths.back().path_origin_names = path_origin_name;

            } else {

                edited_transcript_paths.back().path_origin_names.append("," + path_origin_name);
            }
        }

        bool is_partial = false;

        for (size_t exon_idx = 0; exon_idx < exon_start_node_path_steps.size(); ++exon_idx) {

            auto haplotype_path_start_it = exon_start_node_path_steps.at(exon_idx).find(path_steps_start.first);
            auto haplotype_path_end_it = exon_end_node_path_steps.at(exon_idx).find(path_steps_start.first);

            // Transcript paths are partial if either the start or end exon path 
            // step is empty. Partial transcripts are currently not supported.
            // TODO: Add support for partial transcript paths.
            if ((haplotype_path_start_it == exon_start_node_path_steps.at(exon_idx).end()) || haplotype_path_end_it == exon_end_node_path_steps.at(exon_idx).end()) {

                is_partial = true;
                break;
            }

            // Get path step at exon start if exon start node is in the current path.
            auto haplotype_path_start_step = haplotype_path_start_it->second;

            // Get path mapping at exon end if exon end node is in the current path.
            auto haplotype_path_end_step = haplotype_path_end_it->second;

            bool is_path_origin_reverse = _splice_graph->get_is_reverse(_splice_graph->get_handle_of_step(haplotype_path_start_step));

            // Transcript paths are partial if the orientation of start 
            // and end exon path step is different. Partial transcripts 
            // are currently not supported.            
            // TODO: Add support for this really special case.
            if (is_path_origin_reverse != _splice_graph->get_is_reverse(_splice_graph->get_handle_of_step(haplotype_path_end_step))) {

                is_partial = true;
                break;
            }

            bool is_first_step = true;

            while (true) {

                auto node_length = _splice_graph->get_length(_splice_graph->get_handle_of_step(haplotype_path_start_step));
                int32_t offset = 0;

                // Adjust start position from exon border (last position in upstream intron)
                // to first position in exon.
                if (is_first_step) {

                    if (cur_transcript.exon_border_nodes.at(exon_idx).first.offset() + 1 == node_length) {

                        assert(haplotype_path_start_step != haplotype_path_end_step);

                        if (is_path_origin_reverse) {

                            haplotype_path_start_step = _splice_graph->get_previous_step(haplotype_path_start_step);

                        } else {

                            haplotype_path_start_step = _splice_graph->get_next_step(haplotype_path_start_step);
                        }

                        is_first_step = false;
                        continue;
                    
                    } else {

                        offset = cur_transcript.exon_border_nodes.at(exon_idx).first.offset() + 1;
                    }
                }

                int32_t edit_length = node_length - offset;

                // Adjust end position from exon border (first position in downstream intron)
                // to last position in exon.
                if (haplotype_path_start_step == haplotype_path_end_step) {

                    if (cur_transcript.exon_border_nodes.at(exon_idx).second.offset() == 0) {

                        break;

                    } else {

                        edit_length = cur_transcript.exon_border_nodes.at(exon_idx).second.offset() - offset;
                    }
                }

                assert(0 <= offset && offset < node_length);
                assert(0 < edit_length && edit_length <= node_length);

                // Add new mapping in forward direction. Later the whole path will
                // be reverse complemented if transcript is on the '-' strand.
                auto new_mapping = edited_transcript_paths.back().path.add_mapping();
                new_mapping->set_rank(edited_transcript_paths.back().path.mapping_size());

                new_mapping->mutable_position()->set_node_id(_splice_graph->get_id(_splice_graph->get_handle_of_step(haplotype_path_start_step)));
                new_mapping->mutable_position()->set_offset(offset);
                new_mapping->mutable_position()->set_is_reverse(false);

                // Add new edit representing a complete match.
                auto new_edit = new_mapping->add_edit();
                new_edit->set_from_length(edit_length);
                new_edit->set_to_length(edit_length);
                
                if (haplotype_path_start_step == haplotype_path_end_step) { break; }

                if (is_path_origin_reverse) {

                    haplotype_path_start_step = _splice_graph->get_previous_step(haplotype_path_start_step);

                } else {

                    haplotype_path_start_step = _splice_graph->get_next_step(haplotype_path_start_step);
                }

                is_first_step = false;
            }
        }

        if (is_partial) {

            // Delete partial transcript paths.
            edited_transcript_paths.pop_back();
        
        } else if (edited_transcript_paths.back().path.mapping_size() == 0) {

            // Delete empty paths.
            edited_transcript_paths.pop_back();

        } else {

            // Reverse complement transcript paths that are on the '-' strand.
            if (cur_transcript.is_reverse) {

                reverse_complement_path_in_place(&(edited_transcript_paths.back().path), [&](size_t node_id) {return _splice_graph->get_length(_splice_graph->get_handle(node_id, false));});
            } 
        }  
    } 

    return edited_transcript_paths;
}

void Transcriptome::append_transcript_paths(list<CompletedTranscriptPath> * completed_transcript_paths, list<CompletedTranscriptPath> * new_completed_transcript_paths, const bool add_unqiue_paths_only) const {

    // Add only non unique transcript paths.
    if (add_unqiue_paths_only) {

        for (auto & new_completed_transcript_path: *new_completed_transcript_paths) {

            bool new_path_unqiue = true;

            for (auto & completed_transcript_path: *completed_transcript_paths) {

                // Check if two paths are identical.
                if (completed_transcript_path.path == new_completed_transcript_path.path) {

                    if (completed_transcript_path.transcript_origin != new_completed_transcript_path.transcript_origin) {

                        cerr << "[transcriptome] WARNING: Different transcripts collaped (" << completed_transcript_path.transcript_origin << " & " << new_completed_transcript_path.transcript_origin << ")" << endl;
                    }

                    assert(completed_transcript_path.reference_origin == new_completed_transcript_path.reference_origin || completed_transcript_path.reference_origin.empty() || new_completed_transcript_path.reference_origin.empty());

                    // Merge reference origin name.
                    if (completed_transcript_path.reference_origin.empty()) {

                        completed_transcript_path.reference_origin = new_completed_transcript_path.reference_origin;
                    }

                    // Merge haplotype origin ids.
                    completed_transcript_path.haplotype_origin_ids.insert(completed_transcript_path.haplotype_origin_ids.end(), new_completed_transcript_path.haplotype_origin_ids.begin(), new_completed_transcript_path.haplotype_origin_ids.end());

                    // Merge embedded path origin names.
                    if (completed_transcript_path.path_origin_names.empty()) {

                        completed_transcript_path.path_origin_names = new_completed_transcript_path.path_origin_names;

                    } else if (!new_completed_transcript_path.path_origin_names.empty()) {

                        completed_transcript_path.path_origin_names.append("," + new_completed_transcript_path.path_origin_names);
                    }
                    
                    new_path_unqiue = false;
                    break;
                }
            }

            if (new_path_unqiue) {

                completed_transcript_paths->push_back(new_completed_transcript_path);
            } 
        }
    
    } else {

        completed_transcript_paths->splice(completed_transcript_paths->end(), *new_completed_transcript_paths);
    }
}

list<CompletedTranscriptPath> Transcriptome::construct_completed_transcript_paths(const list<EditedTranscriptPath> & edited_transcript_paths) const {

    list<CompletedTranscriptPath> completed_transcript_paths;

    for (auto & transcript_path: edited_transcript_paths) {

        completed_transcript_paths.emplace_back(transcript_path.transcript_origin);

        completed_transcript_paths.back().name = transcript_path.name;
        completed_transcript_paths.back().reference_origin = transcript_path.reference_origin;
        completed_transcript_paths.back().haplotype_origin_ids = transcript_path.haplotype_origin_ids;
        completed_transcript_paths.back().path_origin_names = transcript_path.path_origin_names;
        completed_transcript_paths.back().path = path_to_handles(transcript_path.path);
    }

    return completed_transcript_paths;     
}

vector<handle_t> Transcriptome::path_to_handles(const Path & path) const {

    vector<handle_t> handle_path;
    handle_path.reserve(path.mapping_size());

    for (auto mapping: path.mapping()) {

        auto handle = _splice_graph->get_handle(mapping.position().node_id(), mapping.position().is_reverse());

        // Check that the path only consist of whole nodes (complete).        
        assert(mapping.edit_size() == 1);
        assert(edit_is_match(mapping.edit(0)));
        assert(mapping.position().offset() == 0);
        assert(mapping.edit(0).from_length() == _splice_graph->get_length(handle));

        handle_path.emplace_back(handle);
    }

    return handle_path;
}

bool Transcriptome::has_novel_exon_boundaries(const list<EditedTranscriptPath> & edited_transcript_paths, const bool include_transcript_ends) const {

    for (auto & transcript_path: edited_transcript_paths) {

        for (size_t i = 0; i < transcript_path.path.mapping_size(); i++) {

            auto & cur_mapping = transcript_path.path.mapping(i);
            auto cur_handle = _splice_graph->get_handle(cur_mapping.position().node_id(), cur_mapping.position().is_reverse());
            
            assert(cur_mapping.edit_size() == 1);
            assert(edit_is_match(cur_mapping.edit(0)));

            // Do not check if left boundary of start exon is novel.
            if (!include_transcript_ends && i == 0) {

                if (cur_mapping.position().offset() + cur_mapping.edit(0).from_length() != _splice_graph->get_length(cur_handle)) {

                    return true;
                }

            // Do not check if right boundary of end exon is novel.
            } else if (!include_transcript_ends && i == transcript_path.path.mapping_size() - 1) {

                if (cur_mapping.position().offset() > 0) {

                    return true;
                }

            // Check if both boundaries are novel.
            } else if (cur_mapping.position().offset() > 0 || cur_mapping.edit(0).from_length() != _splice_graph->get_length(cur_handle)) {

                return true;
            }
        }
    }

    return false;
}

void Transcriptome::augment_splice_graph(list<EditedTranscriptPath> * edited_transcript_paths, unique_ptr<gbwt::GBWT> & haplotype_index, const bool break_at_transcript_ends) {

    assert(_transcript_paths.empty());
    _splice_graph_node_updated = true;

    // Move paths to data structure compatible with augment.
    vector<Path> edited_paths;
    edited_paths.reserve(edited_transcript_paths->size());

    for (auto & transcript_path: *edited_transcript_paths) {

        edited_paths.emplace_back(move(transcript_path.path));
    }

    if (haplotype_index->empty()) {

        // Augment splice graph with edited paths. 
        augment(static_cast<MutablePathMutableHandleGraph *>(_splice_graph.get()), edited_paths, nullptr, nullptr, false, break_at_transcript_ends);
      
    } else {

        vector<Translation> translations;

#ifdef transcriptome_debug
    double time_edit_1 = gcsa::readTimer();
    cerr << "\tDEBUG edit start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

        // Augment splice graph with edited paths. 
        augment(static_cast<MutablePathMutableHandleGraph *>(_splice_graph.get()), edited_paths, &translations, nullptr, false, break_at_transcript_ends);

#ifdef transcriptome_debug
    cerr << "\tDEBUG edit end: " << gcsa::readTimer() - time_edit_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

        // Update threads in gbwt index to match new augmented graph.
        update_haplotype_index(haplotype_index, translations);
    }
}

void Transcriptome::update_haplotype_index(unique_ptr<gbwt::GBWT> & haplotype_index, const vector<Translation> & translations) const {

#ifdef transcriptome_debug
    double time_translation_1 = gcsa::readTimer();
    cerr << "\tDEBUG translation start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    unordered_map<gbwt::node_type, vector<pair<int32_t, gbwt::node_type> > > translation_index;

    // Create translation index 
    for (auto & translation: translations) {

        assert(translation.from().mapping_size() == 1);
        assert(translation.to().mapping_size() == 1);

        auto & from_mapping = translation.from().mapping(0);
        auto & to_mapping = translation.to().mapping(0);

        // Only store changes
        if (!google::protobuf::util::MessageDifferencer::Equals(from_mapping, to_mapping)) {

            auto translation_index_it = translation_index.emplace(mapping_to_gbwt(from_mapping), vector<pair<int32_t, gbwt::node_type> >());
            translation_index_it.first->second.emplace_back(from_mapping.position().offset(), mapping_to_gbwt(to_mapping));
        }
    }

    // Sort translation index by offset
    for (auto & translation: translation_index) {

        sort(translation.second.begin(), translation.second.end());
    }

#ifdef transcriptome_debug
    cerr << "\tDEBUG translation end: " << gcsa::readTimer() - time_translation_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_update_1 = gcsa::readTimer();
    cerr << "\tDEBUG update start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    assert(haplotype_index->bidirectional());

    // Silence GBWT index construction. 
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT); 
    gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(_splice_graph->max_node_id(), true)));

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

            auto translation_index_it = translation_index.find(node);

            if (translation_index_it != translation_index.end()) {

                // First node id is the same (new node offset is 
                // larger than 0). 
                if (translation_index_it->second.front().first > 0) {

                    new_gbwt_threads.emplace_back(node);
                }

                // Add new nodes.
                for (auto & new_node: translation_index_it->second) {

                    new_gbwt_threads.emplace_back(new_node.second);
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
    
#ifdef transcriptome_debug
    cerr << "\tDEBUG update end: " << gcsa::readTimer() - time_update_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif     
}   

void Transcriptome::add_splice_junction_edges(const list<EditedTranscriptPath> & edited_transcript_paths) {

    for (auto & transcript_path: edited_transcript_paths) {

        for (size_t i = 1; i < transcript_path.path.mapping_size(); i++) {

            auto & prev_mapping = transcript_path.path.mapping(i - 1);
            auto & cur_mapping = transcript_path.path.mapping(i);

            auto prev_handle = _splice_graph->get_handle(prev_mapping.position().node_id(), prev_mapping.position().is_reverse());
            auto cur_handle = _splice_graph->get_handle(cur_mapping.position().node_id(), cur_mapping.position().is_reverse());
            
            // Ensure the edge exists.
            _splice_graph->create_edge(prev_handle, cur_handle);
        }
    }
}
  
void Transcriptome::add_splice_junction_edges(const vector<CompletedTranscriptPath> & completed_transcript_paths) {

    for (auto & transcript_path: completed_transcript_paths) {

        for (size_t i = 1; i < transcript_path.path.size(); i++) {
            
            // Ensure the edge exists.
            _splice_graph->create_edge(transcript_path.path.at(i - 1), transcript_path.path.at(i));
        }
    }
}

const vector<CompletedTranscriptPath> & Transcriptome::transcript_paths() const {

    return _transcript_paths;
}

int32_t Transcriptome::size() const {

    return _transcript_paths.size();
}

const MutablePathDeletableHandleGraph & Transcriptome::splice_graph() const {

    return *_splice_graph;
}

bool Transcriptome::splice_graph_node_updated() const {

    return _splice_graph_node_updated;
}

void Transcriptome::remove_non_transcribed(const bool new_reference_paths) {

    vector<path_handle_t> path_handles;
    path_handles.reserve(_splice_graph->get_path_count());

    assert(_splice_graph->for_each_path_handle([&](const path_handle_t & path_handle) {

        path_handles.emplace_back(path_handle);
    }));    

    // Remove all paths.
    for (auto & path_handle: path_handles) {

        _splice_graph->destroy_path(path_handle);
    }

    assert(_splice_graph->get_path_count() == 0);

    // Find all nodes that are in a transcript path.
    unordered_set<nid_t> transcribed_nodes;

    for (auto & transcript_path: _transcript_paths) {

        assert(transcript_path.path.size() > 0);
        for (auto & handle: transcript_path.path) {

            transcribed_nodes.emplace(_splice_graph->get_id(handle));
        }    
    }

    vector<handle_t> non_transcribed_handles;
    non_transcribed_handles.reserve(_splice_graph->get_node_count() - transcribed_nodes.size());

    // Collect all nodes that are not in a transcript path.
    assert(_splice_graph->for_each_handle([&](const handle_t & handle) {
        
        if (transcribed_nodes.count(_splice_graph->get_id(handle)) == 0) {

            non_transcribed_handles.emplace_back(handle); 
        } 
    }));

    for (auto & handle: non_transcribed_handles) {

        // Delete node and in/out edges. 
        _splice_graph->destroy_handle(handle);
    }

    assert(_splice_graph->get_node_count() == transcribed_nodes.size());
}

void Transcriptome::compact_ordered() {

    assert(_transcript_paths.empty());
    _splice_graph->apply_ordering(algorithms::topological_order(_splice_graph.get()), true);
}

int32_t Transcriptome::embed_transcript_paths(const bool add_reference_paths, const bool add_non_reference_paths) {

    int32_t num_embedded_paths = 0;

    // Add transcript paths to graph
    for (auto & transcript_path: _transcript_paths) {

        assert(!transcript_path.reference_origin.empty() || !transcript_path.haplotype_origin_ids.empty() || !transcript_path.path_origin_names.empty());

        if ((add_reference_paths && !transcript_path.reference_origin.empty()) || (add_non_reference_paths && (!transcript_path.haplotype_origin_ids.empty() || !transcript_path.path_origin_names.empty()))) {

            ++num_embedded_paths;

            assert(!_splice_graph->has_path(transcript_path.name));
            auto path_handle = _splice_graph->create_path_handle(transcript_path.name);

            for (auto & handle: transcript_path.path) {

                _splice_graph->append_step(path_handle, handle);
            }
        }
    }

    return num_embedded_paths;
}

int32_t Transcriptome::add_transcripts_to_gbwt(gbwt::GBWTBuilder * gbwt_builder, const bool output_reference_transcripts, const bool add_bidirectional) const {

    int32_t num_added_threads = 0;

    vector<string> sample_names;
    sample_names.reserve(size());

    assert(!gbwt_builder->index.hasMetadata());
    gbwt_builder->index.addMetadata();

    for (auto & transcript_path: _transcript_paths) {

        assert(!transcript_path.reference_origin.empty() || !transcript_path.haplotype_origin_ids.empty() || !transcript_path.path_origin_names.empty());

        if (output_reference_transcripts || !transcript_path.haplotype_origin_ids.empty() || !transcript_path.path_origin_names.empty()) {

            ++num_added_threads;

            // Convert transcript path to GBWT thread.
            gbwt::vector_type gbwt_thread(transcript_path.path.size());
            for (size_t i = 0; i < transcript_path.path.size(); i++) {

                gbwt_thread[i] = handle_to_gbwt(*_splice_graph, transcript_path.path.at(i));
            }

            // Insert transcript path as thread into GBWT index.
            gbwt_builder->insert(gbwt_thread, add_bidirectional);

            // Insert transcript path name into GBWT index.
            gbwt_builder->index.metadata.addPath({static_cast<gbwt::PathName::path_name_type>(sample_names.size()), 0, 0, 0});
            sample_names.emplace_back(transcript_path.name);
        }
    }

    // Set number number of haplotypes and transcript path name in metadata.
    gbwt_builder->index.metadata.setHaplotypes(size());
    gbwt_builder->index.metadata.setSamples(sample_names);

    return num_added_threads;
}

int32_t Transcriptome::write_sequences(ostream * fasta_ostream, const bool output_reference_transcripts) const {

    int32_t num_written_sequences = 0;

    for (auto & transcript_path: _transcript_paths) {

        assert(!transcript_path.reference_origin.empty() || !transcript_path.haplotype_origin_ids.empty() || !transcript_path.path_origin_names.empty());

        if (output_reference_transcripts || !transcript_path.haplotype_origin_ids.empty() || !transcript_path.path_origin_names.empty()) {

            ++num_written_sequences;

            // Construct transcript path sequence.
            string transcript_path_sequence = "";
            for (auto & handle: transcript_path.path) {

                transcript_path_sequence += _splice_graph->get_sequence(handle);
            }

            // Write transcript path name and sequence.
            write_fasta_sequence(transcript_path.name, transcript_path_sequence, *fasta_ostream);
        }
    }

    return num_written_sequences;
}

int32_t Transcriptome::write_info(ostream * tsv_ostream, const gbwt::GBWT & haplotype_index, const bool output_reference_transcripts) const {

    int32_t num_written_info = 0;

    *tsv_ostream << "Name\tLength\tTranscript\tReference\tHaplotypes" << endl; 

    for (auto & transcript_path: _transcript_paths) {

        assert(!transcript_path.reference_origin.empty() || !transcript_path.haplotype_origin_ids.empty() || !transcript_path.path_origin_names.empty());

        if (output_reference_transcripts || !transcript_path.haplotype_origin_ids.empty() || !transcript_path.path_origin_names.empty()) {

            ++num_written_info;

            // Get transcript path length.
            int32_t transcript_path_length = 0;
            for (auto & handle: transcript_path.path) {

                transcript_path_length += _splice_graph->get_length(handle);
            }

            *tsv_ostream << transcript_path.name;
            *tsv_ostream << "\t" << transcript_path_length;
            *tsv_ostream << "\t" << transcript_path.transcript_origin;

            if (transcript_path.reference_origin.empty()) {

                *tsv_ostream << "\t-";
            
            } else {

                *tsv_ostream << "\t" << transcript_path.reference_origin;
            }

            if (transcript_path.haplotype_origin_ids.empty() && transcript_path.path_origin_names.empty()) {

                *tsv_ostream << "\t-";            

            } else {

                *tsv_ostream << "\t";
                bool is_first = true;

                for (auto & id: transcript_path.haplotype_origin_ids) {

                    if (!is_first) {

                        *tsv_ostream << ",";
                    } 

                    is_first = false;             
                    *tsv_ostream << thread_name(haplotype_index, id);
                }

                if (!transcript_path.path_origin_names.empty()) {

                    if (!is_first) {

                        *tsv_ostream << ",";
                    } 

                    *tsv_ostream << transcript_path.path_origin_names;
                }
            }

            *tsv_ostream << endl;
        }
    }

    return num_written_info;
}

void Transcriptome::write_splice_graph(ostream * graph_ostream) const {

    vg::io::save_handle_graph(_splice_graph.get(), *graph_ostream);
}
    
}







#include <thread>

#include "../algorithms/topological_sort.hpp"
#include "../algorithms/apply_bulk_modifications.hpp"

#include "transcriptome.hpp"
#include "../gbwt_helper.hpp"

namespace vg {

using namespace std;

// #define transcriptome_debug


Transcriptome::Transcriptome(const string & graph_filename, const bool show_progress) {

    // Load variation graph.
    get_input_file(graph_filename, [&](istream& in) {
        _splice_graph = new VG(in, show_progress);
    });

    if (!_splice_graph) {
        cerr << "[transcriptome] ERROR: Could not load graph." << endl;
        exit(1);
    }
}

Transcriptome::~Transcriptome() {

    delete _splice_graph;
}

void Transcriptome::add_transcripts(istream & transcript_stream, const gbwt::GBWT & haplotype_index) {

#ifdef transcriptome_debug
    double time_parsing_1 = gcsa::readTimer();
    cerr << "DEBUG parsing start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    vector<Transcript> transcripts;

    // Get mean length of nodes in the graph.
    const float mean_node_length = _splice_graph->length() / static_cast<double>(_splice_graph->size());
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

        if (!_splice_graph->paths.has_path(chrom)) {
        
            cerr << "[transcriptome] ERROR: Chromomsome path \"" << chrom << "\" not found in graph (line " << line_number << ")." << endl;
            exit(1);

        } else if (chrom_path_index.first != chrom) {

            delete chrom_path_index.second;
            chrom_path_index.first = chrom;

            // Construct path index for chromosome/contig.
            chrom_path_index.second = new PathIndex(*_splice_graph, chrom);
        }

        assert(chrom_path_index.second);

        transcript_stream.ignore(numeric_limits<streamsize>::max(), '\t');         
        getline(transcript_stream, feature, '\t');

        // Select only relevant feature types.
        if (feature != feature_type && !feature_type.empty()) {

            transcript_stream.ignore(numeric_limits<streamsize>::max(), '\n');  
            continue;
        }

        // Parse start and end exon position and convert to 0-base.
        getline(transcript_stream, pos, '\t');
        int32_t spos = stoi(pos) - 1;
        getline(transcript_stream, pos, '\t');
        int32_t epos = stoi(pos) - 1;

        assert(spos <= epos);

        // Skip score column.
        transcript_stream.ignore(numeric_limits<streamsize>::max(), '\t');  
        
        // Parse strand and set whether it is reverse.
        getline(transcript_stream, strand, '\t');
        assert(strand == "+" || strand == "-");
        bool is_reverse = (strand == "-") ? true : false;

        // Skip frame column.
        transcript_stream.ignore(numeric_limits<streamsize>::max(), '\t');  

        getline(transcript_stream, attributes, '\n');

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

#ifdef transcriptome_debug
    double time_parsing_2 = gcsa::readTimer();
    cerr << "DEBUG parsing end: " << time_parsing_2 - time_parsing_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_project_1 = gcsa::readTimer();
    cerr << "DEBUG project start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    // Construct transcript paths from transcripts.
    auto proj_transcript_paths = project_transcripts(transcripts, haplotype_index, mean_node_length);

#ifdef transcriptome_debug
    double time_project_2 = gcsa::readTimer();
    cerr << "DEBUG project end: " << time_project_2 - time_project_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

#ifdef transcriptome_debug
    double time_add_1 = gcsa::readTimer();
    cerr << "DEBUG add start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 

    // Add projected transcript paths to transcriptome and
    // augment splice graph with new splice-junctions.
    add_paths_to_transcriptome(&proj_transcript_paths);

#ifdef transcriptome_debug
    double time_add_2 = gcsa::readTimer();
    cerr << "DEBUG add end: " << time_add_2 - time_add_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif 
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

list<TranscriptPath> Transcriptome::project_transcripts(const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const float mean_node_length) const {

    list<TranscriptPath> proj_transcript_paths;
    mutex proj_transcript_paths_mutex;

    vector<thread> projection_threads;
    projection_threads.reserve(num_threads);

    // Spawn projection threads.
    for (int32_t thread_idx = 0; thread_idx < num_threads; thread_idx++) {

        projection_threads.push_back(thread(&Transcriptome::project_transcripts_callback, this, &proj_transcript_paths, &proj_transcript_paths_mutex, thread_idx, ref(transcripts), ref(haplotype_index), mean_node_length));
    }

    // Join projection threads.   
    for (auto & thread: projection_threads) {
        
        thread.join();
    }

    return proj_transcript_paths;
}

void Transcriptome::project_transcripts_callback(list<TranscriptPath> * proj_transcript_paths, mutex * proj_transcript_paths_mutex, const int32_t thread_idx, const vector<Transcript> & transcripts, const gbwt::GBWT & haplotype_index, const float mean_node_length) const {

    list<TranscriptPath> thread_transcript_paths;

    int32_t transcripts_idx = thread_idx;

    while (transcripts_idx < transcripts.size()) {

        // Get next transcript belonging to current thread.
        const Transcript & transcript = transcripts.at(transcripts_idx);

        list<TranscriptPath> cur_transcript_paths;

        if (!haplotype_index.empty()) { 

            // Project transcript onto haplotypes in GBWT index.
            cur_transcript_paths = project_transcript_gbwt(transcript, haplotype_index, mean_node_length); 
        }

        if (use_embedded_paths || use_reference_paths) { 

            // Project transcript onto embedded paths and add new 
            // transcript paths to current set.
            auto new_transcript_paths = project_transcript_embedded(transcript);
            append_transcript_paths(&cur_transcript_paths, &new_transcript_paths, collapse_transcript_paths);
        }

        auto cur_transcript_paths_it = cur_transcript_paths.begin();
        int32_t transcript_path_idx = 1;

        while (cur_transcript_paths_it != cur_transcript_paths.end()) {

            // Set transcript path name. The name contains the original transcript name/id 
            // and a unique index for each copy of the transcript.
            cur_transcript_paths_it->name = cur_transcript_paths_it->transcript_origin + "_" + to_string(transcript_path_idx);
            ++transcript_path_idx;

            ++cur_transcript_paths_it;
        }

        thread_transcript_paths.splice(thread_transcript_paths.end(), cur_transcript_paths);
        transcripts_idx += num_threads;
    }

    // Add transcript paths to transcriptome.
    lock_guard<mutex> trancriptome_lock(*proj_transcript_paths_mutex);
    proj_transcript_paths->splice(proj_transcript_paths->end(), thread_transcript_paths);
}

list<TranscriptPath> Transcriptome::project_transcript_gbwt(const Transcript & cur_transcript, const gbwt::GBWT & haplotype_index, const float mean_node_length) const {

    list<TranscriptPath> cur_transcript_paths;

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
        cur_transcript_paths.emplace_back(cur_transcript.name);
        cur_transcript_paths.back().haplotype_origins.reserve(haplotype.second.size());

        // Add haplotype names as origins.
        for (auto & thread_id: haplotype.second) {

            // Convert bidirectional path id before finding name. 
            cur_transcript_paths.back().haplotype_origins.emplace_back(thread_name(haplotype_index, gbwt::Path::id(thread_id)));
        }

        for (size_t exon_idx = 0; exon_idx < cur_transcript.exons.size(); ++exon_idx) {

            auto exon_border_start_node = cur_transcript.exon_border_nodes.at(exon_idx).first;
            auto exon_border_end_nodes = cur_transcript.exon_border_nodes.at(exon_idx).second;

            assert(gbwt::Node::id(haplotype.first.at(exon_idx).front()) == exon_border_start_node.node_id());
            assert(gbwt::Node::id(haplotype.first.at(exon_idx).back()) == exon_border_end_nodes.node_id());

            for (auto & exon_node: haplotype.first.at(exon_idx)) {

                auto node_id = gbwt::Node::id(exon_node);
                auto node_length = _splice_graph->get_node(node_id)->sequence().size();

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
                auto new_mapping = cur_transcript_paths.back().path.add_mapping();
                new_mapping->set_rank(cur_transcript_paths.back().path.mapping_size());

                new_mapping->mutable_position()->set_node_id(node_id);
                new_mapping->mutable_position()->set_offset(offset);
                new_mapping->mutable_position()->set_is_reverse(false);

                // Add new edit representing a complete match.
                auto new_edit = new_mapping->add_edit();
                new_edit->set_from_length(edit_length);
                new_edit->set_to_length(edit_length);
            }
        }

        assert(cur_transcript_paths.back().path.mapping_size() > 0);

        if (cur_transcript.is_reverse) {

            // Reverse complement transcript paths that are on the '-' strand.
            reverse_complement_path_in_place(&(cur_transcript_paths.back().path), [&](size_t node_id) {return _splice_graph->get_node(node_id)->sequence().size();});
        }

        // Copy paths if collapse of identical transcript paths is not wanted.
        if (!collapse_transcript_paths && cur_transcript_paths.back().haplotype_origins.size() > 1) {

            auto all_haplotype_origins = cur_transcript_paths.back().haplotype_origins;
            cur_transcript_paths.back().haplotype_origins = {cur_transcript_paths.back().haplotype_origins.front()};

            // Create identical copies of all haplotype origins.
            for (size_t i = 1; i < all_haplotype_origins.size(); ++i) {

                cur_transcript_paths.emplace_back(cur_transcript_paths.back());
                cur_transcript_paths.back().haplotype_origins.front() = {all_haplotype_origins.at(i)};
            }
        }
    }

    return cur_transcript_paths; 
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

            auto extended_search = haplotype_index.extend(cur_exon_haplotype.second, out_edges_it->first);

            // Add new extension to queue if not empty (haplotypes found).
            if (!extended_search.empty()) { 

                exon_haplotype_queue.push(make_pair(cur_exon_haplotype.first, extended_search));
                exon_haplotype_queue.back().first.emplace_back(out_edges_it->first);
            }

            ++out_edges_it;
        }

        cur_exon_haplotype.first.emplace_back(out_edges.begin()->first);
        cur_exon_haplotype.second = haplotype_index.extend(cur_exon_haplotype.second, out_edges.begin()->first);        

        // End current extension if empty (no haplotypes found). 
        if (cur_exon_haplotype.second.empty()) { exon_haplotype_queue.pop(); }
    }

    return exon_haplotypes;
}

list<TranscriptPath> Transcriptome::project_transcript_embedded(const Transcript & cur_transcript) const {

    vector<map<int64_t, set<mapping_t*> > *> exon_start_node_mappings;
    vector<map<int64_t, set<mapping_t*> > *> exon_end_node_mappings;

    exon_start_node_mappings.reserve(cur_transcript.exon_border_nodes.size());
    exon_end_node_mappings.reserve(cur_transcript.exon_border_nodes.size());

    // Get embedded path ids and node mappings for all exon border nodes in transcript.
    for (auto & exon_node: cur_transcript.exon_border_nodes) {

        exon_start_node_mappings.emplace_back(&_splice_graph->paths.get_node_mapping(exon_node.first.node_id()));
        exon_end_node_mappings.emplace_back(&_splice_graph->paths.get_node_mapping(exon_node.second.node_id()));
    }

    list<TranscriptPath> cur_transcript_paths;

    // Loop over all paths that contain the transcript start node.
    for (auto & path_mapping_start: *exon_start_node_mappings.front()) {

        // Skip path if transcript end node is not in the current path.
        if (exon_end_node_mappings.back()->find(path_mapping_start.first) == exon_end_node_mappings.back()->end()) {

            continue;
        }

        // Skip alternative allele paths (_alt).
        if (Paths::is_alt(_splice_graph->paths.get_path_name(path_mapping_start.first))) {

            continue;
        }

        const auto path_origin_name = _splice_graph->paths.get_path_name(path_mapping_start.first);

        // Only construct transcript paths originating from a reference chromosome/contig.
        if (path_origin_name != cur_transcript.chrom && !use_embedded_paths && use_reference_paths) {

            continue;
        }

        // Construct transcript path and set transcript origin name.
        cur_transcript_paths.emplace_back(cur_transcript.name);

        // Does transcript path originate from a reference chromosome/contig.
        if (path_origin_name == cur_transcript.chrom) {

            cur_transcript_paths.back().reference_origin = path_origin_name;

        } else {

            cur_transcript_paths.back().haplotype_origins.emplace_back(path_origin_name);
        }

        bool is_partial = false;

        mapping_t * haplotype_path_start_map = nullptr;
        mapping_t * haplotype_path_end_map = nullptr;

        for (size_t exon_idx = 0; exon_idx < exon_start_node_mappings.size(); ++exon_idx) {

            auto haplotype_path_start_it = exon_start_node_mappings.at(exon_idx)->find(path_mapping_start.first);
            auto haplotype_path_end_it = exon_end_node_mappings.at(exon_idx)->find(path_mapping_start.first);

            // Get path mapping at exon start if exon start node is in the current path.
            if (haplotype_path_start_it != exon_start_node_mappings.at(exon_idx)->end()) {

                assert(haplotype_path_start_it->second.size() == 1);
                haplotype_path_start_map = *haplotype_path_start_it->second.begin();
            }

            // Get path mapping at exon end if exon end node is in the current path.
            if (haplotype_path_end_it != exon_end_node_mappings.at(exon_idx)->end()) {

                assert(haplotype_path_end_it->second.size() == 1);
                haplotype_path_end_map = *haplotype_path_end_it->second.begin();
            }

            // Transcript paths are partial if either the start or end exon path 
            // mapping is empty. Partial transcripts are currently not supported.
            // TODO: Add support for partial transcript paths.
            if (!haplotype_path_start_map || !haplotype_path_end_map) {

                is_partial = true;
                break;
            }

            bool is_first_mapping = true;

            while (true) {

                auto cur_node_id = haplotype_path_start_map->node_id();
                auto node_length = _splice_graph->get_node(cur_node_id)->sequence().size();
                assert(node_length == haplotype_path_start_map->length);

                int32_t offset = 0;

                // Adjust start position from exon border (last position in upstream intron)
                // to first position in exon.
                if (is_first_mapping) {

                    if (cur_transcript.exon_border_nodes.at(exon_idx).first.offset() + 1 == node_length) {

                        assert(haplotype_path_start_map != haplotype_path_end_map);

                        haplotype_path_start_map = _splice_graph->paths.traverse_right(haplotype_path_start_map);
                        assert(haplotype_path_start_map);

                        is_first_mapping = false;
                        continue;
                    
                    } else {

                        offset = cur_transcript.exon_border_nodes.at(exon_idx).first.offset() + 1;
                    }
                }

                int32_t edit_length = node_length - offset;

                // Adjust end position from exon border (first position in downstream intron)
                // to last position in exon.
                if (haplotype_path_start_map == haplotype_path_end_map) {

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
                auto new_mapping = cur_transcript_paths.back().path.add_mapping();
                new_mapping->set_rank(cur_transcript_paths.back().path.mapping_size());

                new_mapping->mutable_position()->set_node_id(cur_node_id);
                new_mapping->mutable_position()->set_offset(offset);
                new_mapping->mutable_position()->set_is_reverse(false);

                // Add new edit representing a complete match.
                auto new_edit = new_mapping->add_edit();
                new_edit->set_from_length(edit_length);
                new_edit->set_to_length(edit_length);
                
                if (haplotype_path_start_map == haplotype_path_end_map) { break; }

                haplotype_path_start_map = _splice_graph->paths.traverse_right(haplotype_path_start_map);
                assert(haplotype_path_start_map);
                
                is_first_mapping = false;
            }

            haplotype_path_start_map = nullptr;
            haplotype_path_end_map = nullptr;
        }

        if (is_partial) {

            // Delete partial transcript paths.
            cur_transcript_paths.pop_back();
        
        } else {

            assert(cur_transcript_paths.back().path.mapping_size() > 0);

            // Reverse complement transcript paths that are on the '-' strand.
            if (cur_transcript.is_reverse) {

                reverse_complement_path_in_place(&(cur_transcript_paths.back().path), [&](size_t node_id) {return _splice_graph->get_node(node_id)->sequence().size();});
            } 
        }  
    } 

    return cur_transcript_paths;
}

void Transcriptome::append_transcript_paths(list<TranscriptPath> * cur_transcript_paths, list<TranscriptPath> * new_transcript_paths, const bool add_unqiue_paths_only) const {

    // Add only non unique transcript paths.
    if (add_unqiue_paths_only) {

        for (auto & new_transcript_path: *new_transcript_paths) {

            bool new_path_unqiue = true;

            for (auto & cur_transcript_path: *cur_transcript_paths) {

                // Check if two path mappings are identical.
                if (cur_transcript_path.path.mapping_size() == new_transcript_path.path.mapping_size() && equal(cur_transcript_path.path.mapping().begin(), cur_transcript_path.path.mapping().end(), new_transcript_path.path.mapping().begin(), [](const Mapping & m1, const Mapping & m2) { return google::protobuf::util::MessageDifferencer::Equals(m1, m2); })) {

                    if (cur_transcript_path.transcript_origin != new_transcript_path.transcript_origin) {

                        cerr << "[transcriptome] WARNING: Different transcripts collaped (" << cur_transcript_path.transcript_origin << " & " << new_transcript_path.transcript_origin << ")" << endl;
                    }

                    assert(cur_transcript_path.reference_origin == new_transcript_path.reference_origin || cur_transcript_path.reference_origin.empty() || new_transcript_path.reference_origin.empty());

                    // Merge reference origin name.
                    if (cur_transcript_path.reference_origin.empty()) {

                        cur_transcript_path.reference_origin = new_transcript_path.reference_origin;
                    }

                    // Merge haplotype origin names.
                    cur_transcript_path.haplotype_origins.insert(cur_transcript_path.haplotype_origins.end(), new_transcript_path.haplotype_origins.begin(), new_transcript_path.haplotype_origins.end());
                    
                    new_path_unqiue = false;
                    break;
                }
            }

            if (new_path_unqiue) {

                cur_transcript_paths->push_back(new_transcript_path);
            } 
        }
    
    } else {

        cur_transcript_paths->splice(cur_transcript_paths->end(), *new_transcript_paths);
    }
}

bool Transcriptome::add_novel_transcript_junctions(const list<TranscriptPath> & cur_transcript_paths) {

    bool all_junctions_added = true;

    for (auto & transcript_path: cur_transcript_paths) {

        for (size_t i = 0; i < transcript_path.path.mapping_size(); i++) {

            auto & cur_mapping = transcript_path.path.mapping(i);

            if (cur_mapping.position().offset() > 0 || cur_mapping.edit_size() > 1 || !edit_is_match(cur_mapping.edit(0)) || !_splice_graph->has_node(cur_mapping.position().node_id()) || _splice_graph->get_node(cur_mapping.position().node_id())->sequence().size() != cur_mapping.edit(0).from_length()) {

                all_junctions_added = false;
                i++;
                continue;
            } 

            if (i > 0) {

                auto & prev_mapping = transcript_path.path.mapping(i - 1);

                auto prev_node_side = NodeSide(prev_mapping.position().node_id(), (prev_mapping.position().is_reverse() ? false : true));
                auto cur_node_side = NodeSide(cur_mapping.position().node_id(), (cur_mapping.position().is_reverse() ? true : false));
                
                // Ensure the edge exists.
                _splice_graph->create_edge(prev_node_side, cur_node_side);
            }
        }
    }

    return all_junctions_added;
}

void Transcriptome::add_paths_to_transcriptome(list<TranscriptPath> * new_transcript_paths) {

#ifdef transcriptome_debug
    double time_novel_1 = gcsa::readTimer();
    cerr << "DEBUG novel start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    // Add novel splice-junctions in transcript paths
    // to splice graph which does not require node splitting.
    bool all_junctions_added = add_novel_transcript_junctions(*new_transcript_paths);

#ifdef transcriptome_debug
    double time_novel_2 = gcsa::readTimer();
    cerr << "DEBUG novel end: " << time_novel_2 - time_novel_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

    // Augment splice graph with splice-junctions (including node splitting) 
    // in transcript paths.
    if (!all_junctions_added) {

        // Move paths to data structure compatible with edit.
        vector<Path> edit_paths;
        edit_paths.reserve(new_transcript_paths->size());

        for (auto & transcript_path: *new_transcript_paths) {

            edit_paths.emplace_back(move(transcript_path.path));
        }

#ifdef transcriptome_debug
    double time_edit_1 = gcsa::readTimer();
    cerr << "DEBUG edit start: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

        // Edit splice graph with projected transcript paths and
        // update path traversals to match the augmented graph. 
        _splice_graph->edit(edit_paths, nullptr, false, true, true);

#ifdef transcriptome_debug
    double time_edit_2 = gcsa::readTimer();
    cerr << "DEBUG edit end: " << time_edit_2 - time_edit_1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
#endif

        // Update projected transcript paths with new path traversals. 
        assert(edit_paths.size() == new_transcript_paths->size());

        auto new_transcript_paths_it = new_transcript_paths->begin();
        auto edit_paths_it = edit_paths.begin();

        while (new_transcript_paths_it != new_transcript_paths->end()) {

            new_transcript_paths_it->path = move(*edit_paths_it);

            ++new_transcript_paths_it;
            ++edit_paths_it;
        } 

        assert(edit_paths_it == edit_paths.end());
    }

    _transcript_paths.reserve(_transcript_paths.size() + new_transcript_paths->size());

    // Add projected transcript paths to transcriptome.
    for (auto & transcript_path: *new_transcript_paths) {

        _transcript_paths.emplace_back(move(transcript_path));
    }
}

const vector<TranscriptPath> & Transcriptome::transcript_paths() const {

    return _transcript_paths;
}

int32_t Transcriptome::size() const {

    return _transcript_paths.size();
}

const VG & Transcriptome::splice_graph() const {

    return *_splice_graph;
}

void Transcriptome::remove_non_transcribed(const bool new_reference_paths) {

    // Save copy of embedded reference paths
    Paths reference_paths;
    if (new_reference_paths) {

        reference_paths = _splice_graph->paths;
    }

    // Remove non transcript paths.
    _splice_graph->clear_paths();

    // Find all nodes and edges that are in a transcript path.
    unordered_set<vg::id_t> transcribed_nodes;
    unordered_set<pair<vg::id_t, vg::id_t> > transcribed_edges;

    for (auto & transcript_path: _transcript_paths) {

        assert(transcript_path.path.mapping_size() > 0);
        transcribed_nodes.emplace(transcript_path.path.mapping(0).position().node_id());

        for (size_t i = 1; i < transcript_path.path.mapping_size(); i++) {

            transcribed_nodes.emplace(transcript_path.path.mapping(i).position().node_id());
            transcribed_edges.emplace(transcript_path.path.mapping(i-1).position().node_id(), transcript_path.path.mapping(i).position().node_id());
        }    
    }

    // Find all nodes that are not in a transcript path.
    vector<vg::id_t> non_transcribed_nodes;
    _splice_graph->for_each_node([&](const Node * node) {
        
        if (transcribed_nodes.count(node->id()) == 0) {

             non_transcribed_nodes.emplace_back(node->id());   
        }
    });

    for (auto & node: non_transcribed_nodes) {

        // Delete node and in/out edges. 
        _splice_graph->destroy_node(node);
    }

    // Create new reference paths that only include trancribed nodes and edges.
    if (new_reference_paths) {

        reference_paths.for_each([&](const Path & path) {

            if (!Paths::is_alt(path.name())) {

                vector<Path> new_paths; 

                new_paths.emplace_back(Path());
                new_paths.back().set_name(path.name() + "_" + to_string(new_paths.size() - 1));

                for (auto & mapping: path.mapping()) {

                    auto cur_node_id = mapping.position().node_id();

                    if (new_paths.back().mapping_size() == 0) {

                        if (transcribed_nodes.count(cur_node_id) > 0) {

                            auto new_mapping = new_paths.back().add_mapping();
                            *new_mapping = mapping;
                            new_mapping->set_rank(new_paths.back().mapping_size()); 
                        }               
                    
                    } else {

                        auto prev_node_id = new_paths.back().mapping(new_paths.back().mapping_size() - 1).position().node_id();

                        // Extend new path, if transcribed edge (forward or reverse) exist between 
                        // this and the previous node in the path.
                        if (transcribed_edges.count(make_pair(prev_node_id, cur_node_id)) > 0 || transcribed_edges.count(make_pair(cur_node_id, prev_node_id)) > 0) {

                            auto new_mapping = new_paths.back().add_mapping();
                            *new_mapping = mapping;
                            new_mapping->set_rank(new_paths.back().mapping_size());   

                        } else {

                            new_paths.emplace_back(Path());
                            new_paths.back().set_name(path.name() + "_" + to_string(new_paths.size() - 1));
                        }
                    }
                }

                // Add new reference paths to graph without rebuild paths indexes.
                _splice_graph->paths.extend(new_paths, false, false); 
            }
        });

        // Rebuild paths indexes.
        _splice_graph->paths.compact_ranks();
    }
}

void Transcriptome::compact_ordered() {

    // Find and apply topological ordering 
    auto topological_ordering = algorithms::topological_order(_splice_graph);
    _splice_graph->apply_ordering(topological_ordering);

    // Compact node ids and update embedded paths.
    hash_map<id_t, id_t> compacted_nodes;
    _splice_graph->compact_ids(compacted_nodes);

    // Update transcript paths with compacted node ids
    for (auto & transcript_path: _transcript_paths) {
    
        for (auto & mapping: *transcript_path.path.mutable_mapping()) {

            mapping.mutable_position()->set_node_id(compacted_nodes.at(mapping.position().node_id()));
        }
    }
}

void Transcriptome::embed_transcript_paths(const bool add_reference_paths, const bool add_non_reference_paths, const bool rebuild_indexes) {

    // Add transcript paths to graph
    for (auto & transcript_path: _transcript_paths) {

        assert(!transcript_path.haplotype_origins.empty() || !transcript_path.reference_origin.empty());
        transcript_path.path.set_name(transcript_path.name);

        if (add_reference_paths && !transcript_path.reference_origin.empty()) {

            _splice_graph->paths.extend(transcript_path.path, false, false);
        
        } else if (add_non_reference_paths && !transcript_path.haplotype_origins.empty()) {

            _splice_graph->paths.extend(transcript_path.path, false, false);            
        }   

        transcript_path.path.set_name("");
    }

    // Rebuild paths indexes.
    if (rebuild_indexes) {

        _splice_graph->paths.compact_ranks();
    }
}

void Transcriptome::construct_gbwt(gbwt::GBWTBuilder * gbwt_builder, const bool output_reference_transcripts, const bool add_bidirectional) const {

    vector<string> sample_names;
    sample_names.reserve(size());

    assert(!gbwt_builder->index.hasMetadata());
    gbwt_builder->index.addMetadata();

    for (auto & transcript_path: _transcript_paths) {

        assert(!transcript_path.haplotype_origins.empty() || !transcript_path.reference_origin.empty());
        if (!transcript_path.haplotype_origins.empty() || output_reference_transcripts) {

            // Convert transcript path to GBWT thread.
            gbwt::vector_type gbwt_thread(transcript_path.path.mapping_size());
            for (size_t i = 0; i < transcript_path.path.mapping_size(); i++) {

                assert(transcript_path.path.mapping(i).edit_size() == 1);
                gbwt_thread[i] = mapping_to_gbwt(transcript_path.path.mapping(i));
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
}

void Transcriptome::write_alignments(ostream * gam_ostream, const bool output_reference_transcripts) const {

    vg::io::ProtobufEmitter<Alignment> emitter(*gam_ostream);

    for (auto & transcript_path: _transcript_paths) {

        assert(!transcript_path.haplotype_origins.empty() || !transcript_path.reference_origin.empty());
        if (!transcript_path.haplotype_origins.empty() || output_reference_transcripts) {

            // Write transcript path as alignment 
            Alignment alignment;
            alignment.set_name(transcript_path.name);
            *alignment.mutable_path() = transcript_path.path;
            emitter.write(std::move(alignment));
        }
    }
}

void Transcriptome::write_sequences(ostream * fasta_ostream, const bool output_reference_transcripts) {

    for (auto & transcript_path: _transcript_paths) {

        assert(!transcript_path.haplotype_origins.empty() || !transcript_path.reference_origin.empty());
        if (!transcript_path.haplotype_origins.empty() || output_reference_transcripts) {

            // Write transcript path name and sequence.
            *fasta_ostream << ">" << transcript_path.name << endl;
            *fasta_ostream << _splice_graph->path_sequence(transcript_path.path) << endl;
        }
    }
}

void Transcriptome::write_info(ostream * tsv_ostream, const bool output_reference_transcripts) const {

    *tsv_ostream << "Name\tLength\tTranscript\tReference\tHaplotypes" << endl; 

    for (auto & transcript_path: _transcript_paths) {

        assert(!transcript_path.haplotype_origins.empty() || !transcript_path.reference_origin.empty());
        if (!transcript_path.haplotype_origins.empty() || output_reference_transcripts) {

            *tsv_ostream << transcript_path.name;
            *tsv_ostream << "\t" << path_to_length(transcript_path.path);
            *tsv_ostream << "\t" << transcript_path.transcript_origin;

            if (transcript_path.reference_origin.empty()) {

                *tsv_ostream << "\t-";
            
            } else {

                *tsv_ostream << "\t" << transcript_path.reference_origin;
            }

            if (transcript_path.haplotype_origins.empty()) {

                *tsv_ostream << "\t-";            

            } else {

                auto haplotype_origins_it = transcript_path.haplotype_origins.begin();

                *tsv_ostream << "\t" << *haplotype_origins_it;
                ++haplotype_origins_it;

                while (haplotype_origins_it != transcript_path.haplotype_origins.end()) {

                    *tsv_ostream << "," << *haplotype_origins_it;
                    ++haplotype_origins_it;
                }
            }

            *tsv_ostream << endl;
        }
    }
}

void Transcriptome::write_splice_graph(ostream * vg_ostream) {

    _splice_graph->serialize_to_ostream(*vg_ostream);
}
    
}






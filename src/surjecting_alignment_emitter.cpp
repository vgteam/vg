/**
 * \file surjecting_alignment_emitter.cpp
 * Implementation for SurjectingAlignmentEmitter
 */


#include "surjecting_alignment_emitter.hpp"
#include "hts_alignment_emitter.hpp"

#include <map>

namespace vg {

using namespace std;

unique_ptr<AlignmentEmitter> get_alignment_emitter_with_surjection(const string& filename, const string& format, 
                                                                   const vector<path_handle_t> paths_in_dict_order, size_t max_threads,
                                                                   const HandleGraph* graph) {
    
    // Build a path length map by name
    // TODO: be able to pass along dict order!
    vector<pair<string, int64_t>> path_names_and_lengths;
    if (format == "SAM" || format == "BAM" || format == "CRAM") {
        // Make sure we actually have a PathPositionalHandleGraph
        const PathPositionHandleGraph* path_graph = dynamic_cast<const PathPositionHandleGraph*>(graph);
        if (path_graph == nullptr) {
            cerr << "error[vg::get_alignment_emitter_with_surjection]: Graph does not contain paths to surject into." << endl;
            exit(1);
        }
    
        // We will actually use it, so fill it in.
        for (auto& path_handle : paths_in_dict_order) {
            path_names_and_lengths.emplace_back(path_graph->get_path_name(path_handle), path_graph->get_path_length(path_handle));
        }
    }
    
    // Get the non-surjecting emitter
    unique_ptr<AlignmentEmitter> emitter = get_alignment_emitter(filename, format, path_names_and_lengths, max_threads, graph);
    
    if (format == "SAM" || format == "BAM" || format == "CRAM") {
        // Need to surject
        
        // Make sure (again) we actually have a PathPositionalHandleGraph
        const PathPositionHandleGraph* path_graph = dynamic_cast<const PathPositionHandleGraph*>(graph);
        assert(path_graph != nullptr);
        
        // Make a set of the path handles to surject into
        unordered_set<path_handle_t> target_paths(paths_in_dict_order.begin(), paths_in_dict_order.end());
        // Interpose a surjecting AlignmentEmitter
        emitter = make_unique<SurjectingAlignmentEmitter>(path_graph, target_paths, std::move(emitter));
    }
    
    return emitter;
}

vector<handle_t> get_sequence_dictionary_handles(const string& filename, const PathPositionHandleGraph& graph) {
    
    assert(graph != nullptr);
    
    // We fill in the "dictionary" (which is what SAM calls it; it's not a mapping for us)
    vector<path_handle_t> dictionary;
    
    if (!filename.empty()) {
        // TODO: As of right now HTSLib doesn't really let you iterate the sequence dictionary when you use its parser. So we use our own parser.
        get_input_file(filename, [&](istream& in) {
            for (string line; getline(in, line);) {
                // Each line will produce a sequence name and a handle
                string sequence_name = "";
                path_handle_t path;
            
                if (line.size() == 0) {
                    // Unless it is empty
                    continue;
                }
            
                // See if each line starts with @SQ or if we have to handle it as a name.
                if (starts_with(line, "@SQ")) {
                    // If it is SAM, split on tabs
                    auto parts = split_delims(line, "\t");
                    for (size_t i = 1; i < parts.size(); i++) {
                        if (starts_with(parts[i], "SN:")) {
                            // The rest of this field is the name
                            sequence_name = parts[i].substr(3);
                        } else if (starts_with(parts[i], "LN:")) {
                            // The rest of this field is a length number
                            length = stoll(parts[i].substr(3));
                        }
                    }
                    
                    if (sequence_name == "") {
                        cerr << "error:[vg::get_sequence_dictionary_handles] No sequence name for @SQ line " << line << endl;
                        exit(1);
                    }
                    if (length < 0) {
                        cerr << "error:[vg::get_sequence_dictionary_handles] Unacceptable sequence length " << length << " for sequence " << sequence_name << endl;
                        exit(1);
                    }
                    
                    // Check the sequence against the graph
                    if (!graph.has_path(sequence_name)) {
                        // Name doesn't exist
                        cerr << "error:[vg mpmap] Graph does not have a path named " << line << ", which was indicated in " << filename << endl;
                        exit(1);
                    }
                    path = graph.get_path_handle(sequence_name);
                    size_t graph_path_length = graph.get_path_length(path);
                    if (graph_path_length != length) {
                        // Length doesn't match
                        cerr << "error:[vg mpmap] Graph contains a path " << sequence_name << " of length " << graph_path_length
                            << " but sequence dictionary in " << filename << " indicates a length of " << length << endl;
                        exit(1);
                    }
                    
                } else {
                    // Get the name from the line and the sequence from the graph
                    sequence_name = line;
                    if (graph == nullptr) {
                        cerr << "error:[vg::get_sequence_dictionary_handles] No graph available to get length for sequence " << sequence_name << endl;
                        exit(1);
                    }
                    if (!graph.has_path(sequence_name)) {
                        cerr << "error:[vg mpmap] Graph does not have a path named " << line << ", which was indicated in " << filename << endl;
                        exit(1);
                    }
                    path = graph.get_path_handle(sequence_name);
                }
                
                // Save the path handle
                dictionary.push_back(handle);
            }
        });
        
        if (dictionary.empty()) {
            // There were no entries in the file
            cerr << "error:[vg::get_sequence_dictionary_handles] No sequence dictionary available in file: " << filename << endl;
            exit(1);
        }
    } else {
        // We need to look for non-alt paths in the graph
        if (graph == nullptr) {
            cerr << "error:[vg::get_sequence_dictionary_handles] Graph does not have embedded paths to treat as reference sequences. Cannot get sequence dictionary to produce HTSlib output formats (SAM/BAM/CRAM)." << endl;
            exit(1);
        }
        
        graph.for_each_path_handle([&](const path_handle_t& path_handle) {
            string sequence_name = graph.get_path_name(path_handle);
            if (!Paths::is_alt(sequence_name)) {
                // This isn't an alt allele path, so we want it.
                dictionary.push_back(path_handle);
            }
        });
        
        if (dictionary.empty()) {
            cerr << "error:[vg::get_sequence_dictionary_handles] No non-alt-allele paths available in the graph!" << endl;
            exit(1);
        }
    }
    
    return dictionary;
}

SurjectingAlignmentEmitter::SurjectingAlignmentEmitter(const PathPositionHandleGraph* graph, unordered_set<path_handle_t> paths,
    unique_ptr<AlignmentEmitter>&& backing) : surjector(graph), paths(paths), backing(std::move(backing)) {
    
    // Nothing to do!
    
}

void SurjectingAlignmentEmitter::surject_alignments_in_place(vector<Alignment>& alns) const {
    for (auto& aln : alns) {
        // Surject each alignment and annotate with surjected path position
        aln = surjector.surject(aln, paths, surject_subpath_global);
    }
}

void SurjectingAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln_batch_caught(aln_batch);
    // Surject it in place
    surject_alignments_in_place(aln_batch_caught);
    // Forward it along
    backing->emit_singles(std::move(aln_batch_caught));
}

void SurjectingAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns_batch_caught(alns_batch);
    for (auto& mappings : alns_batch_caught) {
        // Surject all mappings in place
        surject_alignments_in_place(mappings);
    }
    // Forward it along
    backing->emit_mapped_singles(std::move(alns_batch_caught));
}

void SurjectingAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch, vector<Alignment>&& aln2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<Alignment> aln1_batch_caught(aln1_batch);
    vector<Alignment> aln2_batch_caught(aln2_batch);
    // Surject it in place
    surject_alignments_in_place(aln1_batch_caught);
    surject_alignments_in_place(aln2_batch_caught);
    // Forward it along
    backing->emit_pairs(std::move(aln1_batch_caught), std::move(aln2_batch_caught), std::move(tlen_limit_batch));
}

void SurjectingAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch, vector<vector<Alignment>>&& alns2_batch, vector<int64_t>&& tlen_limit_batch) {
    // Intercept the batch on its way
    vector<vector<Alignment>> alns1_batch_caught(alns1_batch);
    vector<vector<Alignment>> alns2_batch_caught(alns2_batch);
    for (auto& mappings : alns1_batch_caught) {
        // Surject all mappings in place
        surject_alignments_in_place(mappings);
    }
    for (auto& mappings : alns2_batch_caught) {
        // Surject all mappings in place
        surject_alignments_in_place(mappings);
    }
    // Forward it along
    backing->emit_mapped_pairs(std::move(alns1_batch_caught), std::move(alns2_batch_caught), std::move(tlen_limit_batch));
}

}

/**
 * \file alignment_emitter.cpp
 *
 * Implements a system for emitting alignments and groups of alignments in multiple formats.
 */

#include "hts_alignment_emitter.hpp"
#include "surjecting_alignment_emitter.hpp"
#include "back_translating_alignment_emitter.hpp"
#include "alignment.hpp"
#include "vg/io/json2pb.h"
#include "algorithms/find_translation.hpp"
#include <vg/io/hfile_cppstream.hpp>
#include <vg/io/stream.hpp>

#include <sstream>

//#define debug

namespace vg {
using namespace std;

unique_ptr<AlignmentEmitter> get_alignment_emitter(const string& filename, const string& format,
                                                   const vector<tuple<path_handle_t, size_t, size_t>>& paths, size_t max_threads,
                                                   const HandleGraph* graph, int flags) {

    
    unique_ptr<AlignmentEmitter> emitter;
    
    if (format == "SAM" || format == "BAM" || format == "CRAM") {
        // We are doing linear HTSLib output
        
        // Make sure we actually have a PathPositionalHandleGraph
        const PathPositionHandleGraph* path_graph = dynamic_cast<const PathPositionHandleGraph*>(graph);
        if (path_graph == nullptr) {
            cerr << "error[vg::get_alignment_emitter]: No graph available supporting path length queries needed for " << format << " output." << endl;
            exit(1);
        }
        
        // Build a path name and length list from the handles (this is for the sequence dictionary and reflects *base* paths
        // as opposed to subpaths in the graph -- if there are no subpaths there is no distinction)
        vector<pair<string, int64_t>> path_names_and_lengths;
        // Remember the actual path lengths (this is for coordinate transformations)        
        unordered_map<string, int64_t> subpath_to_length;
        std::tie(path_names_and_lengths, subpath_to_length) = extract_path_metadata(paths, *path_graph, true);
    
        if (flags & ALIGNMENT_EMITTER_FLAG_HTS_SPLICED) {
            // Use a splicing emitter as the final emitter
            emitter = make_unique<SplicedHTSAlignmentEmitter>(filename, format, path_names_and_lengths, subpath_to_length, *path_graph, max_threads);
        } else {
            // Use a normal emitter
            emitter = make_unique<HTSAlignmentEmitter>(filename, format, path_names_and_lengths, subpath_to_length, max_threads);
        }
        
        if (!(flags & ALIGNMENT_EMITTER_FLAG_HTS_RAW)) {
            // Need to surject
            
            // Make a set of the path handles to surject into
            unordered_set<path_handle_t> target_paths;
            for (const auto& path_info : paths) {
                target_paths.insert(get<0>(path_info));
            }
            // Interpose a surjecting AlignmentEmitter
            emitter = make_unique<SurjectingAlignmentEmitter>(path_graph, target_paths, std::move(emitter),
                flags & ALIGNMENT_EMITTER_FLAG_HTS_PRUNE_SUSPICIOUS_ANCHORS);
        }
    
    } else {
        // The non-HTSlib formats don't actually use the path name and length info.
        // See https://github.com/vgteam/libvgio/issues/34
        
        const NamedNodeBackTranslation* translation = nullptr;
        if (flags & ALIGNMENT_EMITTER_FLAG_VG_USE_SEGMENT_NAMES) {
            // Need to translate from node IDs to segment names
            
            translation = vg::algorithms::find_translation(graph);
            if (translation == nullptr) {
                cerr << "error[vg::get_alignment_emitter]: No graph available supporting translation to named-segment space" << endl;
                exit(1);
            }
        }
        
        // TODO: Push some logic here into libvgio? Or move this top function out of hts_alignment_emitter.cpp?
        // TODO: Only GAF actually handles the translation in the emitter right now.
        // TODO: Move BackTranslatingAlignmentEmitter to libvgio so they all can and we don't have to sniff format here.
        emitter = get_non_hts_alignment_emitter(filename, format, {}, max_threads, graph, translation);
        if (translation && format != "GAF") {
            // Need to translate from node IDs to segment names beforehand.
            // Interpose a translating AlignmentEmitter
            emitter = make_unique<BackTranslatingAlignmentEmitter>(translation, std::move(emitter));
        }
    }
    
    return emitter;
}

/// Run the given iteratee for each path that is either the path with the given
/// name (if present), or a subrange of a path with the given name as the base
/// name (otherwise).
///
/// If a path and subpaths both exist, only look at the full path.
///
/// If the name describes a subpath, look only at that subpath.
///
/// Iteratee returns false to stop.
///
/// Returns true if we reached the end, and false if asked to stop.
static bool for_each_subpath_of(const PathPositionHandleGraph& graph, const string& path_name, const std::function<bool(const path_handle_t& path)>& iteratee) {
    if (graph.has_path(path_name)) {
        // Just look at the full path.
        return iteratee(graph.get_path_handle(path_name));
    }
    
    // Parse out the metadata of the thing we want subpaths of
    PathSense& sense;
    string sample;
    string locus;
    size_t haplotype;
    size_t phase_block;
    subrange_t subrange;
    PathMetadata::parse_path_name(path_name,
                                  sense,
                                  sample,
                                  locus,
                                  haplotype,
                                  phase_block,
                                  subrange);
                                  
    if (subrange != PathMetadata::NO_SUBRANGE) {
        // The path name described a subpath, and we didn't find it.
        // Don't call the iteratee.
        return true;
    }
    
    // Look at every subpath on it
    return graph.for_each_path_matching({sense}, {sample}, {locus}, [&](const path_handle_t& match) {
        // TODO: There's no way to search by haplotype and phase block, we have to scan
        if (graph.get_haplotype(match) != haplotype) {
            // Skip this haplotype
            return true;
        }
        if (graph.get_phase_block(match) != phase_block) {
            // Skip this phase block
            return true;
        }
        // Don't need to check subrange, we know we don't have one and this candidate does.
        return iteratee(match);    
    });
}

/// Returns the base path name for this path (i.e. the path's name without any subrange).
static string get_path_base_name(const PathPositionHandleGraph& graph, const path_handle_t& path) {
    if (graph.get_subrange(path) == PathMetadata::NO_SUBRANGE) {
        // This is a full path
        return graph.get_path_name(path);
    } else {
        // This is a subpath, so remember what it's a subpath of, and use that.
        return PathMetadata::create_path_name(graph.get_path_sense(path),
                                              graph.get_sample_name(path),
                                              graph.get_locus_name(path),
                                              graph.get_haplotype(path),
                                              graph.get_phase_block(path),
                                              PathMetadata::NO_SUBRANGE);
    }
}

pair<vector<pair<string, int64_t>>, unordered_map<string, int64_t>> extract_path_metadata(
    const vector<tuple<path_handle_t, size_t, size_t>>& paths,  const PathPositionHandleGraph& graph,
    bool subpath_support) {

    // Build a path name and length list from the handles (this is for the sequence dictionary and reflects *base* paths
    // as opposed to subpaths in the graph -- if there are no subpaths there is no distinction)
    vector<pair<string, int64_t>> path_names_and_lengths;
    unordered_set<string> base_path_set;    
    // Remember the actual path lengths (this is for coordinate transformations)        
    unordered_map<string, int64_t> subpath_to_length;
    for (const auto& path_info : paths) {
        auto& path = get<0>(path_info);
        auto& own_length = get<1>(path_info);
        auto& base_length = get<2>(path_info);
        string base_path_name = subpath_support ? get_path_base_name(graph, path) : graph.get_path_name(path);
        if (!base_path_set.count(base_path_name)) {
            path_names_and_lengths.push_back(make_pair(base_path_name, base_length));
            base_path_set.insert(base_path_name);
        }
        subpath_to_length[graph.get_path_name(path)] = own_length;
    }

    return make_pair(path_names_and_lengths, subpath_to_length);
}

vector<tuple<path_handle_t, size_t, size_t>> get_sequence_dictionary(const string& filename, const PathPositionHandleGraph& graph) {

    // TODO: We assume we're using the one true default path metadata <-> name mapping
    
    // Subpath support: map to paths from their base name without a subrange.
    // Includes full-length paths if they exost, and subrange-bearing paths otherwise.
    unordered_map<string, vector<path_handle_t>> base_path_to_subpaths;
    
    // Parse the input into this list.  If length was unspecified (ie in regular text file with one column) then it will be -1
    // and filled in later
    vector<pair<string, int64_t>> input_names_lengths;
    
    if (!filename.empty()) {
        // TODO: As of right now HTSLib doesn't really let you iterate the sequence dictionary when you use its parser. So we use our own parser.
        get_input_file(filename, [&](istream& in) {
            for (string line; getline(in, line);) {
                // Each line will produce a sequence name and a handle
                string sequence_name = "";
                int64_t length = -1;
                bool missing_length = false;
            
                // Trim leading and trailing whitespace
                line.erase(line.begin(), find_if(line.begin(), line.end(), [](char ch) {return !isspace(ch);}));
                line.erase(find_if(line.rbegin(), line.rend(), [](char ch) {return !isspace(ch);}).base(), line.end());
            
                if (line.empty()) {
                    // Unless it is empty
                    continue;
                }
            
                // See if each line starts with @SQ and we have to parse it, or @HD and we have to drop it, or if we have to handle it as a name.
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
                                        
                } else if (starts_with(line, "@HD")) {
                    // SAM header line also found in dict files. Drop it.
                    // TODO: Hope nobody named a sequence "@HD"-something
                    continue;
                } else {
                    // Get the name from the line and the sequence from the graph
                    vector<string> toks = split_delims(line, "\t");
                    sequence_name = toks[0];
                    if (toks.size() > 1) {
                        length = std::stol(toks[1]);
                    } else {
                        missing_length = true;
                    }
                }

                if (sequence_name == "") {
                    cerr << "error:[vg::get_sequence_dictionary] No sequence name for line " << line << endl;
                    exit(1);
                }

                if (!missing_length && length < 0) {
                    cerr << "error:[vg::get_sequence_dictionary] Unacceptable sequence length " << length << " for sequence " << sequence_name << endl;
                    exit(1);
                }
                
                if (graph.has_path(sequence_name)) {
                    // If the graph does have a path by this exact name, do a length check.
                    path_handle_t path = graph.get_path_handle(sequence_name);
                    size_t graph_path_length = graph.get_path_length(path);
                    if (missing_length) {
                        // We need to infer the length
                        length = graph_path_length;
                    } else if (graph_path_length != length) {
                        // Length was given but doesn't match
                        cerr << "error:[vg mpmap] Graph contains a path " << sequence_name << " of length " << graph_path_length
                             << " but sequence dictionary in " << filename << " indicates a length of " << length << endl;
                        exit(1);
                    }
                    // Remember the path
                    base_path_to_subpaths[sequence_name].push_back(path);
                } else {
                    // The graph doesn't have this exact path; does it have any subregions of a full path with this name?
                    // If so, remember and use those.
                    for_each_subpath_of(graph, sequence_name, [&](const path_handle_t) {
                        // We know this can't be an exact match, since we already checked for one. It must be a subrange.
                        // We found a subpath we're looking for.
                        base_path_to_subpaths[sequence_name].push_back(match);
                        // Keep looking for more.
                        return true;
                    });
                    if (!base_path_to_subpaths.count(sequence_name)) {
                        // We didn't find any subpaths for this path as a base path either.
                        cerr << "error:[vg mpmap] Graph does not have a path named " << line << ", which was indicated in " << filename << endl;
                        exit(1);
                    }
                    // The length may still be missing.
                } 

                input_names_lengths.push_back(make_pair(sequence_name, length));

            }
        });
        
        if (input_names_lengths.empty()) {
            // There were no entries in the file
            cerr << "error:[vg::get_sequence_dictionary] No sequence dictionary available in file: " << filename << endl;
            exit(1);
        }
    } else {
        // We got no paths in so we need to guess them.
        // We will deduplicate their base names, without subpath info.
        unordered_set<string> base_names;
        
        // When we find a path or subpath we want, we will keep it.
        auto keep_path_or_subpath = [&](const path_handle_t& path) {
            string base_name = get_path_base_name(graph, path);
            int64_t base_length = -1;
            if (graph.get_subrange(path) == PathMetadata::NO_SUBRANGE) {
                // This is a full path so we can determine base length now.
                base_length = graph.get_path_length(path);
            }
            if (!base_names.count(base_name)) {
                // This is the first time we have seen something on this path.
                // Remember we need a length for it.
                // TODO: make this a map so we can just max in here instead of doing another pass.
                input_names_lengths.push_back(make_pair(base_name, base_length));
                // And remember we are using it.
                base_names.insert(base_name);
            }
            // Remember this path as belonging to the right base name.
            base_path_to_subpaths[base_name].push_back(path);
        };
        
        // First look for reference sense paths and their subpaths
        graph.for_each_path_of_sense(PathSense::REFERENCE, keep_path_or_subpath);
        
        if (input_names_lengths.empty()) {
            // If none of those exist, try generic sense paths and their subpaths
            cerr << "warning:[vg::get_sequence_dictionary] No reference-sense paths available in the graph; falling back to generic paths." << endl;
            graph.for_each_path_of_sense(PathSense::GENERIC, [&](const path_handle_t& path) {
                if (Paths::is_alt(graph.get_path_name(path))) {
                    // Skip this path because it is a stored allele.
                    return;
                }
                // Otherwise, keep it.
                keep_path_or_subpath(path);
            });
        }
        
        if (input_names_lengths.empty()) {
            // No non-alt generic paths either
            cerr << "error:[vg::get_sequence_dictionary] No reference or non-alt-allele generic paths available in the graph!" << endl;
            exit(1);
        }
    }
    
    // We fill in the "dictionary" (which is what SAM calls it; it's not a mapping for us)
    // we also store the path length (from the graph) along with the base path length (from the user if specified, from paths otherwise)
    vector<tuple<path_handle_t, size_t, size_t>> dictionary;

    for (auto& base_name_and_length : input_names_lengths) {
        // For every base path name we have stuff on
        
        if (base_name_and_length.second == -1) {
            // We need the overall length of this base path still.
            for_each_subpath_of(graph, base_name_and_length.first, [&](const path_handle_t& path) {
                // So scan it and all its subpaths
                
                subrange_t subrange = graph.get_subrange(path);
                if (subrange == PathMetadata::NO_SUBRANGE) {
                    // Full path, use its length.
                    // TODO: probably we should have seen the full path's length by now if it existed.
                    base_name_and_length.second = graph.get_path_length(path);
                } else {
                    // Subpath, work out where it ends and max it in
                    auto end_offset = subrange.second;
                    if (end_offset == PathMetadata::NO_END_POSITION) {
                        // If there was a stored end we would just trust it, but without it we have to compute it.
                        end_offset = subrange.first + graph.get_path_length(path);
                    }
                    // max the end offset in against all the end offsets we have seen so far.
                    base_name_and_length.second = std::max(base_name_and_length.second, (int64_t) end_offset);
                }
            });
        }
        
        // Now we have the length, so we can fill in the dictionary.
        for (auto& path : base_path_to_subpaths[base_name_and_length.first]) {
            // For every subpath we found on the base path (possibly just
            // itself), remember the subpath, the subpath's length, and the
            // base path's length
            dictionary.push_back(make_tuple(path, graph.get_path_length(path), (size_t)base_name_and_length.second));
        }
    }

    return dictionary;
}

// Give the footer length for rewriting BGZF EOF markers.
const size_t HTSWriter::BGZF_FOOTER_LENGTH = 28;

HTSWriter::HTSWriter(const string& filename, const string& format,
    const vector<pair<string, int64_t>>& path_order_and_length,
    const unordered_map<string, int64_t>& subpath_to_length,
    size_t max_threads) :
    out_file(filename == "-" ? nullptr : new ofstream(filename)),
    multiplexer(out_file.get() != nullptr ? *out_file : cout, max_threads),
    format(format), path_order_and_length(path_order_and_length), subpath_to_length(subpath_to_length),
    backing_files(max_threads, nullptr), sam_files(max_threads, nullptr),
    atomic_header(nullptr), sam_header(), header_mutex(), output_is_bgzf(format != "SAM"),
    hts_mode() {
    
    // We can't work with no streams to multiplex, because we need to be able
    // to write BGZF EOF blocks throught he multiplexer at destruction.
    assert(max_threads > 0);
    
    if (out_file.get() != nullptr && !*out_file) {
        // Make sure we opened a file if we aren't writing to standard output
        cerr << "[vg::HTSWriter] failed to open " << filename << " for writing" << endl;
        exit(1);
    }
    
    // Make sure we have an HTS format
    assert(format == "SAM" || format == "BAM" || format == "CRAM");
    
    // Compute the file mode to send to HTSlib depending on output format
    char out_mode[5];
    string out_format = "";
    strcpy(out_mode, "w");
    if (format == "BAM") {
        out_format = "b";
    } else if (format == "CRAM") {
        out_format = "c";
    } else {
        // Must be SAM
        out_format = "";
    }
    strcat(out_mode, out_format.c_str());
    int compress_level = 9; // TODO: support other compression levels
    if (compress_level >= 0) {
        char tmp[2];
        tmp[0] = compress_level + '0'; tmp[1] = '\0';
        strcat(out_mode, tmp);
    }
    // Save to a C++ string that we will use later.
    hts_mode = out_mode;

    if (this->subpath_to_length.empty()) {
        // no subpath support: just use lengths from path_order_and_length
        for (const auto& pl : path_order_and_length) {
            this->subpath_to_length[pl.first] = pl.second;
        }
    }
    
    // Each thread will lazily open its samFile*, once it has a header ready
}

HTSWriter::~HTSWriter() {
    // Note that the destructor runs in only one thread, and only when
    // destruction is safe. No need to lock the header.
    if (atomic_header.load() != nullptr) {
        // Delete the header
        bam_hdr_destroy(atomic_header.load());
    }
    
    for (size_t thread_number = 0; thread_number < sam_files.size(); thread_number++) {
        // For each thread, find its samFile*
        auto& sam_file = sam_files.at(thread_number);
    
        if (sam_file != nullptr) {
            // Close out all the open samFile*s and flush their data before the
            // multiplexer destructs
            sam_close(sam_file);
            
            if (output_is_bgzf) {
                // Discard all the BGZF EOF marker blocks
                multiplexer.discard_bytes(thread_number, BGZF_FOOTER_LENGTH);
                
                // Put a barrier so subsequent writes come later.
                multiplexer.register_barrier(thread_number);
            }
        }
    }
    
    if (output_is_bgzf) {
        // Now put one BGZF EOF marker in thread 0's stream.
        // It will be the last thing, after all the barriers, and close the file.
        vg::io::finish(multiplexer.get_thread_stream(0), true);
    }
    
}

bam_hdr_t* HTSWriter::ensure_header(const string& read_group,
                                    const string& sample_name,
                                    size_t thread_number) {
    bam_hdr_t* header = atomic_header.load();
    if (header == nullptr) {
        // The header does not exist.
        
        // Lock the header mutex so we have exclusive control of the header
        lock_guard<mutex> header_lock(header_mutex);
        
        // Load into the enclosing scope header. Don't shadow, because the
        // enclosing scope variable is what will get returned.
        header = atomic_header.load();
        if (header == nullptr) {
            // It is our turn to make the header.
        
            // Sniff out the read group and sample, and map from RG to sample
            map<string, string> rg_sample;
            if (!read_group.empty() && !sample_name.empty()) {
                // We have a sample and a read group
                rg_sample[read_group] = sample_name;
            }
            
            // Make the header
            header = hts_string_header(sam_header, path_order_and_length, rg_sample);
            
            // Initialize the SAM file for this thread and actually keep the header
            // we write, since we are the first thread.
            initialize_sam_file(header, thread_number, true);
            
            // Save back to the atomic only after the header has been written and
            // it is safe for other threads to use it.
            atomic_header.store(header);
            
            // We made the header and the SAM file.
            return header;
        }
    }
    
    // Otherwise, someone else beat us to creating the header.
    // Header is ready. We just need to create the samFile* for this thread with it if it doesn't exist.
    
    if (sam_files[thread_number] == nullptr) {
        // The header has been created and written, but hasn't been used to initialize our samFile* yet.
        initialize_sam_file(header, thread_number);
    }
    
    return header;
}


void HTSWriter::save_records(bam_hdr_t* header, vector<bam1_t*>& records, size_t thread_number) {
    // We need a header and an extant samFile*
    assert(header != nullptr);
    assert(sam_files[thread_number] != nullptr);
    
    for (auto& b : records) {
        // Emit each record
        
        if (sam_write1(sam_files[thread_number], header, b) == 0) {
            cerr << "[vg::HTSWriter] error: writing to output file failed" << endl;
            exit(1);
        }
    }
    
    for (auto& b : records) {
        // Deallocate all the records
        bam_destroy1(b);
    }
    
    if (multiplexer.want_breakpoint(thread_number)) {
        // We have written enough that we ought to give the multiplexer a chance to multiplex soon.
        // There's no way to do this without closing and re-opening the HTS file.
        // So just tear down and reamke the samFile* for this thread.
        initialize_sam_file(header, thread_number);
    }
}

void HTSWriter::initialize_sam_file(bam_hdr_t* header, size_t thread_number, bool keep_header) {
    if (sam_files[thread_number] != nullptr) {
        // A samFile* has been created already. Clear it out.
        // Closing the samFile* flushes and destroys the BGZF (if any) and
        // hFILE* backing it.
        sam_close(sam_files[thread_number]);
        
        if (output_is_bgzf) {
            // Now we know there's a closing empty BGZF block that htslib puts to
            // mark EOF. We don't want that in the middle of our stream because it
            // is weird and we aren't actually at EOF.
            // We know how long it is, so we will trim it off.
            multiplexer.discard_bytes(thread_number, BGZF_FOOTER_LENGTH);
        }
        
        // Now place a breakpoint right where we were before that empty block.
        multiplexer.register_breakpoint(thread_number);
    }
    
    // Create a new samFile* for this thread
    // hts_mode was filled in when the header was.
    // hts_hopen demands a filename, but appears to just store it, and
    // doesn't document how required it it.
    backing_files[thread_number] = vg::io::hfile_wrap(multiplexer.get_thread_stream(thread_number));
    sam_files[thread_number] = hts_hopen(backing_files[thread_number], "-", hts_mode.c_str());
    
    if (sam_files[thread_number] == nullptr) {
        // We couldn't open the output samFile*
        cerr << "[vg::HTSWriter] failed to open internal stream for writing " << format << " output" << endl;
        exit(1);
    }
    
    // Write the header again, which is the only way to re-initialize htslib's internals.
    // Remember that sam_hdr_write flushes the BGZF to the hFILE*, but does not flush the hFILE*.
    if (sam_hdr_write(sam_files[thread_number], header) != 0) {
        cerr << "[vg::HTSWriter] error: failed to write the SAM header" << endl;
        exit(1);
    }
    
    // Now flush it out of the hFILE* buffer into the backing C++ stream
    if (hflush(backing_files[thread_number]) != 0) {
        cerr << "[vg::HTSWriter] error: failed to flush the SAM header" << endl;
        exit(1);
    }
    
    if (keep_header) {
        // We are the first thread to write a header, so we actually want it.
        // Place a barrier which is also a breakpoint, so all subsequent writes come later.
        multiplexer.register_barrier(thread_number);
    } else {
        // Discard the header so it won't be in the resulting file again
        multiplexer.discard_to_breakpoint(thread_number);
    }
}

HTSAlignmentEmitter::HTSAlignmentEmitter(const string& filename, const string& format,
                                         const vector<pair<string, int64_t>>& path_order_and_length,
                                         const unordered_map<string, int64_t>& subpath_to_length,
                                         size_t max_threads)
    : HTSWriter(filename, format, path_order_and_length, subpath_to_length, max_threads)
{
    // nothing else to do
}

void HTSAlignmentEmitter::convert_alignment(const Alignment& aln, vector<pair<int, char>>& cigar, bool& pos_rev, int64_t& pos, string& path_name) const {
    
    // We assume the position is available in refpos(0)
    assert(aln.refpos_size() == 1);
    path_name = aln.refpos(0).name();
    size_t path_len = 0;
    if (path_name != "") {
        path_len = subpath_to_length.at(path_name);
    }
    // Extract the position so that it could be adjusted by cigar_against_path if we decided to supperss softclips. Which we don't.
    // TODO: Separate out softclip suppression.
    pos = aln.refpos(0).offset();
    pos_rev = aln.refpos(0).is_reverse();
    cigar = cigar_against_path(aln, pos_rev, pos, path_len, 0);

    // Resolve subpath naming / offset
    auto subpath_info = Paths::parse_subpath_name(path_name);
    if (get<0>(subpath_info)) {
        path_name = get<1>(subpath_info);
        pos += get<2>(subpath_info);
    }
}

void HTSAlignmentEmitter::convert_unpaired(Alignment& aln, bam_hdr_t* header, vector<bam1_t*>& dest) {
    // Look up the stuff we need from the Alignment to express it in BAM.
    vector<pair<int, char>> cigar;
    bool pos_rev;
    int64_t pos;
    string path_name;
    convert_alignment(aln, cigar, pos_rev, pos, path_name);
    
    dest.emplace_back(alignment_to_bam(header,
                                       aln,
                                       path_name,
                                       pos,
                                       pos_rev,
                                       cigar));
}

void HTSAlignmentEmitter::convert_paired(Alignment& aln1, Alignment& aln2, bam_hdr_t* header, int64_t tlen_limit,
                                         vector<bam1_t*>& dest) {
    // Look up the stuff we need from the Alignment to express it in BAM.
    
    
    vector<pair<int, char>> cigar1, cigar2;
    bool pos_rev1, pos_rev2;
    int64_t pos1, pos2;
    string path_name1, path_name2;
    convert_alignment(aln1, cigar1, pos_rev1, pos1, path_name1);
    convert_alignment(aln2, cigar2, pos_rev2, pos2, path_name2);
    
    // Determine the TLEN for each read.
    auto tlens = compute_template_lengths(pos1, cigar1, pos2, cigar2);
        
    dest.emplace_back(alignment_to_bam(header,
                                       aln1,
                                       path_name1,
                                       pos1,
                                       pos_rev1,
                                       cigar1,
                                       path_name2,
                                       pos2,
                                       pos_rev2,
                                       tlens.first,
                                       tlen_limit));
    dest.emplace_back(alignment_to_bam(header,
                                       aln2,
                                       path_name2,
                                       pos2,
                                       pos_rev2,
                                       cigar2,
                                       path_name1,
                                       pos1,
                                       pos_rev1,
                                       tlens.second,
                                       tlen_limit));
    
}

void HTSAlignmentEmitter::emit_singles(vector<Alignment>&& aln_batch) {
    if (aln_batch.empty()) {
        // Nothing to do
        return;
    }
    
    // Work out what thread we are
    size_t thread_number = omp_get_thread_num();
    
    // Make sure header exists
    bam_hdr_t* header = ensure_header(aln_batch.front().read_group(),
                                      aln_batch.front().sample_name(), thread_number);
    assert(header != nullptr);
    assert(sam_files[thread_number] != nullptr);
    
    vector<bam1_t*> records;
    records.reserve(aln_batch.size());
    
    for (auto& aln : aln_batch) {
        // Convert each alignment to HTS format
        convert_unpaired(aln, header, records);
    }
    
    // Save to the stream for this thread.
    save_records(header, records, thread_number);
}

void HTSAlignmentEmitter::emit_mapped_singles(vector<vector<Alignment>>&& alns_batch) {
    // Count the total alignments to do
    size_t count = 0;
    // And find an alignment to base the header on
    Alignment* sniff = nullptr;
    for (auto& alns : alns_batch) {
        count += alns.size();
        if (!alns.empty() && sniff == nullptr) {
            sniff = &alns.front();
        }
    }
    
    if (count == 0) {
        // Nothing to do
        return;
    }
    
    // Work out what thread we are
    size_t thread_number = omp_get_thread_num();
    
    // Make sure header exists
    assert(sniff != nullptr);
    bam_hdr_t* header = ensure_header(sniff->read_group(), sniff->sample_name(),
                                      thread_number);
    assert(header != nullptr);
    assert(sam_files[thread_number] != nullptr);
    
    vector<bam1_t*> records;
    records.reserve(count);
    
    for (auto& alns : alns_batch) {
        for (auto& aln : alns) {
            // Convert each alignment to HTS format
            convert_unpaired(aln, header, records);
        }
    }
    
    // Save to the stream for this thread.
    save_records(header, records, thread_number);

}


void HTSAlignmentEmitter::emit_pairs(vector<Alignment>&& aln1_batch,
                                     vector<Alignment>&& aln2_batch,
                                     vector<int64_t>&& tlen_limit_batch) {
    assert(aln1_batch.size() == aln2_batch.size());
    assert(aln1_batch.size() == tlen_limit_batch.size());
    
    
    if (aln1_batch.empty()) {
        // Nothing to do
        return;
    }
    
    // Work out what thread we are
    size_t thread_number = omp_get_thread_num();
    
    // Make sure header exists
    bam_hdr_t* header = ensure_header(aln1_batch.front().read_group(),
                                      aln1_batch.front().sample_name(), thread_number);
    assert(header != nullptr);
    assert(sam_files[thread_number] != nullptr);
    
    vector<bam1_t*> records;
    records.reserve(aln1_batch.size() * 2);
    
    for (size_t i = 0; i < aln1_batch.size(); i++) {
        // Convert each alignment pair to HTS format
        convert_paired(aln1_batch[i], aln2_batch[i], header, tlen_limit_batch[i], records);
    }
    
    // Save to the stream for this thread.
    save_records(header, records, thread_number);
}
    
    
void HTSAlignmentEmitter::emit_mapped_pairs(vector<vector<Alignment>>&& alns1_batch, 
                                            vector<vector<Alignment>>&& alns2_batch,
                                            vector<int64_t>&& tlen_limit_batch) {
                                            
    assert(alns1_batch.size() == alns2_batch.size());
    assert(alns1_batch.size() == tlen_limit_batch.size());
    
    // Count the total alignments to do
    size_t count = 0;
    // And find an alignment to base the header on
    Alignment* sniff = nullptr;
    for (size_t i = 0; i < alns1_batch.size(); i++) {
        // Go through all pairs
        
        // Make sure each end of the pair has the same number of mappings
        assert(alns1_batch[i].size() == alns2_batch[i].size());
        
        count += alns1_batch[i].size() * 2;
        if (!alns1_batch[i].empty() && sniff == nullptr) {
            sniff = &alns1_batch[i].front();
        }
    }
    
    if (count == 0) {
        // Nothing to do
        return;
    }
    
    // Work out what thread we are
    size_t thread_number = omp_get_thread_num();
    
    // Make sure header exists
    assert(sniff != nullptr);
    bam_hdr_t* header = ensure_header(sniff->read_group(), sniff->sample_name(),
                                      thread_number);
    assert(header != nullptr);
    assert(sam_files[thread_number] != nullptr);
    
    vector<bam1_t*> records;
    records.reserve(count);
    
    for (size_t i = 0; i < alns1_batch.size(); i++) {
        for (size_t j = 0; j < alns1_batch[i].size(); j++) {
            // Convert each alignment pair to HTS format
            convert_paired(alns1_batch[i][j], alns2_batch[i][j], header, tlen_limit_batch[i], records);
        }
    }
    
    // Save to the stream for this thread.
    save_records(header, records, thread_number);
}

SplicedHTSAlignmentEmitter::SplicedHTSAlignmentEmitter(const string& filename, const string& format,
                                                       const vector<pair<string, int64_t>>& path_order_and_length,
                                                       const unordered_map<string, int64_t>& subpath_to_length,
                                                       const PathPositionHandleGraph& graph,
                                                       size_t max_threads) :
    HTSAlignmentEmitter(filename, format, path_order_and_length, subpath_to_length, max_threads), graph(graph) {
    
    // nothing else to do
}

void SplicedHTSAlignmentEmitter::convert_alignment(const Alignment& aln, vector<pair<int, char>>& cigar,
                                                   bool& pos_rev, int64_t& pos, string& path_name) const {
    
    // We assume the position is available in refpos(0)
    assert(aln.refpos_size() == 1);
    path_name = aln.refpos(0).name();
    pos = aln.refpos(0).offset();
    pos_rev = aln.refpos(0).is_reverse();
    
    // Convert to a cigar with spliced deletions
    cigar = spliced_cigar_against_path(aln, path_name, pos, pos_rev);
}

vector<pair<int, char>> SplicedHTSAlignmentEmitter::spliced_cigar_against_path(const Alignment& aln,
                                                                               const string& path_name,
                                                                               int64_t pos, bool rev) const {
    // the return value
    vector<pair<int, char>> cigar;
    
    if (aln.has_path() && aln.path().mapping_size() > 0) {
        // the read is aligned to the path
        
        path_handle_t path_handle = graph.get_path_handle(path_name);
        step_handle_t step = graph.get_step_at_position(path_handle, pos);
        
        // to indicate whether we've found the edit that corresponds to the BAM position
        bool found_pos = false;
        
        const Path& path = aln.path();
        for (size_t i = 0; i < path.mapping_size(); ++i) {
            
            // we traverse backwards on a reverse strand mapping
            const Mapping& mapping = path.mapping(rev ? path.mapping_size() - 1 - i : i);
            
            for (size_t j = 0; j < mapping.edit_size(); ++j) {
                                
                // we traverse backwards on a reverse strand mapping
                const Edit& edit = mapping.edit(rev ? mapping.edit_size() - 1 - j : j);
                
                if (!found_pos) {
                    // we may still be searching through an initial softclip to find
                    // the edit that corresponds to the BAM position
                    if (edit.to_length() > 0 && edit.from_length() == 0) {
                        append_cigar_operation(edit.to_length(), 'S', cigar);
                        // skip the main block where we assign cigar operations
                        continue;
                    }
                    else {
                        found_pos = true;
                    }
                }
                
                // identify the cigar operation
                char cigar_code;
                int length;
                if (edit.from_length() == edit.to_length()) {
                    cigar_code = 'M';
                    length = edit.from_length();
                }
                else if (edit.from_length() > 0 && edit.to_length() == 0) {
                    cigar_code = 'D';
                    length = edit.from_length();
                }
                else if (edit.to_length() > 0 && edit.from_length() == 0) {
                    cigar_code = 'I';
                    length = edit.to_length();
                }
                else {
                    throw std::runtime_error("Spliced CIGAR construction can only convert simple edits");
                }
                
                append_cigar_operation(length, cigar_code, cigar);
            } // close loop over edits
            
            if (found_pos && i + 1 < path.mapping_size()) {
                // we're anchored on the path by the annotated position, and we're transitioning between
                // two mappings, so we should check for a deletion/splice edge
                
                step_handle_t next_step = graph.get_next_step(step);
                
                handle_t next_handle = graph.get_handle_of_step(next_step);
                const Position& next_pos = path.mapping(rev ? path.mapping_size() - 2 - i : i + 1).position();
                if (graph.get_id(next_handle) != next_pos.node_id()
                    || (graph.get_is_reverse(next_handle) != next_pos.is_reverse()) != rev) {
                    
                    // the next mapping in the alignment is not the next mapping on the path, so we must have
                    // taken a deletion
                    
                    // find the closest step that is further along the path than the current one
                    // and matches the next position (usually there will only be one)
                    size_t curr_offset = graph.get_position_of_step(step);
                    size_t nearest_offset = numeric_limits<size_t>::max();
                    graph.for_each_step_on_handle(graph.get_handle(next_pos.node_id()),
                                                  [&](const step_handle_t& candidate) {
                        
                        if (graph.get_path_handle_of_step(candidate) == path_handle) {
                            size_t candidate_offset = graph.get_position_of_step(candidate);
                            if (candidate_offset < nearest_offset && candidate_offset > curr_offset) {
                                nearest_offset = candidate_offset;
                                next_step = candidate;
                            }
                        }
                    });
                    
                    if (nearest_offset == numeric_limits<size_t>::max()) {
                        throw std::runtime_error("Spliced BAM conversion could not find path steps that match alignment");
                    }
                    
                    // the gap between the current step and the next one along the path
                    size_t deletion_length = (nearest_offset - curr_offset -
                                              graph.get_length(graph.get_handle_of_step(step)));
                    
                    // add to the cigar
                    if (deletion_length >= min_splice_length) {
                        // long enough to be a splice
                        append_cigar_operation(deletion_length, 'N', cigar);
                    }
                    else if (deletion_length) {
                        // create or extend a deletion
                        append_cigar_operation(deletion_length, 'D', cigar);
                    }
                }
                
                // iterate along the path
                step = next_step;
            }
        } // close loop over mappings
        
        if (cigar.back().second == 'I') {
            // the final insertion is actually a softclip
            cigar.back().second = 'S';
        }
    }
    
    consolidate_ID_runs(cigar);
    
    return cigar;
}


}

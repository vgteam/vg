#include "gfa.hpp"
#include "utility.hpp"
#include "path.hpp"
#include "rgfa.hpp"
#include <sstream>

#include <gbwtgraph/utils.h>

namespace vg {

using namespace std;

/// Determine if a path should be written as a GFA W line or a GFA P line.
static bool should_write_as_w_line(const PathHandleGraph* graph, path_handle_t path_handle);
/// Write out a W line for a path. Uses a map to keep track of fake offset
/// ranges used to distinguish multiple phase blocks on a haplotype, since GFA
/// doesn't support them.
static void write_w_line(const PathHandleGraph* graph, ostream& out, path_handle_t path_handle, unordered_map<tuple<string, int64_t, string>, size_t>& last_phase_block_end);

void graph_to_gfa(const PathHandleGraph* graph, ostream& out, const set<string>& rgfa_paths,
                  bool rgfa_pline, bool use_w_lines) {

    RGFACover rgfa_cover;
    if (!rgfa_paths.empty()) {
        // index the rgfa cover in the graph (combination of rank-0 rgfa_paths and rank-1 paths with rgfa sample name)
        unordered_set<path_handle_t> rgfa_path_handles;
        for (const string& rgfa_path_name : rgfa_paths) {
            rgfa_path_handles.insert(graph->get_path_handle(rgfa_path_name));
        }
        rgfa_cover.load(graph, rgfa_path_handles);
    }
    // TODO: Support sorting nodes, paths, and/or edges for canonical output
    // TODO: Use a NamedNodeBackTranslation (or forward translation?) to properly round-trip GFA that has had to be chopped.
    
    // Compute reference-sense sample header tags
    unordered_set<string> reference_samples;
    graph->for_each_path_matching({PathSense::REFERENCE}, {}, {}, [&](const path_handle_t& h) {
            if (!rgfa_paths.count(graph->get_path_name(h)) || rgfa_pline) {
                // If it is going to be something other than an rGFA path,
                // we'll have to convey its reference-ness another way.
                reference_samples.insert(graph->get_sample_name(h));
            }
        });
    
    // Start with the header for a GFA1.1 file
    out << "H\tVN:Z:1.1";
    if (!reference_samples.empty()) {
        // Include a reference sample name tag if we have reference paths.
        out << "\t" << gbwtgraph::REFERENCE_SAMPLE_LIST_GFA_TAG << ":Z:" << gbwtgraph::compose_reference_samples_tag(reference_samples);
    }
    out << "\n";

    //Compute the rGFA tags of given paths.  These paths are the rank-0 reference paths (passed in rgfa_paths set)
    //along with paths that have the special rgfa sample name (rank>0) paths.
    unordered_map<nid_t, pair<path_handle_t, size_t>> node_offsets;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        string path_name = graph->get_path_name(path_handle);
        if (rgfa_paths.count(path_name) || RGFACover::is_rgfa_path_name(path_name)) {
            size_t offset = 0;
            subrange_t path_subrange = graph->get_subrange(path_handle);
            if (path_subrange != PathMetadata::NO_SUBRANGE) {
                offset = path_subrange.first;
            }
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                handle_t handle = graph->get_handle_of_step(step_handle);
                nid_t node_id = graph->get_id(handle);
                if (graph->get_is_reverse(handle)) {
                    stringstream ss;
                    ss << "error [gfa]: unable to write rGFA tags for path " << path_name << " because node "
                       << node_id << " is traversed on its reverse strand.  rGFA only supports the forward strand." << endl;
                    throw runtime_error(ss.str());
                }
                if (node_offsets.count(node_id)) {
                    cerr << "warning [gfa]: multiple selected rgfa paths found on node " << node_id << ": keeping tags for "
                         << graph->get_path_name(node_offsets[node_id].first) << " and ignoring those for " << path_name << endl;
                } else {
                    node_offsets[node_id] = make_pair(path_handle, offset);
                }
                offset += graph->get_length(handle);
            });
        }
    });
  
    //Go through each node in the graph
    graph->for_each_handle([&](const handle_t& h) {
        out << "S\t";
        nid_t node_id = graph->get_id(h);
        out << node_id << "\t";
        out << graph->get_sequence(h);
        auto it = node_offsets.find(node_id);
        if (it != node_offsets.end()) {
            // add rGFA tags
            string sample_name = graph->get_sample_name(it->second.first);
            string locus_name = graph->get_locus_name(it->second.first);
            int64_t haplotype = graph->get_haplotype(it->second.first);
            if (RGFACover::is_rgfa_path_name(graph->get_path_name(it->second.first))) {
                std::tie(sample_name, locus_name) = RGFACover::parse_rgfa_locus_name(locus_name);
            }
            string rgfa_sn = PathMetadata::create_path_name(sample_name.empty() ? PathSense::GENERIC : PathSense::REFERENCE,
                                                            sample_name, locus_name,
                                                            sample_name.empty() ? PathMetadata::NO_HAPLOTYPE : haplotype,
                                                            PathMetadata::NO_PHASE_BLOCK,
                                                            PathMetadata::NO_SUBRANGE);
            out << "\t" << "SN:Z:" << rgfa_sn
                << "\t" << "SO:i:" << it->second.second
                << "\t" << "SR:i:" << rgfa_cover.get_rank(node_id);
        }
        out << "\n"; // Writing `std::endl` would flush the buffer.
        return true;
    });
    
    // Sort the paths by name, making sure to treat subpath coordinates numerically
    vector<path_handle_t> path_handles;
    graph->for_each_path_matching(nullptr, nullptr, nullptr, [&](const path_handle_t& h) {
            path_handles.push_back(h);
        });
    std::sort(path_handles.begin(), path_handles.end(), [&](const path_handle_t& p1, const path_handle_t& p2) {
            string n1 = graph->get_path_name(p1);
            string n2 = graph->get_path_name(p2);
            subrange_t subrange1;
            subrange_t subrange2;
            n1 = Paths::strip_subrange(n1, &subrange1);
            n2 = Paths::strip_subrange(n2, &subrange2);
            if (n1 < n2) {
                return true;
            } else if (n1 == n2) {
                return subrange1 < subrange2;
            }
            return false;
        });

    vector<path_handle_t> w_line_paths;

    // Paths as P-lines
    for (const path_handle_t& h : path_handles) {
        auto path_name = graph->get_path_name(h);
        if (rgfa_pline || !rgfa_paths.count(path_name)) {
            if (graph->get_sense(h) != PathSense::REFERENCE && reference_samples.count(graph->get_sample_name(h))) {
                // We have a mix of reference and non-reference paths on the same sample which GFA can't handle.
                cerr << "warning [gfa]: path " << path_name << " will be interpreted as reference sense "
                     << "because reference paths exist on its sample" << endl;
            }
        
            if (use_w_lines && should_write_as_w_line(graph, h)) {
                w_line_paths.push_back(h);
            } else {
                out << "P\t";
                out << path_name << "\t";
                
                bool first = true;
                graph->for_each_step_in_path(h, [&](const step_handle_t& ph) {
                    handle_t step_handle = graph->get_handle_of_step(ph);
                    
                    if (!first) {
                        out << ',';
                    }
                    out << graph->get_id(step_handle);
                    out << (graph->get_is_reverse(step_handle) ? '-' : '+');
                    first = false;
                    return true;
                });
                
                out << "\t*" << "\n";
            }
        }
    }
    
    // Paths as W-lines
    {
        unordered_map<tuple<string, int64_t, string>, size_t> last_phase_block_end;
        for (const path_handle_t& h : w_line_paths) {
            write_w_line(graph, out, h, last_phase_block_end);
        }
    }

    graph->for_each_edge([&](const edge_t& h) {
        
        nid_t from_id = graph->get_id(h.first);
        bool from_is_reverse = graph->get_is_reverse(h.first);
        nid_t to_id = graph->get_id(h.second);
        bool to_is_reverse = graph->get_is_reverse(h.second);
    
        if (from_is_reverse && (to_is_reverse || to_id < from_id)) {
            // Canonicalize edges to be + orientation first if possible, and
            // then low-ID to high-ID if possible, for testability. This edge
            // needs to flip.
            
            // Swap the nodes
            std::swap(from_id, to_id);
            // Swap the orientations
            std::swap(from_is_reverse, to_is_reverse);
            // Reverse the orientations
            from_is_reverse = !from_is_reverse;
            to_is_reverse = !to_is_reverse;
        }
        
        out << "L\t" << from_id << "\t" << (from_is_reverse ? '-' : '+')
            << "\t" << to_id << "\t" << (to_is_reverse ? '-' : '+') << "\t0M\n"; // Writing `std::endl` would flush the buffer.
        return true;
    }, false);
}

bool should_write_as_w_line(const PathHandleGraph* graph, path_handle_t path_handle) {
    // Until we can change the tests, default to sending reference and
    // haplotype paths as W lines, and generic paths as P lines. 
    return graph->get_sense(path_handle) != PathSense::GENERIC;
}

void write_w_line(const PathHandleGraph* graph, ostream& out, path_handle_t path_handle, unordered_map<tuple<string, int64_t, string>, size_t>& last_phase_block_end) {
    // Extract the path metadata
    string sample = graph->get_sample_name(path_handle);
    string contig = graph->get_locus_name(path_handle);
    int64_t hap_index = graph->get_haplotype(path_handle);
    int64_t phase_block = graph->get_phase_block(path_handle);
    auto subrange = graph->get_subrange(path_handle);
    size_t start_offset = 0;
    size_t end_offset = 0;
    if (subrange != PathMetadata::NO_SUBRANGE) {
        start_offset = subrange.first;
        if (subrange.second != PathMetadata::NO_END_POSITION) {
            end_offset = subrange.second;
        }
    }
    
    if (sample == PathMetadata::NO_SAMPLE_NAME) {
        // Represent an elided sample name with "*";
        sample = "*";
    }
    
    if (hap_index == PathMetadata::NO_HAPLOTYPE) {
        // No haplotype is actually assigned here.
        // We probably won't have paths with it assigned and not assigned but
        // the same sample and contig, so assign it 0 and make the sample
        // haploid.
        // TODO: check for collisions somehow?
        hap_index = 0;
    }
     
    // Get the path length.
    // TODO: sniff if the graph has this cached somehow?
    size_t path_length = 0 ;
    graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
            path_length += graph->get_length(graph->get_handle_of_step(step_handle));
        });

    if (end_offset != 0 && start_offset + path_length != end_offset) {
        cerr << "[gfa] warning: incorrect end offset (" << end_offset << ") extracted from from path name " << graph->get_path_name(path_handle)
             << ", using " << (start_offset + path_length) << " instead" << endl;
    }
    
    // See if we need to bump along the start offset to avoid collisions of phase blocks
    auto key = std::tuple<string, int64_t, string>(sample, hap_index, contig);
    auto& phase_block_end_cursor = last_phase_block_end[key];
    if (phase_block_end_cursor != 0) {
        if (start_offset != 0) {
            // TODO: Work out a way to support phase blocks and subranges at the same time.
            cerr << "[gfa] error: cannot write multiple phase blocks on a sample, haplotyope, and contig in GFA format"
                 << " when paths already have subranges. Fix path " << graph->get_path_name(path_handle) << endl;
            exit(1);
        }
        // Budge us to after the last thing and budge the cursor to after us.
        // TODO: GBWTGraph algorithm just uses phase block number as start
        // position so it can roudn trip. Settle on a way to round trip the
        // small phase block numbers somehow?
        start_offset += phase_block_end_cursor;
        phase_block_end_cursor += path_length;
    }

    out << "W\t" << sample << "\t" << hap_index << "\t" << contig << "\t" << start_offset << "\t" << (start_offset + path_length) << "\t";

    graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
            handle_t handle = graph->get_handle_of_step(step_handle);
            out << (graph->get_is_reverse(handle) ? "<" : ">") << graph->get_id(handle);
        });
    out << "\n";
}

}

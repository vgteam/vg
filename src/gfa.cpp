#include "gfa.hpp"
#include "utility.hpp"
#include "path.hpp"
#include <sstream>
#include <algorithm>

#include <gbwtgraph/utils.h>

//#define debug

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

    //Compute the rGFA tags of given paths 
    unordered_map<nid_t, pair<path_handle_t, int64_t>> node_offsets;
    for (const string& path_name : rgfa_paths) {
        path_handle_t path_handle = graph->get_path_handle(path_name);
        size_t offset = 0;
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
  
    //Go through each node in the graph
    graph->for_each_handle([&](const handle_t& h) {
        out << "S\t";
        nid_t node_id = graph->get_id(h);
        out << node_id << "\t";
        out << graph->get_sequence(h);
        auto it = node_offsets.find(node_id);
        if (it != node_offsets.end()) {
            // hack off the rgfa stuff from the path, and move into tags
            string path_name = graph->get_path_name(get<0>(it->second));
            string rgfa_tags;
            string base_name;
            if (get_rgfa_rank(path_name) >= 0) {
                base_name = parse_rgfa_name_into_tags(graph->get_path_name(it->second.first), rgfa_tags);
            } else {
                // no tags, must be the rank-0 reference
                base_name = path_name;
                rgfa_tags = "SR:i:0";
            }                
            // hack off the subrange offset (and add it to SO)
            PathSense sense;
            string sample, locus;
            size_t haplotype, phase_block;
            subrange_t subrange;
            PathMetadata::parse_path_name(base_name, sense, sample, locus, haplotype, phase_block, subrange);
            int64_t base_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
            base_name = PathMetadata::create_path_name(sense, sample, locus, haplotype, phase_block,
                                                       PathMetadata::NO_SUBRANGE);
            // add rGFA tags            
            out << "\t" << "SN:Z:" << base_name
                << "\t" << "SO:i:" << (base_offset + it->second.second);
            if (subrange != PathMetadata::NO_SUBRANGE) {
                // todo : assert this doesn't happen (at least off-reference)
                assert(subrange.second > subrange.first);
                out << "\t" << "SL:i" << to_string(subrange.second - subrange.first);
            }
            out << "\t" << rgfa_tags;
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
        if (get_rgfa_rank(path_name) > 0) {
            if (!rgfa_paths.empty()) {
                // the path was put into tags, no reason to deal with it here
                continue;
            }
        }
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
            string path_name = graph->get_path_name(h);
            if (get_rgfa_rank(path_name) > 0) {
                if (!rgfa_paths.empty()) {
                    // the path was put into tags, no reason to deal with it here
                    continue;
                }
            }
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
            cerr << "[gfa] error: cannot write multiple phase blocks on a sample, haplotype, and contig in GFA format"
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

int get_rgfa_rank(const string& path_name) {
    int rank = -1;

    PathSense path_sense;
    string path_sample;
    string path_locus;
    size_t path_haplotype;
    size_t path_phase_block;
    subrange_t path_subrange;
    PathMetadata::parse_path_name(path_name, path_sense, path_sample, path_locus,
                                  path_haplotype, path_phase_block, path_subrange);
    
    size_t pos = path_locus.rfind(":SR:i:");
    if (pos != string::npos && path_locus.length() - pos >= 6) {
        pos += 6;
        size_t pos2 = path_locus.find(":", pos);
        size_t len = pos2 == string::npos ? pos2 : pos2 - pos;
        string rank_string = path_locus.substr(pos, len);
        rank = parse<int>(rank_string);
    }
    return rank;
}

string create_rgfa_path_name(const string& path_name, int rgfa_rank, const subrange_t& rgfa_subrange,
                             const string& rgfa_sample,
                             const string& start_parent_path_name,
                             int64_t start_parent_offset,
                             int64_t start_parent_node_id,
                             bool start_parent_node_reversed,
                             const string& end_parent_path_name,
                             int64_t end_parent_offset,
                             int64_t end_parent_node_id,
                             bool end_parent_node_reversed) {
    PathSense path_sense;
    string path_sample;
    string path_locus;
    size_t path_haplotype;
    size_t path_phase_block;
    subrange_t path_subrange;
    PathMetadata::parse_path_name(path_name, path_sense, path_sample, path_locus,
                                  path_haplotype, path_phase_block, path_subrange);

    // we're going to store absolutely everything in the locus name, including the sample name (behind SN:Z:)
    assert(path_locus != PathMetadata::NO_LOCUS_NAME);
    string base_name;
    if (path_sample != PathMetadata::NO_SAMPLE_NAME) {
        base_name = "SN:Z:" + path_sample;
    }
    // the contig name will be behind SC...
    base_name += ":SC:Z:" + path_locus;

    // and we also load in our rGFA rank
    base_name += ":SR:i:" + std::to_string(rgfa_rank);

    // start parent parent path (swap # for % as hack to not much up parent paths formatting)
    // turns out we can't embed []'s either, so swap those too
    if (!start_parent_path_name.empty()) {
        string spn = start_parent_path_name;
        std::replace(spn.begin(), spn.end(), '#', '%');
        std::replace(spn.begin(), spn.end(), '[', '(');
        std::replace(spn.begin(), spn.end(), ']', ')');                   
        base_name += ":SPP:Z:" + spn;
        if (start_parent_offset >= 0) {
            base_name += ":SPO:i:" + to_string(start_parent_offset);
        }
        base_name += ":SPN:Z:" + string(start_parent_node_reversed ? "<" : ">") + to_string(start_parent_node_id);
    }
    // end parent parent path (swap # for % as hack to not much up parent paths formatting)
    // turns out we can't embed []'s either, so swap those too    
    if (!end_parent_path_name.empty()) {
        string spn = end_parent_path_name;
        std::replace(spn.begin(), spn.end(), '#', '%');
        std::replace(spn.begin(), spn.end(), '[', '(');
        std::replace(spn.begin(), spn.end(), ']', ')');                           
        base_name += ":EPP:Z:" + spn;
        if (end_parent_offset >= 0) {
            base_name += ":EPO:i:" + to_string(end_parent_offset);
        }
        base_name += ":EPN:Z:" + string(end_parent_node_reversed ? "<" : ">") + to_string(end_parent_node_id);        
    }    

    // and return the final path, with sample/locus/rgfa-rank embedded in locus
    // (as it's a reference path, we alsos strip out the phase block)
    return PathMetadata::create_path_name(PathSense::REFERENCE, rgfa_sample, base_name, path_haplotype,
                                          PathMetadata::NO_PHASE_BLOCK, rgfa_subrange);
}

string parse_rgfa_path_name(const string& path_name, int* rgfa_rank,
                            string* rgfa_sample,
                            string* start_parent_path_name,
                            int64_t* start_parent_offset,
                            int64_t* start_parent_node_id,
                            bool* start_parent_node_reversed,
                            string* end_parent_path_name,
                            int64_t* end_parent_offset,
                            int64_t* end_parent_node_id,
                            bool* end_parent_node_reversed) {


    // begin by parsing the # structure of the path
    // note that all the rgfa stuff we're pulling out here will be embedded in the "locus"
    PathSense path_sense;
    string path_sample;
    string path_locus;
    size_t path_haplotype;
    size_t path_phase_block;
    subrange_t path_subrange;
    PathMetadata::parse_path_name(path_name, path_sense, path_sample, path_locus,
                                  path_haplotype, path_phase_block, path_subrange);

    if (rgfa_sample) {
        *rgfa_sample = path_sample;
    }

    string sample_name;
    string contig_name;
    
    // now we parse the locus which is just a list of tags sepearted by :
    // todo (again), this is too wordy and should be compacted
    vector<string> toks = split_delims(path_locus, ":");
    assert(toks.size() % 3 == 0);
    for (int64_t i = 0; i < toks.size(); i+=3) {
        if (toks[i] == "SN") {
            assert(toks[i+1] == "Z");
            sample_name = toks[i+2];
        } else if (toks[i] == "SC") {
            assert(toks[i+1] == "Z");
            contig_name = toks[i+2];
        } else if (toks[i] == "SR") {
            assert(toks[i+1] == "i");
            if (rgfa_rank) {
                *rgfa_rank = parse<int64_t>(toks[i+2]);
            }
        } else if (toks[i] == "SPP") {
            assert(toks[i+1] == "Z");
            if (start_parent_path_name) {
                *start_parent_path_name = toks[i+2];
                std::replace(start_parent_path_name->begin(), start_parent_path_name->end(), '%', '#');
                std::replace(start_parent_path_name->begin(), start_parent_path_name->end(), '(', '[');
                std::replace(start_parent_path_name->begin(), start_parent_path_name->end(), ')', ']');                
            }
        } else if (toks[i] == "SPO") {
            assert(toks[i+1] == "i");
            if (start_parent_offset) {
                *start_parent_offset = parse<int64_t>(toks[i+2]);
            }
        } else if (toks[i] == "SPN") {
            assert(toks[i+1] == "Z");
            if (start_parent_node_id) {
                *start_parent_node_id = parse<int64_t>(toks[i+2].substr(1));
            }
            if (start_parent_node_reversed) {
                *start_parent_node_reversed = toks[i+2][0] == '<';
            }
        } else if (toks[i] == "EPP") {
            assert(toks[i+1] == "Z");
            if (end_parent_path_name) {
                *end_parent_path_name = toks[i+2];
                std::replace(end_parent_path_name->begin(), end_parent_path_name->end(), '%', '#');
                std::replace(end_parent_path_name->begin(), end_parent_path_name->end(), '(', '[');
                std::replace(end_parent_path_name->begin(), end_parent_path_name->end(), ')', ']');                
            }
        } else if (toks[i] == "EPO") {
            assert(toks[i+1] == "i");
            if (end_parent_offset) {
                *end_parent_offset = parse<int64_t>(toks[i+2]);
            }
        } else if (toks[i] == "EPN") {
            assert(toks[i+1] == "Z");
            if (end_parent_node_id) {
                *end_parent_node_id = parse<int64_t>(toks[i+2].substr(1));
            }
            if (end_parent_node_reversed) {
                *end_parent_node_reversed = toks[i+2][0] == '<';
            }            
        } else {
            assert(false);
        }
    }

    if (path_sense == PathSense::REFERENCE && sample_name.empty()) {
        // the embedded path never actually had a sample, so we need to flip back to generic
        path_sense = PathSense::GENERIC;
        path_haplotype = PathMetadata::NO_HAPLOTYPE;
        path_phase_block = PathMetadata::NO_PHASE_BLOCK;
    }
        
    // reconstruct the vg path with rgfa stuf stripped from locus
    return PathMetadata::create_path_name(path_sense,
                                          sample_name,
                                          contig_name,
                                          path_haplotype,
                                          path_phase_block,
                                          path_subrange);

}

string parse_rgfa_name_into_tags(const string& path_name, string& rgfa_tags) {
    int rgfa_rank;
    string rgfa_sample; 
    string start_parent_rgfa_path_name;
    int64_t start_parent_offset;
    int64_t start_parent_node_id;
    bool start_parent_node_reversed;
    string end_parent_rgfa_path_name;
    int64_t end_parent_offset;
    int64_t end_parent_node_id;
    bool end_parent_node_reversed;
    string vg_path_name = parse_rgfa_path_name(path_name, &rgfa_rank,
                                               &rgfa_sample,
                                               &start_parent_rgfa_path_name,
                                               &start_parent_offset,
                                               &start_parent_node_id,
                                               &start_parent_node_reversed,
                                               &end_parent_rgfa_path_name,
                                               &end_parent_offset,
                                               &end_parent_node_id,
                                               &end_parent_node_reversed);    
    rgfa_tags += "SR:i:" + to_string(rgfa_rank);
    if (!start_parent_rgfa_path_name.empty()) {
        rgfa_tags += "\tSPP:Z:" + start_parent_rgfa_path_name +
            "\tSPO:i:" + to_string(start_parent_offset) +
            "\tSPN:Z:" + (start_parent_node_reversed ? "<" : ">") + to_string(start_parent_node_id);
    }
    if (!end_parent_rgfa_path_name.empty()) {
        rgfa_tags += "\tEPP:Z:" + end_parent_rgfa_path_name +
            "\tEPO:i:" + to_string(end_parent_offset) +
            "\tEPN:Z:" + (end_parent_node_reversed ? "<" : ">") + to_string(end_parent_node_id);
    }    
    
    return vg_path_name;
}



string strip_rgfa_path_name(const string& path_name) {

    PathSense path_sense;
    string path_sample;
    string path_locus;
    size_t path_haplotype;
    size_t path_phase_block;
    subrange_t path_subrange;
    PathMetadata::parse_path_name(path_name, path_sense, path_sample, path_locus,
                                  path_haplotype, path_phase_block, path_subrange);

    assert(path_locus != PathMetadata::NO_LOCUS_NAME);

    size_t sr_pos = path_locus.rfind(":SR:i:");
    if (sr_pos != string::npos && path_locus.length() - sr_pos >= 6) {
        size_t sn_pos = path_locus.rfind("SN:Z:", sr_pos - 1);
        assert(sn_pos != string::npos);

        string orig_sample;
        if (sn_pos > 0) {
            orig_sample = path_locus.substr(0, sn_pos - 1);
        }
        string orig_locus = path_locus.substr(sn_pos + 5, sr_pos - sn_pos - 5);

        // todo: recover path sense  / haploblock?
        if (orig_sample.empty()) {
            path_sense = PathSense::GENERIC;
            path_haplotype = PathMetadata::NO_HAPLOTYPE;
        }
        return PathMetadata::create_path_name(path_sense, orig_sample, orig_locus,
                                              path_haplotype, path_phase_block, path_subrange);
    }
    return path_name;
}

void clamp_path_subrange(string& path_name, int64_t start, int64_t end) {    
    PathSense path_sense;
    string path_sample;
    string path_locus;
    size_t path_haplotype;
    size_t path_phase_block;
    subrange_t path_subrange;
    PathMetadata::parse_path_name(path_name, path_sense, path_sample, path_locus,
                                  path_haplotype, path_phase_block, path_subrange);

    if (path_subrange == PathMetadata::NO_SUBRANGE) {
        path_subrange.first = start;
        path_subrange.second = end;
    } else {
        assert(end <= path_subrange.second - path_subrange.first && end > start);
        path_subrange.first += start;
        path_subrange.second = path_subrange.first + (end - start);
    }
    path_name = PathMetadata::create_path_name(path_sense, path_sample, path_locus,
                                               path_haplotype, path_phase_block, path_subrange);
}

void rgfa_graph_cover(MutablePathMutableHandleGraph* graph,
                      SnarlManager* snarl_manager,
                      const unordered_set<path_handle_t>& reference_paths,
                      int64_t minimum_length,
                      const string& rgfa_sample_name,
                      const unordered_map<string, vector<pair<int64_t, int64_t>>>& preferred_intervals){

    // for sanity's sake, we don't want to ever support multiple rgfa covers, so start by
    // deleting all existing rgfa fragments (except for rank 0 reference paths, of course)
    vector<path_handle_t> existing_cover;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        if (get_rgfa_rank(graph->get_path_name(path_handle)) > 0) {
            existing_cover.push_back(path_handle);
            assert(!reference_paths.count(path_handle));
        }
    });
    for (path_handle_t path_handle : existing_cover) {
        graph->destroy_path(path_handle);
    }    
    
    // we use the path traversal finder for everything
    // (even gbwt haplotypes, as we're using the path handle interface)
    PathTraversalFinder path_trav_finder(*graph, *snarl_manager);
        
    // we collect the rgfa cover in parallel as a list of path fragments
    size_t thread_count = get_thread_count();
    vector<vector<RGFAFragment>> thread_covers(thread_count);
    
    // we process top-level snarls in parallel
    snarl_manager->for_each_top_level_snarl_parallel([&](const Snarl* snarl) {
        // per-thread output
        // each fragment is a rank and vector of steps, the cover is a list of fragments
        // TODO: we can store just a first step and count instead of every fragment
        // The last two numbers are the indexes of the start and end parent fragments (ie values in cover_not_to_fragment)
        vector<RGFAFragment>& cover_fragments = thread_covers.at(omp_get_thread_num());
        // we also index the fragments by their node ids for fast lookups of what's covered by what
        // the value here is an index in the above vector
        unordered_map<nid_t, int64_t> cover_node_to_fragment;
        
        vector<const Snarl*> queue = {snarl}; 

        while(!queue.empty()) {
            const Snarl* cur_snarl = queue.back();
            queue.pop_back();

            // get the snarl cover, writing to cover_nodes and cover_fragments
            // note that we are single-threaded per top-level snarl, at least for now
            // this is because parent snarls and child snarls can potentially cover the
            // sname things
            rgfa_snarl_cover(graph,
                             *cur_snarl,
                             path_trav_finder,
                             reference_paths,
                             minimum_length,
                             cover_fragments,
                             cover_node_to_fragment,
                             preferred_intervals);

            // recurse on the children
            const vector<const Snarl*>& children = snarl_manager->children_of(cur_snarl);
            for (const Snarl* child_snarl : children) {
                queue.push_back(child_snarl);
            }
        }
    });

    // merge up the thread covers
    vector<RGFAFragment> rgfa_fragments = std::move(thread_covers.at(0));
    for (size_t t = 1; t < thread_count; ++t) {
        rgfa_fragments.reserve(rgfa_fragments.size() + thread_covers.at(t).size());
        // adjust the offsets into each vector
        for (auto& other_frag : thread_covers.at(t)) {
            other_frag.start_parent_idx += rgfa_fragments.size();
            other_frag.end_parent_idx += rgfa_fragments.size();
        }
        std::move(thread_covers.at(t).begin(), thread_covers.at(t).end(), std::back_inserter(rgfa_fragments));
    }
    thread_covers.clear();

    
    // we don't have a path position interface, and even if we did we probably wouldn't have it on every path
    // so to keep running time linear, we need to index the fragments so their offsets can be computed in one scan
    // begin by sorting by path
    unordered_map<path_handle_t, vector<int64_t>> path_to_fragments;
    for (size_t i = 0; i <rgfa_fragments.size(); ++i) {
        const auto& rgfa_fragment = rgfa_fragments[i];
        path_handle_t path_handle = graph->get_path_handle_of_step(rgfa_fragment.steps.front());
        path_to_fragments[path_handle].push_back(i);
    }

    unordered_map<step_handle_t, int64_t> step_to_pos;

    for (const auto& path_fragments : path_to_fragments) {
        const path_handle_t& path_handle = path_fragments.first;
        string path_name = graph->get_path_name(path_handle);

        const vector<int64_t>& fragments = get<1>(path_fragments);

        // for each path, start by finding the positional offset of all relevant steps in the path by brute-force scann
        int64_t set_count = 0;
        for (const int64_t& frag_idx : fragments) {
            const vector<step_handle_t>& rgfa_fragment = rgfa_fragments.at(frag_idx).steps;
            for (const step_handle_t& step: rgfa_fragment) {
                step_to_pos[step] = -1;
                assert(graph->get_path_handle_of_step(step) == path_handle);
                ++set_count;
            }
        }
        size_t pos = 0;
        graph->for_each_step_in_path(path_handle, [&](const step_handle_t& step_handle) {
            if (step_to_pos.count(step_handle)) {
                step_to_pos[step_handle] = pos;
                --set_count;
            }                
            handle_t handle = graph->get_handle_of_step(step_handle);                
            pos += graph->get_length(handle);
            return true;
        });
        //assert(set_count == 0);
    }

    for (const auto& path_fragments : path_to_fragments) {
        const path_handle_t& path_handle = path_fragments.first;
        string path_name = graph->get_path_name(path_handle);

        const vector<int64_t>& fragments = get<1>(path_fragments);
        
        // second pass to make the path fragments, now that we know the positional offsets of their endpoints
        for (const int64_t frag_idx : fragments) {
            const RGFAFragment& rgfa_fragment = rgfa_fragments.at(frag_idx);
            // don't do anything for rank-0, which we only kept to help with metadata for other fragments
            // todo: can we upstream this check?
            if (rgfa_fragment.rank == 0) {
                continue;
            }
            
            size_t rgfa_frag_pos = step_to_pos[rgfa_fragment.steps.front()];
            size_t rgfa_frag_length = 0;
            for (const step_handle_t& step : rgfa_fragment.steps) {
                rgfa_frag_length += graph->get_length(graph->get_handle_of_step(step));
            }
            subrange_t rgfa_frag_subrange = graph->get_subrange(path_handle);
            rgfa_frag_subrange.first = rgfa_frag_pos + (rgfa_frag_subrange != PathMetadata::NO_SUBRANGE ? rgfa_frag_subrange.first : 0);
            rgfa_frag_subrange.second = rgfa_frag_subrange.first + rgfa_frag_length;


            string start_parent_rgfa_path_name = "";
            int64_t start_parent_offset = -1;                   
            int64_t start_parent_node_id = -1;
            bool start_parent_node_reversed = false;
            
            // get information about the start parent
            if (rgfa_fragment.start_parent_idx >= 0) {
                const RGFAFragment& start_parent_fragment = rgfa_fragments.at(rgfa_fragment.start_parent_idx);
                path_handle_t start_parent_path_handle = graph->get_path_handle_of_step(start_parent_fragment.steps.front());
                start_parent_rgfa_path_name = graph->get_path_name(start_parent_path_handle);
                start_parent_offset = step_to_pos.at(rgfa_fragment.start_parent_step);
                subrange_t start_parent_subrange = PathMetadata::parse_subrange(start_parent_rgfa_path_name);
                if (start_parent_subrange != PathMetadata::NO_SUBRANGE) {
                    // now we subset the path to just the fragment interval (which I think it more useful)                
                    int64_t start_parent_range_start = step_to_pos.at(start_parent_fragment.steps.front());
                    int64_t start_parent_range_end = start_parent_range_start + start_parent_fragment.length;
                    clamp_path_subrange(start_parent_rgfa_path_name, start_parent_range_start, start_parent_range_end);
                    start_parent_offset -= start_parent_range_start;
                }
                start_parent_node_id = graph->get_id(graph->get_handle_of_step(start_parent_fragment.start_parent_step));
                start_parent_node_reversed = graph->get_is_reverse(graph->get_handle_of_step(start_parent_fragment.start_parent_step));
            }

            string end_parent_rgfa_path_name = "";
            int64_t end_parent_offset = -1;                   
            int64_t end_parent_node_id = -1;
            bool end_parent_node_reversed = false;
            
            // get information about the end parent
            if (rgfa_fragment.end_parent_idx >= 0) {
                const RGFAFragment& end_parent_fragment = rgfa_fragments.at(rgfa_fragment.end_parent_idx);
                path_handle_t end_parent_path_handle = graph->get_path_handle_of_step(end_parent_fragment.steps.front());
                end_parent_rgfa_path_name = graph->get_path_name(end_parent_path_handle);
                assert(graph->get_path_handle_of_step(rgfa_fragment.end_parent_step) == end_parent_path_handle);                
                end_parent_offset = step_to_pos.at(rgfa_fragment.end_parent_step);
                subrange_t end_parent_subrange = PathMetadata::parse_subrange(end_parent_rgfa_path_name);
                if (end_parent_subrange != PathMetadata::NO_SUBRANGE) {
                    // now we subset the path to just the fragment interval (which I think it more useful)                
                    int64_t end_parent_range_start = step_to_pos.at(end_parent_fragment.steps.front());
                    int64_t end_parent_range_end = end_parent_range_start + end_parent_fragment.length;
                    clamp_path_subrange(end_parent_rgfa_path_name, end_parent_range_start, end_parent_range_end);
                    end_parent_offset -= end_parent_range_start;
                }                
                end_parent_node_id = graph->get_id(graph->get_handle_of_step(end_parent_fragment.end_parent_step));
                end_parent_node_reversed = graph->get_is_reverse(graph->get_handle_of_step(end_parent_fragment.end_parent_step));
            }
            
            string rgfa_frag_name = create_rgfa_path_name(path_name, rgfa_fragment.rank, rgfa_frag_subrange, rgfa_sample_name,
                                                          start_parent_rgfa_path_name,
                                                          start_parent_offset,
                                                          start_parent_node_id, start_parent_node_reversed,
                                                          end_parent_rgfa_path_name,
                                                          end_parent_offset,
                                                          end_parent_node_id, end_parent_node_reversed);

#ifdef debug
#pragma omp critical(cerr)
            cerr << "making new rgfa fragment with name " << rgfa_frag_name << " and " << rgfa_fragment.steps.size() << " steps. subrange "
                 << rgfa_frag_subrange.first << "," << rgfa_frag_subrange.second << endl;
#endif
            path_handle_t rgfa_fragment_handle = graph->create_path_handle(rgfa_frag_name);
            for (const step_handle_t& step : rgfa_fragment.steps) {
                graph->append_step(rgfa_fragment_handle, graph->get_handle_of_step(step));
            }            
        }
    }

    // forwardize the graph
    rgfa_forwardize_paths(graph, reference_paths);
}

void rgfa_snarl_cover(const PathHandleGraph* graph,
                      const Snarl& snarl,
                      PathTraversalFinder& path_trav_finder,
                      const unordered_set<path_handle_t>& reference_paths,
                      int64_t minimum_length,
                      vector<RGFAFragment>& cover_fragments,
                      unordered_map<nid_t, int64_t>& cover_node_to_fragment,
                      const unordered_map<string, vector<pair<int64_t, int64_t>>>& preferred_intervals) {

#ifdef debug
#pragma omp critical(cerr)
    cerr << "calling rgfa_snarl_cover on " << pb2json(snarl) << endl;
#endif
    
    // // start by finding the path traversals through the snarl
    vector<vector<step_handle_t>> travs;
    {
        pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > path_travs = path_trav_finder.find_path_traversals(snarl);
        travs.reserve(path_travs.first.size());
        
        // reduce protobuf usage by going back to vector of steps instead of keeping SnarlTraversals around
        for (int64_t i = 0; i < path_travs.first.size(); ++i) {
            string trav_path_name = graph->get_path_name(graph->get_path_handle_of_step(path_travs.second[i].first));
            if (get_rgfa_rank(trav_path_name) > 0) {
                // we ignore existing (off-reference) rGFA paths
                // todo: shoulgd there be better error handling?                
                cerr << "Warning : skipping existing rgfa traversal " << trav_path_name << endl;
                continue;
            }
            bool reversed = false;
            if (graph->get_is_reverse(graph->get_handle_of_step(path_travs.second[i].first)) != snarl.start().backward()) {
                reversed = true;
            }                
            assert((graph->get_is_reverse(graph->get_handle_of_step(path_travs.second[i].second)) != snarl.end().backward()) == reversed);
            vector<step_handle_t> trav;
            trav.reserve(path_travs.first[i].visit_size());
            bool done = false;
            function<step_handle_t(step_handle_t)> visit_next_step = [&graph,&reversed](step_handle_t step_handle) {
                return reversed ? graph->get_previous_step(step_handle) : graph->get_next_step(step_handle);
            };
            for (step_handle_t step_handle = path_travs.second[i].first; !done; step_handle = visit_next_step(step_handle)) {
                trav.push_back(step_handle);
                if (step_handle == path_travs.second[i].second) {
                    done = true;
                }
            }
            if (reversed) {
                std::reverse(trav.begin(), trav.end());
            }
            travs.push_back(trav);
        }
    }    

    // find all reference paths through the snarl
    map<string, int64_t> ref_paths;    
    for (int64_t i = 0; i < travs.size(); ++i) {
        path_handle_t trav_path = graph->get_path_handle_of_step(travs[i].front());
        if (reference_paths.count(trav_path)) {
            ref_paths[graph->get_path_name(trav_path)] = i;
        }
    }

    // note: checking both snarl endpoint here is actually necessary: if the reference path doesn't end in a tip,
    // you can end up with a trivial snarl at its end which will crash on this test. 
    if (ref_paths.empty() && (!cover_node_to_fragment.count(snarl.start().node_id()) || !cover_node_to_fragment.count(snarl.end().node_id()))) {      
        // we're not nested in a reference snarl, and we have no reference path
        // by the current logic, there's nothing to be done.
        cerr << "[rgfa] warning: No referene path through snarl " 
             << pb2json(snarl) << ": unable to process for rGFA cover" << endl;
        return;
    }

    if (ref_paths.size() > 1) {
        // And we could probably cope with this... but don't for now
        cerr << "[rgfa] error: Mutiple reference path traversals found through snarl " 
             << pb2json(snarl) << endl;
    }

    if (!ref_paths.empty()) {
        // update the cover: note that we won't actually make a path out of this
        // since we leave rank-0 paths the way they are, but it's handy to have
        // a consistent representation for the cover while linking up the
        // various nesting metadata.
        int64_t ref_trav_length = 0;
        vector<step_handle_t>& ref_trav = travs.at(ref_paths.begin()->second);
        for (step_handle_t ref_step_handle : ref_trav) {
            nid_t node_id = graph->get_id(graph->get_handle_of_step(ref_step_handle));
            cover_node_to_fragment[node_id] = cover_fragments.size();
            ref_trav_length += graph->get_length(graph->get_handle_of_step(ref_step_handle));
        }

        // todo: check endpoints (and remove from ref_trav??)
        RGFAFragment frag = {            
            0, ref_trav, ref_trav_length, -1, -1, ref_trav.front(), ref_trav.back()
        };
        cover_fragments.push_back(frag);
    }

#ifdef debug
#pragma omp critical(cerr)
    cerr << "found " << travs.size() << " traversals including " << ref_paths.size() << " reference traversals" << endl;
#endif


    // find all intervals within a snarl traversal that are completely uncovered.
    // the returned intervals are open-ended.
    function<vector<pair<int64_t, int64_t>>(const vector<step_handle_t>&)> get_uncovered_intervals = [&](const vector<step_handle_t>& trav) {
        vector<pair<int64_t, int64_t>> intervals;
        int64_t start = -1;
        for (size_t i = 0; i < trav.size(); ++i) {
            bool covered = cover_node_to_fragment.count(graph->get_id(graph->get_handle_of_step(trav[i])));
            if (covered) {
                if (start != -1) {
                    intervals.push_back(make_pair(start, i));
                }
                start = -1;
            } else {
                if (start == -1) {
                    start = i;
                }
            }
        }
        if (start != -1) {
            intervals.push_back(make_pair(start, trav.size()));
        }
        return intervals;
    };

    // build an initial ranked list of candidate traversal fragments
    vector<pair<tuple<int64_t, int64_t, int64_t>, pair<int64_t, pair<int64_t, int64_t>>>> ranked_trav_fragments;
    for (int64_t trav_idx = 0; trav_idx < travs.size(); ++trav_idx) {
        // todo: this map seems backwards?  note really a big deal since
        // we're only allowing one element
        bool is_ref = false;
        for (const auto& ref_trav : ref_paths) {
            if (ref_trav.second == trav_idx) {
                is_ref = true;
                break;
            }
        }
        if (is_ref) {
            continue;
        }
        const vector<step_handle_t>& trav = travs.at(trav_idx);
        vector<pair<int64_t, int64_t>> uncovered_intervals = get_uncovered_intervals(trav);

        for (const auto& uncovered_interval : uncovered_intervals) {
            unordered_set<nid_t> cycle_check;
            bool cyclic = false;
            int64_t interval_length = 0;            
            for (int64_t i = uncovered_interval.first; i < uncovered_interval.second && !cyclic; ++i) {
                handle_t handle = graph->get_handle_of_step(trav[i]);
                interval_length += graph->get_length(handle);
                nid_t node_id = graph->get_id(handle);
                if (cycle_check.count(node_id)) {
                    cyclic = true;
                } else {
                    cycle_check.insert(node_id);
                }                
            }
            if (!cyclic && interval_length >= minimum_length) {
                auto trav_stats = rgfa_traversal_stats(graph, trav, uncovered_interval);
                ranked_trav_fragments.push_back(make_pair(trav_stats, make_pair(trav_idx, uncovered_interval)));
            }
        }
    }

    // todo: typedef!
    function<bool(const pair<tuple<int64_t, int64_t, int64_t>, pair<int64_t, pair<int64_t, int64_t>>>& s1,
                  const pair<tuple<int64_t, int64_t, int64_t>, pair<int64_t, pair<int64_t, int64_t>>>& s2)> heap_comp =
        [](const pair<tuple<int64_t, int64_t, int64_t>, pair<int64_t, pair<int64_t, int64_t>>>& s1,
           const pair<tuple<int64_t, int64_t, int64_t>, pair<int64_t, pair<int64_t, int64_t>>>& s2) {
            return rgfa_traversal_stats_less(s1.first, s2.first);
        };

    // put the fragments into a max heap
    std::make_heap(ranked_trav_fragments.begin(), ranked_trav_fragments.end(), heap_comp);

    // now greedily pull out traversal intervals from the ranked list until none are left
    while (!ranked_trav_fragments.empty()) {

        // get the best scoring (max) fragment from heap
        auto best_stats_fragment = ranked_trav_fragments.front();
        std::pop_heap(ranked_trav_fragments.begin(), ranked_trav_fragments.end(), heap_comp);
        ranked_trav_fragments.pop_back();
        
        const vector<step_handle_t>& trav = travs.at(best_stats_fragment.second.first);
        const pair<int64_t, int64_t>& uncovered_interval = best_stats_fragment.second.second;

        // our traversal may have been partially covered by a different iteration, if so, we need to break it up
        // and continue
        vector<pair<int64_t, int64_t>> chopped_intervals;
        int64_t cur_start = -1;
        bool chopped = false;
        for (int64_t i = uncovered_interval.first; i < uncovered_interval.second; ++i) {
            bool covered = cover_node_to_fragment.count(graph->get_id(graph->get_handle_of_step(trav[i])));
            if (!covered && cur_start == -1) {
                cur_start = i;
            } else if (covered) {
                chopped = true;
                if (cur_start != -1) {
                    chopped_intervals.push_back(make_pair(cur_start, i));
                    cur_start = -1;
                }
            }
        }
        if (cur_start != -1) {
            chopped_intervals.push_back(make_pair(cur_start, uncovered_interval.second));
        }
        if (chopped) {
            for (const pair<int64_t, int64_t>& chopped_interval : chopped_intervals) {
                int64_t chopped_trav_length = 0;
                for (int64_t i = chopped_interval.first; i < chopped_interval.second; ++i) {
                    chopped_trav_length += graph->get_length(graph->get_handle_of_step(trav[i]));
                }
                if (chopped_trav_length >= minimum_length) {
                    auto chopped_stats = rgfa_traversal_stats(graph, trav, chopped_interval);                
                    ranked_trav_fragments.push_back(make_pair(chopped_stats, make_pair(best_stats_fragment.second.first, chopped_interval)));
                    std::push_heap(ranked_trav_fragments.begin(), ranked_trav_fragments.end(), heap_comp);
                }
            }
            continue;
        }

        // we check the "backbone" interval that this interval is coming off
        // since we've already covered the reference, then we know that this interval
        // doesn't span the whole snarl including endpoints, so we can always afford
        // to look back and ahead one
        assert(uncovered_interval.first > 0 && uncovered_interval.second < trav.size());
        int64_t prev_frag_idx = cover_node_to_fragment.at(graph->get_id(graph->get_handle_of_step(trav[uncovered_interval.first - 1])));
        int64_t next_frag_idx = cover_node_to_fragment.at(graph->get_id(graph->get_handle_of_step(trav[uncovered_interval.second])));
        // todo: i'm not sure if/how minigraph treats these cases, where the anchors connect to different ranks
        //       also, can these be avoided entirely?
        int64_t min_frag_idx = std::min(prev_frag_idx, next_frag_idx);
        int64_t fragment_rank;
        if (min_frag_idx == -1) {
            fragment_rank = 1;
        } else {
            fragment_rank = cover_fragments.at(min_frag_idx).rank + 1;
        }

        // now we need to find the steps on the parent path, in order to link back to its position
        // todo: can we avoid this enumeration by keeping track of the parent directly from the beginning?         
        handle_t start_handle = graph->get_handle_of_step(trav[uncovered_interval.first - 1]);
        vector<step_handle_t> start_handle_steps = graph->steps_of_handle(start_handle);
        path_handle_t start_parent_path = graph->get_path_handle_of_step(cover_fragments.at(prev_frag_idx).steps.front());
        vector<int64_t> start_parent_step_indexes;
        for (int64_t i = 0; i < start_handle_steps.size(); ++i) {
            if (graph->get_path_handle_of_step(start_handle_steps[i]) == start_parent_path) {
                start_parent_step_indexes.push_back(i);
            }
        }
        assert(start_parent_step_indexes.size() == 1);
        step_handle_t start_parent_step = start_handle_steps[start_parent_step_indexes[0]];

#ifdef debug
#pragma omp critical(cerr)        
        cerr << "adding step " << graph->get_path_name(start_parent_path) << " " << graph->get_id(graph->get_handle_of_step(start_parent_step)) << " as start parent for " << graph->get_path_name(graph->get_path_handle_of_step(trav.front())) << " starting at " << graph->get_id(graph->get_handle_of_step(trav.at(uncovered_interval.first))) << endl;
#endif
 
        handle_t end_handle = graph->get_handle_of_step(trav[uncovered_interval.second]);
        vector<step_handle_t> end_handle_steps = graph->steps_of_handle(end_handle);
        path_handle_t end_parent_path = graph->get_path_handle_of_step(cover_fragments.at(prev_frag_idx).steps.front());
        vector<int64_t> end_parent_step_indexes;
        for (int64_t i = 0; i < end_handle_steps.size(); ++i) {
            if (graph->get_path_handle_of_step(end_handle_steps[i]) == end_parent_path) {
                end_parent_step_indexes.push_back(i);
            }
        }
        assert(end_parent_step_indexes.size() == 1);
        step_handle_t end_parent_step = end_handle_steps[end_parent_step_indexes[0]];                    

        // update the cover
        vector<step_handle_t> interval;
        int64_t interval_length = 0;
        interval.reserve(uncovered_interval.second - uncovered_interval.first);
        for (int64_t i = uncovered_interval.first; i < uncovered_interval.second; ++i) {
            interval.push_back(trav[i]);
            interval_length += graph->get_length(graph->get_handle_of_step(trav[i]));
            cover_node_to_fragment[graph->get_id(graph->get_handle_of_step(trav[i]))] = cover_fragments.size();
        }        
        RGFAFragment frag = {
            fragment_rank, std::move(interval), interval_length,
            prev_frag_idx, next_frag_idx,
            start_parent_step, end_parent_step
        };

#ifdef debug
#pragma omp critical(cerr)
{
        cerr << "adding fragment path=" << graph->get_path_name(graph->get_path_handle_of_step(frag.steps.front())) << " rank=" << fragment_rank << " steps = ";
        for (auto xx : frag.steps) {
            cerr << graph->get_id(graph->get_handle_of_step(xx)) <<":" << graph->get_is_reverse(graph->get_handle_of_step(xx)) <<",";
        }
        cerr << " sps " <<  graph->get_id(graph->get_handle_of_step(start_parent_step)) <<":" << graph->get_is_reverse(graph->get_handle_of_step(start_parent_step));
        cerr << " eps " <<  graph->get_id(graph->get_handle_of_step(end_parent_step)) <<":" << graph->get_is_reverse(graph->get_handle_of_step(end_parent_step));
        cerr << endl;
}
#endif
        cover_fragments.push_back(frag);
    }
}

tuple<int64_t, int64_t, int64_t> rgfa_traversal_stats(const PathHandleGraph* graph,
                                                      const vector<step_handle_t>& trav,
                                                      const pair<int64_t, int64_t>& trav_fragment) {
    path_handle_t path_handle = graph->get_path_handle_of_step(trav.front());
    int64_t support = 0;
    int64_t reversed_steps = 0;
    int64_t dupes = 0;

    for (int64_t i = trav_fragment.first; i < trav_fragment.second; ++i) {
        const step_handle_t& step = trav[i];
        handle_t handle = graph->get_handle_of_step(step);
        vector<step_handle_t> all_steps = graph->steps_of_handle(handle);
        int64_t length = graph->get_length(handle);
        support += length;
        int64_t self_count = 0;
        for (const step_handle_t& other_step : all_steps) {
            path_handle_t step_path_handle = graph->get_path_handle_of_step(other_step);
            if (step_path_handle == path_handle) {
                ++self_count;
            } else {
                support += length;
            }
        }
        if (self_count > 1) {
            dupes += length * (self_count - 1);
        }
        if (i > 0 && graph->get_is_reverse(handle)) {
            ++reversed_steps;
        }
    }

    return std::make_tuple(support, reversed_steps, dupes);
}

bool rgfa_traversal_stats_less(const tuple<int64_t, int64_t, int64_t>& s1, const tuple<int64_t, int64_t, int64_t>& s2) {
    // duplicates are a deal breaker, if one traversal has no duplicates and the other does, the former wins
    if (get<2>(s1) > 0 && get<2>(s2) == 0) {
        return true;
    } else if (get<2>(s1) == 0 && get<2>(s2) > 0) {
        return false;
    }

    // then support
    if (get<0>(s1) < get<0>(s2)) {
        return true;
    } else if (get<0>(s1) > get<0>(s2)) {
        return false;
    }

    // then duplicates (by value)
    if (get<2>(s1) > get<2>(s2)) {
        return true;
    } else if (get<2>(s1) < get<2>(s2)) {
        return false;
    }

    // then switches
    return get<1>(s1) > get<1>(s2);
}

// copied pretty much verbatem from
// https://github.com/ComparativeGenomicsToolkit/hal2vg/blob/v1.1.2/clip-vg.cpp#L809-L880
void rgfa_forwardize_paths(MutablePathMutableHandleGraph* graph,
                           const unordered_set<path_handle_t>& reference_paths) {

    unordered_map<nid_t, nid_t> id_map;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            if (reference_paths.count(path_handle) || get_rgfa_rank(path_name) >= 0) {
                size_t fw_count = 0;
                size_t total_steps = 0;
                graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        handle_t handle = graph->get_handle_of_step(step_handle);
                        if (graph->get_is_reverse(handle)) {
                            handle_t flipped_handle = graph->create_handle(graph->get_sequence(handle));
                            id_map[graph->get_id(flipped_handle)] = graph->get_id(handle);
			    graph->follow_edges(handle, true, [&](handle_t prev_handle) {
                                    if (graph->get_id(prev_handle) != graph->get_id(handle)) {
                                        graph->create_edge(prev_handle, flipped_handle);
                                    }
			      });
			    graph->follow_edges(handle, false, [&](handle_t next_handle) {
                                    if (graph->get_id(handle) != graph->get_id(next_handle)) {
                                        graph->create_edge(flipped_handle, next_handle);
                                    }
			      });
                            // self-loop cases we punted on above:
                            if (graph->has_edge(handle, handle)) {
                                graph->create_edge(flipped_handle, flipped_handle);
                            }
                            if (graph->has_edge(handle, graph->flip(handle))) {
                                graph->create_edge(flipped_handle, graph->flip(flipped_handle));                                
                            }
                            if (graph->has_edge(graph->flip(handle), handle)) {
                                graph->create_edge(graph->flip(flipped_handle), flipped_handle);
                            }
                            vector<step_handle_t> steps = graph->steps_of_handle(handle);
                            size_t ref_count = 0;
                            for (step_handle_t step : steps) {
                                if (graph->get_path_handle_of_step(step) == path_handle) {
                                    ++ref_count;
                                }
                                step_handle_t next_step = graph->get_next_step(step);
                                handle_t new_handle = graph->get_is_reverse(graph->get_handle_of_step(step)) ? flipped_handle :
                                    graph->flip(flipped_handle);
                                graph->rewrite_segment(step, next_step, {new_handle});
                            }
                            if (ref_count > 1) {
                                cerr << "[rGFA] error: Cycle detected in rGFA path " << path_name << " at node " << graph->get_id(handle) << endl;
                                exit(1);
                            }
                            ++fw_count;
                            assert(graph->steps_of_handle(handle).empty());
                            dynamic_cast<DeletableHandleGraph*>(graph)->destroy_handle(handle);
                        }
                        ++total_steps;
                    });
            }
        });

    // rename all the ids back to what they were (so nodes keep their ids, just get flipped around)
    graph->reassign_node_ids([&id_map](nid_t new_id) {
        return id_map.count(new_id) ? id_map[new_id] : new_id;
    });

    // do a check just to be sure
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            if (reference_paths.count(path_handle) || get_rgfa_rank(path_name) >= 0) {
                graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        handle_t handle = graph->get_handle_of_step(step_handle);
                        if (graph->get_is_reverse(handle)) {
                            cerr << "[rGFA] error: Failed to fowardize node " << graph->get_id(handle) << " in path " << path_name << endl;
                            exit(1);
                        }
                    });
            }
        });
}


}

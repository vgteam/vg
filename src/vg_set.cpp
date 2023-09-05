#include "vg_set.hpp"
#include <vg/io/stream.hpp>
#include "source_sink_overlay.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "io/save_handle_graph.hpp"

namespace vg {
// sets of MutablePathMutableHandleGraphs on disk

void VGset::transform(std::function<void(MutableHandleGraph*)> lambda) {
    // TODO: add a way to cache graphs here for multiple scans
    for (auto& name : filenames) {
        // load
        unique_ptr<MutableHandleGraph> g;
        get_input_file(name, [&](istream& in) {
            // Note: I would have liked to just load a MutableHandleGraph here but the resulting pointer
            // is broken (tested: VG and PackedGraph)
            g = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(in);
            });
        // legacy:
        VG* vg_g = dynamic_cast<VG*>(g.get());
        if (vg_g != nullptr) {
            vg_g->name = name;
        }
        // apply
        lambda(g.get());
        // write to the same file
        vg::io::save_handle_graph(g.get(), name);
    }
}

void VGset::for_each(std::function<void(HandleGraph*)> lambda) {
    // TODO: add a way to cache graphs here for multiple scans
    for (auto& name : filenames) {
        // load
        unique_ptr<HandleGraph> g;
        get_input_file(name, [&](istream& in) {
                g = vg::io::VPKG::load_one<HandleGraph>(in);
            });
        // legacy:
        VG* vg_g = dynamic_cast<VG*>(g.get());
        if (vg_g != nullptr) {
            vg_g->name = name;
        }        
        // apply
        lambda(g.get());
    }
}

id_t VGset::max_node_id(void) {
    id_t max_id = 0;
    for_each([&](HandleGraph* graph) {
            max_id = max(graph->max_node_id(), max_id);
        });
    return max_id;
}

int64_t VGset::merge_id_space(void) {
    int64_t max_node_id = 0;
    auto lambda = [&max_node_id](MutableHandleGraph* g) {
        int64_t delta = max_node_id - g->min_node_id();
        if (delta >= 0) {
            g->increment_node_ids(delta + 1);
        }
        max_node_id = g->max_node_id();
    };
    transform(lambda);
    return max_node_id;
}

void VGset::to_xg(xg::XG& index) {
    // Send a predicate to match nothing
    to_xg(index, [](const string& ignored) {
        return false;
    });
}

void VGset::to_xg(xg::XG& index, const function<bool(const string&)>& paths_to_remove, map<string, Path>* removed_paths) {

    // We make multiple passes through the input graphs, loading each and going through its nodes/edges/paths.
    // TODO: This is going to go load all the graphs multiple times, even if there's only one graph!
    // TODO: streaming HandleGraph API?
    // TODO: push-driven XG build?
    // TODO: XG::from_path_handle_graphs?
    
    // We need to recostruct full removed paths from fragmentary paths encountered in each chunk.
    // This maps from path name to all the Mappings in the path in the order we encountered them
    auto for_each_sequence = [&](const std::function<void(const std::string& seq, const nid_t& node_id)>& lambda) {
        for_each([&](HandleGraph* graph) {
            // For each graph in the set
            
            // ID-sort its handles. This is the order that we have historically
            // used, and may be required for e.g. vg map to work correctly.
            // TODO: Compute this only once and not each time we want to visit
            // all the nodes during XG construction?
            vector<handle_t> order;
            // TODO: reserve if counting handles will be efficient. Can we know?
            graph->for_each_handle([&](const handle_t& h) {
                order.push_back(h);
            });
            std::sort(order.begin(), order.end(), [&](const handle_t& a, const handle_t& b) {
                // Return true if a must come first
                // We know everything is locally forward already.
                return graph->get_id(a) < graph->get_id(b);
            });
            
            for (const handle_t& h : order) {
                // For each node in the graph, tell the XG about it.
                // Assume it is locally forward.
#ifdef debug
                cerr << "Yield node " << graph->get_id(h) << " sequence " << graph->get_sequence(h) << endl;
#endif
                lambda(graph->get_sequence(h), graph->get_id(h));
            }
        });
    };

    auto for_each_edge = [&](const std::function<void(const nid_t& from, const bool& from_rev, const nid_t& to, const bool& to_rev)>& lambda) {
        for_each([&](HandleGraph* graph) {
            // For each graph in the set
            graph->for_each_edge([&](const edge_t& e) {
                // For each edge in the graph, tell the XG about it
                
#ifdef debug
                cerr << "Yield edge " << graph->get_id(e.first) << " " << graph->get_is_reverse(e.first)
                    << " " << graph->get_id(e.second) << " " << graph->get_is_reverse(e.second) << endl;
#endif
                
                lambda(graph->get_id(e.first), graph->get_is_reverse(e.first), graph->get_id(e.second), graph->get_is_reverse(e.second));
            });
        });
    };
    
    // We no longer need to reconstitute paths ourselves; we require that each
    // input file has all the parts of its paths.
    
    auto for_each_path_element = [&](const std::function<void(const std::string& path_name,
                                     const nid_t& node_id, const bool& is_rev,
                                     const std::string& cigar,
                                     bool is_empty, bool is_circular)>& lambda) {
                                 
        for_each([&](HandleGraph* graph) {
            // Look at each graph and see if it has path support
            // TODO: do we need a different for_each to push this down to the loader?
            PathHandleGraph* path_graph = dynamic_cast<PathHandleGraph*>(graph);
            if (path_graph != nullptr) {
                // The graph we loaded actually can have paths, so it can have visits
                path_graph->for_each_path_handle([&](const path_handle_t& p) {
                    // For each path
                    
                    // Get its metadata
                    string path_name = path_graph->get_path_name(p);
                    bool is_circular = path_graph->get_is_circular(p);
                    
                    if(paths_to_remove(path_name)) {
                        // We want to filter out this path
                        if (removed_paths != nullptr) {
                            // When we filter out a path, we need to send our caller a Protobuf version of it.
                            Path proto_path;
                            proto_path.set_name(path_name);
                            proto_path.set_is_circular(is_circular);
                            size_t rank = 1;
                            path_graph->for_each_step_in_path(p, [&](const step_handle_t& s) {
                                handle_t stepped_on = path_graph->get_handle_of_step(s);
                                
                                Mapping* mapping = proto_path.add_mapping();
                                mapping->mutable_position()->set_node_id(path_graph->get_id(stepped_on));
                                mapping->mutable_position()->set_is_reverse(path_graph->get_is_reverse(stepped_on));
                                mapping->set_rank(rank++);
                            });
                            removed_paths->emplace(std::move(path_name), std::move(proto_path));
                        }
                    } else {
                        // We want to leave this path in
                    
                        // Assume it is empty
                        bool is_empty = true;
                        
                        path_graph->for_each_step_in_path(p, [&](const step_handle_t& s) {
                            // For each visit on the path
                            handle_t stepped_on = path_graph->get_handle_of_step(s);
                            // The path can't be empty
                            is_empty = false;
                            // Tell the XG about it
#ifdef debug
                            cerr << "Yield path " << path_name << " visit to "
                                << path_graph->get_id(stepped_on) << " " << path_graph->get_is_reverse(stepped_on) << endl;
#endif
                            lambda(path_name, path_graph->get_id(stepped_on), path_graph->get_is_reverse(stepped_on), "", is_empty, is_circular); 
                        });
                        
                        if (is_empty) {
                            // If the path is empty, tell the XG that.
                            // It still could be circular.
                            
#ifdef debug
                            cerr << "Yield empty path " << path_name << endl;
#endif
                            
                            lambda(path_name, 0, false, "", is_empty, is_circular);
                        }
                    }
                });
            }
        });
    };
    
    // Now build the xg graph, looping over all our graphs in our set whenever we want anything.
    index.from_enumerators(for_each_sequence, for_each_edge, for_each_path_element, false);
}

void VGset::for_each_kmer_parallel(size_t kmer_size, const function<void(const kmer_t&)>& lambda) {
    for_each([&lambda, kmer_size, this](HandleGraph* g) {
            // legacy
            VG* vg_g = dynamic_cast<VG*>(g);
            if (vg_g != nullptr) {
                vg_g->show_progress = show_progress & progress_bars;
                vg_g->preload_progress("processing kmers of " + vg_g->name);
            }
        //g->for_each_kmer_parallel(kmer_size, path_only, edge_max, lambda, stride, allow_dups, allow_negatives);
        for_each_kmer(*g, kmer_size, lambda);
    });
}

void VGset::write_gcsa_kmers_ascii(ostream& out, int kmer_size,
                                   nid_t head_id, nid_t tail_id) {
    if (filenames.size() > 1 && (head_id == 0 || tail_id == 0)) {
        // Detect head and tail IDs in advance if we have multiple graphs
        nid_t max_id = max_node_id(); // expensive, as we'll stream through all the files
        head_id = max_id + 1;
        tail_id = max_id + 2;
    }

    // When we're sure we know what this kmer instance looks like, we'll write
    // it out exactly once. We need the start_end_id actually used in order to
    // go to the correct place when we don't go anywhere (i.e. at the far end of
    // the start/end node.
    auto write_kmer = [&head_id, &tail_id, &out](const kmer_t& kp){
#pragma omp critical (out)
        out << kp << endl;
    };

    for_each([&](HandleGraph* g) {
        // Make an overlay for each graph, without modifying it. Break into tip-less cycle components.
        // Make sure to use a consistent head and tail ID across all graphs in the set.
        SourceSinkOverlay overlay(g, kmer_size, head_id, tail_id);
        
        // Read back the head and tail IDs in case we have only one graph and we just detected them now.
        head_id = overlay.get_id(overlay.get_source_handle());
        tail_id = overlay.get_id(overlay.get_sink_handle());
        
        // Now get the kmers in the graph that pretends to have single head and tail nodes
        for_each_kmer(overlay, kmer_size, write_kmer, head_id, tail_id);
    });
}

// writes to a specific output stream
void VGset::write_gcsa_kmers_binary(ostream& out, int kmer_size, size_t& size_limit,
                                    nid_t head_id, nid_t tail_id) {
    if (filenames.size() > 1 && (head_id == 0 || tail_id == 0)) {
        // Detect head and tail IDs in advance if we have multiple graphs
        nid_t max_id = max_node_id(); // expensive, as we'll stream through all the files
        head_id = max_id + 1;
        tail_id = max_id + 2;
    }

    size_t total_size = 0;
    for_each([&](HandleGraph* g) {
        // Make an overlay for each graph, without modifying it. Break into tip-less cycle components.
        // Make sure to use a consistent head and tail ID across all graphs in the set.
        SourceSinkOverlay overlay(g, kmer_size, head_id, tail_id);
        
        // Read back the head and tail IDs in case we have only one graph and we just detected them now.
        head_id = overlay.get_id(overlay.get_source_handle());
        tail_id = overlay.get_id(overlay.get_sink_handle());
        
        size_t current_bytes = size_limit - total_size;
        write_gcsa_kmers(overlay, kmer_size, out, current_bytes, head_id, tail_id);
        total_size += current_bytes;
    });
    size_limit = total_size;
}

// writes to a set of temp files and returns their names
vector<string> VGset::write_gcsa_kmers_binary(int kmer_size, size_t& size_limit,
                                              nid_t head_id, nid_t tail_id) {
    if (filenames.size() > 1 && (head_id == 0 || tail_id == 0)) {
        // Detect head and tail IDs in advance if we have multiple graphs
        nid_t max_id = max_node_id(); // expensive, as we'll stream through all the files
        head_id = max_id + 1;
        tail_id = max_id + 2;
    }

    vector<string> tmpnames;
    size_t total_size = 0;
    for_each([&](HandleGraph* g) {
        // Make an overlay for each graph, without modifying it. Break into tip-less cycle components.
        // Make sure to use a consistent head and tail ID across all graphs in the set.
        SourceSinkOverlay overlay(g, kmer_size, head_id, tail_id);
        
        // Read back the head and tail IDs in case we have only one graph and we just detected them now.
        head_id = overlay.get_id(overlay.get_source_handle());
        tail_id = overlay.get_id(overlay.get_sink_handle());
        
        size_t current_bytes = size_limit - total_size;
        try {
            tmpnames.push_back(write_gcsa_kmers_to_tmpfile(overlay, kmer_size, current_bytes, head_id, tail_id));
        }
        catch (SizeLimitExceededException& ex) {
            // clean up the temporary files before continuing to throw the exception
            for (const auto& tmpname : tmpnames) {
                temp_file::remove(tmpname);
            }
            throw ex;
        }
        total_size += current_bytes;
    });
    size_limit = total_size;
    return tmpnames;
}

}

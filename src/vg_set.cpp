#include "vg_set.hpp"
#include <vg/io/stream.hpp>
#include "source_sink_overlay.hpp"

namespace vg {
// sets of VGs on disk

void VGset::transform(std::function<void(VG*)> lambda) {
    for (auto& name : filenames) {
        // load
        VG* g = NULL;
        if (name == "-") {
            g = new VG(std::cin, show_progress & progress_bars);
        } else {
            ifstream in(name.c_str());
            if (!in) throw ifstream::failure("failed to open " + name);
            g = new VG(in, show_progress & progress_bars);
            in.close();
        }
        g->name = name;
        // apply
        lambda(g);
        // write to the same file
        ofstream out(name.c_str());
        g->serialize_to_ostream(out);
        out.close();
        delete g;
    }
}

void VGset::for_each(std::function<void(VG*)> lambda) {
    for (auto& name : filenames) {
        // load
        VG* g = NULL;
        if (name == "-") {
            g = new VG(std::cin, show_progress & progress_bars);
        } else {
            ifstream in(name.c_str());
            if (!in) throw ifstream::failure("failed to open " + name);
            g = new VG(in, show_progress & progress_bars);
            in.close();
        }
        g->name = name;
        // apply
        lambda(g);
        delete g;
    }
}

void VGset::for_each_graph_chunk(std::function<void(Graph&)> lamda) {
    for (auto& name : filenames) {
        ifstream in(name.c_str());
        vg::io::for_each(in, lamda);
    }
}

id_t VGset::max_node_id(void) {
    id_t max_id = 0;
    for_each_graph_chunk([&](const Graph& graph) {
            for (size_t i = 0; i < graph.node_size(); ++i) {
                max_id = max(graph.node(i).id(), max_id);
            }
        });
    return max_id;
}

int64_t VGset::merge_id_space(void) {
    int64_t max_node_id = 0;
    auto lambda = [&max_node_id](VG* g) {
        if (max_node_id > 0) g->increment_node_ids(max_node_id);
        max_node_id = g->max_node_id();
    };
    transform(lambda);
    return max_node_id;
}

void VGset::to_xg(xg::XG& index, bool store_threads) {
    // Send a predicate to match nothing
    to_xg(index, store_threads, [](const string& ignored) {
        return false;
    });
}

void VGset::to_xg(xg::XG& index, bool store_threads, const regex& paths_to_take, map<string, Path>* removed_paths) {
    to_xg(index, store_threads, [&](const string& path_name) -> bool {
        // Take paths that match the regex.
        return std::regex_match(path_name, paths_to_take);
    }, removed_paths);
}

void VGset::to_xg(xg::XG& index, bool store_threads, const function<bool(const string&)>& paths_to_take, map<string, Path>* removed_paths) {
    
    // We need to recostruct full removed paths from fragmentary paths encountered in each chunk.
    // This maps from path name to all the Mappings in the path in the order we encountered them
    auto for_each_sequence = [&](const std::function<void(const std::string& seq, const nid_t& node_id)>& lambda) {
        for (auto& name : filenames) {
            std::ifstream in(name);
            vg::io::for_each(in, (function<void(Graph&)>)[&](Graph& graph) {
                    for (uint64_t i = 0; i < graph.node_size(); ++i) {
                        auto& node = graph.node(i);
                        lambda(node.sequence(), node.id());
                    }
                });
        }
    };

    auto for_each_edge = [&](const std::function<void(const nid_t& from, const bool& from_rev, const nid_t& to, const bool& to_rev)>& lambda) {
        for (auto& name : filenames) {
            std::ifstream in(name);
            vg::io::for_each(in, (function<void(Graph&)>)[&](Graph& graph) {
                    for (uint64_t i = 0; i < graph.edge_size(); ++i) {
                        auto& edge = graph.edge(i);
                        lambda(edge.from(), edge.from_start(), edge.to(), edge.to_end());
                    }
                });
        }
    };

    // thanks to the silliness that is vg's graph chunking, we have to reconstitute our paths here
    // to make sure that they are constructed with the correct order
    map<string, vector<pair<nid_t, bool>>> paths;
    for (auto& name : filenames) {
        std::ifstream in(name);
        vg::io::for_each(in, (function<void(Graph&)>)[&](Graph& graph) {
                for (uint64_t i = 0; i < graph.path_size(); ++i) {
                    auto& path = graph.path(i);
                    auto& steps = paths[path.name()];
                    for (auto& mapping : path.mapping()) {
                        if (mapping.rank() > steps.size()) {
                            steps.resize(mapping.rank());
                        }
                        steps[mapping.rank()-1] = make_pair(mapping.position().node_id(),
                                                            mapping.position().is_reverse());
                    }
                }
            });
    }

    auto for_each_path_element = [&](
        const std::function<void(const std::string& path_name,
                                 const nid_t& node_id, const bool& is_rev,
                                 const std::string& cigar)>& lambda) {
        for (auto& path : paths) {
            auto& path_name = path.first;
            auto& path_steps = path.second;
            for (auto& step : path_steps) {
                if (step.first) { // ids are > 0, so we skip anything that's not filled
                    lambda(path_name, step.first, step.second, "");
                }
            }
        }
    };

    index.from_enumerators(for_each_sequence, for_each_edge, for_each_path_element, false);

}

void VGset::for_each_kmer_parallel(size_t kmer_size, const function<void(const kmer_t&)>& lambda) {
    for_each([&lambda, kmer_size, this](VG* g) {
        g->show_progress = show_progress & progress_bars;
        g->preload_progress("processing kmers of " + g->name);
        //g->for_each_kmer_parallel(kmer_size, path_only, edge_max, lambda, stride, allow_dups, allow_negatives);
        for_each_kmer(*g, kmer_size, lambda);
    });
}

void VGset::write_gcsa_kmers_ascii(ostream& out, int kmer_size,
                                   id_t head_id, id_t tail_id) {
    if (filenames.size() > 1 && (head_id == 0 || tail_id == 0)) {
        // Detect head and tail IDs in advance if we have multiple graphs
        id_t max_id = max_node_id(); // expensive, as we'll stream through all the files
        head_id = max_id + 1;
        tail_id = max_id + 2;
    }

    // When we're sure we know what this kmer instance looks like, we'll write
    // it out exactly once. We need the start_end_id actually used in order to
    // go to the correct place when we don't go anywhere (i.e. at the far end of
    // the start/end node.
    auto write_kmer = [&head_id, &tail_id](const kmer_t& kp){
#pragma omp critical (cout)
        cout << kp << endl;
    };

    for_each([&](VG* g) {
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
                                    id_t head_id, id_t tail_id) {
    if (filenames.size() > 1 && (head_id == 0 || tail_id == 0)) {
        // Detect head and tail IDs in advance if we have multiple graphs
        id_t max_id = max_node_id(); // expensive, as we'll stream through all the files
        head_id = max_id + 1;
        tail_id = max_id + 2;
    }

    size_t total_size = 0;
    for_each([&](VG* g) {
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
                                              id_t head_id, id_t tail_id) {
    if (filenames.size() > 1 && (head_id == 0 || tail_id == 0)) {
        // Detect head and tail IDs in advance if we have multiple graphs
        id_t max_id = max_node_id(); // expensive, as we'll stream through all the files
        head_id = max_id + 1;
        tail_id = max_id + 2;
    }

    vector<string> tmpnames;
    size_t total_size = 0;
    for_each([&](VG* g) {
        // Make an overlay for each graph, without modifying it. Break into tip-less cycle components.
        // Make sure to use a consistent head and tail ID across all graphs in the set.
        SourceSinkOverlay overlay(g, kmer_size, head_id, tail_id);
        
        // Read back the head and tail IDs in case we have only one graph and we just detected them now.
        head_id = overlay.get_id(overlay.get_source_handle());
        tail_id = overlay.get_id(overlay.get_sink_handle());
        
        size_t current_bytes = size_limit - total_size;
        tmpnames.push_back(write_gcsa_kmers_to_tmpfile(overlay, kmer_size, current_bytes, head_id, tail_id));
        total_size += current_bytes;
    });
    size_limit = total_size;
    return tmpnames;
}

}

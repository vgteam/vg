#include "vg_set.hpp"
#include "stream.hpp"
#include "source_sink_overlay.hpp"

namespace vg {
// sets of VGs on disk

void VGset::transform(std::function<void(VG*)> lambda) {
    for (auto& name : filenames) {
        // load
        VG* g = NULL;
        if (name == "-") {
            g = new VG(std::cin, show_progress);
        } else {
            ifstream in(name.c_str());
            if (!in) throw ifstream::failure("failed to open " + name);
            g = new VG(in, show_progress);
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
            g = new VG(std::cin, show_progress);
        } else {
            ifstream in(name.c_str());
            if (!in) throw ifstream::failure("failed to open " + name);
            g = new VG(in, show_progress);
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
        stream::for_each(in, lamda);
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
    // Nothing matches the default-constructed regex, so nothing will ever be
    // sent to the map.
    to_xg(index, store_threads, regex());
}

void VGset::to_xg(xg::XG& index, bool store_threads, const regex& paths_to_take, map<string, Path>* removed_paths) {
    
    // We need to recostruct full removed paths from fragmentary paths encountered in each chunk.
    // This maps from path name to all the Mappings in the path in the order we encountered them
    map<string, list<Mapping>> mappings;
    
    // Set up an XG index
    index.from_callback([&](function<void(Graph&)> callback) {
        for (auto& name : filenames) {
#ifdef debug
            cerr << "Loading chunks from " << name << endl;
#endif
            // Load chunks from all the files and pass them into XG.
            std::ifstream in(name);
            
            if (name == "-"){
                if (!in) throw ifstream::failure("vg_set: cannot read from stdin. Failed to open " + name);
            }
            
            if (!in) throw ifstream::failure("failed to open " + name);
            
            function<void(Graph&)> handle_graph = [&](Graph& graph) {
#ifdef debug
                cerr << "Got chunk of " << name << "!" << endl;
#endif

                // We'll move all the path fragments into one of these (if removed_paths is not null)
                std::list<Path> paths_taken;

                // Remove the matching paths.
                remove_paths(graph, paths_to_take, removed_paths ? &paths_taken : nullptr);

                for (auto& path : paths_taken) {
                    // Copy all the mappings from each path into the collection of mappings in order encountered
                    std::copy(path.mapping().begin(), path.mapping().end(), std::back_inserter(mappings[path.name()]));
                }

                // Ship out the corrected graph
                callback(graph);
            };
            
            stream::for_each(in, handle_graph);
            
            // Now that we got all the chunks, reconstitute any siphoned-off paths into Path objects and return them.
            // We have to handle chunks being encountered in any order, if ranks are set, or in path-forward order, if ranks are missing.
            for(auto& kv : mappings) {
                // We'll fill in this Path object
                Path path;
                path.set_name(kv.first);

                // This will hold mappings by rank.
                // If they don't have ranks assigned, we will assign them in the order encountered.
                map<int64_t, Mapping> mappings_by_rank;

                for (auto& mapping : kv.second) {
                    if (mapping.rank() != 0) {
                        // It has a rank already.
                        // Make sure we didn't try to assign something to its rank
                        assert(!mappings_by_rank.count(mapping.rank()));

                        mappings_by_rank[mapping.rank()] = mapping;
                    } else {
                        // Assign it the next available rank in its path
                        int64_t next_rank;
                        if (!mappings_by_rank.empty()) {
                            next_rank = mappings_by_rank.rbegin()->first + 1;
                        } else {
                            next_rank = 1;
                        }

                        mapping.set_rank(next_rank);

                        assert(!mappings_by_rank.count(mapping.rank()));
                        mappings_by_rank[mapping.rank()] = mapping;
                    }
                }

                for(auto& rank_and_mapping : mappings_by_rank) {
                    // Put in all the mappings. Ignore the rank since thay're already marked with and sorted by rank.
                    *path.add_mapping() = rank_and_mapping.second;
                }
                
                // Now the Path is rebuilt; stick it in the big output map.
                (*removed_paths)[path.name()] = path;
            }
            
#ifdef debug
            cerr << "Got all chunks; building XG index" << endl;
#endif
        }
    });
}

void VGset::for_each_kmer_parallel(int kmer_size, const function<void(const kmer_t&)>& lambda) {
    for_each([&lambda, kmer_size, this](VG* g) {
        g->show_progress = show_progress;
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

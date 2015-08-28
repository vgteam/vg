#include "vg_set.hpp"
#include "stream.hpp"

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

int64_t VGset::merge_id_space(void) {
    int64_t max_node_id = 0;
    int64_t max_path_id = 0;
    auto lambda = [&max_node_id, &max_path_id](VG* g) {
        if (max_node_id > 0) g->increment_node_ids(max_node_id);
        max_node_id = g->max_node_id();
    };
    transform(lambda);
    return max_node_id;
}

void VGset::to_xg(const string& xg_db_name) {
    // get a temporary graph
    string tmp_graph = xg_db_name + ".tmp_vg";
    ofstream os(tmp_graph);
    // concatenate the graphs we've been passed
    auto lambda = [&os](VG* g) {
        g->serialize_to_ostream(os);
    };
    for_each(lambda);
    os.close();
    // and load them into xg
    ifstream is(tmp_graph);
    xg::XG index;
    index.from_vg(is);
    is.close();
    remove(tmp_graph.c_str());
    // save the xg version to the file name we've been given
    ofstream db_out(xg_db_name.c_str());
    index.serialize(db_out);
    db_out.close();
}

void VGset::store_in_index(Index& index) {
    for_each([&index, this](VG* g) {
        g->show_progress = show_progress;
        index.load_graph(*g);
    });
}

void VGset::store_paths_in_index(Index& index) {
    for_each([&index, this](VG* g) {
        g->show_progress = show_progress;
        index.load_paths(*g);
    });
}

// stores kmers of size kmer_size with stride over paths in graphs in the index
void VGset::index_kmers(Index& index, int kmer_size, int edge_max, int stride, bool allow_negatives) {

    // create a vector of output files
    // as many as there are threads
    for_each([&index, kmer_size, edge_max, stride, allow_negatives, this](VG* g) {

        int thread_count;
#pragma omp parallel
        {
#pragma omp master
            thread_count = omp_get_num_threads();
        }

        // these are indexed by thread
        vector<vector<KmerMatch> > buffer;
        for (int i = 0; i < thread_count; ++i) {
            buffer.emplace_back();
        }
        // how many kmer entries to hold onto
        uint64_t buffer_max_size = 100000; // 100k

        // this may need a guard
        auto write_buffer = [&index](int tid, vector<KmerMatch>& buf) {
            rocksdb::WriteBatch batch;
            function<void(KmerMatch&)> keep_kmer = [&index, &batch](KmerMatch& k) {
                index.batch_kmer(k.sequence(), k.node_id(), k.position(), batch);
            };
            std::for_each(buf.begin(), buf.end(), keep_kmer);
            rocksdb::Status s = index.db->Write(rocksdb::WriteOptions(), &batch);
        };

        auto cache_kmer = [&buffer, &buffer_max_size, &write_buffer,
                           this](string& kmer, list<NodeTraversal>::iterator n, int p, list<NodeTraversal>& path, VG& graph) {
            if (allATGC(kmer)) {
                int tid = omp_get_thread_num();
                // note that we don't need to guard this
                // each thread has its own buffer!
                auto& buf = buffer[tid];
                KmerMatch k;
                k.set_sequence(kmer); k.set_node_id((*n).node->id()); k.set_position(p); k.set_backward((*n).backward);
                buf.push_back(k);
                if (buf.size() > buffer_max_size) {
                    write_buffer(tid, buf);
                    buf.clear();
                }
            }
        };

        g->create_progress("indexing kmers of " + g->name, buffer.size());
        g->for_each_kmer_parallel(kmer_size, edge_max, cache_kmer, stride, false, allow_negatives);
        g->destroy_progress();

        g->create_progress("flushing kmer buffers " + g->name, g->size());
        int tid = 0;
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < buffer.size(); ++i) {
            auto& buf = buffer[i];
            write_buffer(i, buf);
            g->update_progress(tid);
        }
        buffer.clear();
        g->destroy_progress();
    });

    index.remember_kmer_size(kmer_size);

}

void VGset::for_each_kmer_parallel(function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)>& lambda,
                                   int kmer_size, int edge_max, int stride, bool allow_dups, bool allow_negatives) {
    for_each([&lambda, kmer_size, edge_max, stride, allow_dups, allow_negatives, this](VG* g) {
        g->show_progress = show_progress;
        g->progress_message = "processing kmers of " + g->name;
        g->for_each_kmer_parallel(kmer_size, edge_max, lambda, stride, allow_dups, allow_negatives);
    });
}

void VGset::write_gcsa_out(ostream& out, int kmer_size, int edge_max, int stride, bool forward_only,
                           int64_t start_id, int64_t end_id) {

    // When we're sure we know what this kmer instance looks like, we'll write
    // it out exactly once. We need the start_end_id actually used in order to
    // go to the correct place when we don't go anywhere (i.e. at the far end of
    // the start/end node.
    auto write_kmer = [&start_id, &end_id](KmerPosition& kp){
        // We're going to write out every KmerPosition
        stringstream line;
        // Columns 1 and 2 are the kmer string and the node id:offset start position.
        line << kp.kmer << '\t' << kp.pos << '\t';
        // Column 3 is the comma-separated preceeding character options for this kmer instance.
        for (auto c : kp.prev_chars) line << c << ',';
        // If there are previous characters, kill the last comma. Otherwise, say "$" is the only previous character.
        if (!kp.prev_chars.empty()) { line.seekp(-1, line.cur);
        } else { line << '$'; }
        line << '\t';
        // Column 4 is the next character options from this kmer instance. Works just like column 3.
        for (auto c : kp.next_chars) line << c << ',';
        if (!kp.next_chars.empty()) { line.seekp(-1, line.cur);
        } else { line << '#'; }
        line << '\t';
        // Column 5 is the node id:offset positions of the places we can go
        // from here. They all start immediately after the last character of
        // this kmer.
        for (auto& p : kp.next_positions) line << p << ',';
        string rec = line.str();
        // handle origin marker
        // Go to the start/end node in forward orientation.
        if (kp.next_positions.empty()) { line << start_id << ":0"; rec = line.str(); }
        else { rec.pop_back(); }
#pragma omp critical (cout)
        {
            cout << rec << endl;
        }
    };

    // Run on each KmerPosition
    for_each_gcsa_kmer_position_parallel(kmer_size, edge_max, stride,
                                         forward_only,
                                         start_id, end_id,
                                         write_kmer);
    
}

void VGset::for_each_gcsa_kmer_position_parallel(int kmer_size, int edge_max, int stride,
                                                 bool forward_only,
                                                 int64_t& head_id, int64_t& tail_id,
                                                 function<void(KmerPosition&)> lambda) {

    // We have pointers to our single start/end node, and we will own it.
    // None of the VG graphs can own it since they get destroyed during the
    // for_each. TODO: the next free ID in the first graph (which creates this
    // node) must be free in all the graphs.
    Node* head_node = nullptr;
    Node* tail_node = nullptr;

    auto handle_node_in_graph = [&kmer_size, &edge_max, &stride, &forward_only, &head_node, &tail_node, &lambda, this]
        (VG* graph, Node* node) {
        // Go through all the kpaths of this node, and produce the GCSA2 kmers on both strands that start in this node.
#ifdef debug
        cerr << "Visiting node " << node->id() << endl;
#endif
        // This function runs in only one thread on a given node, so we can keep
        // our cache here. We gradually fill in each KmerPosition with all the
        // next positions and characters reachable with its string from its
        // orientation and offset along that strand in this node.
        map<tuple<string, bool, int32_t>, KmerPosition> cache;

        // We're going to visit every kmer of the node and run this:
        function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)>
            visit_kmer = [&cache, &kmer_size, &edge_max, &node, &forward_only, &head_node, &tail_node, this]
                         (string& kmer, list<NodeTraversal>::iterator start_node, int start_pos, 
                          list<NodeTraversal>& path, VG& graph) {
                         
            // We should never see negative offset kmers; _for_each_kmer ought to
            // have turned them around for positive offsets on the opposite strand.
            assert(start_pos >= 0);

            // todo, handle edge bounding
            // we need to check if the previous or next kmer will be excluded based on
            // edge bounding
            // if so, we should connect to the source or sink node
            
            // Get the information from the graph about what's before and after
            // this kmer, and where it ends.
            list<NodeTraversal>::iterator end_node;
            int32_t end_pos; // This counts in from the right of end_node.

            set<tuple<char, int64_t, bool, int32_t>> prev_positions;
            set<tuple<char, int64_t, bool, int32_t>> next_positions;

            // Fill in prev_chars, next_chars, prev_positions, and next_positions for the kmer by walking the path.
            graph.kmer_context(kmer,
                               kmer_size,
                               edge_max,
                               forward_only,
                               path,
                               start_node,
                               start_pos,
                               end_node,
                               end_pos,
                               prev_positions,
                               next_positions);

            if(start_node->node == node) {
                if (forward_only && start_node->backward) return;
                // This kmer starts on the node we're currently processing.
                // Store the information about it's forward orientation.
                
                // Get the KmerPosition to fill, creating it if it doesn't exist already.
                auto cache_key = make_tuple(kmer, start_node->backward, start_pos);
#ifdef debug
                if(cache.count(cache_key)) {
                    cerr << "F: Adding to " << kmer << " at " << start_node->node->id() << " " << start_node->backward
                         << " offset " << start_pos << endl;
                } else {
                    cerr << "F: Creating " << kmer << " at " << start_node->node->id() << " " << start_node->backward
                         << " offset " << start_pos << endl;
                }
#endif
                KmerPosition& forward_kmer = cache[cache_key];

                // fix up the context by swapping reverse-complements of the head and tail node
                // with their forward versions
                auto next_positions_copy = next_positions;
                next_positions.clear();
                for (auto p : next_positions_copy) {
                    char c = get<0>(p);
                    int64_t node_id = get<1>(p);
                    bool is_backward = get<2>(p);
                    int32_t pos = get<3>(p);
                    if (node_id == tail_node->id() && is_backward) {
                        // fixe the char as well...
                        c = head_node->sequence()[0]; // let's hope it has sequence
                        node_id = head_node->id();
                        is_backward = false;
                    } else if (node_id == head_node->id() && is_backward) {
                        c = tail_node->sequence()[0]; // let's hope it has sequence
                        node_id = tail_node->id();
                        is_backward = false;
                    }
                    next_positions.insert(make_tuple(c, node_id, is_backward, pos));
                }

                // Add in the kmer string
                if (forward_kmer.kmer.empty()) forward_kmer.kmer = kmer;
                
                // Add in the start position
                if (forward_kmer.pos.empty()) {
                    // And the distance from the end of the kmer to the end of its ending node.
                    stringstream ps;
                    if (start_node->node->id() == tail_node->id() && start_node->backward) {
                        ps << head_node->id() << ":" << start_pos;
                    } else if (start_node->node->id() == head_node->id() && start_node->backward) {
                        ps << tail_node->id() << ":" << start_pos;
                    } else {
                        ps << start_node->node->id() << ":" << (start_node->backward?"-":"") << start_pos;
                    }
                    forward_kmer.pos = ps.str();
                }
                
                // Add in the prev and next characters.
                for (auto& t : prev_positions) {
                    char c = get<0>(t);
                    forward_kmer.prev_chars.insert(c);
                }
                for (auto& t : next_positions) {
                    char c = get<0>(t);
                    forward_kmer.next_chars.insert(c);
                }
                
                // Add in the next positions
                for (auto& p : next_positions) {
                    // Figure out if the forward kmer should go next to the forward or reverse copy of the next node.
                    int64_t target_node = get<1>(p);
                    bool target_node_backward = get<2>(p);
                    int32_t target_off = get<3>(p);
                    // Say we go to it at the correct offset
                    stringstream ps;
                    ps << target_node << ":" << (target_node_backward?"-":"") << target_off;
                    forward_kmer.next_positions.insert(ps.str());
                }
            }

            if(end_node->node == node && !forward_only) {
                // This kmer ends on the node we're currently processing.
                // *OR* this kmer starts on the head node (in which case, we'd never see it from the other end)
                // Store the information about it's reverse orientation.

                // Get the KmerPosition to fill, creating it if it doesn't exist
                // already. We flip the backwardness because we look at the kmer
                // the other way, but since end_pos already counts from the end,
                // we don't touch it.
                auto cache_key = make_tuple(reverse_complement(kmer), !end_node->backward, end_pos);
#ifdef debug
                if(cache.count(cache_key)) {
                    cerr << "R: Adding to " << reverse_complement(kmer) << " at " << end_node->node->id() 
                         << " " << !(*end_node).backward << " offset " << end_pos << endl;
                } else {
                    cerr << "R: Creating " << reverse_complement(kmer) << " at " << end_node->node->id()
                         << " " << !(*end_node).backward << " offset " << end_pos << endl;
                }
#endif
                KmerPosition& reverse_kmer = cache[cache_key];

                // fix up the context by swapping reverse-complements of the head and tail node
                // with their forward versions
                auto prev_positions_copy = prev_positions;
                prev_positions.clear();
                for (auto p : prev_positions_copy) {
                    char c = get<0>(p);
                    int64_t node_id = get<1>(p);
                    bool is_backward = get<2>(p);
                    int32_t pos = get<3>(p);
                    if (node_id == tail_node->id() && !is_backward) {
                        node_id = head_node->id();
                        is_backward = true;
                        pos = graph.get_node(node_id)->sequence().size() - pos - 1;
                    } else if (node_id == head_node->id() && !is_backward) {
                        node_id = tail_node->id();
                        is_backward = true;
                        pos = tail_node->sequence().size() - pos - 1;
                    } else {
                        pos = graph.get_node(node_id)->sequence().size() - pos - 1;
                    }
                    prev_positions.insert(make_tuple(c, node_id, is_backward, pos));
                }

                // Add in the kmer string
                if (reverse_kmer.kmer.empty()) reverse_kmer.kmer = reverse_complement(kmer);
                
                // Add in the start position
                if (reverse_kmer.pos.empty()) {
                    // Use the other node ID, facing the other way
                    stringstream ps;
                    if (end_node->node->id() == tail_node->id() && !end_node->backward) {
                        ps << head_node->id() << ":" << end_pos;
                    } else if (end_node->node->id() == head_node->id() && !end_node->backward) {
                        ps << tail_node->id() << ":" << end_pos;
                    } else {
                        ps << (*end_node).node->id() << ":" << (!end_node->backward?"-":"") << end_pos;
                    }
                    // And the distance from the end of the kmer to the end of its ending node.
                    reverse_kmer.pos = ps.str();
                }
                    
                // Add in the prev and next characters.
                // fixme ... reverse complements things that should be translated to the head or tail node
                for (auto& t : prev_positions) {
                    char c = get<0>(t);
                    reverse_kmer.next_chars.insert(reverse_complement(c));
                }
                for (auto& t : next_positions) {
                    char c = get<0>(t);
                    reverse_kmer.prev_chars.insert(reverse_complement(c));
                }
                
                // Add in the next positions (using the prev positions since we're reversing)
                for (auto& p : prev_positions) {
                    // Figure out if the reverse kmer should go next to the forward or reverse copy of the next node.
                    int64_t target_node = get<1>(p);
                    bool target_node_backward = get<2>(p);
                    int32_t off = get<3>(p);

                    // Say we go to it at the correct offset
                    stringstream ps;
                    ps << target_node << ":" << (!target_node_backward?"-":"") << off;
                    reverse_kmer.next_positions.insert(ps.str());
                }
            }
        };
        
        // Now we visit every kmer of this node and fill in the cache. Don't
        // allow negative offsets; force them to be converted to positive
        // offsets on the reverse strand. But do allow different paths that
        // produce the same kmer, since GCSA2 needs those.
        graph->for_each_kmer_of_node(node, kmer_size, edge_max, visit_kmer, stride, true, false);
        
        // Now that the cache is full and correct, containing each kmer starting
        // on either strand of this node, send out all its entries.
        for(auto& kv : cache) {
            lambda(kv.second);
        }
        
    };

    // For every graph in our set (in serial), visit all the nodes in parallel and handle them.
    for_each([&handle_node_in_graph, kmer_size, edge_max, stride,
              &head_node, &tail_node, &head_id, &tail_id, this](VG* g) {
        g->show_progress = show_progress;
        g->progress_message = "processing kmers of " + g->name;
        
        if(head_node == nullptr) {
            assert(tail_node == nullptr); // they should be only set together
            // This is the first graph.
            // Add the start/end node, but make our own copy before we destroy the graph.
            g->add_start_end_markers(kmer_size, '#', '$', head_node, tail_node, head_id, tail_id);
            head_node = new Node(*head_node);
            tail_node = new Node(*tail_node);
            
            // Save its ID
            head_id = head_node->id();
            tail_id = tail_node->id();
        } else {
            // Add the existing start/end node
            int64_t maxid = g->max_node_id();
            if(head_node->id() <= maxid || tail_node->id() <= maxid) {
                // If the ID we got for the node when we made it in the
                // first graph is too small, we have to complain. It would be
                // nice if we could make a path through all the graphs, get the
                // max ID, and then use that to determine the new node ID.
                cerr << "error:[for_each_gcsa_kmer_position_parallel] created a start/end "
                     << "node in first graph with id used by later graph " << g->name 
                     << ". Put the graph with the largest node id first and try again." << endl;
                exit(1);
            }
            g->add_start_end_markers(kmer_size, '#', '$', head_node, tail_node, head_id, tail_id);
        }
        
        // Process all the kmers in the graph on a node-by-node basis. Make sure
        // to know about the graph when we do so (so we can get the kmers).
        g->for_each_node_parallel([&](Node* node) {
            handle_node_in_graph(g, node);
        });

        ofstream out("pre-gcsa.vg");
        g->serialize_to_ostream(out);
        out.close();
    });
    
    // delete the head and tail nodes
    if(head_node != nullptr) {
        delete head_node;
    }
    if(tail_node != nullptr) {
        delete tail_node;
    }
}

void VGset::get_gcsa_kmers(int kmer_size, int edge_max, int stride,
                           bool forward_only,
                           vector<gcsa::KMer>& kmers_out,
                           int64_t head_id, int64_t tail_id) {

    // TODO: This function goes through an internal string format that should
    // really be replaced by making some API changes to gcsa2.

    // We need an alphabet to parse the internal string format
    const gcsa::Alphabet alpha;
    
    // Each thread is going to make its own KMers, then we'll concatenate these all together at the end.
    vector<vector<gcsa::KMer>> thread_outputs;
    
#pragma omp parallel
    {
#pragma omp single
        {
            // Become parallel, get our number of threads, and make one of them make the per-thread outputs big enough.
            thread_outputs.resize(omp_get_num_threads());
        }
    }
    
    auto convert_kmer = [&thread_outputs, &alpha, &head_id, &tail_id](KmerPosition& kp) {
        // Convert this KmerPosition to several gcsa::Kmers, and save them in thread_outputs
                               
        // We need to make this kmer into a series of tokens
        vector<string> tokens;
        
        // First the kmer
        tokens.push_back(kp.kmer);
        
        // Then the node id:offset
        tokens.push_back(kp.pos);
        
        // Then the comma-separated preceeding characters. See <http://stackoverflow.com/a/18427254/402891>
        stringstream preceeding;
        copy(kp.prev_chars.begin(), kp.prev_chars.end(), ostream_iterator<char>(preceeding, ","));
        if(kp.prev_chars.empty()) {
            // If we don't have any previous characters, we come from "$"
            preceeding << "$";
        }
        tokens.push_back(preceeding.str());
        
        // And the comma-separated subsequent characters.
        stringstream subsequent;
        copy(kp.next_chars.begin(), kp.next_chars.end(), ostream_iterator<char>(subsequent, ","));
        if(kp.next_chars.empty()) {
            // If we don't have any next characters, we go to "#"
            subsequent << "#";
        }
        tokens.push_back(subsequent.str());
        
        // Finally, each of the node id:offset positions you can go to next (the successors).
        tokens.insert(tokens.end(), kp.next_positions.begin(), kp.next_positions.end());
        
        if (kp.next_positions.empty()) {
            // If we didn't have any successors, we have to say we go to the start of the start node
            tokens.push_back(to_string(tail_id) + ":0");
        }    
        
        for(size_t successor_index = 4; successor_index < tokens.size(); successor_index++) {
            // Now make a GCSA KMer for each of those successors, by passing the
            // tokens, the alphabet, and the index in the tokens of the
            // successor.
            
            thread_outputs[omp_get_thread_num()].emplace_back(tokens, alpha, successor_index);
            
            // Kmers that go to the sink/have stop characters still need to be marked as sorted.
            if(kp.kmer.rfind('$') != string::npos) {
                //(*(thread_outputs[omp_get_thread_num()].rbegin())).makeSorted();
#pragma omp critical
                {
                   // cout << "Marked " << *(thread_outputs[omp_get_thread_num()].rbegin()) << " as sorted early" << endl;
                }
            }
        }
        
    };
    
    // Run on each KmerPosition. This populates start_end_id, if it was 0, before calling convert_kmer.
    for_each_gcsa_kmer_position_parallel(kmer_size, edge_max, stride,
                                         forward_only,
                                         head_id, tail_id, convert_kmer);
                                         
    
    for(auto& thread_output : thread_outputs) {
        // Now throw everything into the output vector
        kmers_out.insert(kmers_out.end(), make_move_iterator(thread_output.begin()), make_move_iterator(thread_output.end()));
    }
    
}

}

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
                           this](string& kmer, Node* n, int p, list<Node*>& path, VG& graph) {
            if (allATGC(kmer)) {
                int tid = omp_get_thread_num();
                // note that we don't need to guard this
                // each thread has its own buffer!
                auto& buf = buffer[tid];
                KmerMatch k;
                k.set_sequence(kmer); k.set_node_id(n->id()); k.set_position(p);
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

void VGset::for_each_kmer_parallel(function<void(string&, Node*, int, list<Node*>&, VG&)>& lambda,
                                   int kmer_size, int edge_max, int stride, bool allow_dups, bool allow_negatives) {
    for_each([&lambda, kmer_size, edge_max, stride, allow_dups, allow_negatives, this](VG* g) {
        g->show_progress = show_progress;
        g->progress_message = "processing kmers of " + g->name;
        g->for_each_kmer_parallel(kmer_size, edge_max, lambda, stride, allow_dups, allow_negatives);
    });
}

void VGset::write_gcsa_out(ostream& out, int kmer_size, int edge_max, int stride, bool allow_dups) {

    // When we're sure we know what this kmer instance looks like, we'll write it out exactly once.
    auto write_kmer = [](KmerPosition& kp){
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
        if (kp.next_positions.empty()) { line << "0:0"; rec = line.str(); }
        else { rec.pop_back(); }
#pragma omp critical (cout)
        {
            cout << rec << endl;
        }
    };
    
    // Run on each KmerPosition
    for_each_gcsa_kmer_position_parallel(kmer_size, edge_max, stride, write_kmer, allow_dups);
    
}

void VGset::for_each_gcsa_kmer_position_parallel(int kmer_size, int edge_max, int stride, 
                                                 function<void(KmerPosition&)> lambda, bool allow_dups) {

    // This maps from thread number to current node ID, and a map from kmer
    // string and offset in the node to the KmerPosition describing that kmer
    // instance. Each node belongs to only one thread, and threads handle nodes
    // one at a time, so threads don't care about each others' caches or
    // information for other nodes.
    map<int, pair<int64_t, map<pair<string, int32_t>, KmerPosition>>> thread_caches;
#pragma omp parallel
    {
#pragma omp single
        for (int i = 0; i < omp_get_num_threads(); ++i) {
            // And we need to make sure that every thread has an empty set #0
            thread_caches[i].first = 0;
        }
    }

    // When we're done filling in our caches, we will process them with this function.
    auto process_cache = [&lambda](map<pair<string, int32_t>, KmerPosition >& cache){
        for (auto& k : cache) {
            // We're going to process every KmerPosition
            auto& kp = k.second;
            lambda(kp);
        }
    };

    // We're going to visit every kmer in each graph.
    function<void(string&, Node*, int, list<Node*>&, VG&)>
        visit_kmer = [&process_cache, &thread_caches, &kmer_size, &edge_max]
                     (string& kmer, Node* node, int pos, list<Node*>& path, VG& graph) {
        if (pos >= 0) {
            //kmer, starting position = (node id, offset), previous characters, successive characters, successive positions
            // todo, handle edge bounding
            // we need to check if the previous or next kmer will be excluded based on
            // edge bounding
            // if so, we should connect to the source or sink node
            set<char> prev_chars;
            set<char> next_chars;
            set<pair<int64_t, int32_t> > next_positions;
            // Fill in prev_chars, next_chars, and next_positions for the kmer by walking the path.
            graph.kmer_context(kmer,
                               kmer_size,
                               edge_max,
                               path,
                               node,
                               pos,
                               prev_chars,
                               next_chars,
                               next_positions);
            // We should put the kmer in our cache
            auto& cache = thread_caches[omp_get_thread_num()];
            if (cache.first != node->id()) {
                // If the last thing we put in the cache was from a different
                // node, process all those KmerPositions and start on the ones
                // for this node.
                process_cache(cache.second);
                cache.second.clear();
                cache.first = node->id();
            }
            // Is there something in the cache already for this kmer at this position?
            // If not we will get a default-constructed empty record and we can fill it in.
            auto& other = cache.second[make_pair(kmer, pos)];
            if (other.kmer.empty()) other.kmer = kmer;
            if (other.pos.empty()) {
                stringstream ps; ps << node->id() << ":" << pos;
                other.pos = ps.str();
            }
            // Add more previous and next characters and positions to the kmer.
            // This should actually find all of them the first time, no matter
            // what path we come by, since we just ask the graph about its
            // topology in kmer_context.
            for (auto c : prev_chars) other.prev_chars.insert(c);
            for (auto c : next_chars) other.next_chars.insert(c);
            for (auto p : next_positions) {
                stringstream ps; ps << p.first << ":" << p.second;
                other.next_positions.insert(ps.str());
            }
        }
    };

    // We have pointers to our own head and tail nodes, and we will own them.
    // None of the VG graphs can own them since they get destroyed during the
    // for_each. TODO: the next free ID in the first graph (which creates these
    // nodes) must be free in all the graphs.
    Node* head_node = nullptr;
    Node* tail_node = nullptr;

    // For every graph in our set, visit all the kmers in parallel and make and process KmerPositions for them.
    for_each([&visit_kmer, kmer_size, edge_max, stride, allow_dups, &head_node, &tail_node, this](VG* g) {
        g->show_progress = show_progress;
        g->progress_message = "processing kmers of " + g->name;
        
        if(head_node == nullptr) {
            // This is the first graph. Add the head and tail nodes, but make our own copies before we destroy the graph.
            g->add_start_and_end_markers(kmer_size, '#', '$', head_node, tail_node);
            head_node = new Node(*head_node);
            tail_node = new Node(*tail_node);
        } else {
            // Add the existing head and tail nodes
            if(min(head_node->id(), tail_node->id()) <= g->max_node_id()) {
                // If the IDs we got for these nodes when we made them in the
                // first graph are too small, we have to complain. It would be
                // nice if we could make a path through all the graphs, get the
                // max ID, and then use that to determine the head and tail node
                // IDs, but we can't because we ahve to be able to stream the
                // graphs.
                cerr << "error:[for_each_gcsa_kmer_position_parallel] created a head or "
                     << "tail node in first graph with id used by later graph " << g->name 
                     << ". Put the graph with the largest ids first and try again." << endl;
                exit(1);
            }
            g->add_start_and_end_markers(kmer_size, '#', '$', head_node, tail_node);
        }
        // Process all the kmers in the graph
        g->for_each_kmer_parallel(kmer_size, edge_max, visit_kmer, stride, allow_dups);
    });
    
// clean up caches
#pragma omp parallel
    {
        auto& cache = thread_caches[omp_get_thread_num()];
        process_cache(cache.second);
    }
    
    // delete the head and tail nodes
    if(head_node != nullptr) {
        delete head_node;
    }
    if(tail_node != nullptr) {
        delete tail_node;
    }    
}

void VGset::get_gcsa_kmers(int kmer_size, int edge_max, int stride, vector<gcsa::KMer>& kmers_out, bool allow_dups) {

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
    
    auto convert_kmer = [&thread_outputs, &alpha](KmerPosition& kp) {
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
            // If we didn't have any successors, we have to say we go to 0:0.
            tokens.push_back("0:0");
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
    
    // Run on each KmerPosition
    for_each_gcsa_kmer_position_parallel(kmer_size, edge_max, stride, convert_kmer, allow_dups);
    
    for(auto& thread_output : thread_outputs) {
        // Now throw everything into the output vector
        kmers_out.insert(kmers_out.end(), make_move_iterator(thread_output.begin()), make_move_iterator(thread_output.end()));
    }
    
}

}

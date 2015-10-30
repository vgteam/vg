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
    {
        stringstream cmd;
        cmd << "cat ";
        for (auto& name : filenames) {
            cmd << name << " ";
        }
        cmd << ">" << tmp_graph;
        if (system(cmd.str().c_str())) {
            cerr << "[vg::map] could not concatenate graphs" << endl;
            exit(1);
        }
    }
    // and load them into xg
    ifstream is(tmp_graph);
    xg::XG index;
    index.from_stream(is);
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

void VGset::for_each_kmer_parallel(
    const function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG&)>& lambda,
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

    // For every graph in our set (in serial), visit all the nodes in parallel and handle them.
    for_each([kmer_size, edge_max,stride, forward_only,
              &head_id, &tail_id, lambda](VG* g) {
                 g->for_each_gcsa_kmer_position_parallel(kmer_size, edge_max, stride,
                                                         forward_only,
                                                         head_id, tail_id,
                                                         lambda);
             });
}

void VGset::get_gcsa_kmers(int kmer_size, int edge_max, int stride,
                           bool forward_only,
                           vector<gcsa::KMer>& kmers_out,
                           int64_t head_id, int64_t tail_id) {
    for_each([kmer_size, edge_max,stride, forward_only, &kmers_out,
              &head_id, &tail_id](VG* g) {
                 g->get_gcsa_kmers(kmer_size, edge_max, stride,
                                   forward_only,
                                   kmers_out,
                                   head_id, tail_id);
             });
}

}

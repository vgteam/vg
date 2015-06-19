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

    struct KmerPosition {
        string kmer;
        string pos;
        set<char> prev_chars;
        set<char> next_chars;
        set<string> next_positions;
    };

    map<int, pair<int64_t, map<pair<string, int32_t>, KmerPosition > > > output_cache;
#pragma omp parallel
    {
#pragma omp single
        for (int i = 0; i < omp_get_num_threads(); ++i) {
            output_cache[i].first = 0;
        }
    }

    auto write_cache = [](map<pair<string, int32_t>, KmerPosition >& cache){
        for (auto& k : cache) {
            auto& kp = k.second;
            stringstream line;
            line << kp.kmer << '\t' << kp.pos << '\t';
            for (auto c : kp.prev_chars) line << c << ',';
            if (!kp.prev_chars.empty()) { line.seekp(-1, line.cur);
            } else { line << '$'; }
            line << '\t';
            for (auto c : kp.next_chars) line << c << ',';
            if (!kp.next_chars.empty()) { line.seekp(-1, line.cur);
            } else { line << '#'; }
            line << '\t';
            for (auto& p : kp.next_positions) line << p << ',';
            string rec = line.str();
            // handle origin marker
            if (kp.next_positions.empty()) { line << "0:0"; rec = line.str(); }
            else { rec.pop_back(); }
#pragma omp critical (cout)
            {
                cout << rec << endl;
            }
        }
    };

    Node* head_node=NULL;
    Node* tail_node=NULL;

    function<void(string&, Node*, int, list<Node*>&, VG&)>
        lambda = [&write_cache, &output_cache, &kmer_size, &edge_max, &head_node, &tail_node]
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
            graph.kmer_context(kmer,
                               kmer_size,
                               edge_max,
                               path,
                               node,
                               pos,
                               prev_chars,
                               next_chars,
                               next_positions);

            auto& cache = output_cache[omp_get_thread_num()];
            if (cache.first != node->id()) {
                write_cache(cache.second);
                cache.second.clear();
                cache.first = node->id();
            }
            auto& other = cache.second[make_pair(kmer, pos)];
            if (other.kmer.empty()) other.kmer = kmer;
            if (other.pos.empty()) {
                stringstream ps; ps << node->id() << ":" << pos;
                other.pos = ps.str();
            }
            for (auto c : prev_chars) other.prev_chars.insert(c);
            for (auto c : next_chars) other.next_chars.insert(c);
            for (auto p : next_positions) {
                stringstream ps; ps << p.first << ":" << p.second;
                other.next_positions.insert(ps.str());
            }
        }
    };

    for_each([&lambda, kmer_size, edge_max, stride, allow_dups, &head_node, &tail_node, this](VG* g) {
        g->show_progress = show_progress;
        g->progress_message = "processing kmers of " + g->name;
        // add in start and end markers that are required by GCSA
        g->add_start_and_end_markers(kmer_size, '#', '$', head_node, tail_node);
        g->for_each_kmer_parallel(kmer_size, edge_max, lambda, stride, allow_dups);
    });

// clean up caches
#pragma omp parallel
    {
        auto& cache = output_cache[omp_get_thread_num()];
        write_cache(cache.second);
    }
}

}

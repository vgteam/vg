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

// stores kmers of size kmer_size with stride over paths in graphs in the index
void VGset::index_kmers(Index& index, int kmer_size, int edge_max, int stride) {

    // create a vector of output files
    // as many as there are threads
    for_each([&index, kmer_size, edge_max, stride, this](VG* g) {

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

        auto write_buffer = [&index](int tid, vector<KmerMatch>& buf) {
            rocksdb::WriteBatch batch;
            function<void(KmerMatch&)> keep_kmer = [&index, &batch](KmerMatch& k) {
                index.batch_kmer(k.sequence(), k.node_id(), k.position(), batch);
            };
            std::for_each(buf.begin(), buf.end(), keep_kmer);
            rocksdb::Status s = index.db->Write(rocksdb::WriteOptions(), &batch);
        };

        auto cache_kmer = [&buffer, &buffer_max_size, &write_buffer,
                           this](string& kmer, Node* n, int p) {
            if (allATGC(kmer)) {
                int tid = omp_get_thread_num();
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
        g->for_each_kmer_parallel(kmer_size, edge_max, cache_kmer, stride);
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

void VGset::for_each_kmer_parallel(function<void(string&, Node*, int)>& lambda,
                                   int kmer_size, int edge_max, int stride) {
    for_each([&lambda, kmer_size, edge_max, stride, this](VG* g) {
        g->show_progress = show_progress;
        g->progress_message = "processing kmers of " + g->name;
        g->for_each_kmer_parallel(kmer_size, edge_max, lambda, stride);
    });
}


}

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
    int64_t max_id = 0;
    auto lambda = [&max_id](VG* g) {
        if (max_id > 0) g->increment_node_ids(max_id);
        max_id = g->max_node_id();
    };
    transform(lambda);
    return max_id;
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
        vector<int> written_buffer_count;
        for (int i = 0; i < thread_count; ++i) {
            buffer.emplace_back();
            written_buffer_count.push_back(0);
        }
        // how many kmer entries to hold onto
        uint64_t buffer_max_size = 1000000;

        auto write_buffer = [&written_buffer_count, &index](int tid, vector<KmerMatch>& buf) {
            stringstream file_name;
            file_name << index.name << "/" << tid << "." << written_buffer_count[tid]++;
            ofstream out(file_name.str());
            function<KmerMatch(uint64_t)> write_kmer =
                [&buf](uint64_t i) {
                return buf.at(i);
            };
            stream::write(out, buf.size(), write_kmer);
            out.close();
        };

        auto cache_kmer = [&buffer, &written_buffer_count,
                           &buffer_max_size, &write_buffer,
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

        g->create_progress("caching kmers of " + g->name, buffer.size());
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

        int total_buffers = 0;
        for (int count : written_buffer_count) {
            total_buffers += count;
        }

        g->create_progress("indexing kmers " + g->name, total_buffers);
        int written_buffers = 0;
//#pragma omp parallel for schedule(dynamic)
        for (int tid = 0; tid < written_buffer_count.size(); ++tid) {
            int count = written_buffer_count[tid];
            for (int i = 0; i < count; ++ i) {
                stringstream file_name;
                file_name << index.name << "/" << tid << "." << i;
                ifstream in(file_name.str());
                function<void(KmerMatch&)> keep_kmer = [&index, this](KmerMatch& k) {
                    index.put_kmer(k.sequence(), k.node_id(), k.position());
                };
                stream::for_each(in, keep_kmer);
                int r = system((string("rm ") + file_name.str()).c_str());
                g->update_progress(++written_buffers);
            }
        }
        index.remember_kmer_size(kmer_size);
        g->destroy_progress();

    });

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

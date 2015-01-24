#include "vg_set.hpp"

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
    index.close();
    index.prepare_for_bulk_load();
    index.open();

    for_each([&index, this](VG* g) {
        g->show_progress = show_progress;
        index.load_graph(*g);
    });
    // clean up after bulk load
    index.flush();
    index.close();
    index.reset_options();
    index.open();
    index.compact();
}

// stores kmers of size kmer_size with stride over paths in graphs in the index
void VGset::index_kmers(Index& index, int kmer_size, int stride) {

    index.close();
    index.prepare_for_bulk_load();
    index.open();

    for_each([&index, kmer_size, stride, this](VG* g) {

        // set up an index per process
        // we merge them at the end of this graph
        vector<Index*> indexes;
        vector<int64_t> counts;
        int thread_count = 1;
#pragma omp parallel
        {
#pragma omp master
            {
                thread_count = omp_get_num_threads();
                for (int i = 0; i < thread_count; ++i) {
                    stringstream s;
                    s << index.name << "." << i;
                    string n = s.str();
                    Index* idx = new Index;
                    idx->prepare_for_bulk_load();
                    idx->open(n);
                    indexes.push_back(idx);
                    counts.push_back(0);
                }
            }
        }

        auto keep_kmer = [&indexes, &counts, this](string& kmer, Node* n, int p) {
            if (allATGC(kmer)) {
                counts[omp_get_thread_num()]++;
                indexes[omp_get_thread_num()]->put_kmer(kmer, n->id(), p);
                //index->put_kmer(kmer, n->id(), p);
            }
        };

        g->create_progress("indexing kmers of " + g->name, g->size());
        g->for_each_kmer_parallel(kmer_size, keep_kmer, stride);
        g->destroy_progress();

        int64_t total_kmers = 0;
        for (auto i : counts) total_kmers += i;
        int64_t count = 0;

        // merge results
        index.remember_kmer_size(kmer_size);
        g->create_progress("merging kmers of " + g->name, total_kmers);

//#pragma omp parallel for schedule(static, 1)
        for (int i = 0; i < thread_count; ++i) {
            auto* idx = indexes[i];
            idx->for_all([g, &index, &count](string& k, string& v) {
#pragma omp atomic
                ++count;
                if (count % 10000) g->update_progress(count);
                index.db->Put(index.write_options, k, v);
            });
            string dbname = idx->name;
            delete idx;
            int r = system((string("rm -r ") + dbname).c_str());
        }
        g->destroy_progress();
    });

    // clean up after bulk load
    index.flush();
    index.close();
    index.reset_options();
    index.open();
    index.compact();
}

void VGset::for_each_kmer_parallel(function<void(string&, Node*, int)>& lambda,
                                   int kmer_size, int stride) {
    for_each([&lambda, kmer_size, stride, this](VG* g) {
        g->show_progress = show_progress;
        g->progress_message = "processing kmers of " + g->name;
        g->for_each_kmer_parallel(kmer_size, lambda, stride);
    });
}


}

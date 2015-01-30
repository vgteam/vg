#include "mapper.hpp"

namespace vg {

Mapper::Mapper(Index* idex)
    : index(idex)
    , best_clusters(0) {
    kmer_sizes = index->stored_kmer_sizes();
    if (kmer_sizes.empty()) {
        cerr << "error:[vg::Mapper] the index (" 
             << index->name << ") does not include kmers" << endl;
        exit(1);
    }
}

Mapper::~Mapper(void) {
    // noop
}

Alignment Mapper::align(string& sequence, int stride) {
    Alignment alignment;
    alignment.set_sequence(sequence);
    return align(alignment, stride);
}

Alignment& Mapper::align(Alignment& alignment, int stride) {

    if (index == NULL) {
        cerr << "error:[vg::Mapper] no index loaded, cannot map alignment!" << endl;
        exit(1);
    }

    // establish kmers

    const string& sequence = alignment.sequence();
    auto kmers = kmers_of(sequence, stride);

    vector<map<int64_t, set<int32_t> > > positions(kmers.size());
    int i = 0;
    for (auto& k : kmers) {
        index->get_kmer_positions(k, positions.at(i++));
    }

    // make threads
    // these start whenever we have a kmer match which is outside of
    // one of the last positions (for the previous kmer) + the kmer stride % wobble (hmm)

    map<pair<int64_t, int32_t>, vector<int64_t> > position_threads;
    map<int64_t, vector<int64_t> > node_threads;
    int node_wobble = 2;
    int position_wobble = 2;
    int kmer_size = *kmer_sizes.begin(); // just use the first for now

    i = 0;
    for (auto& p : positions) {
        auto& kmer = kmers.at(i++);
        for (auto& x : p) {
            int64_t id = x.first;
            set<int32_t>& pos = x.second;
            for (auto& y : pos) {
                //cout << kmer << "\t" << i << "\t" << id << "\t" << y << endl;
                // thread rules
                // if we find the previous position
                int m = 0;
                vector<int64_t> thread;
                for (int j = 0; j < 2*position_wobble + 1; ++j) {
                    if (j == 0) { // on point
                    } else if (j % 2 == 0) { // subtract
                        m *= -1;
                    } else { // add
                        m *= -1; ++m;
                    }
                    //cout << "checking " << id << " " << y << " - " << kmer_size << " + " << m << endl;
                    if (position_threads.find(make_pair(id, y - stride + m)) != position_threads.end()) {
                        //length = position_threads[make_pair(id, y - stride + m)] + 1;
                        thread = position_threads[make_pair(id, y - stride + m)];
                        //cout << "thread is " << thread.size() << " long" << endl;
                        break;
                    }
                }
                // if length == 1, maybe we should look at a neighboring (previous) node
                // so we wobble around looking for one
                for (int j = 1; thread.empty() && j < node_wobble; ++j) {
                    //cout << "checking " << id << " " << y << " - " << kmer_size << " + " << m << endl;
                    if (node_threads.find(id - j) != node_threads.end()) {
                        thread = node_threads[id - j];
                        // isn't this guaranteed?
                        //cout << "thread is " << thread.size() << " long" << endl;
                        break;
                    }
                }

                thread.push_back(id);
                position_threads[make_pair(id, y)] = thread;
                node_threads[id] = thread;

            }
        }
    }

    map<int, vector<vector<int64_t> > > threads_by_length;
    for (auto& t : node_threads) {
        auto& thread = t.second;
        auto& threads = threads_by_length[thread.size()];
        threads.push_back(thread);
    }

    // for debugging
    /*
    for (auto& t : threads_by_length) {
        auto& length = t.first;
        auto& threads = t.second;
        cerr << length << ":" << endl;
        for (auto& thread : threads) {
            cerr << "\t";
            for (auto& id : thread) {
                cerr << id << " ";
            }
            cerr << endl;
        }
        cerr << endl;
    }
    */

    // collect the nodes from the best N threads by length
    // and expand subgraphs as before
    set<int64_t> nodes;
    VG* graph = new VG;
    map<int, vector<vector<int64_t> > >::reverse_iterator tl = threads_by_length.rbegin();
    for (int i = 0; tl != threads_by_length.rend() && (best_clusters == 0 || i < best_clusters); ++i, ++tl) {
        VG g;
        auto& threads = tl->second;
        for (auto& thread : threads) {
            for (auto& id : thread) {
                if (!nodes.count(id)) {
                    nodes.insert(id);
                    index->get_context(id, g);
                    index->get_connected_nodes(g);
                    graph->extend(g);
                }
            }
        }
    }

    /*
    ofstream f("vg_align.vg");
    graph->serialize_to_ostream(f);
    f.close();
    */

    int max_iter = 10;
    int iter = 0;
    int context_step = 1;
    int64_t max_subgraph_size = 0;

    auto get_max_subgraph_size = [this, &max_subgraph_size, &graph]() {
        list<VG> subgraphs;
        graph->disjoint_subgraphs(subgraphs);
        for (auto& subgraph : subgraphs) {
            max_subgraph_size = max(subgraph.total_length_of_nodes(), max_subgraph_size);
        }
    };

    get_max_subgraph_size();

    while (max_subgraph_size < sequence.size()*2 && iter < max_iter) {
        index->expand_context(*graph, context_step);
        get_max_subgraph_size();
        ++iter;
    }

    return graph->align(alignment);

}

vector<string> Mapper::kmers_of(const string& seq, const int stride) {
    vector<string> kmers;
    if (!seq.empty()) {
        for (int kmer_size : kmer_sizes) {
            for (int i = 0; i < seq.size()-kmer_size; i+=stride) {
                kmers.push_back(seq.substr(i,kmer_size));
            }
        }
    }
    return kmers;
}

}

#include "mapper.hpp"

namespace vg {

Mapper::Mapper(Index* idex)
    : index(idex)
    , best_clusters(0)
    , hit_max(100) {
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
    return align_threaded(alignment, stride);
}

Alignment& Mapper::align_threaded(Alignment& alignment, int stride) {

    // parameters, some of which should probably be modifiable
    // TODO -- move to Mapper object

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
        index->get_kmer_positions(k, positions.at(i));
        if (positions.at(i).size() > hit_max) positions.at(i).clear();
        ++i;
    }

    // make threads
    // these start whenever we have a kmer match which is outside of
    // one of the last positions (for the previous kmer) + the kmer stride % wobble (hmm)

    map<int64_t, vector<int> > node_kmer_order;
    map<pair<int64_t, int32_t>, vector<int64_t> > position_threads;
    map<int64_t, vector<int64_t> > node_threads;
    //int node_wobble = 0; // turned off...
    int position_wobble = 2;
    int kmer_size = *kmer_sizes.begin(); // for simplicity, use the first available kmer size

    int max_iter = sequence.size();
    int iter = 0;
    int context_step = 1;
    int64_t max_subgraph_size = 0;
    int max_thread_gap = 30; // counted in nodes

    i = 0;
    for (auto& p : positions) {
        auto& kmer = kmers.at(i++);
        for (auto& x : p) {
            int64_t id = x.first;
            set<int32_t>& pos = x.second;
            node_kmer_order[id].push_back(i-1);
            for (auto& y : pos) {
                //cerr << kmer << "\t" << i << "\t" << id << "\t" << y << endl;
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
                    //cerr << "checking " << id << " " << y << " - " << kmer_size << " + " << m << endl;
                    auto previous = position_threads.find(make_pair(id, y - stride + m));
                    if (previous != position_threads.end()) {
                        //length = position_threads[make_pair(id, y - stride + m)] + 1;
                        thread = previous->second;
                        position_threads.erase(previous);
                        //cerr << "thread is " << thread.size() << " long" << endl;
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

    // now sort the threads and re-cluster them

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

    // sort threads by ids
    // if they are at least 2 hits long
    set<vector<int64_t> > sorted_threads;
    auto tl = threads_by_length.rbegin();
    for (auto& t : node_threads) {
        auto& thread = t.second;
        // if you use very short kmers, this may become necessary
        //if (thread.size() > 1) {
        sorted_threads.insert(thread);
        //}
    }
    // clean up
    threads_by_length.clear();

    // go back through and combine closely-linked threads
    // ... but only if their kmer order is proper
    map<int64_t, vector<int64_t> > threads_by_last;
    for (auto& thread : sorted_threads) {
        //cerr << thread.front() << "-" << thread.back() << endl;
        auto prev = threads_by_last.upper_bound(thread.front()-max_thread_gap);
        //if (prev != threads_by_last.begin()) --prev;
        // now we should be at the highest thread within the bounds
        //cerr << prev->first << " " << thread.front() << endl;
        // todo: it may also make sense to check that the kmer order makes sense
        // what does this mean? it means that the previous 
        if (prev != threads_by_last.end()
            && prev->first > thread.front() - max_thread_gap) {
            vector<int64_t> new_thread;
            auto& prev_thread = prev->second;
            new_thread.reserve(prev_thread.size() + thread.size());
            new_thread.insert(new_thread.end(), prev_thread.begin(), prev_thread.end());
            new_thread.insert(new_thread.end(), thread.begin(), thread.end());
            threads_by_last.erase(prev);
            // this will clobber... not good
            // maybe overwrite only if longer?
            threads_by_last[new_thread.back()] = new_thread;
        } else {
            threads_by_last[thread.back()] = thread;
        }
    }

    // debugging
    /*
    for (auto& t : threads_by_last) {
        auto& thread = t.second;
        cerr << t.first << "\t";
        for (auto& id : thread) {
            cerr << id << " ";
        }
        cerr << endl;
    }
    */

    // rebuild our threads_by_length set
    for (auto& t : threads_by_last) {
        auto& thread = t.second;
        auto& threads = threads_by_length[thread.size()];
        threads.push_back(thread);
    }

    // collect the nodes from the best N threads by length
    // and expand subgraphs as before
    VG* graph = new VG;
    tl = threads_by_length.rbegin();
    for (int i = 0; tl != threads_by_length.rend() && (best_clusters == 0 || i < best_clusters); ++i, ++tl) {
        auto& threads = tl->second;
        for (auto& thread : threads) {
            VG g;
            set<int64_t> nodes;
            for (auto& id : thread) {
                if (!nodes.count(id)) {
                    nodes.insert(id);
                    index->get_context(id, g);
                }
            }
            // by definition, our thread should construct a contiguous graph
            // however, it might not due to complexities of thread determination
            // so, enforce this here
            list<VG> subgraphs;
            int c = 0;
            do {
                subgraphs.clear();
                index->expand_context(g, 1);
                // to find the disjoint subgraphs we need a complete graph (no orphan edges)
                // so we make a temporary copy of the graph
                VG gprime = g;
                index->get_connected_nodes(gprime);
                gprime.disjoint_subgraphs(subgraphs);
            } while(subgraphs.size() > 1 && c++ < max_thread_gap);
            // now that we have a completed graph for our thread, save it
            graph->extend(g);
        }
    }

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
        index->get_connected_nodes(*graph);
        get_max_subgraph_size();
        ++iter;
    }
    // ensure we have a complete graph prior to alignment
    index->get_connected_nodes(*graph);

    // align, check for soft clipping, and expand context against the target if needed

    graph->align(alignment);

    // check if we start or end with soft clips
    // if so, try to expand the graph until we don't have any more (or we hit a threshold)

    int sc_start = softclip_start(alignment);
    int sc_end = softclip_end(alignment);
    //cerr << alignment.score() << " " << sc_start << " " << sc_end << endl;

    // NB it will probably be faster here to not query the DB again,
    // but instead to grow the matching graph in the right direction
    if (sc_start > 0 || sc_end > 0) {
        // get the target graph
        delete graph;
        graph = new VG;
        // nodes
        Path* path = alignment.mutable_path();
        for (int i = 0; i < path->mapping_size(); ++i) {
            index->get_context(path->mutable_mapping(i)->node_id(), *graph);
        }
        while (graph->total_length_of_nodes() < sequence.size() * 2) {
            index->expand_context(*graph, context_step); // expand faster here
        }
        index->get_connected_nodes(*graph);

        /*
        cerr << "graph " << graph->size() << endl;
        ofstream f("vg_align.vg");
        graph->serialize_to_ostream(f);
        f.close();
        */

        alignment.clear_path();
        graph->align(alignment);

    }


    return alignment;

}

int Mapper::softclip_start(Alignment& alignment) {
    if (alignment.mutable_path()->mapping_size() > 0) {
        Path* path = alignment.mutable_path();
        Mapping* first_mapping = path->mutable_mapping(0);
        Edit* first_edit = first_mapping->mutable_edit(0);
        if (first_edit->type() == Edit_Type_SOFTCLIP) {
            return first_edit->length();
        }
    }
    return 0;
}

int Mapper::softclip_end(Alignment& alignment) {
    if (alignment.mutable_path()->mapping_size() > 0) {
        Path* path = alignment.mutable_path();
        Mapping* last_mapping = path->mutable_mapping(path->mapping_size()-1);
        Edit* last_edit = last_mapping->mutable_edit(last_mapping->edit_size()-1);
        if (last_edit->type() == Edit_Type_SOFTCLIP) {
            return last_edit->length();
        }
    }
    return 0;
}

Alignment& Mapper::align_simple(Alignment& alignment, int stride) {

    if (index == NULL) {
        cerr << "error:[vg::Mapper] no index loaded, cannot map alignment!" << endl;
        exit(1);
    }

    // establish kmers

    const string& sequence = alignment.sequence();
    auto kmers = kmers_of(sequence, stride);

    map<string, int32_t> kmer_counts;
    vector<map<int64_t, set<int32_t> > > positions(kmers.size());
    int i = 0;
    for (auto& k : kmers) {
        index->get_kmer_positions(k, positions.at(i++));
        kmer_counts[k] = positions.at(i-1).size();
    }
    positions.clear();
    VG* graph = new VG;
    for (auto& c : kmer_counts) {
        if (c.second < hit_max) {
            index->get_kmer_subgraph(c.first, *graph);
        }
    }

    int max_iter = sequence.size();
    int iter = 0;
    int context_step = 1;
    int64_t max_subgraph_size = 0;

    // use kmers which are informative
    // and build up the graph

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
        index->get_connected_nodes(*graph);
        get_max_subgraph_size();
        ++iter;
    }
    // ensure we have a complete graph prior to alignment
    index->get_connected_nodes(*graph);

    /*
    ofstream f("vg_align.vg");
    graph->serialize_to_ostream(f);
    f.close();
    */

    graph->align(alignment);

    return alignment;

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

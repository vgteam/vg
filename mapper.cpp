#include "mapper.hpp"

namespace vg {

Mapper::Mapper(Index* idex)
    : index(idex)
    , best_clusters(0)
    , hit_max(100)
    , hit_size_threshold(0)
    , kmer_min(21)
    , thread_extension(10)
    , thread_extension_max(80)
{
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

Alignment Mapper::align(string& sequence, int kmer_size, int stride) {

    std::chrono::time_point<std::chrono::system_clock> start_both, end_both;
    start_both = std::chrono::system_clock::now();

    // forward
    Alignment alignment_f;
    alignment_f.set_sequence(sequence);

    {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        align_threaded(alignment_f, kmer_size, stride);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        cerr << elapsed_seconds.count() << "\t" << "+" << "\t" << alignment_f.sequence() << endl;
    }

    // reverse
    Alignment alignment_r;
    alignment_r.set_sequence(reverse_complement(sequence));

    {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        align_threaded(alignment_r, kmer_size, stride);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;

        cerr << elapsed_seconds.count() << "\t" << "-" << "\t" << alignment_r.sequence() << endl;
    }

    end_both = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_both = end_both-start_both;

    cerr << elapsed_seconds_both.count() << "\t" << "b" << "\t" << sequence << endl;

    if (alignment_r.score() > alignment_f.score()) {
        return alignment_r;
    } else {
        return alignment_f;
    }
}

Alignment& Mapper::align_threaded(Alignment& alignment, int kmer_size, int stride, int attempt, int hit_count) {

    // parameters, some of which should probably be modifiable
    // TODO -- move to Mapper object

    if (index == NULL) {
        cerr << "error:[vg::Mapper] no index loaded, cannot map alignment!" << endl;
        exit(1);
    }

    // establish kmers

    const string& sequence = alignment.sequence();

    // if kmer size is not specified, pick it up from the index
    // for simplicity, use the first available kmer size; this could change
    if (kmer_size == 0) kmer_size = *kmer_sizes.begin();
    // and start with stride such that we barely cover the read with kmers
    if (stride == 0) stride = sequence.size() / ceil((double)sequence.size() / kmer_size);

    int kmer_hit_count = 0;
    int kept_kmer_count = 0;

    auto align_with_increased_sensitivity = [this,
                                             &kmer_size,
                                             &stride,
                                             &sequence,
                                             &alignment,
                                             &attempt,
                                             &kmer_hit_count]() -> Alignment& {
        if ((double)stride/kmer_size < 0.3 && kmer_size -2 >= kmer_min) {
            kmer_size -= 2;
            stride = sequence.size() / ceil((double)sequence.size() / kmer_size);
            cerr << "realigning with " << kmer_size << " " << stride << endl;
            alignment.clear_path();
            return align_threaded(alignment, kmer_size, stride, ++attempt, kmer_hit_count);
        } else if ((double)stride/kmer_size >= 0.3 && kmer_size >= kmer_min) {
            stride = max(1, stride/2);
            cerr << "realigning with " << kmer_size << " " << stride << endl;
            alignment.clear_path();
            return align_threaded(alignment, kmer_size, stride, ++attempt, kmer_hit_count);
        }
    };

    auto kmers = balanced_kmers(sequence, kmer_size, stride);

    vector<map<int64_t, vector<int32_t> > > positions(kmers.size());
    int i = 0;
    for (auto& k : kmers) {
        cerr << k << "\t" << index->approx_size_of_kmer_matches(k) << endl;
        // if we have more than one block worth of kmers on disk, consider this kmer non-informative
        if (index->approx_size_of_kmer_matches(k) > hit_size_threshold) {
            continue;
        }
        index->get_kmer_positions(k, positions.at(i));
        kmer_hit_count += positions.at(i).size();
        if (positions.at(i).size() > hit_max) positions.at(i).clear();
        kept_kmer_count += positions.at(i).size();
        ++i;
    }

    cerr << "kmer hits " << kmer_hit_count << endl;
    cerr << "kept kmer hits " << kept_kmer_count << endl;

    //if (kept_kmer_count == 0 && attempt == 0) {
    //return align_with_increased_sensitivity();
//}
/*
    if (kmer_hit_count > 0 && kept_kmer_count == 0) { // || kept_kmer_count > 5* alignment.sequence().size()) {
        cerr << "bailout" << endl;
        //return alignment;
    }
*/


    // make threads
    // these start whenever we have a kmer match which is outside of
    // one of the last positions (for the previous kmer) + the kmer stride % wobble (hmm)

    map<int64_t, vector<int> > node_kmer_order;
    map<pair<int64_t, int32_t>, vector<int64_t> > position_threads;
    map<int64_t, vector<int64_t> > node_threads;
    //int node_wobble = 0; // turned off...
    int position_wobble = 2;

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
            vector<int32_t>& pos = x.second;
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
        if (thread.size() > 1) {
            sorted_threads.insert(thread);
        }
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

    int thread_ex = thread_extension;
    while (alignment.score() == 0 && thread_ex <= thread_extension_max) {
        // collect the nodes from the best N threads by length
        // and expand subgraphs as before
        VG* graph = new VG;
        tl = threads_by_length.rbegin();
        for (int i = 0; tl != threads_by_length.rend() && (best_clusters == 0 || i < best_clusters); ++i, ++tl) {
            auto& threads = tl->second;
            // by definition, our thread should construct a contiguous graph
            for (auto& thread : threads) {
                int64_t first = *thread.begin() - thread_ex;
                int64_t last = *thread.rbegin() + thread_ex;
                // so we can pick it up efficiently from the index by pulling the range from first to last
                index->get_range(first, last, *graph);
            }
        }

        // by default, expand the graph a bit so we are likely to map
        index->expand_context(*graph, 1);
        index->get_connected_nodes(*graph);

        // align
        alignment.clear_path();
        graph->align(alignment);
        delete graph;
        if (alignment.score() == 0) {
            thread_ex *= 2;
        }
    }

    // did we still fail to align?
    // if so, decrease the stride; if we are already at decreased stride, decrease the kmer size
    if (alignment.score() == 0 && (kmer_hit_count > 0 || hit_count > 0)) {
        return align_with_increased_sensitivity();
    }

    // check if we start or end with soft clips
    // if so, try to expand the graph until we don't have any more (or we hit a threshold)

    int sc_start = softclip_start(alignment);
    int sc_end = softclip_end(alignment);
    //cerr << alignment.score() << " " << sc_start << " " << sc_end << endl;

    // NB it will probably be faster here to not fully query the DB again,
    // but instead to grow the matching graph in the right direction
    if (sc_start > 0 || sc_end > 0) {
        // get the target graph
        VG* graph = new VG;
        // nodes
        Path* path = alignment.mutable_path();
        for (int i = 0; i < path->mapping_size(); ++i) {
            index->get_context(path->mutable_mapping(i)->node_id(), *graph);
        }
        while (graph->total_length_of_nodes() < sequence.size() * 3) {
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
        delete graph;

    }

    if (alignment.score() == 0) cerr << "failed alignment" << endl;

    return alignment;

}

int softclip_start(Alignment& alignment) {
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

int softclip_end(Alignment& alignment) {
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

string reverse_complement(const string& seq) {
    string rc;
    rc.assign(seq.rbegin(), seq.rend());
    for (auto& c : rc) {
        switch (c) {
        case 'A': c = 'T'; break;
        case 'T': c = 'A'; break;
        case 'G': c = 'C'; break;
        case 'C': c = 'G'; break;
        case 'N': c = 'N'; break;
        default: break;
        }
    }
    return rc;
}

Alignment& Mapper::align_simple(Alignment& alignment, int kmer_size, int stride) {

    if (index == NULL) {
        cerr << "error:[vg::Mapper] no index loaded, cannot map alignment!" << endl;
        exit(1);
    }

    // establish kmers
    const string& sequence = alignment.sequence();
    //  
    auto kmers = balanced_kmers(sequence, kmer_size, stride);

    map<string, int32_t> kmer_counts;
    vector<map<int64_t, vector<int32_t> > > positions(kmers.size());
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

const int balanced_stride(int read_length, int kmer_size, int stride) {
    double r = read_length;
    double k = kmer_size;
    double j = stride;
    return round((r-k)/round((r-k)/j));
}

const vector<string> balanced_kmers(const string& seq, const int kmer_size, const int stride) {
    // choose the closest stride that will generate balanced kmers
    vector<string> kmers;
    int b = balanced_stride(seq.size(), kmer_size, stride);
    if (!seq.empty()) {
        for (int i = 0; i < seq.size()-kmer_size; i+=b) {
            kmers.push_back(seq.substr(i,kmer_size));
        }
    }
    return kmers;
}

}

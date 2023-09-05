#include "kmer.hpp"

#include <atomic>

//#define debug

namespace vg {


const string SizeLimitExceededException::msg = "error: exceeded limit of size on disk";

const char* SizeLimitExceededException::what() const throw() {
    return msg.c_str();
}

void for_each_kmer(const HandleGraph& graph, size_t k,
                   const function<void(const kmer_t&)>& lambda,
                   id_t head_id, id_t tail_id, atomic<int>* stop_flag) {
    // for each position on the forward and reverse of the graph
    // TODO -- add parallel interface in handlegraph
    bool using_head_tail = head_id + tail_id > 0;
#ifdef debug
    cerr << "Looping over kmers" << endl;
#endif
    graph.for_each_handle([&](const handle_t& h) {
#ifdef debug
        cerr << "Process handle " << graph.get_id(h) << endl;
#endif
        // for the forward and reverse of this handle
        // walk k bases from the end, so that any kmer starting on the node will be represented in the tree we build
        for (auto handle_is_rev : { false, true }) {
            handle_t handle = handle_is_rev ? graph.flip(h) : h;
            list<kmer_t> kmers;
            // for each position in the node, set up a kmer with that start position and the node end or kmer length as the end position
            // determine next positions
            id_t handle_id = graph.get_id(handle);
            size_t handle_length = graph.get_length(handle);
            string handle_seq = graph.get_sequence(handle);
            for (size_t i = 0; i < handle_length; ++i) {
                pos_t begin = make_pos_t(handle_id, handle_is_rev, i);
                pos_t end = make_pos_t(handle_id, handle_is_rev, min(handle_length, i+k));
                kmer_t kmer = kmer_t(handle_seq.substr(offset(begin), offset(end)-offset(begin)), begin, end, handle);
                // determine previous context
                // if we are running with head/tail nodes, we'll need to do some trickery to eliminate the reverse complement versions of both
                if (i == 0) {
                    // look at the previous nodes
                    graph.follow_edges(handle, true, [&](const handle_t& prev) {
                        size_t prev_length = graph.get_length(prev);
                        kmer.prev_pos.emplace_back(graph.get_id(prev), graph.get_is_reverse(prev), prev_length-1);
                        kmer.prev_char.emplace_back(graph.get_sequence(prev).substr(prev_length-1, 1)[0]);
                        if (stop_flag) {
                            // stop if it evaluates to true
                            return !stop_flag->load();
                        }
                        else {
                            // always keep going
                            return true;
                        }
                    });
                    // if we're on the forward head or reverse tail, we need to point to the end of the opposite node
                    if (kmer.prev_pos.empty() && using_head_tail) {
                        if (id(begin) == head_id) {
                            kmer.prev_pos.emplace_back(tail_id, false, 0);
                            kmer.prev_char.emplace_back(graph.get_sequence(graph.get_handle(tail_id, false))[0]);
                        } else if (id(begin) == tail_id) {
                            kmer.prev_pos.emplace_back(head_id, true, 0);
                            kmer.prev_char.emplace_back(graph.get_sequence(graph.get_handle(head_id, true))[0]);
                        }
                    }
                } else {
                    // the previous is in this node
                    kmer.prev_pos.emplace_back(handle_id, handle_is_rev, i-1);
                    kmer.prev_char.emplace_back(handle_seq[i-1]);
                }
                if (kmer.seq.size() < k) {
                    kmer.seq.reserve(k); // may reduce allocation costs
                                         // follow edges if we haven't completed the kmer here
                    graph.follow_edges(kmer.curr, false, [&](const handle_t& next) {
                        kmers.push_back(kmer);
                        auto& todo = kmers.back();
                        todo.curr = next;
                        if (stop_flag) {
                            // stop if it evaluates to true
                            return !stop_flag->load();
                        }
                        else {
                            // always keep going
                            return true;
                        }
                    });
                } else {
                    kmers.push_back(kmer);
                }
                
                if (stop_flag && stop_flag->load()) {
                    break;
                }
            }
            
            // now expand the kmers until they reach k
            while (!kmers.empty()) {
                // first we check which ones have reached length k in the current handle; for each of these we run lambda and remove them from our list
                auto kmers_end = kmers.end();
                for (list<kmer_t>::iterator q = kmers.begin(); q != kmers_end; ++q) {
                    auto& kmer = *q;
                    // did we reach our target length?
                    if (kmer.seq.size() == k) {
                        // TODO here check if we are at the beginning of the reverse head or the beginning of the forward tail and would need special handling
                        // establish the context
                        handle_t end_handle = graph.get_handle(id(kmer.end), is_rev(kmer.end));
                        size_t end_length = graph.get_length(end_handle);
                        if (offset(kmer.end) == end_length) {
                            // have to check which nodes are next
                            graph.follow_edges(kmer.curr, false, [&](const handle_t& next) {
                                kmer.next_pos.emplace_back(graph.get_id(next), graph.get_is_reverse(next), 0);
                                kmer.next_char.emplace_back(graph.get_sequence(next)[0]);
                            });
                            if (kmer.next_pos.empty() && using_head_tail) {
                                if (id(kmer.begin) == head_id) {
                                    kmer.next_pos.emplace_back(tail_id, true, 0);
                                    kmer.next_char.emplace_back(graph.get_sequence(graph.get_handle(tail_id, true))[0]);
                                } else if (id(kmer.begin) == tail_id) {
                                    kmer.next_pos.emplace_back(head_id, false, 0);
                                    kmer.next_char.emplace_back(graph.get_sequence(graph.get_handle(head_id, false))[0]);
                                }
                                //cerr << "done head or tail" << endl;
                            }
                        } else {
                            // on node
                            kmer.next_pos.push_back(kmer.end);
                            kmer.next_char.push_back(graph.get_sequence(end_handle)[offset(kmer.end)]);
                        }
                        // if we have head and tail ids set, iterate through our positions and do the flip
                        if (using_head_tail) {
                            // flip the beginning
                            if (id(kmer.begin) == head_id && is_rev(kmer.begin)) {
                                get_id(kmer.begin) = tail_id;
                                get_is_rev(kmer.begin) = false;
                            } else if (id(kmer.begin) == tail_id && is_rev(kmer.begin)) {
                                get_id(kmer.begin) = head_id;
                                get_is_rev(kmer.begin) = false;
                            }
                            // flip the nexts
                            for (auto& pos : kmer.next_pos) {
                                if (id(pos) == head_id && is_rev(pos)) {
                                    get_id(pos) = tail_id;
                                    get_is_rev(pos) = false;
                                } else if (id(pos) == tail_id && is_rev(pos)) {
                                    get_id(pos) = head_id;
                                    get_is_rev(pos) = false;
                                }
                            }
                            // if we aren't both from and to a head/tail node, emit
                            /*
                             if (!((offset(kmer.begin) == 0
                             && id(kmer.begin) == head_id
                             && kmer.next_pos.size() == 1
                             && id(kmer.next_pos.front()) == tail_id)
                             || (offset(kmer.begin) == 0
                             && id(kmer.begin) == tail_id
                             && kmer.next_pos.size() == 1
                             && id(kmer.next_pos.front()) == head_id))) {
                             lambda(kmer);
                             }
                             */
                            if (kmer.prev_pos.size() == 1 && kmer.next_pos.size() == 1
                                && (offset(kmer.begin) == 0)
                                && (id(kmer.begin) == head_id || id(kmer.begin) == tail_id)
                                && (id(kmer.prev_pos.front()) == head_id || id(kmer.prev_pos.front()) == tail_id)
                                && (id(kmer.next_pos.front()) == head_id || id(kmer.next_pos.front()) == tail_id)) {
                                // skip
                            } else {
                                lambda(kmer);
                            }
                        } else {
                            // now pass the kmer and its context to our callback
                            lambda(kmer);
                        }
                        q = kmers.erase(q);
                    } else {
                        // do we finish in the current node?
                        id_t curr_id = graph.get_id(kmer.curr);
                        size_t curr_length = graph.get_length(kmer.curr);
                        bool curr_is_rev = graph.get_is_reverse(kmer.curr);
                        string curr_seq = graph.get_sequence(kmer.curr);
                        size_t take = min(curr_length, k-kmer.seq.size());
                        kmer.end = make_pos_t(curr_id, curr_is_rev, take);
                        kmer.seq.append(curr_seq.substr(0,take));
                        if (kmer.seq.size() < k) {
                            // if not, we need to expand through the node then follow on
                            graph.follow_edges(kmer.curr, false, [&](const handle_t& next) {
                                kmers.push_back(kmer);
                                auto& todo = kmers.back();
                                todo.curr = next;
                            });
                            q = kmers.erase(q);
                        } else {
                            if (kmer.seq.size() > k) {
                                assert(kmer.seq.size() <= k);
                            }
                        }
                    }
                }
                if (stop_flag && stop_flag->load()) {
                    break;
                }
            }
            if (stop_flag && stop_flag->load()) {
                break;
            }
        }
        if (stop_flag) {
            // stop if it evaluates to true
            return !stop_flag->load();
        }
        else {
            // always keep going
            return true;
        }
    }, true);
}

ostream& operator<<(ostream& out, const kmer_t& kmer) {
    out << kmer.seq << "\t"
        << id(kmer.begin) << ":" << (is_rev(kmer.begin) ? "-":"") << offset(kmer.begin) << "\t";
    for (size_t i = 0; i < kmer.prev_char.size(); ++i) {
        if (i != 0) out << ",";
        out << kmer.prev_char[i];
    }
    out << "\t";
    for (size_t i = 0; i < kmer.next_char.size(); ++i) {
        if (i != 0) out << ",";
        out << kmer.next_char[i];
    }
    out << "\t";
    for (size_t i = 0; i < kmer.next_pos.size(); ++i) {
        if (i != 0) out << ",";
        auto& pos = kmer.next_pos[i];
        out << id(pos) << ":" << (is_rev(pos) ? "-":"") << offset(pos);
    }
    return out;
}

void kmer_to_gcsa_kmers(const kmer_t& kmer, const gcsa::Alphabet& alpha, const function<void(const gcsa::KMer&)>& lambda) {
    assert(kmer.next_pos.size());
    gcsa::byte_type predecessors = encode_chars(kmer.prev_char, alpha);
    gcsa::byte_type successors = encode_chars(kmer.next_char, alpha);
    gcsa::KMer k;
    k.key = gcsa::Key::encode(alpha, kmer.seq, predecessors, successors);
    if (offset(kmer.begin) >= 1024) {
#pragma omp critical (error)
        {
            cerr << "error: Found kmer with offset >= 1024. GCSA2 cannot handle nodes greater than 1024 bases long. "
                 << "To enable indexing, modify your graph using `vg mod -X 256 x.vg >y.vg`. "
                 << kmer << endl;
            exit(1);
        }
    }
    k.from = gcsa::Node::encode(id(kmer.begin), offset(kmer.begin), is_rev(kmer.begin));
    for (auto& pos : kmer.next_pos) {
        k.to = gcsa::Node::encode(id(pos), offset(pos), is_rev(pos));
        lambda(k);
    }
}

gcsa::byte_type encode_chars(const vector<char>& chars, const gcsa::Alphabet& alpha) {
    gcsa::byte_type val = 0;
    for (char c : chars) val |= 1 << alpha.char2comp[c];
    return val;
}

void write_gcsa_kmers(const HandleGraph& graph, int kmer_size, ostream& out, size_t& size_limit, id_t head_id, id_t tail_id) {

    // We need an alphabet to parse the internal string format
    const gcsa::Alphabet alpha;
    // Each thread is going to make its own KMers, then we'll concatenate these all together at the end.
    vector<vector<gcsa::KMer> > thread_outputs;
#pragma omp parallel
    {
#pragma omp single
        {
            // Set up our write buffers at the given parallelism we expect
            thread_outputs.resize(omp_get_num_threads());
        }
    }
    
    // we can't throw from within an OMP block, so instead we have to use some machinery to flag when
    // we need to throw
    atomic<int> size_limit_exceeded(0);
    
    // This handles the buffered writing for each thread
    size_t buffer_limit = 1e5; // max 100k kmers per buffer
    size_t total_bytes = 0;
    auto handle_kmers = [&](vector<gcsa::KMer>& kmers, bool more) {
        if (!more || kmers.size() > buffer_limit) {
            size_t bytes_required = kmers.size() * sizeof(gcsa::KMer) + sizeof(gcsa::GraphFileHeader);
#pragma omp critical
            {
                if (!size_limit_exceeded.load()) {
                    // we didn't exceed the size limit while waiting for the critical block
                    if (total_bytes + bytes_required > size_limit) {
                        cerr << "error: [write_gcsa_kmers()] size limit of " << size_limit << " bytes exceeded" << endl;
                        size_limit_exceeded.store(1);
                    }
                    else {
                        gcsa::writeBinary(out, kmers, kmer_size);
                        total_bytes += bytes_required;
                    }
                }
            }
            kmers.clear();
        }
    };
    // Here we convert our kmer_t to gcsa::KMer
    auto convert_kmer = [&thread_outputs, &alpha, &head_id, &tail_id, &handle_kmers](const kmer_t& kmer) {
        // Convert this KmerPosition to several gcsa::KMers, and save them in thread_outputs
        vector<gcsa::KMer>& thread_output = thread_outputs[omp_get_thread_num()];
        kmer_to_gcsa_kmers(kmer, alpha, [&thread_output](const gcsa::KMer& k) { thread_output.push_back(k); });
        // Handle kmer buffered writes, indicating we're not yet done
        handle_kmers(thread_output, true);
    };
    // Run on each KmerPosition. This populates start_end_id, if it was 0, before calling convert_kmer.
    for_each_kmer(graph, kmer_size, convert_kmer, head_id, tail_id, &size_limit_exceeded);
    for(auto& thread_output : thread_outputs) {
        // Flush our buffers
        handle_kmers(thread_output, false);
    }
    // did we end execution because we hit the size limit
    if (size_limit_exceeded.load()) {
        throw SizeLimitExceededException();
    }
    
    // FIXME: we seem to use this behavior in VGset, but this is not good semantics
    size_limit = total_bytes;
}

string write_gcsa_kmers_to_tmpfile(const HandleGraph& graph, int kmer_size, size_t& size_limit, id_t head_id, id_t tail_id,
                                   const string& base_file_name) {
    // open a temporary file for the kmers
    string tmpfile = temp_file::create(base_file_name);
    ofstream out(tmpfile);
    // write the kmers to the temporary file
    try {
        write_gcsa_kmers(graph, kmer_size, out, size_limit, head_id, tail_id);
    }
    catch (SizeLimitExceededException& ex) {
        out.close();
        temp_file::remove(tmpfile);
        throw ex;
    }
    out.close();
    return tmpfile;
}



}

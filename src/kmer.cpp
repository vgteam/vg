#include "kmer.hpp"

namespace vg {

void for_each_kmer(const HandleGraph& graph, size_t k,
                   const function<void(const kmer_t&)>& lambda,
                   id_t head_id, id_t tail_id) {
    // for each position on the forward and reverse of the graph
    // TODO -- add parallel interface in handlegraph
    bool using_head_tail = head_id + tail_id > 0;
    graph.for_each_handle([&](const handle_t& h) {
            // for the forward and reverse of this handle
            // walk k bases from the end, so that any kmer starting on the node will be represented in the tree we build
            for (auto handle_is_rev : { false, true }) {
                //cerr << "###########################################" << endl;
                handle_t handle = handle_is_rev ? graph.flip(h) : h;
                list<kmer_t> kmers;
                // for each position in the node, set up a kmer with that start position and the node end or kmer length as the end position
                // determine next positions
                id_t handle_id = graph.get_id(handle);
                size_t handle_length = graph.get_length(handle);
                string handle_seq = graph.get_sequence(handle);
                for (size_t i = 0; i < handle_length;  ++i) {
                    pos_t begin = make_pos_t(handle_id, handle_is_rev, i);
                    pos_t end = make_pos_t(handle_id, handle_is_rev, min(handle_length, i+k));
                    kmer_t kmer = kmer_t(handle_seq.substr(offset(begin), offset(end)-offset(begin)), begin, end, handle);
                    // determine previous context
                    if (i == 0) {
                        // look at the previous nodes
                        graph.follow_edges(handle, true, [&](const handle_t& prev) {
                                size_t prev_length = graph.get_length(prev);
                                kmer.prev_pos.emplace_back(graph.get_id(prev), graph.get_is_reverse(prev), prev_length-1);
                                kmer.prev_char.emplace_back(graph.get_sequence(prev).substr(prev_length-1, 1)[0]);
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
                            });
                    } else {
                        kmers.push_back(kmer);
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
                            // now pass the kmer and its context to our callback
                            lambda(kmer);
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
                }
            }
        });
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

}

#include "kmer.hpp"

namespace vg {

namespace algorithms {

void for_each_kmer(const HandleGraph& graph, size_t k, size_t edge_max,
                   const std::function<void(const kmer_t&)>& lambda) {
    graph.for_each_handle([&](const handle_t& h) {
            // for the forward and reverse of this handle
            // walk k bases from the end, so that any kmer starting on the node will be represented in the tree we build
            for (auto handle_is_rev : { false, true }) {
                handle_t handle = handle_is_rev ? graph.flip(h) : h;
                std::list<kmer_t> kmers;
                // for each position in the node, set up a kmer with that start position and the node end or kmer length as the end position
                // determine next positions
                nid_t handle_id = graph.get_id(handle);
                size_t handle_length = graph.get_length(handle);
                std::string handle_seq = graph.get_sequence(handle);
                for (size_t i = 0; i < handle_length;  ++i) {
                    pos_t begin = make_pos_t(handle_id, handle_is_rev, i);
                    pos_t end = make_pos_t(handle_id, handle_is_rev, std::min(handle_length, i+k));
                    kmer_t kmer = kmer_t(handle_seq.substr(offset(begin), offset(end)-offset(begin)), begin, end, handle);
                    if (kmer.seq.size() < k) {
                        size_t next_count = 0;
                        if (edge_max) graph.follow_edges(kmer.curr, false, [&](const handle_t& next) { ++next_count; return next_count <= 1; });
                        //kmer.seq.reserve(k); // may reduce allocation costs
                        // follow edges if we haven't completed the kmer here
                        if (next_count > 1 && (edge_max && edge_max == kmer.forks)) {
                        } else {
                            graph.follow_edges(kmer.curr, false, [&](const handle_t& next) {
                                    kmers.push_back(kmer);
                                    auto& todo = kmers.back();
                                    todo.curr = next;
                                    if (next_count > 1) {
                                        ++todo.forks;
                                    }
                            });
                        }
                    } else {
                        kmers.push_back(kmer);
                    }
                }

                // now expand the kmers until they reach k
                while (!kmers.empty()) {
                    // first we check which ones have reached length k in the current handle; for each of these we run lambda and remove them from our list
                    auto kmers_end = kmers.end();
                    for (std::list<kmer_t>::iterator q = kmers.begin(); q != kmers_end; ++q) {
                        auto& kmer = *q;
                        // did we reach our target length?
                        if (kmer.seq.size() == k) {
                            // now pass the kmer to our callback
                            lambda(kmer);
                            q = kmers.erase(q);
                        } else {
                            // do we finish in the current node?
                            nid_t curr_id = graph.get_id(kmer.curr);
                            size_t curr_length = graph.get_length(kmer.curr);
                            bool curr_is_rev = graph.get_is_reverse(kmer.curr);
                            std::string curr_seq = graph.get_sequence(kmer.curr);
                            size_t take = std::min(curr_length, k-kmer.seq.size());
                            kmer.end = make_pos_t(curr_id, curr_is_rev, take);
                            kmer.seq.append(curr_seq.substr(0,take));
                            if (kmer.seq.size() < k) {
                                size_t next_count = 0;
                                if (edge_max) graph.follow_edges(kmer.curr, false, [&](const handle_t& next) { ++next_count; return next_count <= 1; });
                                //kmer.seq.reserve(k); // may reduce allocation costs
                                // follow edges if we haven't completed the kmer here
                                if (next_count > 1 && (edge_max && edge_max == kmer.forks)) {
                                } else {
                                    graph.follow_edges(kmer.curr, false, [&](const handle_t& next) {
                                            kmers.push_back(kmer);
                                            auto& todo = kmers.back();
                                            todo.curr = next;
                                            if (next_count > 1) {
                                                ++todo.forks;
                                            }
                                        });
                                }
                                // if not, we need to expand through the node then follow on
                                /*
                                graph.follow_edges(kmer.curr, false, [&](const handle_t& next) {
                                        kmers.push_back(kmer);
                                        auto& todo = kmers.back();
                                        todo.curr = next;
                                    });
                                */
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
        }, true);
}

std::ostream& operator<<(std::ostream& out, const kmer_t& kmer) {
    out << kmer.seq << "\t"
        << id(kmer.begin) << ":" << (is_rev(kmer.begin) ? "-":"") << offset(kmer.begin) << "\t";
    return out;
}

}

}

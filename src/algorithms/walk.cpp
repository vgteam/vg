#include "walk.hpp"

#include <gbwtgraph/utils.h>

namespace vg {

namespace algorithms {

void for_each_walk(const HandleGraph& graph, size_t k, size_t edge_max,
                   const std::function<void(const walk_t&)>& lambda) {
    graph.for_each_handle([&](const handle_t& h) {
            // for the forward and reverse of this handle
            // walk k bases from the end, so that any walk starting on the node will be represented in the tree we build
            for (auto handle_is_rev : { false, true }) {
                handle_t handle = handle_is_rev ? graph.flip(h) : h;
                std::list<walk_t> walks;
                // for each position in the node, set up a walk with that start position and the node end or walk length as the end position
                // determine next positions
                nid_t handle_id = graph.get_id(handle);
                size_t handle_length = graph.get_length(handle);
                std::string handle_seq = graph.get_sequence(handle);
                for (size_t i = 0; i < handle_length;  ++i) {
                    pos_t begin = make_pos_t(handle_id, handle_is_rev, i);
                    pos_t end = make_pos_t(handle_id, handle_is_rev, std::min(handle_length, i+k));
                    walk_t walk = walk_t(handle_seq.substr(offset(begin), offset(end)-offset(begin)), begin, end, handle);
                    if (walk.seq.size() < k) {
                        size_t next_count = 0;
                        if (edge_max) graph.follow_edges(walk.curr, false, [&](const handle_t& next) { ++next_count; return next_count <= 1; });
                        //walk.seq.reserve(k); // may reduce allocation costs
                        // follow edges if we haven't completed the walk here
                        if (next_count > 1 && (edge_max && edge_max == walk.forks)) {
                        } else {
                            graph.follow_edges(walk.curr, false, [&](const handle_t& next) {
                                    walks.push_back(walk);
                                    auto& todo = walks.back();
                                    todo.curr = next;
                                    if (next_count > 1) {
                                        ++todo.forks;
                                    }
                            });
                        }
                    } else {
                        walks.push_back(walk);
                    }
                }

                // now expand the walks until they reach k
                while (!walks.empty()) {
                    // first we check which ones have reached length k in the current handle; for each of these we run lambda and remove them from our list
                    auto walks_end = walks.end();
                    for (std::list<walk_t>::iterator q = walks.begin(); q != walks_end; ++q) {
                        auto& walk = *q;
                        // did we reach our target length?
                        if (walk.seq.size() == k) {
                            // TODO here check if we are at the beginning of the reverse head or the beginning of the forward tail and would need special handling
                            // establish the context
                            handle_t end_handle = graph.get_handle(id(walk.end), is_rev(walk.end));
                            size_t end_length = graph.get_length(end_handle);
                            // now pass the walk to our callback
                            lambda(walk);
                            q = walks.erase(q);
                        } else {
                            // do we finish in the current node?
                            nid_t curr_id = graph.get_id(walk.curr);
                            size_t curr_length = graph.get_length(walk.curr);
                            bool curr_is_rev = graph.get_is_reverse(walk.curr);
                            std::string curr_seq = graph.get_sequence(walk.curr);
                            size_t take = std::min(curr_length, k-walk.seq.size());
                            walk.end = make_pos_t(curr_id, curr_is_rev, take);
                            walk.seq.append(curr_seq.substr(0,take));
			    walk.path.push_back(walk.curr);
                            if (walk.seq.size() < k) {
                                size_t next_count = 0;
                                if (edge_max) graph.follow_edges(walk.curr, false, [&](const handle_t& next) { ++next_count; return next_count <= 1; });
                                //walk.seq.reserve(k); // may reduce allocation costs
                                // follow edges if we haven't completed the walk here
                                if (next_count > 1 && (edge_max && edge_max == walk.forks)) {
                                } else {
                                    graph.follow_edges(walk.curr, false, [&](const handle_t& next) {
                                            walks.push_back(walk);
                                            auto& todo = walks.back();
                                            todo.curr = next;
                                            if (next_count > 1) {
                                                ++todo.forks;
                                            }
                                        });
                                }
                                // if not, we need to expand through the node then follow on
                                /*
                                graph.follow_edges(walk.curr, false, [&](const handle_t& next) {
                                        walks.push_back(walk);
                                        auto& todo = walks.back();
                                        todo.curr = next;
                                    });
                                */
                                q = walks.erase(q);
                            } else {
                                if (walk.seq.size() > k) {
                                    assert(walk.seq.size() <= k);
                                }
                            }
                        }
                    }
                }
            }
        }, true);
}

std::ostream& operator<<(std::ostream& out, const walk_t& walk) {
    out << walk.seq << "\t"
        << id(walk.begin) << ":" << (is_rev(walk.begin) ? "-":"") << offset(walk.begin) << "\t";
    return out;
}

uint64_t walk_haplotype_frequency(const HandleGraph& graph,
                                  const gbwt::GBWT& haplotypes,
                                  const walk_t& walk) {
    if (walk.path.empty()) {
        return 0;
    }
    auto& first_step = walk.path.front();
    gbwt::node_type start_node = gbwt::Node::encode(graph.get_id(first_step), graph.get_is_reverse(first_step));
    gbwt::SearchState search_state = haplotypes.find(start_node);
    for (uint64_t i = 1; i < walk.path.size(); ++i) {
        auto& next = walk.path[i];
        gbwt::node_type next_node = gbwt::Node::encode(graph.get_id(next), graph.get_is_reverse(next));
        search_state = haplotypes.extend(search_state, next_node);
        if (search_state.empty()) {
            break;
        }
    }
    return search_state.size();
}

std::vector<std::string> walk_haplotype_names(const HandleGraph& graph,
                                              const gbwt::GBWT& haplotypes,
                                              const walk_t& walk) {
    std::vector<std::string> names;
    if (walk.path.empty()) {
        return names;
    }
    auto& first_step = walk.path.front();
    gbwt::node_type start_node = gbwt::Node::encode(graph.get_id(first_step), graph.get_is_reverse(first_step));
    gbwt::SearchState search_state = haplotypes.find(start_node);
    for (uint64_t i = 1; i < walk.path.size(); ++i) {
        auto& next = walk.path[i];
        gbwt::node_type next_node = gbwt::Node::encode(graph.get_id(next), graph.get_is_reverse(next));
        search_state = haplotypes.extend(search_state, next_node);
        if (search_state.empty()) {
            break;
        }
    }
    assert(haplotypes.hasMetadata() && haplotypes.metadata.hasSampleNames());
    // Pre-parse the reference samples.
    // TODO: Can we pass this down?
    auto reference_samples = gbwtgraph::parse_reference_samples_tag(haplotypes);
    for (auto& thread : haplotypes.locate(search_state)) {
        // For each match
        auto id = gbwt::Path::id(thread);
        
        // Figure out what kind of path it is
        PathSense sense = gbwtgraph::get_path_sense(haplotypes, id, reference_samples);
        
        // Figure out its haplotype number
        auto haplotype = gbwtgraph::get_path_haplotype(haplotypes, id, sense);
        if (haplotype == PathMetadata::NO_HAPLOTYPE) {
            // If no haplotype is applicable, use 0.
            haplotype = 0;
        }
        
        // Compose a name for it
        std::stringstream ss;
        ss << gbwtgraph::get_path_sample_name(haplotypes, id, sense) << "#" << haplotype;
        names.push_back(ss.str());
        // TODO: should we do something special for generic sense or NO_SAMPLE_NAME?
    }
    return names;
}


}

}

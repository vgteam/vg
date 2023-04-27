/**
 * \file load_proto_to_graph.cpp
 * Implementation for backend-agnostic Protobuf loading.
 */

#include "load_proto_to_graph.hpp"

#include "../hash_map.hpp"
#include "../handle.hpp"
#include "../crash.hpp"

#include "load_proto_to_graph.hpp"

#include "vg/io/json2pb.h"

#include <vg/vg.pb.h>
#include <vg/io/registry.hpp>
#include <vg/io/protobuf_iterator.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>

#include <unordered_map>
#include <mutex>

//#define debug

namespace vg {

namespace io {

using namespace std;
using namespace vg;

void load_proto_to_graph(vg::MutablePathMutableHandleGraph* destination, const vg::io::message_sender_function_t& for_each_message) {
    
    load_proto_to_graph(destination, [&](const function<void(Graph&)>& process_chunk) {
        // Now we can give all the deserialized Graph chunks to process_chunk.
        for_each_message([&](const string& serialized_graph) {
            // For each Graph, unpack it
            Graph g;
            if (!ProtobufIterator<Graph>::parse_from_string(g, serialized_graph)) {
                // TODO: make this an exception if we ever want to be allowed to continue from this.
                cerr << "error[load_proto_to_graph]: invalid Graph message" << endl;
                exit(1);
            }
            
            // And send it along.
            process_chunk(g);
        });
    });
}

void load_proto_to_graph(vg::MutablePathMutableHandleGraph* destination, const function<void(const function<void(Graph&)>&)>& chunk_sender) {
    
    // This holds edges we couldn't make when we saw them because both end nodes didn't exist.
    // We organize them by the ID of the node they were waiting for.
    hash_map<nid_t, vector<tuple<nid_t, bool, nid_t, bool>>> deferred_edges;
    
    // This represents a path in progress
    struct path_record_t {
        size_t min_rank_added = 0;
        size_t max_rank_added = 0;
        map<size_t, pair<nid_t, bool>> earlier_visits_by_rank;
        map<size_t, pair<nid_t, bool>> later_visits_by_rank;
    };
    
    // This holds, for each path, the min and max ranks added to the
    // graph, and the trees of everything not added yet above and below those
    // ranks.
    hash_map<path_handle_t, path_record_t> paths_in_progress;
    
    // We can only deal with one chunk coming in at a time, but we have
    // ambitions for a parallel Constructor. So make sure to only handle one
    // chunk at a time.
    mutex chunk_mutex;
    
    chunk_sender([&](Graph& g) {
        // For each Graph chunk...
        
        // Handle one at a time
        lock_guard<mutex> chunk_guard(chunk_mutex);
        
        // Within this chunk, we keep a node to handle cache
        unordered_map<nid_t, handle_t> node_to_handle;
        
        // Define a way to get a handle from a cached handle, or the graph, or to fail.
        auto get_handle = [&](nid_t id, bool is_reverse, handle_t& dest) {
             auto handle_iter = node_to_handle.find(id);
             if (handle_iter == node_to_handle.end()) {
                // Handle is not cached.
                if (destination->has_node(id)) {
                    // But we do have the node in the destination graph already.
                    dest = destination->get_handle(id, is_reverse);
                    crash_unless(destination->get_id(dest) == id);
                    crash_unless(destination->get_is_reverse(dest) == is_reverse);
                    return true;
                } else {
                    // The node doesn't exist yet.
                    return false;
                }
            } else {
                // Handle is cached. All we have to do is get the correct orientation.
                dest = is_reverse ? destination->flip(handle_iter->second) : handle_iter->second;
                crash_unless(destination->get_id(dest) == id);
                crash_unless(destination->get_is_reverse(dest) == is_reverse);
                return true;
            }
        };
        
        for (auto& n : g.node()) {
            // Create all the nodes
            
#ifdef debug
                cerr << "Found node " << n.id() << endl;
#endif
            
            node_to_handle[n.id()] = destination->create_handle(n.sequence(), n.id());
            
            auto edges_waiting = deferred_edges.find(n.id());
            if (edges_waiting != deferred_edges.end()) {
            
#ifdef debug
                cerr << "Node has deferred edges waiting on it" << endl;
#endif
            
                // There were edges waiting for this node. Make them and get them out of our memory.
                for (auto& edge : edges_waiting->second) {
                
#ifdef debug
                    cerr << "Make deferred edge " << get<0>(edge) << " " << get<1>(edge)
                        << " " << get<2>(edge) << " " << get<3>(edge) << endl;
#endif
                
                    // For each edge that only lacked this node, get the handles.
                    handle_t from_handle;
                    handle_t to_handle;
                    // We can assert because if we did our bookkepping right both nodes exist now.
                    crash_unless(get_handle(get<0>(edge), get<1>(edge), from_handle));
                    crash_unless(get_handle(get<2>(edge), get<3>(edge), to_handle));
                    // Make the edge
                    destination->create_edge(from_handle, to_handle);
                }
                
                // Forget all those edges we just made.
                deferred_edges.erase(edges_waiting);
            }
        }
        
        for (auto& e : g.edge()) {
            // For each edge
            
            // See if we have the handle for each end cached
            auto from_handle_iter = node_to_handle.find(e.from());
            auto to_handle_iter = node_to_handle.find(e.to());
            
            // Get the actual handle for each end, from the cache or graph, if available.
            // If not available, defer the edge.
            handle_t from_handle;
            handle_t to_handle;
            if (!get_handle(e.from(), e.from_start(), from_handle)) {
                // From node isn't made yet. Defer the edge.
                
#ifdef debug
                cerr << "Defer edge " << e.from() << " " << e.from_start() << " "
                     << e.to() << " " << e.to_end() << " on " << e.from() << endl;
#endif
                
                deferred_edges[e.from()].emplace_back(e.from(), e.from_start(), e.to(), e.to_end());
                // Don't do the edge now.
                continue;
            } else if (!get_handle(e.to(), e.to_end(), to_handle)) {
                // To node isn't made yet. Defer the edge.
                
#ifdef debug
                cerr << "Defer edge " << e.from() << " " << e.from_start() << " "
                     << e.to() << " " << e.to_end() << " on " << e.to() << endl;
#endif
                
                deferred_edges[e.to()].emplace_back(e.from(), e.from_start(), e.to(), e.to_end());
                // Don't do the edge now.
                continue;
            }
            
#ifdef debug
            cerr << "Make edge " << e.from() << " " << e.from_start() << " "
                << e.to() << " " << e.to_end() << endl;
#endif
            
            // Make the edge
            destination->create_edge(from_handle, to_handle);
        }
        
        // For paths, we have append_step and prepend_step, so we can add path
        // stuff to the graph as long as we have a continuous rank range.
        // If we get stuff not in a contiguous rank range, we collate it with a map.
        
        // Path visits can only occur in or after chunks their nodes and edges are in.
        
        for (auto& p : g.path()) {
            // We're going to find or place our path in our index of records.
            path_handle_t path;
            if (!destination->has_path(p.name())) {
                // We need to create a new path!
#ifdef debug
                cerr << "Found new path " << p.name() << endl;
#endif
                path = destination->create_path_handle(p.name(), p.is_circular());
            } else {
                // We need to extend an existing path
#ifdef debug
                cerr << "Found existing path " << p.name() << endl;
#endif
                path = destination->get_path_handle(p.name());
                
                if (p.is_circular() && !destination->get_is_circular(path)) {
                    // If it ever shows up as circular, make it curcular
                    destination->set_circularity(path, true);
                }
            }
            
            // Find or make the record for the path in progress.
            //
            // We only need to read records for paths that are split across
            // chunks, but we always need to write records in case we need to
            // read them.
            auto& record = paths_in_progress[path];
                
            for (auto& m : p.mapping()) {
                // For each mapping in the path part for this chunk, in order
                if (m.rank() != 0) {
                    
#ifdef debug
                    cerr << "Found ranked mapping (" << m.rank() << ")" << endl;
#endif
                
                    // If it has a rank
                    if (record.min_rank_added == 0) {
                        // If it is the first rank, put it in.
                        
#ifdef debug
                        cerr << "Add as first mapping" << endl;
#endif
                        
                        handle_t visited;
                        crash_unless(get_handle(m.position().node_id(), m.position().is_reverse(), visited));
                        destination->append_step(path, visited);
                        
                        // And save its rank as the only one added
                        record.min_rank_added = m.rank();
                        record.max_rank_added = record.min_rank_added;
                    } else if(record.max_rank_added + 1 == m.rank()) {
                        // If it is adjacent to the previous rank on the high side (most common case)
                        
#ifdef debug
                        cerr << "Add as adjacent mapping" << endl;
#endif
                        
                        // Add it
                        handle_t visited;
                        if (!get_handle(m.position().node_id(), m.position().is_reverse(), visited)) {
                            throw std::runtime_error("Could not find node " + std::to_string(m.position().node_id()) + " " + (m.position().is_reverse() ? "-" : "+") + " to add mapping " + pb2json(m) + " adjacent to existing mapping at rank " + std::to_string(record.max_rank_added));
                        }
                        destination->append_step(path, visited);
                        
                        // And update the ranks
                        record.max_rank_added = m.rank();
                    } else {
                        // If it isn't adjacent on the high side, stick it in the appropriate ordered map.
                        // We can resolve these later.
                        
#ifdef debug
                        cerr << "Save for later" << endl;
#endif
                        
                        if (m.rank() >= record.min_rank_added && m.rank() <= record.max_rank_added) {
                            // Prohibit duplicate ranks in our contiguous, added region.
                            // Note that we may miss them if they are not in the contiguous region when we see them.
                            cerr << "error[load_proto_to_graph]: duplicate rank " << m.rank() << " in path " << p.name() << endl;
                            exit(1);
                        }
                        
                        // Decide on which side of the contiguous region we go
                        auto& dest_map = (m.rank() < record.min_rank_added) ? record.earlier_visits_by_rank : record.later_visits_by_rank;
                        // Add the visit in
                        dest_map.emplace(m.rank(), make_pair(m.position().node_id(), m.position().is_reverse())); 
                    }
                } else {
                    // If it has no rank, just stick it into the path when it occurs in the file.
                    // Mixing ranked and unranked in the same path is undefined behavior.
                    
#ifdef debug
                    cerr << "Found unranked ranked mapping" << endl;
#endif
                    
                    // Make sure we have the node we are visiting
                    handle_t visited;
                    crash_unless(get_handle(m.position().node_id(), m.position().is_reverse(), visited));
                    
                    // Make a step to it.
                    destination->append_step(path, visited);
                }
            }
            
            // Now resolve any out-of-order ranked mappings that we can. Maybe we filled in a gap.
            {
                auto it = record.earlier_visits_by_rank.rbegin();
                while(it != record.earlier_visits_by_rank.rend() && it->first + 1 == record.min_rank_added) {
                    // We have something to prepend
                    
#ifdef debug
                    cerr << "Resolve earlier mapping (" << it->first << ")" << endl;
#endif
                    
                    // Get the handle
                    handle_t visited;
                    crash_unless(get_handle(it->second.first, it->second.second, visited));
                    
                    // Prepend it
                    destination->prepend_step(path, visited);
                    
                    // Update rank
                    record.min_rank_added = it->first;
                    
                    // Drop the visit. Because we have a reverse iterator, we
                    // need to get the next thing (possibly rend) and convert
                    // to a forward iterator. See <https://stackoverflow.com/a/1830240>
                    record.earlier_visits_by_rank.erase(std::next(it).base());
                    
                    // Find the next one
                    it = record.earlier_visits_by_rank.rbegin();
                }
            }
            {
                auto it = record.later_visits_by_rank.begin();
                while(it != record.later_visits_by_rank.end() && it->first == record.max_rank_added + 1) {
                    // We have something to append
                    
#ifdef debug
                    cerr << "Resolve later mapping (" << it->first << ")" << endl;
#endif
                    
                    // Get the handle
                    handle_t visited;
                    crash_unless(get_handle(it->second.first, it->second.second, visited));
                    
                    // Prepend it
                    destination->append_step(path, visited);
                    
                    // Update rank
                    record.max_rank_added = it->first;
                    
                    // Drop the visit
                    record.later_visits_by_rank.erase(it);
                    
                    // Find the next one
                    it = record.later_visits_by_rank.begin();
                }
            }
        }
    });
    
    if (!deferred_edges.empty()) {
        // If there are any deferred edges left, we are missing a node.
        // Sometimes we have to deal with dangling edges. We just remove them.
        cerr << "warning[load_proto_to_graph]: dangling edges on missing node " << deferred_edges.begin()->first
            << " and " << (deferred_edges.size() - 1) << " other missing nodes removed" << endl;
    }
    
    // Now make all the path steps we didn't make yet, allowing for rank gaps
    
    for (auto& path_and_record : paths_in_progress) {
        // For each path and record that we touched.
        // TODO: this could be a lot of them when we have alt paths!
        auto& path = path_and_record.first;
        auto& record = path_and_record.second;
        
        for (auto it = record.earlier_visits_by_rank.rbegin(); it != record.earlier_visits_by_rank.rend(); ++it) {
            // For each earlier thing we still have, add it
#ifdef debug
            cerr << "Resolve final earlier mapping (" << it->first << ")" << endl;
#endif
            handle_t visited = destination->get_handle(it->second.first, it->second.second);
            destination->prepend_step(path, visited);
        }
        record.earlier_visits_by_rank.clear();
        
        for (auto it = record.later_visits_by_rank.begin(); it != record.later_visits_by_rank.end(); ++it) {
            // For each earlier thing we still have, add it
#ifdef debug
            cerr << "Resolve final later mapping (" << it->first << ")" << endl;
#endif
            handle_t visited = destination->get_handle(it->second.first, it->second.second);
            destination->append_step(path, visited);
        }
        record.later_visits_by_rank.clear();
    }
    
    // Now we're done!
}

}

}

//#define debug_distance_indexing

#include "snarl_distance_index.hpp"

using namespace std;
using namespace handlegraph;
namespace vg {

bool check_regularity(const SnarlDistanceIndex::TemporaryDistanceIndex& temp_index, const SnarlDistanceIndex::temp_record_ref_t& snarl_index, const SnarlDistanceIndex::TemporaryDistanceIndex::TemporarySnarlRecord& temp_snarl_record, const vector<SnarlDistanceIndex::temp_record_ref_t>& all_children, const HandleGraph* graph) {
#ifdef debug_distance_indexing
    std::cerr << "Check if snarl " << temp_snarl_record.start_node_id << " to " << temp_snarl_record.end_node_id << " with " << all_children.size() << " children is regular" << std::endl;
#endif

    if (temp_snarl_record.is_root_snarl) {
        // Roots can't be regular.
#ifdef debug_distance_indexing
        std::cerr << "Snarl is not regular because it is a root snarl." << std::endl;
#endif
        return false;
    }
    if (temp_snarl_record.is_simple) {
        // Simple snarls are always also regular.
#ifdef debug_distance_indexing
        std::cerr << "Snarl is regular because it is simple." << std::endl;
#endif
        return true;
    }

    // Get the snarl boundary nodes, facing out
    handle_t start_out = graph->get_handle(temp_snarl_record.start_node_id, !temp_snarl_record.start_node_rev);
    handle_t end_out = graph->get_handle(temp_snarl_record.end_node_id, temp_snarl_record.end_node_rev);

    // Define accessors to get bounding graph handles for children, facing out.
    auto child_start_out = [&](const SnarlDistanceIndex::temp_record_ref_t& child_index) {
        return child_index.first == SnarlDistanceIndex::TEMP_NODE ? 
            graph->get_handle(child_index.second, true) :
            graph->get_handle(
                temp_index.get_chain(child_index).start_node_id,
                !temp_index.get_chain(child_index).start_node_rev
            );
    };
    auto child_end_out = [&](const SnarlDistanceIndex::temp_record_ref_t& child_index) {
        return child_index.first == SnarlDistanceIndex::TEMP_NODE ? 
            graph->get_handle(child_index.second, false) :
            graph->get_handle(
                temp_index.get_chain(child_index).end_node_id,
                temp_index.get_chain(child_index).end_node_rev
            );
    };

    for (const SnarlDistanceIndex::temp_record_ref_t& child_index : all_children) {
        // We should only have nodes and chains as children
        assert(child_index.first == SnarlDistanceIndex::TEMP_NODE
            || child_index.first == SnarlDistanceIndex::TEMP_CHAIN);
        if (child_index.first == SnarlDistanceIndex::TEMP_NODE
            && (child_index.second == temp_snarl_record.start_node_id
                || child_index.second == temp_snarl_record.end_node_id)) {
            // Don't think about children for the snarl bounds now; we handle the bounds later.
            continue;
        }

        // Have we seen the snarl start?
        bool saw_start = false;
        // Have we seen the snarl end?
        bool saw_end = false;
        // Have we seen anything else, or a duplicate snarl boundary?
        bool saw_other = false;

        auto handle_destination = [&](const handle_t& next_handle) {
#ifdef debug_distance_indexing
            std::cerr << "\tConnects to " << graph->get_id(next_handle) << (graph->get_is_reverse(next_handle) ? "-" : "+") << std::endl;
#endif

            // Every edge out the end the child must go to a snarl boundary out
            // that hasn't been reached yet.
            if (next_handle == start_out && !saw_start) {
                saw_start = true;
#ifdef debug_distance_indexing
                std::cerr << "\t\tThis is a new connection to snarl start" << std::endl;
#endif
                return true;
            } else if (next_handle == end_out && !saw_end) {
                saw_end = true;
#ifdef debug_distance_indexing
                std::cerr << "\t\tThis is a new connection to snarl end" << std::endl;
#endif
                return true;
            } else {
                saw_other = true;
                // We don't care if we have an edge going the right way because
                // we found an edge going the wrong way.
#ifdef debug_distance_indexing
                std::cerr << "\t\tThis is an unwanted connection!" << std::endl;
#endif
                return false;
            }
        };
        
        // Check the edges off the child start
        handle_t here = child_start_out(child_index);
#ifdef debug_distance_indexing
            std::cerr << "Look right from " << graph->get_id(here) << (graph->get_is_reverse(here) ? "-" : "+") << std::endl;
#endif
        graph->follow_edges(here, false, handle_destination);

        if (saw_other || !(saw_start != saw_end)) {
            // We have an edge we shouldn't, or we don't connect to exactly one boundary.
#ifdef debug_distance_indexing
            std::cerr << "\tWe must not be regular" << std::endl;
#endif
            return false;
        }
        
        // Check the edges off the child end
        here = child_end_out(child_index);
#ifdef debug_distance_indexing
            std::cerr << "Look right from " << graph->get_id(here) << (graph->get_is_reverse(here) ? "-" : "+") << std::endl;
#endif
        graph->follow_edges(here, false, handle_destination);

        if (saw_other || !saw_start || !saw_end) {
            // We have an edge we shouldn't, or we haven't reached both
            // boundaries exactly once across the two ends of the child.
#ifdef debug_distance_indexing
            std::cerr << "\tWe must not be regular" << std::endl;
#endif
            return false;
        }

        if (child_index.first == SnarlDistanceIndex::TEMP_CHAIN) {
            // If a child is a chain, check it for loops
#ifdef debug_distance_indexing
            std::cerr << "Check child chain for loops." << std::endl;
#endif
            const SnarlDistanceIndex::TemporaryDistanceIndex::TemporaryChainRecord& temp_chain_record = temp_index.get_chain(child_index);
#ifdef debug_distance_indexing
            std::cerr << "Forward loops:";
            for (auto& l : temp_chain_record.forward_loops) {
                std::cerr << " " << l;
            }
            std::cerr << std::endl;
#endif

            if (!temp_chain_record.forward_loops.empty() && temp_chain_record.forward_loops.front() != std::numeric_limits<size_t>::max()) {
                // There's a forward loop in this child chain, so the snarl's not regular.
#ifdef debug_distance_indexing
                std::cerr << "We are not regular because there's a forward loop in this child chain." << std::endl;
#endif
                return false;
            }

#ifdef debug_distance_indexing
            std::cerr << "Backward loops:";
            for (auto& l : temp_chain_record.backward_loops) {
                std::cerr << " " << l;
            }
            std::cerr << std::endl;
#endif

            if (!temp_chain_record.backward_loops.empty() && temp_chain_record.backward_loops.back() != std::numeric_limits<size_t>::max()) {
                // There's a backward loop in this child chain, so the snarl's not regular.
#ifdef debug_distance_indexing
                std::cerr << "We are not regular because there's a backward loop in this child chain." << std::endl;
#endif
                return false;
            }
        }
    }

    // Now we know the children are fine; check for disallowed edges between
    // the sentinels.

    handle_t start_in = graph->flip(start_out);
    if (graph->has_edge(start_in, start_out)) {
#ifdef debug_distance_indexing
        std::cerr << "We are not regular because we have a start-start loop." << std::endl;
#endif
        return false;
    }

    handle_t end_in = graph->flip(end_out);
    if (graph->has_edge(end_in, end_out)) {
#ifdef debug_distance_indexing
        std::cerr << "We are not regular because we have an end-end loop." << std::endl;
#endif
        return false;
    }

    // If we don't have any disallowed edges, and we don't have any children
    // without the exact right connectivity, we must be regular.

    // We don't make sure we actually had any children.
    
#ifdef debug_distance_indexing
    std::cerr << "We are a regular snarl." << std::endl;
#endif

    return true;
}

} // namespace vg

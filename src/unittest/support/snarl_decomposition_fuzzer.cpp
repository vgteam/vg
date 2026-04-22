#include "snarl_decomposition_fuzzer.hpp"

#include <cassert>
#include <stack>

namespace vg {
namespace unittest {

using ET = DecompositionEventType;

// SnarlDecompositionFuzzer deterministic constructor
SnarlDecompositionFuzzer::SnarlDecompositionFuzzer(
    const HandleGraph* graph,
    const HandleGraphSnarlFinder* finder,
    const std::unordered_set<nid_t>& chains_to_flip)
    : HandleGraphSnarlFinder(graph), wrapped(finder)
{

    should_flip = [chains_to_flip, graph](nid_t node_id) -> bool {
        return chains_to_flip.count(node_id);
    };
}

void SnarlDecompositionFuzzer::traverse_decomposition(
    const function<void(handle_t)>& begin_chain,
    const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl,
    const function<void(handle_t)>& end_snarl) const
{
    // Step 1: Capture all events from the wrapped finder.
    std::vector<DecompositionHandleEvent> events = capture_events(*wrapped);

    if (events.empty()) {
        return;
    }

    // Step 2: Build pairing vector mapping each begin to its matching end
    // and vice versa, using separate stacks for chains and snarls.
    std::vector<size_t> other_bound(events.size());
    {
        stack<size_t> chain_stack, snarl_stack;
        for (size_t i = 0; i < events.size(); i++) {
            switch (events[i].type) {
            case ET::BEGIN_CHAIN:
                chain_stack.push(i);
                break;
            case ET::END_CHAIN:
                assert(!chain_stack.empty());
                other_bound[i] = chain_stack.top();
                other_bound[chain_stack.top()] = i;
                chain_stack.pop();
                break;
            case ET::BEGIN_SNARL:
                snarl_stack.push(i);
                break;
            case ET::END_SNARL:
                assert(!snarl_stack.empty());
                other_bound[i] = snarl_stack.top();
                other_bound[snarl_stack.top()] = i;
                snarl_stack.pop();
                break;
            }
        }
    }

    // Step 3: Walk through events with a cursor, flipping chains as needed.
    // When we flip a chain, we jump to the other end and reverse direction,
    // pushing the entry point onto a stack. When the cursor reaches a stack
    // entry point, we jump back to the far end and restore direction.
    struct FlipEntry {
        size_t entry_index;
        bool original_reverse;
    };
    std::stack<FlipEntry> flip_stack;

    auto emitter = event_emitter(begin_chain, end_chain, begin_snarl, end_snarl);

    bool reverse = false;
    for (size_t cursor = 0; cursor != events.size(); cursor += reverse ? -1 : 1) {
        // We know if we're entering a chain, we can't be at a stack pop point.
        // So we can handle those cases separately.

        if (events[cursor].type == (reverse ? ET::END_CHAIN : ET::BEGIN_CHAIN) && 
            should_flip(graph->get_id(events[cursor].handle))) {
            
            // We're entering a chain, and this is a chain we want to flip. So
            // flip before emitting anything.

            // Flip: remember where we entered, jump to the other end,
            // reverse direction, emit the entry event there.
            flip_stack.push({cursor, reverse});
            cursor = other_bound[cursor];
            reverse = !reverse;
        }
        
        // Emit the event here
        emitter(reverse ? flip(events[cursor], graph) : events[cursor]);

        if (!flip_stack.empty() && cursor == flip_stack.top().entry_index) {
            // We've returned to the entry point of a flipped chain, so after
            // emitting, go back to the entry orientation and jump to the other
            // side, so we can advance out of it. 
            
            FlipEntry entry = flip_stack.top();
            flip_stack.pop();
            cursor = other_bound[entry.entry_index];
            reverse = entry.original_reverse;
        }
    }
}

// ReplaySnarlFinder implementation

ReplaySnarlFinder::ReplaySnarlFinder(const HandleGraph* graph, const std::vector<DecompositionEvent>& events) : HandleGraphSnarlFinder(graph) {
    this->events.reserve(events.size());
    for (const DecompositionEvent& e : events) {
        // Translate input events into handles
        this->events.emplace_back(e.type, graph->get_handle(e.id, e.is_reverse));
    }
}

void ReplaySnarlFinder::traverse_decomposition(
    const std::function<void(handle_t)>& begin_chain,
    const std::function<void(handle_t)>& end_chain,
    const std::function<void(handle_t)>& begin_snarl,
    const std::function<void(handle_t)>& end_snarl) const
{
    auto emitter = event_emitter(begin_chain, end_chain, begin_snarl, end_snarl);
    for (auto& event : events) {
        emitter(event);
    }
}

std::function<void(const DecompositionHandleEvent&)> event_emitter(
    const std::function<void(handle_t)>& begin_chain,
    const std::function<void(handle_t)>& end_chain,
    const std::function<void(handle_t)>& begin_snarl,
    const std::function<void(handle_t)>& end_snarl
) {
    return [&](const DecompositionHandleEvent& event) { 
        switch (event.type) {
        case ET::BEGIN_CHAIN:
            begin_chain(event.handle);
            break;
        case ET::END_CHAIN:
            end_chain(event.handle);
            break;
        case ET::BEGIN_SNARL:
            begin_snarl(event.handle);
            break;
        case ET::END_SNARL:
            end_snarl(event.handle);
            break;
        }
    };
}

std::vector<DecompositionEvent> capture_events(const HandleGraphSnarlFinder& finder, const HandleGraph& graph) {
    // Get all the events in terms of handles
    std::vector<DecompositionHandleEvent> handle_result = capture_events(finder);
    // And translate them to IDs and orientations
    std::vector<DecompositionEvent> result;
    result.reserve(handle_result.size());
    for (DecompositionHandleEvent& e : handle_result) {
        result.emplace_back(e.type, graph.get_id(e.handle), graph.get_is_reverse(e.handle));
    }
    return result;
}

std::vector<DecompositionHandleEvent> capture_events(const HandleGraphSnarlFinder& finder) {
    std::vector<DecompositionHandleEvent> result;
    // Mint out functions that push events of different types.
    auto event_pusher = [&result](ET event) {
        return [event,&result](const handle_t& h) {
            result.push_back({event, h});
        };
    };
    finder.traverse_decomposition(
        event_pusher(ET::BEGIN_CHAIN),
        event_pusher(ET::END_CHAIN),
        event_pusher(ET::BEGIN_SNARL),
        event_pusher(ET::END_SNARL)
    );
    return result;
}

} // namespace unittest
} // namespace vg

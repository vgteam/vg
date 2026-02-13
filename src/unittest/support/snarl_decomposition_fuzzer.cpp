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
    // This vector is indexable with the event types.
    std::vector<const std::function<void(handle_t)>*> handlers {&begin_chain, &end_chain, &begin_snarl, &end_snarl};

    // Make a helper to emit an event, transforming it based on direction.
    // Forward: emit as-is.
    // Backward: swap begin/end types and flip handles.
    std::function<void(const DecompositionEvent&, bool)> emit_event = [&](const DecompositionEvent& event, bool reverse) {
        if (reverse) {
            // Flip the event around if needed
            emit_event(flip(event, graph), false);
        } else {
            // Call the right handler on the event's handle.
            (*handlers.at((int)event.type))(event.handle);
        }
    };

    // Step 1: Capture all events from the wrapped finder.
    std::vector<DecompositionEvent> events;
    wrapped->traverse_decomposition(
        [&](handle_t h) { events.push_back({ET::BEGIN_CHAIN, h}); },
        [&](handle_t h) { events.push_back({ET::END_CHAIN, h}); },
        [&](handle_t h) { events.push_back({ET::BEGIN_SNARL, h}); },
        [&](handle_t h) { events.push_back({ET::END_SNARL, h}); }
    );

    if (events.empty()) {
        return;
    }

    // Step 2: Build pairing vector mapping each begin to its matching end
    // and vice versa, using separate stacks for chains and snarls.
    std::vector<size_t> pair_of(events.size());
    {
        stack<size_t> chain_stack, snarl_stack;
        for (size_t i = 0; i < events.size(); i++) {
            switch (events[i].type) {
            case ET::BEGIN_CHAIN:
                chain_stack.push(i);
                break;
            case ET::END_CHAIN:
                assert(!chain_stack.empty());
                pair_of[i] = chain_stack.top();
                pair_of[chain_stack.top()] = i;
                chain_stack.pop();
                break;
            case ET::BEGIN_SNARL:
                snarl_stack.push(i);
                break;
            case ET::END_SNARL:
                assert(!snarl_stack.empty());
                pair_of[i] = snarl_stack.top();
                pair_of[snarl_stack.top()] = i;
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
            cursor = pair_of[cursor];
            reverse = !reverse;
        }
        
        // Emit the event here
        emit_event(events[cursor], reverse);

        if (!flip_stack.empty() && cursor == flip_stack.top().entry_index) {
            // We've returned to the entry point of a flipped chain, so after
            // emitting, go back to the entry orientation and jump to the other
            // side, so we can advance out of it. 
            
            FlipEntry entry = flip_stack.top();
            flip_stack.pop();
            cursor = pair_of[entry.entry_index];
            reverse = entry.original_reverse;
        }
    }
}

// ReplaySnarlFinder implementation

ReplaySnarlFinder::ReplaySnarlFinder(const std::vector<ReplaySnarlFinder::Event>& events)
    : HandleGraphSnarlFinder(nullptr), events(events)
{
}

void ReplaySnarlFinder::traverse_decomposition(
    const std::function<void(handle_t)>& begin_chain,
    const std::function<void(handle_t)>& end_chain,
    const std::function<void(handle_t)>& begin_snarl,
    const std::function<void(handle_t)>& end_snarl) const
{
    for (const auto& event : events) {
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
    }
}

} // namespace unittest
} // namespace vg

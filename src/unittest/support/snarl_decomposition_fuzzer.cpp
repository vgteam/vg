#include "snarl_decomposition_fuzzer.hpp"

#include <cassert>
#include <stack>

namespace vg {
namespace unittest {

using namespace std;
using ET = DecompositionEventType;

// SnarlDecompositionFuzzer deterministic constructor
SnarlDecompositionFuzzer::SnarlDecompositionFuzzer(
    const HandleGraph* graph,
    const HandleGraphSnarlFinder* finder,
    const set<pair<handle_t, handle_t>>& chains_to_flip)
    : HandleGraphSnarlFinder(graph), wrapped(finder)
{
    should_flip = [chains_to_flip](handle_t begin, handle_t end) -> bool {
        return chains_to_flip.count({begin, end}) > 0;
    };
}

void SnarlDecompositionFuzzer::emit_event(
    const DecompositionEvent& event,
    bool forward,
    const function<void(handle_t)>& begin_chain,
    const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl,
    const function<void(handle_t)>& end_snarl) const
{
    if (forward) {
        switch (event.type) {
        case ET::BEGIN_CHAIN: begin_chain(event.handle); break;
        case ET::END_CHAIN:   end_chain(event.handle); break;
        case ET::BEGIN_SNARL: begin_snarl(event.handle); break;
        case ET::END_SNARL:   end_snarl(event.handle); break;
        }
    } else {
        handle_t flipped = graph->flip(event.handle);
        switch (event.type) {
        case ET::BEGIN_CHAIN: end_chain(flipped); break;
        case ET::END_CHAIN:   begin_chain(flipped); break;
        case ET::BEGIN_SNARL: end_snarl(flipped); break;
        case ET::END_SNARL:   begin_snarl(flipped); break;
        }
    }
}

void SnarlDecompositionFuzzer::traverse_decomposition(
    const function<void(handle_t)>& begin_chain,
    const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl,
    const function<void(handle_t)>& end_snarl) const
{
    // Step 1: Capture all events from the wrapped finder.
    vector<DecompositionEvent> events;
    wrapped->traverse_decomposition(
        [&](handle_t h) { events.push_back({ET::BEGIN_CHAIN, h}); },
        [&](handle_t h) { events.push_back({ET::END_CHAIN, h}); },
        [&](handle_t h) { events.push_back({ET::BEGIN_SNARL, h}); },
        [&](handle_t h) { events.push_back({ET::END_SNARL, h}); }
    );

    if (events.empty()) return;

    // Step 2: Build pairing vector mapping each begin to its matching end
    // and vice versa, using separate stacks for chains and snarls.
    size_t n = events.size();
    vector<size_t> pair_of(n);
    {
        stack<size_t> chain_stack, snarl_stack;
        for (size_t i = 0; i < n; i++) {
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
        bool original_forward;
    };
    stack<FlipEntry> flip_stack;

    int64_t cursor = 0;
    bool forward = true;

    while (cursor >= 0 && cursor < (int64_t)n) {
        size_t idx = (size_t)cursor;

        // Check if we've returned to the entry point of a flipped chain.
        if (!flip_stack.empty() && idx == flip_stack.top().entry_index) {
            emit_event(events[idx], forward,
                       begin_chain, end_chain, begin_snarl, end_snarl);
            FlipEntry entry = flip_stack.top();
            flip_stack.pop();
            cursor = (int64_t)pair_of[entry.entry_index];
            forward = entry.original_forward;
            cursor += forward ? 1 : -1;
            continue;
        }

        // Check if we're entering a chain.
        bool is_chain_entry =
            (forward && events[idx].type == ET::BEGIN_CHAIN) ||
            (!forward && events[idx].type == ET::END_CHAIN);

        if (is_chain_entry) {
            size_t begin_idx = forward ? idx : pair_of[idx];
            size_t end_idx = forward ? pair_of[idx] : idx;
            handle_t begin_handle = events[begin_idx].handle;
            handle_t end_handle = events[end_idx].handle;

            if (should_flip(begin_handle, end_handle)) {
                // Flip: remember where we entered, jump to the other end,
                // reverse direction, emit the entry event there.
                flip_stack.push({idx, forward});
                cursor = (int64_t)pair_of[idx];
                forward = !forward;
                emit_event(events[(size_t)cursor], forward,
                           begin_chain, end_chain, begin_snarl, end_snarl);
                cursor += forward ? 1 : -1;
                continue;
            }
        }

        // Normal event: emit and advance.
        emit_event(events[idx], forward,
                   begin_chain, end_chain, begin_snarl, end_snarl);
        cursor += forward ? 1 : -1;
    }
}

// ReplaySnarlFinder implementation

ReplaySnarlFinder::ReplaySnarlFinder(const vector<ReplaySnarlFinder::Event>& events)
    : HandleGraphSnarlFinder(nullptr), events(events)
{
}

void ReplaySnarlFinder::traverse_decomposition(
    const function<void(handle_t)>& begin_chain,
    const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl,
    const function<void(handle_t)>& end_snarl) const
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

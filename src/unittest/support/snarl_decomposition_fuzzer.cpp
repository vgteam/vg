#include "snarl_decomposition_fuzzer.hpp"

#include <cassert>

namespace vg {
namespace unittest {

using namespace std;

// An event captured from the wrapped finder's decomposition
struct CapturedEvent {
    enum Type { BEGIN_CHAIN, END_CHAIN, BEGIN_SNARL, END_SNARL };
    Type type;
    handle_t handle;
};

/// Find the matching end event for a begin event at position `start`.
/// For BEGIN_CHAIN finds matching END_CHAIN; for BEGIN_SNARL finds END_SNARL.
/// Handles nesting correctly.
static size_t find_matching_end(const vector<CapturedEvent>& events, size_t start) {
    bool is_chain = (events[start].type == CapturedEvent::BEGIN_CHAIN);
    CapturedEvent::Type begin_type = is_chain ? CapturedEvent::BEGIN_CHAIN : CapturedEvent::BEGIN_SNARL;
    CapturedEvent::Type end_type = is_chain ? CapturedEvent::END_CHAIN : CapturedEvent::END_SNARL;

    int depth = 0;
    for (size_t i = start; i < events.size(); i++) {
        if (events[i].type == begin_type) {
            depth++;
        } else if (events[i].type == end_type) {
            depth--;
            if (depth == 0) {
                return i;
            }
        }
    }
    assert(false);
    return events.size();
}

// Forward declaration.
// do_flip: whether this chain should be emitted reversed. A chain's effective
// flip is the XOR of all flip decisions from itself and its ancestors:
// flipping a parent flips all its descendants, and if a descendant is also
// in the flip set it gets flipped again (canceling the parent's flip).
static void process_chain(
    const vector<CapturedEvent>& events,
    size_t chain_start, size_t chain_end,
    const HandleGraph& graph,
    const function<bool(handle_t, handle_t)>& should_flip,
    bool do_flip,
    const function<void(handle_t)>& begin_chain,
    const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl,
    const function<void(handle_t)>& end_snarl);

/// Process events in the range [start, end), recursively applying chain flipping.
/// This handles top-level chains in the root snarl.
static void process_events(
    const vector<CapturedEvent>& events,
    size_t start, size_t end,
    const HandleGraph& graph,
    const function<bool(handle_t, handle_t)>& should_flip,
    const function<void(handle_t)>& begin_chain,
    const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl,
    const function<void(handle_t)>& end_snarl)
{
    size_t i = start;
    while (i < end) {
        if (events[i].type == CapturedEvent::BEGIN_CHAIN) {
            size_t chain_end_idx = find_matching_end(events, i);
            handle_t orig_begin = events[i].handle;
            handle_t orig_end = events[chain_end_idx].handle;
            bool do_flip = should_flip(orig_begin, orig_end);
            process_chain(events, i, chain_end_idx, graph, should_flip,
                         do_flip, begin_chain, end_chain, begin_snarl, end_snarl);
            i = chain_end_idx + 1;
        } else {
            switch (events[i].type) {
            case CapturedEvent::BEGIN_SNARL:
                begin_snarl(events[i].handle);
                break;
            case CapturedEvent::END_SNARL:
                end_snarl(events[i].handle);
                break;
            default:
                break;
            }
            i++;
        }
    }
}

/// Collect the index spans of child chains within a snarl's interior
/// (between interior_start and interior_end, exclusive of snarl boundaries).
static vector<pair<size_t, size_t>> collect_child_chains(
    const vector<CapturedEvent>& events,
    size_t interior_start, size_t interior_end)
{
    vector<pair<size_t, size_t>> children;
    size_t i = interior_start;
    while (i < interior_end) {
        if (events[i].type == CapturedEvent::BEGIN_CHAIN) {
            size_t chain_end_idx = find_matching_end(events, i);
            children.push_back({i, chain_end_idx});
            i = chain_end_idx + 1;
        } else {
            i++;
        }
    }
    return children;
}

/// Compute the effective flip for a child chain: the XOR of the parent's
/// effective flip and the child's own should_flip decision (always evaluated
/// on the child's original handles from the event stream).
static bool child_effective_flip(
    bool parent_do_flip,
    const function<bool(handle_t, handle_t)>& should_flip,
    handle_t child_orig_begin, handle_t child_orig_end)
{
    return parent_do_flip != should_flip(child_orig_begin, child_orig_end);
}

/// Process a single chain (events[chain_start] to events[chain_end], inclusive).
/// do_flip indicates whether this chain should be emitted reversed.
static void process_chain(
    const vector<CapturedEvent>& events,
    size_t chain_start, size_t chain_end,
    const HandleGraph& graph,
    const function<bool(handle_t, handle_t)>& should_flip,
    bool do_flip,
    const function<void(handle_t)>& begin_chain,
    const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl,
    const function<void(handle_t)>& end_snarl)
{
    handle_t orig_begin = events[chain_start].handle;
    handle_t orig_end = events[chain_end].handle;

    if (do_flip) {
        // Flip this chain: reverse snarl order, flip all boundary handles,
        // and recursively process children with their effective flip state.

        // Collect the snarls inside this chain.
        struct SnarlSpan {
            size_t begin_idx;
            size_t end_idx;
        };
        vector<SnarlSpan> snarls;

        size_t i = chain_start + 1;
        while (i < chain_end) {
            if (events[i].type == CapturedEvent::BEGIN_SNARL) {
                size_t snarl_end = find_matching_end(events, i);
                snarls.push_back({i, snarl_end});
                i = snarl_end + 1;
            } else {
                i++;
            }
        }

        begin_chain(graph.flip(orig_end));

        // Emit snarls in reverse order with flipped boundaries
        for (int s = (int)snarls.size() - 1; s >= 0; s--) {
            const auto& snarl = snarls[s];
            handle_t snarl_begin_h = events[snarl.begin_idx].handle;
            handle_t snarl_end_h = events[snarl.end_idx].handle;

            begin_snarl(graph.flip(snarl_end_h));

            // Collect and emit child chains in reverse order
            auto children = collect_child_chains(events, snarl.begin_idx + 1, snarl.end_idx);
            for (int c = (int)children.size() - 1; c >= 0; c--) {
                handle_t child_begin = events[children[c].first].handle;
                handle_t child_end = events[children[c].second].handle;
                bool child_flip = child_effective_flip(do_flip, should_flip, child_begin, child_end);
                process_chain(events, children[c].first, children[c].second,
                             graph, should_flip, child_flip,
                             begin_chain, end_chain, begin_snarl, end_snarl);
            }

            end_snarl(graph.flip(snarl_begin_h));
        }

        end_chain(graph.flip(orig_begin));
    } else {
        // Don't flip this chain, but still recursively process nested chains
        // which may have their own flip decisions.
        begin_chain(orig_begin);

        size_t i = chain_start + 1;
        while (i < chain_end) {
            if (events[i].type == CapturedEvent::BEGIN_SNARL) {
                size_t snarl_end = find_matching_end(events, i);
                begin_snarl(events[i].handle);

                auto children = collect_child_chains(events, i + 1, snarl_end);
                for (auto& child : children) {
                    handle_t child_begin = events[child.first].handle;
                    handle_t child_end = events[child.second].handle;
                    bool child_flip = child_effective_flip(do_flip, should_flip, child_begin, child_end);
                    process_chain(events, child.first, child.second,
                                 graph, should_flip, child_flip,
                                 begin_chain, end_chain, begin_snarl, end_snarl);
                }

                end_snarl(events[snarl_end].handle);
                i = snarl_end + 1;
            } else {
                i++;
            }
        }

        end_chain(orig_end);
    }
}

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

void SnarlDecompositionFuzzer::traverse_decomposition(
    const function<void(handle_t)>& begin_chain,
    const function<void(handle_t)>& end_chain,
    const function<void(handle_t)>& begin_snarl,
    const function<void(handle_t)>& end_snarl) const
{
    // Step 1: Capture all events from the wrapped finder
    vector<CapturedEvent> events;
    wrapped->traverse_decomposition(
        [&](handle_t h) { events.push_back({CapturedEvent::BEGIN_CHAIN, h}); },
        [&](handle_t h) { events.push_back({CapturedEvent::END_CHAIN, h}); },
        [&](handle_t h) { events.push_back({CapturedEvent::BEGIN_SNARL, h}); },
        [&](handle_t h) { events.push_back({CapturedEvent::END_SNARL, h}); }
    );

    // Step 2: Process events, flipping chains as needed
    process_events(events, 0, events.size(), *graph, should_flip,
                  begin_chain, end_chain, begin_snarl, end_snarl);
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
        case EventType::BEGIN_CHAIN:
            begin_chain(event.handle);
            break;
        case EventType::END_CHAIN:
            end_chain(event.handle);
            break;
        case EventType::BEGIN_SNARL:
            begin_snarl(event.handle);
            break;
        case EventType::END_SNARL:
            end_snarl(event.handle);
            break;
        }
    }
}

} // namespace unittest
} // namespace vg

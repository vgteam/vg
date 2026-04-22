#ifndef VG_UNITTEST_SNARL_DECOMPOSITION_FUZZER_HPP_INCLUDED
#define VG_UNITTEST_SNARL_DECOMPOSITION_FUZZER_HPP_INCLUDED

/**
 * \file snarl_decomposition_fuzzer.hpp
 * Provides SnarlDecompositionFuzzer, which wraps a HandleGraphSnarlFinder and
 * randomly flips chains in the snarl decomposition, and ReplaySnarlFinder,
 * which replays a scripted sequence of decomposition events.
 */

#include <functional>
#include <random>
#include <vector>
#include <set>
#include <utility>
#include "snarls.hpp"
#include "handle.hpp"

namespace vg {
namespace unittest {

/// Event types for snarl decomposition traversal.
enum class DecompositionEventType {
    BEGIN_CHAIN = 0,
    END_CHAIN,
    BEGIN_SNARL,
    END_SNARL
};

inline std::ostream& operator<<(std::ostream& out, const DecompositionEventType& t) {
    int bits = (int)t;
    return out << (bits & 1 ? "END" : "BEGIN") << "_" << (bits & 2 ? "SNARL" : "CHAIN");
}

/// Flip the polatiry of an event type (start vs. end)
inline DecompositionEventType flip(const DecompositionEventType& t) {
    // We can flip by toggling the low bit.
    return (DecompositionEventType)((int) t ^ 1);
}

/// A single event in a snarl decomposition traversal.
/// This is in terms of IDs and orientations because those are easier to write in test code.
struct DecompositionEvent {
    DecompositionEventType type;
    nid_t id;
    bool is_reverse;

    inline bool operator==(const DecompositionEvent& other) const {
        return type == other.type && id == other.id && is_reverse == other.is_reverse;
    }

    inline bool operator!=(const DecompositionEvent& other) const {
        return ! (*this == other);
    }
};

inline std::ostream& operator<<(std::ostream& out, const DecompositionEvent& e) {
    return out << e.type << "(" << e.id << (e.is_reverse ? "-" : "+") << ")";
}

/// A single event in a snarl decomposition traversal.
/// This is in terms of handles because those are easier to work with internally.
struct DecompositionHandleEvent {
    DecompositionEventType type;
    handle_t handle;
};

/// Flip the polarity of a whole event (event type between begin and end, and handle orientation)
inline DecompositionHandleEvent flip(const DecompositionHandleEvent& e, const HandleGraph* g) {
    return {flip(e.type), g->flip(e.handle)};
}

/// Turn begin and end functions to call into a function that emits an event by
/// type. The provided functions must outlive the returned function.
std::function<void(const DecompositionHandleEvent&)> event_emitter(
    const std::function<void(handle_t)>& begin_chain,
    const std::function<void(handle_t)>& end_chain,
    const std::function<void(handle_t)>& begin_snarl,
    const std::function<void(handle_t)>& end_snarl
);

/// Capture all events emitted by a snarl finder, in terms of IDs and orientations.
std::vector<DecompositionEvent> capture_events(const HandleGraphSnarlFinder& finder, const HandleGraph& graph);

/// Capture all events emitted by a snarl finder, in terms of handles.
std::vector<DecompositionHandleEvent> capture_events(const HandleGraphSnarlFinder& finder);

/**
 * A HandleGraphSnarlFinder that wraps another HandleGraphSnarlFinder and
 * randomly flips chains in the snarl decomposition. Flipping a chain reverses
 * the entire chain including all children; if a child chain is also selected
 * for flipping, it gets flipped again (canceling the parent's flip for that
 * child).
 *
 * For deterministic testing, chains_to_flip can be provided, which is a set
 * of (begin_handle, end_handle) pairs identifying chains to flip. When
 * provided, p_flip and the generator are ignored.
 */
class SnarlDecompositionFuzzer : public HandleGraphSnarlFinder {
public:
    /**
     * Construct a fuzzer wrapping the given finder, flipping chains with
     * probability p_flip using the given random generator.
     * The graph pointer is needed to flip handles.
     */
    template<typename URNG>
    SnarlDecompositionFuzzer(const HandleGraph* graph,
                             const HandleGraphSnarlFinder* finder,
                             double p_flip, URNG& generator);

    /**
     * Construct a fuzzer wrapping the given finder, flipping the chains
     * bounded by the given node IDs.
     *
     * You should provide both bounding IDs for each chain, but only the one
     * that the chain is actually arrived at through during the traversal will
     * really get used.
     *
     * Note that a node can bound at most one chain.
     *
     * This is mostly for testing the fuzzer itself.
     */
    SnarlDecompositionFuzzer(const HandleGraph* graph,
                             const HandleGraphSnarlFinder* finder,
                             const std::unordered_set<nid_t>& chains_to_flip);

    virtual ~SnarlDecompositionFuzzer() = default;

    /**
     * Traverse the snarl decomposition, flipping selected chains.
     */
    virtual void traverse_decomposition(
        const std::function<void(handle_t)>& begin_chain,
        const std::function<void(handle_t)>& end_chain,
        const std::function<void(handle_t)>& begin_snarl,
        const std::function<void(handle_t)>& end_snarl
    ) const override;

private:
    /// The wrapped snarl finder
    const HandleGraphSnarlFinder* wrapped;

    /// Function that decides whether to flip a chain, given either of its
    /// bounding node IDs. May be nondeterministic.
    std::function<bool(nid_t)> should_flip;
};

/**
 * A HandleGraphSnarlFinder that replays a scripted sequence of decomposition
 * events. Useful for testing SnarlDecompositionFuzzer without needing a real
 * graph or snarl finder.
 */
class ReplaySnarlFinder : public HandleGraphSnarlFinder {
public:
    /**
     * Construct a replay finder that will emit the given events.
     */
    ReplaySnarlFinder(const HandleGraph* graph, const std::vector<DecompositionEvent>& events);

    virtual ~ReplaySnarlFinder() = default;

    /**
     * Replay the scripted events.
     */
    virtual void traverse_decomposition(
        const std::function<void(handle_t)>& begin_chain,
        const std::function<void(handle_t)>& end_chain,
        const std::function<void(handle_t)>& begin_snarl,
        const std::function<void(handle_t)>& end_snarl
    ) const override;

private:

    using EventType = DecompositionEventType;
    using Event = DecompositionHandleEvent;
    
    /// This stores events we are going to replay.
    std::vector<Event> events;
};

// Template implementation

template<typename URNG>
SnarlDecompositionFuzzer::SnarlDecompositionFuzzer(
    const HandleGraph* graph,
    const HandleGraphSnarlFinder* finder,
    double p_flip, URNG& generator)
    : HandleGraphSnarlFinder(graph), wrapped(finder)
{
    should_flip = [&generator, p_flip](nid_t ignored) -> bool {
        return std::uniform_real_distribution<double>(0.0, 1.0)(generator) < p_flip;
    };
}

} // namespace unittest
} // namespace vg

#endif

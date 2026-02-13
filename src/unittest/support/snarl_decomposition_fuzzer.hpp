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
    BEGIN_CHAIN,
    END_CHAIN,
    BEGIN_SNARL,
    END_SNARL
};

/// A single event in a snarl decomposition traversal.
struct DecompositionEvent {
    DecompositionEventType type;
    handle_t handle;
};

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
     * Construct a fuzzer wrapping the given finder, flipping exactly the
     * chains identified by the given set of (begin_handle, end_handle) pairs.
     * The handles should be the inward-facing begin and outward-facing end
     * handles as originally emitted by the wrapped finder.
     */
    SnarlDecompositionFuzzer(const HandleGraph* graph,
                             const HandleGraphSnarlFinder* finder,
                             const std::set<std::pair<handle_t, handle_t>>& chains_to_flip);

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

    /// Function that decides whether to flip a chain given its begin and end handles
    std::function<bool(handle_t, handle_t)> should_flip;
};

/**
 * A HandleGraphSnarlFinder that replays a scripted sequence of decomposition
 * events. Useful for testing SnarlDecompositionFuzzer without needing a real
 * graph or snarl finder.
 */
class ReplaySnarlFinder : public HandleGraphSnarlFinder {
public:
    /// Alias for shared event types
    using EventType = DecompositionEventType;
    using Event = DecompositionEvent;

    /**
     * Construct a replay finder that will emit the given events.
     * The graph pointer can be null since we never actually use it.
     */
    ReplaySnarlFinder(const std::vector<Event>& events);

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
    auto gen = std::make_shared<URNG>(generator);
    auto dist = std::make_shared<std::uniform_real_distribution<double>>(0.0, 1.0);
    should_flip = [gen, dist, p_flip](handle_t, handle_t) -> bool {
        return (*dist)(*gen) < p_flip;
    };
}

} // namespace unittest
} // namespace vg

#endif

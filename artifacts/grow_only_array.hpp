#ifndef VG_GROW_ONLY_ARRAY_HPP_INCLUDED
#define VG_GROW_ONLY_ARRAY_HPP_INCLUDED

/** \file
 * An array that threads can add to and read from at once, without locks.
 */

#include <array>
#include <atomic>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace vg {

/**
 * An array that only grows, which threads can add to and read from at once, without locks.
 *
 * Each element gets an index when it is added, and keeps it and its address until the array is
 * destroyed. Elements cannot be changed or removed. Reading an element by index takes constant
 * time.
 *
 * The elements live in blocks that are never moved. Block `b` holds `2^(FirstBlockBits + b)`
 * elements, so the blocks double in size, and few are needed. Adding an element takes the next
 * index from an atomic counter, so threads adding at once each get their own. The first thread
 * to need a block publishes it with a compare-and-swap, and the others use that block.
 *
 * A thread may read an element once the `emplace_back()` that added it has returned, and the
 * index has reached the reading thread through something that orders the two threads, such as a
 * mutex or the thread that added it handing it over. Reading an index that was never returned is
 * undefined.
 *
 * @tparam T The element type. It must be constructible from the arguments given to
 * `emplace_back()` without throwing, so that every index handed out holds an element.
 * @tparam FirstBlockBits The base-2 logarithm of the number of elements in the first block.
 */
template<typename T, size_t FirstBlockBits = 6>
class GrowOnlyArray {
public:
    GrowOnlyArray() = default;

    /// Destroy the elements. No thread may be using the array.
    ~GrowOnlyArray() {
        const size_t count = next_index.load();
        for (size_t i = 0; i < count; ++i) {
            const auto [block, offset] = locate(i);
            std::destroy_at(blocks[block].load(std::memory_order_relaxed) + offset);
        }
        std::allocator<T> allocator;
        for (size_t block = 0; block < BLOCK_COUNT; ++block) {
            if (T* storage = blocks[block].load(std::memory_order_relaxed); storage != nullptr) {
                allocator.deallocate(storage, block_size(block));
            }
        }
    }

    // Threads hold the addresses of elements, so the array cannot move.
    GrowOnlyArray(const GrowOnlyArray&) = delete;
    GrowOnlyArray& operator=(const GrowOnlyArray&) = delete;

    /**
     * Add an element, made from `args`, and return its index.
     *
     * Safe to call from several threads at once, and while other threads read.
     */
    template<typename... Args>
    size_t emplace_back(Args&&... args) {
        static_assert(std::is_nothrow_constructible_v<T, Args...>,
                      "GrowOnlyArray elements must be constructible without throwing");
        const size_t index = next_index.fetch_add(1);
        if (index > SIZE_MAX - block_size(0)) {
            throw std::length_error("GrowOnlyArray: too many elements");
        }
        const auto [block, offset] = locate(index);
        T* storage = blocks[block].load(std::memory_order_acquire);
        if (storage == nullptr) {
            // Make the block, and publish it unless another thread has already.
            std::allocator<T> allocator;
            T* fresh = allocator.allocate(block_size(block));
            if (blocks[block].compare_exchange_strong(storage, fresh, std::memory_order_acq_rel,
                                                      std::memory_order_acquire)) {
                storage = fresh;
            } else {
                // `storage` now holds the other thread's block.
                allocator.deallocate(fresh, block_size(block));
            }
        }
        std::construct_at(storage + offset, std::forward<Args>(args)...);
        return index;
    }

    /// Get the element at `index`.
    const T& operator[](size_t index) const {
        const auto [block, offset] = locate(index);
        return blocks[block].load(std::memory_order_acquire)[offset];
    }

private:
    /// Enough blocks to hold every index a `size_t` can count.
    static constexpr size_t BLOCK_COUNT = sizeof(size_t) * 8 - FirstBlockBits;

    /// Get the number of elements block `block` holds.
    static constexpr size_t block_size(size_t block) {
        return size_t(1) << (FirstBlockBits + block);
    }

    /**
     * Get the block holding element `index`, and its offset in the block.
     *
     * Counting from `block_size(0)` instead of 0, block `b` starts at `2^(FirstBlockBits + b)`,
     * so the block is the position of the highest set bit, less `FirstBlockBits`.
     */
    static std::pair<size_t, size_t> locate(size_t index) noexcept {
        const size_t shifted = index + block_size(0);
        const size_t block = std::bit_width(shifted) - 1 - FirstBlockBits;
        return {block, shifted - block_size(block)};
    }

    /// The blocks, each null until an element in it is added.
    std::array<std::atomic<T*>, BLOCK_COUNT> blocks {};
    /// The index the next element added gets.
    std::atomic<size_t> next_index {0};
};

}

#endif

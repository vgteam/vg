#include "haplotype_minimizer_index.hpp"

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/utils.h>

#include <omp.h>

#include <algorithm>
#include <compare>
#include <cstdlib>
#include <deque>
#include <iostream>

namespace vg {

namespace {

using code_type = HaplotypeMinimizerIndex::code_type;
using minimizer_type = HaplotypeMinimizerIndex::minimizer_type;
using finder_type = gbwtgraph::MinimizerIndex<HaplotypeMinimizerIndex::key_type>;

/// The number of bits in a bucket number.
constexpr size_t BUCKET_BITS = 10;
constexpr size_t BUCKETS = size_t(1) << BUCKET_BITS;

/**
 * Return the bucket of an encoded key
 *
 * A minimizer is the kmer with the smallest hash in its window, so the top bits of the hash of a
 * minimizer key are mostly zero. The bucket is the top bits of a different hash, the finalizer of
 * MurmurHash3, which spreads minimizer keys evenly over the buckets.
 */
size_t bucket_of(code_type key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key >> (64 - BUCKET_BITS);
}

/// A minimizer key and an encoded gbwtgraph::Position where it occurs.
struct Occurrence {
    code_type key;
    code_type position;

    auto operator<=>(const Occurrence&) const = default;
};

/**
 * A thread's store of the minimizer occurrences it found, with one bin for each bucket
 *
 * Each thread adds to its own store at every occurrence, so the stores are aligned to cache lines,
 * which keeps one thread's writes from invalidating the store of another.
 *
 * Occurrences wait in a buffer until it fills. Consecutive windows of one thread overlap and share
 * many occurrences, so the buffer is sorted and its duplicates removed before the occurrences go to
 * their bins. A bin is a deque, which grows in small blocks without copying.
 */
class alignas(gbwtgraph::CACHE_LINE_SIZE) OccurrenceBins {
public:
    OccurrenceBins() : bins(BUCKETS) {}

    /// Add an occurrence.
    void add(Occurrence occurrence) {
        this->buffer.push_back(occurrence);
    }

    /// Move the buffered occurrences to their bins if the buffer is full.
    void flush_if_full() {
        if (this->buffer.size() >= BUFFER_SIZE) {
            this->flush();
        }
    }

    /// Move the buffered occurrences to their bins.
    void flush() {
        gbwt::removeDuplicates(this->buffer, false);
        for (const Occurrence& occurrence : this->buffer) {
            this->bins[bucket_of(occurrence.key)].push_back(occurrence);
        }
        this->buffer.clear();
    }

    /// Return the number of occurrences in the bin of the bucket.
    size_t bin_size(size_t bucket) const {
        return this->bins[bucket].size();
    }

    /// Append the occurrences in the bin of the bucket to the vector, and free the bin.
    void take(size_t bucket, std::vector<Occurrence>& out) {
        std::deque<Occurrence>& bin = this->bins[bucket];
        out.insert(out.end(), bin.begin(), bin.end());
        std::deque<Occurrence>().swap(bin);
    }

private:
    /// The number of occurrences at which the buffer is flushed.
    constexpr static size_t BUFFER_SIZE = 65536;

    std::vector<Occurrence> buffer;
    std::vector<std::deque<Occurrence>> bins;
};

/**
 * Add the occurrences of the minimizers in a haplotype window to the bins
 *
 * The position of a minimizer is the graph position of its first base, read in the orientation of
 * the minimizer. For a reverse complement, that base is the last one of the minimizer in the
 * window. We find the positions as gbwtgraph::index_haplotypes() does, so that both indexes hold
 * the same occurrences.
 *
 * @param traversal The handles that the window visits.
 * @param sequence The sequence of the window. It starts at the start of the first handle.
 */
void add_window_occurrences(
    const gbwtgraph::GBWTGraph& graph, const finder_type& finder,
    const std::vector<handlegraph::handle_t>& traversal, const std::string& sequence, OccurrenceBins& bins
) {
    auto iter = traversal.begin();
    size_t node_start = 0;
    for (const minimizer_type& minimizer : finder.minimizers(sequence)) {
        if (minimizer.empty()) {
            continue;
        }
        // Find the node that holds the base at the offset of the minimizer.
        size_t node_length = graph.get_length(*iter);
        while (node_start + node_length <= minimizer.offset) {
            node_start += node_length;
            ++iter;
            node_length = graph.get_length(*iter);
        }
        gbwtgraph::pos_t pos = gbwtgraph::make_pos_t(
            graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start
        );
        if (minimizer.is_reverse) {
            pos = gbwtgraph::reverse_base_pos(pos, node_length);
        }
        if (!gbwtgraph::Position::valid_offset(pos)) {
            #pragma omp critical (cerr)
            {
                std::cerr << "error: [HaplotypeMinimizerIndex] node offset " << gbwtgraph::offset(pos)
                          << " is too large" << std::endl;
            }
            std::exit(EXIT_FAILURE);
        }
        bins.add({ minimizer.key.get_key(), gbwtgraph::Position(pos).value });
    }
    bins.flush_if_full();
}

/**
 * Find the minimizer occurrences in all haplotype windows of the graph
 *
 * The windows are the haplotype-consistent walks that gbwtgraph::for_each_haplotype_window()
 * reports.
 *
 * @returns The occurrences found by each thread, indexed by OpenMP thread number.
 */
std::vector<OccurrenceBins> collect_occurrences(const gbwtgraph::GBWTGraph& graph, const finder_type& finder) {
    int threads = omp_get_max_threads();
    std::vector<OccurrenceBins> result(static_cast<size_t>(threads));
    gbwtgraph::for_each_haplotype_window(graph, finder.window_bp(),
        [&graph, &finder, &result](const std::vector<handlegraph::handle_t>& traversal, const std::string& sequence) {
            add_window_occurrences(graph, finder, traversal, sequence, result[omp_get_thread_num()]);
        }, threads > 1);
    for (OccurrenceBins& thread_bins : result) {
        thread_bins.flush();
    }
    return result;
}

/// Take the occurrences in a bucket from the bins of every thread. Return them sorted and distinct.
std::vector<Occurrence> take_bucket(std::vector<OccurrenceBins>& bins, size_t bucket) {
    size_t total = 0;
    for (const OccurrenceBins& thread_bins : bins) {
        total += thread_bins.bin_size(bucket);
    }
    std::vector<Occurrence> result;
    result.reserve(total);
    for (OccurrenceBins& thread_bins : bins) {
        thread_bins.take(bucket, result);
    }
    gbwt::removeDuplicates(result, false);
    return result;
}

/// Return the number of distinct keys among occurrences sorted by key.
size_t distinct_keys(const std::vector<Occurrence>& occurrences) {
    size_t result = 0;
    for (size_t i = 0; i < occurrences.size(); i++) {
        if (i == 0 || occurrences[i].key != occurrences[i - 1].key) {
            result++;
        }
    }
    return result;
}

} // anonymous namespace

//------------------------------------------------------------------------------

HaplotypeMinimizerIndex::HaplotypeMinimizerIndex(const gbwtgraph::GBZ& gbz, size_t k, size_t w) :
    finder(k, w, 0, false), buckets(BUCKETS) {
    std::vector<OccurrenceBins> bins = collect_occurrences(gbz.graph, this->finder);

    // A thread that builds a bucket reads only the bins of that bucket, so the buckets are built in
    // parallel without a lock.
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t bucket_id = 0; bucket_id < BUCKETS; bucket_id++) {
        std::vector<Occurrence> occurrences = take_bucket(bins, bucket_id);
        // The bucket is built in a local object and moved into place, so threads building
        // neighbouring buckets do not write to the same cache line.
        Bucket bucket;
        bucket.keys.reserve(distinct_keys(occurrences));
        bucket.positions.reserve(occurrences.size());
        for (const Occurrence& occurrence : occurrences) {
            bucket.append(occurrence.key, occurrence.position);
        }
        this->buckets[bucket_id] = std::move(bucket);
    }
}

std::span<const code_type> HaplotypeMinimizerIndex::Bucket::find(code_type key) const {
    auto iter = std::ranges::lower_bound(this->keys, key, {}, &Entry::key);
    if (iter == this->keys.end() || iter->key != key) {
        return {};
    }
    size_t start = (iter == this->keys.begin() ? 0 : (iter - 1)->end);
    return std::span<const code_type>(this->positions.data() + start, iter->end - start);
}

std::span<const code_type> HaplotypeMinimizerIndex::positions(code_type key) const {
    return this->buckets[bucket_of(key)].find(key);
}

size_t HaplotypeMinimizerIndex::size() const {
    size_t result = 0;
    for (const Bucket& bucket : this->buckets) {
        result += bucket.keys.size();
    }
    return result;
}

size_t HaplotypeMinimizerIndex::number_of_values() const {
    size_t result = 0;
    for (const Bucket& bucket : this->buckets) {
        result += bucket.positions.size();
    }
    return result;
}

size_t HaplotypeMinimizerIndex::unique_keys() const {
    size_t result = 0;
    for (const Bucket& bucket : this->buckets) {
        size_t start = 0;
        for (const Bucket::Entry& entry : bucket.keys) {
            if (entry.end - start == 1) {
                result++;
            }
            start = entry.end;
        }
    }
    return result;
}

}

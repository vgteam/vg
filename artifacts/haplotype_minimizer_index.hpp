#ifndef VG_HAPLOTYPE_MINIMIZER_INDEX_HPP_INCLUDED
#define VG_HAPLOTYPE_MINIMIZER_INDEX_HPP_INCLUDED

/**
 * \file haplotype_minimizer_index.hpp
 * An index of the minimizers in the haplotypes of a GBZ graph, built on all threads.
 */

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>

#include <compare>
#include <cstddef>
#include <span>
#include <string>
#include <vector>

namespace vg {

/**
 * A table from each minimizer in the haplotypes of a GBZ graph to the graph positions where it
 * occurs.
 *
 * It holds the keys and positions that gbwtgraph::index_haplotypes() puts into a
 * gbwtgraph::MinimizerIndex<Key64> with no payload, unweighted minimizers and no syncmers, and it
 * looks them up by key.
 *
 * The gbwtgraph index inserts every occurrence under one lock, so only one thread builds it at a
 * time. This index takes no lock. It divides the keys into buckets by a hash, and each thread
 * keeps the occurrences it finds in bins of its own, one for each bucket. Each bucket is then
 * gathered from the bins of all threads, sorted and cleared of duplicates on a thread of its own.
 *
 * The positions of each key are distinct and sorted. Each key has at least one position.
 */
class HaplotypeMinimizerIndex {
public:
    using key_type = gbwtgraph::Key64;
    using minimizer_type = gbwtgraph::Kmer<key_type>;
    using code_type = gbwtgraph::KmerEncoding::code_type;

    /**
     * Index the minimizers in all haplotypes of the graph
     *
     * Runs on all OpenMP threads.
     *
     * @param k The length of a kmer, at most 31.
     * @param w The number of consecutive kmers in a minimizer window.
     */
    HaplotypeMinimizerIndex(const gbwtgraph::GBZ& gbz, size_t k, size_t w);

    /// Return the length of a kmer.
    size_t k() const { return this->finder.k(); }

    /// Return the number of consecutive kmers in a minimizer window.
    size_t w() const { return this->finder.w(); }

    /// Return the minimizers in the sequence, sorted by offset.
    std::vector<minimizer_type> minimizers(const std::string& sequence) const {
        return this->finder.minimizers(sequence);
    }

    /// A minimizer key and a graph position where it occurs, encoded as a gbwtgraph::Position.
    struct Occurrence {
        code_type key;
        code_type position;

        auto operator<=>(const Occurrence&) const = default;
    };

    /// Return the number of graph positions where the kmer with the key occurs.
    size_t count(code_type key) const { return this->occurrences(key).size(); }

    /// Return the occurrences of the kmer with the key, sorted by position. The list is empty when
    /// the kmer does not occur.
    std::span<const Occurrence> occurrences(code_type key) const;

    /// Return the number of distinct keys.
    size_t size() const;

    /// Return the number of (key, position) occurrences.
    size_t number_of_values() const;

    /// Return the number of keys with a single position.
    size_t unique_keys() const;

private:
    /// Finds the minimizers in a sequence. It holds no keys, and only its parameters are used.
    gbwtgraph::MinimizerIndex<key_type> finder;

    /// The occurrences in each bucket, sorted and distinct. A bucket is found by the top bits of
    /// a hash of the key.
    std::vector<std::vector<Occurrence>> buckets;
};

}

#endif

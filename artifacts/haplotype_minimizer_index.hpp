#ifndef VG_HAPLOTYPE_MINIMIZER_INDEX_HPP_INCLUDED
#define VG_HAPLOTYPE_MINIMIZER_INDEX_HPP_INCLUDED

/**
 * \file haplotype_minimizer_index.hpp
 * An index of the minimizers in the haplotypes of a GBZ graph, built on all threads.
 */

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>

#include <cstddef>
#include <string>
#include <vector>

namespace vg {

/**
 * A table from each minimizer in the haplotypes of a GBZ graph to the graph positions where it
 * occurs.
 *
 * It holds the keys and positions that gbwtgraph::index_haplotypes() puts into a
 * gbwtgraph::MinimizerIndex<Key64> with no payload, unweighted minimizers and no syncmers. Its
 * lookup methods have the same names, types and answers as that index's, so it can replace that
 * index where nothing else is needed.
 *
 * The gbwtgraph index inserts every occurrence under one lock, so only one thread builds it at a
 * time. This index takes no lock. It divides the keys into buckets by a hash. Each thread keeps
 * the occurrences it finds in bins of its own, one for each bucket. Then each bucket is gathered
 * from the bins of all threads, sorted and cleared of duplicates on a thread of its own. Each
 * bucket is sized by the occurrences it holds. The graph is traversed twice, for half of the
 * buckets each time, so that the collected occurrences take no more memory than the finished
 * index.
 *
 * The positions of each key are distinct and sorted. Each key has at least one position.
 */
class HaplotypeMinimizerIndex {
public:
    using key_type = gbwtgraph::Key64;
    using minimizer_type = gbwtgraph::Kmer<key_type>;
    using code_type = gbwtgraph::KmerEncoding::code_type;
    using value_type = gbwtgraph::KmerEncoding::value_type;
    using multi_value_type = gbwtgraph::KmerEncoding::multi_value_type;

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

    /// Return the number of graph positions where the minimizer occurs.
    size_t count(const minimizer_type& minimizer) const { return this->find(minimizer).second; }

    /**
     * Return the graph positions where the minimizer occurs
     *
     * @returns A pointer to the sorted, encoded positions and their number. The pointer is null
     * when the minimizer does not occur. Decode the positions with get_value().
     */
    multi_value_type find(const minimizer_type& minimizer) const;

    /**
     * Return the position at the index in a list of positions from find()
     *
     * The payload pointer of the value is null, as this index has no payloads. An index past the end
     * of the list gives gbwtgraph::Position::no_pos().
     */
    value_type get_value(multi_value_type values, size_t index) const;

    /// Return the number of distinct keys.
    size_t size() const;

    /// Return the number of (key, position) occurrences.
    size_t number_of_values() const;

    /// Return the number of keys with a single position.
    size_t unique_keys() const;

private:
    /**
     * A table of the keys in one bucket, each with its positions
     *
     * `keys` is sorted by key. The positions of the key at rank i are
     * `positions[keys[i - 1].end, keys[i].end)`, starting at 0 for rank 0.
     */
    struct Bucket {
        /// A key and the end of its positions in `positions`.
        struct Entry {
            code_type key;
            size_t end;
        };

        std::vector<Entry> keys;
        /// Encoded gbwtgraph::Position values, sorted within each key.
        std::vector<code_type> positions;

        /// Add a position of the last key, or of a new last key. Keys and positions must come in
        /// sorted order, without duplicates.
        void append(code_type key, code_type position) {
            if (this->keys.empty() || this->keys.back().key != key) {
                this->keys.push_back({ key, 0 });
            }
            this->positions.push_back(position);
            this->keys.back().end = this->positions.size();
        }

        /// Return the positions of the key, or an empty list if the bucket does not hold it.
        multi_value_type find(code_type key) const;
    };

    /// Finds the minimizers in a sequence. It holds no keys, and only its parameters are used.
    gbwtgraph::MinimizerIndex<key_type> finder;

    /// The buckets, indexed by the top bits of a hash of the key.
    std::vector<Bucket> buckets;
};

}

#endif

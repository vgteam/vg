#include "diploid_pair_scorer.hpp"

#include <algorithm>
#include <bit>

namespace vg {

//------------------------------------------------------------------------------

DiploidPairScorer::DiploidPairScorer(
    const Haplotypes::Subchain& subchain, const std::vector<size_t>& sequences,
    const std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>>& kmer_types
) :
    words((kmer_types.size() + 63) / 64),
    classes(this->words),
    rows(sequences.size() * this->words, 0)
{
    for (size_t kmer_id = 0; kmer_id < kmer_types.size(); kmer_id++) {
        std::uint64_t bit = std::uint64_t(1) << (kmer_id % 64);
        ClassMasks& masks = this->classes[kmer_id / 64];
        switch (kmer_types[kmer_id].first) {
        case KmerPresenceMatrix::absent:
            masks.absent |= bit;
            break;
        case KmerPresenceMatrix::heterozygous:
            masks.heterozygous |= bit;
            break;
        case KmerPresenceMatrix::present:
            masks.present |= bit;
            break;
        default:
            break; // Frequent and ignored kmers score 0.
        }
    }

    std::vector<std::uint64_t> row;
    for (size_t i = 0; i < sequences.size(); i++) {
        subchain.presence_row(sequences[i], row);
        std::copy_n(row.begin(), std::min(row.size(), this->words), this->rows.begin() + i * this->words);
    }
}

std::int64_t DiploidPairScorer::score(size_t first, size_t second) const {
    const std::uint64_t* first_row = this->rows.data() + first * this->words;
    const std::uint64_t* second_row = this->rows.data() + second * this->words;
    std::int64_t result = 0;
    for (size_t i = 0; i < this->words; i++) {
        // The kmers with two, one and no copies in the pair. The bits past the
        // last kmer are in no class, so they do not count.
        std::uint64_t two = first_row[i] & second_row[i];
        std::uint64_t one = first_row[i] ^ second_row[i];
        std::uint64_t none = ~(first_row[i] | second_row[i]);
        const ClassMasks& masks = this->classes[i];
        result += std::popcount(masks.absent & none) - std::popcount(masks.absent & two);
        result += std::popcount(masks.present & two) - std::popcount(masks.present & none);
        result += std::popcount(masks.heterozygous & one);
    }
    return result;
}

//------------------------------------------------------------------------------

} // namespace vg

#ifndef VG_RYMER_INDEX_HPP_INCLUDED
#define VG_RYMER_INDEX_HPP_INCLUDED

/** \file
 * The RYmer index, which SAFARI uses to find seeds for ancient DNA reads.
 *
 * An RYmer is a k-mer written in a two-letter alphabet: purines (A and G) become R and pyrimidines
 * (C and T) become Y. Deamination, the commonest damage in ancient DNA, turns C into T and, on the
 * other strand, G into A. Neither change alters the RYmer, so a k-mer with damage still has the
 * RYmer of the k-mer it came from.
 *
 * The RYmer index is a minimizer index over the haplotypes of a GBZ graph. Minimizers are chosen
 * from the DNA sequence exactly as in an ordinary minimizer index with the same k and w. Each
 * chosen minimizer is then stored under the key of its RYmer, with its graph position, and with
 * a payload of one word that holds the DNA key of the minimizer. Looking up the RYmer of a read's
 * minimizer therefore gives the graph k-mers that the read's k-mer could be a damaged copy of.
 *
 * An RYmer key packs one bit per base, purine 0 and pyrimidine 1, with the first base in the most
 * significant bit used. These are the low bits of the two-bit DNA codes A=0, C=1, G=2, T=3.
 */

#include <string>
#include <vector>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>

namespace vg {

using namespace std;

/// The tag that marks a minimizer index as an RYmer index, and its value.
extern const string RYMER_INDEX_TAG; // "rymer_index"
extern const string RYMER_INDEX_TAG_VALUE; // "1"

/// The number of payload words in an RYmer index: one, for the DNA key of the minimizer.
constexpr size_t RYMER_PAYLOAD_SIZE = 1;

/// Returns the RYmer key of the k-mer whose DNA key is `dna_key`.
gbwtgraph::Key64 rymer_key(gbwtgraph::Key64 dna_key, size_t k);

/// Returns the minimizer that has the RYmer key of `dna_minimizer` in place of its DNA key, and
/// the hash of that key. Its offset and orientation are those of `dna_minimizer`.
gbwtgraph::DefaultMinimizerIndex::minimizer_type
to_rymer_minimizer(const gbwtgraph::DefaultMinimizerIndex::minimizer_type& dna_minimizer, size_t k);

/// Returns the RYmer of a DNA sequence, as a string of 'R' and 'Y', with 'N' for any character
/// that is not A, C, G or T in either case.
string rymer_string(const string& sequence);

/// Adds the minimizers of all haplotypes in `gbz` to `index` under their RYmer keys, as the file
/// comment describes, and tags the index as an RYmer index. The index must have a payload size of
/// `RYMER_PAYLOAD_SIZE` words and must not use syncmers.
void index_rymers(const gbwtgraph::GBZ& gbz, gbwtgraph::DefaultMinimizerIndex& index);

/// Returns true if the index has the tag that `index_rymers()` sets.
bool is_rymer_index(const gbwtgraph::DefaultMinimizerIndex& index);

/// Returns the DNA keys of the graph k-mers stored under the RYmer key of `rymer_minimizer`, in
/// the orientation in which the graph k-mers were chosen as minimizers. A key appears once for
/// each graph position it was found at.
vector<gbwtgraph::Key64> graph_kmers_for_rymer(const gbwtgraph::DefaultMinimizerIndex& index,
                                               const gbwtgraph::DefaultMinimizerIndex::minimizer_type& rymer_minimizer);

}

#endif

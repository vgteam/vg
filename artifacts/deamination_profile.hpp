#ifndef VG_DEAMINATION_PROFILE_HPP_INCLUDED
#define VG_DEAMINATION_PROFILE_HPP_INCLUDED

/** \file
 * Substitution rates along an ancient DNA fragment, and the test SAFARI uses to decide whether a
 * read k-mer could be a damaged copy of a graph k-mer.
 *
 * Ancient DNA molecules are damaged most near their ends. A profile gives, for each distance from
 * the 5' end and from the 3' end of a fragment, the probability of reading each base as each
 * other base. The profiles are 12-column `.prof` tables: a header line, then one row per
 * distance from the end, with the columns A>C A>G A>T C>A C>G C>T G>A G>C G>T T>A T>C T>G
 * separated by tabs. The first row is the outermost base. Positions beyond the last row take the
 * last row's rates.
 */

#include <array>
#include <string>
#include <vector>

namespace vg {

using namespace std;

/// The probability of reading each base as each other base at one position of a fragment,
/// indexed [true base][read base], each in the order A, C, G, T.
using SubstitutionMatrix = array<array<double, 4>, 4>;

class DeaminationProfile {
public:
    /// Reads the profile of the 5' end and the profile of the 3' end from two `.prof` files,
    /// which may be gzip-compressed. Throws std::runtime_error if a file cannot be read or is not
    /// a valid profile.
    DeaminationProfile(const string& five_prime_filename, const string& three_prime_filename);

    /// Makes a profile from substitution matrices for the positions counted from each end, the
    /// outermost first. Throws std::invalid_argument if either list is empty.
    DeaminationProfile(vector<SubstitutionMatrix> five_prime, vector<SubstitutionMatrix> three_prime);

    /// Returns the probability of reading `read_base` where the fragment has `true_base`, at
    /// 0-based `position` from the 5' end of a fragment `fragment_length` bases long. Where 'N'
    /// is either base the probability is 0.25. Other characters count as 'A'.
    ///
    /// For each true base this uses the substitution rates of whichever end, 5' or 3', is more
    /// likely to change that base at this position; that is, the end whose rates give the true
    /// base the lower probability of being read unchanged, the 5' end on a tie.
    double substitution_probability(size_t fragment_length, size_t position,
                                    char true_base, char read_base) const;

    /// Parses the text of a `.prof` file into its substitution matrices. `source` names the text
    /// in error messages. Throws std::runtime_error if the text is not a valid profile.
    static vector<SubstitutionMatrix> parse_profile(const string& text, const string& source);

private:
    /// Substitution rates by distance from the 5' end and from the 3' end, the outermost first.
    vector<SubstitutionMatrix> five_prime;
    vector<SubstitutionMatrix> three_prime;
};

/// Settings for `damage_posterior()`.
struct DamageModelParameters {
    /// The prior probability that a graph k-mer with the same RYmer as a read k-mer is not where
    /// the read comes from.
    double spurious_prior = 0.5;
};

/// Returns the posterior probability that `read_kmer`, which starts at 0-based `kmer_start` in a
/// read `read_length` bases long, is a copy of `graph_kmer` with damage, rather than a k-mer from
/// somewhere else that happens to share its RYmer. Both k-mers are in read orientation.
///
/// Under the damage model each read base comes from the graph base at the same offset with the
/// probability the profile gives at the position of the k-mer's first base, `kmer_start`; SAFARI
/// uses that one position for every base of the k-mer. Under the spurious model the number of
/// mismatches is binomial with a rate falling with the k-mer length k as 1.0014 * k^-0.6628, for
/// k-mers shorter than 22 bases, and has likelihood 0.01 for longer ones.
///
/// Returns 1 if the k-mers are equal and 0 if the damage model gives the read k-mer probability 0.
double damage_posterior(const string& graph_kmer, const string& read_kmer, size_t read_length,
                        size_t kmer_start, const DeaminationProfile& profile,
                        const DamageModelParameters& parameters);

}

#endif

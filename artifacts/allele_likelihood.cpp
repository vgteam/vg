#include "allele_likelihood.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "../src/alignment.hpp"
#include "off_panel.hpp"
#include "../src/path.hpp"
#include "../src/statistics.hpp"

namespace vg {

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// DepthContext
////////////////////////////////////////////////////////////////////////////////

double DepthContext::expected_reads(const vector<int>& genotype) const {
    double total = 0.0;
    for (int allele : genotype) {
        double len = (allele >= 0 && (size_t)allele < lengths.size())
                         ? (double)lengths[allele] : 0.0;
        total += len + read_length - 1.0;
    }
    return rate * total;
}

double DepthContext::ratio(const vector<int>& genotype) const {
    double expected = expected_reads(genotype);
    return expected > 0.0 ? observed / expected : -1.0;
}

////////////////////////////////////////////////////////////////////////////////
// AlleleReadLikelihoods
////////////////////////////////////////////////////////////////////////////////

void AlleleReadLikelihoods::set_contents(size_t reads, size_t alleles, vector<double>&& rows,
                                         vector<double>&& mismap, vector<double>&& best_ln,
                                         vector<string>&& names, size_t dropped) {
    this->n_reads = reads;
    this->n_alleles = alleles;
    this->matrix = std::move(rows);
    this->read_mismap_prob = std::move(mismap);
    this->read_best_ln = std::move(best_ln);
    this->read_names = std::move(names);
    this->unplaceable = dropped;
    this->effective_read_total = 0.0;
    for (double e : this->read_mismap_prob) {
        this->effective_read_total += 1.0 - e;
    }
    count_observed();
}

void AlleleReadLikelihoods::count_observed() {
    depth.observed = depth_effective ? effective_read_total : (double)n_reads;
}

double AlleleReadLikelihoods::rel(size_t r, size_t a) const {
    return matrix.at(r * n_alleles + a);
}

double AlleleReadLikelihoods::mismap_prob(size_t r) const {
    return read_mismap_prob.at(r);
}

double AlleleReadLikelihoods::best_ln_likelihood(size_t r) const {
    return read_best_ln.at(r);
}


/// Return n ln(lambda) - lambda - ln Gamma(n + 1), the log likelihood of lambda expected reads
/// when n reads are seen.
///
/// The observation n is `sum_r (1 - e_r)` rather than a row count, so it is usually not a whole
/// number, and ln Gamma(n + 1) stands in for ln n!. As a function of n >= 0 this is, up to a
/// normalising constant, the log density of a continuous analogue of the Poisson distribution.
/// The ln Gamma(n + 1) term does not depend on lambda, so it is the same for every genotype at a
/// site.
static double ln_poisson_pmf(double n, double lambda) {
    if (lambda <= 0.0) {
        return -numeric_limits<double>::infinity();
    }
    if (n <= 0.0) {
        return -lambda;
    }
    return n * log(lambda) - lambda - lgamma(n + 1.0);
}

/// Return ln Z_beta(lambda) by the normal approximation, where Z_beta(lambda) is the integral over
/// n >= 0 of f(n; lambda)^beta, and f is the density of `ln_poisson_pmf`.
///
/// For large lambda, f is close to a normal density with mean and variance lambda. So f^beta is
/// (2 pi lambda)^(-beta/2) times a normal density with variance lambda / beta, scaled by
/// sqrt(2 pi lambda / beta), and
///
///     ln Z_beta(lambda) ~= ((1 - beta) / 2) ln(2 pi lambda) - (1/2) ln beta.
///
/// Below lambda = 2 the approximation falls away from the integral, so we hold lambda at 2 there.
static double ln_tempered_poisson_normaliser(double lambda, double beta) {
    const double held = max(lambda, 2.0);
    return 0.5 * (1.0 - beta) * log(2.0 * M_PI * held) - 0.5 * log(beta);
}

double AlleleReadLikelihoods::weighted_length(const vector<int>& genotype, size_t i) const {
    const int allele = genotype[i];
    if (!unique_lengths.empty() && allele >= 0 && (size_t)allele < unique_lengths.size()) {
        // Take the sequence of this allele that no other member of the genotype carries, since
        // reads outside it cannot tell the two apart.
        size_t u = numeric_limits<size_t>::max();
        for (size_t j = 0; j < genotype.size(); ++j) {
            const int other = genotype[j];
            if (j == i || other < 0 || (size_t)other >= unique_lengths[allele].size()) {
                continue;
            }
            u = min(u, unique_lengths[allele][other]);
        }
        return (u == numeric_limits<size_t>::max()) ? 0.0 : (double)u;
    }
    if (allele >= 0 && (size_t)allele < allele_lengths.size()) {
        return (double)allele_lengths[allele];
    }
    return 0.0;
}

vector<double> AlleleReadLikelihoods::mixture_weights(const vector<int>& genotype) const {
    if (genotype.empty()) {
        return {};
    }
    double flat = 1.0 / (double)genotype.size();

    vector weights(genotype.size(), flat);
    if (uses_length_weights()) {
        double sum = 0.0;
        for (size_t i = 0; i < genotype.size(); ++i) {
            weights[i] = weighted_length(genotype, i) + mean_read_length - 1.0;
            sum += weights[i];
        }
        if (sum > 0.0) {
            for (double& w : weights) {
                w /= sum;
            }
        } else {
            weights.assign(genotype.size(), flat);
        }
    }
    return weights;
}

double AlleleReadLikelihoods::genotype_likelihood(const vector<int>& genotype) const {
    if (genotype.empty()) {
        return 0.0;
    }

    double total = 0.0;
    vector<double> weights = mixture_weights(genotype);

    for (size_t r = 0; r < n_reads; ++r) {
        // Marginalise over which haplotype of the genotype produced this read.
        double mixture = 0.0;
        for (size_t i = 0; i < genotype.size(); ++i) {
            int allele = genotype[i];
            // Skip the negative markers that the VCF layer uses for star and missing alleles.
            if (allele < 0 || (size_t)allele >= n_alleles) {
                continue;
            }
            mixture += weights[i] * rel(r, (size_t)allele);
        }

        // Mix in the chance that the read is mismapped, in which case its likelihood is the row
        // maximum, 1. The result lies in [e_r, 1], so its log is finite.
        double e_r = read_mismap_prob[r];
        total += log((1.0 - e_r) * mixture + e_r);
    }

    if (uses_depth_term()) {
        // Add ln of the tempered density f^beta / Z_beta at the observed count, with beta the
        // depth weight. Z_beta grows with lambda when beta < 1, so without it we would favour
        // the genotypes that expect more reads.
        double lambda = expected_reads(genotype);
        total += depth_weight * ln_poisson_pmf(observed_reads(), lambda)
                 - ln_tempered_poisson_normaliser(lambda, depth_weight);
    }

    return total;
}

double AlleleReadLikelihoods::achievable_gap(const vector<int>& called,
                                              const vector<int>& runner_up) const {
    if (called.empty() || runner_up.empty() || n_reads == 0) {
        return 0.0;
    }

    vector<double> w_called = mixture_weights(called);
    vector<double> w_runner = mixture_weights(runner_up);

    // An ideal read from allele a has rel 1 for a and 0 elsewhere, so the mixture of a genotype
    // for it is the total weight the genotype places on a. A homozygote sees 1, and a
    // heterozygote sees the weight of the read's own haplotype.
    //
    // Each haplotype of `called` gets a slot, with its share of the ideal reads and the mixture
    // each genotype shows for a read from it.
    struct Slot {
        double share;
        double mix_called;
        double mix_runner;
    };
    vector<Slot> slots;
    slots.reserve(called.size());
    for (size_t i = 0; i < called.size(); ++i) {
        int a = called[i];
        if (a < 0 || (size_t)a >= n_alleles) {
            // Star and missing alleles carry no sequence to tell apart.
            continue;
        }
        double mix_called = 0.0;
        double mix_runner = 0.0;
        for (size_t j = 0; j < called.size(); ++j) {
            if (called[j] == a) {
                mix_called += w_called[j];
            }
        }
        for (size_t j = 0; j < runner_up.size(); ++j) {
            if (runner_up[j] == a) {
                mix_runner += w_runner[j];
            }
        }
        slots.push_back({w_called[i], mix_called, mix_runner});
    }
    if (slots.empty()) {
        return 0.0;
    }

    // Give every ideal read e_r at the floor, whatever e_r the site's own reads have. With their
    // own e_r, a poorly mapped site would get a small denominator, and a weak call there would
    // look strong.
    double e = mismap_floor;
    double per_read = 0.0;
    for (const Slot& s : slots) {
        per_read += s.share * (log((1.0 - e) * s.mix_called + e)
                               - log((1.0 - e) * s.mix_runner + e));
    }
    double total = per_read * (double)n_reads;

    // Identical genotypes give 0, and rounding can put the sum a hair below it.
    return max(total, 0.0);
}

vector<vector<int>> AlleleReadLikelihoods::enumerate_genotypes(size_t num_alleles, int ploidy) {
    vector<vector<int>> genotypes;
    if (num_alleles == 0 || ploidy <= 0) {
        return genotypes;
    }

    // VCF puts a genotype with alleles a_1 <= ... <= a_P at an index given by a recursion over
    // its largest allele. Generating non-decreasing tuples in colexicographic order gives the
    // same order. For diploid this is (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), ...
    vector current(ploidy, 0);

    function<void(int, int)> recurse = [&](int depth, int max_allele) {
        if (depth < 0) {
            genotypes.push_back(current);
            return;
        }
        for (int allele = 0; allele <= max_allele; ++allele) {
            current[depth] = allele;
            recurse(depth - 1, allele);
        }
    };

    for (int top = 0; top < (int)num_alleles; ++top) {
        current[ploidy - 1] = top;
        recurse(ploidy - 2, top);
    }

    return genotypes;
}

vector<pair<vector<int>, double>> AlleleReadLikelihoods::score_genotypes(int ploidy) const {
    vector<pair<vector<int>, double>> scored;
    for (auto& genotype : enumerate_genotypes(n_alleles, ploidy)) {
        double ll = genotype_likelihood(genotype);
        scored.emplace_back(genotype, ll);
    }
    return scored;
}

void AlleleReadLikelihoods::dump(ostream& out, const string& site_name) const {
    out << "#site\tread\tmismap_prob\tbest_ln";
    for (size_t a = 0; a < n_alleles; ++a) {
        out << "\tallele_" << a;
    }
    out << endl;
    for (size_t r = 0; r < n_reads; ++r) {
        out << site_name << "\t" << (read_names.empty() ? string(".") : read_names[r]) << "\t"
            << read_mismap_prob[r] << "\t" << read_best_ln[r];
        for (size_t a = 0; a < n_alleles; ++a) {
            out << "\t" << rel(r, a);
        }
        out << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////
// AlleleReadLikelihoodsBuilder
////////////////////////////////////////////////////////////////////////////////

AlleleReadLikelihoodsBuilder::AlleleReadLikelihoodsBuilder(size_t num_alleles, double min_mismap,
                                                           double max_mismap)
    : n_alleles(num_alleles), min_mismap(min_mismap), max_mismap(max_mismap) {
}

bool AlleleReadLikelihoodsBuilder::add_read(const vector<double>& raw_ln_likelihood,
                                            double mismap_prob, const string& name,
                                            size_t read_length, const Position* start) {
    assert(raw_ln_likelihood.size() == n_alleles);

    // Divide the row by the best fit of the read over all alleles at the site. The divisor does
    // not depend on the genotype, so it cancels from every genotype comparison.
    double best = -numeric_limits<double>::infinity();
    for (double ll : raw_ln_likelihood) {
        best = max(best, ll);
    }

    if (!(best > -numeric_limits<double>::infinity())) {
        // The read placed on no allele, so there is no row maximum to divide by. Drop it and
        // count it.
        ++unplaceable;
        return false;
    }

    rows.emplace_back();
    rows.back().reserve(n_alleles);
    for (double ll : raw_ln_likelihood) {
        // Get a value in [0, 1], with exactly 0 for -inf.
        rows.back().push_back(exp(ll - best));
    }

    if (read_length > 0) {
        read_length_total += (double)read_length;
        ++read_length_count;
    }
    mismap_probs.push_back(min(max(mismap_prob, min_mismap), max_mismap));
    best_lns.push_back(best);
    names.push_back(name);
    if (start != nullptr) {
        starts.emplace_back(start->node_id(), start->offset(), start->is_reverse());
    } else {
        starts.emplace_back(0, 0, false);
    }
    return true;
}

AlleleReadLikelihoods AlleleReadLikelihoodsBuilder::build() {
    size_t n_reads = rows.size();

    // Put the rows in a canonical order, so that a floating-point sum over the reads gives the
    // same answer whatever order they arrived in. The name and the start identify an alignment.
    // The values of the row break any remaining tie, between rows that are then interchangeable.
    order.resize(n_reads);
    for (size_t i = 0; i < n_reads; ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (names[a] != names[b]) {
            return names[a] < names[b];
        }
        if (starts[a] != starts[b]) {
            return starts[a] < starts[b];
        }
        if (best_lns[a] != best_lns[b]) {
            return best_lns[a] < best_lns[b];
        }
        if (mismap_probs[a] != mismap_probs[b]) {
            return mismap_probs[a] < mismap_probs[b];
        }
        return rows[a] < rows[b];
    });

    vector<double> matrix;
    matrix.reserve(n_reads * n_alleles);
    vector<double> sorted_mismap;
    vector<double> sorted_best;
    vector<string> sorted_names;
    sorted_mismap.reserve(n_reads);
    sorted_best.reserve(n_reads);
    sorted_names.reserve(n_reads);
    for (size_t i : order) {
        matrix.insert(matrix.end(), rows[i].begin(), rows[i].end());
        sorted_mismap.push_back(mismap_probs[i]);
        sorted_best.push_back(best_lns[i]);
        sorted_names.push_back(std::move(names[i]));
    }

    AlleleReadLikelihoods result;
    result.set_contents(n_reads, n_alleles, std::move(matrix), std::move(sorted_mismap),
                        std::move(sorted_best), std::move(sorted_names), unplaceable);
    result.set_mismap_floor(min_mismap);
    if (read_length_count > 0) {
        // Set R whichever mixture is in use, since the depth term's lambda = rate * (L + R - 1)
        // needs it too.
        result.set_mean_read_length(read_length_total / (double)read_length_count);
    }
    if (!allele_lengths.empty() && read_length_count > 0) {
        result.set_length_weights(allele_lengths, read_length_total / (double)read_length_count);
        result.set_unique_lengths(std::move(unique_lengths));
    }
    return result;
}

////////////////////////////////////////////////////////////////////////////////
// GraphAlignedAlleleLikelihoodCalculator
////////////////////////////////////////////////////////////////////////////////

GraphAlignedAlleleLikelihoodCalculator::GraphAlignedAlleleLikelihoodCalculator(
    const PathHandleGraph& graph, const SiteReadSource& read_source,
    const EditAlignmentScorer& qual_scorer, const EditAlignmentScorer& plain_scorer, const Params& params)
    : graph(graph), read_source(read_source),
      qual_scorer(qual_scorer), plain_scorer(plain_scorer), params(params) {
}

vector<GraphAlignedAlleleLikelihoodCalculator::AlleleStep>
GraphAlignedAlleleLikelihoodCalculator::get_allele_steps(const Traversal& walk,
                                                         const StepSequences* own) const {
    vector<AlleleStep> steps;
    steps.reserve(walk.size());
    for (size_t i = 0; i < walk.size(); ++i) {
        AlleleStep step;
        step.node_id = graph.get_id(walk[i]);
        step.backward = graph.get_is_reverse(walk[i]);
        step.sequence = graph.get_sequence(walk[i]);
        if (own != nullptr && own->at(i) != step.sequence) {
            step.sequence = own->at(i);
            step.own_sequence = true;
        }
        steps.push_back(std::move(step));
    }
    return steps;
}

size_t GraphAlignedAlleleLikelihoodCalculator::spelled_length(const vector<AlleleStep>& steps,
                                                              size_t from, size_t to) {
    size_t length = 0;
    for (size_t j = from; j < to; ++j) {
        length += steps[j].sequence.size();
    }
    return length;
}

/// Return the stretch of an alignment that lies in a site.
///
/// The stretch runs from the first mapping on a node of `site_nodes` to the last, with the read
/// sequence and qualities cut to match. Flipping only this stretch of a reverse-strand read keeps
/// the cost of flipping proportional to the site rather than to the read, which matters for long
/// reads.
///
/// Flipping the stretch gives the same scores as flipping the whole read. A step at read offset
/// `off` with length `len` moves to `to - off - len` in the flipped stretch and to
/// `total - off - len` in the flipped read. The two differ by the constant `total - to`, and so do
/// the two flipped sequences, so each step still indexes the same bases.
static Alignment site_span_of(const SiteRead& read, const unordered_set<nid_t>& site_nodes) {
    const Alignment& aln = *read.aln;
    const Path& path = aln.path();
    int64_t first = -1;
    int64_t last = -1;
    size_t from = 0;
    size_t to = 0;
    if (read.indexed()) {
        // Take both ends from the mappings the source found, and the read offset before each,
        // without walking the rest of the read.
        for (size_t k = 0; k < read.mapping_count; ++k) {
            auto i = (int64_t)read.mappings[k];
            if (!site_nodes.contains(path.mapping(i).position().node_id())) {
                continue;
            }
            if (first < 0) {
                first = i;
                from = read.read_offsets[i];
            }
            last = i;
            to = (size_t)read.read_offsets[i] + (size_t)mapping_to_length(path.mapping(i));
        }
    } else {
        size_t offset = 0;
        for (int64_t i = 0; i < path.mapping_size(); ++i) {
            auto to_length = (size_t)mapping_to_length(path.mapping(i));
            if (site_nodes.contains(path.mapping(i).position().node_id())) {
                if (first < 0) {
                    first = i;
                    from = offset;
                }
                last = i;
                to = offset + to_length;
            }
            offset += to_length;
        }
    }
    if (first < 0) {
        // The caller has already found that the read visits the site, so this cannot happen.
        // Flipping the whole alignment would still give the same scores.
        return aln;
    }

    Alignment span;
    span.set_name(aln.name());
    span.set_mapping_quality(aln.mapping_quality());
    to = min(to, aln.sequence().size());
    if (from < to) {
        span.set_sequence(aln.sequence().substr(from, to - from));
        if (aln.quality().size() >= to) {
            span.set_quality(aln.quality().substr(from, to - from));
        }
    }
    Path* out = span.mutable_path();
    for (int64_t i = first; i <= last; ++i) {
        *out->add_mapping() = path.mapping(i);
    }
    return span;
}

bool GraphAlignedAlleleLikelihoodCalculator::get_read_steps(
    const SiteRead& read, const unordered_set<nid_t>& site_nodes,
    const unordered_set<nid_t>& boundary_nodes, vector<ReadStep>& steps_out) const {

    steps_out.clear();
    const Alignment& aln = *read.aln;
    const Path& path = aln.path();

    bool touches_interior = false;

    auto take = [&](int64_t i, size_t read_offset) {
        const Mapping& mapping = path.mapping(i);
        nid_t node_id = mapping.position().node_id();
        if (!site_nodes.contains(node_id)) {
            return;
        }
        ReadStep step;
        step.node_id = node_id;
        step.backward = mapping.position().is_reverse();
        step.read_offset = read_offset;
        step.read_length = (size_t)mapping_to_length(mapping);
        step.mapping = &mapping;
        steps_out.push_back(step);

        if (!boundary_nodes.contains(node_id)) {
            touches_interior = true;
        }
    };

    if (read.indexed()) {
        // Take the mappings the source found inside the queried ranges, which it hands over in
        // ascending order, with their read offsets.
        for (size_t k = 0; k < read.mapping_count; ++k) {
            take((int64_t)read.mappings[k], read.read_offsets[read.mappings[k]]);
        }
    } else {
        // Count the read offset across every mapping, including those outside the site, so that
        // the offsets of the ones inside are right.
        size_t read_offset = 0;
        for (int64_t i = 0; i < path.mapping_size(); ++i) {
            take(i, read_offset);
            read_offset += (size_t)mapping_to_length(path.mapping(i));
        }
    }

    if (steps_out.empty()) {
        return false;
    }

    // A read can tell the alleles apart if it touches an interior node, or if it moves between
    // two different nodes of the site. The second case includes a read that crosses straight from
    // one boundary node to the other, which supports an allele that deletes the site's interior.
    bool uses_internal_edge = false;
    for (size_t i = 1; i < steps_out.size(); ++i) {
        if (steps_out[i].node_id != steps_out[i - 1].node_id) {
            uses_internal_edge = true;
            break;
        }
    }

    if (!touches_interior && !uses_internal_edge) {
        // The read stays inside one boundary node, and fits every allele equally.
        steps_out.clear();
        return false;
    }

    return true;
}

bool GraphAlignedAlleleLikelihoodCalculator::read_is_reverse_of_alleles(
    const vector<ReadStep>& read_steps,
    const unordered_map<nid_t, bool>& allele_orientations) const {

    // Vote over the read's visits, so that one node the alleles visit in different orientations
    // cannot decide it.
    size_t agree = 0;
    size_t disagree = 0;
    for (const ReadStep& step : read_steps) {
        auto found = allele_orientations.find(step.node_id);
        if (found == allele_orientations.end()) {
            continue;
        }
        if (found->second == step.backward) {
            ++agree;
        } else {
            ++disagree;
        }
    }
    return disagree > agree;
}

int32_t GraphAlignedAlleleLikelihoodCalculator::score_shared_node(
    const Alignment& aln, const ReadStep& step, const EditAlignmentScorer& read_scorer,
    double& nat_adjust, const AlleleStep* allele_step) const {

    int32_t score = 0;
    const string& seq = aln.sequence();
    const string& qual = aln.quality();
    size_t read_pos = step.read_offset;
    const string* own = (allele_step != nullptr && allele_step->own_sequence)
                            ? &allele_step->sequence : nullptr;
    auto node_pos = (size_t)step.mapping->position().offset();

    for (int64_t i = 0; i < step.mapping->edit_size(); ++i) {
        const Edit& edit = step.mapping->edit(i);
        auto to_len = (size_t)edit.to_length();
        auto from_len = (size_t)edit.from_length();

        if (read_pos + to_len > seq.size()) {
            // The alignment is malformed. Do not read past the sequence.
            break;
        }

        if (from_len == to_len && own != nullptr) {
            score += score_substitution(aln, read_pos, to_len, *own, node_pos, read_scorer);
        } else if (from_len == to_len) {
            auto begin = seq.begin() + read_pos;
            auto end = begin + to_len;
            auto qual_begin = qual.empty() ? seq.begin() : qual.begin() + read_pos;
            if (edit.sequence().empty()) {
                score += read_scorer.score_exact_match(begin, end, qual_begin);
            } else {
                score += read_scorer.score_mismatch(begin, end, qual_begin);
            }
        } else if (from_len > to_len) {
            // The read lacks bases of the node.
            score += read_scorer.score_gap(from_len - to_len);
        } else {
            // The read carries bases the node lacks.
            score += read_scorer.score_gap(to_len - from_len);
            nat_adjust += params.insertion_gap_nats;
        }

        read_pos += to_len;
        node_pos += from_len;
    }

    return score;
}

int32_t GraphAlignedAlleleLikelihoodCalculator::score_substitution(
    const Alignment& aln, size_t read_offset, size_t length, const string& allele_bases,
    size_t allele_offset, const EditAlignmentScorer& read_scorer) const {

    int32_t score = 0;
    const string& seq = aln.sequence();
    const string& qual = aln.quality();

    // Walk the two sequences together, and score each run of equal bases, or of differing
    // bases, in one call.
    size_t i = 0;
    while (i < length) {
        if (read_offset + i >= seq.size() || allele_offset + i >= allele_bases.size()) {
            break;
        }
        bool matching = seq[read_offset + i] == allele_bases[allele_offset + i];
        size_t run = 1;
        while (i + run < length && read_offset + i + run < seq.size() &&
               allele_offset + i + run < allele_bases.size() &&
               (seq[read_offset + i + run] == allele_bases[allele_offset + i + run]) == matching) {
            ++run;
        }

        auto begin = seq.begin() + read_offset + i;
        auto end = begin + run;
        auto qual_begin = qual.empty() ? seq.begin() : qual.begin() + read_offset + i;
        score += matching ? read_scorer.score_exact_match(begin, end, qual_begin)
                          : read_scorer.score_mismatch(begin, end, qual_begin);
        i += run;
    }

    return score;
}

int32_t GraphAlignedAlleleLikelihoodCalculator::score_node_substitution(
    const Alignment& aln, const ReadStep& read_step, const AlleleStep& allele_step,
    const EditAlignmentScorer& read_scorer, double& nat_adjust) const {

    // Score the overlapping bases one by one, so that an equal-length pair of nodes, as at a SNP,
    // costs mismatches rather than two gaps.
    const size_t shared = min(read_step.read_length, allele_step.sequence.size());
    int32_t score = score_substitution(aln, read_step.read_offset, shared, allele_step.sequence, 0,
                                       read_scorer);

    // Charge the difference in length as a gap.
    const size_t difference = read_step.read_length > allele_step.sequence.size()
                                  ? read_step.read_length - allele_step.sequence.size()
                                  : allele_step.sequence.size() - read_step.read_length;
    if (difference > 0) {
        score += read_scorer.score_gap(difference);
        if (read_step.read_length > allele_step.sequence.size()) {
            nat_adjust += params.insertion_gap_nats;
        }
    }
    return score;
}

/// Return a visit, a node ID and an orientation, packed into one integer.
///
/// Sorted vectors of these keys answer membership by binary search, which is cheaper than a hash
/// lookup for the short traversals of most sites.
static inline int64_t step_key(nid_t node, bool backward) {
    return ((int64_t)node << 1) | (int64_t)backward;
}

void GraphAlignedAlleleLikelihoodCalculator::prepare_read_scratch(
    const Alignment& aln, const vector<ReadStep>& read_steps,
    const EditAlignmentScorer& read_scorer, ReadScratch& scratch) const {

    const size_t m = read_steps.size();
    scratch.own.assign(m, 0);
    scratch.own_nats.assign(m, 0.0);
    scratch.keys.clear();
    scratch.keys.reserve(m);
    for (size_t i = 0; i < m; ++i) {
        double nats = 0.0;
        scratch.own[i] = score_shared_node(aln, read_steps[i], read_scorer, nats);
        scratch.own_nats[i] = nats;
        scratch.keys.push_back(step_key(read_steps[i].node_id, read_steps[i].backward));
    }
    std::ranges::sort(scratch.keys);
}

vector<int64_t> GraphAlignedAlleleLikelihoodCalculator::sorted_allele_keys(
    const vector<AlleleStep>& allele_steps) {

    vector<int64_t> keys;
    keys.reserve(allele_steps.size());
    for (const AlleleStep& a : allele_steps) {
        keys.push_back(step_key(a.node_id, a.backward));
    }
    std::ranges::sort(keys);
    return keys;
}

GraphAlignedAlleleLikelihoodCalculator::AlleleStepPositions
GraphAlignedAlleleLikelihoodCalculator::index_allele_steps(
    const vector<AlleleStep>& allele_steps) {
    AlleleStepPositions positions;
    positions.reserve(allele_steps.size() * 2);
    for (size_t j = 0; j < allele_steps.size(); ++j) {
        positions[step_key(allele_steps[j].node_id, allele_steps[j].backward)]
            .push_back((uint32_t)j);
    }
    return positions;
}

/// Return the first index at or after `from` at which `positions` lists the visit `key`, or the
/// largest size_t if there is none.
static size_t next_occurrence(const unordered_map<int64_t, vector<uint32_t>>& positions,
                              int64_t key, size_t from) {
    auto found = positions.find(key);
    if (found == positions.end()) {
        return numeric_limits<size_t>::max();
    }
    auto at = std::ranges::lower_bound(found->second, (uint32_t)from);
    return at == found->second.end() ? numeric_limits<size_t>::max() : *at;
}

int32_t GraphAlignedAlleleLikelihoodCalculator::score_by_greedy_pairing(
    const Alignment& aln, const vector<ReadStep>& read_steps,
    const vector<AlleleStep>& allele_steps, const AlleleStepPositions& allele_positions,
    const EditAlignmentScorer& read_scorer, double& nat_adjust) const {

    int32_t score = 0;
    size_t allele_index = 0;
    bool have_anchor = false;
    size_t bases_accounted = 0;

    // Make a run of unpaired read visits after the first same-visit pair one gap, as optimal
    // pairing does, so that a pairing scores the same whichever search found it. The first visit
    // of the run opens the gap and the others extend it. Before that pair, each unpaired read
    // visit is a gap of its own. A visit with no read bases is no gap, and neither starts nor ends
    // a run.
    const int32_t extend_per_base = read_scorer.score_gap(2) - read_scorer.score_gap(1);
    bool in_insertion = false;
    auto leave_unpaired = [&](const ReadStep& step) {
        if (step.read_length == 0) {
            return;
        }
        if (in_insertion) {
            score += (int32_t)step.read_length * extend_per_base;
        } else {
            score += read_scorer.score_gap(step.read_length);
            in_insertion = have_anchor;
        }
        nat_adjust += params.insertion_gap_nats;
        bases_accounted += step.read_length;
    };

    // Find the last read visit of each kind, to tell whether the read makes an allele visit later.
    unordered_map<int64_t, size_t> read_last_visit;
    read_last_visit.reserve(read_steps.size() * 2);
    for (size_t i = 0; i < read_steps.size(); ++i) {
        read_last_visit[step_key(read_steps[i].node_id, read_steps[i].backward)] = i;
    }

    for (size_t read_index = 0; read_index < read_steps.size(); ++read_index) {
        const ReadStep& read_step = read_steps[read_index];
        const int64_t read_key = step_key(read_step.node_id, read_step.backward);
        const size_t found = next_occurrence(allele_positions, read_key, allele_index);

        if (found != numeric_limits<size_t>::max()) {
            // The allele makes this read visit. Charge the allele visits skipped since the
            // previous pair as a deletion. Allele visits before the first pair are outside the
            // read and cost nothing, so that a read is not penalised for being short.
            const size_t deleted = have_anchor ? spelled_length(allele_steps, allele_index, found)
                                               : 0;
            if (deleted > 0) {
                score += read_scorer.score_gap(deleted);
            }
            score += score_shared_node(aln, read_step, read_scorer, nat_adjust,
                                       &allele_steps[found]);
            bases_accounted += read_step.read_length;
            have_anchor = true;
            in_insertion = false;
            allele_index = found + 1;
            continue;
        }

        if (allele_index >= allele_steps.size()) {
            // The allele is used up but the read goes on. Charge the rest of the read as an
            // insertion, so that every allele is scored over the same read bases.
            leave_unpaired(read_step);
            continue;
        }

        // A visit that the read and the allele both make may not be substituted for another, as
        // in optimal pairing. If the read makes the current allele visit anywhere, or the allele
        // made this read visit before the current allele visit, leave the read visit unpaired.
        // The allele visit stays to be paired later, or to be charged as a deletion at the next
        // pair.
        const AlleleStep& allele_step = allele_steps[allele_index];
        if (read_last_visit.contains(step_key(allele_step.node_id, allele_step.backward))
            || allele_positions.contains(read_key)) {
            leave_unpaired(read_step);
            continue;
        }

        // The read and the allele put different nodes here.
        score += score_node_substitution(aln, read_step, allele_step, read_scorer, nat_adjust);
        bases_accounted += read_step.read_length;
        in_insertion = false;
        ++allele_index;
    }

    // Check that every read base inside the site was scored, whatever the allele.
    size_t window_bases = 0;
    for (const ReadStep& read_step : read_steps) {
        window_bases += read_step.read_length;
    }
    assert(bases_accounted == window_bases);

    return score;
}

/**
 * A dynamic-programming grid that finds the optimal pairing of the visits of one read with those
 * of one allele.
 *
 * Base-level edits are still read off the mapper's alignment. Only the pairing of visits is
 * searched. Row i of the grid holds the best scores after the first i read visits, and column j
 * after the first j allele visits. Each cell has four states:
 *
 *   P  No shared visit is paired yet, in the leading flank. Allele visits cost nothing here, since
 *      allele sequence before the read is not the read's to explain. Read visits are still
 *      charged. Pairing a shared visit leaves P for M, and a substitution does not.
 *   M  This read visit is paired with this allele visit.
 *   I  We are inside a run of read visits the allele lacks.
 *   D  We are inside a run of allele visits the read lacks.
 *
 * Pairing the same visit scores the read's own edits in it. Pairing two visits that each appear
 * nowhere in the other sequence is a substitution. A visit the read and the allele share may not
 * be paired with a different visit, though it may be left unpaired.
 *
 * I and D may follow each other directly. Base-level alignment usually forbids that, as a
 * duplicate of a substitution. Here a substitution scores two different nodes base by base, which
 * gives a different score from a deletion next to an insertion.
 *
 * A run of inserted read visits after the first match is one gap, as in greedy pairing. A read
 * visit with no bases is no gap when left unpaired, and its row passes every state through.
 *
 * Each cell also carries the `insertion_gap_nats` total of its own best path.
 */
class GraphAlignedAlleleLikelihoodCalculator::PairingGrid {
public:
    PairingGrid(const GraphAlignedAlleleLikelihoodCalculator& calculator, const Alignment& aln,
                const vector<ReadStep>& read_steps, const vector<AlleleStep>& allele_steps,
                const ReadScratch& scratch, const vector<int64_t>& allele_keys,
                const EditAlignmentScorer& read_scorer)
        : calculator(calculator), aln(aln), read_steps(read_steps), allele_steps(allele_steps),
          scratch(scratch), allele_keys(allele_keys), read_scorer(read_scorer),
          extend_per_base(read_scorer.score_gap(2) - read_scorer.score_gap(1)) {
    }

    /**
     * Fill the grid and return the score of the best pairing.
     *
     * The best pairing ends in state P, M or I in the last row. It may not end in D, which has
     * charged allele visits after the read. It may end in P, so that a read that pairs nothing
     * still gets a finite score.
     *
     * @param nat_adjust Gets the `insertion_gap_nats` of the best pairing added to it.
     */
    int32_t best_score(double& nat_adjust);

private:
    /// A score in one state of one cell, with the `insertion_gap_nats` total of its path.
    struct Cell {
        int32_t score;
        double nats;
    };

    /// One row of the grid, with a cell for each column in each state.
    struct Row {
        vector<Cell> P;
        vector<Cell> M;
        vector<Cell> I;
        vector<Cell> D;

        /// Fill the row with `width` unreachable cells.
        void reset(size_t width);
        /// Make the cells in columns [from, to) unreachable.
        void clear(size_t from, size_t to);
    };

    /// The score of an unreachable state, low enough that adding to it cannot overflow.
    static constexpr int32_t NEG = numeric_limits<int32_t>::min() / 4;
    static constexpr Cell NONE{NEG, 0.0};

    /// The number of cells above which only a band of the grid is filled.
    static constexpr size_t BAND_THRESHOLD = 20000;
    /// How far the band extends on each side of the shared visits.
    static constexpr size_t BAND_WIDTH = 64;

    /// Return the higher-scoring cell, or `a` on a tie.
    static Cell better(const Cell& a, const Cell& b) { return a.score >= b.score ? a : b; }

    /// Return `c` with `d` added to its score and `dn` to its nats, or an unreachable cell.
    static Cell plus(const Cell& c, int32_t d, double dn) {
        return c.score == NEG ? Cell{NEG, 0.0} : Cell{c.score + d, c.nats + dn};
    }

    /// Fill `band_from` and `band_to` with the first and last column to fill in each row.
    void find_band();

    /// Return the score of pairing read visit `i` (from 1) with allele visit `as`, or NEG where
    /// the pairing is forbidden.
    int32_t pair_score(size_t i, const AlleleStep& as, bool is_match, double& pair_nats) const;

    /// Fill row `i` (from 1) of the grid from the row before it.
    void fill_row(size_t i, const Row& previous, Row& current) const;

    const GraphAlignedAlleleLikelihoodCalculator& calculator;
    const Alignment& aln;
    const vector<ReadStep>& read_steps;
    const vector<AlleleStep>& allele_steps;
    const ReadScratch& scratch;
    const vector<int64_t>& allele_keys;
    const EditAlignmentScorer& read_scorer;
    /// The cost of extending an open gap by one more base.
    ///
    /// We take it as the difference between a two-base and a one-base gap, because the scorer
    /// does not expose its gap extension score.
    const int32_t extend_per_base;
    /// Whether only a band of each row is filled.
    bool banded = false;
    /// The first column filled in each row, when banded.
    vector<size_t> band_from;
    /// The last column filled in each row, when banded.
    vector<size_t> band_to;
};

void GraphAlignedAlleleLikelihoodCalculator::PairingGrid::Row::reset(size_t width) {
    P.assign(width, NONE);
    M.assign(width, NONE);
    I.assign(width, NONE);
    D.assign(width, NONE);
}

void GraphAlignedAlleleLikelihoodCalculator::PairingGrid::Row::clear(size_t from, size_t to) {
    for (size_t j = from; j < to; ++j) {
        P[j] = NONE;
        M[j] = NONE;
        I[j] = NONE;
        D[j] = NONE;
    }
}

void GraphAlignedAlleleLikelihoodCalculator::PairingGrid::find_band() {
    // A read and an allele usually differ at a few nodes, so the best pairing stays near the path
    // through their shared visits. The match cell of a shared visit at read step i_s and allele
    // step j_s is (i_s + 1, j_s + 1). The band of row i spans the projections, along the
    // diagonal, of the last shared match cell in row i or earlier and of the first one after row
    // i, widened by BAND_WIDTH on each side. Between two shared visits the band therefore holds
    // every cell of a deletion or an insertion that joins them, however long. Where they share a
    // diagonal, the band is centred on it.
    const size_t m = read_steps.size();
    const size_t n = allele_steps.size();
    band_from.assign(m + 1, 0);
    band_to.assign(m + 1, 0);
    size_t a_index = 0;
    vector<pair<size_t, size_t>> shared;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = a_index; j < n; ++j) {
            if (allele_steps[j].node_id == read_steps[i].node_id &&
                allele_steps[j].backward == read_steps[i].backward) {
                shared.emplace_back(i + 1, j + 1);
                a_index = j + 1;
                break;
            }
        }
    }
    // shared[s - 1] is the last match cell in row i or earlier, and shared[s] the first after it.
    size_t s = 0;
    for (size_t i = 0; i <= m; ++i) {
        while (s < shared.size() && shared[s].first <= i) {
            ++s;
        }
        // With no shared visit at all, use the diagonal i == j.
        size_t from = i;
        size_t to = i;
        if (s > 0) {
            from = shared[s - 1].second + (i - shared[s - 1].first);
            to = from;
        }
        if (s < shared.size()) {
            const size_t ahead = shared[s].second - min(shared[s].second, shared[s].first - i);
            from = s > 0 ? min(from, ahead) : ahead;
            to = s > 0 ? max(to, ahead) : ahead;
        }
        from = min(from, n);
        to = min(to, n);
        band_from[i] = from > BAND_WIDTH ? from - BAND_WIDTH : 1;
        band_to[i] = min(n, to + BAND_WIDTH);
    }
}

int32_t GraphAlignedAlleleLikelihoodCalculator::PairingGrid::pair_score(
    size_t i, const AlleleStep& as, bool is_match, double& pair_nats) const {
    const ReadStep& rs = read_steps[i - 1];
    if (is_match && as.own_sequence) {
        // Score against the allele's own sequence, which the read's stored score of the node
        // does not know.
        return calculator.score_shared_node(aln, rs, read_scorer, pair_nats, &as);
    }
    if (is_match) {
        pair_nats = scratch.own_nats[i - 1];
        return scratch.own[i - 1];
    }
    // Membership is an exact test of sharing while neither sequence repeats a node. With a repeat
    // it forbids some substitutions it need not, so the result can be less than optimal, but is
    // still a valid pairing.
    if (!std::ranges::binary_search(allele_keys, step_key(rs.node_id, rs.backward)) &&
        !std::ranges::binary_search(scratch.keys, step_key(as.node_id, as.backward))) {
        return calculator.score_node_substitution(aln, rs, as, read_scorer, pair_nats);
    }
    return NEG;
}

void GraphAlignedAlleleLikelihoodCalculator::PairingGrid::fill_row(size_t i, const Row& previous,
                                                                   Row& current) const {
    const size_t n = allele_steps.size();
    const ReadStep& rs = read_steps[i - 1];
    const auto rlen = (int32_t)rs.read_length;
    const int32_t rgap = read_scorer.score_gap(rs.read_length);
    const double rnat = calculator.params.insertion_gap_nats;
    const bool no_bases = rs.read_length == 0;
    const Row& p = previous;
    Row& c = current;

    c.P[0] = no_bases ? p.P[0] : plus(p.P[0], rgap, rnat);
    c.M[0] = NONE;
    c.I[0] = no_bases ? p.I[0]
                      : better(plus(better(p.M[0], p.D[0]), rgap, rnat),
                               plus(p.I[0], rlen * extend_per_base, rnat));
    c.D[0] = NONE;

    size_t j_lo = 1;
    size_t j_hi = n;
    if (banded) {
        j_lo = band_from[i];
        j_hi = band_to[i];
        // Clear the cells outside the band, which hold values from two rows ago.
        c.clear(1, j_lo);
        c.clear(j_hi + 1, n + 1);
    }
    for (size_t j = j_lo; j <= j_hi; ++j) {
        const AlleleStep& as = allele_steps[j - 1];
        const auto alen = (int32_t)as.sequence.size();
        const int32_t agap = read_scorer.score_gap(as.sequence.size());
        const bool is_match = (as.node_id == rs.node_id && as.backward == rs.backward);

        double pair_nats = 0.0;
        const int32_t pair = pair_score(i, as, is_match, pair_nats);
        if (pair == NEG) {
            c.M[j] = NONE;
            c.P[j] = NONE;
        } else if (is_match) {
            // A match closes the flank, so it may be entered from P as well.
            c.M[j] = plus(better(better(p.M[j - 1], p.P[j - 1]), better(p.I[j - 1], p.D[j - 1])),
                          pair, pair_nats);
            c.P[j] = NONE;
        } else {
            c.M[j] = plus(better(p.M[j - 1], better(p.I[j - 1], p.D[j - 1])), pair, pair_nats);
            c.P[j] = plus(p.P[j - 1], pair, pair_nats);
        }
        // Delete leading allele visits for free.
        c.P[j] = better(c.P[j], c.P[j - 1]);

        if (no_bases) {
            // Leaving a visit with no bases unpaired keeps the state it follows.
            c.M[j] = better(c.M[j], p.M[j]);
            c.P[j] = better(c.P[j], p.P[j]);
            c.I[j] = p.I[j];
        } else {
            // Insert the read visit before any anchor.
            c.P[j] = better(c.P[j], plus(p.P[j], rgap, rnat));
            // A shared visit may be left unpaired, though it may not be substituted. A gap can
            // explain a read with poor edits inside a shared node better.
            c.I[j] = better(plus(better(p.M[j], p.D[j]), rgap, rnat),
                            plus(p.I[j], rlen * extend_per_base, rnat));
        }
        c.D[j] = better(plus(better(c.M[j - 1], c.I[j - 1]), agap, 0.0),
                        plus(c.D[j - 1], alen * extend_per_base, 0.0));
        if (no_bases) {
            c.D[j] = better(c.D[j], p.D[j]);
        }
    }
}

int32_t GraphAlignedAlleleLikelihoodCalculator::PairingGrid::best_score(double& nat_adjust) {
    const size_t m = read_steps.size();
    const size_t n = allele_steps.size();

    // Fill only a band of the grid at large sites. The band can miss the best pairing, so there
    // the search is an approximation.
    banded = m * n > BAND_THRESHOLD;
    if (banded) {
        find_band();
    }

    thread_local Row previous;
    thread_local Row current;
    previous.reset(n + 1);
    current.reset(n + 1);
    for (size_t j = 0; j <= n; ++j) {
        // Any allele prefix may be skipped for free.
        previous.P[j] = Cell{0, 0.0};
    }

    for (size_t i = 1; i <= m; ++i) {
        fill_row(i, previous, current);
        std::swap(previous, current);
    }

    // Each read visit is taken by exactly one transition, so every read base inside the site is
    // scored, whatever the allele.
    Cell best = NONE;
    for (size_t j = 0; j <= n; ++j) {
        best = better(best, better(previous.P[j], better(previous.M[j], previous.I[j])));
    }
    assert(best.score != NEG);
    nat_adjust += best.nats;
    return best.score;
}

int32_t GraphAlignedAlleleLikelihoodCalculator::score_by_optimal_pairing(
    const Alignment& aln, const vector<ReadStep>& read_steps,
    const vector<AlleleStep>& allele_steps, const ReadScratch& scratch,
    const vector<int64_t>& allele_keys, const EditAlignmentScorer& read_scorer,
    double& nat_adjust) const {
    PairingGrid grid(*this, aln, read_steps, allele_steps, scratch, allele_keys, read_scorer);
    return grid.best_score(nat_adjust);
}

void GraphAlignedAlleleLikelihoodCalculator::set_rate_reference(
    const PathPositionHandleGraph* position_graph, const vector<path_handle_t>& reference_paths) {
    rate_graph = position_graph;
    rate_paths = reference_paths;
    std::scoped_lock guard(window_bp_mutex);
    ref_bucket_counts.clear();
    ref_window_rate.clear();
}

bool GraphAlignedAlleleLikelihoodCalculator::rate_position(const MaterializedSite& site,
                                                          const vector<nid_t>& fallback_nodes,
                                                          size_t& path_index,
                                                          int64_t& position) const {
    if (rate_graph == nullptr || rate_paths.empty()) {
        return false;
    }
    // Place a boundary node at the earliest step on the first listed reference path that visits
    // it. Path handles, not node IDs, choose among the steps, so the choice survives renumbering.
    auto place = [&](nid_t node_id) {
        if (!rate_graph->has_node(node_id)) {
            return false;
        }
        bool found = false;
        rate_graph->for_each_step_on_handle(rate_graph->get_handle(node_id), [&](const step_handle_t& step) {
            path_handle_t path = rate_graph->get_path_handle_of_step(step);
            for (size_t i = 0; i < rate_paths.size(); ++i) {
                if (rate_paths[i] != path) {
                    continue;
                }
                auto pos = (int64_t)rate_graph->get_position_of_step(step);
                if (!found || i < path_index || (i == path_index && pos < position)) {
                    path_index = i;
                    position = pos;
                    found = true;
                }
                break;
            }
        });
        return found;
    };
    if (place(graph.get_id(site.start)) || place(graph.get_id(site.end))) {
        return true;
    }
    return std::ranges::any_of(fallback_nodes, place);
}

GraphAlignedAlleleLikelihoodCalculator::WindowReadStats
GraphAlignedAlleleLikelihoodCalculator::StartCounts::stats() const {
    WindowReadStats result;
    result.start_rate = (reads <= 0.0 || bp == 0) ? 0.0 : reads / (double)bp;
    result.mean_read_length = length_count > 0 ? length_total / (double)length_count : 0.0;
    return result;
}

GraphAlignedAlleleLikelihoodCalculator::StartCounts
GraphAlignedAlleleLikelihoodCalculator::count_starts(const vector<nid_t>& nodes, size_t bp) const {
    StartCounts counts;
    counts.bp = bp;
    vector<nid_t> sorted(nodes);
    std::ranges::sort(sorted);
    const auto duplicates = std::ranges::unique(sorted);
    sorted.erase(duplicates.begin(), duplicates.end());
    if (sorted.empty()) {
        return counts;
    }

    // Count reads as they are counted at a site. Under `depth_effective_reads` each counts as
    // 1 - e_r, with the same clamps and the same `use_mismap_term` switch. Counting them
    // differently would put a constant factor between N and lambda.
    //
    // 1 - e_r depends on the MAPQ of the read alone, so we tally reads by MAPQ and sum the tallies
    // in MAPQ order. The sum then does not depend on the order the read source returns reads in,
    // which can change with its fetch window. The length total is a sum of whole numbers, exact
    // in a double, so its order does not change it.
    //
    // The read source counts the starts itself, so a source that fetches reads by node-ID window
    // can answer from tallies it made while fetching windows for the sites.
    double& length_total = counts.length_total;
    size_t& length_count = counts.length_count;
    map<int32_t, size_t> by_mapq;
    read_source.for_each_read_start(sorted, [&](int32_t mapq, size_t reads, size_t bases) {
        length_total += (double)bases;
        length_count += reads;
        by_mapq[mapq] += reads;
    });
    for (const auto& [mapq, reads] : by_mapq) {
        if (!params.depth_effective_reads) {
            counts.reads += (double)reads;
            continue;
        }
        double mismap = params.use_mismap_term
                            ? phred_to_prob((double)mapq)
                            : params.min_mismap_prob;
        counts.reads += (double)reads
                        * (1.0 - min(max(mismap, params.min_mismap_prob), params.max_mismap_prob));
    }
    return counts;
}

GraphAlignedAlleleLikelihoodCalculator::StartCounts
GraphAlignedAlleleLikelihoodCalculator::bucket_counts(size_t path_index, int64_t bucket) const {
    path_handle_t path = rate_paths[path_index];
    int64_t lo = bucket * RATE_BUCKET;
    int64_t hi = min((int64_t)rate_graph->get_path_length(path), lo + RATE_BUCKET);
    if (bucket < 0 || lo >= hi) {
        return StartCounts();
    }
    pair<size_t, int64_t> key(path_index, bucket);
    {
        std::scoped_lock guard(window_bp_mutex);
        auto found = ref_bucket_counts.find(key);
        if (found != ref_bucket_counts.end()) {
            return found->second;
        }
    }

    // Take the reference nodes whose step begins in [lo, hi), each once however often the path
    // visits it.
    vector<nid_t> nodes;
    size_t bp = 0;
    unordered_set<nid_t> seen;
    step_handle_t step = rate_graph->get_step_at_position(path, (size_t)lo);
    while (step != rate_graph->path_end(path)) {
        auto start = (int64_t)rate_graph->get_position_of_step(step);
        if (start >= hi) {
            break;
        }
        if (start >= lo) {
            handle_t handle = rate_graph->get_handle_of_step(step);
            nid_t id = rate_graph->get_id(handle);
            if (seen.insert(id).second) {
                nodes.push_back(id);
                bp += rate_graph->get_length(handle);
            }
        }
        step = rate_graph->get_next_step(step);
    }

    StartCounts counts = count_starts(nodes, bp);
    std::scoped_lock guard(window_bp_mutex);
    ref_bucket_counts[key] = counts;
    return counts;
}

GraphAlignedAlleleLikelihoodCalculator::WindowReadStats
GraphAlignedAlleleLikelihoodCalculator::local_read_stats(
    const MaterializedSite& site, const vector<nid_t>& fallback_nodes,
    const vector<pair<nid_t, nid_t>>& site_ranges) const {

    if (site_ranges.empty() || params.depth_ploidy <= 0) {
        return WindowReadStats();
    }
    size_t path_index = 0;
    int64_t position = 0;
    if (!rate_position(site, fallback_nodes, path_index, position)) {
        ++id_fallbacks;
        return id_window_read_stats(site_ranges);
    }
    int64_t bucket = position / RATE_BUCKET;
    pair<size_t, int64_t> key(path_index, bucket);
    {
        std::scoped_lock guard(window_bp_mutex);
        auto found = ref_window_rate.find(key);
        if (found != ref_window_rate.end()) {
            return found->second;
        }
    }

    // Take the site's bucket and one on each side, so that the window reaches at least a bucket
    // past the site in both directions.
    StartCounts counts;
    for (int64_t b = bucket - 1; b <= bucket + 1; ++b) {
        counts.add(bucket_counts(path_index, b));
    }
    WindowReadStats stats = counts.stats();
    std::scoped_lock guard(window_bp_mutex);
    ref_window_rate[key] = stats;
    return stats;
}

GraphAlignedAlleleLikelihoodCalculator::WindowReadStats
GraphAlignedAlleleLikelihoodCalculator::id_window_read_stats(
    const vector<pair<nid_t, nid_t>>& site_ranges) const {

    nid_t lo = site_ranges.front().first;
    for (const auto& [range_start, range_end] : site_ranges) {
        lo = min(lo, range_start);
    }
    auto window_index = (size_t)(lo / RATE_ID_WINDOW);
    {
        std::scoped_lock guard(window_bp_mutex);
        auto found = window_rate.find(window_index);
        if (found != window_rate.end()) {
            return found->second;
        }
    }

    nid_t first = (nid_t)window_index * RATE_ID_WINDOW;
    nid_t last = first + RATE_ID_WINDOW - 1;

    // Ask the graph which IDs are nodes, since node IDs need not be dense.
    vector<nid_t> nodes;
    size_t bp = 0;
    for (nid_t id = first; id <= last; ++id) {
        if (graph.has_node(id)) {
            nodes.push_back(id);
            bp += graph.get_length(graph.get_handle(id));
        }
    }

    WindowReadStats stats = count_starts(nodes, bp).stats();
    std::scoped_lock guard(window_bp_mutex);
    window_rate[window_index] = stats;
    return stats;
}

/// Return the lengths of the longest prefix and the longest suffix that `a` and `b` share, which
/// do not overlap in either string.
static pair<size_t, size_t> shared_ends(const string& a, const string& b) {
    const size_t limit = min(a.size(), b.size());
    size_t prefix = 0;
    while (prefix < limit && a[prefix] == b[prefix]) {
        ++prefix;
    }
    size_t suffix = 0;
    while (suffix < limit - prefix && a[a.size() - 1 - suffix] == b[b.size() - 1 - suffix]) {
        ++suffix;
    }
    return {prefix, suffix};
}

/// Return the number of bases of `a` outside the longest prefix and suffix it shares with `b`.
static size_t bases_not_shared(const string& a, const string& b) {
    const auto [prefix, suffix] = shared_ends(a, b);
    return a.size() - prefix - suffix;
}

/// Return what it costs to align all of `read` to all of `allele`.
///
/// A mismatch costs the match score less the mismatch score of `scorer`, and each base of a gap
/// costs its one-base gap penalty.
static int32_t edit_cost(const string& read, const string& allele,
                         const EditAlignmentScorer& scorer) {
    const int32_t substitution = (int32_t)scorer.match + (int32_t)scorer.mismatch;
    const int32_t gap = -scorer.score_gap(1);
    vector<int32_t> cost(allele.size() + 1);
    for (size_t j = 0; j <= allele.size(); ++j) {
        cost[j] = gap * (int32_t)j;
    }
    for (size_t i = 1; i <= read.size(); ++i) {
        int32_t diagonal = cost[0];
        cost[0] = gap * (int32_t)i;
        for (size_t j = 1; j <= allele.size(); ++j) {
            const int32_t up = cost[j];
            cost[j] = std::min({diagonal + (read[i - 1] == allele[j - 1] ? 0 : substitution),
                                up + gap, cost[j - 1] + gap});
            diagonal = up;
        }
    }
    return cost[allele.size()];
}

/**
 * A record of how the alleles of a site inside nodes change the one walk they all take.
 *
 * The record compares the step sequences of each allele with the spelling of the walk. An allele
 * whose step sequences differ from the walk by one indel is an *indel allele*. It is scored as
 * its walk, and then each read's row is set by what the read says about the indel (see
 * `observe_indel`). We cannot compare its spelling base by base with a read aligned to the walk,
 * once the two differ in length. And the aligner may have put the read's own copy of the indel
 * anywhere in its repeat tract, on another node too.
 *
 * A site with an allele that both deletes and inserts, or with several alleles with step
 * sequences of which one changes the length, is *window scored* instead. Every allele is scored
 * as the walk. Then each allele is judged by what it costs to align the read's bases in the
 * window, between the unchanged bases around the change, to the allele's own bases there.
 *
 * For any other site the record is empty, and changes nothing.
 */
struct GraphAlignedAlleleLikelihoodCalculator::WalkChanges {
    /// An indel allele, with the rates of the read errors that would show or hide its indel.
    struct IndelAllele {
        size_t allele;
        WalkIndel indel;
        /// The end of the repeat tract of the indel, from `tract_end`.
        int64_t after;
        /// The ln of the rate of the error that turns the walk's sequence into the read's, for a
        /// read with the indel.
        double ln_shown;
        /// The ln of the rate of the error that turns the allele's sequence into the read's, for a
        /// read without the indel.
        double ln_hidden;
        /// Whether the indel has `LARGE_INDEL_MIN_LENGTH` bases or more, and is matched loosely.
        bool large;
    };

    /// The spelling of the walk.
    string walk_spelling;
    /// Where each step of the walk starts in `walk_spelling`, followed by the spelling's length.
    vector<int64_t> walk_starts;
    /// The index of each step of the walk, by `step_key`.
    unordered_map<int64_t, size_t> walk_index;
    /// The first walk base that the alleles change, counting the whole repeat tract of an indel.
    int64_t changed_lo = numeric_limits<int64_t>::max();
    /// The end of the walk bases that the alleles change, exclusive.
    ///
    /// In a window-scored site it is at most the last walk base, so that a read always has an
    /// aligned base at it to end the window.
    int64_t changed_hi = -1;
    /// The indel alleles, which are empty in a window-scored site.
    vector<IndelAllele> indel_alleles;
    /// Whether each allele is an indel allele.
    vector<bool> is_indel_allele;
    /// Whether the site is window scored.
    bool window_scored = false;
    /// The bases of each allele in the window, in a window-scored site.
    vector<string> window;

    /// Compare the alleles of `site` to their walk.
    WalkChanges(const HandleGraph& graph, const MaterializedSite& site, const Params& params);

    /// Return whether reads at the site need their rows changed by `adjust_row`.
    bool changes_reads() const { return !indel_alleles.empty() || window_scored; }

    /// Return whether allele `allele` is scored as the walk, rather than as its step sequences.
    bool scored_as_walk(size_t allele) const {
        return window_scored || is_indel_allele[allele];
    }

    /// Return a function that gives each base of the walk spelling, by offset.
    ///
    /// The function refers to this record, and must not outlive it.
    function<char(int64_t)> base_lookup() const {
        return [this](int64_t x) { return walk_spelling[x]; };
    }

    /**
     * Return the pin at walk point `point`, on the step holding walk base `base`.
     *
     * The pin is the offset of the point along the forward strand of the step's node.
     *
     * @param handle Gets the step's handle.
     */
    int64_t pin_at(const HandleGraph& graph, const Traversal& walk, int64_t point, int64_t base,
                   handle_t& handle) const;

    /**
     * Change the row of a read that follows the walk, by what it says about the changes.
     *
     * @param read_steps The read's visits inside the site.
     * @param sequence The read's bases, in the direction of the walk.
     * @param read_scorer The scorer that scored the row.
     * @param plain_scorer The scorer that charges alignments to window bases.
     * @param insertion_gap_nats The nats added for each gap of read bases an allele lacks.
     */
    void adjust_row(const vector<ReadStep>& read_steps, const string& sequence,
                    const EditAlignmentScorer& read_scorer,
                    const EditAlignmentScorer& plain_scorer, double insertion_gap_nats,
                    vector<double>& row) const;

private:
    /// Widen the changed bases to every place an aligner could put the length change of an
    /// allele that replaces the `deleted` walk bases from `prefix` with `inserted` of its own.
    void widen_for_replacement(const string& spelled, size_t prefix, size_t deleted,
                               size_t inserted);

    /// Record allele `allele`, which deletes `deleted` walk bases from `prefix` or inserts
    /// `inserted` of its own there, as an indel allele.
    void add_indel_allele(size_t allele, const string& spelled, size_t prefix, size_t deleted,
                          size_t inserted, const Params& params);

    /// Fill `window` with the bases of each allele in the window.
    void cut_windows(const HandleGraph& graph, const vector<Traversal>& traversals,
                     const AlleleSequences* sequences);

    /// Find the mappings of a read whose visits follow consecutive steps of the walk, and where
    /// each of their steps starts in the walk spelling.
    ///
    /// @return False if the read leaves the walk.
    bool follow_walk(const vector<ReadStep>& read_steps, vector<const Mapping*>& mappings,
                     vector<int64_t>& starts, size_t& first, size_t& last) const;

    /// Charge each allele in `row` with what it costs to align the read's window bases to the
    /// allele's own.
    void charge_window(const WalkAlignment& walked, const string& sequence,
                       const EditAlignmentScorer& plain_scorer, vector<double>& row) const;

    /// Set the row entries of each indel allele, and of the alleles scored as the walk, by
    /// whether the read shows the indel.
    ///
    /// The read covers the walk bases [reach_begin, reach_end), and runs on past either end that
    /// is a large sentinel.
    void observe_indels(const WalkAlignment& walked, int64_t reach_begin, int64_t reach_end,
                        const EditAlignmentScorer& read_scorer, double insertion_gap_nats,
                        vector<double>& row) const;
};

GraphAlignedAlleleLikelihoodCalculator::WalkChanges::WalkChanges(
    const HandleGraph& graph, const MaterializedSite& site, const Params& params)
    : is_indel_allele(site.alleles.walks.size(), false) {
    const vector<Traversal>& traversals = site.alleles.walks;
    if (!site.bounds_inside_nodes(graph) || site.alleles.edits.empty()
        || !std::ranges::all_of(traversals,
                                [&](const Traversal& t) { return t == traversals.front(); })) {
        return;
    }
    for (size_t i = 0; i < traversals.front().size(); ++i) {
        const handle_t& h = traversals.front()[i];
        walk_starts.push_back((int64_t)walk_spelling.size());
        walk_index[step_key(graph.get_id(h), graph.get_is_reverse(h))] = i;
        walk_spelling += graph.get_sequence(h);
    }
    walk_starts.push_back((int64_t)walk_spelling.size());

    bool resized = false;
    size_t own_alleles = 0;
    for (size_t a = 0; a < traversals.size(); ++a) {
        if (site.alleles.spelled_by_walk(a)) {
            continue;
        }
        const string spelled = site.alleles.sequence(graph, a);
        if (spelled == walk_spelling) {
            continue;
        }
        const auto [prefix, suffix] = shared_ends(spelled, walk_spelling);
        changed_lo = std::min(changed_lo, (int64_t)prefix);
        changed_hi = std::max(changed_hi, (int64_t)(walk_spelling.size() - suffix));
        ++own_alleles;
        if (spelled.size() == walk_spelling.size()) {
            continue;
        }
        resized = true;
        const size_t deleted = walk_spelling.size() - prefix - suffix;
        const size_t inserted = spelled.size() - prefix - suffix;
        if (deleted > 0 && inserted > 0) {
            window_scored = true;
            widen_for_replacement(spelled, prefix, deleted, inserted);
        } else {
            add_indel_allele(a, spelled, prefix, deleted, inserted, params);
        }
    }
    window_scored = window_scored || (resized && own_alleles > 1);
    if (window_scored) {
        const AlleleSequences sequences = site.alleles.all_step_sequences(graph);
        cut_windows(graph, traversals, &sequences);
    }
}

void GraphAlignedAlleleLikelihoodCalculator::WalkChanges::widen_for_replacement(
    const string& spelled, size_t prefix, size_t deleted, size_t inserted) {
    // The length change lies in the repeat at one end of the replaced bases or the other,
    // wherever a read's aligner put it, so cover the repeat at both ends.
    const function<char(int64_t)> base = base_lookup();
    const size_t shift = std::min(deleted, inserted);
    const size_t change = std::max(deleted, inserted) - shift;
    for (const size_t at : {prefix, prefix + shift}) {
        WalkIndel moved;
        moved.at = (int64_t)at;
        moved.deleted = deleted > inserted ? (uint32_t)change : 0;
        moved.inserted = deleted > inserted ? "" : spelled.substr(at, change);
        left_align(moved, 1, base);
        changed_lo = std::min(changed_lo, moved.at);
        changed_hi = std::max(changed_hi, tract_end(moved, (int64_t)walk_spelling.size(), base));
    }
}

void GraphAlignedAlleleLikelihoodCalculator::WalkChanges::add_indel_allele(
    size_t allele, const string& spelled, size_t prefix, size_t deleted, size_t inserted,
    const Params& params) {
    const function<char(int64_t)> base = base_lookup();
    IndelAllele entry;
    entry.allele = allele;
    entry.indel.at = (int64_t)prefix;
    entry.indel.deleted = (uint32_t)deleted;
    entry.indel.inserted = spelled.substr(prefix, inserted);
    left_align(entry.indel, 1, base);
    entry.after = tract_end(entry.indel, (int64_t)walk_spelling.size(), base);
    if (entry.after >= (int64_t)walk_spelling.size()) {
        // No read can cover a walk base after the repeat tract, which `observe_indel` needs, so
        // leave the allele scored as its step sequences.
        return;
    }
    changed_lo = std::min(changed_lo, entry.indel.at);
    changed_hi = std::max(changed_hi, entry.after);
    // Take the rate of a large indel as that of an error the size of its repeat tract, not of its
    // length.
    entry.large = entry.indel.length() >= LARGE_INDEL_MIN_LENGTH;
    const auto tract = (double)std::max<int64_t>(
        1, entry.after - entry.indel.at - (entry.large ? (int64_t)deleted : 0));
    const double insertion = std::min(0.5, params.indel_insertion_rate * tract);
    const double deletion = std::min(0.5, params.indel_deletion_rate * tract);
    entry.ln_shown = log(deleted > 0 ? deletion : insertion);
    entry.ln_hidden = log(deleted > 0 ? insertion : deletion);
    // The allele is scored as its walk, and measured as its walk too. A read must cover the tract
    // on either allele to tell them apart, so neither has bases of its own that only its reads
    // cover, and the mixture weights are even.
    is_indel_allele[allele] = true;
    indel_alleles.push_back(std::move(entry));
}

void GraphAlignedAlleleLikelihoodCalculator::WalkChanges::cut_windows(
    const HandleGraph& graph, const vector<Traversal>& traversals,
    const AlleleSequences* sequences) {
    indel_alleles.clear();
    changed_hi = std::min(changed_hi, (int64_t)walk_spelling.size() - 1);
    for (size_t a = 0; a < traversals.size(); ++a) {
        const string spelled = allele_sequence(graph, traversals[a], own_sequences(sequences, a));
        // Every allele agrees with the walk outside the window.
        window.push_back(spelled.substr(
            changed_lo, spelled.size() - (walk_spelling.size() - changed_hi) - changed_lo));
    }
}

int64_t GraphAlignedAlleleLikelihoodCalculator::WalkChanges::pin_at(
    const HandleGraph& graph, const Traversal& walk, int64_t point, int64_t base,
    handle_t& handle) const {
    const size_t step = std::ranges::upper_bound(walk_starts, base) - walk_starts.begin() - 1;
    handle = walk[step];
    const int64_t k = point - walk_starts[step];
    return graph.get_is_reverse(handle) ? (int64_t)graph.get_length(handle) - k : k;
}

bool GraphAlignedAlleleLikelihoodCalculator::WalkChanges::follow_walk(
    const vector<ReadStep>& read_steps, vector<const Mapping*>& mappings,
    vector<int64_t>& starts, size_t& first, size_t& last) const {
    for (size_t j = 0; j < read_steps.size(); ++j) {
        auto found = walk_index.find(step_key(read_steps[j].node_id, read_steps[j].backward));
        if (found == walk_index.end() || (j > 0 && found->second != last + 1)) {
            return false;
        }
        if (j == 0) {
            first = found->second;
        }
        last = found->second;
        mappings.push_back(read_steps[j].mapping);
        starts.push_back(walk_starts[found->second]);
    }
    return !mappings.empty();
}

void GraphAlignedAlleleLikelihoodCalculator::WalkChanges::charge_window(
    const WalkAlignment& walked, const string& sequence, const EditAlignmentScorer& plain_scorer,
    vector<double>& row) const {
    const int64_t r0 = walked.read_offset_at(changed_lo - 1);
    const int64_t r1 = walked.read_offset_at(changed_hi);
    if (r0 < 0 || r1 <= r0) {
        // The read does not cover the window and a base on each side.
        return;
    }
    string bases = sequence.substr(r0 + 1, r1 - r0 - 1);
    for (char& b : bases) {
        b = (char)toupper(b);
    }
    // Take the cost from the read's score on the walk, which every allele shares here.
    for (size_t b = 0; b < row.size(); ++b) {
        row[b] -= plain_scorer.get_log_base() * (double)edit_cost(bases, window[b], plain_scorer);
    }
}

void GraphAlignedAlleleLikelihoodCalculator::WalkChanges::observe_indels(
    const WalkAlignment& walked, int64_t reach_begin, int64_t reach_end,
    const EditAlignmentScorer& read_scorer, double insertion_gap_nats,
    vector<double>& row) const {
    const double log_base = read_scorer.get_log_base();
    for (const IndelAllele& entry : indel_alleles) {
        if (!(row[entry.allele] > -numeric_limits<double>::infinity())) {
            continue;
        }
        int clip = 0;
        const IndelObservation seen =
            entry.large ? observe_large_indel(walked, entry.indel, entry.after, LargeIndelParams(),
                                              reach_begin, reach_end, &clip)
                        : observe_indel(walked, entry.indel, entry.after);
        if (seen == IndelObservation::ALT) {
            // The allele explains the read without the gap of its indel, and the walk explains
            // it only by an error that shows the indel. A read clipped at the end of the indel
            // paid no gap.
            double fit = row[entry.allele];
            if (clip == 0) {
                const double gap_nats = entry.indel.deleted == 0 ? insertion_gap_nats : 0.0;
                fit = row[entry.allele]
                      - log_base * (double)read_scorer.score_gap(entry.indel.length())
                      - gap_nats;
            }
            row[entry.allele] = fit;
            for (size_t b = 0; b < row.size(); ++b) {
                if (!is_indel_allele[b] && row[b] > -numeric_limits<double>::infinity()) {
                    row[b] = fit + entry.ln_shown;
                }
            }
        } else if (seen == IndelObservation::REF) {
            // The allele explains the read only by an error that hides the indel.
            row[entry.allele] += entry.ln_hidden;
        }
    }
}

void GraphAlignedAlleleLikelihoodCalculator::WalkChanges::adjust_row(
    const vector<ReadStep>& read_steps, const string& sequence,
    const EditAlignmentScorer& read_scorer, const EditAlignmentScorer& plain_scorer,
    double insertion_gap_nats, vector<double>& row) const {
    vector<const Mapping*> mappings;
    vector<int64_t> starts;
    size_t first = 0;
    size_t last = 0;
    if (!follow_walk(read_steps, mappings, starts, first, last)) {
        return;
    }
    WalkAlignment walked = walk_alignment(mappings, starts, read_steps.front().read_offset,
                                          sequence);
    left_align(walked, base_lookup());
    if (window_scored) {
        charge_window(walked, sequence, plain_scorer, row);
    }
    // The alignment's clips are real only at the ends of the read itself. Where the read goes on
    // past an end of the walk, it covers everything beyond that end.
    const ReadStep& final_step = read_steps.back();
    const bool starts_read = read_steps.front().read_offset == 0;
    const bool ends_read = final_step.read_offset + final_step.read_length >= sequence.size();
    if (!starts_read) {
        walked.clip_begin = 0;
    }
    if (!ends_read) {
        walked.clip_end = 0;
    }
    const int64_t reach_begin = !starts_read && first == 0 && walked.begin == 0
                                    ? numeric_limits<int64_t>::min() / 4
                                    : walked.begin;
    const int64_t reach_end = !ends_read && last + 1 == walk_starts.size() - 1
                                      && walked.end == (int64_t)walk_spelling.size()
                                  ? numeric_limits<int64_t>::max() / 4
                                  : walked.end;
    observe_indels(walked, reach_begin, reach_end, read_scorer, insertion_gap_nats, row);
}

/**
 * The two pins of a site, where the anchors place each read in the site's bounds.
 *
 * Each pin is an offset along the forward strand of its node. The start pin follows the start
 * node, and the end pin precedes the end node, in the site's direction. A site inside nodes puts
 * its start pin after the walk base before the bases its alleles change. It puts its end pin
 * before the base after them, or after the last of them where they end the walk. Where they
 * change no bases, its two pins are one.
 */
struct GraphAlignedAlleleLikelihoodCalculator::SitePins {
    handle_t start_handle;
    handle_t end_handle;
    int64_t start_pin;
    int64_t end_pin;
    /// Whether the end pin is resolved as an exit pin, against the base before it.
    ///
    /// This is so for a site inside nodes whose changed bases end the walk, leaving no base after
    /// them.
    bool end_exit;
    bool inside_node;

    /// Place the pins of `site`, whose alleles make the changes `changes`.
    SitePins(const HandleGraph& graph, const MaterializedSite& site, const WalkChanges& changes);

    /// Describe the pins in the anchor evidence of the site.
    void describe(const HandleGraph& graph, AnchorSiteEvidence& evidence) const;

    /// Return the anchor record of a read, with its mismapping probability `mismap`.
    AnchorRead anchor_read(const SiteRead& read, const HandleGraph& graph, double mismap,
                           AnchorCounters& counters) const;
};

GraphAlignedAlleleLikelihoodCalculator::SitePins::SitePins(const HandleGraph& graph,
                                                            const MaterializedSite& site,
                                                            const WalkChanges& changes)
    : start_handle(site.start), end_handle(site.end),
      start_pin(graph.get_is_reverse(site.start) ? 0 : (int64_t)graph.get_length(site.start)),
      end_pin(graph.get_is_reverse(site.end) ? (int64_t)graph.get_length(site.end) : 0),
      end_exit(site.bounds_inside_nodes(graph)
               && changes.changed_hi == (int64_t)changes.walk_spelling.size()),
      inside_node(site.bounds_inside_nodes(graph)) {
    if (!inside_node || changes.changed_lo > changes.changed_hi) {
        return;
    }
    const Traversal& walk = site.alleles.walks.front();
    const int64_t lo = changes.changed_lo;
    const int64_t hi = changes.changed_hi;
    start_pin = changes.pin_at(graph, walk, lo, std::max<int64_t>(lo - 1, 0), start_handle);
    end_handle = start_handle;
    if (lo != hi) {
        end_pin = changes.pin_at(graph, walk, hi, hi - (end_exit ? 1 : 0), end_handle);
    } else {
        end_pin = start_pin;
    }
}

void GraphAlignedAlleleLikelihoodCalculator::SitePins::describe(
    const HandleGraph& graph, AnchorSiteEvidence& evidence) const {
    evidence.start_node = graph.get_id(start_handle);
    evidence.end_node = graph.get_id(end_handle);
    evidence.inside_node = inside_node;
    evidence.start_pin = (uint32_t)start_pin;
    evidence.end_pin = (uint32_t)end_pin;
}

AnchorRead GraphAlignedAlleleLikelihoodCalculator::SitePins::anchor_read(
    const SiteRead& read, const HandleGraph& graph, double mismap,
    AnchorCounters& counters) const {
    AnchorRead record;
    record.read = read_names().intern(read.aln->name());
    record.mismap = (float)mismap;
    record.start_pin = resolve_anchor_pin(read, graph, graph.get_id(start_handle),
                                          graph.get_is_reverse(start_handle), true, counters,
                                          inside_node ? start_pin : -1);
    record.end_pin = resolve_anchor_pin(read, graph, graph.get_id(end_handle),
                                        graph.get_is_reverse(end_handle), end_exit, counters,
                                        inside_node ? end_pin : -1);
    return record;
}

/**
 * A table of the alleles of one site, as node visits with sequences, ready for scoring reads.
 */
struct GraphAlignedAlleleLikelihoodCalculator::SiteAlleles {
    /// The visits of each allele.
    vector<vector<AlleleStep>> steps;
    /// The sorted step keys of each allele.
    vector<vector<int64_t>> keys;
    /// Where each step key occurs in each allele, for greedy pairing only.
    vector<AlleleStepPositions> positions;
    /// The orientation in which every allele visits each node, for the nodes where they agree.
    ///
    /// A read aligned to the opposite strand is flipped into the alleles' reading direction by a
    /// vote of these nodes.
    unordered_map<nid_t, bool> orientations;

    /// Fill the keys, positions and orientations from `steps`.
    void index(bool greedy);

    /// Return the length that each allele spells.
    vector<size_t> spelled_lengths() const;

    /**
     * Return the length of each allele without the site's boundary nodes, T_a in the depth term.
     *
     * A read that lies inside one boundary node yields no row, so boundary sequence yields no
     * rows. An allele that deletes the whole interior gets length 0, and lambda then counts only
     * the R - 1 reads that span its junction.
     */
    vector<size_t> interior_lengths(const unordered_set<nid_t>& boundary_nodes) const;

    /**
     * Return the number of bases of each allele that each other allele lacks.
     *
     * A node visited more than once by an allele counts once. A node both alleles visit, but
     * spelled differently by an allele's own sequence, counts the bases outside the longest
     * prefix and suffix the two spellings share.
     */
    vector<vector<size_t>> unique_lengths() const;
};

void GraphAlignedAlleleLikelihoodCalculator::SiteAlleles::index(bool greedy) {
    keys.reserve(steps.size());
    for (const auto& allele : steps) {
        keys.push_back(sorted_allele_keys(allele));
    }
    if (greedy) {
        positions.reserve(steps.size());
        for (const auto& allele : steps) {
            positions.push_back(index_allele_steps(allele));
        }
    }
    unordered_set<nid_t> ambiguous;
    for (const auto& allele : steps) {
        for (const AlleleStep& step : allele) {
            auto found = orientations.find(step.node_id);
            if (found == orientations.end()) {
                orientations[step.node_id] = step.backward;
            } else if (found->second != step.backward) {
                ambiguous.insert(step.node_id);
            }
        }
    }
    for (nid_t node_id : ambiguous) {
        orientations.erase(node_id);
    }
}

vector<size_t> GraphAlignedAlleleLikelihoodCalculator::SiteAlleles::spelled_lengths() const {
    vector<size_t> lengths;
    lengths.reserve(steps.size());
    for (const auto& allele : steps) {
        lengths.push_back(spelled_length(allele, 0, allele.size()));
    }
    return lengths;
}

vector<size_t> GraphAlignedAlleleLikelihoodCalculator::SiteAlleles::interior_lengths(
    const unordered_set<nid_t>& boundary_nodes) const {
    vector<size_t> lengths;
    lengths.reserve(steps.size());
    for (const auto& allele : steps) {
        size_t len = 0;
        for (const AlleleStep& step : allele) {
            if (!boundary_nodes.contains(step.node_id)) {
                len += step.sequence.size();
            }
        }
        lengths.push_back(len);
    }
    return lengths;
}

vector<vector<size_t>> GraphAlignedAlleleLikelihoodCalculator::SiteAlleles::unique_lengths() const {
    vector<map<nid_t, const AlleleStep*>> content(steps.size());
    for (size_t a = 0; a < steps.size(); ++a) {
        for (const AlleleStep& step : steps[a]) {
            content[a][step.node_id] = &step;
        }
    }
    vector lengths(steps.size(), vector<size_t>(steps.size(), 0));
    for (size_t a = 0; a < content.size(); ++a) {
        for (size_t b = 0; b < content.size(); ++b) {
            if (a == b) {
                continue;
            }
            size_t total = 0;
            for (const auto& [node_id, step] : content[a]) {
                auto other = content[b].find(node_id);
                if (other == content[b].end()) {
                    total += step->sequence.size();
                } else if ((step->own_sequence || other->second->own_sequence)
                           && step->sequence != other->second->sequence) {
                    total += bases_not_shared(step->sequence, other->second->sequence);
                }
            }
            lengths[a][b] = total;
        }
    }
    return lengths;
}

const Alignment* GraphAlignedAlleleLikelihoodCalculator::orient_read(
    const SiteRead& read, const unordered_set<nid_t>& site_nodes,
    const unordered_set<nid_t>& boundary_nodes, const SiteAlleles& alleles,
    vector<ReadStep>& read_steps, Alignment& flipped) const {
    if (!get_read_steps(read, site_nodes, boundary_nodes, read_steps)) {
        return nullptr;
    }
    if (!read_is_reverse_of_alleles(read_steps, alleles.orientations)) {
        return read.aln;
    }
    // Flip only the stretch of the read inside the site.
    function<int64_t(nid_t)> node_length = [&](nid_t node_id) {
        return (int64_t)graph.get_length(graph.get_handle(node_id));
    };
    flipped = reverse_complement_alignment(site_span_of(read, site_nodes), node_length);
    // The stretch holds only the site's mappings, so it needs no index.
    SiteRead flipped_read;
    flipped_read.aln = &flipped;
    if (!get_read_steps(flipped_read, site_nodes, boundary_nodes, read_steps)) {
        // Flipping keeps the nodes the read visits, so this cannot happen.
        return nullptr;
    }
    return &flipped;
}

void GraphAlignedAlleleLikelihoodCalculator::score_row(
    const Alignment& aln, const vector<ReadStep>& read_steps,
    const EditAlignmentScorer& read_scorer, const SiteAlleles& alleles, ReadScratch& scratch,
    vector<double>& row) const {
    const double log_base = read_scorer.get_log_base();
    if (params.optimal_pairing) {
        prepare_read_scratch(aln, read_steps, read_scorer, scratch);
    }
    for (size_t a = 0; a < alleles.steps.size(); ++a) {
        const vector<AlleleStep>& steps = alleles.steps[a];
        if (steps.empty()) {
            // An allele with no visits, such as the empty traversal of a star allele, cannot place
            // the read.
            row[a] = -numeric_limits<double>::infinity();
            continue;
        }
        double nat_adjust = 0.0;
        const int32_t score =
            params.optimal_pairing
                ? score_by_optimal_pairing(aln, read_steps, steps, scratch, alleles.keys[a],
                                           read_scorer, nat_adjust)
                : score_by_greedy_pairing(aln, read_steps, steps, alleles.positions[a],
                                          read_scorer, nat_adjust);
        row[a] = log_base * (double)score + nat_adjust;
    }
}

/// Collect the nodes of `site`, including its boundary nodes, and its boundary nodes alone.
///
/// A site inside nodes is made of the nodes its alleles walk, and has no boundary node. Its bounds
/// are inside its first and last nodes, so every read on them may tell its alleles apart, and
/// their bases all count towards the depth term.
static void collect_site_nodes(const HandleGraph& graph, const MaterializedSite& site,
                               unordered_set<nid_t>& site_nodes,
                               unordered_set<nid_t>& boundary_nodes) {
    if (site.bounds_inside_nodes(graph)) {
        site_nodes.insert(graph.get_id(site.start));
        for (const Traversal& walk : site.alleles.walks) {
            for (const handle_t& h : walk) {
                site_nodes.insert(graph.get_id(h));
            }
        }
    } else {
        const auto [nodes, edges] = site_contents(graph, site.start, site.end, true);
        site_nodes.insert(nodes.begin(), nodes.end());
        boundary_nodes = {graph.get_id(site.start), graph.get_id(site.end)};
    }
}

/// Return node IDs merged into ranges of consecutive IDs, in ascending order.
static vector<pair<nid_t, nid_t>> id_ranges(const unordered_set<nid_t>& ids) {
    vector<nid_t> sorted_ids(ids.begin(), ids.end());
    std::ranges::sort(sorted_ids);
    vector<pair<nid_t, nid_t>> ranges;
    for (nid_t id : sorted_ids) {
        if (!ranges.empty() && id == ranges.back().second + 1) {
            ranges.back().second = id;
        } else {
            ranges.emplace_back(id, id);
        }
    }
    return ranges;
}

/// Put `kept`, which holds an entry for each read that the builder kept, in the order of the
/// rows of the matrix it built.
template<typename T>
static void put_in_row_order(vector<T>& kept, const vector<size_t>& order) {
    vector<T> reordered;
    reordered.reserve(order.size());
    for (size_t i : order) {
        reordered.push_back(std::move(kept[i]));
    }
    kept = std::move(reordered);
}

/// Make the read-phasing evidence of a site, before any read is added.
static unique_ptr<PhaseReadEvidence> new_phase_evidence(const vector<size_t>& spelled_lengths,
                                                        bool length_weighted) {
    auto evidence = make_unique<PhaseReadEvidence>();
    evidence->n_alleles = spelled_lengths.size();
    evidence->length_weighted = length_weighted;
    for (size_t length : spelled_lengths) {
        evidence->allele_length.push_back((uint32_t)length);
    }
    return evidence;
}

/// Make the anchor evidence of a site, before any read or pin is added.
static unique_ptr<AnchorSiteEvidence> new_anchor_evidence(const vector<size_t>& spelled_lengths,
                                                          bool length_weighted) {
    auto evidence = make_unique<AnchorSiteEvidence>();
    evidence->n_alleles = spelled_lengths.size();
    evidence->length_weighted = length_weighted;
    for (size_t length : spelled_lengths) {
        evidence->allele_length.push_back((uint32_t)length);
    }
    return evidence;
}

/// Copy the rows of `matrix` into the anchor evidence of its site, whose reads are in row order.
static void copy_rows(const AlleleReadLikelihoods& matrix, AnchorSiteEvidence& evidence) {
    evidence.mean_read_length = (float)matrix.mean_read_length_estimate();
    evidence.rel.resize(matrix.num_reads() * matrix.num_alleles());
    for (size_t r = 0; r < matrix.num_reads(); ++r) {
        evidence.reads[r].mismap = (float)matrix.mismap_prob(r);
        for (size_t a = 0; a < matrix.num_alleles(); ++a) {
            evidence.rel[r * matrix.num_alleles() + a] = (float)matrix.rel(r, a);
        }
    }
}

/// Copy the rows of `matrix` into the phasing evidence of its site, whose keys are in row order.
static void copy_rows(const AlleleReadLikelihoods& matrix, PhaseReadEvidence& evidence) {
    evidence.mean_read_length = (float)matrix.mean_read_length_estimate();
    evidence.mismap.resize(matrix.num_reads());
    evidence.rel.resize(matrix.num_reads() * matrix.num_alleles());
    for (size_t r = 0; r < matrix.num_reads(); ++r) {
        evidence.mismap[r] = (float)matrix.mismap_prob(r);
        for (size_t a = 0; a < matrix.num_alleles(); ++a) {
            evidence.rel[r * matrix.num_alleles() + a] = (float)matrix.rel(r, a);
        }
    }
}

void AlleleReadLikelihoods::attach_anchor_evidence(unique_ptr<AnchorSiteEvidence> evidence) {
    copy_rows(*this, *evidence);
    anchor_evidence = std::move(evidence);
}

void AlleleReadLikelihoods::attach_phase_evidence(unique_ptr<PhaseReadEvidence> evidence) {
    copy_rows(*this, *evidence);
    phase_evidence = std::move(evidence);
}

/**
 * Give `matrix` the evidence collected alongside its rows, in its row order.
 *
 * @param order The row order of the builder that made `matrix`.
 * @param anchor The anchor evidence, or null.
 * @param phase The read-phasing evidence, or null.
 */
static void attach_evidence(AlleleReadLikelihoods& matrix, const vector<size_t>& order,
                            unique_ptr<AnchorSiteEvidence> anchor,
                            unique_ptr<PhaseReadEvidence> phase) {
    if (anchor != nullptr) {
        if (anchor->reads.size() == order.size()) {
            put_in_row_order(anchor->reads, order);
        }
        matrix.attach_anchor_evidence(std::move(anchor));
    }
    if (phase != nullptr && phase->read_key.size() == order.size()) {
        put_in_row_order(phase->read_key, order);
    }
    if (phase != nullptr && phase->read_key.size() == matrix.num_reads()) {
        matrix.attach_phase_evidence(std::move(phase));
    }
    // Otherwise drop the phasing evidence. Padding the keys would give every padded read the same
    // key, so they would be taken for one fragment at every site.
}

double GraphAlignedAlleleLikelihoodCalculator::mismap_probability(const Alignment& aln) const {
    // vg's mappers build MAPQ from distinct graph placements, so it measures whether the read
    // could be at another locus, rather than whether another allele of this site fits too. That
    // is the meaning the mismapping term needs.
    return params.use_mismap_term ? phred_to_prob((double)aln.mapping_quality())
                                  : params.min_mismap_prob;
}

AlleleReadLikelihoods GraphAlignedAlleleLikelihoodCalculator::compute(
    const MaterializedSite& site, const vector<nid_t>& fallback_nodes, int region_ploidy) {
    const vector<Traversal>& traversals = site.alleles.walks;

    AlleleReadLikelihoodsBuilder builder(traversals.size(), params.min_mismap_prob,
                                        params.max_mismap_prob);
    if (traversals.empty()) {
        return builder.build();
    }

    unordered_set<nid_t> site_nodes;
    unordered_set<nid_t> boundary_nodes;
    collect_site_nodes(graph, site, site_nodes, boundary_nodes);
    if (site_nodes.empty()) {
        return builder.build();
    }
    const vector<pair<nid_t, nid_t>> ranges = id_ranges(site_nodes);

    const WalkChanges changes(graph, site, params);
    SiteAlleles alleles;
    alleles.steps.reserve(traversals.size());
    for (size_t a = 0; a < traversals.size(); ++a) {
        const StepSequences own = changes.scored_as_walk(a) ? StepSequences()
                                                            : site.alleles.step_sequences(graph, a);
        alleles.steps.push_back(get_allele_steps(traversals[a], own.empty() ? nullptr : &own));
    }
    alleles.index(!params.optimal_pairing);

    if (params.length_weighted_mixture) {
        builder.set_allele_lengths(alleles.spelled_lengths());
        builder.set_unique_lengths(alleles.unique_lengths());
    }
    // These lengths ask how much sequence yields reads, while the unique lengths of the mixture
    // weights ask which sequence tells the alleles apart.
    const vector<size_t> depth_lengths = alleles.interior_lengths(boundary_nodes);
    const SitePins pins(graph, site, changes);

    // Collect the anchor evidence, if anchors are being written, with the reads in the order
    // `add_read` kept them. Give it the spelled lengths whatever the mixture setting, because the
    // anchors weight their own scores.
    unique_ptr<AnchorSiteEvidence> anchor_evidence;
    if (params.collect_anchors) {
        anchor_evidence = new_anchor_evidence(alleles.spelled_lengths(),
                                              params.length_weighted_mixture);
        pins.describe(graph, *anchor_evidence);
    }
    // Collect the `rel` rows for read phasing, if it is on and anchors are not being written,
    // since the anchor evidence already holds the rows.
    unique_ptr<PhaseReadEvidence> phase_evidence;
    if (params.collect_read_phasing && !params.collect_anchors) {
        phase_evidence = new_phase_evidence(alleles.spelled_lengths(),
                                            params.length_weighted_mixture);
    }

    // Count anchor pin resolutions into the run's counters where there are some, and into this
    // local otherwise, as when a unit test drives the calculator directly.
    AnchorCounters unowned_counters;
    AnchorCounters& anchor_pin_counters =
        params.anchor_counters != nullptr ? *params.anchor_counters : unowned_counters;

    vector<ReadStep> read_steps;
    ReadScratch read_scratch;
    vector<double> row(traversals.size());

    read_source.for_each_read(ranges, [&](const SiteRead& read) {
        Alignment flipped;
        const Alignment* scored_aln = orient_read(read, site_nodes, boundary_nodes, alleles,
                                                  read_steps, flipped);
        if (scored_aln == nullptr) {
            return;
        }
        // Score a read with no base qualities with the scorer that does not need them.
        const EditAlignmentScorer& read_scorer =
            scored_aln->quality().empty() ? plain_scorer : qual_scorer;
        score_row(*scored_aln, read_steps, read_scorer, alleles, read_scratch, row);
        if (changes.changes_reads()) {
            changes.adjust_row(read_steps, scored_aln->sequence(), read_scorer, plain_scorer,
                               params.insertion_gap_nats, row);
        }

        const Alignment& aln = *read.aln;
        const double mismap = mismap_probability(aln);
        const Position* start =
            aln.path().mapping_size() > 0 ? &aln.path().mapping(0).position() : nullptr;
        if (!builder.add_read(row, mismap, aln.name(), aln.sequence().size(), start)) {
            return;
        }

        if (phase_evidence != nullptr) {
            // Keep only the hash of the name, because the names make the anchor evidence too
            // heavy to keep for phasing.
            phase_evidence->read_key.push_back((uint64_t)std::hash<string>{}(aln.name()));
        }
        if (anchor_evidence != nullptr) {
            // Resolve the pins against the read's own alignment, whose offsets are in the read
            // file. The record's strand field gives the orientation.
            anchor_evidence->reads.push_back(
                pins.anchor_read(read, graph, mismap, anchor_pin_counters));
        }
    });

    AlleleReadLikelihoods result = builder.build();

    // Set R, the mean read length, for the depth term's lambda, the mixture weights and the anchor
    // slot weights. The builder's R is the mean over this site's reads, which over-represents long
    // reads, because a long read overlaps more sites. The rate window's mean counts each read once,
    // where it begins.
    WindowReadStats stats = local_read_stats(site, fallback_nodes, ranges);
    if (stats.mean_read_length > 0.0) {
        result.set_mean_read_length(stats.mean_read_length);
    }

    attach_evidence(result, builder.row_order(), std::move(anchor_evidence),
                    std::move(phase_evidence));
    if (!depth_lengths.empty()) {
        // Set the depth context whether or not the depth term is on, so that DR is always
        // written. The rate is per haplotype, and the window counts reads from every haplotype
        // in the region, not only those crossing the site. Where no read began in the window, as
        // in a small window at the end of a contig, take R from the site's own reads.
        int haplotypes = region_ploidy > 0 ? region_ploidy : params.depth_ploidy;
        result.set_depth_context(depth_lengths,
                                 stats.start_rate / (double)haplotypes,
                                 stats.mean_read_length > 0.0 ? stats.mean_read_length
                                                              : result.mean_read_length_estimate(),
                                 params.depth_weight, params.depth_effective_reads);
    }
    return result;
}

}

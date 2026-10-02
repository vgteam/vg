#include "allele_likelihood.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "alignment.hpp"
#include "path.hpp"
#include "statistics.hpp"

namespace vg {

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// AlleleReadLikelihoods
////////////////////////////////////////////////////////////////////////////////

void AlleleReadLikelihoods::set_contents(size_t n_reads, size_t n_alleles, vector<double>&& matrix,
                                         vector<double>&& mismap, vector<double>&& best_ln,
                                         vector<string>&& names, size_t unplaceable) {
    this->n_reads = n_reads;
    this->n_alleles = n_alleles;
    this->matrix = std::move(matrix);
    this->read_mismap_prob = std::move(mismap);
    this->read_best_ln = std::move(best_ln);
    this->read_names = std::move(names);
    this->unplaceable = unplaceable;
    // sum_r (1 - e_r), the expected number of these reads that came from this locus.
    // Cached because the depth term asks for it once per genotype.
    this->effective_read_total = 0.0;
    for (double e : this->read_mismap_prob) {
        this->effective_read_total += 1.0 - e;
    }
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

double AlleleReadLikelihoods::expected_reads(const vector<int>& genotype) const {
    return expected_reads_from(traversal_lengths, depth_rate, depth_read_length, genotype);
}

double AlleleReadLikelihoods::expected_reads_from(const vector<size_t>& traversal_lengths,
                                                  double rate, double read_length,
                                                  const vector<int>& genotype) {
    double total = 0.0;
    for (int allele : genotype) {
        double len = (allele >= 0 && (size_t)allele < traversal_lengths.size())
                         ? (double)traversal_lengths[allele] : 0.0;
        total += len + read_length - 1.0;
    }
    return rate * total;
}

double AlleleReadLikelihoods::observed_reads() const {
    return depth_effective ? effective_read_total : (double)n_reads;
}

double AlleleReadLikelihoods::depth_ratio(const vector<int>& genotype) const {
    double expected = expected_reads(genotype);
    return expected > 0.0 ? observed_reads() / expected : -1.0;
}

/// n ln(lambda) - lambda - ln Gamma(n + 1), with n the observed read count and lambda
/// the expected one. The observation is `sum_r (1 - e_r)` rather than a row count, so n
/// is usually not a whole number, and ln Gamma(n + 1) stands in for ln n!. As a function
/// of n >= 0 this is, up to a normalising constant, the log density of a continuous
/// analogue of the Poisson distribution over read counts. Here it is evaluated at the
/// observed n as a log likelihood of lambda. The depth term raises it to the power beta and
/// divides by the normaliser, ln_tempered_poisson_normaliser. The ln Gamma(n + 1) term does
/// not depend on lambda, so it is the same for every genotype at a site.
static double ln_poisson_pmf(double n, double lambda) {
    if (lambda <= 0.0) {
        return -numeric_limits<double>::infinity();
    }
    if (n <= 0.0) {
        return -lambda;
    }
    return n * log(lambda) - lambda - lgamma(n + 1.0);
}

/// ln Z_beta(lambda), where Z_beta(lambda) is the integral over n >= 0 of f(n; lambda)^beta and
/// f is the density of ln_poisson_pmf, by the normal approximation. For large lambda, f is close
/// to a normal density with mean and variance lambda, so f^beta is (2 pi lambda)^(-beta/2) times
/// a normal density with variance lambda / beta, scaled by sqrt(2 pi lambda / beta), and
///
///     ln Z_beta(lambda) ~= ((1 - beta) / 2) ln(2 pi lambda) - (1/2) ln beta.
///
/// Against numerical integration at beta = 0.1, 0.5 and 1, it is within 0.13 nats for every
/// lambda >= 2 and within 0.02 for lambda >= 30. Below lambda = 2 the approximation falls away
/// from the integral, so lambda is held at 2 there.
static double ln_tempered_poisson_normaliser(double lambda, double beta) {
    const double held = max(lambda, 2.0);
    return 0.5 * (1.0 - beta) * log(2.0 * M_PI * held) - 0.5 * log(beta);
}

vector<double> AlleleReadLikelihoods::mixture_weights(const vector<int>& genotype) const {
    if (genotype.empty()) {
        return {};
    }
    double flat = 1.0 / (double)genotype.size();

    // Expected share of this site's reads for each haplotype of the genotype. Flat
    // 1/|G| unless lengths were supplied; see set_length_weights.
    vector<double> weights(genotype.size(), flat);
    if (uses_length_weights()) {
        double sum = 0.0;
        for (size_t i = 0; i < genotype.size(); ++i) {
            int allele = genotype[i];
            double own;
            if (!unique_lengths.empty() && allele >= 0
                && (size_t)allele < unique_lengths.size()) {
                // Sequence in this allele that no other member of the genotype
                // carries. Reads outside it cannot separate the two.
                size_t u = numeric_limits<size_t>::max();
                for (size_t j = 0; j < genotype.size(); ++j) {
                    int other = genotype[j];
                    if (j == i || other < 0
                        || (size_t)other >= unique_lengths[allele].size()) {
                        continue;
                    }
                    u = min(u, unique_lengths[allele][other]);
                }
                own = (u == numeric_limits<size_t>::max()) ? 0.0 : (double)u;
            } else if (allele >= 0 && (size_t)allele < allele_lengths.size()) {
                own = (double)allele_lengths[allele];
            } else {
                own = 0.0;
            }
            weights[i] = own + mean_read_length - 1.0;
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
            // The VCF layer uses negative sentinels for star and missing
            // alleles, so be defensive rather than reading out of bounds.
            if (allele < 0 || (size_t)allele >= n_alleles) {
                continue;
            }
            mixture += weights[i] * rel(r, (size_t)allele);
        }

        // Mix in the chance that the read is mismapped, in which case its likelihood
        // is the row maximum, 1. The result lies in [e_r, 1], so its log is finite.
        double e_r = read_mismap_prob[r];
        total += log((1.0 - e_r) * mixture + e_r);
    }

    if (uses_depth_term()) {
        // ln of the tempered density f^beta / Z_beta at the observed count, with beta the
        // depth weight. Z_beta grows with lambda when beta < 1, so leaving it out would favour
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

    // Per distinct haplotype slot of `called`: what share of an ideal pileup it
    // contributes, and what mixture each genotype would then show. An ideal read from
    // allele a has rel 1 for a and 0 elsewhere, so the mixture collapses to the total
    // weight the genotype places on a -- which is why a homozygote sees 1 and a
    // heterozygote sees its own strand's weight.
    struct Slot { double share, mix_called, mix_runner; };
    vector<Slot> slots;
    slots.reserve(called.size());
    for (size_t i = 0; i < called.size(); ++i) {
        int a = called[i];
        if (a < 0 || (size_t)a >= n_alleles) {
            // Star and missing alleles carry no sequence to discriminate on.
            continue;
        }
        double mix_called = 0.0, mix_runner = 0.0;
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

    // An ideal read is well mapped, so every ideal read has e_r at the floor, whatever
    // e_r the site's own reads have. With their own e_r, a poorly mapped site would get a
    // small denominator, and a weak call there would look strong.
    double e = mismap_floor;
    double per_read = 0.0;
    for (const Slot& s : slots) {
        per_read += s.share * (log((1.0 - e) * s.mix_called + e)
                               - log((1.0 - e) * s.mix_runner + e));
    }
    double total = per_read * (double)n_reads;

    // Identical genotypes give exactly 0, and floating point can drift a hair below.
    return max(total, 0.0);
}

vector<vector<int>> AlleleReadLikelihoods::enumerate_genotypes(size_t num_alleles, int ploidy) {
    vector<vector<int>> genotypes;
    if (num_alleles == 0 || ploidy <= 0) {
        return genotypes;
    }

    // VCF orders genotypes so that a genotype with alleles a_1 <= ... <= a_P sits
    // at an index given by a recursion over its largest allele. Generating
    // non-decreasing tuples in colexicographic order reproduces that ordering, so a
    // genotype's position here is its GL index. For diploid this is (0,0), (0,1),
    // (1,1), (0,2), (1,2), (2,2), ...
    vector<int> current(ploidy, 0);

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

    // The row's divisor is the read's best fit over all alleles at the site, so it
    // does not depend on the genotype and cancels from every genotype comparison.
    double best = -numeric_limits<double>::infinity();
    for (double ll : raw_ln_likelihood) {
        best = max(best, ll);
    }

    if (!(best > -numeric_limits<double>::infinity())) {
        // This read placed on no allele, so there is no row maximum to divide by and
        // normalising would give NaN. Drop it and count it.
        ++unplaceable;
        return false;
    }

    rows.emplace_back();
    rows.back().reserve(n_alleles);
    for (double ll : raw_ln_likelihood) {
        // exp of a non-positive number, so in [0,1]; -inf gives exactly 0.
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

    // A canonical row order, so that a floating-point sum over the reads gives the same answer
    // whatever order they arrived in. The name and the start identify an alignment; the row's
    // values break any remaining tie, between rows that are then interchangeable.
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
    vector<double> sorted_mismap, sorted_best;
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
        // Set R whichever mixture is in use, since the depth term's
        // lambda = rate * (L + R - 1) needs it too.
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
    const PathHandleGraph& graph, SnarlManager& snarl_manager, const SiteReadSource& read_source,
    const EditAlignmentScorer& qual_scorer, const EditAlignmentScorer& plain_scorer, const Params& params)
    : graph(graph), snarl_manager(snarl_manager), read_source(read_source),
      qual_scorer(qual_scorer), plain_scorer(plain_scorer), params(params) {
}

vector<GraphAlignedAlleleLikelihoodCalculator::AlleleStep>
GraphAlignedAlleleLikelihoodCalculator::get_allele_steps(const SnarlTraversal& traversal) const {
    vector<AlleleStep> steps;
    steps.reserve(traversal.visit_size());
    for (int64_t i = 0; i < traversal.visit_size(); ++i) {
        const Visit& visit = traversal.visit(i);
        if (visit.node_id() == 0) {
            // A visit to a child snarl rather than a node: the traversal has not
            // been fully expanded, so its sequence cannot be materialised here.
            // Skip it; the flanking nodes still anchor the comparison.
            continue;
        }
        AlleleStep step;
        step.node_id = visit.node_id();
        step.backward = visit.backward();
        step.sequence = graph.get_sequence(graph.get_handle(visit.node_id(), visit.backward()));
        steps.push_back(std::move(step));
    }
    return steps;
}

/// The stretch of an alignment that lies in this site: mappings from the first to the last that
/// touches a site node, with the read sequence and qualities sliced to match.
///
/// A reverse-strand read is flipped once per site it reaches, and flipping copies every mapping,
/// edit and base. Flipping only this stretch keeps the cost proportional to the site rather than
/// to the read, which matters for long reads.
///
/// Flipping the slice gives the same scores as flipping the whole read. A step at read offset
/// `off` with length `len` moves to `to - off - len` in the flipped slice and to
/// `total - off - len` in the flipped read. The two differ by the constant `total - to`, and so do
/// the two flipped sequences, so each step still indexes the same bases.
static Alignment site_span_of(const SiteRead& read, const unordered_set<nid_t>& site_nodes) {
    const Alignment& aln = *read.aln;
    const Path& path = aln.path();
    int64_t first = -1, last = -1;
    size_t from = 0, to = 0;
    if (read.indexed()) {
        // The source has already found the site's mappings, and the read offset before
        // each, so both ends come from those without walking the rest of the read. The
        // queried ranges are the site's nodes, so no second membership test is needed.
        for (size_t k = 0; k < read.mapping_count; ++k) {
            int64_t i = (int64_t)read.mappings[k];
            if (!site_nodes.count(path.mapping(i).position().node_id())) {
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
            size_t to_length = (size_t)mapping_to_length(path.mapping(i));
            if (site_nodes.count(path.mapping(i).position().node_id())) {
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
        // The caller has already established that the read touches the site, so this cannot happen;
        // returning the whole alignment keeps it correct rather than empty if it ever does.
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
        if (!site_nodes.count(node_id)) {
            return;
        }
        ReadStep step;
        step.node_id = node_id;
        step.backward = mapping.position().is_reverse();
        step.read_offset = read_offset;
        step.read_length = (size_t)mapping_to_length(mapping);
        step.mapping = &mapping;
        steps_out.push_back(step);

        if (!boundary_nodes.count(node_id)) {
            touches_interior = true;
        }
    };

    if (read.indexed()) {
        // Only the mappings the source found inside the queried ranges, with their read
        // offsets already summed. Read order is preserved: the index hands them over in
        // ascending mapping order.
        for (size_t k = 0; k < read.mapping_count; ++k) {
            take((int64_t)read.mappings[k], read.read_offsets[read.mappings[k]]);
        }
    } else {
        // Track the read offset across every mapping, including those outside the
        // site: otherwise the offsets of the ones inside would be wrong.
        size_t read_offset = 0;
        for (int64_t i = 0; i < path.mapping_size(); ++i) {
            take(i, read_offset);
            read_offset += (size_t)mapping_to_length(path.mapping(i));
        }
    }

    if (steps_out.empty()) {
        return false;
    }

    // A read can tell the alleles apart if it touches an interior node, or if it moves
    // between two different nodes of the site. The second case includes a read that
    // crosses straight from one boundary node to the other, which supports an allele
    // that deletes the site's interior. A read that stays inside one boundary node fits
    // every allele equally.
    bool uses_internal_edge = false;
    for (size_t i = 1; i < steps_out.size(); ++i) {
        if (steps_out[i].node_id != steps_out[i - 1].node_id) {
            uses_internal_edge = true;
            break;
        }
    }

    if (!touches_interior && !uses_internal_edge) {
        // The read would add the same constant to every allele. A read that fails to
        // place on some allele is different: it is kept, as evidence against that
        // allele.
        steps_out.clear();
        return false;
    }

    return true;
}

bool GraphAlignedAlleleLikelihoodCalculator::read_is_reverse_of_alleles(
    const vector<ReadStep>& read_steps,
    const unordered_map<nid_t, bool>& allele_orientations) const {

    // Alleles all run from the snarl's start to its end. A read aligned to the opposite
    // strand visits the same nodes in the opposite orientation, and must be flipped
    // before its visits can be paired with an allele's. We decide by a vote over the
    // read's visits, so that one node the alleles visit in different orientations
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
    double& nat_adjust) const {

    int32_t score = 0;
    const string& seq = aln.sequence();
    const string& qual = aln.quality();
    size_t read_pos = step.read_offset;

    for (int64_t i = 0; i < step.mapping->edit_size(); ++i) {
        const Edit& edit = step.mapping->edit(i);
        size_t to_len = (size_t)edit.to_length();
        size_t from_len = (size_t)edit.from_length();

        if (read_pos + to_len > seq.size()) {
            // Malformed alignment; refuse to read past the sequence.
            break;
        }

        if (from_len == to_len) {
            auto begin = seq.begin() + read_pos;
            auto end = begin + to_len;
            auto qual_begin = qual.empty() ? seq.begin() : qual.begin() + read_pos;
            if (edit.sequence().empty()) {
                // Match run.
                score += read_scorer.score_exact_match(begin, end, qual_begin);
            } else {
                // Substitution run: each mismatched base charged its own quality.
                score += read_scorer.score_mismatch(begin, end, qual_begin);
            }
        } else if (from_len > to_len) {
            // Deletion relative to the graph.
            score += read_scorer.score_gap(from_len - to_len);
        } else {
            // Insertion relative to the graph: the read carries bases the allele lacks.
            score += read_scorer.score_gap(to_len - from_len);
            nat_adjust += params.insertion_gap_nats;
        }

        read_pos += to_len;
    }

    return score;
}

int32_t GraphAlignedAlleleLikelihoodCalculator::score_substitution(
    const Alignment& aln, size_t read_offset, size_t length, const string& allele_bases,
    size_t allele_offset, const EditAlignmentScorer& read_scorer) const {

    int32_t score = 0;
    const string& seq = aln.sequence();
    const string& qual = aln.quality();

    // Walk the two sequences together, batching consecutive equal bases and
    // consecutive differing bases so each run is scored in one call, with each
    // mismatch charged at its own base quality.
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

// (node, orientation) packed into one integer, so membership is a binary search over a
// sorted vector rather than a hash lookup, which is cheaper for the short traversals of
// most sites.
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
    std::sort(scratch.keys.begin(), scratch.keys.end());
}

vector<int64_t> GraphAlignedAlleleLikelihoodCalculator::sorted_allele_keys(
    const vector<AlleleStep>& allele_steps) {

    vector<int64_t> keys;
    keys.reserve(allele_steps.size());
    for (const AlleleStep& a : allele_steps) {
        keys.push_back(step_key(a.node_id, a.backward));
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

GraphAlignedAlleleLikelihoodCalculator::AlleleStepPositions
GraphAlignedAlleleLikelihoodCalculator::index_allele_steps(
    const vector<AlleleStep>& allele_steps) {
    AlleleStepPositions positions;
    positions.reserve(allele_steps.size() * 2);
    for (size_t j = 0; j < allele_steps.size(); ++j) {
        positions[((int64_t)allele_steps[j].node_id << 1) | (int64_t)allele_steps[j].backward]
            .push_back((uint32_t)j);
    }
    return positions;
}

// Greedy pairing: one left-to-right pass that pairs each read visit with the next occurrence of
// the same visit in the allele, and never revises a pair.
int32_t GraphAlignedAlleleLikelihoodCalculator::score_by_greedy_pairing(
    const Alignment& aln, const vector<ReadStep>& read_steps,
    const vector<AlleleStep>& allele_steps, const AlleleStepPositions& allele_positions,
    const EditAlignmentScorer& read_scorer,
    bool& placed_out, double& nat_adjust) const {

    placed_out = !allele_steps.empty();
    if (!placed_out) {
        return 0;
    }

    int32_t score = 0;
    size_t allele_index = 0;
    bool have_anchor = false;
    size_t bases_accounted = 0;

    // A run of consecutive unpaired read visits after the first same-visit pair is one gap, as
    // in optimal pairing, so that a pairing scores the same whichever search found it: the run's
    // first visit opens the gap and the others extend it. Before that pair, each unpaired read
    // visit is a gap of its own. Each unpaired visit adds --insertion-nats, as in optimal pairing.
    // A visit with no read bases, where the read deletes its whole node, is no gap: it adds
    // nothing and neither starts nor ends a run.
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

    // Last read position visiting each (node, orientation). The pass below has to ask
    // "is the allele's current node still to come in the read?", which is the mirror of
    // the anchor search's "is the read's node still to come in the allele?". Both
    // sequences are topologically ordered, so a later visit means the allele node is not
    // this read node's counterpart.
    unordered_map<int64_t, size_t> read_last_visit;
    read_last_visit.reserve(read_steps.size() * 2);
    for (size_t i = 0; i < read_steps.size(); ++i) {
        read_last_visit[((int64_t)read_steps[i].node_id << 1) | (int64_t)read_steps[i].backward] = i;
    }

    for (size_t read_index = 0; read_index < read_steps.size(); ++read_index) {
        const ReadStep& read_step = read_steps[read_index];
        // Find the first allele step at or after allele_index that makes this read visit.
        // Pairing on shared visits takes the correspondence from the mapper's alignment. The
        // positions are sorted, so lower_bound finds the step without scanning the allele,
        // which matters at sites with long alleles.
        size_t found = numeric_limits<size_t>::max();
        auto positions = allele_positions.find(
            ((int64_t)read_step.node_id << 1) | (int64_t)read_step.backward);
        if (positions != allele_positions.end()) {
            auto at = std::lower_bound(positions->second.begin(), positions->second.end(),
                                       (uint32_t)allele_index);
            if (at != positions->second.end()) {
                found = *at;
            }
        }

        if (found != numeric_limits<size_t>::max()) {
            if (have_anchor) {
                // Allele nodes between the previous pair and this one are sequence the
                // allele has and the read skipped: a deletion. Allele sequence before the
                // first pair or after the last is outside the read's window and is not
                // charged, so that a read is not penalised for being short.
                size_t deleted = 0;
                for (size_t j = allele_index; j < found; ++j) {
                    deleted += allele_steps[j].sequence.size();
                }
                if (deleted > 0) {
                    score += read_scorer.score_gap(deleted);
                }
            }

            score += score_shared_node(aln, read_step, read_scorer, nat_adjust);
            bases_accounted += read_step.read_length;
            have_anchor = true;
            in_insertion = false;
            allele_index = found + 1;
            continue;
        }

        // The allele has no matching node from here on. Either the read took a
        // node this allele lacks, or the two substituted nodes for one another.
        if (allele_index < allele_steps.size()) {
            const AlleleStep& allele_step = allele_steps[allele_index];

            // If the read visits this allele node later on, the current read visit is one
            // the allele lacks: an insertion. We charge it and leave allele_index where it
            // is, so the allele node is still there to pair with the read's later visit.
            auto later = read_last_visit.find(((int64_t)allele_step.node_id << 1)
                                              | (int64_t)allele_step.backward);
            if (later != read_last_visit.end() && later->second > read_index) {
                leave_unpaired(read_step);
                continue;
            }

            // Optimal pairing never substitutes a visit that the other sequence makes elsewhere:
            // here, a read visit the allele made before the last pair, or an allele visit the
            // read made before this one. The read visit is left unpaired instead, and the allele
            // visit stays where it is, to be charged as a deletion at the next pair or left
            // uncharged after the last.
            if (positions != allele_positions.end() || later != read_last_visit.end()) {
                leave_unpaired(read_step);
                continue;
            }

            size_t shared = min(read_step.read_length, allele_step.sequence.size());

            // Score the overlapping extent base by base, so an equal-length
            // substituted node (the common SNP case) is charged as mismatches
            // rather than as a pair of gaps.
            score += score_substitution(aln, read_step.read_offset, shared, allele_step.sequence, 0,
                                        read_scorer);

            // Whatever length the two disagree by is an indel.
            size_t difference = read_step.read_length > allele_step.sequence.size()
                                    ? read_step.read_length - allele_step.sequence.size()
                                    : allele_step.sequence.size() - read_step.read_length;
            if (difference > 0) {
                score += read_scorer.score_gap(difference);
                if (read_step.read_length > allele_step.sequence.size()) {
                    nat_adjust += params.insertion_gap_nats;
                }
            }

            bases_accounted += read_step.read_length;
            in_insertion = false;
            ++allele_index;
        } else {
            // The allele is used up but the read continues. Those read bases cannot be
            // placed on this allele, so they are charged as an insertion, which keeps
            // every allele scored over the same read bases.
            leave_unpaired(read_step);
        }
    }

    // Every read base inside the site must have been scored, whatever the allele, so
    // that all alleles are scored over the same read bases.
    size_t window_bases = 0;
    for (const ReadStep& read_step : read_steps) {
        window_bases += read_step.read_length;
    }
    assert(bases_accounted == window_bases);

    return score;
}

// Optimal pairing (--optimal-pairing): finds the highest-scoring pairing of the read's node visits with
// the allele's by dynamic programming. Base-level edits are still read off the mapper's alignment;
// only the pairing of visits is searched.
//
// The dynamic program fills an m x n grid, m read visits against n allele visits, in four states:
//
//   P  no shared visit paired yet, the leading flank. Allele visits cost nothing here, since
//      allele sequence before the read's window is not the read's to explain; read visits are
//      still charged. Pairing a shared visit leaves P for M; a substitution does not.
//   M  this read visit paired with this allele visit.
//   I  inside a run of read visits the allele lacks.
//   D  inside a run of allele visits the read lacks.
//
// Pairing the same visit scores the read's own edits in it. Pairing two visits that each appear
// nowhere in the other sequence is a substitution. Any other pair is forbidden: a visit the read
// and the allele share may not be paired with a different visit, though it may be left unpaired.
//
// A run of k inserted read visits after the first match is one gap, as in greedy pairing. A read
// visit with no bases, where the read deletes its whole node, is no gap when left unpaired: its
// row passes every state through unchanged. The answer is the best final-row cell in P, M or I. D is excluded because it has charged
// allele visits after the read's window, and P is allowed so that a read that pairs nothing still
// gets a finite score.
int32_t GraphAlignedAlleleLikelihoodCalculator::score_by_optimal_pairing(
    const Alignment& aln, const vector<ReadStep>& read_steps,
    const vector<AlleleStep>& allele_steps, const ReadScratch& scratch,
    const vector<int64_t>& allele_keys, const EditAlignmentScorer& read_scorer,
    bool& placed_out, double& nat_adjust) const {

    placed_out = !allele_steps.empty();
    if (!placed_out) {
        return 0;
    }

    const size_t m = read_steps.size(), n = allele_steps.size();

    const int32_t NEG = numeric_limits<int32_t>::min() / 4;
    // Extending an open gap by one more base, as the difference between a two-base and a
    // one-base gap. Equals gap_extension without assuming the scorer exposes it.
    const int32_t extend_per_base = read_scorer.score_gap(2) - read_scorer.score_gap(1);

    struct Cell { int32_t score; double nats; };

    // A node the read and the allele share may not be paired with a different node, since
    // two walks through one node read the same bases and the mapper has already aligned the
    // read there. Otherwise an insertion allele could explain reads that lack the insertion.
    // A shared visit may still be left unpaired, which can explain a read with poor edits in
    // that node better.
    //
    // "Shared" is tested by membership, which is exact while neither sequence repeats a
    // node. With a repeat it only forbids some substitutions it need not, so the result can
    // be less than optimal but is still a valid pairing.
    const vector<int64_t>& read_keys = scratch.keys;
    auto holds = [](const vector<int64_t>& v, int64_t k) {
        return std::binary_search(v.begin(), v.end(), k);
    };

    // The states P, M, I and D are described above the function. Each cell carries the
    // --insertion-nats total of its own best path in `nats`, beside `score`.
    //
    // I and D may follow each other directly. Base-level alignment usually forbids that, as
    // a duplicate of a substitution, but here a substitution scores two different nodes base
    // by base, which gives a different score from a deletion next to an insertion.
    const Cell NONE{NEG, 0.0};
    auto better = [](const Cell& a, const Cell& b) { return a.score >= b.score ? a : b; };
    auto plus = [NEG](const Cell& c, int32_t d, double dn) {
        return c.score == NEG ? Cell{NEG, 0.0} : Cell{c.score + d, c.nats + dn};
    };

    thread_local vector<Cell> pP, pM, pI, pD, P, M, I, D;
    pP.assign(n + 1, NONE); pM.assign(n + 1, NONE);
    pI.assign(n + 1, NONE); pD.assign(n + 1, NONE);
    P.assign(n + 1, NONE); M.assign(n + 1, NONE);
    I.assign(n + 1, NONE); D.assign(n + 1, NONE);
    for (size_t j = 0; j <= n; ++j) {
        pP[j] = Cell{0, 0.0};       // any allele prefix may be consumed free
    }

    // On large sites, only fill a band of the grid. A read and an allele usually differ at a
    // few nodes, so the best pairing stays near the path through their shared visits. The match
    // cell of a shared visit at read step i_s and allele step j_s is (i_s + 1, j_s + 1). Row i's
    // band spans the projections, along the diagonal, of the last shared visit's match cell in
    // row i or earlier and of the first one after row i, widened by `band` on each side. Between
    // two shared visits the band therefore holds every cell of a deletion or an insertion that
    // joins them, however long; where they share a diagonal it is centred on it. Small sites are filled
    // in full. The band can miss the best pairing, so on large sites this search is an
    // approximation.
    const bool banded = m * n > 20000;
    const size_t band = 64;
    vector<size_t> band_from, band_to;
    if (banded) {
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
        // shared[s - 1] is the last match cell in row i or earlier, shared[s] the first after it.
        size_t s = 0;
        for (size_t i = 0; i <= m; ++i) {
            while (s < shared.size() && shared[s].first <= i) {
                ++s;
            }
            size_t from = i, to = i;   // the i == j diagonal, with no shared visit at all
            if (s > 0) {
                from = to = shared[s - 1].second + (i - shared[s - 1].first);
            }
            if (s < shared.size()) {
                const size_t ahead =
                    shared[s].second - min(shared[s].second, shared[s].first - i);
                from = s > 0 ? min(from, ahead) : ahead;
                to = s > 0 ? max(to, ahead) : ahead;
            }
            from = min(from, n);
            to = min(to, n);
            band_from[i] = from > band ? from - band : 1;
            band_to[i] = min(n, to + band);
        }
    }

    for (size_t i = 1; i <= m; ++i) {
        const ReadStep& rs = read_steps[i - 1];
        const int32_t rlen = (int32_t)rs.read_length;
        const int32_t rgap = read_scorer.score_gap(rs.read_length);
        const double rnat = params.insertion_gap_nats;
        const int64_t rkey = step_key(rs.node_id, rs.backward);
        const bool no_bases = rs.read_length == 0;

        P[0] = no_bases ? pP[0] : plus(pP[0], rgap, rnat);
        M[0] = NONE;
        I[0] = no_bases ? pI[0]
                        : better(plus(better(pM[0], pD[0]), rgap, rnat),
                                 plus(pI[0], rlen * extend_per_base, rnat));
        D[0] = NONE;

        size_t j_lo = 1, j_hi = n;
        if (banded) {
            j_lo = band_from[i];
            j_hi = band_to[i];
            // Cells outside the band this row must not carry a stale value from two rows ago.
            for (size_t j = 1; j < j_lo; ++j) { P[j] = M[j] = I[j] = D[j] = NONE; }
            for (size_t j = j_hi + 1; j <= n; ++j) { P[j] = M[j] = I[j] = D[j] = NONE; }
        }
        for (size_t j = j_lo; j <= j_hi; ++j) {
            const AlleleStep& as = allele_steps[j - 1];
            const int32_t alen = (int32_t)as.sequence.size();
            const int32_t agap = read_scorer.score_gap(as.sequence.size());
            const bool is_match = (as.node_id == rs.node_id && as.backward == rs.backward);

            int32_t pair = NEG;
            double pair_nats = 0.0;
            if (is_match) {
                pair = scratch.own[i - 1];
                pair_nats = scratch.own_nats[i - 1];
            } else if (!holds(allele_keys, rkey) &&
                       !holds(read_keys, step_key(as.node_id, as.backward))) {
                // Neither node anchors elsewhere, so this really is a substitution. The
                // overlapping extent is scored base by base, which keeps an equal-length pair
                // -- the SNP case -- a mismatch rather than two gaps.
                const size_t shared = min(rs.read_length, as.sequence.size());
                pair = score_substitution(aln, rs.read_offset, shared, as.sequence, 0, read_scorer);
                const size_t difference = rs.read_length > as.sequence.size()
                                              ? rs.read_length - as.sequence.size()
                                              : as.sequence.size() - rs.read_length;
                if (difference > 0) {
                    pair += read_scorer.score_gap(difference);
                    if (rs.read_length > as.sequence.size()) {
                        pair_nats = params.insertion_gap_nats;
                    }
                }
            }

            if (pair == NEG) {
                M[j] = NONE;
                P[j] = NONE;
            } else if (is_match) {
                // A match closes the flank, so it may be entered from P as well.
                M[j] = plus(better(better(pM[j - 1], pP[j - 1]), better(pI[j - 1], pD[j - 1])),
                            pair, pair_nats);
                P[j] = NONE;
            } else {
                M[j] = plus(better(pM[j - 1], better(pI[j - 1], pD[j - 1])), pair, pair_nats);
                P[j] = plus(pP[j - 1], pair, pair_nats);
            }
            // Free leading allele deletion.
            P[j] = better(P[j], P[j - 1]);

            if (no_bases) {
                // Leaving a visit with no bases unpaired keeps the state it follows.
                M[j] = better(M[j], pM[j]);
                P[j] = better(P[j], pP[j]);
                I[j] = pI[j];
            } else {
                // A read insertion before any anchor.
                P[j] = better(P[j], plus(pP[j], rgap, rnat));
                // A shared visit may be left unpaired, though it may not be substituted. A read
                // with poor edits inside a shared node can be explained better by a gap.
                I[j] = better(plus(better(pM[j], pD[j]), rgap, rnat),
                              plus(pI[j], rlen * extend_per_base, rnat));
            }
            D[j] = better(plus(better(M[j - 1], I[j - 1]), agap, 0.0),
                          plus(D[j - 1], alen * extend_per_base, 0.0));
            if (no_bases) {
                D[j] = better(D[j], pD[j]);
            }
        }
        pP.swap(P); pM.swap(M); pI.swap(I); pD.swap(D);
    }

    // Every read step is consumed by exactly one transition, so every read base inside the
    // site is scored whatever the allele. Allele nodes after the last pair are outside the
    // window and cost nothing, so the pairing may not end in D, which has charged them.
    Cell best = NONE;
    for (size_t j = 0; j <= n; ++j) {
        best = better(best, better(pP[j], better(pM[j], pI[j])));
    }
    assert(best.score != NEG);
    nat_adjust += best.nats;
    return best.score;
}

void GraphAlignedAlleleLikelihoodCalculator::set_rate_reference(
    const PathPositionHandleGraph* position_graph, const vector<path_handle_t>& reference_paths) {
    rate_graph = position_graph;
    rate_paths = reference_paths;
    lock_guard<std::mutex> guard(window_bp_mutex);
    ref_bucket_counts.clear();
    ref_window_rate.clear();
}

bool GraphAlignedAlleleLikelihoodCalculator::rate_position(const Snarl& snarl, size_t& path_index,
                                                          int64_t& position) const {
    if (rate_graph == nullptr || rate_paths.empty()) {
        return false;
    }
    // Where a boundary node lies on the reference: the earliest step on the first listed
    // reference path that visits it. Path handles, not node IDs, choose among the steps, so
    // the choice survives renumbering.
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
                int64_t pos = (int64_t)rate_graph->get_position_of_step(step);
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
    // A nested site's boundaries are often off the reference. Its ancestors' are not, and an
    // ancestor lies within a window's width of the site unless it is very large.
    const Snarl* current = &snarl;
    while (current != nullptr) {
        if (place(current->start().node_id()) || place(current->end().node_id())) {
            return true;
        }
        const Snarl* managed = snarl_manager.into_which_snarl(current->start().node_id(),
                                                              current->start().backward());
        current = managed == nullptr ? nullptr : snarl_manager.parent_of(managed);
    }
    return false;
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
    // The nodes as sorted, coalesced ID ranges: the read source visits each read once however
    // many ranges it touches.
    vector<nid_t> sorted(nodes);
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    vector<pair<nid_t, nid_t>> ranges;
    for (nid_t id : sorted) {
        if (!ranges.empty() && ranges.back().second + 1 == id) {
            ranges.back().second = id;
        } else {
            ranges.emplace_back(id, id);
        }
    }
    if (ranges.empty()) {
        return counts;
    }

    // Reads are counted as they are at a site: under `depth_effective_reads` each counts as
    // 1 - e_r, with the same clamps and the same `use_mismap_term` switch. Counting them
    // differently would put a constant factor between N and lambda.
    //
    // 1 - e_r depends on the read's MAPQ alone, so reads are tallied by MAPQ and the tallies
    // summed in MAPQ order. The sum then does not depend on the order the read source returns
    // reads in, which can change with its fetch window.
    double& length_total = counts.length_total;
    size_t& length_count = counts.length_count;
    map<int32_t, size_t> by_mapq;
    read_source.for_each_alignment(ranges, [&](const Alignment& aln) {
        // Count only the reads that begin on one of the nodes; the fetch also returns reads
        // that only pass through them.
        const Path& path = aln.path();
        if (path.mapping_size() == 0) {
            return;
        }
        nid_t start_node = path.mapping(0).position().node_id();
        if (!std::binary_search(sorted.begin(), sorted.end(), start_node)) {
            return;
        }
        length_total += (double)aln.sequence().size();
        ++length_count;
        ++by_mapq[aln.mapping_quality()];
    });
    for (const auto& tally : by_mapq) {
        if (!params.depth_effective_reads) {
            counts.reads += (double)tally.second;
            continue;
        }
        double mismap = params.use_mismap_term
                            ? phred_to_prob((double)tally.first)
                            : params.min_mismap_prob;
        counts.reads += (double)tally.second
                        * (1.0 - min(max(mismap, params.min_mismap_prob), params.max_mismap_prob));
    }
    return counts;
}

GraphAlignedAlleleLikelihoodCalculator::StartCounts
GraphAlignedAlleleLikelihoodCalculator::bucket_counts(size_t path_index, int64_t bucket) const {
    path_handle_t path = rate_paths[path_index];
    int64_t lo = bucket * RATE_BUCKET;
    int64_t hi = min<int64_t>((int64_t)rate_graph->get_path_length(path), lo + RATE_BUCKET);
    if (bucket < 0 || lo >= hi) {
        return StartCounts();
    }
    pair<size_t, int64_t> key(path_index, bucket);
    {
        lock_guard<std::mutex> guard(window_bp_mutex);
        auto found = ref_bucket_counts.find(key);
        if (found != ref_bucket_counts.end()) {
            return found->second;
        }
    }

    // The reference nodes whose step begins in [lo, hi), each counted once however often the
    // path visits it.
    vector<nid_t> nodes;
    size_t bp = 0;
    unordered_set<nid_t> seen;
    step_handle_t step = rate_graph->get_step_at_position(path, (size_t)lo);
    while (step != rate_graph->path_end(path)) {
        int64_t start = (int64_t)rate_graph->get_position_of_step(step);
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
    lock_guard<std::mutex> guard(window_bp_mutex);
    ref_bucket_counts[key] = counts;
    return counts;
}

GraphAlignedAlleleLikelihoodCalculator::WindowReadStats
GraphAlignedAlleleLikelihoodCalculator::local_read_stats(
    const Snarl& snarl, const vector<pair<nid_t, nid_t>>& site_ranges) const {

    // The window is fixed rather than taken from the read source's fetch window, so that DR
    // and the depth term do not depend on how the reads were supplied.
    if (site_ranges.empty() || params.depth_ploidy <= 0) {
        return WindowReadStats();
    }
    size_t path_index = 0;
    int64_t position = 0;
    if (!rate_position(snarl, path_index, position)) {
        ++id_fallbacks;
        return id_window_read_stats(site_ranges);
    }
    int64_t bucket = position / RATE_BUCKET;
    pair<size_t, int64_t> key(path_index, bucket);
    {
        lock_guard<std::mutex> guard(window_bp_mutex);
        auto found = ref_window_rate.find(key);
        if (found != ref_window_rate.end()) {
            return found->second;
        }
    }

    // Computed once per bucket and shared by the sites in it, since counting a window's reads
    // costs far more than scoring a site. The window is the site's bucket and one on each
    // side, so that it reaches at least a bucket past the site in both directions.
    StartCounts counts;
    for (int64_t b = bucket - 1; b <= bucket + 1; ++b) {
        counts.add(bucket_counts(path_index, b));
    }
    WindowReadStats stats = counts.stats();
    lock_guard<std::mutex> guard(window_bp_mutex);
    ref_window_rate[key] = stats;
    return stats;
}

GraphAlignedAlleleLikelihoodCalculator::WindowReadStats
GraphAlignedAlleleLikelihoodCalculator::id_window_read_stats(
    const vector<pair<nid_t, nid_t>>& site_ranges) const {

    // The rate window is the block of RATE_ID_WINDOW IDs containing the site's lowest node ID.
    nid_t lo = site_ranges.front().first;
    for (const auto& r : site_ranges) {
        lo = min(lo, r.first);
    }
    size_t window_index = (size_t)(lo / RATE_ID_WINDOW);
    {
        lock_guard<std::mutex> guard(window_bp_mutex);
        auto found = window_rate.find(window_index);
        if (found != window_rate.end()) {
            return found->second;
        }
    }

    nid_t first = (nid_t)window_index * RATE_ID_WINDOW;
    nid_t last = first + RATE_ID_WINDOW - 1;

    // Node IDs are dense in a GBZ but not guaranteed to be, so ask the graph.
    vector<nid_t> nodes;
    size_t bp = 0;
    for (nid_t id = first; id <= last; ++id) {
        if (graph.has_node(id)) {
            nodes.push_back(id);
            bp += graph.get_length(graph.get_handle(id));
        }
    }

    WindowReadStats stats = count_starts(nodes, bp).stats();
    lock_guard<std::mutex> guard(window_bp_mutex);
    window_rate[window_index] = stats;
    return stats;
}

AlleleReadLikelihoods GraphAlignedAlleleLikelihoodCalculator::compute(
    const Snarl& snarl, const vector<SnarlTraversal>& traversals, int region_ploidy) {


    AlleleReadLikelihoodsBuilder builder(traversals.size(), params.min_mismap_prob,
                                        params.max_mismap_prob);
    if (traversals.empty()) {
        return builder.build();
    }

    // Nodes making up the site, including its boundaries.
    auto contents = snarl_manager.deep_contents(&snarl, graph, true);
    unordered_set<nid_t> site_nodes(contents.first.begin(), contents.first.end());
    if (site_nodes.empty()) {
        return builder.build();
    }
    unordered_set<nid_t> boundary_nodes{snarl.start().node_id(), snarl.end().node_id()};

    // Merge the site's node IDs into ranges for the read source to query.
    vector<nid_t> sorted_ids(site_nodes.begin(), site_nodes.end());
    sort(sorted_ids.begin(), sorted_ids.end());
    vector<pair<nid_t, nid_t>> ranges;
    for (nid_t id : sorted_ids) {
        if (!ranges.empty() && id == ranges.back().second + 1) {
            ranges.back().second = id;
        } else {
            ranges.emplace_back(id, id);
        }
    }

    // Materialise each allele's node sequences once per site. This is per allele,
    // not per (read, allele), so it stays off the hot path.
    vector<vector<AlleleStep>> allele_steps;
    allele_steps.reserve(traversals.size());
    for (const SnarlTraversal& traversal : traversals) {
        allele_steps.push_back(get_allele_steps(traversal));
    }
    // Sorted once per allele, not once per (read, allele): the keys do not mention the read.
    vector<vector<int64_t>> allele_keys;
    allele_keys.reserve(allele_steps.size());
    for (const auto& steps : allele_steps) {
        allele_keys.push_back(sorted_allele_keys(steps));
    }
    // Likewise once per allele, and only for greedy pairing, which is its only consumer.
    vector<AlleleStepPositions> allele_positions;
    if (!params.optimal_pairing) {
        allele_positions.reserve(allele_steps.size());
        for (const auto& steps : allele_steps) {
            allele_positions.push_back(index_allele_steps(steps));
        }
    }

    // Each allele's length for the depth term's lambda, T_a: the traversal without the
    // site's two boundary nodes. A read that lies inside one boundary node is dropped by
    // get_read_steps, so boundary sequence yields no rows. An allele with no interior, which
    // deletes the whole interior, gets length 0, and lambda then counts only the R - 1
    // reads that span its junction.
    //
    // This differs from the unique lengths the mixture weights use: lambda asks how much
    // sequence yields reads, and the weights ask which sequence tells the alleles apart.
    vector<size_t> depth_lengths;
    depth_lengths.reserve(allele_steps.size());
    for (const auto& steps : allele_steps) {
        size_t len = 0;
        for (const AlleleStep& step : steps) {
            if (!boundary_nodes.count(step.node_id)) {
                len += step.sequence.size();
            }
        }
        depth_lengths.push_back(len);
    }

    if (params.length_weighted_mixture) {
        // Length of each allele as it is actually spelled by its traversal, which
        // is the quantity the mixture weight needs. Computed from the same steps
        // the scorer uses, so it cannot drift from what the reads are scored
        // against.
        vector<size_t> allele_lengths;
        allele_lengths.reserve(allele_steps.size());
        for (const auto& steps : allele_steps) {
            size_t len = 0;
            for (const AlleleStep& step : steps) {
                len += step.sequence.size();
            }
            allele_lengths.push_back(len);
        }
        builder.set_allele_lengths(allele_lengths);

        {
            // Per-allele node content, then pairwise set differences. Computed once
            // per site off the hot path; a node visited more than once by an allele
            // counts its sequence once, which is what "does this allele carry this
            // sequence" means.
            vector<map<nid_t, size_t>> content(allele_steps.size());
            for (size_t a = 0; a < allele_steps.size(); ++a) {
                for (const AlleleStep& step : allele_steps[a]) {
                    content[a][step.node_id] = step.sequence.size();
                }
            }
            vector<vector<size_t>> unique_lengths(
                allele_steps.size(), vector<size_t>(allele_steps.size(), 0));
            for (size_t a = 0; a < content.size(); ++a) {
                for (size_t b = 0; b < content.size(); ++b) {
                    if (a == b) {
                        continue;
                    }
                    size_t total = 0;
                    for (const auto& entry : content[a]) {
                        if (!content[b].count(entry.first)) {
                            total += entry.second;
                        }
                    }
                    unique_lengths[a][b] = total;
                }
            }
            builder.set_unique_lengths(std::move(unique_lengths));
        }
    }

    // The orientation each allele visits each node in, so a read aligned to the
    // opposite strand can be flipped into the alleles' reading direction before
    // being compared to them. A node different alleles disagree about is left out
    // rather than allowed to cast a vote.
    unordered_map<nid_t, bool> allele_orientations;
    unordered_set<nid_t> ambiguous_orientation;
    for (const auto& steps : allele_steps) {
        for (const AlleleStep& step : steps) {
            auto found = allele_orientations.find(step.node_id);
            if (found == allele_orientations.end()) {
                allele_orientations[step.node_id] = step.backward;
            } else if (found->second != step.backward) {
                ambiguous_orientation.insert(step.node_id);
            }
        }
    }
    for (nid_t node_id : ambiguous_orientation) {
        allele_orientations.erase(node_id);
    }

    // Node lengths for reverse-complementing a read's alignment.
    function<int64_t(nid_t)> node_length = [&](nid_t node_id) {
        return (int64_t)graph.get_length(graph.get_handle(node_id));
    };

    // Anchor evidence, when anchors are armed. Filled per read alongside the matrix rows, and in
    // the same order: `add_read` drops a read that placed on nothing, so the evidence follows what
    // it kept rather than what was offered.
    unique_ptr<AnchorSiteEvidence> anchor_evidence;
    // Read phasing needs only the `rel` rows, without read names or positions. Built only when
    // anchors are not being written, since the anchor evidence already contains the rows.
    unique_ptr<PhaseReadEvidence> phase_evidence;
    if (params.collect_read_phasing && !params.collect_anchors) {
        phase_evidence = make_unique<PhaseReadEvidence>();
        phase_evidence->n_alleles = traversals.size();
        phase_evidence->length_weighted = params.length_weighted_mixture;
    }
    if (params.collect_anchors) {
        anchor_evidence = make_unique<AnchorSiteEvidence>();
        anchor_evidence->n_alleles = traversals.size();
        anchor_evidence->start_node = snarl.start().node_id();
        anchor_evidence->end_node = snarl.end().node_id();
        anchor_evidence->length_weighted = params.length_weighted_mixture;
        // The alleles as spelled, for the mixture weights the per-read score uses. Computed here
        // whatever the mixture setting, because the score's weighting is its own decision.
        anchor_evidence->allele_length.reserve(allele_steps.size());
        for (const auto& steps : allele_steps) {
            size_t len = 0;
            for (const AlleleStep& step : steps) {
                len += step.sequence.size();
            }
            anchor_evidence->allele_length.push_back((uint32_t)len);
        }
    }
    if (phase_evidence != nullptr) {
        phase_evidence->allele_length.reserve(allele_steps.size());
        for (const auto& steps : allele_steps) {
            size_t len = 0;
            for (const AlleleStep& step : steps) {
                len += step.sequence.size();
            }
            phase_evidence->allele_length.push_back((uint32_t)len);
        }
    }

    // Counted into the run's counters where there are some, and into this local otherwise, as
    // in a unit test that drives the calculator directly. Declared outside the per-read
    // callback so that it is constructed once per site.
    AnchorCounters unowned_counters;
    AnchorCounters& anchor_pin_counters =
        params.anchor_counters != nullptr ? *params.anchor_counters : unowned_counters;

    vector<ReadStep> read_steps;
    ReadScratch read_scratch;
    vector<double> row(traversals.size());

    read_source.for_each_read(ranges, [&](const SiteRead& read) {
        const Alignment& aln = *read.aln;

        if (!get_read_steps(read, site_nodes, boundary_nodes, read_steps)) {
            return;
        }

        // Flip a reverse-strand read into the alleles' reading direction, so that its
        // visits can be paired with theirs. Only the stretch inside the site is flipped; see
        // site_span_of.
        Alignment flipped;
        const Alignment* scored_aln = &aln;
        if (read_is_reverse_of_alleles(read_steps, allele_orientations)) {
            flipped = reverse_complement_alignment(site_span_of(read, site_nodes), node_length);
            scored_aln = &flipped;
            // The slice is a fresh alignment holding only the site's mappings, so there is
            // no index for it and none is wanted: walking it is already walking the site.
            SiteRead flipped_read;
            flipped_read.aln = &flipped;
            if (!get_read_steps(flipped_read, site_nodes, boundary_nodes, read_steps)) {
                // Should not happen: flipping preserves which nodes are touched.
                return;
            }
        }

        // Pick the scorer this read can actually support. A read with no base
        // qualities cannot be scored by a quality-adjusted scorer without
        // inventing data.
        const EditAlignmentScorer& read_scorer =
            scored_aln->quality().empty() ? plain_scorer : qual_scorer;
        double log_base = read_scorer.get_log_base();

        // Depends on the read alone, so it is computed here rather than per allele. Only
        // optimal pairing consumes it, and building it is not free, so skip it for greedy pairing.
        if (params.optimal_pairing) {
            prepare_read_scratch(*scored_aln, read_steps, read_scorer, read_scratch);
        }

        for (size_t a = 0; a < traversals.size(); ++a) {
            bool placed = false;
            double nat_adjust = 0.0;
            int32_t score =
                params.optimal_pairing
                    ? score_by_optimal_pairing(*scored_aln, read_steps, allele_steps[a],
                                               read_scratch, allele_keys[a], read_scorer,
                                               placed, nat_adjust)
                    : score_by_greedy_pairing(*scored_aln, read_steps, allele_steps[a],
                                              allele_positions[a], read_scorer,
                                              placed, nat_adjust);
            // nat_adjust carries corrections the int32 score cannot express; see
            // AlleleLikelihoodParams::insertion_gap_nats. It is zero by default.
            row[a] = placed ? log_base * (double)score + nat_adjust
                            : -numeric_limits<double>::infinity();
        }

        // MAPQ is a phred probability that the read is in the wrong place. vg's
        // mappers build the score vector it comes from out of distinct graph
        // placements, so it measures "this read could be at another locus"
        // rather than "another allele of this snarl fits too" -- which is what
        // this term needs it to mean.
        double mismap = params.use_mismap_term
                            ? phred_to_prob((double)aln.mapping_quality())
                            : params.min_mismap_prob;

        const Position* start =
            aln.path().mapping_size() > 0 ? &aln.path().mapping(0).position() : nullptr;
        if (!builder.add_read(row, mismap, aln.name(), aln.sequence().size(), start)) {
            return;
        }

        if (phase_evidence != nullptr) {
            // Hashed here, where the name is in hand, and never stored: the string is the whole
            // reason the anchor evidence is too heavy to retain for phasing.
            phase_evidence->read_key.push_back((uint64_t)std::hash<string>{}(aln.name()));
        }
        if (anchor_evidence != nullptr) {
            // Resolved against `aln`, not `scored_aln`. The flipped copy's offsets are in the
            // reverse-complemented read, which the read file does not contain; the record's
            // strand field gives the orientation instead.
            AnchorRead record;
            record.name = aln.name();
            record.mismap = (float)mismap;
            record.start_pin = resolve_anchor_pin(read, graph, snarl.start().node_id(),
                                                  snarl.start().backward(), true,
                                                  anchor_pin_counters);
            record.end_pin = resolve_anchor_pin(read, graph, snarl.end().node_id(),
                                                snarl.end().backward(), false,
                                                anchor_pin_counters);
            anchor_evidence->reads.push_back(std::move(record));
        }
    });

    AlleleReadLikelihoods result = builder.build();

    // The builder put the rows in its canonical order; put what was kept per read alongside
    // them in the same order, so that row r is still reads[r] and read_key[r].
    const vector<size_t>& order = builder.row_order();
    if (anchor_evidence != nullptr && anchor_evidence->reads.size() == order.size()) {
        vector<AnchorRead> reordered;
        reordered.reserve(order.size());
        for (size_t i : order) {
            reordered.push_back(std::move(anchor_evidence->reads[i]));
        }
        anchor_evidence->reads = std::move(reordered);
    }
    if (phase_evidence != nullptr && phase_evidence->read_key.size() == order.size()) {
        vector<uint64_t> reordered;
        reordered.reserve(order.size());
        for (size_t i : order) {
            reordered.push_back(phase_evidence->read_key[i]);
        }
        phase_evidence->read_key = std::move(reordered);
    }

    // R, the mean read length, for all of its users: the depth term's lambda, the mixture
    // weights and the anchor slot weights. The builder's R is the mean over this site's reads,
    // which over-represents long reads because a long read reaches more sites. The rate
    // window's mean counts each read once, where it begins, so we use it instead.
    WindowReadStats stats = local_read_stats(snarl, ranges);
    if (stats.mean_read_length > 0.0) {
        result.set_mean_read_length(stats.mean_read_length);
    }

    if (anchor_evidence != nullptr) {
        // The rows in the builder's order, which `reads` was put in above, so row r is reads[r].
        anchor_evidence->mean_read_length = (float)result.mean_read_length_estimate();
        anchor_evidence->rel.resize(result.num_reads() * result.num_alleles());
        for (size_t r = 0; r < result.num_reads(); ++r) {
            anchor_evidence->reads[r].mismap = (float)result.mismap_prob(r);
            for (size_t a = 0; a < result.num_alleles(); ++a) {
                anchor_evidence->rel[r * result.num_alleles() + a] = (float)result.rel(r, a);
            }
        }
        result.anchor_evidence = std::move(anchor_evidence);
    }
    if (phase_evidence != nullptr) {
        // Same rows in the same order, so row r is read_key[r]: the keys are pushed inside the
        // callback, after `add_read` has accepted the read, so a read it drops contributes neither,
        // and were put in the builder's order above.
        phase_evidence->mean_read_length = (float)result.mean_read_length_estimate();
        // If the counts disagree, drop the site's phasing evidence. Padding the keys instead
        // would give every padded read the same key, so they would be taken for one fragment
        // at every site.
        if (phase_evidence->read_key.size() != result.num_reads()) {
            phase_evidence.reset();
        }
    }
    if (phase_evidence != nullptr) {
        phase_evidence->mismap.resize(result.num_reads());
        phase_evidence->rel.resize(result.num_reads() * result.num_alleles());
        for (size_t r = 0; r < result.num_reads(); ++r) {
            phase_evidence->mismap[r] = (float)result.mismap_prob(r);
            for (size_t a = 0; a < result.num_alleles(); ++a) {
                phase_evidence->rel[r * result.num_alleles() + a] = (float)result.rel(r, a);
            }
        }
        result.phase_evidence = std::move(phase_evidence);
    }
    if (!depth_lengths.empty()) {
        // Set whether or not the depth term is on, so that DR is always written; a zero
        // weight leaves the likelihood unchanged. The rate is per haplotype, and the window
        // counts reads from every haplotype in the region, not only those crossing the site.
        int haplotypes = region_ploidy > 0 ? region_ploidy : params.depth_ploidy;
        result.set_depth_context(depth_lengths,
                                 stats.start_rate / (double)haplotypes,
                                 // The window's population mean, not the site's size-biased one.
                                 // Falls back to the site mean where no read began in the window,
                                 // which is a small window at the end of a contig rather than a
                                 // condition worth a branch elsewhere.
                                 stats.mean_read_length > 0.0 ? stats.mean_read_length
                                                              : result.mean_read_length_estimate(),
                                 params.depth_weight, params.depth_effective_reads);
    }
    return result;
}

}

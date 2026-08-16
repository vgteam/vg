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
    // sum_r (1 - e_r), the expected number of these reads that are genuinely from
    // this locus. Cached because the depth term asks for it once per genotype.
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

const string& AlleleReadLikelihoods::read_name(size_t r) const {
    return read_names.at(r);
}

double AlleleReadLikelihoods::expected_reads(const vector<int>& genotype) const {
    double total = 0.0;
    for (int allele : genotype) {
        double len = (allele >= 0 && (size_t)allele < traversal_lengths.size())
                         ? (double)traversal_lengths[allele] : 0.0;
        total += max(len + depth_read_length - 1.0, 1.0);
    }
    return depth_rate * total;
}

double AlleleReadLikelihoods::observed_reads() const {
    return depth_effective ? effective_read_total : (double)n_reads;
}

double AlleleReadLikelihoods::depth_ratio(const vector<int>& genotype) const {
    double expected = expected_reads(genotype);
    return expected > 0.0 ? observed_reads() / expected : -1.0;
}

/// ln of a Poisson pmf, continued to real `n` through lgamma. The observation is
/// `sum_r (1 - e_r)` rather than a row count, so it is fractional by construction.
/// The `-ln n!` normaliser is the same for every genotype at a site and cancels in
/// every comparison; it is carried anyway because GL is reported, not just ranked.
static double ln_poisson_pmf(double n, double lambda) {
    if (lambda <= 0.0) {
        return -numeric_limits<double>::infinity();
    }
    if (n <= 0.0) {
        return -lambda;
    }
    return n * log(lambda) - lambda - lgamma(n + 1.0);
}

vector<double> AlleleReadLikelihoods::mixture_weights(const vector<int>& genotype) const {
    if (genotype.empty()) {
        return {};
    }
    double flat = 1.0 / (double)genotype.size();

    // Expected share of this site's reads per haplotype of the genotype. Flat
    // 1/|G| unless lengths were supplied; see set_length_weights for why the flat
    // weight is wrong whenever the alleles differ in length.
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
            double eff = own + mean_read_length - 1.0;
            // A traversal shorter than one read still admits reads spanning it,
            // so the effective length can never fall to zero.
            weights[i] = max(eff, 1.0);
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

        // Fold in "this read did not come from this site at all". Because the
        // rows are normalised the background is exactly 1, so the bracket lies
        // in [e_r, 1] and its log is always finite: no logsumexp needed, and no
        // single read can penalise a genotype without bound.
        double e_r = read_mismap_prob[r];
        total += log((1.0 - e_r) * mixture + e_r);
    }

    if (uses_depth_term()) {
        total += depth_weight * ln_poisson_pmf(observed_reads(), expected_reads(genotype));
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

    // The reads' *own* e_r is deliberately not used here, and this is the one thing
    // about this function that is easy to get wrong -- the first version used it.
    //
    // An ideal read is well-fitting *and* well-mapped, so the denominator is built at
    // the mismap floor. Using each read's own e_r instead folds the site's unreliability
    // into both sides of the ratio, where it cancels: a window of MAPQ-0 reads has
    // -ln(0.7) = 0.36 of achievable gap per read against -ln(0.02) = 3.91, so a badly
    // mapped site is scored against a denominator small enough to make a weak call look
    // strong. Measured over the titration, the per-read-e_r version scored 0.427 against
    // 0.347 for raw GQ -- worse than no normalisation at all -- while the floor version
    // scores 0.260.
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
    // non-decreasing tuples in colexicographic order reproduces exactly that
    // ordering, so a genotype's position here is its GL index. For diploid this
    // is (0,0), (0,1), (1,1), (0,2), (1,2), (2,2), ... Note this differs from the
    // order PoissonSupportSnarlCaller emits, which iterates the low allele in the
    // outer loop and so disagrees with the spec once there are 3+ alleles.
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

void AlleleReadLikelihoodsBuilder::add_read(const vector<double>& raw_ln_likelihood,
                                            double mismap_prob, const string& name,
                                            size_t read_length) {
    assert(raw_ln_likelihood.size() == n_alleles);

    // The row's divisor is the read's best fit over ALL alleles at the site, not
    // just those in some genotype. That keeps it genotype-independent, which is
    // what lets it drop out of every genotype comparison.
    double best = -numeric_limits<double>::infinity();
    for (double ll : raw_ln_likelihood) {
        best = max(best, ll);
    }

    if (!(best > -numeric_limits<double>::infinity())) {
        // This read placed on nothing at all, so there is no row maximum to
        // divide by. Normalising would give NaN and quietly poison every
        // genotype at the site. Drop it and count it: a rising count means the
        // read source is over-fetching, or the reads and graph do not match.
        ++unplaceable;
        return;
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
}

AlleleReadLikelihoods AlleleReadLikelihoodsBuilder::build() {
    size_t n_reads = rows.size();
    vector<double> matrix;
    matrix.reserve(n_reads * n_alleles);
    for (auto& row : rows) {
        matrix.insert(matrix.end(), row.begin(), row.end());
    }

    AlleleReadLikelihoods result;
    result.set_contents(n_reads, n_alleles, std::move(matrix), std::move(mismap_probs),
                        std::move(best_lns), std::move(names), unplaceable);
    result.set_mismap_floor(min_mismap);
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

bool GraphAlignedAlleleLikelihoodCalculator::get_read_steps(
    const Alignment& aln, const unordered_set<nid_t>& site_nodes,
    const unordered_set<nid_t>& boundary_nodes, vector<ReadStep>& steps_out) const {

    steps_out.clear();
    const Path& path = aln.path();

    // Track the read offset across every mapping, including those outside the
    // site: otherwise the offsets of the ones inside would be wrong.
    size_t read_offset = 0;
    bool touches_interior = false;

    for (int64_t i = 0; i < path.mapping_size(); ++i) {
        const Mapping& mapping = path.mapping(i);
        size_t to_length = (size_t)mapping_to_length(mapping);
        nid_t node_id = mapping.position().node_id();

        if (site_nodes.count(node_id)) {
            ReadStep step;
            step.node_id = node_id;
            step.backward = mapping.position().is_reverse();
            step.read_offset = read_offset;
            step.read_length = to_length;
            step.mapping = &mapping;
            steps_out.push_back(step);

            if (!boundary_nodes.count(node_id)) {
                touches_interior = true;
            }
        }

        read_offset += to_length;
    }

    if (steps_out.empty()) {
        return false;
    }

    // A read is informative if it can discriminate between alleles at all.
    //
    // Touching an interior node does it, but that is NOT the only way, and
    // assuming it was silently destroyed deletion genotyping: a read that
    // traverses straight from one boundary node to the other uses the deletion
    // edge, touches no interior node, and is the *only* direct evidence the
    // deletion allele ever gets. The discriminating signal there is in the edge,
    // not in the node set.
    //
    // So: informative if it touches an interior node, or if it moves between two
    // distinct nodes inside the site (which for a boundary-to-boundary read means
    // it used an edge no reference traversal has). A read sitting entirely within
    // one boundary node has neither and genuinely cannot discriminate.
    bool uses_internal_edge = false;
    for (size_t i = 1; i < steps_out.size(); ++i) {
        if (steps_out[i].node_id != steps_out[i - 1].node_id) {
            uses_internal_edge = true;
            break;
        }
    }

    if (!touches_interior && !uses_internal_edge) {
        // Genuinely uninformative: contributes an identical constant to every
        // allele. Note this is deliberately not the same as failing to place on
        // some allele, which is informative and must be kept.
        steps_out.clear();
        return false;
    }

    return true;
}

bool GraphAlignedAlleleLikelihoodCalculator::read_is_reverse_of_alleles(
    const vector<ReadStep>& read_steps,
    const unordered_map<nid_t, bool>& allele_orientations) const {

    // Alleles all run from the snarl's start to its end, so they impose a reading
    // direction on the site. A read aligned to the opposite strand visits the same
    // nodes with the opposite orientation flag, and would fail to anchor against
    // any of them -- which silently mis-scored every reverse-strand read, roughly
    // half of all reads, against the wrong allele.
    //
    // Decide by vote rather than from a single step, so a node that different
    // alleles visit in different orientations cannot flip the whole read.
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
    const Alignment& aln, const ReadStep& step, const EditAlignmentScorer& read_scorer) const {

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
            // Insertion relative to the graph.
            score += read_scorer.score_gap(to_len - from_len);
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
    // consecutive differing bases so each run is scored in one call. Charging
    // per base like this is exactly what scoring from the graph's implied
    // alignment buys over a length-averaged DP score.
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

int32_t GraphAlignedAlleleLikelihoodCalculator::score_read_against_allele(
    const Alignment& aln, const vector<ReadStep>& read_steps,
    const vector<AlleleStep>& allele_steps, const EditAlignmentScorer& read_scorer,
    bool& placed_out) const {

    placed_out = !allele_steps.empty();
    if (!placed_out) {
        return 0;
    }

    int32_t score = 0;
    size_t allele_index = 0;
    bool have_anchor = false;
    size_t bases_accounted = 0;

    for (const ReadStep& read_step : read_steps) {
        // Look for this read node ahead in the allele. Anchoring on shared node
        // visits is what makes this a read-off of the alignment the graph already
        // asserts rather than an alignment we invent.
        size_t found = numeric_limits<size_t>::max();
        for (size_t j = allele_index; j < allele_steps.size(); ++j) {
            if (allele_steps[j].node_id == read_step.node_id &&
                allele_steps[j].backward == read_step.backward) {
                found = j;
                break;
            }
        }

        if (found != numeric_limits<size_t>::max()) {
            if (have_anchor) {
                // Allele nodes between the previous anchor and this one are
                // sequence the allele has and the read skipped: a deletion.
                //
                // Only *internal* skips count. Allele sequence before the first
                // anchor or after the last is simply outside the read's window,
                // and charging for it would penalise a read for being short --
                // the very length artefact the window invariant exists to stop.
                size_t deleted = 0;
                for (size_t j = allele_index; j < found; ++j) {
                    deleted += allele_steps[j].sequence.size();
                }
                if (deleted > 0) {
                    score += read_scorer.score_gap(deleted);
                }
            }

            score += score_shared_node(aln, read_step, read_scorer);
            bases_accounted += read_step.read_length;
            have_anchor = true;
            allele_index = found + 1;
            continue;
        }

        // The allele has no matching node from here on. Either the read took a
        // node this allele lacks, or the two substituted nodes for one another.
        if (allele_index < allele_steps.size()) {
            const AlleleStep& allele_step = allele_steps[allele_index];
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
            }

            bases_accounted += read_step.read_length;
            ++allele_index;
        } else {
            // The allele is exhausted but the read continues. Those read bases
            // cannot be placed on this allele, so they are charged as an
            // insertion. They are NOT dropped: omitting them would score this
            // allele over fewer read bases than its competitors and fabricate a
            // likelihood ratio out of the length difference alone.
            score += read_scorer.score_gap(read_step.read_length);
            bases_accounted += read_step.read_length;
        }
    }

    // The window invariant: every read base inside the site was accounted for,
    // whatever the allele. If this ever fires, some allele is being scored over a
    // different span than its competitors and the likelihoods are miscalibrated
    // in a way that still produces plausible-looking VCF.
    size_t window_bases = 0;
    for (const ReadStep& read_step : read_steps) {
        window_bases += read_step.read_length;
    }
    assert(bases_accounted == window_bases);

    return score;
}

double GraphAlignedAlleleLikelihoodCalculator::local_read_rate(
    const vector<pair<nid_t, nid_t>>& site_ranges) const {

    // The neighbourhood the rate is measured over, in node IDs. Prefer the source's own
    // fetch window so the query lands on a cache entry the site already populated.
    //
    // The fallback is not a tuning knob and was demoted from one: it exists so that a
    // source with no window of its own -- an in-memory source answers every range
    // exactly -- still produces the same `DR` as an indexed source over the same reads.
    // Without it the emitted FORMAT fields depended on how the reads were supplied, which
    // an in-tree test caught immediately. Measured through --read-window on real data, the
    // width saturates well below the 4096 an indexed source uses: 1024 loses 0.004
    // structural-variant F1, 16384 gains 0.001. There is nothing here to tune, and every
    // indexed source supplies a window anyway, so this constant is only ever reached by
    // in-memory sources and tests.
    static const size_t FALLBACK_RATE_WINDOW = 4096;
    size_t span = read_source.get_window_span();
    if (span == 0) {
        span = FALLBACK_RATE_WINDOW;
    }
    if (site_ranges.empty() || params.depth_ploidy <= 0) {
        return 0.0;
    }
    // The window the source would have fetched to answer this site's own query.
    nid_t lo = site_ranges.front().first;
    for (const auto& r : site_ranges) {
        lo = min(lo, r.first);
    }
    size_t window_index = (size_t)(lo / (nid_t)span);
    {
        lock_guard<std::mutex> guard(window_bp_mutex);
        auto found = window_rate.find(window_index);
        if (found != window_rate.end()) {
            return found->second;
        }
    }

    nid_t first = (nid_t)window_index * (nid_t)span;
    nid_t last = first + (nid_t)span - 1;

    // Memoised per window, not per site. Counting a window's reads costs O(reads in
    // window), which is far more than a snarl, and neighbouring sites share windows;
    // without this the diagnostic would cost more than the genotyping.
    //
    // Reads are counted here exactly as they are counted at a site: under
    // `depth_effective_reads` each contributes `1 - e_r` rather than 1, using the same
    // MAPQ, the same clamps and the same `use_mismap_term` switch. Weighting one side
    // and not the other would put a constant scale factor between N and lambda and
    // bias every DR in the same direction, which is not a signal.
    double reads = 0.0;
    read_source.for_each_read({{first, last}}, [&](const Alignment& aln) {
        if (!params.depth_effective_reads) {
            reads += 1.0;
            return;
        }
        double mismap = params.use_mismap_term
                            ? phred_to_prob((double)aln.mapping_quality())
                            : params.min_mismap_prob;
        reads += 1.0 - min(max(mismap, params.min_mismap_prob), params.max_mismap_prob);
    });

    // Node IDs are dense in a GBZ but not guaranteed to be, so ask the graph.
    size_t bp = 0;
    for (nid_t id = first; id <= last; ++id) {
        if (graph.has_node(id)) {
            bp += graph.get_length(graph.get_handle(id));
        }
    }

    double rate = (reads <= 0.0 || bp == 0) ? 0.0 : reads / (double)bp;
    lock_guard<std::mutex> guard(window_bp_mutex);
    window_rate[window_index] = rate;
    return rate;
}

AlleleReadLikelihoods GraphAlignedAlleleLikelihoodCalculator::compute(
    const Snarl& snarl, const vector<SnarlTraversal>& traversals, int ploidy) {

    last_site_reads = 0;
    last_site_uninformative = 0;

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

    // Per-allele length for the depth term's lambda: the sequence over which a read
    // can become a row of this matrix, which is the traversal's **interior** only.
    //
    // Not the whole traversal, and the difference is not cosmetic. A SnarlTraversal
    // runs from the snarl's start visit to its end visit inclusive, so its length
    // includes both boundary nodes -- but get_read_steps drops a read that sits
    // entirely inside one boundary node as uninformative, so those bases recruit no
    // reads. Counting them made lambda too large by roughly the two anchors' length,
    // a constant per site, which put the median DR at 0.59 instead of 1 and diluted
    // exactly the contrast between genotypes the term exists to see. It showed as DR
    // rising with event size -- 0.58 at SNVs, 0.87 above 1 kb -- because a fixed
    // overhead matters less the longer the allele.
    //
    // A traversal with no interior at all is the deletion edge: no interior node, so
    // only a junction-spanning read can be a row, and max(len + R - 1, 1) gives the
    // R - 1 junction positions, which is right.
    //
    // Distinct from the unique content the mixture weights use: lambda asks how much
    // sequence generates reads, the weights ask which sequence tells alleles apart.
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

    vector<ReadStep> read_steps;
    vector<double> row(traversals.size());

    read_source.for_each_read(ranges, [&](const Alignment& aln) {
        ++last_site_reads;

        if (!get_read_steps(aln, site_nodes, boundary_nodes, read_steps)) {
            ++last_site_uninformative;
            return;
        }

        // Flip a reverse-strand read into the alleles' reading direction. Without
        // this it anchors on nothing, falls through to the substitution path, and
        // is scored against the wrong allele entirely.
        Alignment flipped;
        const Alignment* scored_aln = &aln;
        if (read_is_reverse_of_alleles(read_steps, allele_orientations)) {
            flipped = reverse_complement_alignment(aln, node_length);
            scored_aln = &flipped;
            if (!get_read_steps(flipped, site_nodes, boundary_nodes, read_steps)) {
                // Should not happen: flipping preserves which nodes are touched.
                ++last_site_uninformative;
                return;
            }
        }

        // Pick the scorer this read can actually support. A read with no base
        // qualities cannot be scored by a quality-adjusted scorer without
        // inventing data.
        const EditAlignmentScorer& read_scorer =
            scored_aln->quality().empty() ? plain_scorer : qual_scorer;
        double log_base = read_scorer.get_log_base();

        for (size_t a = 0; a < traversals.size(); ++a) {
            bool placed = false;
            int32_t score = score_read_against_allele(*scored_aln, read_steps, allele_steps[a],
                                                      read_scorer, placed);
            row[a] = placed ? log_base * (double)score : -numeric_limits<double>::infinity();
        }

        // MAPQ is a phred probability that the read is in the wrong place. vg's
        // mappers build the score vector it comes from out of distinct graph
        // placements, so it measures "this read could be at another locus"
        // rather than "another allele of this snarl fits too" -- which is what
        // this term needs it to mean.
        double mismap = params.use_mismap_term
                            ? phred_to_prob((double)aln.mapping_quality())
                            : params.min_mismap_prob;

        builder.add_read(row, mismap, aln.name(), aln.sequence().size());
    });

    AlleleReadLikelihoods result = builder.build();
    if (!depth_lengths.empty()) {
        // Set unconditionally so `DR` is emitted whether or not the term is armed:
        // the observable should be measurable as a ranking signal before the model
        // is allowed to act on it. A zero weight leaves the likelihood untouched.
        // Per haplotype, from the site's own ploidy rather than an assumed one.
        int effective_ploidy = ploidy > 0 ? ploidy : params.depth_ploidy;
        result.set_depth_context(depth_lengths,
                                 local_read_rate(ranges) / (double)effective_ploidy,
                                 result.mean_read_length_estimate(),
                                 params.depth_weight, params.depth_effective_reads);
    }
    return result;
}

}

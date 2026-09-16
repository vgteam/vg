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

bool AlleleReadLikelihoodsBuilder::add_read(const vector<double>& raw_ln_likelihood,
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
    return true;
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
    if (read_length_count > 0) {
        // Always, not only under the length-weighted mixture: the depth term's
        // lambda = rate * (L + R - 1) needs R whichever mixture is in use, and gating it on
        // allele_lengths left R at 0 under --flat-mixture while the term stayed armed.
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
/// Exists so that flipping a reverse-strand read costs the site rather than the whole read.
/// `reverse_complement_alignment` deep-copies every mapping, every edit and the sequence, and it runs
/// once per (read, site). For a 150 bp read that is nothing; for a 19 kb ONT read delivered to the
/// ~160 sites it spans it was **67.6% of the entire scoring cost** in a profile of chr20.
///
/// Slicing is exact, not an approximation. Flipping the slice puts the step that was at `off` at
/// `to - off - len`, where the full flip puts it at `total - off - len`; the two differ by the
/// constant `total - to`, and the sliced sequence is shifted by exactly the same constant. Base for
/// base, `revcomp(read)[full] == revcomp(read[from:to])[sliced]`. Only read_offset and read_length
/// are ever used to index the sequence, so nothing else can see the difference.
static Alignment site_span_of(const SiteRead& read, const unordered_set<nid_t>& site_nodes) {
    const Alignment& aln = *read.aln;
    const Path& path = aln.path();
    int64_t first = -1, last = -1;
    size_t from = 0, to = 0;
    if (read.indexed()) {
        // The source has already found the site's mappings, and the read offset before
        // each. Both ends come straight out of that, without touching the rest of the
        // read -- which for an ONT alignment is thousands of mappings the site does not
        // want. The index lists mappings in the queried ranges, which are exactly the
        // site's nodes, so no second membership test is needed.
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
            // This is the dominant gap path on ONT -- a homopolymer stutter inside a node
            // both the read and the allele visit.
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

// (node, orientation) packed into one integer, so membership is a binary search over a
// sorted vector rather than a hash lookup. The median traversal is three nodes, where a
// hash's constant costs more than the scan it replaces.
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

// The greedy walk: one left-to-right pass, anchoring on shared node visits and never revising a
// pairing once made. This is the DEFAULT, and it is the right default for short reads.
//
// The exact walk below buys +0.0042 chr20 / +0.0035 chr6 ONT indel F1 for +17% CPU, which is a
// good trade on long reads. On short reads the same machinery buys +0.0006 indel -- on BOTH
// contigs -- and nothing on SNVs, for 3.10x the CPU, because a 150 bp read barely diverges from an
// allele and there is almost nothing for an optimal correspondence to resolve. So --realign is off
// unless asked for, and --preset ont asks for it.
int32_t GraphAlignedAlleleLikelihoodCalculator::score_read_against_allele_greedy(
    const Alignment& aln, const vector<ReadStep>& read_steps,
    const vector<AlleleStep>& allele_steps, const EditAlignmentScorer& read_scorer,
    bool& placed_out, double& nat_adjust) const {

    placed_out = !allele_steps.empty();
    if (!placed_out) {
        return 0;
    }

    int32_t score = 0;
    size_t allele_index = 0;
    bool have_anchor = false;
    size_t bases_accounted = 0;

    // Last read position visiting each (node, orientation). The walk below has to ask
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

            score += score_shared_node(aln, read_step, read_scorer, nat_adjust);
            bases_accounted += read_step.read_length;
            have_anchor = true;
            allele_index = found + 1;
            continue;
        }

        // The allele has no matching node from here on. Either the read took a
        // node this allele lacks, or the two substituted nodes for one another.
        if (allele_index < allele_steps.size()) {
            const AlleleStep& allele_step = allele_steps[allele_index];

            // Does the read visit this allele node later on? Then it is not this read
            // node's counterpart -- the read took a node the allele simply lacks, which is
            // an insertion, not a substitution. Consuming the allele node here would burn
            // the anchor the read is about to need: its own visit would find the allele
            // exhausted and be charged a second time, so ONE inserted base cost a
            // substitution plus two gaps. The same event costs a single gap when it is the
            // allele that carries the extra node, and an indel must cost the same from
            // either side. Charge the insertion and leave allele_index where it is.
            auto later = read_last_visit.find(((int64_t)allele_step.node_id << 1)
                                              | (int64_t)allele_step.backward);
            if (later != read_last_visit.end() && later->second > read_index) {
                score += read_scorer.score_gap(read_step.read_length);
                nat_adjust += params.insertion_gap_nats;
                bases_accounted += read_step.read_length;
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
            ++allele_index;
        } else {
            // The allele is exhausted but the read continues. Those read bases
            // cannot be placed on this allele, so they are charged as an
            // insertion. They are NOT dropped: omitting them would score this
            // allele over fewer read bases than its competitors and fabricate a
            // likelihood ratio out of the length difference alone.
            score += read_scorer.score_gap(read_step.read_length);
            nat_adjust += params.insertion_gap_nats;
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

// The optimal walk (--realign, and --preset ont). Chooses the correspondence between the read's
// node visits and the allele's by dynamic programming rather than in one greedy pass. It does NOT
// re-align bases: base-level edits are still read off the mapper's alignment, and what is searched
// for is which read VISIT pairs with which allele VISIT. That distinction is why this is
// affordable where base-level realignment is not -- WFA over the same windows was measured at 24x
// this walk's CPU without finishing chr20 ONT. See doc/read-likelihood-genotyping.md,
// "The two walks".
//
// Four states over the m x n grid of read visits against allele visits, m and n being VISITS, not
// bases:
//
//   P  no shared visit matched yet -- the leading flank. Allele visits are FREE here, because
//      allele sequence outside the read's window is not the read's to explain, and charging it
//      would penalise a read for being short. Read visits are still charged. A match closes the
//      flank (P -> M); a leading substitution does not.
//   M  this read visit paired with this allele visit.
//   I  inside a run of read visits the allele lacks.
//   D  inside a run of allele visits the read lacks.
//
// Pairing cost, three cases: the same node visit costs the read's own edits inside it; two nodes
// where NEITHER appears anywhere in the other sequence is a genuine substitution; anything else is
// forbidden, which is the anchor rule -- a visit the two share may not be paired with a different
// node, though it may still be gapped. Both halves of that were measured: allowing the
// substitution costs 0.063 indel F1, and forbidding the gap costs a further 0.0071.
//
// Where this differs from the greedy walk in COST, not merely in correspondence: greedy charges a
// fresh gap open per inserted read visit, while I opens once and extends, so a k-visit insertion
// run differs by (k-1)*(gap_open - gap_extend). That is exactly zero at --preset ont, where both
// are 1, and non-zero at the short-read defaults of 6 and 1.
//
// The answer is the best final-row cell over P, M or I, never D -- D has charged trailing allele
// bases that lie outside the read's window. P is a legal terminal state, so a read that matches
// nothing still scores finitely.
int32_t GraphAlignedAlleleLikelihoodCalculator::score_read_against_allele(
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

    // A node the read and the allele SHARE may not be paired with a different node. In a
    // pangenome graph two paths through one node traverse identically the same bases, so a
    // shared visit is the alignment the mapper already asserted rather than a guess -- and
    // letting the correspondence decline it and substitute elsewhere costs 0.063 of indel F1
    // on chr20 ONT, because an insertion allele can then explain reads that do not carry the
    // insertion.
    //
    // It may still be GAPPED. That distinction is load-bearing and was measured: FORCING every
    // shared visit to match, rather than merely forbidding it to substitute, is 0.0060 of
    // indel F1 worse than this and slightly worse than the greedy walk this replaces. A read
    // whose own edits inside a shared node are bad -- an ONT homopolymer run -- is sometimes
    // better explained by gapping it.
    //
    // The predicate is membership, which is exact only while neither sequence repeats a node.
    // Real traversals do not: 0 repeats in 234,001 chr20 allele traversals and in 19,999 reads
    // averaging 936 nodes. Under a repeat it is merely over-restrictive -- it forbids a
    // substitution on account of an occurrence already consumed -- so it degrades optimality,
    // never correctness.
    const vector<int64_t>& read_keys = scratch.keys;
    auto holds = [](const vector<int64_t>& v, int64_t k) {
        return std::binary_search(v.begin(), v.end(), k);
    };

    // Four states. `nats` rides with `score` so the chosen path's --insertion-nats bookkeeping
    // is its own and not some other path's.
    //
    //   P  no shared visit matched yet. Allele nodes are FREE here: sequence before the first
    //      anchor is outside the read's window, and charging it would penalise a read for
    //      being short. A leading SUBSTITUTION does not close the flank.
    //   M  read step paired with allele step
    //   I  inside a run of read nodes the allele lacks
    //   D  inside a run of allele nodes the read lacks
    //
    // I and D are each reachable from the other. In base-level alignment that adjacency is
    // conventionally forbidden as a duplicate of a substitution, but here "substitution" means
    // base-level scoring of two different nodes, which is a different quantity -- forbidding it
    // lost real optima on 7,887 short-read cells.
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

    // Band the DP on big sites. A read and an allele differ over a handful of nodes, so the
    // correspondence hugs the diagonal the shared visits define; exploring the whole m x n
    // rectangle is wasted on a traversal of thousands of nodes. Unbanded this walk cost 1.93x the
    // CPU on short reads -- measured with the per-read hoisting held constant on both arms --
    // nearly all of it in that tail.
    //
    // The centre for read row i is projected back along the diagonal from the NEXT shared visit,
    // falling forward from the last one once none remains ahead -- so the band follows the anchors
    // rather than the i == j diagonal, which an indel would immediately push it off. Small sites -- the overwhelming majority, median
    // three nodes -- are left exhaustive, so the common case is untouched.
    //
    // It is an approximation, and measured as one: against the unbanded walk it moves a single
    // chr20 record of 115,255, leaves indel F1 identical at 0.86659 and SNV F1 marginally
    // better. Forcing perfect-match pairings instead -- a cheaper bound, and exact in a
    // standard alignment -- moves 43 records, because this walk is not a standard alignment:
    // allele sequence outside the read's window is free, so closing that flank early to take a
    // match can cost more than the match earns.
    const bool banded = m * n > 20000;
    const size_t band = 64;
    vector<size_t> centre;
    if (banded) {
        centre.assign(m + 1, 0);
        size_t a_index = 0, prev_i = 0, prev_j = 0;
        vector<pair<size_t, size_t>> shared;
        for (size_t i = 0; i < m; ++i) {
            for (size_t j = a_index; j < n; ++j) {
                if (allele_steps[j].node_id == read_steps[i].node_id &&
                    allele_steps[j].backward == read_steps[i].backward) {
                    shared.emplace_back(i, j);
                    a_index = j + 1;
                    break;
                }
            }
        }
        size_t s = 0;
        for (size_t i = 0; i <= m; ++i) {
            while (s < shared.size() && shared[s].first < i) {
                prev_i = shared[s].first;
                prev_j = shared[s].second;
                ++s;
            }
            centre[i] = s < shared.size()
                            ? shared[s].second - min(shared[s].second, shared[s].first - i)
                            : prev_j + (i - min(i, prev_i));
            centre[i] = min(centre[i], n);
        }
    }

    for (size_t i = 1; i <= m; ++i) {
        const ReadStep& rs = read_steps[i - 1];
        const int32_t rlen = (int32_t)rs.read_length;
        const int32_t rgap = read_scorer.score_gap(rs.read_length);
        const double rnat = params.insertion_gap_nats;
        const int64_t rkey = step_key(rs.node_id, rs.backward);

        P[0] = plus(pP[0], rgap, rnat);
        M[0] = NONE;
        I[0] = better(plus(better(pM[0], pD[0]), rgap, rnat),
                      plus(pI[0], rlen * extend_per_base, rnat));
        D[0] = NONE;

        size_t j_lo = 1, j_hi = n;
        if (banded) {
            const size_t mid = centre[i];
            j_lo = mid > band ? mid - band : 1;
            j_hi = min(n, mid + band);
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
            // Free leading allele deletion, and a read insertion before any anchor.
            P[j] = better(P[j], P[j - 1]);
            P[j] = better(P[j], plus(pP[j], rgap, rnat));

            // A shared visit may be GAPPED even though it may not be SUBSTITUTED away.
            // Forbidding the gap too -- forcing every shared visit to match, which is what a
            // pure partition of the two paths does -- costs 0.0071 of chr20 ONT indel F1,
            // because a read whose own bases inside a shared node are bad, typically an ONT
            // homopolymer run, is sometimes better explained by gapping the node than by
            // paying for those edits.
            I[j] = better(plus(better(pM[j], pD[j]), rgap, rnat),
                          plus(pI[j], rlen * extend_per_base, rnat));
            D[j] = better(plus(better(M[j - 1], I[j - 1]), agap, 0.0),
                          plus(D[j - 1], alen * extend_per_base, 0.0));
        }
        pP.swap(P); pM.swap(M); pI.swap(I); pD.swap(D);
    }

    // Every read step is consumed by exactly one transition, so every read base inside the
    // site is accounted for whatever the allele -- the invariant that stops one allele being
    // scored over a different span than its competitors.
    //
    // Allele nodes after the last match are outside the window and free, so simply stop:
    // never end in D, which has charged them.
    Cell best = NONE;
    for (size_t j = 0; j <= n; ++j) {
        best = better(best, better(pP[j], better(pM[j], pI[j])));
    }
    assert(best.score != NEG);
    nat_adjust += best.nats;
    return best.score;
}

GraphAlignedAlleleLikelihoodCalculator::WindowReadStats
GraphAlignedAlleleLikelihoodCalculator::local_read_stats(
    const vector<pair<nid_t, nid_t>>& site_ranges) const {

    // The neighbourhood the rate is measured over, in node IDs -- one constant for every read
    // source, deliberately decoupled from the source's fetch window. Taking the fetch window made
    // DR (and any --depth-quality GQ) depend on how the same reads were supplied: an in-memory
    // source fell back to 4096 while an indexed GAM's default fetch window is 256, sixteen times
    // narrower with differently quantized boundaries -- despite an older comment here claiming the
    // two matched. (The in-tree in-memory-vs-indexed test could not see the difference: its graph
    // fits inside one 256-ID window, where the two spans cover identical reads.) Measured through
    // --read-window on real data the width is not worth tuning: 1024 loses 0.004 structural-variant
    // F1, 16384 gains 0.001.
    static const size_t RATE_WINDOW = 4096;
    size_t span = RATE_WINDOW;
    if (site_ranges.empty() || params.depth_ploidy <= 0) {
        return WindowReadStats();
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
    double length_total = 0.0;
    size_t length_count = 0;
    read_source.for_each_alignment({{first, last}}, [&](const Alignment& aln) {
        // Only reads that BEGIN here. The fetch hands over everything overlapping the window, and
        // counting all of it is an overlap rate where the geometry wants a start rate; it is also
        // what makes the length mean size-biased. One test fixes both.
        const Path& path = aln.path();
        if (path.mapping_size() == 0) {
            return;
        }
        nid_t start_node = path.mapping(0).position().node_id();
        if (start_node < first || start_node > last) {
            return;
        }
        length_total += (double)aln.sequence().size();
        ++length_count;
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

    WindowReadStats stats;
    stats.start_rate = (reads <= 0.0 || bp == 0) ? 0.0 : reads / (double)bp;
    stats.mean_read_length = length_count > 0 ? length_total / (double)length_count : 0.0;
    lock_guard<std::mutex> guard(window_bp_mutex);
    window_rate[window_index] = stats;
    return stats;
}

AlleleReadLikelihoods GraphAlignedAlleleLikelihoodCalculator::compute(
    const Snarl& snarl, const vector<SnarlTraversal>& traversals, int ploidy) {


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

    // Anchor evidence, when anchors are armed. Filled per read alongside the matrix rows, and in
    // the same order: `add_read` drops a read that placed on nothing, so the evidence follows what
    // it kept rather than what was offered.
    unique_ptr<AnchorSiteEvidence> anchor_evidence;
    // Read phasing needs the same `rel` rows and nothing else -- no name, no pin. Built only when
    // anchors are NOT armed, because the anchor evidence already carries everything it wants.
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

    // Counted into the run's instance where there is one. A local stands in otherwise -- a unit
    // test driving the calculator directly has no run to count into, and a dropped count is better
    // than a singleton nobody can see. Hoisted out of the per-read callback: it is 20 atomics, and
    // constructing them once per read to throw them away is pure waste.
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

        // Flip a reverse-strand read into the alleles' reading direction. Without
        // this it anchors on nothing, falls through to the substitution path, and
        // is scored against the wrong allele entirely.
        //
        // Only the stretch inside the site is flipped -- see site_span_of. Flipping the whole
        // alignment is the same answer at, on long reads, many times the cost.
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

        // Depends on the read alone, so it is computed here rather than per allele. Only the
        // exact walk consumes it, and building it is not free, so skip it for the greedy one.
        if (params.realign) {
            prepare_read_scratch(*scored_aln, read_steps, read_scorer, read_scratch);
        }

        for (size_t a = 0; a < traversals.size(); ++a) {
            bool placed = false;
            double nat_adjust = 0.0;
            int32_t score =
                params.realign
                    ? score_read_against_allele(*scored_aln, read_steps, allele_steps[a],
                                                read_scratch, allele_keys[a], read_scorer,
                                                placed, nat_adjust)
                    : score_read_against_allele_greedy(*scored_aln, read_steps, allele_steps[a],
                                                       read_scorer, placed, nat_adjust);
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

        if (!builder.add_read(row, mismap, aln.name(), aln.sequence().size())) {
            return;
        }

        if (phase_evidence != nullptr) {
            // Hashed here, where the name is in hand, and never stored: the string is the whole
            // reason the anchor evidence is too heavy to retain for phasing.
            phase_evidence->read_key.push_back((uint64_t)std::hash<string>{}(aln.name()));
        }
        if (anchor_evidence != nullptr) {
            // Resolved against `aln`, NOT `scored_aln`. The flipped copy exists so the read can be
            // compared to alleles that read the other way; its sequence is reverse-complemented and
            // its offsets are in that frame, so pinning against it would report positions in a read
            // the read file does not contain. The strand field is what expresses orientation here.
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

    // R, for every consumer of it, from one place. The builder can only offer the mean over the
    // reads delivered to THIS site, and that estimator is size-biased -- a long read overlaps more
    // sites, so it is sampled more often, and the site mean estimates E[L^2]/E[L] rather than E[L]
    // (78,165 against 33,449 bp on a 33 kb ONT set; equal at fixed read length, which is why this
    // was invisible on 150 bp reads). Every use of R here is Lander-Waterman geometry -- how much
    // sequence can produce a read that overlaps something -- and all of it wants the population
    // mean. Three consumers: the depth term's lambda below, the mixture weights
    // (`set_length_weights`, via this same field), and the anchor slot weights, which read it back
    // through `mean_read_length_estimate` a few lines down. Fixing only the depth term would leave
    // the three disagreeing about one quantity.
    WindowReadStats stats = local_read_stats(ranges);
    if (stats.mean_read_length > 0.0) {
        result.set_mean_read_length(stats.mean_read_length);
    }

    if (anchor_evidence != nullptr) {
        // The rows the builder kept, in the order it kept them, so row r is reads[r].
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
        // callback, after `add_read` has accepted the read, so a read it drops contributes neither.
        phase_evidence->mean_read_length = (float)result.mean_read_length_estimate();
        // Fail CLOSED on a row-count mismatch rather than resizing into one. `resize` here would
        // truncate or zero-pad, and a zero key is a key -- every padded read would collide and be
        // silently treated as the same fragment at every site. Dropping the site's evidence costs
        // one site's phase; mis-keying costs a wrong haplotype with nothing to show for it.
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
        // Set unconditionally so `DR` is emitted whether or not the term is armed:
        // the observable should be measurable as a ranking signal before the model
        // is allowed to act on it. A zero weight leaves the likelihood untouched.
        // Per haplotype, from the site's own ploidy rather than an assumed one.
        int effective_ploidy = ploidy > 0 ? ploidy : params.depth_ploidy;
        result.set_depth_context(depth_lengths,
                                 stats.start_rate / (double)effective_ploidy,
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

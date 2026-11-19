#include "recombinator.hpp"

#include "../artifacts/diploid_pair_scorer.hpp"
#include "gref.hpp"
#include "statistics.hpp"
#include "algorithms/component.hpp"
#include "algorithms/extract_subchain.hpp"
#include "hash_map.hpp"
#include "utility.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

//#define debug_selected_haplotypes

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t HaplotypePartitioner::APPROXIMATE_JOBS;

constexpr size_t Recombinator::NUM_HAPLOTYPES;
constexpr size_t Recombinator::NUM_CANDIDATES;
constexpr size_t Recombinator::COVERAGE;
constexpr double Recombinator::PRESENT_DISCOUNT;
constexpr double Recombinator::HET_ADJUSTMENT;
constexpr double Recombinator::ABSENT_SCORE;
constexpr size_t Recombinator::LONG_READ_HOMOPOLYMERS;
constexpr size_t Recombinator::MIN_GREF_LENGTH;
constexpr double Recombinator::MODEL_TAU;
constexpr double Recombinator::MODEL_RHO;

// Other static class variables.

const std::vector<std::pair<Recombinator::Parameters::preset_t, std::string>> Recombinator::Parameters::PRESETS = {
    { Recombinator::Parameters::preset_default, "default" },
    { Recombinator::Parameters::preset_haploid, "haploid" },
    { Recombinator::Parameters::preset_diploid, "diploid" },
    { Recombinator::Parameters::preset_hifi, "hifi" },
    { Recombinator::Parameters::preset_ont, "ont" },
    { Recombinator::Parameters::preset_model, "model" }
};

//------------------------------------------------------------------------------

// Returns true if the path belongs to a gref cover: a copy of a reference path
// (`gref_CHM13#0#chr1`, or `gref_x` for a generic path) or a fragment of one.
bool is_gref_path(const gbwt::FullPathName& path_name) {
    if (path_name.sample_name == gbwtgraph::GENERIC_PATH_SAMPLE_NAME) {
        return GrefCover::is_gref_derived(path_name.contig_name);
    }
    return GrefCover::is_gref_derived(path_name.sample_name);
}

//------------------------------------------------------------------------------

HaplotypePartitioner::HaplotypePartitioner(
    const gbwtgraph::GBZ& gbz,
    const gbwt::FastLocate& r_index,
    const SnarlDistanceIndex& distance_index,
    const minimizer_index_type& minimizer_index,
    Verbosity verbosity
) :
    gbz(gbz), fragment_map(gbz.index.metadata, verbosity >= Haplotypes::verbosity_extra_debug), r_index(r_index),
    distance_index(distance_index), minimizer_index(minimizer_index),
    verbosity(verbosity)
{
    if (this->verbosity >= Haplotypes::verbosity_detailed) {
        std::cerr << "HaplotypePartitioner: " << this->gbz.index.metadata.paths() << " fragments for " << this->fragment_map.size() << " haplotype sequences" << std::endl;
    }
}

void HaplotypePartitioner::Parameters::print(std::ostream& out) const {
    out << "Partitioning parameters:" << std::endl;
    out << "- target length " << this->subchain_length << " bp" << std::endl;
    if (this->linear_structure) {
        out << "- strictly linear structure" << std::endl;
    }
    if (this->subchain_kmers) {
        out << "- subchain kmers with copies" << std::endl;
    }
    out << "- " << this->approximate_jobs << " jobs" << std::endl;
}

//------------------------------------------------------------------------------

Haplotypes HaplotypePartitioner::partition_haplotypes(const Parameters& parameters) const {

    // Sanity checks.
    if (parameters.subchain_length == 0) {
        std::string msg = "HaplotypePartitioner::partition_haplotypes(): subchain length cannot be 0";
        throw std::runtime_error(msg);
    }
    if (parameters.approximate_jobs == 0) {
        std::string msg = "HaplotypePartitioner::partition_haplotypes(): number of jobs cannot be 0";
        throw std::runtime_error(msg);
    }
    if (this->verbosity >= Haplotypes::verbosity_detailed) {
        parameters.print(std::cerr);
    }

    // TODO: Should this be a constructor?
    Haplotypes result;
    result.header.k = this->minimizer_index.k();
    result.set_graph_name(this->gbz);

    // Determine GBWT construction jobs.
    double start = gbwt::readTimer();
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Determining construction jobs" << std::endl;
    }
    size_t size_bound = this->gbz.graph.get_node_count() / parameters.approximate_jobs;
    gbwtgraph::ConstructionJobs jobs = gbwtgraph::gbwt_construction_jobs(this->gbz.graph, size_bound);
    result.header.construction_jobs = jobs.size();

    // Determine the number of top-level chains and fill in basic information.
    size_t total_chains = 0;
    this->distance_index.for_each_child(this->distance_index.get_root(), [&](const handlegraph::net_handle_t& chain) {
        if (this->distance_index.is_looping_chain(chain)) {
            std::string msg = "HaplotypePartitioner::partition_haplotypes(): top-level chain " +
                std::to_string(total_chains) + " is a loop; haplotype sampling cannot be used with this graph\n." +
                "This can sometimes be resolved by using the vg index -P option to specify a reference backbone " +
                "when computing the distance index.";
            throw std::runtime_error(msg);
        }
        total_chains++;
    });
    if (jobs.components() != total_chains) {
        // TODO: Could we instead identify the components with multiple top-level chains
        // and skip them?
        std::string msg = "HaplotypePartitioner::partition_haplotypes(): there are "
            + std::to_string(total_chains) + " top-level chains and " + std::to_string(jobs.components())
            + " weakly connected components; haplotype sampling cannot be used with this graph";
        throw std::runtime_error(msg);
    }
    result.header.top_level_chains = jobs.components();
    result.chains.resize(result.components());
    for (size_t chain_id = 0; chain_id < result.components(); chain_id++) {
        result.chains[chain_id].offset = chain_id;
        result.chains[chain_id].job_id = result.jobs(); // Not assigned to any job yet.
        result.chains[chain_id].placeholder = 0; // Placeholder for the default scoring model.
        result.chains[chain_id].contig_name = "component_" + std::to_string(chain_id); // Placeholder name.
    }

    // Assign chains to jobs and fill in contig names.
    auto chains_by_job = gbwtgraph::partition_chains(this->distance_index, this->gbz.graph, jobs);
    // Chain names become the contigs of the generated haplotypes, so a gref cover must not
    // change them. Skip gref fragments (`chr1_7_alt`) and loci in the reserved gref namespace
    // (`gref_x`, a non-PanSN copy); a PanSN gref copy (`gref_CHM13#0#chr1`) has the base locus.
    auto contig_names = jobs.contig_names(this->gbz.graph, [&](const path_handle_t& path) -> bool {
        std::string path_name = this->gbz.graph.get_path_name(path);
        return !GrefCover::is_gref_derived(this->gbz.graph.get_locus_name(path)) && !GrefCover::is_gref_fragment(path_name);
    });
    for (size_t job_id = 0; job_id < result.jobs(); job_id++) {
        for (auto& chain : chains_by_job[job_id]) {
            result.chains[chain.offset].job_id = job_id;
            nid_t node_id = this->gbz.graph.get_id(chain.handle);
            size_t component_id = jobs.component(node_id);
            if (component_id >= contig_names.size()) {
                std::string msg = "HaplotypePartitioner::partition_haplotypes(): cannot map top-level chain " + std::to_string(chain.offset) + " to a component";
                throw std::runtime_error(msg);
            }
            result.chains[chain.offset].contig_name = contig_names[component_id];
            if (this->verbosity >= Haplotypes::verbosity_debug) {
                std::cerr << "Using contig name " << contig_names[component_id] << " for chain " << chain.offset << std::endl;
            }
        }
    }

    // Assign paths to construction jobs. We do this for all paths, as reference status
    // is not a fixed property of a path.
    size_t total_paths = this->gbz.index.metadata.paths();
    result.jobs_for_paths = std::vector<size_t>(total_paths, result.jobs());
    for (size_t path_id = 0; path_id < total_paths; path_id++) {
        gbwt::size_type gbwt_seq_id = gbwt::Path::encode(path_id, false);
        gbwt::edge_type start = this->gbz.index.start(gbwt_seq_id);
        nid_t node_id = gbwt::Node::id(start.first);
        result.jobs_for_paths[path_id] = jobs.job(node_id);
    }

    jobs = gbwtgraph::ConstructionJobs(); // Save memory.
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Partitioned " << result.components() << " components into " << result.jobs() << " jobs in " << seconds << " seconds" << std::endl;
    }

    // Determine the subchains and sequences for each top-level chain.
    if (verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Running " << omp_get_max_threads() << " jobs in parallel" << std::endl;
    }
    start = gbwt::readTimer();
    // With a single job, building subchain kmers can use the threads within the chains.
    #pragma omp parallel for schedule(dynamic, 1) if(chains_by_job.size() > 1)
    for (size_t job = 0; job < chains_by_job.size(); job++) {
        const std::vector<gbwtgraph::TopLevelChain>& chains = chains_by_job[job];
        size_t total_subchains = 0, total_kmers = 0;
        size_t total_parents = 0, total_sequences = 0;
        for (auto& chain : chains) {
            try {
                this->build_subchains(chain, result.chains[chain.offset], parameters);
            } catch (const std::runtime_error& e) {
                std::cerr << "error: [job " << job << "]: " << e.what() << std::endl;
                std::exit(EXIT_FAILURE);
            }
            total_subchains += result.chains[chain.offset].subchains.size();
            for (auto& subchain : result.chains[chain.offset].subchains) {
                total_kmers += subchain.num_kmers();
                total_parents += subchain.kmers_present.get_num_parents();
                total_sequences += subchain.num_sequences();
            }
        }
        #pragma omp critical
        {
            result.header.total_subchains += total_subchains;
            result.header.total_kmers += total_kmers;
            if (this->verbosity >= Haplotypes::verbosity_detailed) {
                std::cerr << "Finished job " << job << " with " << chains.size() << " chains (";
                for (size_t i = 0; i < chains.size(); i++) {
                    if (i > 0) {
                        std::cerr << ", ";
                    }
                    std::cerr << result.chains[chains[i].offset].contig_name;
                }
                std::cerr << "), " << total_subchains << " subchains, " << total_kmers << " kmers, and "
                    << total_parents << " parents for " << total_sequences << " sequences" << std::endl;
            }
        }
    }
    if (verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Finished the jobs in " << seconds << " seconds" << std::endl;
    }

    return result;
}

//------------------------------------------------------------------------------

size_t HaplotypePartitioner::get_distance(handle_t from, handle_t to) const {
    return this->distance_index.minimum_distance(
        this->gbz.graph.get_id(from), this->gbz.graph.get_is_reverse(from), this->gbz.graph.get_length(from) - 1,
        this->gbz.graph.get_id(to), this->gbz.graph.get_is_reverse(to), 0,
        false, &this->gbz.graph
    );
}

bool HaplotypePartitioner::contains_reversals(handle_t handle) const {
    gbwt::node_type forward = gbwtgraph::GBWTGraph::handle_to_node(handle);
    std::vector<gbwt::size_type> forward_da = this->r_index.decompressDA(forward);
    std::sort(forward_da.begin(), forward_da.end());

    gbwt::node_type reverse = gbwt::Node::reverse(forward);
    std::vector<gbwt::size_type> reverse_da = this->r_index.decompressDA(reverse);
    std::sort(reverse_da.begin(), reverse_da.end());

    auto fw_iter = forward_da.begin();
    auto rv_iter = reverse_da.begin();
    while (fw_iter != forward_da.end() && rv_iter != reverse_da.end()) {
        if (*fw_iter == *rv_iter) {
            return true;
        } else if (*fw_iter < *rv_iter) {
            ++fw_iter;
        } else {
            ++rv_iter;
        }
    }
    return false;
}

std::vector<HaplotypePartitioner::Subchain>
HaplotypePartitioner::get_subchains(const gbwtgraph::TopLevelChain& chain, const Parameters& parameters) const {
    std::vector<Subchain> result;

    // First pass: take all snarls as subchains.
    std::vector<Subchain> snarls;
    handle_t snarl_start = empty_gbwtgraph_handle();
    bool has_start = false;
    bool was_snarl = false;
    // Closed interval of net handles, from the start bound to the end bound.
    // We list the children of the chain instead of stepping with
    // follow_net_edges(), which gives only traversals that a walk can take and
    // so cannot step over a snarl that no walk crosses. Such snarls are split
    // below.
    net_handle_t chain_end = this->distance_index.get_bound(chain.chain, true, false);
    this->distance_index.for_each_child(chain.chain, [&](const net_handle_t& curr) -> bool {
        if (this->distance_index.is_node(curr)) {
            handle_t handle = this->distance_index.get_handle(curr, &this->gbz.graph);
            if (was_snarl) {
                if (!has_start) {
                    // If the chain starts with a snarl, we take it as a prefix.
                    snarls.push_back({ Haplotypes::Subchain::prefix, empty_gbwtgraph_handle(), handle, 0, 0 });
                } else {
                    size_t distance = this->get_distance(snarl_start, handle);
                    if (distance < std::numeric_limits<size_t>::max()) {
                        // Normal snarl with two boundary nodes.
                        snarls.push_back({ Haplotypes::Subchain::normal, snarl_start, handle, 0, 0 });
                    } else {
                        // The snarl is not connected, so we break it into two.
                        snarls.push_back({ Haplotypes::Subchain::suffix, snarl_start, empty_gbwtgraph_handle(), 0, 0 });
                        snarls.push_back({ Haplotypes::Subchain::prefix, empty_gbwtgraph_handle(), handle, 0, 0 });
                    }
                }
            }
            snarl_start = handle;
            has_start = true;
            was_snarl = false;
        } else if (this->distance_index.is_snarl(curr)) {
            was_snarl = true;
        }

        // We check the loop condition here, as we have an inclusive endpoint.
        return curr != chain_end;
    });
    if (was_snarl && has_start) {
        // If the chain ends with a snarl, we take it as a suffix.
        snarls.push_back({ Haplotypes::Subchain::suffix, snarl_start, empty_gbwtgraph_handle(), 0, 0 });
    }

    // Second pass: Combine snarls into subchains.
    size_t head = 0;
    while (head < snarls.size()) {
        if (snarls[head].type != Haplotypes::Subchain::normal) {
            // Prefixes and suffixes should be rare, so we can simply use them directly.
            result.push_back(snarls[head]);
            head++;
            continue;
        }
        size_t tail = head;
        std::uint32_t distance = this->get_distance(snarls[head].start, snarls[tail].end);
        std::uint32_t extra_snarls = 0;
        while (tail + 1 < snarls.size()) {
            if (snarls[tail + 1].type != Haplotypes::Subchain::normal) {
                break;
            }
            size_t candidate = this->get_distance(snarls[head].start, snarls[tail + 1].end);
            if (candidate > parameters.subchain_length) {
                // Including the next snarl would exceed target length. But if a haplotype visits
                // the tail in both orientations, it flips the orientation in a subsequent subchain,
                // returns back, flips again, and eventually continues forward. In such situations,
                // sampling minimal haplotypes within this subchain would lead to sequence loss,
                // while sampling maximal haplotypes could make some kmers specific to the next
                // subchain shared with haplotypes in this subchain. We therefore move forward until
                // we can make the subchain contain the reversals.
                if (parameters.linear_structure && this->contains_reversals(snarls[tail].end)) {
                    extra_snarls++;
                } else {
                    break;
                }
            }
            tail++;
            distance = candidate;
        }
        result.push_back({
            Haplotypes::Subchain::normal,
            snarls[head].start, snarls[tail].end,
            distance, extra_snarls
        });
        head = tail + 1;
    }

    return result;
}

//------------------------------------------------------------------------------

// Returns (SA[i], i) for the handle, ordered by i.
std::vector<HaplotypePartitioner::sequence_type> get_sequence_visits(handle_t handle, const gbwt::FastLocate& r_index) {
    std::vector<gbwt::size_type> sa = r_index.decompressSA(gbwtgraph::GBWTGraph::handle_to_node(handle));
    std::vector<HaplotypePartitioner::sequence_type> result;
    result.reserve(sa.size());
    for (size_t i = 0; i < sa.size(); i++) {
        result.push_back({ sa[i], i });
    }
    return result;
}

std::vector<HaplotypePartitioner::sequence_type> HaplotypePartitioner::get_sequences(handle_t handle) const {
    auto result = get_sequence_visits(handle, this->r_index);
    std::sort(result.begin(), result.end(), [&](sequence_type a, sequence_type b) -> bool {
        gbwt::size_type a_id = this->r_index.seqId(a.first);
        gbwt::size_type a_offset = this->r_index.seqOffset(a.first);
        gbwt::size_type b_id = this->r_index.seqId(b.first);
        gbwt::size_type b_offset = this->r_index.seqOffset(b.first);
        return ((a_id < b_id) || ((a_id == b_id) && (a_offset > b_offset)));
    });
    for (auto& sequence : result) {
        sequence.first = this->r_index.seqId(sequence.first);
    }
    return result;
}

// A sequence visit (SA[i], i) decomposed into a visit in a fragmented haplotype.
// Visits can be sorted by (chain, orientation, position within the chain).
struct FragmentedHaplotypeVisit {
    // (DA[i], i).
    gbwt::size_type sequence_id;
    gbwt::size_type gbwt_node_offset;

    // Information used for ordering visits.
    gbwt::size_type oriented_chain_id;
    gbwt::size_type rank_in_chain;
    gbwt::size_type rank_in_sequence;

    FragmentedHaplotypeVisit(HaplotypePartitioner::sequence_type sequence_visit, const HaplotypePartitioner& partitioner) {
        // Sequence id is encoded in the SA value and node offset is already given.
        this->sequence_id = partitioner.r_index.seqId(sequence_visit.first);
        this->gbwt_node_offset = sequence_visit.second;

        // We get the chain id from the fragment map and the orientation from the sequence id.
        // Then we abuse the path/sequence name encoding to encode the orientation in the id.
        gbwt::size_type path_id = gbwt::Path::id(this->sequence_id);
        bool is_reverse = gbwt::Path::is_reverse(this->sequence_id);
        this->oriented_chain_id = gbwt::Path::encode(partitioner.fragment_map.chain(path_id), is_reverse);

        // Rank in chain is based on the count/fragment/starting offset field in path name.
        // But we must consider the orientation of the chain.
        gbwt::PathName path_name = partitioner.gbz.index.metadata.path(path_id);
        this->rank_in_chain = path_name.count;
        if (is_reverse) {
            this->rank_in_chain = std::numeric_limits<gbwt::size_type>::max() - this->rank_in_chain;
        }

        // Sequence offset in the SA value is the distance to the end of the sequence in nodes.
        this->rank_in_sequence = std::numeric_limits<gbwt::size_type>::max() - partitioner.r_index.seqOffset(sequence_visit.first);
    }

    // Returns true if the visits are in the same chain in the same orientation.
    bool same_chain(const FragmentedHaplotypeVisit& another) const {
        return (this->oriented_chain_id == another.oriented_chain_id);
    }

    // Orders the visits by chain and orientation, and visits within the chain
    // by their position in the chain.
    bool operator<(const FragmentedHaplotypeVisit& another) const {
        if (this->oriented_chain_id != another.oriented_chain_id) {
            return (this->oriented_chain_id < another.oriented_chain_id);
        } else if (this->rank_in_chain != another.rank_in_chain) {
            return (this->rank_in_chain < another.rank_in_chain);
        } else {
            return (this->rank_in_sequence < another.rank_in_sequence);
        }
    }

    // Returns (DA[i], i).
    HaplotypePartitioner::sequence_type to_sequence() const {
        return HaplotypePartitioner::sequence_type(this->sequence_id, this->gbwt_node_offset);
    }
};

// Matches the minimal sequence visits in `from` and `to` by the same oriented chain of haplotype sequences.
// Requires that the fragments remain within the subchain defined by start and end.
// Returns (DA[i], i) for the matching visits in `from`, sorted by rank in the chain.
std::vector<HaplotypePartitioner::sequence_type> match_visits(
    handle_t start, handle_t end,
    const std::vector<HaplotypePartitioner::sequence_type>& from_visits,
    const std::vector<HaplotypePartitioner::sequence_type>& to_visits,
    const HaplotypePartitioner& partitioner
) {
    std::vector<FragmentedHaplotypeVisit> from, to;
    for (auto& visit : from_visits) {
        from.push_back(FragmentedHaplotypeVisit(visit, partitioner));
    }
    std::sort(from.begin(), from.end());
    for (auto& visit : to_visits) {
        to.push_back(FragmentedHaplotypeVisit(visit, partitioner));
    }
    std::sort(to.begin(), to.end());

    // When there are multiple fragments, some of them could be outside the subchain.
    // We want to avoid that to keep the kmers we use specific to the subchain.
    hash_set<nid_t> node_ids = extract_subchain(partitioner.gbz.graph, start, end);
    auto within_subchain = [&](const FragmentedHaplotypeVisit& from_visit, const FragmentedHaplotypeVisit& to_visit) -> bool {
        gbwt::size_type sequence_id = from_visit.sequence_id;
        gbwt::edge_type pos(gbwtgraph::GBWTGraph::handle_to_node(start), from_visit.gbwt_node_offset);
        while (sequence_id != to_visit.sequence_id) {
            while (pos.first != gbwt::ENDMARKER) {
                if (node_ids.find(gbwt::Node::id(pos.first)) == node_ids.end()) {
                    // The fragment is outside the subchain.
                    return false;
                }
                pos = partitioner.gbz.index.LF(pos);
            }
            sequence_id = partitioner.fragment_map.oriented_next(sequence_id);
            if (sequence_id == gbwt::invalid_sequence()) {
                // No more fragments in the chain.
                return false;
            }
            pos = partitioner.gbz.index.start(sequence_id);
            if (pos.first != gbwt::ENDMARKER && node_ids.find(gbwt::Node::id(pos.first)) == node_ids.end()) {
                // The fragment starts outside the subchain.
                return false;
            }
        }

        // We have ensured that this is a minimal end-to-end visit. Once the sequence ids match
        // and we know that the position is within the subchain, we cannot leave the subchain and
        // later reach the end without re-entering the subchain.
        return true;
    };

    std::vector<HaplotypePartitioner::sequence_type> result;
    auto from_iter = from.begin();
    auto to_iter = to.begin();
    while (from_iter != from.end() && to_iter != to.end()) {
        if (*from_iter < *to_iter) {
            if (from_iter->same_chain(*to_iter)) {
                auto peek = from_iter + 1;
                while (peek != from.end() && *peek < *to_iter) {
                    from_iter = peek;
                    ++peek;
                }
                if (within_subchain(*from_iter, *to_iter)) {
                    result.push_back(from_iter->to_sequence());
                }
                ++from_iter; ++to_iter;
            } else {
                ++from_iter;
            }
        } else {
            ++to_iter;
        }
    }

    return result;
}

std::vector<HaplotypePartitioner::sequence_type> HaplotypePartitioner::get_sequences(Subchain subchain) const {
    if (subchain.type == Haplotypes::Subchain::prefix) {
        // NOTE: If a sequence flips after this subchain and returns to it, we will
        // take the initial prefix that exits the subchain multiple times.
        return this->get_sequences(subchain.end);
    }
    if (subchain.type == Haplotypes::Subchain::suffix) {
        // NOTE: If the sequence flips, exits the subchain, and later returns to it,
        // some of the haplotypes will overlap the preceding subchain(s).
        return this->get_sequences(subchain.start);
    }
    auto from = get_sequence_visits(subchain.start, this->r_index);
    auto to = get_sequence_visits(subchain.end, this->r_index);
    return match_visits(subchain.start, subchain.end, from, to, *this);
}

//------------------------------------------------------------------------------

// Generate a haplotype over the closed range from `pos` to `end`.
// The haplotype may consist of multiple fragments.
// Take at most start_max and end_max characters from the initial and the final
// node, respectively.
// Returns an empty haplotype if there is only one node.
// Set `pos.first = gbwt::ENDMARKER` to start from the beginning.
// Set `end = empty_gbwtgraph_handle()` to continue until the end without a final node.
std::vector<std::string> generate_haplotype(
    gbwt::size_type sequence_id,
    gbwt::edge_type pos, handle_t end, size_t start_max, size_t end_max,
    const HaplotypePartitioner& partitioner, gbwt::CachedGBWT& cached_gbwt
) {
    // Determine the start for a prefix.
    bool single_fragment = (pos.first == gbwt::ENDMARKER || end == empty_gbwtgraph_handle());
    if (pos.first == gbwt::ENDMARKER) {
        pos = cached_gbwt.start(sequence_id);
    }

    // Empty haplotype.
    std::vector<std::string> haplotype;
    if (pos == gbwt::invalid_edge() || pos.first == gbwt::ENDMARKER) {
        return haplotype;
    }

    // Handle the initial node.
    handle_t curr = gbwtgraph::GBWTGraph::node_to_handle(pos.first);
    if (curr == end) {
        return haplotype;
    }
    std::string_view view = partitioner.gbz.graph.get_sequence_view(curr);
    size_t offset = (view.size() > start_max ? view.size() - start_max : 0);
    haplotype.emplace_back();
    haplotype.back().append(view.data() + offset, view.size() - offset);

    while (true) {
        pos = cached_gbwt.LF(pos);
        if (pos.first == gbwt::ENDMARKER) {
            if (single_fragment) {
                break;
            }
            // Try to continue with the next non-empty fragment.
            while (pos.first == gbwt::ENDMARKER) {
                sequence_id = partitioner.fragment_map.oriented_next(sequence_id);
                if (sequence_id == gbwt::invalid_sequence()) {
                    break;
                }
                pos = cached_gbwt.start(sequence_id);
            }
            if (sequence_id == gbwt::invalid_sequence()) {
                break;
            }
            haplotype.emplace_back();
        }
        curr = gbwtgraph::GBWTGraph::node_to_handle(pos.first);
        view = partitioner.gbz.graph.get_sequence_view(curr);
        if (curr == end) {
            haplotype.back().append(view.data(), std::min(view.size(), end_max));
            break;
        } else {
            haplotype.back().append(view.data(), view.size());
        }
    }

    return haplotype;
}

// Return the ids of the nodes the sequence visits from node `start`, at the sequence's
// offset there, to node `end`, sorted and without duplicates. Without a start
// (`gbwt::ENDMARKER`), the walk begins at the start of the sequence, and without an end,
// it continues to the end of the sequence. With both, the walk continues with the next
// fragment of the haplotype when a fragment ends, because the generated haplotypes do
// the same.
std::vector<nid_t> subchain_walk_nodes(
    const gbwt::CachedGBWT& index, gbwt::node_type start, gbwt::node_type end,
    Haplotypes::sequence_type sequence, const gbwt::FragmentMap& fragment_map
) {
    gbwt::edge_type pos = (start != gbwt::ENDMARKER ? gbwt::edge_type(start, sequence.second) : index.start(sequence.first));
    std::vector<nid_t> nodes;
    bool normal = (start != gbwt::ENDMARKER && end != gbwt::ENDMARKER);
    gbwt::size_type sequence_id = sequence.first;
    while (pos.first != gbwt::ENDMARKER) {
        nodes.push_back(gbwt::Node::id(pos.first));
        if (pos.first == end) {
            break;
        }
        pos = index.LF(pos);
        while (normal && pos.first == gbwt::ENDMARKER) {
            sequence_id = fragment_map.oriented_next(sequence_id);
            if (sequence_id == gbwt::invalid_sequence()) {
                break;
            }
            pos = index.start(sequence_id);
        }
    }
    gbwt::removeDuplicates(nodes, false);
    return nodes;
}

// Generate the haplotype of the sequence over the subchain. Keep only the last k - 1 bases
// of the start node and the first k - 1 bases of the end node, where k is the kmer length
// of the minimizer index, so that every kmer in the haplotype has a base past a boundary
// node.
std::vector<std::string> subchain_haplotype(
    HaplotypePartitioner::sequence_type sequence, const HaplotypePartitioner::Subchain& subchain,
    const HaplotypePartitioner& partitioner, gbwt::CachedGBWT& cached_gbwt
) {
    size_t k = partitioner.minimizer_index.k();
    gbwt::edge_type pos(gbwt::ENDMARKER, 0);
    size_t start_max = std::numeric_limits<size_t>::max();
    if (subchain.has_start()) {
        pos = gbwt::edge_type(gbwtgraph::GBWTGraph::handle_to_node(subchain.start), sequence.second);
        start_max = k - 1;
    }
    return generate_haplotype(sequence.first, pos, subchain.end, start_max, k - 1, partitioner, cached_gbwt);
}


// Return the sorted set of kmers that are minimizers in the sequences and have a single
// occurrence in the graph.
std::vector<HaplotypePartitioner::kmer_type> take_unique_minimizers(
    const std::vector<std::string>& sequences,
    const HaplotypePartitioner::minimizer_index_type& minimizer_index,
    vg::hash_map<HaplotypePartitioner::kmer_type, size_t>& kmer_cache
) {
    std::vector<HaplotypePartitioner::kmer_type> result;
    for (const std::string& sequence : sequences) {
        auto minimizers = minimizer_index.minimizers(sequence);
        for (auto& minimizer : minimizers) {
            auto iter = kmer_cache.find(minimizer.key.get_key());
            if (iter != kmer_cache.end()) {
                if (iter->second == 1) {
                    result.push_back(minimizer.key.get_key());
                }
            } else {
                size_t count = minimizer_index.count(minimizer);
                kmer_cache[minimizer.key.get_key()] = count;
                if (count == 1) {
                    result.push_back(minimizer.key.get_key());
                }

            }
        }
    }
    gbwt::removeDuplicates(result, false);
    return result;
}

std::vector<HaplotypePartitioner::kmer_type> HaplotypePartitioner::unique_minimizers(
    gbwt::size_type sequence_id,
    gbwt::CachedGBWT& cached_gbwt,
    vg::hash_map<HaplotypePartitioner::kmer_type, size_t>& kmer_cache
) const {
    size_t limit = std::numeric_limits<size_t>::max();
    std::vector<std::string> haplotype = generate_haplotype(
        sequence_id,
        gbwt::edge_type(gbwt::ENDMARKER, 0), empty_gbwtgraph_handle(), limit, limit,
        *this, cached_gbwt
    );
    return take_unique_minimizers(haplotype, this->minimizer_index, kmer_cache);
}

std::vector<HaplotypePartitioner::kmer_type> HaplotypePartitioner::unique_minimizers(
    sequence_type sequence, Subchain subchain, size_t& fragments,
    gbwt::CachedGBWT& cached_gbwt,
    vg::hash_map<HaplotypePartitioner::kmer_type, size_t>& kmer_cache
) const {
    std::vector<std::string> haplotype = subchain_haplotype(sequence, subchain, *this, cached_gbwt);
    fragments = haplotype.size();
    return take_unique_minimizers(haplotype, this->minimizer_index, kmer_cache);
}

//------------------------------------------------------------------------------

/*
  Take a set of sequences defined by the sorted set of kmers that are present, as
  well as their ids. Build the `kmers`, `kmer_counts`, `sequences`, and
  `kmers_present` fields in the subchain, which are assumed to be empty.
*/
void present_kmers(
    const std::vector<std::vector<HaplotypePartitioner::kmer_type>>& kmers_by_sequence,
    const std::vector<Haplotypes::sequence_type>& sequences,
    Haplotypes::Subchain& subchain
) {
    // First pass: Determine the distinct informative kmers and their counts.
    size_t max_count = 0;
    std::vector<std::vector<HaplotypePartitioner::kmer_type>::const_iterator> iterators; 
    for (const auto& kmers : kmers_by_sequence) {
        iterators.push_back(kmers.begin());
    }
    std::vector<size_t> kmer_counts;
    while (true) {
        HaplotypePartitioner::kmer_type min_value = std::numeric_limits<HaplotypePartitioner::kmer_type>::max();
        size_t count = 0;
        for (size_t i = 0; i < kmers_by_sequence.size(); i++) {
            if (iterators[i] == kmers_by_sequence[i].end()) {
                continue;
            }
            if (*iterators[i] < min_value) {
                min_value = *iterators[i];
                count = 1;
            } else if (*iterators[i] == min_value) {
                count++;
            }
        }
        if (count == 0) {
            break;
        }
        if (count < kmers_by_sequence.size()) {
            // The kmer is informative.
            subchain.kmers.push_back(min_value);
            kmer_counts.push_back(count);
            if (count > max_count) {
                max_count = count;
            }
        }
        for (size_t i = 0; i < kmers_by_sequence.size(); i++) {
            if (iterators[i] != kmers_by_sequence[i].end() && *iterators[i] == min_value) {
                ++iterators[i];
            }
        }
    }
    subchain.kmer_counts = sdsl::int_vector<0>(kmer_counts.size(), 0, sdsl::bits::length(max_count));
    for (size_t i = 0; i < kmer_counts.size(); i++) {
        subchain.kmer_counts[i] = kmer_counts[i];
    }

    // Second pass: Build the kmer presence matrix and record the new order of the sequences.
    sdsl::bit_vector buffer(kmers_by_sequence.size() * subchain.kmers.size(), 0);
    for (size_t sequence_id = 0; sequence_id < kmers_by_sequence.size(); sequence_id++) {
        auto iter = kmers_by_sequence[sequence_id].begin();
        size_t start = sequence_id * subchain.kmers.size();
        for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
            while (iter != kmers_by_sequence[sequence_id].end() && *iter < subchain.kmers[kmer_id]) {
                ++iter;
            }
            if (iter != kmers_by_sequence[sequence_id].end() && *iter == subchain.kmers[kmer_id]) {
                buffer[start + kmer_id] = 1;
                ++iter;
            }
        }
    }
    subchain.kmers_present = KmerPresenceMatrix(
        sequences, subchain.sequences,
        subchain.kmers.size(), buffer
    );
}

void HaplotypePartitioner::build_subchains(const gbwtgraph::TopLevelChain& chain, Haplotypes::TopLevelChain& output, const Parameters& parameters) const {

    // Determine subchains and calculate some statistics.
    std::vector<Subchain> subchains = this->get_subchains(chain, parameters);
    if (this->verbosity >= Haplotypes::verbosity_debug) {
        size_t long_subchains = 0, with_extra_snarls = 0, extra_snarls = 0;
        for (const Subchain& subchain : subchains) {
            if (subchain.length > parameters.subchain_length) {
                long_subchains++;
            }
            if (subchain.extra_snarls > 0) {
                with_extra_snarls++;
                extra_snarls += subchain.extra_snarls;
            }
        }
        // The second condition should be redundant.
        if (long_subchains > 0 || with_extra_snarls > 0) {
            #pragma omp critical
            {
                std::cerr << "Chain " << chain.offset << " (" << output.contig_name << "): " << long_subchains << " long subchains ("
                    << with_extra_snarls << " with " << extra_snarls << " additional snarls)" << std::endl;
            }
        }
    }

    // Determine the sequences crossing each subchain.
    std::vector<std::pair<Subchain, std::vector<sequence_type>>> to_process;
    for (const Subchain& subchain : subchains) {
        auto sequences = this->get_sequences(subchain);
        if (sequences.empty()) {
            // There are no haplotypes crossing the subchain, so we break it into
            // a suffix and a prefix.
            // NOTE: See the general get_sequences() for the asymmetry in handling prefixes and suffixes.
            to_process.push_back({
                { Haplotypes::Subchain::suffix, subchain.start, empty_gbwtgraph_handle(), 0, 0 },
                this->get_sequences(subchain.start)
            });
            to_process.push_back({
                { Haplotypes::Subchain::prefix, empty_gbwtgraph_handle(), subchain.end, 0, 0 },
                this->get_sequences(subchain.end)
            });
        } else {
            to_process.push_back({ subchain, std::move(sequences) });
        }
    }

    // Convert the subchains to actual subchains.
    size_t fragmented_subchains = 0;
    size_t fragmented_sequences = 0;
    size_t additional_fragments = 0;
    if (parameters.subchain_kmers) {
        this->build_kmer_subchains(to_process, output);
    } else {
        for (auto iter = to_process.begin(); iter != to_process.end(); ++iter) {
            // We cache GBWT records and kmer counts at subchain level.
            gbwt::CachedGBWT cached_gbwt(this->gbz.index);
            vg::hash_map<kmer_type, size_t> kmer_cache;
            bool fragmented = false;
            output.subchains.push_back({
                iter->first.type,
                gbwtgraph::GBWTGraph::handle_to_node(iter->first.start), gbwtgraph::GBWTGraph::handle_to_node(iter->first.end),
                0, // placeholder
                {}, sdsl::int_vector<0>(0, 0, 64), {}, KmerPresenceMatrix()
            });
            Haplotypes::Subchain& subchain = output.subchains.back();
            std::vector<std::vector<kmer_type>> kmers_by_sequence;
            kmers_by_sequence.reserve(iter->second.size());
            for (sequence_type sequence : iter->second) {
                size_t fragments = 0;
                kmers_by_sequence.emplace_back(this->unique_minimizers(
                    sequence, iter->first, fragments, cached_gbwt, kmer_cache
                ));
                if (fragments > 1) {
                    fragmented = true;
                    fragmented_sequences++;
                    additional_fragments += fragments - 1;
                }
            }
            present_kmers(kmers_by_sequence, iter->second, subchain);
            if (fragmented) {
                fragmented_subchains++;
            }
        }
    }

    if (fragmented_subchains > 0 && this->verbosity >= Haplotypes::verbosity_debug) {
        #pragma omp critical
        {
            std::cerr << "Chain " << chain.offset << " (" << output.contig_name << "): " << fragmented_subchains << " fragmented subchains ("
                << fragmented_sequences << " fragmented sequences with " << additional_fragments << " additional fragments)" << std::endl;
        }
    }

    // Take entire sequences if we could not generate any haplotypes.
    // Note that the kmer sets should be empty, as the sequences should
    // be identical.
    if (subchains.empty()) {
        // We cache GBWT records and kmer counts at subchain level.
        gbwt::CachedGBWT cached_gbwt(this->gbz.index);
        vg::hash_map<kmer_type, size_t> kmer_cache;
        output.subchains.push_back({
            Haplotypes::Subchain::full_haplotype,
            gbwt::ENDMARKER, gbwt::ENDMARKER,
            0, // placeholder
            {}, sdsl::int_vector<0>(0, 0, 64), {}, KmerPresenceMatrix()
        });
        Haplotypes::Subchain& subchain = output.subchains.back();
        gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(chain.handle);
        auto da = this->r_index.decompressDA(node);
        std::vector<sequence_type> sequences;
        sequences.reserve(da.size());
        std::vector<std::vector<kmer_type>> kmers_by_sequence;
        kmers_by_sequence.reserve(da.size());
        for (auto seq_id : da) {
            kmers_by_sequence.emplace_back(this->unique_minimizers(seq_id, cached_gbwt, kmer_cache));
            sequences.emplace_back(seq_id, 0); // No start node, so we don't need a offset in it either.
        }
        present_kmers(kmers_by_sequence, sequences, subchain);
    }
}

hash_map<Haplotypes::Subchain::kmer_type, std::uint32_t> oriented_kmer_positions(
    const std::vector<Haplotypes::Subchain::kmer_type>& kmers, size_t k
) {
    hash_map<Haplotypes::Subchain::kmer_type, std::uint32_t> result;
    result.reserve(2 * kmers.size());
    for (size_t i = 0; i < kmers.size(); i++) {
        auto position = static_cast<std::uint32_t>(i);
        result[kmers[i]] = position;
        result[gbwtgraph::Key64(kmers[i]).reverse_complement(k).get_key()] = position;
    }
    return result;
}

hash_map<std::uint32_t, std::uint32_t> count_oriented_kmers(
    const std::vector<std::string>& pieces, const hash_map<Haplotypes::Subchain::kmer_type, std::uint32_t>& positions, size_t k
) {
    hash_map<std::uint32_t, std::uint32_t> counts;
    for (const std::string& piece : pieces) {
        gbwtgraph::Key64 key;
        size_t valid_chars = 0;
        for (char c : piece) {
            key.forward(k, static_cast<unsigned char>(c), valid_chars);
            if (valid_chars < k) {
                continue;
            }
            if (auto iter = positions.find(key.get_key()); iter != positions.end()) {
                counts[iter->second]++;
            }
        }
    }
    return counts;
}

namespace {

// A tally of the copies of m kmers in n sequences of a subchain, up to
// `Haplotypes::Subchain::MAX_COPIES`. A kmer is informative unless every sequence has
// the same number of copies of it. The tally writes the informative kmers into the
// subchain, numbered in the order of the kmers.
class KmerCopyTally {
public:
    KmerCopyTally(size_t n, size_t m) :
        n(n), m(m), present(n * m, 0), multi(n), present_in(m, 0), min_copies(m, MAX_COPIES), max_copies(m, 0) {
    }

    // Record the occurrences of the kmers in a sequence.
    void add(size_t sequence, const hash_map<std::uint32_t, std::uint32_t>& occurrences) {
        for (const auto& [kmer, count] : occurrences) {
            auto copies = static_cast<std::uint8_t>(std::min(size_t(count), MAX_COPIES));
            this->present[sequence * this->m + kmer] = true;
            this->present_in[kmer]++;
            this->min_copies[kmer] = std::min(this->min_copies[kmer], copies);
            this->max_copies[kmer] = std::max(this->max_copies[kmer], copies);
            if (copies > 1) {
                this->multi[sequence].emplace_back(kmer, copies);
            }
        }
        std::ranges::sort(this->multi[sequence]);
    }

    // Write the informative kmers, their presence in the sequences and their copies
    // into the subchain. The tally can write only once. The kmers are the m kmers of the tally, and the sequences are
    // the n sequences of the tally, which the subchain stores in its own order.
    void write(
        const std::vector<Haplotypes::Subchain::kmer_type>& kmers,
        const std::vector<HaplotypePartitioner::sequence_type>& sequences, Haplotypes::Subchain& output
    ) {
        this->number_informative();
        output.kmers.clear();
        output.kmers.reserve(this->informative);
        output.kmer_counts = sdsl::int_vector<0>(this->informative, 0, static_cast<std::uint8_t>(sdsl::bits::length(this->max_count)));
        for (size_t j = 0; j < this->m; j++) {
            if (this->final_id[j] < this->informative) {
                output.kmer_counts[output.kmers.size()] = this->present_in[j];
                output.kmers.push_back(kmers[j]);
            }
        }
        sdsl::bit_vector matrix(this->n * this->informative, 0);
        for (size_t i = 0; i < this->n; i++) {
            for (size_t j = 0; j < this->m; j++) {
                if (this->present[i * this->m + j] && this->final_id[j] < this->informative) {
                    matrix[i * this->informative + this->final_id[j]] = true;
                }
            }
        }
        this->present = sdsl::bit_vector();
        output.sequences.clear();
        output.kmers_present = KmerPresenceMatrix(sequences, output.sequences, this->informative, matrix);
        this->write_copies(sequences, output);
    }

private:
    constexpr static size_t MAX_COPIES = Haplotypes::Subchain::MAX_COPIES;

    size_t n;
    size_t m;

    // Presence of kmer j in sequence i at i * m + j.
    sdsl::bit_vector present;

    // For each sequence, (kmer, copies) for the kmers with more than one copy, sorted.
    std::vector<std::vector<std::pair<std::uint32_t, std::uint8_t>>> multi;

    // For each kmer: the number of sequences with it, and its fewest and most copies.
    std::vector<std::uint32_t> present_in;
    std::vector<std::uint8_t> min_copies;
    std::vector<std::uint8_t> max_copies;

    // The number of informative kmers, the id of each kmer among them (or the largest
    // value if it is not informative), and the largest number of sequences with one.
    size_t informative = 0;
    std::vector<std::uint32_t> final_id;
    size_t max_count = 0;

    void number_informative() {
        this->final_id.assign(this->m, std::numeric_limits<std::uint32_t>::max());
        for (size_t j = 0; j < this->m; j++) {
            if (this->present_in[j] < this->n || this->min_copies[j] != this->max_copies[j]) {
                this->final_id[j] = static_cast<std::uint32_t>(this->informative++);
                this->max_count = std::max(this->max_count, size_t(this->present_in[j]));
            }
        }
    }

    // Store the copies of the informative kmers in the order of the sequences in the subchain.
    void write_copies(const std::vector<HaplotypePartitioner::sequence_type>& sequences, Haplotypes::Subchain& output) const {
        std::map<std::pair<std::uint32_t, std::uint32_t>, size_t> original;
        for (size_t i = 0; i < this->n; i++) {
            original[{ static_cast<std::uint32_t>(sequences[i].first), static_cast<std::uint32_t>(sequences[i].second) }] = i;
        }
        std::vector<Haplotypes::Subchain::kmer_copies_type> lists(this->n);
        for (size_t r = 0; r < this->n; r++) {
            auto iter = original.find({ output.sequences[r].first, output.sequences[r].second });
            if (iter == original.end()) {
                throw std::runtime_error("HaplotypePartitioner::build_subchain_kmers(): cannot find a reordered sequence");
            }
            for (const auto& [kmer, copies] : this->multi[iter->second]) {
                if (this->final_id[kmer] < this->informative) {
                    lists[r].emplace_back(this->final_id[kmer], copies);
                }
            }
        }
        output.set_copies(lists);
    }
};

} // anonymous namespace

std::vector<Haplotypes::Subchain::kmer_type> HaplotypePartitioner::kmers_on_nodes(
    const hash_map<kmer_type, minimizer_index_type::minimizer_type>& minimizers, const hash_set<nid_t>& nodes
) const {
    std::vector<kmer_type> result;
    for (const auto& [key, minimizer] : minimizers) {
        auto occurrences = this->minimizer_index.find(minimizer);
        bool own = (occurrences.second > 0);
        for (size_t i = 0; own && i < occurrences.second; i++) {
            own = (nodes.find(this->minimizer_index.get_value(occurrences, i).first.id()) != nodes.end());
        }
        if (own) {
            result.push_back(key);
        }
    }
    std::ranges::sort(result);
    return result;
}

void HaplotypePartitioner::build_kmer_subchains(
    const std::vector<std::pair<Subchain, std::vector<sequence_type>>>& to_process, Haplotypes::TopLevelChain& output
) const {
    // Subchains are processed in parallel, unless the chains already are.
    size_t first = output.subchains.size();
    output.subchains.resize(first + to_process.size());
    #pragma omp parallel for schedule(dynamic, 1) if(!omp_in_parallel())
    for (size_t i = 0; i < to_process.size(); i++) {
        Haplotypes::Subchain& subchain = output.subchains[first + i];
        subchain.type = to_process[i].first.type;
        subchain.start = gbwtgraph::GBWTGraph::handle_to_node(to_process[i].first.start);
        subchain.end = gbwtgraph::GBWTGraph::handle_to_node(to_process[i].first.end);
        subchain.placeholder = 0;
        try {
            this->build_subchain_kmers(to_process[i].first, to_process[i].second, subchain);
        } catch (const std::runtime_error& e) {
            std::cerr << "error: [subchain " << (first + i) << "]: " << e.what() << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
}

void HaplotypePartitioner::build_subchain_kmers(
    const Subchain& subchain, const std::vector<sequence_type>& sequences, Haplotypes::Subchain& output
) const {
    gbwt::CachedGBWT cached_gbwt(this->gbz.index);
    size_t k = this->minimizer_index.k();
    gbwt::node_type start = gbwtgraph::GBWTGraph::handle_to_node(subchain.start);
    gbwt::node_type end = gbwtgraph::GBWTGraph::handle_to_node(subchain.end);

    // The nodes the sequences visit and the distinct minimizers in the sequences.
    hash_set<nid_t> nodes;
    hash_map<kmer_type, minimizer_index_type::minimizer_type> distinct;
    for (const sequence_type& sequence : sequences) {
        for (nid_t node : subchain_walk_nodes(cached_gbwt, start, end, sequence, this->fragment_map)) {
            nodes.insert(node);
        }
        for (const std::string& piece : subchain_haplotype(sequence, subchain, *this, cached_gbwt)) {
            for (const auto& minimizer : this->minimizer_index.minimizers(piece)) {
                distinct.emplace(minimizer.key.get_key(), minimizer);
            }
        }
    }
    std::vector<kmer_type> kept = this->kmers_on_nodes(distinct, nodes);
    distinct = hash_map<kmer_type, minimizer_index_type::minimizer_type>();

    // Count the occurrences of the kept kmers in each sequence, in both orientations.
    hash_map<kmer_type, std::uint32_t> positions = oriented_kmer_positions(kept, k);
    KmerCopyTally tally(sequences.size(), kept.size());
    for (size_t i = 0; i < sequences.size(); i++) {
        tally.add(i, count_oriented_kmers(subchain_haplotype(sequences[i], subchain, *this, cached_gbwt), positions, k));
    }
    tally.write(kept, sequences, output);
}

//------------------------------------------------------------------------------

/*
 * A haplotype beging generated as a GBWT path.
 *
 * GBWT metadata will be set as following:
 *
 * * Sample name is "recombination".
 * * Contig name is taken from the top-level chain.
 * * Haplotype identifier is set during construction.
 * * Fragment identifier is set as necessary.
 */
struct RecombinatorHaplotype {
    typedef Recombinator::sequence_type sequence_type;

    // Recombinator instance we are using.
    const Recombinator& recombinator;

    // GBWT builder and metadata builder used by this job.
    gbwt::GBWTBuilder& builder;
    gbwtgraph::MetadataBuilder& metadata;

    // Contig name in GBWT metadata.
    const std::string& contig_name;

    // Whether this haplotype should be wrapped, appending its origin (first)
    // fragment onto its last fragment to create an origin-spanning edge (for
    // circular contigs).
    bool wrap;

    // Nodes of the origin (first) fragment, saved when wrapping so that they
    // can be appended onto the last fragment. Empty until the origin fragment
    // has been inserted.
    gbwt::vector_type origin_fragment;

    // Haplotype identifier in GBWT metadata.
    size_t id;

    // Fragment identifier in GBWT metadata.
    // If no original haplotype crosses a subchain, a new fragment will
    // start after the subchain.
    size_t fragment;

    // GBWT sequence indentifier in the previous subchain, or
    // `gbwt::invalid_sequence()` if there was no such sequence.
    gbwt::size_type sequence_id;

    // GBWT position at the end of the latest `extend()` call.
    // `gbwt::invalid_edge()` otherwise.
    gbwt::edge_type position;

    // The path being generated.
    gbwt::vector_type path;

    // Constructor that starts a new haplotype with the given contig name and identifier.
    RecombinatorHaplotype(
        const Recombinator& recombinator,
        gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata,
        const std::string& contig_name, bool wrap, size_t id
    ) :
        recombinator(recombinator),
        builder(builder), metadata(metadata),
        contig_name(contig_name), wrap(wrap), id(id), fragment(0),
        sequence_id(gbwt::invalid_sequence()), position(gbwt::invalid_edge()) {}

    /*
     * Extends the haplotype over the given subchain by using the given
     * original haplotype. If the haplotype is fragmented, this may
     * finish the current fragment and start a new one.
     *
     * This assumes that the original haplotype crosses the subchain.
     *
     * If `extend()` has been called for this fragment, there must be a
     * unary path connecting the subchains, which will be used in the
     * generated haplotype.
     *
     * If `extend()` has not been called, the generated haplotype will
     * take the prefix of the original original haplotype until the start
     * of the subchain.
     *
     * `is_last_subchain` indicates whether this is the last subchain of the
     * chain. It is needed because a suffix subchain finishes a fragment
     * immediately, and that fragment is the last one only when the suffix is
     * the last subchain (a broken snarl can put a suffix in the middle).
     */
    void extend(sequence_type sequence, const Haplotypes::Subchain& subchain, bool is_last_subchain);

    // Takes an existing haplotype from the GBWT index and inserts it into
    // the builder. This is intended for fragments that do not contain
    // subchains crossed by the original haplotypes. The call will fail if
    // `extend()` has been called.
    void take(gbwt::size_type sequence_id);

    // Finishes the current haplotype fragment and starts a new one.
    // If `with_suffix` is true, the current fragment will be extended
    // until the end of the corresponding path. In that case, the call will
    // fail if `extend()` has not been called for this fragment.
    // If `is_final` is true, this is the last fragment of the haplotype, and
    // the origin fragment will be appended onto it when wrapping.
    void finish(bool with_suffix, bool is_final = true);

private:
    // Extends the haplotype over a unary path from a previous subchain.
    void connect(gbwt::node_type until);

    // Takes a prefix of the current sequence until `until`, assuming that
    // `sequence_id` has been set.
    void prefix(gbwt::node_type until);

    // Extends the haplotype from the previous subchain until the end.
    void suffix();

    // Inserts the current fragment into the builder. If `is_final` is true and
    // wrapping is enabled, the origin fragment is appended onto it first.
    void insert(bool is_final);
};

void RecombinatorHaplotype::extend(sequence_type sequence, const Haplotypes::Subchain& subchain, bool is_last_subchain) {
    if (subchain.type == Haplotypes::Subchain::full_haplotype) {
        throw std::runtime_error("Haplotype::extend(): cannot extend a full haplotype");
    }

    if (sequence.first == gbwt::invalid_sequence()) {
        throw std::runtime_error("Haplotype::extend(): invalid sequence id");
    }
    this->sequence_id = sequence.first;

    if (subchain.type == Haplotypes::Subchain::prefix) {
        if (!this->path.empty())  {
            throw std::runtime_error("Haplotype::extend(): got a prefix subchain after the start of a fragment");
        }
        this->prefix(subchain.end);
        return;
    }

    // Suffixes and normal subchains have a start node, so we must reach it first.
    if (!this->path.empty()) {
        this->connect(subchain.start);
    } else {
        this->prefix(subchain.start);
    }

    const gbwt::GBWT& index = this->recombinator.gbz.index;
    this->position = gbwt::edge_type(subchain.start, sequence.second);
    if (!index.contains(this->position)) {
        throw std::runtime_error("Haplotype::extend(): the GBWT index does not contain position (" + std::to_string(this->position.first) + ", " + std::to_string(this->position.second) + ")");
    }

    if (subchain.type == Haplotypes::Subchain::suffix) {
        this->finish(true, is_last_subchain);
        return;
    }

    // This is a normal subchain.
    while (this->position.first != subchain.end) {
        this->position = index.LF(this->position);
        if (this->position.first == gbwt::ENDMARKER) {
            gbwt::size_type prev = this->sequence_id;
            this->finish(false, false);
            do {
                // Find the next non-empty fragment.
                this->sequence_id = this->recombinator.fragment_map.oriented_next(prev);
                if (this->sequence_id == gbwt::invalid_sequence()) {
                    std::string msg = "Haplotype::extend(): no successor for GBWT sequence " + std::to_string(prev);
                    throw std::runtime_error(msg);
                }
                this->position = index.start(this->sequence_id);
                prev = this->sequence_id;
            } while (this->position.first == gbwt::ENDMARKER);
        }
        this->path.push_back(this->position.first);
    }
}

void RecombinatorHaplotype::take(gbwt::size_type sequence_id) {
    const gbwt::GBWT& index = this->recombinator.gbz.index;
    if (!this->path.empty()) {
        throw std::runtime_error("Haplotype::take(): the current fragment is not empty");
    }
    if (sequence_id >= index.sequences()) {
        throw std::runtime_error("Haplotype::take(): the GBWT index does not contain sequence " + std::to_string(sequence_id));
    }
    this->path = index.extract(sequence_id);
    this->finish(false);
}

void RecombinatorHaplotype::finish(bool with_suffix, bool is_final) {
    if (with_suffix) {
        if (this->position == gbwt::invalid_edge()) {
            throw std::runtime_error("Haplotype::finish(): there is no current position");
        }
        this->suffix();
    }
    this->insert(is_final);
    this->fragment++;
    this->sequence_id = gbwt::invalid_sequence();
    this->position = gbwt::invalid_edge();
    this->path.clear();
}

void RecombinatorHaplotype::connect(gbwt::node_type until) {
    const gbwtgraph::GBWTGraph& graph = this->recombinator.gbz.graph;
    handle_t curr = gbwtgraph::GBWTGraph::node_to_handle(this->position.first);
    handle_t end = gbwtgraph::GBWTGraph::node_to_handle(until);
    this->position = gbwt::invalid_edge();
    hash_set<handle_t> visited;
    while (curr != end) {
        if (visited.find(curr) != visited.end()) {
            throw std::runtime_error("Haplotype::connect(): the path contains a cycle");
        }
        visited.insert(curr);
        handle_t successor = empty_gbwtgraph_handle();
        size_t successors = 0;
        graph.follow_edges(curr, false, [&](const handle_t& next) {
            successor = next;
            successors++;
        });
        if (successors != 1) {
            throw std::runtime_error("Haplotype::connect(): the path is not unary");
        }
        this->path.push_back(gbwtgraph::GBWTGraph::handle_to_node(successor));
        curr = successor;
    }
}

void RecombinatorHaplotype::prefix(gbwt::node_type until) {
    const gbwt::GBWT& index = this->recombinator.gbz.index;
    this->position = gbwt::invalid_edge();
    if (this->sequence_id >= index.sequences()) {
        throw std::runtime_error("Haplotype::prefix(): invalid GBWT sequence id " + std::to_string(sequence_id));
    }
    for (this->position = index.start(this->sequence_id); this->position.first != gbwt::ENDMARKER; this->position = index.LF(this->position)) {
        this->path.push_back(this->position.first);
        if (this->position.first == until) {
            return;
        }
    }
    throw std::runtime_error("Haplotype::prefix(): GBWT sequence " + std::to_string(this->sequence_id) + " did not reach GBWT node " + std::to_string(until));
}

void RecombinatorHaplotype::suffix() {
    const gbwt::GBWT& index = this->recombinator.gbz.index;
    for (gbwt::edge_type curr = index.LF(this->position); curr.first != gbwt::ENDMARKER; curr = index.LF(curr)) {
        this->path.push_back(curr.first);
    }
}

void RecombinatorHaplotype::insert(bool is_final) {
    std::string sample_name = "recombination"; // TODO: Make this a static class variable.
    this->metadata.add_haplotype(sample_name, this->contig_name, this->id, this->fragment);
    if (this->wrap) {
        // Linearizing a circular sequence drops the adjacency joining its end
        // back to its start. We restore it by appending the origin (first)
        // fragment onto the last fragment. We save the origin fragment here,
        // before it is inserted, so that a distinct copy remains available for
        // an unfragmented haplotype where the origin is also the last fragment.
        if (this->fragment == 0) {
            this->origin_fragment = this->path;
        }
        if (is_final) {
            append_wrap_fragment(this->path, this->origin_fragment);
        }
    }
    this->builder.insert(this->path, true);
}

//------------------------------------------------------------------------------

/*
 * An additional haplotype fragment from a bad subchain.
 *
 * GBWT metadata will be set as following:
 *
 * * Sample name is "fragment".
 * * Contig name is "contig_N", with the name taken from the top-level chain and the number from the subchain.
 * * Haplotype identifier is set during construction.
 * * Fragment identifier is used if the fragment itself is fragmented.
 */
struct RecombinatorFragment {
    // Sequence identifier for the start of the fragment.
    gbwt::size_type sequence_id;

    // GBWT starting position (inclusive).
    gbwt::edge_type from;

    // GBWT ending position (inclusive).
    gbwt::node_type to;

    // Identifier of the subchain.
    size_t subchain_id;

    // Haplotype identifier.
    size_t haplotype_id;

    RecombinatorFragment(gbwt::size_type sequence_id, gbwt::edge_type from, gbwt::node_type to, size_t subchain_id, size_t haplotype_id) :
        sequence_id(sequence_id), from(from), to(to), subchain_id(subchain_id), haplotype_id(haplotype_id) {}

    // Generates the GBWT path(s) and the metadata for the fragment.
    // Returns the number of paths generated.
    size_t generate(
        const gbwt::GBWT& index, const gbwt::FragmentMap& fragment_map,
        gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata,
        const std::string& contig_name
    ) const;
};

size_t RecombinatorFragment::generate(
    const gbwt::GBWT& index, const gbwt::FragmentMap& fragment_map,
    gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata,
    const std::string& contig_name
) const {
    gbwt::vector_type path;
    std::string sample = "fragment"; // TODO: Make this a static class variable.
    std::string contig = contig_name + "_" + std::to_string(this->subchain_id);
    size_t fragment = 0;
    auto finish_fragment = [&]() {
        if (!path.empty()) {
            metadata.add_haplotype(sample, contig, this->haplotype_id, fragment);
            builder.insert(path, true);
            fragment++;
            path.clear();
        }
    };

    size_t sequence = this->sequence_id;
    gbwt::edge_type pos = this->from;
    while (true) {
        if (pos.first == gbwt::ENDMARKER) {
            if (this->to == gbwt::ENDMARKER) {
                // We are in a suffix subchain, so we are not interested in subsequent fragments.
                break;
            }
            do {
                finish_fragment();
                sequence = fragment_map.oriented_next(sequence);
                if (sequence == gbwt::invalid_sequence()) {
                    std::string msg = "RecombinatorFragment::generate(): no successor for GBWT sequence " + std::to_string(this->sequence_id);
                    throw std::runtime_error(msg);
                }
                pos = index.start(sequence);
            } while (pos.first == gbwt::ENDMARKER);
        }
        path.push_back(pos.first);
        if (pos.first == this->to) {
            break;
        }
        pos = index.LF(pos);
    }
    finish_fragment();

    return fragment;
}

//------------------------------------------------------------------------------

void Recombinator::Statistics::combine(const Statistics& another) {
    this->chains += another.chains;
    this->subchains += another.subchains;
    this->bad_subchains += another.bad_subchains;
    this->fragments += another.fragments;
    this->full_haplotypes += another.full_haplotypes;
    this->haplotypes = std::max(this->haplotypes, another.haplotypes);
    this->extra_fragments += another.extra_fragments;
    this->connections += another.connections;
    this->ref_paths += another.ref_paths;
    this->copied_paths += another.copied_paths;
    this->kmers += another.kmers;
    this->score += another.score;
}

std::ostream& Recombinator::Statistics::print(std::ostream& out) const {
    out << this->haplotypes << " haplotypes for " << this->chains << " chains ("
        << this->full_haplotypes << " full, " << this->subchains << " subchains, " << this->fragments << " fragments)";
    if (this->bad_subchains > 0) {
        out << "; " << this->bad_subchains << " bad subchains with " << this->extra_fragments << " extra fragments";
    }
    if (this->subchains > 0) {
        double connection_rate = static_cast<double>(this->connections) / (this->subchains * this->haplotypes);
        out << "; connection rate " << connection_rate;
    }
    if (this->ref_paths > 0) {
        out << "; included " << this->ref_paths << " reference paths";
    }
    if (this->copied_paths > 0) {
        out << "; copied " << this->copied_paths << " paths from excluded chains";
    }
    if (this->kmers > 0) {
        double average_score = this->score / (this->kmers * this->haplotypes);
        out << "; used " << this->kmers << " kmers with average score " << average_score;
    }
    return out;
}

//------------------------------------------------------------------------------

Recombinator::Recombinator(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, Verbosity verbosity) :
    gbz(gbz), haplotypes(haplotypes), fragment_map(gbz.index.metadata, verbosity >= Haplotypes::verbosity_extra_debug), verbosity(verbosity)
{
    if (this->verbosity >= Haplotypes::verbosity_detailed) {
        std::cerr << "Recombinator: " << this->gbz.index.metadata.paths() << " fragments for " << this->fragment_map.size() << " haplotype sequences" << std::endl;
    }

    this->jobs_for_cached_paths = this->haplotypes.assign_reference_paths(this->gbz, this->verbosity);
}

//------------------------------------------------------------------------------

std::string Recombinator::Parameters::preset_list(std::string_view separator) {
    std::string result = "{";
    for (const auto& preset : PRESETS) {
        if (result.size() > 1) {
            result += separator;
        }
        result += preset.second;
    }
    result += "}";
    return result;
}

std::optional<Recombinator::Parameters::preset_t> Recombinator::Parameters::parse_preset(std::string_view name) {
    for (const auto& preset : PRESETS) {
        if (preset.second == name) {
            return preset.first;
        }
    }
    return std::nullopt;
}

Recombinator::Parameters::Parameters(preset_t preset) {
    if (preset == preset_haploid) {
        this->haploid_scoring = true;
        this->include_reference = true;
    } else if (preset == preset_diploid) {
        this->num_haplotypes = NUM_CANDIDATES;
        this->diploid_sampling = true;
        this->include_reference = true;
    } else if (preset == preset_hifi) {
        this->num_haplotypes = NUM_CANDIDATES;
        this->diploid_sampling = true;
        this->include_reference = true;
        this->ignore_homopolymers = LONG_READ_HOMOPOLYMERS;
    } else if (preset == preset_ont) {
        this->num_haplotypes = NUM_CANDIDATES;
        this->diploid_sampling = true;
        this->include_reference = true;
        this->ignore_homopolymers = LONG_READ_HOMOPOLYMERS;
    } else if (preset == preset_model) {
        this->num_haplotypes = MODEL_NUM_HAPLOTYPES;
        this->include_reference = true;
        this->continuity_bonus = MODEL_CONTINUITY_BONUS;
        this->max_new_bp = MODEL_MAX_NEW_BP;
        this->bp_escape = MODEL_BP_ESCAPE;
        this->max_n_bp = MODEL_MAX_N_BP;
        this->similar_threading = true;
        this->model_pair = true;
    }
}

void Recombinator::Parameters::print(std::ostream& out) const {
    out << "Sampling parameters:" << std::endl;
    if (this->haploid_scoring) {
        out << "- haploid scoring (absent " << this->absent_score << ", present " << this->present_discount << ")" << std::endl;
    } else {
        out << "- diploid scoring (absent " << this->absent_score << ", het " << this->het_adjustment << ", present " << this->present_discount << ")" << std::endl;
    }
    if (this->coverage > 0) {
        out << "- kmer coverage " << this->coverage << std::endl;
    }
    if (this->diploid_sampling) {
        out << "- diploid sampling (" << this->num_haplotypes << " candidates";
        if (this->extra_fragments) {
            out << ", badness threshold " << this->badness_threshold;
        }
        out << ")" << std::endl;
    } else {
        out << "- heuristic sampling (" << this->num_haplotypes << " haplotypes)" << std::endl;
    }
    if (this->include_reference) {
        out << "- include reference paths (gref fragments clipped, min length " << this->min_gref_length << ")" << std::endl;
    }
    if (this->ignore_homopolymers != std::numeric_limits<size_t>::max()) {
        out << "- ignore homopolymers of length >= " << this->ignore_homopolymers << std::endl;
    }
    if (this->continuity_bonus > 0.0) {
        out << "- continuity bonus " << this->continuity_bonus << std::endl;
    }
    if (this->has_new_bp_limit()) {
        out << "- at most " << this->max_new_bp << " bp of new node sequence";
        if (this->bp_escape > 0.0) {
            out << ", unless a new sample kmer per " << this->bp_escape << " bp";
        }
        out << std::endl;
    }
    if (this->has_n_base_limit()) {
        out << "- at most " << this->max_n_bp << " N bases" << std::endl;
    }
    if (this->similar_threading) {
        out << "- similar threading" << std::endl;
    }
    if (this->model_pair) {
        out << "- model pair" << std::endl;
    }
    if (!this->high_coverage_chains.empty()) {
        out << "- " << this->high_coverage_chains.size() << " high-coverage chains ("
            << this->high_coverage_num_haplotypes << " haplotypes each)" << std::endl;
    }
    if (!this->half_coverage_chains.empty()) {
        out << "- " << this->half_coverage_chains.size() << " half-coverage chains ("
            << this->half_coverage_num_haplotypes << " haplotypes each)" << std::endl;
    }
    if (!this->excluded_chains.empty()) {
        out << "- " << this->excluded_chains.size() << " chains excluded from personalization" << std::endl;
    }
}

//------------------------------------------------------------------------------

void add_path(const gbwt::GBWT& source, gbwt::size_type path_id, gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata) {
    // We know that sufficient metadata exists, because this is a cached path.
    // This keeps the subrange of a generic path. Reference samples will be copied later.
    metadata.add_gbwt_path(source.metadata.fullPath(path_id));
    builder.insert(source.extract(gbwt::Path::encode(path_id, false)), true);
}

//------------------------------------------------------------------------------

std::vector<PathPiece> clip_path_to_index(
    const gbwt::vector_type& path, const gbwt::DynamicGBWT& index,
    const std::function<size_t(gbwt::node_type)>& node_length
) {
    std::vector<PathPiece> result;
    size_t offset = 0;
    for (size_t i = 0; i < path.size(); i++) {
        gbwt::node_type node = path[i];
        size_t length = node_length(node);
        if (!result.empty() && result.back().start + result.back().nodes == i && index.hasEdge(path[i - 1], node)) {
            result.back().nodes++;
            result.back().bp_length += length;
        } else if (index.contains(node) && !index.empty(node)) {
            // As in GBWTGraph, a node within the alphabet is absent if its record is empty.
            result.push_back({ i, 1, offset, length });
        }
        offset += length;
    }
    return result;
}

// Adds the pieces of the given gref fragments that only use nodes and edges already in the
// builder, named as subranges of the fragment in the full graph, and returns their number.
// Call after `builder.finish()`. Everything is clipped before anything is inserted, because
// an insert may flush the builder and start a thread that modifies the index.
size_t add_gref_fragments(
    const gbwtgraph::GBZ& gbz, const std::vector<gbwt::size_type>& path_ids, size_t min_length,
    gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata
) {
    auto node_length = [&](gbwt::node_type node) -> size_t {
        return gbz.graph.get_length(gbwtgraph::GBWTGraph::node_to_handle(node));
    };
    std::vector<gbwt::vector_type> pieces;
    for (gbwt::size_type path_id : path_ids) {
        gbwt::vector_type path = gbz.index.extract(gbwt::Path::encode(path_id, false));
        for (const PathPiece& piece : clip_path_to_index(path, builder.index, node_length)) {
            if (piece.bp_length < min_length) {
                continue;
            }
            // The metadata must be added in the same order as the paths.
            gbwt::FullPathName name = gbz.index.metadata.fullPath(path_id);
            name.offset += piece.bp_offset;
            metadata.add_gbwt_path(name);
            pieces.emplace_back(path.begin() + piece.start, path.begin() + piece.start + piece.nodes);
        }
    }
    for (const gbwt::vector_type& piece : pieces) {
        builder.insert(piece, true);
    }
    return pieces.size();
}

//------------------------------------------------------------------------------

void validate_recombinator_parameters(const Recombinator::Parameters& parameters) {
    if (parameters.num_haplotypes == 0) {
        std::string msg = "validate_recombinator_parameters(): number of haplotypes cannot be 0";
        throw std::runtime_error(msg);
    }
    if (parameters.badness_threshold <= 0.0) {
        std::string msg = "validate_recombinator_parameters(): badness threshold must be positive";
        throw std::runtime_error(msg);
    }
    if (parameters.diploid_sampling && parameters.num_haplotypes < 2) {
        std::string msg = "validate_recombinator_parameters(): diploid sampling requires at least 2 haplotypes";
        throw std::runtime_error(msg);
    }
    if (parameters.present_discount < 0.0 || parameters.present_discount > 1.0) {
        std::string msg = "validate_recombinator_parameters(): present discount must be between 0.0 and 1.0";
        throw std::runtime_error(msg);
    }
    if (parameters.het_adjustment < 0.0) {
        std::string msg = "validate_recombinator_parameters(): het adjustment must be non-negative";
        throw std::runtime_error(msg);
    }
    if (parameters.absent_score < 0.0) {
        std::string msg = "validate_recombinator_parameters(): absent score must be non-negative";
        throw std::runtime_error(msg);
    }
    if (parameters.continuity_bonus < 0.0) {
        std::string msg = "validate_recombinator_parameters(): continuity bonus must be non-negative";
        throw std::runtime_error(msg);
    }
    if (parameters.bp_escape < 0.0) {
        std::string msg = "validate_recombinator_parameters(): bp escape must be non-negative";
        throw std::runtime_error(msg);
    }
    if (parameters.bp_escape > 0.0 && !parameters.has_new_bp_limit()) {
        std::string msg = "validate_recombinator_parameters(): bp escape needs a new-bp limit";
        throw std::runtime_error(msg);
    }
    if (parameters.model_pair && parameters.num_haplotypes < 2) {
        std::string msg = "validate_recombinator_parameters(): the model pair requires at least 2 haplotypes";
        throw std::runtime_error(msg);
    }
    if (parameters.model_pair && parameters.diploid_sampling) {
        std::string msg = "validate_recombinator_parameters(): the model pair cannot be combined with diploid sampling";
        throw std::runtime_error(msg);
    }
    if (parameters.model_pair && parameters.scoring_model != Recombinator::Parameters::standard_scoring) {
        std::string msg = "validate_recombinator_parameters(): the model pair needs the standard scoring model";
        throw std::runtime_error(msg);
    }
    if (parameters.model_pair && parameters.haploid_scoring) {
        std::string msg = "validate_recombinator_parameters(): the model pair cannot be combined with haploid scoring";
        throw std::runtime_error(msg);
    }
    if (!parameters.high_coverage_chains.empty() && parameters.high_coverage_num_haplotypes == 0) {
        std::string msg = "validate_recombinator_parameters(): number of high-coverage haplotypes cannot be 0";
        throw std::runtime_error(msg);
    }
    if (!parameters.half_coverage_chains.empty() && parameters.half_coverage_num_haplotypes == 0) {
        std::string msg = "validate_recombinator_parameters(): number of half-coverage haplotypes cannot be 0";
        throw std::runtime_error(msg);
    }
    // A chain may be assigned to at most one of the excluded / high-coverage /
    // half-coverage sets.
    for (size_t chain_id : parameters.high_coverage_chains) {
        if (parameters.half_coverage_chains.find(chain_id) != parameters.half_coverage_chains.end()) {
            std::string msg = "validate_recombinator_parameters(): a chain cannot be both high-coverage and half-coverage";
            throw std::runtime_error(msg);
        }
        if (parameters.excluded_chains.find(chain_id) != parameters.excluded_chains.end()) {
            std::string msg = "validate_recombinator_parameters(): a chain cannot be both high-coverage and excluded";
            throw std::runtime_error(msg);
        }
    }
    for (size_t chain_id : parameters.half_coverage_chains) {
        if (parameters.excluded_chains.find(chain_id) != parameters.excluded_chains.end()) {
            std::string msg = "validate_recombinator_parameters(): a chain cannot be both half-coverage and excluded";
            throw std::runtime_error(msg);
        }
    }
}

double get_or_estimate_coverage(
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& counts,
    const Recombinator::Parameters& parameters,
    Haplotypes::Verbosity verbosity) {
    if (parameters.coverage > 0) {
        return parameters.coverage;
    }

    double start = gbwt::readTimer();
    if (verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Estimating kmer coverage" << std::endl;
    }
    std::map<size_t, size_t> count_to_frequency;
    // FIXME: parameter to include singletons in the estimation.
    for (auto iter = counts.begin(); iter != counts.end(); ++iter) {
        // We are only interested in kmers with multiple occurrences, as unique
        // kmers are likely sequencing errors.
        if (iter->second > 1) {
            count_to_frequency[iter->second]++;
        }
    }

    // Use mode as the initial estimate for coverage.
    auto statistics = summary_statistics(count_to_frequency);
    if (verbosity >= Haplotypes::verbosity_detailed) {
        std::cerr << "Coverage: median " << statistics.median
            << ", mean " << statistics.mean
            << ", stdev " << statistics.stdev
            << ", mode " << statistics.mode;
    }

    // First attempt: The most common kmer count is the true coverage.
    double coverage = statistics.mode;
    bool reliable = (statistics.mode >= statistics.median);

    // Second attempt: If we do diploid scoring, the most common count may
    // be for heterozygous kmers. We therefore try to find a secondary peak
    // at ~2x mode and use it if it is good enough and above the median.
    if (!reliable && !parameters.haploid_scoring) {
        size_t low = 1.7 * statistics.mode, high = 2.3 * statistics.mode;
        size_t peak = count_to_frequency[coverage];
        size_t best = low, secondary = count_to_frequency[low];
        for (size_t i = low + 1; i <= high; i++) {
            if (count_to_frequency[i] > secondary) {
                best = i; secondary = count_to_frequency[i];
            }
        }
        if (verbosity >= Haplotypes::verbosity_detailed) {
            std::cerr << "; secondary peak at " << best;
        }
        if (best >= size_t(statistics.median) && secondary >= peak / 2) {
            coverage = best;
            reliable = true;
        }
    }

    // Third attempt: Use the median as an unreliable estimate.
    if (!reliable) {
        coverage = statistics.median;
    }

    if (verbosity >= Haplotypes::verbosity_detailed) {
        std::cerr << "; using " << coverage << std::endl;
    }
    if (!reliable) {
        std::cerr << "warning: Kmer coverage estimate is unreliable" << std::endl;
    }

    if (verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Estimated kmer coverage in " << seconds << " seconds" << std::endl;
    }
    return coverage;
}

gbwt::GBWT Recombinator::generate_haplotypes(const std::string& kff_file, const Parameters& parameters) const {

    // Input validation (may throw).
    validate_recombinator_parameters(parameters);

    if (this->verbosity >= Haplotypes::verbosity_detailed) {
        parameters.print(std::cerr);
    }

    // Get kmer counts (may throw) and determine coverage.
    hash_map<Haplotypes::Subchain::kmer_type, size_t> counts = this->haplotypes.kmer_counts(kff_file, this->verbosity);
    double coverage = get_or_estimate_coverage(counts, parameters, this->verbosity);

    double start = gbwt::readTimer();
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Building GBWT" << std::endl;
    }

    // Determine construction jobs.
    std::vector<std::vector<size_t>> jobs(this->haplotypes.jobs());
    for (auto& chain : this->haplotypes.chains) {
        if (chain.job_id < this->haplotypes.jobs()) {
            jobs[chain.job_id].push_back(chain.offset);
        }
    }

    // Paths crossing excluded chains are copied through by copy_chain(), so we
    // must not add them again as reference paths.
    std::unordered_set<gbwt::size_type> excluded_path_ids;
    for (size_t chain_id : parameters.excluded_chains) {
        for (const Haplotypes::Subchain& subchain : this->haplotypes.chains[chain_id].subchains) {
            for (const Haplotypes::compact_sequence_type& sequence : subchain.sequences) {
                excluded_path_ids.insert(gbwt::Path::id(sequence.first));
            }
        }
    }

    // Figure out GBWT path ids for reference paths in each job. Gref fragments are set
    // aside: copying them would keep the whole gref cover, so they are clipped later.
    std::vector<std::vector<gbwt::size_type>> reference_paths(this->haplotypes.jobs());
    std::vector<std::vector<gbwt::size_type>> gref_fragment_paths(this->haplotypes.jobs());
    if (parameters.include_reference) {
        for (size_t i = 0; i < this->gbz.graph.named_paths.size(); i++) {
            size_t job_id = this->jobs_for_cached_paths[i];
            gbwt::size_type path_id = this->gbz.graph.named_paths[i].id;
            if (job_id < this->haplotypes.jobs() && excluded_path_ids.find(path_id) == excluded_path_ids.end()) {
                std::string path_name = this->gbz.graph.get_path_name(handlegraph::as_path_handle(i));
                if (GrefCover::is_gref_fragment(path_name)) {
                    gref_fragment_paths[job_id].push_back(path_id);
                } else {
                    reference_paths[job_id].push_back(path_id);
                }
            }
        }
    }

    // Choose the model pairs.
    std::vector<std::vector<std::pair<size_t, size_t>>> pairs(this->haplotypes.chains.size());
    if (parameters.model_pair) {
        pairs = this->model_pairs(counts, coverage, parameters);
    }

    // Build partial indexes.
    // We use a separate MetadataBuilder for each job, because we don't know in advance
    // how many fragments there will be for each generated haplotype.
    double checkpoint = gbwt::readTimer();
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Running " << omp_get_max_threads() << " GBWT construction jobs in parallel" << std::endl;
    }
    std::vector<gbwt::GBWT> indexes(jobs.size());
    Statistics statistics;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < jobs.size(); job++) {
        gbwt::GBWTBuilder builder(sdsl::bits::length(this->gbz.index.sigma() - 1), parameters.buffer_size);
        gbwtgraph::MetadataBuilder metadata;
        Statistics job_statistics;
        // Add named and reference paths.
        for (auto path_id : reference_paths[job]) {
            add_path(this->gbz.index, path_id, builder, metadata);
            job_statistics.ref_paths++;
        }
        // Add haplotypes for each chain.
        for (auto chain_id : jobs[job]) {
            try {
                Statistics chain_statistics;
                if (parameters.excluded_chains.find(chain_id) != parameters.excluded_chains.end()) {
                    // Excluded chains are copied through verbatim instead of
                    // personalizing them. This takes precedence over the scoring
                    // models.
                    chain_statistics = this->copy_chain(
                        this->haplotypes.chains[chain_id], builder, metadata
                    );
                } else {
                    // High- and half-coverage chains use their own scoring model,
                    // no diploid sampling, and their own haplotype count.
                    const Parameters* chain_parameters = &parameters;
                    Parameters chain_parameters_override;
                    if (parameters.high_coverage_chains.find(chain_id) != parameters.high_coverage_chains.end()) {
                        chain_parameters_override = parameters;
                        chain_parameters_override.scoring_model = Parameters::high_coverage_scoring;
                        chain_parameters_override.diploid_sampling = false;
                        chain_parameters_override.num_haplotypes = parameters.high_coverage_num_haplotypes;
                        chain_parameters = &chain_parameters_override;
                    } else if (parameters.half_coverage_chains.find(chain_id) != parameters.half_coverage_chains.end()) {
                        chain_parameters_override = parameters;
                        chain_parameters_override.scoring_model = Parameters::half_coverage_scoring;
                        chain_parameters_override.diploid_sampling = false;
                        chain_parameters_override.num_haplotypes = parameters.half_coverage_num_haplotypes;
                        chain_parameters = &chain_parameters_override;
                    }
                    chain_statistics = this->generate_haplotypes(
                        this->haplotypes.chains[chain_id], counts, builder, metadata, *chain_parameters, coverage, pairs[chain_id]
                    );
                }
                job_statistics.combine(chain_statistics);
            } catch (const std::runtime_error& e) {
                std::cerr << "error: [job " << job << "]: " << e.what() << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        builder.finish();
        // Gref fragments can only be clipped against the finished topology.
        if (!gref_fragment_paths[job].empty()) {
            job_statistics.ref_paths += add_gref_fragments(this->gbz, gref_fragment_paths[job], parameters.min_gref_length, builder, metadata);
            builder.finish();
        }
        builder.index.addMetadata();
        builder.index.metadata = metadata.get_metadata();
        indexes[job] = builder.index;
        #pragma omp critical
        {
            if (this->verbosity >= Haplotypes::verbosity_detailed) {
                std::cerr << "Job " << job << ": "; job_statistics.print(std::cerr) << std::endl;
            }
            statistics.combine(job_statistics);
        }
    }
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Total: "; statistics.print(std::cerr) << std::endl;
        std::cerr << "Finished the jobs in " << seconds << " seconds" << std::endl;
    }

    // Merge the partial indexes and set the reference samples tag.
    checkpoint = gbwt::readTimer();
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Merging the partial indexes" << std::endl;
    }
    gbwt::GBWT merged(indexes);
    if (parameters.include_reference || !parameters.excluded_chains.empty()) {
        // Excluded chains are copied through verbatim and may contain reference
        // paths, so preserve the reference sample tag in that case as well.
        copy_reference_samples(this->gbz.index, merged);
    }
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Merged the indexes in " << seconds << " seconds" << std::endl;
    }

    if (this->verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Built the GBWT in " << seconds << " seconds" << std::endl;
    }
    return merged;
}

//------------------------------------------------------------------------------

// Return true if the kmer contains a homopolymer of length threshold or more.
// The kmer must be encoded with 2 bits per symbol in the low 2k bits. The high
// bits must be 0.
bool contains_homopolymer(const Haplotypes::Subchain::kmer_type kmer, size_t k, size_t threshold) {
    if (threshold > k) {
        return false;
    }
    if (threshold <= 1) {
        return true;
    }

    // 2-bit field i is 0 iff symbols i and i+1 are the same.
    Haplotypes::Subchain::kmer_type diff = kmer ^ (kmer >> 2);
    // Field i is 1 iff symbols i and i+1 are the same.
    Haplotypes::Subchain::kmer_type equal = ~(diff | (diff >> 1)) & 0x5555555555555555ull;
    // Mask out the last field, as it compares the last symbol to the padding.
    equal &= (Haplotypes::Subchain::kmer_type(1) << (2 * (k - 1))) - 1;

    // A homopolymer of length `threshold` is a run of `threshold - 1` set bits in
    // `equal`, 2 bits apart. We find such runs by shifting and ANDing. After each
    // step, bit i is set when the original `equal` has covered set bits starting
    // at bit i. Each step doubles `covered`, and the last one adds only what is
    // still needed.
    size_t needed = threshold - 1;
    for (size_t covered = 1; covered < needed && equal != 0; ) {
        size_t step = std::min(covered, needed - covered);
        equal &= equal >> (2 * step);
        covered += step;
    }
    return (equal != 0);
}

// Classify the kmers in the subchain according to the kmer counts and the
// coverage. Return a vector of kmer types and their initial scores. Update the
// statistics with the number of non-frequent kmers if necessary.
std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>> classify_kmers(
    const Haplotypes::Subchain& subchain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    double coverage, size_t k,
    Recombinator::Statistics* statistics,
    const Recombinator::Parameters& parameters
) {
    // TODO: What are the proper thresholds?
    double absent_threshold = (parameters.haploid_scoring ? coverage * 0.2 : coverage * 0.1);
    double heterozygous_threshold = (parameters.haploid_scoring ? 0.0 : coverage / std::log(4.0));
    double homozygous_threshold = coverage * 2.5;

    // Determine the type of each kmer in the sample and the score for the kmer.
    // A haplotype with the kmer gets +1.0 * score, while a haplotype without it
    // gets -1.0 * score.
    std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>> kmer_types;
    size_t selected_kmers = 0;
    for (size_t kmer_id = 0; kmer_id < subchain.num_kmers(); kmer_id++) {
        double count = kmer_counts.at(subchain.kmers[kmer_id]);
        if (contains_homopolymer(subchain.kmers[kmer_id], k, parameters.ignore_homopolymers)) {
            kmer_types.push_back({ KmerPresenceMatrix::ignored, 0.0 });
            continue;
        }
        switch (parameters.scoring_model) {
        case Recombinator::Parameters::high_coverage_scoring:
            // High-coverage model: the frequent component is the true signal,
            // while everything below it is contamination. Reward presence of
            // frequent kmers and penalize presence of all others.
            if (count >= homozygous_threshold) {
                kmer_types.push_back({ KmerPresenceMatrix::present, 1.0 });
            } else {
                kmer_types.push_back({ KmerPresenceMatrix::absent, -1.0 * parameters.absent_score });
            }
            selected_kmers++;
            break;
        case Recombinator::Parameters::half_coverage_scoring:
            // Half-coverage model for heterogametic allosomes: the single true
            // copy sits at ~cov/2, i.e. the band the standard model calls
            // heterozygous, so reward that band as present. There is no real
            // heterozygous component outside the PAR. The homozygous (~cov) band
            // and above is paralog / contamination in the body and
            // non-discriminative backbone in the PAR, so it is treated as
            // uninformative (like the frequent band): the selected haplotype
            // still emits that sequence, we just don't let it steer selection.
            if (count < absent_threshold) {
                kmer_types.push_back({ KmerPresenceMatrix::absent, -1.0 * parameters.absent_score });
            } else if (count < heterozygous_threshold) {
                kmer_types.push_back({ KmerPresenceMatrix::present, 1.0 });
            } else {
                kmer_types.push_back({ KmerPresenceMatrix::frequent, 0.0 });
            }
            selected_kmers++;
            break;
        default:
            // Standard scoring model.
            if (count < absent_threshold) {
                kmer_types.push_back({ KmerPresenceMatrix::absent, -1.0 * parameters.absent_score });
                selected_kmers++;
            } else if (count < heterozygous_threshold) {
                kmer_types.push_back({ KmerPresenceMatrix::heterozygous, 0.0 });
                selected_kmers++;
            } else if (count < homozygous_threshold) {
                kmer_types.push_back({ KmerPresenceMatrix::present, 1.0 });
                selected_kmers++;
            } else {
                kmer_types.push_back({ KmerPresenceMatrix::frequent, 0.0 });
            }
            break;
        }
    }
    if (statistics != nullptr) {
        statistics->kmers += selected_kmers;
    }

    return kmer_types;
}

//------------------------------------------------------------------------------

// Select the best pair of haplotypes from the candidates. Each haplotype gets
// +1 for getting a kmer right and -1 for getting it wrong.
// Also returns the remaining non-duplicated haplotypes as extra fragments if
// the appropriate paremeters have been set.
std::vector<std::pair<size_t, double>> select_diploid(
    const gbwtgraph::GBZ& gbz,
    const Haplotypes::Subchain& subchain,
    const std::vector<std::pair<size_t, double>>& candidates,
    const std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>>& kmer_types,
    const Recombinator::Parameters& parameters
) {
    std::int64_t best_score = std::numeric_limits<std::int64_t>::min();
    size_t best_left = 0, best_right = 0;

    // Decode the presence rows of the candidates once, so that scoring a pair
    // takes a few popcounts per 64 kmers.
    std::vector<size_t> sequences;
    sequences.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        sequences.push_back(candidate.first);
    }
    DiploidPairScorer scorer(subchain, sequences, kmer_types);

    // We now consider taking the same haplotype twice. If this is a bad subchain
    // and we sample extra fragments, the number of fragments will depend on
    // whether we take the same haplotype twice.
    for (size_t left = 0; left < candidates.size(); left++) {
        for (size_t right = left; right < candidates.size(); right++) {
            std::int64_t score = scorer.score(left, right);
            if (score > best_score) {
                best_score = score;
                best_left = left;
                best_right = right;
            }
        }
    }

    // If this is a bad subchain, we move the selected haplotypes to the front
    // and return the rest as extra fragments.
    // Otherwise we return only the selected haplotypes.
    if (parameters.extra_fragments) {
        double badness = subchain.badness(gbz);
        if (badness > parameters.badness_threshold) {
            std::vector<std::pair<size_t, double>> result;
            result.push_back(candidates[best_left]);
            result.push_back(candidates[best_right]);
            for (size_t i = 0; i < candidates.size(); i++) {
                if (i != best_left && i != best_right) {
                    result.push_back(candidates[i]);
                }
            }
            return result;
        }
    }

    return { candidates[best_left], candidates[best_right] };
}

//------------------------------------------------------------------------------

namespace {

// A memo of the length and the number of N bases of each node it is asked about, which
// reads the sequence of a node only once.
class NodeBaseCounts {
public:
    explicit NodeBaseCounts(const gbwtgraph::GBZ& gbz) : gbz(gbz) {}

    // Return the length and the number of N bases of the node.
    std::pair<size_t, size_t> operator()(nid_t node) {
        auto [iter, inserted] = this->counts.try_emplace(node);
        if (inserted) {
            std::string_view view = this->gbz.graph.get_sequence_view(this->gbz.graph.get_handle(node));
            iter->second.first = view.size();
            iter->second.second = static_cast<size_t>(std::ranges::count(view, 'N') + std::ranges::count(view, 'n'));
        }
        return iter->second;
    }

    // Return the total length and the total number of N bases of the nodes.
    std::pair<size_t, size_t> totals(const std::vector<nid_t>& nodes) {
        std::pair<size_t, size_t> result(0, 0);
        for (nid_t node : nodes) {
            auto [length, n_bases] = (*this)(node);
            result.first += length;
            result.second += n_bases;
        }
        return result;
    }

private:
    const gbwtgraph::GBZ& gbz;
    std::unordered_map<nid_t, std::pair<size_t, size_t>> counts;
};

} // anonymous namespace

// Return true if the sequence can be selected: it is not from a banned sample or a gref cover.
bool can_select(const gbwtgraph::GBZ& gbz, const Haplotypes::Subchain& subchain, size_t seq_offset, const Recombinator::Parameters& parameters) {
    gbwt::FullPathName path_name = gbz.index.metadata.fullPath(gbwt::Path::id(subchain.sequences[seq_offset].first));
    // A gref cover is never a candidate. Its copy of the reference would take haplotype
    // slots while adding nothing the reference does not, and its fragments are pieces
    // of other haplotypes. Without it, a gref graph samples the same haplotypes as the
    // graph it was computed from.
    return (!parameters.banned_samples.contains(path_name.sample_name) && !is_gref_path(path_name));
}

//------------------------------------------------------------------------------

// Return the sample's count of the kmer for the count model: `pair_model::IGNORED_COUNT`
// for a kmer with a homopolymer that the parameters ignore, and otherwise the count, at
// most `pair_model::CountModel::MAX_COUNT`.
size_t model_count(
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    Haplotypes::Subchain::kmer_type kmer, size_t k, const Recombinator::Parameters& parameters
) {
    if (contains_homopolymer(kmer, k, parameters.ignore_homopolymers)) {
        return pair_model::IGNORED_COUNT;
    }
    auto iter = kmer_counts.find(kmer);
    return std::min((iter == kmer_counts.end() ? size_t(0) : iter->second), pair_model::CountModel::MAX_COUNT);
}

namespace {

// Return which sequences the model uses: those that can be selected, or every sequence if
// none can, minus those with more than `max_n_bp` N bases unless that would remove them
// all. The N bases of each sequence are used only with an N-base limit.
std::vector<bool> usable_sequences(
    const gbwtgraph::GBZ& gbz, const Haplotypes::Subchain& subchain, const std::vector<size_t>& n_bases,
    const Recombinator::Parameters& parameters
) {
    size_t n = subchain.num_sequences();
    std::vector<bool> usable(n, false);
    for (size_t i = 0; i < n; i++) {
        usable[i] = can_select(gbz, subchain, i, parameters);
    }
    auto num_usable = static_cast<size_t>(std::ranges::count(usable, true));
    if (num_usable == 0) {
        usable.assign(n, true);
        num_usable = n;
    }
    if (!parameters.has_n_base_limit()) {
        return usable;
    }
    size_t over_limit = 0;
    for (size_t i = 0; i < n; i++) {
        if (usable[i] && n_bases[i] > parameters.max_n_bp) {
            over_limit++;
        }
    }
    if (over_limit == 0 || over_limit == num_usable) {
        return usable;
    }
    for (size_t i = 0; i < n; i++) {
        if (n_bases[i] > parameters.max_n_bp) {
            usable[i] = false;
        }
    }
    return usable;
}

// Add the counts of the kmers of the subchain to the histogram of counts for the count model.
void add_to_histogram(
    std::vector<double>& histogram, const Haplotypes::Subchain& subchain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts, size_t k, const Recombinator::Parameters& parameters
) {
    for (auto kmer : subchain.kmers) {
        if (size_t count = model_count(kmer_counts, kmer, k, parameters); count != pair_model::IGNORED_COUNT) {
            histogram[count] += 1.0;
        }
    }
}

} // anonymous namespace

pair_model::Candidates model_candidates(
    const gbwtgraph::GBZ& gbz, const Haplotypes::Subchain& subchain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts, const pair_model::CountTable& table,
    const gbwt::FragmentMap& fragment_map, size_t k, const Recombinator::Parameters& parameters
) {
    size_t n = subchain.num_sequences();
    size_t J = subchain.num_kmers();

    // Kmer presence, copies, N bases, length, and haplotype of each sequence. The
    // haplotype of a sequence is the sample and the phase of its path. A path without a
    // phase has the largest phase, whose low 32 bits keep it apart from real phases.
    gbwt::CachedGBWT cache(gbz.index);
    NodeBaseCounts node_bases(gbz);
    const gbwt::Metadata& metadata = gbz.index.metadata;
    std::vector<std::vector<std::uint64_t>> rows(n);
    std::vector<pair_model::KmerCopies> copies(n);
    std::vector<size_t> n_bases(n, 0);
    std::vector<size_t> length(n, 0);
    std::vector<std::uint64_t> haplotype_ids(n);
    for (size_t i = 0; i < n; i++) {
        subchain.presence_row(i, rows[i]);
        pair_model::KmerCopies& sequence_copies = copies[i];
        subchain.for_each_copy(i, [&sequence_copies](size_t kmer_id, size_t value) {
            sequence_copies.emplace_back(kmer_id, value);
        });
        std::tie(length[i], n_bases[i]) = node_bases.totals(
            subchain_walk_nodes(cache, subchain.start, subchain.end, subchain.sequences[i], fragment_map)
        );
        gbwt::PathName name = metadata.path(gbwt::Path::id(subchain.sequences[i].first));
        haplotype_ids[i] = (std::uint64_t(name.sample) << 32) | std::uint64_t(static_cast<std::uint32_t>(name.phase));
    }

    std::vector<bool> usable = usable_sequences(gbz, subchain, n_bases, parameters);

    std::vector<size_t> counts(J);
    for (size_t j = 0; j < J; j++) {
        counts[j] = model_count(kmer_counts, subchain.kmers[j], k, parameters);
    }
    std::ranges::sort(length);
    double median_length = 0.0;
    if (n > 0) {
        median_length = (n % 2 == 1 ? double(length[n / 2]) : 0.5 * (double(length[n / 2 - 1]) + double(length[n / 2])));
    }
    return pair_model::score_pairs(rows, copies, haplotype_ids, usable, counts, table, median_length);
}

std::vector<std::vector<std::pair<size_t, size_t>>> Recombinator::model_pairs(
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts, double coverage, const Parameters& parameters
) const {
    double start = gbwt::readTimer();
    if (!(coverage > 0.0)) {
        throw std::runtime_error("Recombinator::model_pairs(): the count model needs a positive kmer coverage");
    }

    // The model is used in chains sampled with the standard scoring model.
    std::vector<size_t> chains;
    for (const Haplotypes::TopLevelChain& chain : this->haplotypes.chains) {
        bool full = (chain.subchains.size() == 1 && chain.subchains.front().type == Haplotypes::Subchain::full_haplotype);
        bool other = (parameters.excluded_chains.contains(chain.offset) || parameters.high_coverage_chains.contains(chain.offset) ||
            parameters.half_coverage_chains.contains(chain.offset));
        if (!other && !full && !chain.subchains.empty()) {
            chains.push_back(chain.offset);
        }
    }

    // Fit the count model to the counts of the kmers in those chains.
    std::vector<double> histogram(pair_model::CountModel::MAX_COUNT + 1, 0.0);
    for (size_t chain_id : chains) {
        for (const Haplotypes::Subchain& subchain : this->haplotypes.chains[chain_id].subchains) {
            add_to_histogram(histogram, subchain, kmer_counts, this->haplotypes.k(), parameters);
        }
    }
    pair_model::CountModel count_model = pair_model::CountModel::fit(histogram, coverage);
    pair_model::CountTable table = pair_model::count_table(count_model, histogram);
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Count model: mu " << count_model.mu << ", eps " << count_model.eps << ", r1 " << count_model.r1
            << ", r0 " << count_model.r0 << ", pi";
        for (double p : count_model.pi) {
            std::cerr << " " << p;
        }
        std::cerr << " (" << (gbwt::readTimer() - start) << " seconds)" << std::endl;
    }

    // Find the candidates in the subchains of each chain in parallel, and take the Viterbi path.
    std::vector<std::vector<std::pair<size_t, size_t>>> result(this->haplotypes.chains.size());
    size_t total = 0;
    for (size_t chain_id : chains) {
        const Haplotypes::TopLevelChain& chain = this->haplotypes.chains[chain_id];
        std::vector<pair_model::Candidates> candidates(chain.subchains.size());
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
            candidates[subchain_id] = model_candidates(
                this->gbz, chain.subchains[subchain_id], kmer_counts, table, this->fragment_map, this->haplotypes.k(), parameters
            );
        }
        std::vector<std::pair<size_t, size_t>> path = pair_model::viterbi(candidates, MODEL_TAU, MODEL_RHO);
        result[chain_id] = pair_model::choose_sequences(candidates, path);
        total += result[chain_id].size();
    }
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Model pairs for " << total << " subchains in " << (gbwt::readTimer() - start) << " seconds" << std::endl;
    }
    return result;
}

//------------------------------------------------------------------------------

namespace {

// Return the sequences in the subchain that can be selected, as (sequence offset, score 0).
std::vector<std::pair<size_t, double>> selectable_sequences(
    const gbwtgraph::GBZ& gbz, const Haplotypes::Subchain& subchain, const Recombinator::Parameters& parameters
) {
    std::vector<std::pair<size_t, double>> result;
    for (size_t seq_offset = 0; seq_offset < subchain.num_sequences(); seq_offset++) {
        if (can_select(gbz, subchain, seq_offset, parameters)) {
            result.emplace_back(seq_offset, 0.0);
        }
    }
    return result;
}

// Return the sequences in the subchain that get the continuity bonus: those a generated
// haplotype took in the previous subchain.
std::vector<bool> continuing_sequences(
    const Haplotypes::Subchain& subchain, const Recombinator::SelectionContext& context, const Recombinator::Parameters& parameters
) {
    std::vector<bool> result(subchain.num_sequences(), false);
    if (context.sources == nullptr || !(parameters.continuity_bonus > 0.0)) {
        return result;
    }
    hash_set<gbwt::size_type> sources;
    sources.insert(context.sources->begin(), context.sources->end());
    for (size_t seq_offset = 0; seq_offset < subchain.num_sequences(); seq_offset++) {
        result[seq_offset] = (sources.find(subchain.sequences[seq_offset].first) != sources.end());
    }
    return result;
}

// Remove the sequence from the candidates. Return false if it was not a candidate.
bool remove_candidate(std::vector<std::pair<size_t, double>>& candidates, size_t seq_offset) {
    auto found = std::ranges::find_if(candidates, [seq_offset](const std::pair<size_t, double>& candidate) {
        return candidate.first == seq_offset;
    });
    if (found == candidates.end()) {
        return false;
    }
    candidates.erase(found);
    return true;
}

// Update the kmer scores after the sequence has been selected: a heterozygous kmer
// becomes less wanted if the sequence has it and more wanted if not, and a present kmer
// the sequence has is discounted.
void adjust_kmer_scores(
    const Haplotypes::Subchain& subchain, size_t seq_offset,
    std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>>& kmer_types,
    const Recombinator::Parameters& parameters
) {
    subchain.for_each_kmer(seq_offset, [&kmer_types, &parameters](size_t kmer_id, bool is_present) {
        switch (kmer_types[kmer_id].first) {
        case KmerPresenceMatrix::heterozygous:
            kmer_types[kmer_id].second += (is_present ? -1.0 : 1.0) * parameters.het_adjustment;
            break;
        case KmerPresenceMatrix::present:
            if (is_present) {
                kmer_types[kmer_id].second *= parameters.present_discount;
            }
            break;
        default:
            break;
        }
    });
}

// A judge of whether a candidate in a subchain stays within the new-bp and N-base limits
// (`max_new_bp`, `bp_escape` and `max_n_bp` in `Recombinator::Parameters`), given the
// haplotypes already selected there. It is told each selected haplotype, and it finds
// the nodes of a candidate when first asked about it.
class SelectionLimits {
public:
    SelectionLimits(
        const gbwtgraph::GBZ& gbz, const Haplotypes::Subchain& subchain, const gbwt::FragmentMap* fragment_map,
        const std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>>& kmer_types,
        const std::vector<std::pair<size_t, double>>& candidates, const Recombinator::Parameters& parameters
    ) : subchain(subchain), fragment_map(fragment_map), node_bases(gbz), parameters(parameters),
        escape(parameters.has_new_bp_limit() && parameters.bp_escape > 0.0),
        candidate_nodes((parameters.has_new_bp_limit() || parameters.has_n_base_limit()) ? subchain.num_sequences() : 0),
        words((subchain.num_kmers() + 63) / 64) {
        if (parameters.has_new_bp_limit() || parameters.has_n_base_limit()) {
            if (fragment_map == nullptr) {
                throw std::invalid_argument("select_haplotypes(): the new-bp and N-base limits need a fragment map");
            }
            this->cache.emplace(gbz.index);
        }
        if (!this->escape) {
            return;
        }
        this->sample_bits.assign(this->words, 0);
        this->selected_bits.assign(this->words, 0);
        for (size_t kmer_id = 0; kmer_id < kmer_types.size(); kmer_id++) {
            KmerPresenceMatrix::kmer_classification type = kmer_types[kmer_id].first;
            if (type == KmerPresenceMatrix::present || type == KmerPresenceMatrix::heterozygous) {
                this->sample_bits[kmer_id / 64] |= std::uint64_t(1) << (kmer_id % 64);
            }
        }
        this->candidate_bits.resize(subchain.num_sequences());
        for (const auto& [seq_offset, score] : candidates) {
            subchain.presence_row(seq_offset, this->candidate_bits[seq_offset]);
        }
    }

    // Return true if the candidate is within the limits. The new-bp limit does not
    // apply to the first haplotype.
    bool allows(size_t seq_offset, bool first) {
        bool within_bp = (!this->parameters.has_new_bp_limit() || first || this->within_bp_limit(seq_offset));
        return within_bp && (!this->parameters.has_n_base_limit() || this->within_n_limit(seq_offset));
    }

    // Count a selected haplotype among those that candidates are compared with.
    void add_selected(size_t seq_offset) {
        if (this->escape) {
            for (size_t w = 0; w < this->words; w++) {
                this->selected_bits[w] |= this->candidate_bits[seq_offset][w];
            }
        }
        if (this->parameters.has_new_bp_limit()) {
            for (nid_t node : this->nodes_of(seq_offset)) {
                this->selected_nodes.insert(node);
            }
        }
    }

private:
    const Haplotypes::Subchain& subchain;
    const gbwt::FragmentMap* fragment_map;

    // For walking the candidates through the subchain, with a limit.
    std::optional<gbwt::CachedGBWT> cache;

    // Finds the lengths and the N bases of the nodes.
    NodeBaseCounts node_bases;

    const Recombinator::Parameters& parameters;

    // The escape from the new-bp limit applies.
    bool escape;

    // The nodes of each candidate, once found.
    std::vector<std::optional<std::vector<nid_t>>> candidate_nodes;

    // The nodes the selected haplotypes visit.
    hash_set<nid_t> selected_nodes;

    // Kmer presence as 64-bit words: the sample kmers, the kmers of each candidate, and
    // the kmers of the selected haplotypes. Set only with the escape.
    size_t words;
    std::vector<std::uint64_t> sample_bits;
    std::vector<std::vector<std::uint64_t>> candidate_bits;
    std::vector<std::uint64_t> selected_bits;

    const std::vector<nid_t>& nodes_of(size_t seq_offset) {
        std::optional<std::vector<nid_t>>& nodes = this->candidate_nodes[seq_offset];
        if (!nodes) {
            nodes = subchain_walk_nodes(
                *this->cache, this->subchain.start, this->subchain.end, this->subchain.sequences[seq_offset], *this->fragment_map
            );
        }
        return *nodes;
    }

    bool within_bp_limit(size_t seq_offset) {
        size_t new_bp = 0;
        for (nid_t node : this->nodes_of(seq_offset)) {
            if (this->selected_nodes.find(node) == this->selected_nodes.end()) {
                new_bp += this->node_bases(node).first;
            }
        }
        if (new_bp <= this->parameters.max_new_bp) {
            return true;
        }
        if (!this->escape) {
            return false;
        }
        size_t new_sample = 0;
        for (size_t w = 0; w < this->words; w++) {
            new_sample += __builtin_popcountll(this->candidate_bits[seq_offset][w] & ~this->selected_bits[w] & this->sample_bits[w]);
        }
        return (double(new_sample) * this->parameters.bp_escape >= double(new_bp));
    }

    bool within_n_limit(size_t seq_offset) {
        return (this->node_bases.totals(this->nodes_of(seq_offset)).second <= this->parameters.max_n_bp);
    }
};

// Return the index of the highest-scoring candidate within the limits, trying the
// candidates best first, so that the limits find the nodes of a candidate only when
// needed. Ties go to the candidate listed first. If no candidate is within the limits,
// return the best one for the first haplotype, and nothing otherwise.
std::optional<size_t> best_within_limits(
    const std::vector<std::pair<size_t, double>>& candidates, SelectionLimits& limits, bool first
) {
    if (candidates.empty()) {
        return std::nullopt;
    }
    // The best candidate is usually within the limits, so we try it before sorting.
    auto best = static_cast<size_t>(std::ranges::max_element(candidates, {}, &std::pair<size_t, double>::second) - candidates.begin());
    if (limits.allows(candidates[best].first, first)) {
        return best;
    }
    std::vector<size_t> order(candidates.size());
    std::iota(order.begin(), order.end(), size_t(0));
    std::ranges::stable_sort(order, [&candidates](size_t a, size_t b) {
        return (candidates[a].second > candidates[b].second);
    });
    for (size_t i : order) {
        if (i != best && limits.allows(candidates[i].first, first)) {
            return i;
        }
    }
    if (first) {
        return best;
    }
    return std::nullopt;
}

} // anonymous namespace

std::vector<std::pair<size_t, double>> select_haplotypes(
    const gbwtgraph::GBZ& gbz,
    const Haplotypes::Subchain& subchain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    double coverage, size_t k,
    Recombinator::Statistics* statistics,
    std::vector<Recombinator::LocalHaplotype>* local_haplotypes,
    const Recombinator::Parameters& parameters,
    const Recombinator::SelectionContext& context
) {
    // Classify the kmers.
    std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>> kmer_types = classify_kmers(
        subchain, kmer_counts, coverage, k, statistics, parameters
    );

    // Select the haplotypes greedily.
    std::vector<std::pair<size_t, double>> selected_haplotypes;
    std::vector<std::pair<size_t, double>> remaining_haplotypes = selectable_sequences(gbz, subchain, parameters);
    std::vector<double> scores(subchain.num_sequences(), 0.0);
    std::vector<bool> continues = continuing_sequences(subchain, context, parameters);
    SelectionLimits limits(gbz, subchain, context.fragment_map, kmer_types, remaining_haplotypes, parameters);
    auto select = [&](std::pair<size_t, double> pick) {
        selected_haplotypes.push_back(pick);
        limits.add_selected(pick.first);
        adjust_kmer_scores(subchain, pick.first, kmer_types, parameters);
    };

    // The model pair becomes the first two haplotypes. A haplotype may pair with itself.
    // The members get the kmer scores they have before either is selected.
    if (context.pair != nullptr && !remaining_haplotypes.empty()) {
        subchain.score_sequences(kmer_types, [&scores](size_t haplotype_id, double score) {
            scores[haplotype_id] = score;
        });
        for (size_t seq_offset : { context.pair->first, context.pair->second }) {
            if (bool repeated = (!selected_haplotypes.empty() && selected_haplotypes.front().first == seq_offset);
                !remove_candidate(remaining_haplotypes, seq_offset) && !repeated) {
                throw std::runtime_error("select_haplotypes(): the model pair includes a haplotype that cannot be selected");
            }
            select({ seq_offset, scores[seq_offset] });
        }
    }

    while (selected_haplotypes.size() < parameters.num_haplotypes && !remaining_haplotypes.empty()) {
        // Score the haplotypes. We score all of them, as that is cheaper than
        // determining which ones are still relevant.
        subchain.score_sequences(
            kmer_types,
            [&](size_t haplotype_id, double score) {
                scores[haplotype_id] = score;
            }
        );
        for (auto& [seq_offset, score] : remaining_haplotypes) {
            score = scores[seq_offset];
            if (continues[seq_offset]) {
                score += parameters.continuity_bonus;
            }
        }

        // Report ranks and scores for each remaining haplotype.
        if (local_haplotypes != nullptr) {
            auto copy = remaining_haplotypes;
            std::sort(
                copy.begin(), copy.end(),
                [](std::pair<size_t, double> a, std::pair<size_t, double> b) -> bool {
                    return (a.second > b.second);
                }
            );
            for (size_t i = 0; i < copy.size(); i++) {
                size_t sequence_id = copy[i].first;
                double score = copy[i].second;
                (*local_haplotypes)[sequence_id].scores.emplace_back(i, score);
            }
        }

        // Select the highest-scoring haplotype within the limits.
        std::optional<size_t> selected = best_within_limits(remaining_haplotypes, limits, selected_haplotypes.empty());
        if (!selected.has_value()) {
            break;
        }
        // The selected haplotype keeps its kmer score without the continuity bonus.
        size_t pick = remaining_haplotypes[*selected].first;
        remaining_haplotypes.erase(remaining_haplotypes.begin() + *selected);
        select({ pick, scores[pick] });
    }

    if (parameters.diploid_sampling) {
        // Do diploid sampling. If this is a bad subchain, we also return the
        // extra fragments starting from the third haplotype.
        return select_diploid(gbz, subchain, selected_haplotypes, kmer_types, parameters);
    } else {
        // If we did not have enough haplotypes in the subchain, repeat them as necessary.
        size_t original_selected = selected_haplotypes.size();
        for (size_t i = original_selected; i < parameters.num_haplotypes; i++) {
            auto next = selected_haplotypes[i % original_selected];
            selected_haplotypes.push_back(next);
        }
        return selected_haplotypes;
    }
}

void thread_by_similarity(
    const Haplotypes::Subchain& subchain,
    const std::vector<std::pair<size_t, double>>& selected_haplotypes,
    const std::vector<gbwt::size_type>& sources,
    std::vector<size_t>& haplotype_to_selected, sdsl::bit_vector& selected_in_use
) {
    size_t num_haplotypes = sources.size();
    std::map<size_t, std::vector<std::uint64_t>> rows;
    auto row_of = [&](size_t offset) -> const std::vector<std::uint64_t>& {
        auto [iter, inserted] = rows.try_emplace(offset);
        if (inserted) {
            subchain.presence_row(offset, iter->second);
        }
        return iter->second;
    };

    // Unused selections by sequence offset, with the lowest index last.
    std::map<size_t, std::vector<size_t>> free_entries;
    for (size_t selected = num_haplotypes; selected > 0; selected--) {
        if (!selected_in_use[selected - 1]) {
            free_entries[selected_haplotypes[selected - 1].first].push_back(selected - 1);
        }
    }

    // (distance, generated haplotype, sequence offset) for each candidate match.
    std::vector<std::tuple<std::uint64_t, size_t, size_t>> pairs;
    for (size_t haplotype = 0; haplotype < num_haplotypes; haplotype++) {
        if (haplotype_to_selected[haplotype] < num_haplotypes || sources[haplotype] == gbwt::invalid_sequence()) {
            continue;
        }
        size_t source = subchain.sequences.size();
        for (size_t i = 0; i < subchain.sequences.size(); i++) {
            if (subchain.sequences[i].first == sources[haplotype]) {
                source = i;
                break;
            }
        }
        if (source == subchain.sequences.size()) {
            continue;
        }
        const std::vector<std::uint64_t>& a = row_of(source);
        for (const auto& [offset, entries] : free_entries) {
            const std::vector<std::uint64_t>& b = row_of(offset);
            std::uint64_t distance = 0;
            for (size_t w = 0; w < a.size(); w++) {
                distance += __builtin_popcountll(a[w] ^ b[w]);
            }
            pairs.emplace_back(distance, haplotype, offset);
        }
    }
    std::ranges::sort(pairs);
    for (const auto& [distance, haplotype, offset] : pairs) {
        std::vector<size_t>& entries = free_entries[offset];
        if (haplotype_to_selected[haplotype] < num_haplotypes || entries.empty()) {
            continue;
        }
        haplotype_to_selected[haplotype] = entries.back();
        selected_in_use[entries.back()] = true;
        entries.pop_back();
    }
}

namespace {

// Return the GBWT sequence each generated haplotype took in the previous subchain, or
// `gbwt::invalid_sequence()`.
std::vector<gbwt::size_type> current_sequences(const std::vector<RecombinatorHaplotype>& haplotypes) {
    std::vector<gbwt::size_type> result;
    for (const RecombinatorHaplotype& haplotype : haplotypes) {
        result.push_back(haplotype.sequence_id);
    }
    return result;
}

} // anonymous namespace

Recombinator::Statistics Recombinator::generate_haplotypes(const Haplotypes::TopLevelChain& chain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata,
    const Parameters& parameters,
    double coverage,
    const std::vector<std::pair<size_t, size_t>>& pairs
) const {
    size_t final_haplotypes = (parameters.diploid_sampling ? 2 : parameters.num_haplotypes);
    bool wrap = (parameters.wrap_contigs.find(chain.contig_name) != parameters.wrap_contigs.end());
    std::vector<RecombinatorHaplotype> haplotypes;
    for (size_t i = 0; i < final_haplotypes; i++) {
        haplotypes.emplace_back(*this, builder, metadata, chain.contig_name, wrap, i + 1);
    }

    Statistics statistics;
    statistics.chains = 1; statistics.haplotypes = haplotypes.size();
    if (chain.subchains.size() == 1 && chain.subchains.front().type == Haplotypes::Subchain::full_haplotype) {
        // Full haplotypes should all be identical, because there are no snarls.
        // Therefore we do not need kmers.
        auto& subchain = chain.subchains.front();
        for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
            assert(!subchain.sequences.empty());
            size_t seq = haplotype % subchain.num_sequences();
            haplotypes[haplotype].take(subchain.sequences[seq].first);
        }
        statistics.full_haplotypes = 1;
    } else {
        bool have_haplotypes = false;
        std::vector<RecombinatorFragment> extra_fragments;
        for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
            const auto& subchain = chain.subchains[subchain_id];
            if (subchain.type == Haplotypes::Subchain::full_haplotype) {
                throw std::runtime_error("Recombinator::generate_haplotypes(): nontrivial chain " + std::to_string(chain.offset) + " contains a subchain with full haplotypes");
            }
            assert(!subchain.sequences.empty());

            // Select the haplotypes greedily. If we are doing diploid sampling and we get
            // extra fragments, store them for later processing.
            std::vector<gbwt::size_type> sources = current_sequences(haplotypes);
            SelectionContext context;
            context.fragment_map = &this->fragment_map;
            context.sources = &sources;
            context.pair = (pairs.empty() ? nullptr : &pairs[subchain_id]);
            std::vector<std::pair<size_t, double>> selected_haplotypes = select_haplotypes(
                this->gbz, subchain, kmer_counts, coverage, this->haplotypes.k(), &statistics, nullptr, parameters, context
            );
#ifdef debug_selected_haplotypes
            for (size_t i = 0; i < selected_haplotypes.size(); i++) {
                gbwt::size_type sequence_id = subchain.sequences[selected_haplotypes[i].first].first;
                gbwt::size_type path_id = gbwt::Path::id(sequence_id);
                gbwt::FullPathName path_name = this->gbz.index.metadata.fullPath(path_id);
                cerr << "Selected haplotype " << path_name.contig_name
                     << " with score " << selected_haplotypes[i].second << endl;
            }
#endif
            if (parameters.diploid_sampling && selected_haplotypes.size() > 2) {
                for (size_t i = 2; i < selected_haplotypes.size(); i++) {
                    gbwt::size_type sequence_id = subchain.sequences[selected_haplotypes[i].first].first;
                    gbwt::edge_type start;
                    if (subchain.has_start()) {
                        start = gbwt::edge_type(subchain.start, subchain.sequences[selected_haplotypes[i].first].second);
                    } else {
                        // This is a prefix.
                        start = this->gbz.index.start(sequence_id);
                    }
                    extra_fragments.emplace_back(sequence_id, start, subchain.end, subchain_id, i - 1);
                }
                statistics.bad_subchains++;
                selected_haplotypes.resize(2);
            }

            // Try to match the existing haplotypes with the selected sequences based on
            // GBWT sequence id.
            // TODO: There are often many equally good haplotypes, and because we choose
            // a different subset of them in each subchain, we often cannot connect the
            // same haplotype from subchain to subchain, even when it would be one of the
            // top ones in each of them.
            std::vector<size_t> haplotype_to_selected(haplotypes.size(), haplotypes.size());
            sdsl::bit_vector selected_in_use(haplotypes.size(), 0);
            for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                for (size_t selected = 0; selected < haplotypes.size(); selected++) {
                    if (subchain.sequences[selected_haplotypes[selected].first].first == haplotypes[haplotype].sequence_id && !selected_in_use[selected]) {
                        haplotype_to_selected[haplotype] = selected;
                        selected_in_use[selected] = 1;
                        statistics.connections++;
                        break;
                    }
                }
            }
            if (parameters.similar_threading) {
                thread_by_similarity(subchain, selected_haplotypes, sources, haplotype_to_selected, selected_in_use);
            }
            for (size_t haplotype = 0, selected = 0; haplotype < haplotypes.size(); haplotype++) {
                if (haplotype_to_selected[haplotype] < haplotypes.size()) {
                    continue;
                }
                while (selected < haplotypes.size() && selected_in_use[selected]) {
                    selected++;
                }
                assert(selected < haplotypes.size());
                haplotype_to_selected[haplotype] = selected;
                selected_in_use[selected] = 1;
                selected++;
            }

            // Finally extend the haplotypes with the selected and matched sequences.
            bool is_last_subchain = (subchain_id + 1 == chain.subchains.size());
            for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                size_t selected = haplotype_to_selected[haplotype];
                size_t seq_offset = selected_haplotypes[selected].first;
                statistics.score += selected_haplotypes[selected].second;
                haplotypes[haplotype].extend(subchain.sequences[seq_offset], subchain, is_last_subchain);
            }
            have_haplotypes = subchain.has_end();
            statistics.subchains++;
        }
        if (have_haplotypes) {
            for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                haplotypes[haplotype].finish(true);
            }
        }

        // Each haplotype should consist of at least one fragment.
        statistics.fragments = 0;
        for (const RecombinatorHaplotype& haplotype : haplotypes) {
            statistics.fragments += haplotype.fragment;
        }

        // Add the extra fragments as separate paths.
        statistics.extra_fragments = 0;
        for (const RecombinatorFragment& fragment : extra_fragments) {
            statistics.extra_fragments += fragment.generate(this->gbz.index, fragment_map, builder, metadata, chain.contig_name);
        }
    }

    return statistics;
}

//------------------------------------------------------------------------------

Recombinator::Statistics Recombinator::copy_chain(const Haplotypes::TopLevelChain& chain,
    gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata
) const {
    // Collect the distinct paths crossing the chain. Each GBWT sequence in a
    // subchain corresponds to a path; the same path may appear in multiple
    // subchains and in either orientation.
    std::unordered_set<gbwt::size_type> path_ids;
    for (const Haplotypes::Subchain& subchain : chain.subchains) {
        for (const Haplotypes::compact_sequence_type& sequence : subchain.sequences) {
            path_ids.insert(gbwt::Path::id(sequence.first));
        }
    }

    // Copy each path through verbatim, preserving its metadata. These are
    // haplotype-sense paths, not reference paths, so we count them separately.
    Statistics statistics;
    statistics.chains = 1;
    for (gbwt::size_type path_id : path_ids) {
        add_path(this->gbz.index, path_id, builder, metadata);
        statistics.copied_paths++;
    }

    return statistics;
}

//------------------------------------------------------------------------------

} // namespace vg

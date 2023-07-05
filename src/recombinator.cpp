#include "recombinator.hpp"

#include "kff.hpp"
#include "statistics.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr std::uint32_t Haplotypes::Header::MAGIC_NUMBER;
constexpr std::uint32_t Haplotypes::Header::VERSION;
constexpr std::uint64_t Haplotypes::Header::DEFAULT_K;

constexpr size_t HaplotypePartitioner::SUBCHAIN_LENGTH;
constexpr size_t HaplotypePartitioner::APPROXIMATE_JOBS;

constexpr size_t Recombinator::NUM_HAPLOTYPES;
constexpr size_t Recombinator::COVERAGE;
constexpr size_t Recombinator::KFF_BLOCK_SIZE;
constexpr double Recombinator::PRESENT_DISCOUNT;
constexpr double Recombinator::HET_ADJUSTMENT;
constexpr double Recombinator::ABSENT_SCORE;

//------------------------------------------------------------------------------

// Returns a GBWTGraph handle as a string (id, orientation).
std::string to_string(handle_t handle) {
    gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(handle);
    return std::string("(") + std::to_string(gbwt::Node::id(node)) + std::string(", ") + std::to_string(gbwt::Node::is_reverse(node)) + std::string(")");
}

//------------------------------------------------------------------------------

hash_map<Haplotypes::Subchain::kmer_type, size_t>::iterator
find_kmer(hash_map<Haplotypes::Subchain::kmer_type, size_t>& counts, Haplotypes::Subchain::kmer_type kmer, size_t k) {
    Haplotypes::Subchain::kmer_type rc = minimizer_reverse_complement(kmer, k);
    auto forward = counts.find(kmer);
    auto reverse = counts.find(rc);
    return (forward != counts.end() ? forward : reverse);
}

hash_map<Haplotypes::Subchain::kmer_type, size_t> Haplotypes::kmer_counts(const std::string& kff_file, Verbosity verbosity) const {
    double start = gbwt::readTimer();
    if (verbosity >= verbosity_basic) {
        std::cerr << "Reading kmer counts" << std::endl;
    }

    // Open and validate the kmer count file.
    ParallelKFFReader reader(kff_file);

    // Populate the map with the kmers we are interested in.
    double checkpoint = gbwt::readTimer();
    hash_map<Subchain::kmer_type, size_t> result;
    result.reserve(this->header.total_kmers);
    for (size_t chain_id = 0; chain_id < this->chains.size(); chain_id++) {
        const TopLevelChain& chain = this->chains[chain_id];
        for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
            const Subchain& subchain = chain.subchains[subchain_id];
            for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                result[subchain.kmers[kmer_id].first] = 0;
            }
        }
    }
    if (verbosity >= verbosity_detailed) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Initialized the hash map with " << result.size() << " kmers in " << seconds << " seconds" << std::endl;
    }

    // Read the KFF file and add the counts using multiple threads.
    checkpoint = gbwt::readTimer();
    size_t kmer_count = 0;
    #pragma omp parallel
    {
        #pragma omp task
        {
            while (true) {
                std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> block = reader.read(Recombinator::KFF_BLOCK_SIZE);
                if (block.empty()) {
                    break;
                }
                std::vector<std::pair<hash_map<Subchain::kmer_type, size_t>::iterator, size_t>> buffer;
                for (auto kmer : block) {
                    auto iter = find_kmer(result, kmer.first, this->k());
                    if (iter != result.end()) {
                        buffer.push_back({ iter, kmer.second });
                    }
                }
                #pragma omp critical
                {
                    for (auto to_update : buffer) {
                        to_update.first->second += to_update.second;
                    }
                    kmer_count += block.size();
                }
            }
        }
    }
    if (verbosity >= verbosity_detailed) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Read " << kmer_count << " kmers in " << seconds << " seconds" << std::endl;
    }

    if (verbosity >= verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Read the kmer counts in " << seconds << " seconds" << std::endl;
    }
    return result;
}

//------------------------------------------------------------------------------

std::string Haplotypes::Subchain::to_string() const {
    std::string result;
    switch (this->type) {
    case normal:
        result.append("normal");
        break;
    case prefix:
        result.append("prefix");
        break;
    case suffix:
        result.append("suffix");
        break;
    case full_haplotype:
        result.append("full");
        break;
    default:
        result.append("invalid");
        break;
    }

    result.append(" from ");
    result.append(to_string_gbwtgraph(this->start));
    result.append(" to ");
    result.append(to_string_gbwtgraph(this->end));

    return result;
}

void Haplotypes::Subchain::simple_sds_serialize(std::ostream& out) const {
    sdsl::simple_sds::serialize_value<std::uint64_t>(this->type, out);
    sdsl::simple_sds::serialize_value<gbwt::node_type>(this->start, out);
    sdsl::simple_sds::serialize_value<gbwt::node_type>(this->end, out);
    sdsl::simple_sds::serialize_vector(this->kmers, out);
    sdsl::simple_sds::serialize_vector(this->sequences, out);
    this->kmers_present.simple_sds_serialize(out);
}

void Haplotypes::Subchain::simple_sds_load(std::istream& in) {
    std::uint64_t temp = sdsl::simple_sds::load_value<std::uint64_t>(in);
    switch (temp) {
    case normal: // Fall through.
    case prefix: // Fall through.
    case suffix: // Fall through.
    case full_haplotype:
        this->type = static_cast<subchain_t>(temp);
        break;
    default:
        throw sdsl::simple_sds::InvalidData("Invalid subchain type: " + std::to_string(temp));
    }

    this->start = sdsl::simple_sds::load_value<gbwt::node_type>(in);
    this->end = sdsl::simple_sds::load_value<gbwt::node_type>(in);
    bool should_have_start = (this->type == normal || this->type == suffix);
    bool should_have_end = (this->type == normal || this->type == prefix);
    if ((this->start != gbwt::ENDMARKER) != should_have_start) {
        throw sdsl::simple_sds::InvalidData("Subchain start node " + std::to_string(this->start) + " does not match type " + std::to_string(temp));
    }
    if ((this->end != gbwt::ENDMARKER) != should_have_end) {
        throw sdsl::simple_sds::InvalidData("Subchain end node " + std::to_string(this->end) + " does not match type" + std::to_string(temp));
    }

    this->kmers = sdsl::simple_sds::load_vector<std::pair<kmer_type, size_t>>(in);
    this->sequences = sdsl::simple_sds::load_vector<sequence_type>(in);
    this->kmers_present.simple_sds_load(in);
    if (kmers_present.size() != kmers.size() * sequences.size()) {
        throw sdsl::simple_sds::InvalidData("Invalid length for the kmer presence bitvector in subchain from " +
            std::to_string(this->start) + " to " + std::to_string(this->end));
    }
}

size_t Haplotypes::Subchain::simple_sds_size() const {
    size_t result = sdsl::simple_sds::value_size<std::uint64_t>() + 2 * sdsl::simple_sds::value_size<gbwt::node_type>();
    result += sdsl::simple_sds::vector_size(this->kmers);
    result += sdsl::simple_sds::vector_size(this->sequences);
    result += this->kmers_present.simple_sds_size();
    return result;
}

void Haplotypes::TopLevelChain::simple_sds_serialize(std::ostream& out) const {
    sdsl::simple_sds::serialize_value<size_t>(this->offset, out);
    sdsl::simple_sds::serialize_value<size_t>(this->job_id, out);
    sdsl::simple_sds::serialize_value<size_t>(this->subchains.size(), out);
    for (auto& subchain : this->subchains) {
        subchain.simple_sds_serialize(out);
    }
}

void Haplotypes::TopLevelChain::simple_sds_load(std::istream& in) {
    this->offset = sdsl::simple_sds::load_value<size_t>(in);
    this->job_id = sdsl::simple_sds::load_value<size_t>(in);
    size_t subchain_count = sdsl::simple_sds::load_value<size_t>(in);
    this->subchains.resize(subchain_count);
    for (size_t i = 0; i < subchain_count; i++) {
        this->subchains[i].simple_sds_load(in);
    }
}

size_t Haplotypes::TopLevelChain::simple_sds_size() const {
    size_t result = 3 * sdsl::simple_sds::value_size<size_t>();
    for (auto& subchain : this->subchains) {
        result += subchain.simple_sds_size();
    }
    return result;
}

void Haplotypes::simple_sds_serialize(std::ostream& out) const {
    sdsl::simple_sds::serialize_value<Header>(this->header, out);
    sdsl::simple_sds::serialize_vector(this->jobs_for_cached_paths, out);
    for (auto& chain : this->chains) {
        chain.simple_sds_serialize(out);
    }
}

void Haplotypes::simple_sds_load(std::istream& in) {
    this->header = sdsl::simple_sds::load_value<Header>(in);
    if (this->header.magic_number != Header::MAGIC_NUMBER) {
        throw sdsl::simple_sds::InvalidData("Haplotypes::simple_sds_load(): Expected magic number " + std::to_string(Header::MAGIC_NUMBER) +
            ", got " + std::to_string(this->header.magic_number));
    }
    if (this->header.version != Header::VERSION) {
        throw sdsl::simple_sds::InvalidData("Haplotypes::simple_sds_load(): Expected version " + std::to_string(Header::VERSION) +
            ", got " + std::to_string(this->header.version));
    }

    this->jobs_for_cached_paths = sdsl::simple_sds::load_vector<size_t>(in);

    this->chains.resize(this->header.top_level_chains);
    for (auto& chain : this->chains) {
        chain.simple_sds_load(in);
    }
}

size_t Haplotypes::simple_sds_size() const {
    size_t result = sdsl::simple_sds::value_size<Header>();
    result += sdsl::simple_sds::vector_size(this->jobs_for_cached_paths);
    for (auto& chain : this->chains) {
        result += chain.simple_sds_size();
    }
    return result;
}

//------------------------------------------------------------------------------

HaplotypePartitioner::HaplotypePartitioner(const gbwtgraph::GBZ& gbz,
    const gbwt::FastLocate& r_index,
    const SnarlDistanceIndex& distance_index,
    const minimizer_index_type& minimizer_index,
    Verbosity verbosity) :
    gbz(gbz), r_index(r_index), distance_index(distance_index), minimizer_index(minimizer_index),
    verbosity(verbosity)
{
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

    Haplotypes result;
    result.header.k = this->minimizer_index.k();

    // Determine top-level chains and assign them to jobs.
    double start = gbwt::readTimer();
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Determining construction jobs" << std::endl;
    }
    size_t size_bound = this->gbz.graph.get_node_count() / parameters.approximate_jobs;
    gbwtgraph::ConstructionJobs jobs = gbwtgraph::gbwt_construction_jobs(this->gbz.graph, size_bound);
    auto chains_by_job = gbwtgraph::partition_chains(this->distance_index, this->gbz.graph, jobs);
    result.header.top_level_chains = jobs.components;
    result.header.construction_jobs = jobs.size();
    result.chains.resize(result.components());
    for (size_t chain_id = 0; chain_id < result.components(); chain_id++) {
        result.chains[chain_id].offset = chain_id;
        result.chains[chain_id].job_id = result.jobs(); // Not assigned to any job.
    }
    for (size_t job_id = 0; job_id < result.jobs(); job_id++) {
        for (auto& chain : chains_by_job[job_id]) {
            result.chains[chain.offset].job_id = job_id;
        }
    }

    // Assign named and reference paths to jobs.
    result.jobs_for_cached_paths.reserve(this->gbz.graph.named_paths.size());
    for (size_t i = 0; i < this->gbz.graph.named_paths.size(); i++) {
        const gbwtgraph::NamedPath& path = this->gbz.graph.named_paths[i];
        if (path.from == gbwt::invalid_edge()) {
            // Skip empty paths.
            result.jobs_for_cached_paths.push_back(result.jobs());
            continue;
        }
        nid_t node_id = gbwt::Node::id(path.from.first);
        auto iter = jobs.node_to_job.find(node_id);
        if (iter == jobs.node_to_job.end()) {
            std::string msg = "HaplotypePartitioner::partition_haplotypes(): cannot assign node " + std::to_string(node_id) + " to a job";
            throw std::runtime_error(msg);
        }
        result.jobs_for_cached_paths.push_back(iter->second);
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
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < chains_by_job.size(); job++) {
        const std::vector<gbwtgraph::TopLevelChain>& chains = chains_by_job[job];
        size_t total_subchains = 0, total_kmers = 0;
        for (auto& chain : chains) {
            try {
                this->build_subchains(chain, result.chains[chain.offset], parameters);
            } catch (const std::runtime_error& e) {
                std::cerr << "error: [job " << job << "]: " << e.what() << std::endl;
                std::exit(EXIT_FAILURE);
            }
            total_subchains += result.chains[chain.offset].subchains.size();
            for (auto& subchain : result.chains[chain.offset].subchains) {
                total_kmers += subchain.kmers.size();
            }
        }
        #pragma omp critical
        {
            result.header.total_subchains += total_subchains;
            result.header.total_kmers += total_kmers;
            if (this->verbosity >= Haplotypes::verbosity_detailed) {
                std::cerr << "Finished job " << job << " with " << chains.size() << " chains, " << total_subchains << " subchains, and " << total_kmers << " kmers" << std::endl;
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

std::vector<HaplotypePartitioner::Subchain>
HaplotypePartitioner::get_subchains(const gbwtgraph::TopLevelChain& chain, const Parameters& parameters) const {
    std::vector<Subchain> result;

    // First pass: take all snarls as subchains.
    std::vector<Subchain> snarls;
    handle_t snarl_start = empty_gbwtgraph_handle();
    bool has_start = false;
    bool was_snarl = false;
    net_handle_t curr = this->distance_index.get_bound(chain.chain, false, true);
    net_handle_t chain_end = this->distance_index.get_bound(chain.chain, true, false);
    while (curr != chain_end) {
        if (this->distance_index.is_node(curr)) {
            handle_t handle = this->distance_index.get_handle(curr, &this->gbz.graph);
            if (was_snarl) {
                if (!has_start) {
                    // If the chain starts with a snarl, we take it as a prefix.
                    snarls.push_back({ Haplotypes::Subchain::prefix, empty_gbwtgraph_handle(), handle });
                } else {
                    size_t distance = this->get_distance(snarl_start, handle);
                    if (distance < std::numeric_limits<size_t>::max()) {
                        // Normal snarl with two boundary nodes.
                        snarls.push_back({ Haplotypes::Subchain::normal, snarl_start, handle });
                    } else {
                        // The snarl is not connected, so we break it into two.
                        snarls.push_back({ Haplotypes::Subchain::suffix, snarl_start, empty_gbwtgraph_handle() });
                        snarls.push_back({ Haplotypes::Subchain::prefix, empty_gbwtgraph_handle(), handle });
                    }
                }
            }
            snarl_start = handle;
            has_start = true;
            was_snarl = false;
        } else if (this->distance_index.is_snarl(curr)) {
            was_snarl = true;
        }
        net_handle_t next;
        size_t successors = 0;
        this->distance_index.follow_net_edges(curr, &this->gbz.graph, false, [&](const net_handle_t& child) {
            successors++;
            next = child;
        });
        if (successors != 1) {
            throw std::runtime_error("HaplotypePartitioner::get_subchains(): chain " + std::to_string(chain.offset) + " has " + std::to_string(successors) + " successors for a child");
        }
        curr = next;
    }
    if (was_snarl && has_start) {
        // If the chain ends with a snarl, we take it as a suffix.
        snarls.push_back({ Haplotypes::Subchain::suffix, snarl_start, empty_gbwtgraph_handle() });
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
        while (tail + 1 < snarls.size()) {
            if (snarls[tail + 1].type != Haplotypes::Subchain::normal) {
                break;
            }
            size_t candidate = this->get_distance(snarls[head].start, snarls[tail + 1].end);
            if (candidate <= parameters.subchain_length) {
                tail++;
            } else {
                break;
            }
        }
        result.push_back({ Haplotypes::Subchain::normal, snarls[head].start, snarls[tail].end });
        head = tail + 1;
    }

    return result;
}

//------------------------------------------------------------------------------

std::vector<HaplotypePartitioner::sequence_type> HaplotypePartitioner::get_sequence_visits(handle_t handle) const {
    std::vector<gbwt::size_type> sa = this->r_index.decompressSA(gbwtgraph::GBWTGraph::handle_to_node(handle));
    std::vector<sequence_type> result;
    result.reserve(sa.size());
    for (size_t i = 0; i < sa.size(); i++) {
        result.push_back({ sa[i], i });
    }
    std::sort(result.begin(), result.end(), [&](sequence_type a, sequence_type b) -> bool {
        gbwt::size_type a_id = r_index.seqId(a.first);
        gbwt::size_type a_offset = r_index.seqOffset(a.first);
        gbwt::size_type b_id = r_index.seqId(b.first);
        gbwt::size_type b_offset = r_index.seqOffset(b.first);
        return ((a_id < b_id) || ((a_id == b_id) && (a_offset > b_offset)));
    });
    return result;
}

void sa_to_da(std::vector<HaplotypePartitioner::sequence_type>& sequences, const gbwt::FastLocate& r_index) {
    for (auto& sequence : sequences) {
        sequence.first = r_index.seqId(sequence.first);
    }
}

std::vector<HaplotypePartitioner::sequence_type> HaplotypePartitioner::get_sequences(handle_t handle) const {
    auto result = this->get_sequence_visits(handle);
    sa_to_da(result, this->r_index);
    return result;
}

std::vector<HaplotypePartitioner::sequence_type> HaplotypePartitioner::get_sequences(Subchain subchain) const {
    if (subchain.type == Haplotypes::Subchain::prefix) {
        return this->get_sequences(subchain.end);
    }
    if (subchain.type == Haplotypes::Subchain::suffix) {
        return this->get_sequences(subchain.start);
    }
    auto from = this->get_sequence_visits(subchain.start);
    auto to = this->get_sequence_visits(subchain.end);

    auto from_iter = from.begin();
    auto to_iter = to.begin();
    std::vector<sequence_type> result;
    while (from_iter != from.end() && to_iter != to.end()) {
        gbwt::size_type from_id = this->r_index.seqId(from_iter->first);
        gbwt::size_type to_id = this->r_index.seqId(to_iter->first);
        if (from_id == to_id) {
            // If a haplotype crosses the subchain multiple times, we take the last entry before
            // each exit.
            gbwt::size_type to_offset = this->r_index.seqOffset(to_iter->first);
            if (this->r_index.seqOffset(from_iter->first) >= to_offset) {
                auto peek = from_iter +1;
                while (peek != from.end() && this->r_index.seqId(peek->first) == from_id && this->r_index.seqOffset(peek->first) >= to_offset) {
                    from_iter = peek;
                    ++peek;
                }
                result.push_back(*from_iter);
                ++from_iter; ++to_iter;
            } else {
                ++to_iter;
            }
        } else if (from_id < to_id) {
            ++from_iter;
        } else if (from_id > to_id) {
            ++to_iter;
        }
    }

    sa_to_da(result, this->r_index);
    return result;
}

//------------------------------------------------------------------------------

// Generate a haplotype over the closed range from `pos` to `end`.
// Take at most start_max and end_max characters from the initial and the final
// node, respectively
// Returns an empty haplotype if there is only one node.
// Set `end = empty_gbwtgraph_handle()` to continue until the end without a final node.
std::string generate_haplotype(gbwt::edge_type pos, handle_t end, size_t start_max, size_t end_max, const gbwtgraph::GBWTGraph& graph) {
    std::string haplotype;
    if (pos == gbwt::invalid_edge() || pos.first == gbwt::ENDMARKER) {
        return haplotype;
    }

    // Handle the initial node.
    handle_t curr = gbwtgraph::GBWTGraph::node_to_handle(pos.first);
    if (curr == end) {
        return haplotype;
    }
    gbwtgraph::view_type view = graph.get_sequence_view(curr);
    size_t offset = (view.second > start_max ? view.second - start_max : 0);
    haplotype.append(view.first + offset, view.second - offset);

    while (true) {
        pos = graph.index->LF(pos);
        if (pos.first == gbwt::ENDMARKER) {
            break;
        }
        curr = gbwtgraph::GBWTGraph::node_to_handle(pos.first);
        view = graph.get_sequence_view(curr);
        if (curr == end) {
            haplotype.append(view.first, std::min(view.second, end_max));
            break;
        } else {
            haplotype.append(view.first, view.second);
        }
    }

    return haplotype;
}

// Return the sorted set of kmers that are minimizers in the sequence and have a single
// occurrence in the graph.
std::vector<HaplotypePartitioner::kmer_type> take_unique_minimizers(const std::string& sequence, const HaplotypePartitioner::minimizer_index_type& minimizer_index) {
    std::vector<HaplotypePartitioner::kmer_type> result;
    auto minimizers = minimizer_index.minimizers(sequence);
    result.reserve(minimizers.size());
    for (auto& minimizer : minimizers) {
        if (minimizer_index.count(minimizer) == 1) {
            result.push_back(minimizer.key.get_key());
        }
    }
    gbwt::removeDuplicates(result, false);
    return result;
}

std::vector<HaplotypePartitioner::kmer_type> HaplotypePartitioner::unique_minimizers(gbwt::size_type sequence_id) const {
    gbwt::edge_type pos = this->gbz.index.start(sequence_id);
    size_t limit = std::numeric_limits<size_t>::max();
    std::string haplotype = generate_haplotype(pos, empty_gbwtgraph_handle(), limit, limit, this->gbz.graph);
    return take_unique_minimizers(haplotype, this->minimizer_index);
}

std::vector<HaplotypePartitioner::kmer_type> HaplotypePartitioner::unique_minimizers(sequence_type sequence, Subchain subchain) const {
    gbwt::edge_type pos;
    size_t start_max = std::numeric_limits<size_t>::max(), end_max = this->minimizer_index.k() - 1;
    if (subchain.has_start()) {
        pos = gbwt::edge_type(gbwtgraph::GBWTGraph::handle_to_node(subchain.start), sequence.second);
        start_max = this->minimizer_index.k() - 1;
    } else {
        pos = this->gbz.index.start(sequence.first);
    }
    std::string haplotype = generate_haplotype(pos, subchain.end, start_max, end_max, this->gbz.graph);
    return take_unique_minimizers(haplotype, this->minimizer_index);
}

//------------------------------------------------------------------------------

// TODO: Return the number of skipped non-informative kmers.
/*
  Take a set of sequences defined by the sorted set of kmers that are present.
  Output a sorted vector of all kmers and the concatenated bitvectors that mark
  the presence of kmers in the sequences. The output vectors are assumed to be
  empty.
*/
void present_kmers(const std::vector<std::vector<HaplotypePartitioner::kmer_type>>& sequences,
    std::vector<std::pair<HaplotypePartitioner::kmer_type, size_t>>& all_kmers,
    sdsl::bit_vector& kmers_present) {

    // Build a map of distinct kmers. For each kmer, record the largest sequence
    // id containing the kmer and the number of sequences containing it.
    std::map<HaplotypePartitioner::kmer_type, std::pair<size_t, size_t>> present;
    for (size_t sequence_id = 0; sequence_id < sequences.size(); sequence_id++) {
        auto& sequence = sequences[sequence_id];
        for (auto kmer : sequence) {
            auto iter = present.find(kmer);
            if (iter != present.end()) {
                if (iter->second.first < sequence_id) {
                    iter->second.first = sequence_id;
                    iter->second.second++;
                }
            } else {
                present[kmer] = std::pair<size_t, size_t>(sequence_id, 1);
            }
        }
    }

    // Now take those kmers that occur in some but not in all sequences.
    // Use the first field for storing the offset of the kmer in the vector.
    all_kmers.reserve(present.size());
    size_t offset = 0;
    for (auto iter = present.begin(); iter != present.end(); ++iter) {
        if (iter->second.second < sequences.size()) {
            all_kmers.push_back({ iter->first, iter->second.second });
            iter->second.first = offset;
            offset++;
        }
    }

    // Transform the sequences into kmer presence bitvectors.
    kmers_present = sdsl::bit_vector(sequences.size() * all_kmers.size());
    for (size_t i = 0; i < sequences.size(); i++) {
        size_t start = i * all_kmers.size();
        for (auto kmer : sequences[i]) {
            auto iter = present.find(kmer);
            if (iter->second.second < sequences.size()) {
                kmers_present[start + iter->second.first] = 1;
            }
        }
    }
}

void HaplotypePartitioner::build_subchains(const gbwtgraph::TopLevelChain& chain, Haplotypes::TopLevelChain& output, const Parameters& parameters) const {
    std::vector<Subchain> subchains = this->get_subchains(chain, parameters);
    for (const Subchain& subchain : subchains) {
        std::vector<std::pair<Subchain, std::vector<sequence_type>>> to_process;
        auto sequences = this->get_sequences(subchain);
        if (sequences.empty()) {
            // There are no haplotypes crossing the subchain, so we break it into
            // a suffix and a prefix.
            to_process.push_back({ { Haplotypes::Subchain::suffix, subchain.start, empty_gbwtgraph_handle() }, this->get_sequences(subchain.start) });
            to_process.push_back({ { Haplotypes::Subchain::prefix, empty_gbwtgraph_handle(), subchain.end }, this->get_sequences(subchain.end) });
        } else {
            to_process.push_back({ subchain, std::move(sequences) });
        }
        for (auto iter = to_process.begin(); iter != to_process.end(); ++iter) {
            output.subchains.push_back({
                iter->first.type,
                gbwtgraph::GBWTGraph::handle_to_node(iter->first.start), gbwtgraph::GBWTGraph::handle_to_node(iter->first.end),
                {}, {}, sdsl::bit_vector()
            });
            Haplotypes::Subchain& subchain = output.subchains.back();
            std::vector<std::vector<kmer_type>> kmers_by_sequence;
            kmers_by_sequence.reserve(iter->second.size());
            for (sequence_type sequence : iter->second) {
                kmers_by_sequence.emplace_back(this->unique_minimizers(sequence, iter->first));
            }
            present_kmers(kmers_by_sequence, subchain.kmers, subchain.kmers_present);
            subchain.sequences = std::move(iter->second);
        }
    }

    // Take entire sequences if we could not generate any haplotypes.
    // Note that the kmer sets should be empty, as the sequences should
    // be identical.
    if (subchains.empty()) {
        output.subchains.push_back({
            Haplotypes::Subchain::full_haplotype,
            gbwt::ENDMARKER, gbwt::ENDMARKER,
            {}, {}, sdsl::bit_vector()
        });
        Haplotypes::Subchain& subchain = output.subchains.back();
        gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(chain.handle);
        auto sequences = this->r_index.decompressDA(node);
        std::vector<std::vector<kmer_type>> kmers_by_sequence;
        kmers_by_sequence.reserve(sequences.size());
        for (auto seq_id : sequences) {
            kmers_by_sequence.emplace_back(this->unique_minimizers(seq_id));
        }
        present_kmers(kmers_by_sequence, subchain.kmers, subchain.kmers_present);
        subchain.sequences.reserve(sequences.size());
        for (size_t i = 0; i < sequences.size(); i++) {
            subchain.sequences.push_back({ sequences[i], 0 });
        }
    }
}

//------------------------------------------------------------------------------

/*
 * A haplotype beging generated as a GBWT path.
 *
 * GBWT metadata will be set as following:
 *
 * * Sample name is "recombination".
 * * Contig name is "chain_X", where X is the chain identifier.
 * * Haplotype identifier is set during construction.
 * * Fragment identifier is set as necessary.
 */
struct RecombinatorHaplotype {
    typedef Recombinator::sequence_type sequence_type;

    // Contig identifier in GBWT metadata.
    // Offset of the top-level chain in the children of the root snarl.
    size_t chain;

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

    /*
     * Extends the haplotype over the given subchain by using the given
     * original haplotype.
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
     */
    void extend(sequence_type sequence, const Haplotypes::Subchain& subchain, const Recombinator& recombinator, gbwt::GBWTBuilder& builder);

    // Takes an existing haplotype from the GBWT index and inserts it into
    /// the builder. This is intended for fragments that do not contain
    /// subchains crossed by the original haplotypes. The call will fail if
    /// `extend()` has been called.
    void take(gbwt::size_type sequence_id, const Recombinator& recombinator, gbwt::GBWTBuilder& builder);

    // Extends the original haplotype from the latest `extend()` call until
    // the end, inserts it into the builder, and starts a new fragment.
    // The call will fail if `extend()` has not been called for this
    // fragment.
    void finish(const Recombinator& recombinator, gbwt::GBWTBuilder& builder);

private:
    // Extends the haplotype over a unary path from a previous subchain.
    void connect(gbwt::node_type until, const gbwtgraph::GBWTGraph& graph);

    // Takes a prefix of a sequence.
    void prefix(gbwt::size_type sequence_id, gbwt::node_type until, const gbwt::GBWT& index);

    // Extends the haplotype from the previous subchain until the end.
    void suffix(const gbwt::GBWT& index);

    // Inserts the current fragment into the builder.
    void insert(gbwt::GBWTBuilder& builder);
};

void RecombinatorHaplotype::extend(sequence_type sequence, const Haplotypes::Subchain& subchain, const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
    if (subchain.type == Haplotypes::Subchain::full_haplotype) {
        throw std::runtime_error("Haplotype::extend(): cannot extend a full haplotype");
    }

    if (subchain.type == Haplotypes::Subchain::prefix) {
        if (!this->path.empty())  {
            throw std::runtime_error("Haplotype::extend(): got a prefix subchain after the start of a fragment");
        }
        this->prefix(sequence.first, subchain.end, recombinator.gbz.index);
        return;
    }

    // Suffixes and normal subchains have a start node, so we must reach it first.
    if (!this->path.empty()) {
        this->connect(subchain.start, recombinator.gbz.graph);
    } else {
        this->prefix(sequence.first, subchain.start, recombinator.gbz.index);
    }

    gbwt::edge_type curr(subchain.start, sequence.second);
    if (!recombinator.gbz.index.contains(curr)) {
        throw std::runtime_error("Haplotype::extend(): the GBWT index does not contain position (" + std::to_string(curr.first) + ", " + std::to_string(curr.second) + ")");
    }

    if (subchain.type == Haplotypes::Subchain::suffix) {
        this->position = curr;
        this->finish(recombinator, builder);
        return;
    }

    // This is a normal subchain.
    while (curr.first != subchain.end) {
        curr = recombinator.gbz.index.LF(curr);
        if (curr.first == gbwt::ENDMARKER) {
            throw std::runtime_error("Haplotype::extend(): the sequence did not reach the end of the subchain at GBWT node " + std::to_string(subchain.end));
        }
        this->path.push_back(curr.first);
    }
    this->sequence_id = sequence.first;
    this->position = curr;
}

void RecombinatorHaplotype::take(gbwt::size_type sequence_id, const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
    if (!this->path.empty()) {
        throw std::runtime_error("Haplotype::take(): the current fragment is not empty");
    }
    if (sequence_id >= recombinator.gbz.index.sequences()) {
        throw std::runtime_error("Haplotype::take(): the GBWT index does not contain sequence " + std::to_string(sequence_id));
    }
    this->path = recombinator.gbz.index.extract(sequence_id);
    this->insert(builder);
    this->fragment++;
    this->sequence_id = gbwt::invalid_sequence();
    this->position = gbwt::invalid_edge();
    this->path.clear();
}

void RecombinatorHaplotype::finish(const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
    if (this->position == gbwt::invalid_edge()) {
        throw std::runtime_error("Haplotype::finish(): there is no current position");
    }
    this->suffix(recombinator.gbz.index);
    this->insert(builder);
    this->fragment++;
    this->sequence_id = gbwt::invalid_sequence();
    this->position = gbwt::invalid_edge();
    this->path.clear();
}

void RecombinatorHaplotype::connect(gbwt::node_type until, const gbwtgraph::GBWTGraph& graph) {
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

void RecombinatorHaplotype::prefix(gbwt::size_type sequence_id, gbwt::node_type until, const gbwt::GBWT& index) {
    this->position = gbwt::invalid_edge();
    if (sequence_id >= index.sequences()) {
        throw std::runtime_error("Haplotype::prefix(): invalid GBWT sequence id " + std::to_string(sequence_id));
    }
    this->sequence_id = sequence_id;
    for (gbwt::edge_type curr = index.start(sequence_id); curr.first != gbwt::ENDMARKER; curr = index.LF(curr)) {
        this->path.push_back(curr.first);
        if (curr.first == until) {
            this->position = curr;
            return;
        }
    }
    throw std::runtime_error("Haplotype::prefix(): GBWT sequence " + std::to_string(sequence_id) + " did not reach GBWT node " + std::to_string(until));
}

void RecombinatorHaplotype::suffix(const gbwt::GBWT& index) {
    for (gbwt::edge_type curr = index.LF(this->position); curr.first != gbwt::ENDMARKER; curr = index.LF(curr)) {
        this->path.push_back(curr.first);
    }
}

void RecombinatorHaplotype::insert(gbwt::GBWTBuilder& builder) {
    std::string sample_name = "recombination";
    gbwt::size_type sample_id = builder.index.metadata.sample(sample_name);
    if (sample_id >= builder.index.metadata.samples()) {
        builder.index.metadata.addSamples({ sample_name });
    }

    std::string contig_name = "chain_" + std::to_string(this->chain);
    gbwt::size_type contig_id = builder.index.metadata.contig(contig_name);
    if (contig_id >= builder.index.metadata.contigs()) {
        builder.index.metadata.addContigs({ contig_name });
    }

    builder.index.metadata.addPath(sample_id, contig_id, this->id, this->fragment);
    builder.insert(this->path, true);
}

//------------------------------------------------------------------------------

void Recombinator::Statistics::combine(const Statistics& another) {
    this->chains += another.chains;
    this->subchains += another.subchains;
    this->fragments += another.fragments;
    this->full_haplotypes += another.full_haplotypes;
    this->haplotypes = std::max(this->haplotypes, another.haplotypes);
    this->connections += another.connections;
    this->ref_paths += another.ref_paths;
    this->kmers += another.kmers;
    this->score += another.score;
}

std::ostream& Recombinator::Statistics::print(std::ostream& out) const {
    out << this->haplotypes << " haplotypes for " << this->chains << " chains ("
        << this->full_haplotypes << " full, " << this->subchains << " subchains, " << this->fragments << " fragments)";
    if (this->subchains > 0) {
        double connection_rate = static_cast<double>(this->connections) / (this->subchains * this->haplotypes);
        out << "; connection rate " << connection_rate;
    }
    if (this->ref_paths > 0) {
        out << "; included " << this->ref_paths << " reference paths";
    }
    if (this->kmers > 0) {
        double average_score = this->score / (this->kmers * this->haplotypes);
        out << "; used " << this->kmers << " kmers with average score " << average_score;
    }
    return out;
}

//------------------------------------------------------------------------------

Recombinator::Recombinator(const gbwtgraph::GBZ& gbz, Verbosity verbosity) :
    gbz(gbz), verbosity(verbosity)
{
}

//------------------------------------------------------------------------------

void add_path(const gbwt::GBWT& source, gbwt::size_type path_id, gbwt::GBWTBuilder& builder) {
    // We know that sufficient metadata exists, because this is a cached path.
    gbwt::PathName path_name = source.metadata.path(path_id);

    std::string sample_name = source.metadata.sample(path_name.sample);
    path_name.sample = builder.index.metadata.sample(sample_name);
    if (path_name.sample >= builder.index.metadata.samples()) {
        builder.index.metadata.addSamples({ sample_name });
    }

    std::string contig_name = source.metadata.contig(path_name.contig);
    path_name.contig = builder.index.metadata.contig(contig_name);
    if (path_name.contig >= builder.index.metadata.contigs()) {
        builder.index.metadata.addContigs({ contig_name });
    }

    builder.index.metadata.addPath(path_name);

    gbwt::vector_type path = source.extract(gbwt::Path::encode(path_id, false));
    builder.insert(path, true);
}

//------------------------------------------------------------------------------

void recombinator_sanity_checks(const Recombinator::Parameters& parameters) {
    if (parameters.num_haplotypes == 0) {
        std::string msg = "recombinator_sanity_checks(): number of haplotypes cannot be 0";
        throw std::runtime_error(msg);
    }
    if (parameters.present_discount < 0.0 || parameters.present_discount > 1.0) {
        std::string msg = "recombinator_sanity_checks(): present discount must be between 0.0 and 1.0";
        throw std::runtime_error(msg);
    }
    if (parameters.het_adjustment < 0.0) {
        std::string msg = "recombinator_sanity_checks(): het adjustment must be non-negative";
        throw std::runtime_error(msg);
    }
    if (parameters.absent_score < 0.0) {
        std::string msg = "recombinator_sanity_checks(): absent score must be non-negative";
        throw std::runtime_error(msg);
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
    for (auto iter = counts.begin(); iter != counts.end(); ++iter) {
        // We are only interested in kmers with multiple occurrences, as unique
        // kmers are likely sequencing errors.
        if (iter->second > 1) {
            count_to_frequency[iter->second]++;
        }
    }
    auto statistics = summary_statistics(count_to_frequency);
    double coverage = std::max(statistics.median, statistics.mode);
    if (verbosity >= Haplotypes::verbosity_detailed) {
        std::cerr << "Coverage: median " << statistics.median
            << ", mean " << statistics.mean
            << ", stdev " << statistics.stdev
            << ", mode " << statistics.mode
            << "; using " << coverage << std::endl;
    }
    if (verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Estimated kmer coverage in " << seconds << " seconds" << std::endl;
    }
    return coverage;
}

gbwt::GBWT Recombinator::generate_haplotypes(const Haplotypes& haplotypes, const std::string& kff_file, const Parameters& parameters) const {

    // Sanity checks (may throw).
    recombinator_sanity_checks(parameters);

    // Get kmer counts (may throw) and determine coverage.
    hash_map<Haplotypes::Subchain::kmer_type, size_t> counts = haplotypes.kmer_counts(kff_file, this->verbosity);
    double coverage = get_or_estimate_coverage(counts, parameters, this->verbosity);

    double start = gbwt::readTimer();
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        if (parameters.random_sampling) {
            std::cerr << "Building GBWT (random sampling)" << std::endl;
        } else {
            std::cerr << "Building GBWT" << std::endl;
        }
    }

    // Determine construction jobs.
    std::vector<std::vector<size_t>> jobs(haplotypes.jobs());
    for (auto& chain : haplotypes.chains) {
        if (chain.job_id < haplotypes.jobs()) {
            jobs[chain.job_id].push_back(chain.offset);
        }
    }

    // Figure out GBWT path ids for reference paths in each job.
    std::vector<std::vector<gbwt::size_type>> reference_paths(haplotypes.jobs());
    if (parameters.include_reference) {
        for (size_t i = 0; i < this->gbz.graph.named_paths.size(); i++) {
            size_t job_id = haplotypes.jobs_for_cached_paths[i];
            if (job_id < haplotypes.jobs()) {
                reference_paths[job_id].push_back(this->gbz.graph.named_paths[i].id);
            }
        }
    }

    // Build partial indexes.
    double checkpoint = gbwt::readTimer();
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Running " << omp_get_max_threads() << " GBWT construction jobs in parallel" << std::endl;
    }
    std::vector<gbwt::GBWT> indexes(jobs.size());
    Statistics statistics;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < jobs.size(); job++) {
        gbwt::GBWTBuilder builder(sdsl::bits::length(this->gbz.index.sigma() - 1), parameters.buffer_size);
        builder.index.addMetadata();
        Statistics job_statistics;
        // Add haplotypes for each chain.
        for (auto chain_id : jobs[job]) {
            try {
                Statistics chain_statistics = this->generate_haplotypes(haplotypes.chains[chain_id], counts, builder, parameters, coverage);
                job_statistics.combine(chain_statistics);
            } catch (const std::runtime_error& e) {
                std::cerr << "error: [job " << job << "]: " << e.what() << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        // Add named and reference paths.
        for (auto path_id : reference_paths[job]) {
            add_path(this->gbz.index, path_id, builder);
            job_statistics.ref_paths++;
        }
        builder.finish();
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

    // Merge the partial indexes.
    checkpoint = gbwt::readTimer();
    if (this->verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Merging the partial indexes" << std::endl;
    }
    gbwt::GBWT merged(indexes);
    if (parameters.include_reference) {
        // If we included reference paths, set the same samples as references in the output GBWT.
        std::string reference_samples = this->gbz.index.tags.get(gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG);
        if (!reference_samples.empty()) {
            merged.tags.set(gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG, reference_samples);
        }
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

Recombinator::Statistics Recombinator::generate_haplotypes(const Haplotypes::TopLevelChain& chain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    gbwt::GBWTBuilder& builder,
    const Parameters& parameters,
    double coverage) const {

    // TODO: What are the proper thresholds?
    enum kmer_presence { absent, heterozygous, present, ignore };
    double absent_threshold = coverage * 0.1;
    double heterozygous_threshold = coverage / std::log(4.0);
    double homozygous_threshold = coverage * 2.5;

    std::vector<RecombinatorHaplotype> haplotypes;
    for (size_t i = 0; i < parameters.num_haplotypes; i++) {
        haplotypes.push_back({ chain.offset, i, 0, gbwt::invalid_sequence(), gbwt::invalid_edge(), {} });
    }

    // TODO: Random seed.
    std::mt19937_64 rng(0xACDCABBA);

    Statistics statistics;
    statistics.chains = 1; statistics.haplotypes = parameters.num_haplotypes;
    if (chain.subchains.size() == 1 && chain.subchains.front().type == Haplotypes::Subchain::full_haplotype) {
        // Full haplotypes should all be identical, because there are no snarls.
        // Therefore we do not need kmers.
        auto& subchain = chain.subchains.front();
        for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
            assert(!subchain.sequences.empty());
            size_t seq = haplotype % subchain.sequences.size();
            haplotypes[haplotype].take(subchain.sequences[seq].first, *this, builder);
        }
        statistics.full_haplotypes = 1;
    } else {
        bool have_haplotypes = false;
        for (auto& subchain : chain.subchains) {
            if (subchain.type == Haplotypes::Subchain::full_haplotype) {
                throw std::runtime_error("Recombinator::generate_haplotypes(): nontrivial chain " + std::to_string(chain.offset) + " contains a subchain with full haplotypes");
            }
            assert(!subchain.sequences.empty());

            // TODO: -log prob may be the right score once we have enough haplotypes, but
            // right now +1 works better, because we don't have haplotypes with the right
            // combination of rare kmers.
            // Determine the type of each kmer in the sample and the score for the kmer.
            // A haplotype with the kmer gets +1.0 * score, while a haplotype without it
            // gets -1.0 * score.
            std::vector<std::pair<kmer_presence, double>> kmer_types;
            size_t selected_kmers = 0;
            for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                double count = kmer_counts.at(subchain.kmers[kmer_id].first);
                if (count < absent_threshold) {
                    kmer_types.push_back({ absent, -1.0 * parameters.absent_score });
                    selected_kmers++;
                } else if (count < heterozygous_threshold) {
                    kmer_types.push_back({ heterozygous, 0.0 });
                    selected_kmers++;
                } else if (count < homozygous_threshold) {
                    kmer_types.push_back({ present, 1.0 });
                    selected_kmers++;
                } else {
                    kmer_types.push_back({ ignore, 0.0 });
                }
            }
            statistics.kmers += selected_kmers;

            // Select the haplotypes greedily.
            std::vector<std::pair<size_t, double>> selected_haplotypes;
            std::vector<std::pair<size_t, double>> remaining_haplotypes;
            for (size_t seq_offset = 0; seq_offset < subchain.sequences.size(); seq_offset++) {
                remaining_haplotypes.push_back( { seq_offset, 0.0 });
            }
            while (selected_haplotypes.size() < haplotypes.size() && !remaining_haplotypes.empty()) {
                // Score the remaining haplotypes.
                for (size_t i = 0; i < remaining_haplotypes.size(); i++) {
                    size_t offset = remaining_haplotypes[i].first * subchain.kmers.size();
                    double score = 0.0;
                    for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                        double multiplier = -1.0 + 2.0 * subchain.kmers_present[offset + kmer_id];
                        score += multiplier * kmer_types[kmer_id].second;
                    }
                    remaining_haplotypes[i].second = score;
                }

                // Select the haplotype to use.
                size_t selected = 0;
                if (parameters.random_sampling) {
                    selected = rng() % remaining_haplotypes.size();
                } else {
                    for (size_t i = 1; i < remaining_haplotypes.size(); i++) {
                        if (remaining_haplotypes[i].second > remaining_haplotypes[selected].second) {
                            selected = i;
                        }
                    }
                }
                selected_haplotypes.push_back(remaining_haplotypes[selected]);
                remaining_haplotypes.erase(remaining_haplotypes.begin() + selected);

                // Adjust kmer scores based on the selected haplotype.
                size_t offset = selected_haplotypes.back().first * subchain.kmers.size();
                for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                    switch (kmer_types[kmer_id].first) {
                    case heterozygous:
                        kmer_types[kmer_id].second += (subchain.kmers_present[offset + kmer_id] ? -1.0 : 1.0) * parameters.het_adjustment;
                        break;
                    case present:
                        if (subchain.kmers_present[offset + kmer_id]) {
                            kmer_types[kmer_id].second *= parameters.present_discount;
                        }
                        break;
                    default:
                        break;
                    }
                }
            }

            // If we did not have enough haplotypes in the subchain, repeat them as necessary.
            size_t original_selected = selected_haplotypes.size();
            for (size_t i = original_selected; i < haplotypes.size(); i++) {
                selected_haplotypes.push_back(selected_haplotypes[i % original_selected]);
            }

            // Try to match the existing haplotypes with the selected sequences based on
            // GBWT sequence id.
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
            for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                size_t selected = haplotype_to_selected[haplotype];
                size_t seq_offset = selected_haplotypes[selected].first;
                statistics.score += selected_haplotypes[selected].second;
                haplotypes[haplotype].extend(subchain.sequences[seq_offset], subchain, *this, builder);
            }
            have_haplotypes = subchain.has_end();
            statistics.subchains++;
        }
        if (have_haplotypes) {
            for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                haplotypes[haplotype].finish(*this, builder);
            }
        }
        statistics.fragments = haplotypes.front().fragment;
    }

    return statistics;
}

//------------------------------------------------------------------------------

std::vector<Recombinator::LocalHaplotype> Recombinator::extract_sequences(const Haplotypes& haplotypes,
    const std::string& kff_file,
    size_t chain_id, size_t subchain_id) const {

    if (chain_id >= haplotypes.chains.size()) {
        std::string msg = "Recombinator::extract_sequences(): invalid chain id " + std::to_string(chain_id);
        throw std::runtime_error(msg);
    }
    if (subchain_id >= haplotypes.chains[chain_id].subchains.size()) {
        std::string msg = "Recombinator::extract_sequences(): invalid subchain id " + std::to_string(subchain_id) +
            " in chain " + std::to_string(chain_id);
        throw std::runtime_error(msg);
    }

    // Extract the haplotypes.
    const Haplotypes::Subchain& subchain = haplotypes.chains[chain_id].subchains[subchain_id];
    std::vector<LocalHaplotype> result(subchain.sequences.size());
    for (size_t i = 0; i < subchain.sequences.size(); i++) {
        size_t path_id = gbwt::Path::id(subchain.sequences[i].first);
        path_handle_t path_handle = this->gbz.graph.path_to_handle(path_id);
        result[i].name = this->gbz.graph.get_path_name(path_handle);

        gbwt::edge_type pos;
        if (subchain.has_start()) {
            pos = gbwt::edge_type(subchain.start, subchain.sequences[i].second);
        } else {
            pos = this->gbz.index.start(subchain.sequences[i].first);
        }
        handle_t until = gbwtgraph::GBWTGraph::node_to_handle(subchain.end);
        size_t limit = std::numeric_limits<size_t>::max();
        result[i].sequence = generate_haplotype(pos, until, limit, limit, this->gbz.graph);
    }

    // FIXME: Score the haplotypes.
    // 1. read kmer counts
    // 2. estimate coverage
    // 3. "Select the haplotypes greedily" loop

    return result;
}

//------------------------------------------------------------------------------

} // namespace vg

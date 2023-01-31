#include "recombinator.hpp"

#include "kff.hpp"

#include <algorithm>
#include <map>
#include <unordered_set>

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

//------------------------------------------------------------------------------

// Returns a GBWTGraph handle as a string (id, orientation).
std::string to_string(handle_t handle) {
    gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(handle);
    return std::string("(") + std::to_string(gbwt::Node::id(node)) + std::string(", ") + std::to_string(gbwt::Node::is_reverse(node)) + std::string(")");
}

//------------------------------------------------------------------------------

// Returns the big-endian interger representation of the kmer in the canonical orientation.
std::uint64_t kff_to_key(const std::uint8_t* kmer, size_t k, size_t bytes, const std::uint8_t* encoding) {
    std::vector<std::uint8_t> rc = kff_reverse_complement(kmer, k, encoding);
    return std::min(kff_parse(kmer, bytes), kff_parse(rc.data(), bytes));
}

std::unordered_map<Haplotypes::Subchain::kmer_type, size_t> Haplotypes::kmer_counts(const std::string& kff_file) const {
    // Open and validate the kmer count file.
    Kff_reader reader(kff_file);
    std::uint64_t kff_k = reader.get_var("k");
    if (kff_k != this->k()) {
        throw std::runtime_error("Haplotypes::kmer_counts(): expected " + std::to_string(this->k()) +
            "-mers but KFF file " + kff_file + " contains " + std::to_string(kff_k) + "-mers");
    }
    const std::uint8_t* encoding = reader.get_encoding();
    std::string decoding = kff_invert(encoding);
    size_t bytes = kff_bytes(this->k());
    size_t data_bytes = reader.get_var("data_size");

    // Populate the map with the kmers we are interested in.
    // FIXME: Faster hash table?
    std::unordered_map<Subchain::kmer_type, size_t> result;
    for (size_t chain_id = 0; chain_id < this->chains.size(); chain_id++) {
        const TopLevelChain& chain = this->chains[chain_id];
        for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
            const Subchain& subchain = chain.subchains[subchain_id];
            for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                result[subchain.kmers[kmer_id]] = 0;
            }
        }
    }

    // Add the counts.
    std::uint8_t* kmer;
    std::uint8_t* data;
    while (reader.has_next()) {
        reader.next_kmer(kmer, data);
        Subchain::kmer_type recoded = kff_recode(kmer, this->k(), decoding);
        auto iter = result.find(recoded);
        if (iter != result.end()) {
            iter->second += kff_parse(data, data_bytes);
        } else {
            Subchain::kmer_type rc = minimizer_reverse_complement(recoded, this->k());
            iter = result.find(rc);
            if (iter != result.end()) {
                iter->second += kff_parse(data, data_bytes);
            }
        }
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

    this->kmers = sdsl::simple_sds::load_vector<kmer_type>(in);
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

    this->chains.resize(this->header.top_level_chains);
    for (auto& chain : this->chains) {
        chain.simple_sds_load(in);
    }
}

size_t Haplotypes::simple_sds_size() const {
    size_t result = sdsl::simple_sds::value_size<Header>();
    for (auto& chain : this->chains) {
        result += chain.simple_sds_size();
    }
    return result;
}

//------------------------------------------------------------------------------

HaplotypePartitioner::HaplotypePartitioner(const gbwtgraph::GBZ& gbz,
    const gbwt::FastLocate& r_index,
    const SnarlDistanceIndex& distance_index,
    const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
    Verbosity verbosity) :
    gbz(gbz), r_index(r_index), distance_index(distance_index), minimizer_index(minimizer_index),
    verbosity(verbosity)
{
}

//------------------------------------------------------------------------------

Haplotypes HaplotypePartitioner::partition_haplotypes(const Parameters& parameters) const {
    // FIXME parameter sanity checks

    Haplotypes result;
    result.header.k = this->minimizer_index.k();

    // Determine top-level chains and assign them to jobs.
    double start = gbwt::readTimer();
    if (this->verbosity >= verbosity_basic) {
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
    jobs = gbwtgraph::ConstructionJobs(); // Save memory.
    if (this->verbosity >= verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Partitioned " << result.components() << " components into " << result.jobs() << " jobs in " << seconds << " seconds" << std::endl;
    }

    // Determine the subchains and sequences for each top-level chain.
    if (verbosity >= verbosity_basic) {
        std::cerr << "Running " << omp_get_max_threads() << " jobs in parallel" << std::endl;
    }
    start = gbwt::readTimer();
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < chains_by_job.size(); job++) {
        const std::vector<gbwtgraph::TopLevelChain>& chains = chains_by_job[job];
        for (auto& chain : chains) {
            try {
                this->build_subchains(chain, result.chains[chain.offset], parameters);
            } catch (const std::runtime_error& e) {
                std::cerr << "error: [job " << job << "]: " << e.what() << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        if (this->verbosity >= verbosity_detailed) {
            #pragma omp critical
            {
                std::cerr << "Finished job " << job << " with " << chains.size() << " chains" << std::endl;
            }
        }
    }
    if (verbosity >= verbosity_basic) {
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
            if (candidate <= SUBCHAIN_LENGTH) {
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
std::vector<HaplotypePartitioner::kmer_type> take_unique_minimizers(const std::string& sequence, const gbwtgraph::DefaultMinimizerIndex& minimizer_index) {
    std::vector<HaplotypePartitioner::kmer_type> result;
    auto minimizers = minimizer_index.minimizers(sequence);
    result.reserve(minimizers.size());
    for (auto& minimizer : minimizers) {
        if (minimizer_index.count(minimizer) == 1) {
            result.push_back(minimizer.key.get_key());
        }
    }
    std::sort(result.begin(), result.end());
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

/*
  Take a set of sequences defined by the sorted set of kmers that are present.
  Output a sorted vector of all kmers and the concatenated bitvectors that mark
  the presence of kmers in the sequences. The output vectors are assumed to be
  empty.
*/
void present_kmers(std::vector<std::vector<HaplotypePartitioner::kmer_type>>& sequences,
    std::vector<HaplotypePartitioner::kmer_type>& all_kmers,
    sdsl::bit_vector& kmers_present) {

    // Determine the distinct kmers and their ranks.
    std::map<HaplotypePartitioner::kmer_type, size_t> present;
    for (auto& sequence : sequences) {
        for (auto kmer : sequence) {
            present[kmer] = 0;
        }
    }
    all_kmers.reserve(present.size());
    size_t offset = 0;
    for (auto iter = present.begin(); iter != present.end(); ++iter) {
        all_kmers.push_back(iter->first);
        iter->second = offset;
        offset++;
    }

    // Transform the sequences into kmer presence bitvectors.
    kmers_present = sdsl::bit_vector(sequences.size() * present.size());
    for (size_t i = 0; i < sequences.size(); i++) {
        size_t start = i * present.size();
        for (auto kmer : sequences[i]) {
            kmers_present[start + present[kmer]] = 1;
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
            // TODO: It will eventually be faster to find all unique minimizers in the relevant subgraph,
            // align them to the haplotypes, and use the r-index to determine which haplotypes contain
            // that minimizer. On the other hand, that will cause some issues with haplotypes that visit
            // the same subchain multiple times.
            output.subchains.push_back({
                iter->first.type,
                gbwtgraph::GBWTGraph::handle_to_node(iter->first.start), gbwtgraph::GBWTGraph::handle_to_node(iter->first.end),
                {}, {}, sdsl::bit_vector()
            });
            Haplotypes::Subchain& subchain = output.subchains.back();
            std::vector<std::vector<kmer_type>> kmers_by_sequence;
            kmers_by_sequence.reserve(iter->second.size());
            for (sequence_type sequence : iter->second) {
                kmers_by_sequence.push_back(this->unique_minimizers(sequence, iter->first));
            }
            present_kmers(kmers_by_sequence, subchain.kmers, subchain.kmers_present);
            subchain.sequences = std::move(iter->second);
        }
    }

    // Take entire sequences if we could not generate any haplotypes.
    // TODO: Maybe we don't need kmers here, as the sequences should be identical.
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
            kmers_by_sequence.push_back(this->unique_minimizers(seq_id));
        }
        present_kmers(kmers_by_sequence, subchain.kmers, subchain.kmers_present);
        subchain.sequences.reserve(sequences.size());
        for (size_t i = 0; i < sequences.size(); i++) {
            subchain.sequences.push_back({ sequences[i], 0 });
        }
    }
}

//------------------------------------------------------------------------------

void Recombinator::Haplotype::extend(sequence_type sequence, const Haplotypes::Subchain& subchain, const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
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
    this->position = curr;
}

void Recombinator::Haplotype::take(gbwt::size_type sequence_id, const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
    if (!this->path.empty()) {
        throw std::runtime_error("Haplotype::take(): the current fragment is not empty");
    }
    if (sequence_id >= recombinator.gbz.index.sequences()) {
        throw std::runtime_error("Haplotype::take(): the GBWT index does not contain sequence " + std::to_string(sequence_id));
    }
    this->path = recombinator.gbz.index.extract(sequence_id);
    this->insert(builder);
    this->fragment++;
    this->position = gbwt::invalid_edge();
    this->path.clear();
}

void Recombinator::Haplotype::finish(const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
    if (this->position == gbwt::invalid_edge()) {
        throw std::runtime_error("Haplotype::finish(): there is no current position");
    }
    this->suffix(recombinator.gbz.index);
    this->insert(builder);
    this->fragment++;
    this->position = gbwt::invalid_edge();
    this->path.clear();
}

void Recombinator::Haplotype::connect(gbwt::node_type until, const gbwtgraph::GBWTGraph& graph) {
    handle_t curr = gbwtgraph::GBWTGraph::node_to_handle(this->position.first);
    handle_t end = gbwtgraph::GBWTGraph::node_to_handle(until);
    this->position = gbwt::invalid_edge();
    std::unordered_set<handle_t> visited;
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

void Recombinator::Haplotype::prefix(gbwt::size_type sequence_id, gbwt::node_type until, const gbwt::GBWT& index) {
    this->position = gbwt::invalid_edge();
    if (sequence_id >= index.sequences()) {
        throw std::runtime_error("Haplotype::prefix(): invalid GBWT sequence id " + std::to_string(sequence_id));
    }
    for (gbwt::edge_type curr = index.start(sequence_id); curr.first != gbwt::ENDMARKER; curr = index.LF(curr)) {
        this->path.push_back(curr.first);
        if (curr.first == until) {
            this->position = curr;
            return;
        }
    }
    throw std::runtime_error("Haplotype::prefix(): GBWT sequence " + std::to_string(sequence_id) + " did not reach GBWT node " + std::to_string(until));
}

void Recombinator::Haplotype::suffix(const gbwt::GBWT& index) {
    for (gbwt::edge_type curr = index.LF(this->position); curr.first != gbwt::ENDMARKER; curr = index.LF(curr)) {
        this->path.push_back(curr.first);
    }
    this->position = gbwt::invalid_edge();
}

void Recombinator::Haplotype::insert(gbwt::GBWTBuilder& builder) {
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
    this->kmers += another.kmers;
    this->score += another.score;
}

std::ostream& Recombinator::Statistics::print(std::ostream& out) const {
    size_t normal_chains = this->chains - this->full_haplotypes;
    out << normal_chains << " chains with " << this->subchains << " subchains and " << this->fragments << " fragments; "
        << this->full_haplotypes << " chains with full haplotypes; generated "
        << this->haplotypes << " haplotypes";
    if (this->kmers > 0) {
        double average_score = static_cast<double>(this->score) / (this->kmers * this->haplotypes);
        out << " with " << this->kmers << " kmers and average score " << average_score;
    }
    return out;
}

//------------------------------------------------------------------------------

Recombinator::Recombinator(const gbwtgraph::GBZ& gbz, HaplotypePartitioner::Verbosity verbosity) :
    gbz(gbz), verbosity(verbosity)
{
}

//------------------------------------------------------------------------------

gbwt::GBWT Recombinator::generate_haplotypes(const Haplotypes& haplotypes, const std::string& kff_file, const Parameters& parameters) const {

    double start = gbwt::readTimer();

    // FIXME sanity checks for parameters

    // Get kmer counts.
    double checkpoint = gbwt::readTimer();
    if (this->verbosity >= HaplotypePartitioner::verbosity_basic) {
        std::cerr << "Reading kmer counts" << std::endl;
    }
    std::unordered_map<Haplotypes::Subchain::kmer_type, size_t> counts;
    try {
        counts = haplotypes.kmer_counts(kff_file);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (this->verbosity >= HaplotypePartitioner::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Read the kmer counts in " << seconds << " seconds";
    }

    if (this->verbosity >= HaplotypePartitioner::verbosity_basic) {
        std::cerr << "Building GBWT" << std::endl;
    }

    // Determine construction jobs.
    std::vector<std::vector<size_t>> jobs(haplotypes.jobs());
    for (auto& chain : haplotypes.chains) {
        if (chain.job_id < haplotypes.jobs()) {
            jobs[chain.job_id].push_back(chain.offset);
        }
    }

    // Build partial indexes.
    checkpoint = gbwt::readTimer();
    if (this->verbosity >= HaplotypePartitioner::verbosity_basic) {
        std::cerr << "Running " << omp_get_max_threads() << " GBWT construction jobs in parallel" << std::endl;
    }
    std::vector<gbwt::GBWT> indexes(jobs.size());
    Statistics statistics;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < jobs.size(); job++) {
        gbwt::GBWTBuilder builder(sdsl::bits::length(this->gbz.index.sigma() - 1), parameters.buffer_size);
        builder.index.addMetadata();
        Statistics job_statistics;
        for (auto chain_id : jobs[job]) {
            try {
                Statistics chain_statistics = this->generate_haplotypes(haplotypes.chains[chain_id], counts, builder, parameters);
                job_statistics.combine(chain_statistics);
            } catch (const std::runtime_error& e) {
                std::cerr << "error: [job " << job << "]: " << e.what() << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        builder.finish();
        indexes[job] = builder.index;
        #pragma omp critical
        {
            if (this->verbosity >= HaplotypePartitioner::verbosity_detailed) {
                std::cerr << "Job " << job << ": "; job_statistics.print(std::cerr) << std::endl;
            }
            statistics.combine(job_statistics);
        }
    }
    if (this->verbosity >= HaplotypePartitioner::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Processed "; statistics.print(std::cerr) << std::endl;
        std::cerr << "Finished the jobs in " << seconds << " seconds" << std::endl;
    }

    // Merge the partial indexes.
    checkpoint = gbwt::readTimer();
    if (this->verbosity >= HaplotypePartitioner::verbosity_basic) {
        std::cerr << "Merging the partial indexes" << std::endl;
    }
    gbwt::GBWT merged(indexes);
    if (this->verbosity >= HaplotypePartitioner::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Merged the indexes in " << seconds << " seconds" << std::endl;
    }

    if (this->verbosity >= HaplotypePartitioner::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Built the GBWT in " << seconds << " seconds" << std::endl;
    }
    return merged;
}

Recombinator::Statistics Recombinator::generate_haplotypes(const Haplotypes::TopLevelChain& chain,
    const std::unordered_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    gbwt::GBWTBuilder& builder,
    const Parameters& parameters) const {

    // TODO: This is a placeholder.
    enum kmer_presence { absent, present, ignore };
    double absent_threshold = parameters.coverage * 0.2;
    double heterozygous_threshold = parameters.coverage * 0.7;
    double homozygous_threshold = parameters.coverage * 1.5;

    std::vector<Haplotype> haplotypes;
    for (size_t i = 0; i < parameters.num_haplotypes; i++) {
        haplotypes.push_back({ chain.offset, i, 0, gbwt::invalid_edge(), {} });
    }

    Statistics statistics;
    statistics.chains = 1; statistics.haplotypes = parameters.num_haplotypes;
    if (chain.subchains.size() == 1 && chain.subchains.front().type == Haplotypes::Subchain::full_haplotype) {
        // TODO: Full haplotypes should all be identical, because there are no snarls. Therefore we do not need kmers.
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

            // TODO: This is a placeholder.
            // Determine the type of each kmer in the sample.
            std::vector<kmer_presence> kmer_types;
            for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                double count = kmer_counts.at(subchain.kmers[kmer_id]);
                if (count < absent_threshold) {
                    kmer_types.push_back(absent);
                } else if (count < heterozygous_threshold) {
                    kmer_types.push_back(ignore);
                } else if (count < homozygous_threshold) {
                    kmer_types.push_back(present);
                } else {
                    kmer_types.push_back(ignore);
                }
            }

            // Score the sequences by kmer presence and sort them by score in descending order.
            std::vector<std::pair<size_t, std::int64_t>> sequence_scores;
            for (size_t sequence_id = 0; sequence_id < subchain.sequences.size(); sequence_id++) {
                size_t offset = sequence_id * subchain.kmers.size();
                std::int64_t score = 0.0;
                for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                    switch (kmer_types[kmer_id]) {
                    case present:
                        score += (subchain.kmers_present[offset + kmer_id] ? 1 : -1);
                        statistics.kmers++;
                        break;
                    case absent:
                        score += (subchain.kmers_present[offset + kmer_id] ? -1 : 1);
                        statistics.kmers++;
                        break;
                    default:
                        break;
                    }
                }
                sequence_scores.push_back({ sequence_id, score });
            }
            std::sort(sequence_scores.begin(), sequence_scores.end(), [](std::pair<size_t, std::int64_t> a, std::pair<size_t, std::int64_t> b) -> bool {
                return (a.second > b.second);
            });

            // Extend the haplotypes with the highest-scoring sequences.
            for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                size_t seq_id = sequence_scores[haplotype % sequence_scores.size()].first;
                haplotypes[haplotype].extend(subchain.sequences[seq_id], subchain, *this, builder);
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

} // namespace vg

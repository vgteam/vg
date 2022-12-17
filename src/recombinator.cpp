#include "recombinator.hpp"

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

//------------------------------------------------------------------------------

// An empty handle.
handle_t empty_handle() {
    return gbwtgraph::GBWTGraph::node_to_handle(0);
}

// Returns a GBWTGraph handle as a string (id, orientation).
std::string to_string(handle_t handle) {
    gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(handle);
    return std::string("(") + std::to_string(gbwt::Node::id(node)) + std::string(", ") + std::to_string(gbwt::Node::is_reverse(node)) + std::string(")");
}

//------------------------------------------------------------------------------

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

Haplotypes HaplotypePartitioner::partition_haplotypes() const {
    Haplotypes result;
    result.header.k = this->minimizer_index.k();

    // Determine top-level chains and assign them to jobs.
    double start = gbwt::readTimer();
    if (this->verbosity >= verbosity_basic) {
        std::cerr << "Determining construction jobs" << std::endl;
    }
    size_t size_bound = this->gbz.graph.get_node_count() / APPROXIMATE_JOBS;
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

    // FIXME number of construction jobs
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
                this->build_subchains(chain, result.chains[chain.offset]);
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
HaplotypePartitioner::get_subchains(const gbwtgraph::TopLevelChain& chain) const {
    std::vector<Subchain> result;

    // First pass: take all snarls as subchains.
    std::vector<Subchain> snarls;
    handle_t snarl_start = empty_handle();
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
                    snarls.push_back({ Haplotypes::Subchain::prefix, empty_handle(), handle });
                } else {
                    size_t distance = this->get_distance(snarl_start, handle);
                    if (distance < std::numeric_limits<size_t>::max()) {
                        // Normal snarl with two boundary nodes.
                        snarls.push_back({ Haplotypes::Subchain::normal, snarl_start, handle });
                    } else {
                        // The snarl is not connected, so we break it into two.
                        snarls.push_back({ Haplotypes::Subchain::suffix, snarl_start, empty_handle() });
                        snarls.push_back({ Haplotypes::Subchain::prefix, empty_handle(), handle });
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
        snarls.push_back({ Haplotypes::Subchain::suffix, snarl_start, empty_handle() });
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
            // To avoid weird situations, we only consider the first time a haplotype crosses this
            // subchain. We combine the first exit from the subchain and the last entry before it.
            // These situations should be rare.
            // TODO: Should we consider all minimal crossings?
            bool multiple_visits = false;
            gbwt::size_type to_offset = this->r_index.seqOffset(to_iter->first);
            if (this->r_index.seqOffset(from_iter->first) >= this->r_index.seqOffset(to_iter->first)) {
                auto peek = from_iter + 1;
                while (peek != from.end() && this->r_index.seqId(peek->first) == from_id && this->r_index.seqOffset(peek->first) >= to_offset) {
                    from_iter = peek;
                    ++peek;
                    multiple_visits = true;
                }
                result.push_back(*from_iter);
            }
            ++from_iter;
            while (from_iter != from.end() && this->r_index.seqId(from_iter->first) == from_id) {
                ++from_iter;
                multiple_visits = true;
            }
            ++to_iter;
            while (to_iter != to.end() && this->r_index.seqId(to_iter->first) == to_id) {
                ++to_iter;
                multiple_visits = true;
            }
            if (multiple_visits && this->verbosity >= verbosity_debug) {
                #pragma omp critical
                {
                    std::cerr << "Taking the first visit of sequence " << from_id << " to subchain " << to_string(subchain.start) << " to " << to_string(subchain.end) << std::endl;
                }
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

void append_to(std::string& haplotype, handle_t handle, const gbwtgraph::GBWTGraph& graph) {
    gbwtgraph::view_type view = graph.get_sequence_view(handle);
    haplotype.append(view.first, view.second);
}

// Generate a haplotype over the closed range from `pos` to `end`.
// Set `end = empty_handle()` to continue until the end.
std::string generate_haplotype(gbwt::edge_type pos, handle_t end, const gbwtgraph::GBWTGraph& graph) {
    std::string haplotype;
    if (pos == gbwt::invalid_edge() || pos.first == gbwt::ENDMARKER) {
        return haplotype;
    }

    handle_t curr = gbwtgraph::GBWTGraph::node_to_handle(pos.first);
    append_to(haplotype, curr, graph);
    while (curr != end) {
        pos = graph.index->LF(pos);
        if (pos.first == gbwt::ENDMARKER) {
            break;
        }
        curr = gbwtgraph::GBWTGraph::node_to_handle(pos.first);
        append_to(haplotype, curr, graph);
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
    std::string haplotype = generate_haplotype(pos, empty_handle(), this->gbz.graph);
    return take_unique_minimizers(haplotype, this->minimizer_index);
}

std::vector<HaplotypePartitioner::kmer_type> HaplotypePartitioner::unique_minimizers(sequence_type sequence, Subchain subchain) const {
    gbwt::edge_type pos;
    if (subchain.has_start()) {
        pos = gbwt::edge_type(gbwtgraph::GBWTGraph::handle_to_node(subchain.start), sequence.second);
    } else {
        pos = this->gbz.index.start(sequence.first);
    }
    std::string haplotype = generate_haplotype(pos, subchain.end, this->gbz.graph);
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

void HaplotypePartitioner::build_subchains(const gbwtgraph::TopLevelChain& chain, Haplotypes::TopLevelChain& output) const {
    std::vector<Subchain> subchains = this->get_subchains(chain);
    for (const Subchain& subchain : subchains) {
        std::vector<std::pair<Subchain, std::vector<sequence_type>>> to_process;
        auto sequences = this->get_sequences(subchain);
        if (sequences.empty()) {
            // There are no haplotypes crossing the subchain, so we break it into
            // a suffix and a prefix.
            to_process.push_back({ { Haplotypes::Subchain::suffix, subchain.start, empty_handle() }, this->get_sequences(subchain.start) });
            to_process.push_back({ { Haplotypes::Subchain::prefix, empty_handle(), subchain.end }, this->get_sequences(subchain.end) });
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
        handle_t successor = empty_handle();
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
}

std::ostream& Recombinator::Statistics::print(std::ostream& out) const {
    out << (this->chains - this->full_haplotypes) << " chains with " << this->subchains << " subchains and " << this->fragments << " fragments; "
        << this->full_haplotypes << " chains with full haplotypes";
    return out;
}

//------------------------------------------------------------------------------

Recombinator::Recombinator(const gbwtgraph::GBZ& gbz, HaplotypePartitioner::Verbosity verbosity) :
    gbz(gbz), verbosity(verbosity)
{
}

//------------------------------------------------------------------------------

gbwt::GBWT Recombinator::generate_haplotypes(const Haplotypes& haplotypes) const {
    double start = gbwt::readTimer();
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
    double checkpoint = gbwt::readTimer();
    if (this->verbosity >= HaplotypePartitioner::verbosity_basic) {
        std::cerr << "Running " << omp_get_max_threads() << " GBWT construction jobs in parallel" << std::endl;
    }
    std::vector<gbwt::GBWT> indexes(jobs.size());
    Statistics statistics;
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t job = 0; job < jobs.size(); job++) {
        // FIXME buffer size?
        gbwt::GBWTBuilder builder(sdsl::bits::length(this->gbz.index.sigma() - 1));
        builder.index.addMetadata();
        Statistics job_statistics;
        for (auto chain_id : jobs[job]) {
            try {
                Statistics chain_statistics = this->generate_haplotypes(haplotypes.chains[chain_id], builder);
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

Recombinator::Statistics Recombinator::generate_haplotypes(const Haplotypes::TopLevelChain& chain, gbwt::GBWTBuilder& builder) const {
    std::vector<Haplotype> haplotypes;
    for (size_t i = 0; i < NUM_HAPLOTYPES; i++) {
        haplotypes.push_back({ chain.offset, i, 0, gbwt::invalid_edge(), {} });
    }

    Statistics statistics; statistics.chains = 1;
    if (chain.subchains.size() == 1 && chain.subchains.front().type == Haplotypes::Subchain::full_haplotype) {
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
            for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                assert(!subchain.sequences.empty());
                size_t seq = haplotype % subchain.sequences.size();
                haplotypes[haplotype].extend(subchain.sequences[seq], subchain, *this, builder);
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

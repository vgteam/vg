#include "recombinator.hpp"

#include "kff.hpp"
#include "statistics.hpp"
#include "algorithms/component.hpp"
#include "algorithms/extract_subchain.hpp"

#include <algorithm>
#include <cmath>
#include <map>

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr std::uint32_t Haplotypes::Header::MAGIC_NUMBER;
constexpr std::uint32_t Haplotypes::Header::VERSION;
constexpr std::uint32_t Haplotypes::Header::MIN_VERSION;
constexpr std::uint64_t Haplotypes::Header::DEFAULT_K;

constexpr size_t HaplotypePartitioner::SUBCHAIN_LENGTH;
constexpr size_t HaplotypePartitioner::APPROXIMATE_JOBS;

constexpr size_t Recombinator::NUM_HAPLOTYPES;
constexpr size_t Recombinator::NUM_CANDIDATES;
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
                result[subchain.kmers[kmer_id]] = 0;
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

size_t Haplotypes::Subchain::distance(const gbwtgraph::GBZ& gbz, size_t i) const {
    if (this->type != normal || i >= this->sequences.size()) {
        return 0;
    }

    size_t result = 1;
    gbwt::edge_type curr(this->start, this->sequences[i].second);
    while (true) {
        curr = gbz.index.LF(curr);
        if (curr.first == gbwt::ENDMARKER || curr.first == this->end) {
            break;
        }
        result += gbz.graph.get_length(gbwtgraph::GBWTGraph::node_to_handle(curr.first));
    }

    return result;
}

// TODO: What is the right formula?
double Haplotypes::Subchain::badness(const gbwtgraph::GBZ& gbz) const {
    double result = 0.0;

    // Factor 1: Subchain length, ideally over a reference path.
    if (this->type == normal) {
        size_t expected_length = HaplotypePartitioner::SUBCHAIN_LENGTH;
        size_t selected = 0;
        for (size_t i = 0; i < this->sequences.size(); i++) {
            gbwt::size_type path_id = gbwt::Path::id(this->sequences[i].first);
            auto found = gbz.graph.id_to_path.find(path_id);
            if (found != gbz.graph.id_to_path.end()) {
                selected = i; break;
            }
        }
        size_t length = this->distance(gbz, selected);
        if (length > expected_length) {
            result += std::log(static_cast<double>(length) / static_cast<double>(expected_length));
        }
    }

    // Factor 2: Number of haplotypes relative to the expected number.
    size_t expected_haplotypes = gbz.index.metadata.haplotypes();
    size_t haplotypes = this->sequences.size();
    if (haplotypes < expected_haplotypes) {
        result += std::log(static_cast<double>(expected_haplotypes) / static_cast<double>(haplotypes));
    }

    // Factor 3: Information content of the kmers.
    // Disabled for the moment.
    /* double expected_entropy = 4.0 * std::log(static_cast<double>(haplotypes));
    double entropy = 0.0;
    for (size_t i = 0; i < this->kmer_counts.size(); i++) {
        double p = static_cast<double>(this->kmer_counts[i]) / static_cast<double>(haplotypes);
        entropy -= p * std::log(p);
    }
    if (entropy < expected_entropy) {
        result += expected_entropy - entropy;
    } */

    return result;
}

void Haplotypes::Subchain::simple_sds_serialize(std::ostream& out) const {
    sdsl::simple_sds::serialize_value<std::uint64_t>(this->type, out);
    sdsl::simple_sds::serialize_value<gbwt::node_type>(this->start, out);
    sdsl::simple_sds::serialize_value<gbwt::node_type>(this->end, out);
    sdsl::simple_sds::serialize_vector(this->kmers, out);
    this->kmer_counts.simple_sds_serialize(out);
    sdsl::simple_sds::serialize_vector(this->sequences, out);
    this->kmers_present.simple_sds_serialize(out);
}

void load_subchain_header(Haplotypes::Subchain& subchain, std::istream& in) {
    std::uint64_t temp = sdsl::simple_sds::load_value<std::uint64_t>(in);
    switch (temp) {
    case Haplotypes::Subchain::normal: // Fall through.
    case Haplotypes::Subchain::prefix: // Fall through.
    case Haplotypes::Subchain::suffix: // Fall through.
    case Haplotypes::Subchain::full_haplotype:
        subchain.type = static_cast<Haplotypes::Subchain::subchain_t>(temp);
        break;
    default:
        throw sdsl::simple_sds::InvalidData("Invalid subchain type: " + std::to_string(temp));
    }

    subchain.start = sdsl::simple_sds::load_value<gbwt::node_type>(in);
    subchain.end = sdsl::simple_sds::load_value<gbwt::node_type>(in);
    bool should_have_start = (subchain.type == Haplotypes::Subchain::normal || subchain.type == Haplotypes::Subchain::suffix);
    bool should_have_end = (subchain.type == Haplotypes::Subchain::normal || subchain.type == Haplotypes::Subchain::prefix);
    if ((subchain.start != gbwt::ENDMARKER) != should_have_start) {
        throw sdsl::simple_sds::InvalidData("Subchain start node " + std::to_string(subchain.start) + " does not match type " + std::to_string(temp));
    }
    if ((subchain.end != gbwt::ENDMARKER) != should_have_end) {
        throw sdsl::simple_sds::InvalidData("Subchain end node " + std::to_string(subchain.end) + " does not match type" + std::to_string(temp));
    }
}

void load_subchain_kmers_present(Haplotypes::Subchain& subchain, std::istream& in) {
    subchain.kmers_present.simple_sds_load(in);
    if (subchain.kmers_present.size() != subchain.kmers.size() * subchain.sequences.size()) {
        throw sdsl::simple_sds::InvalidData("Invalid length for the kmer presence bitvector in subchain from " +
            std::to_string(subchain.start) + " to " + std::to_string(subchain.end));
    }
}

void Haplotypes::Subchain::load_v1(std::istream& in) {
    load_subchain_header(*this, in);

    // Kmer and sequence information must be converted to a more compact format.
    auto kmers_counts = sdsl::simple_sds::load_vector<std::pair<kmer_type, size_t>>(in);
    auto seqs = sdsl::simple_sds::load_vector<sequence_type>(in);
    this->kmers = std::vector<kmer_type>(kmers_counts.size());
    this->kmer_counts = sdsl::int_vector<0>(kmers_counts.size(), 0, sdsl::bits::length(this->sequences.size()));
    for (size_t i = 0; i < kmers_counts.size(); i++) {
        this->kmers[i] = kmers_counts[i].first;
        this->kmer_counts[i] = kmers_counts[i].second;
    }
    this->sequences = std::vector<compact_sequence_type>(seqs.size());
    for (size_t i = 0; i < seqs.size(); i++) {
        this->sequences[i].first = seqs[i].first;
        this->sequences[i].second = seqs[i].second;
    }

    load_subchain_kmers_present(*this, in);
}

void Haplotypes::Subchain::simple_sds_load(std::istream& in) {
    load_subchain_header(*this, in);

    this->kmers = sdsl::simple_sds::load_vector<kmer_type>(in);
    this->kmer_counts.simple_sds_load(in);
    this->sequences = sdsl::simple_sds::load_vector<compact_sequence_type>(in);

    load_subchain_kmers_present(*this, in);
}

size_t Haplotypes::Subchain::simple_sds_size() const {
    size_t result = sdsl::simple_sds::value_size<std::uint64_t>() + 2 * sdsl::simple_sds::value_size<gbwt::node_type>();
    result += sdsl::simple_sds::vector_size(this->kmers);
    result += this->kmer_counts.simple_sds_size();
    result += sdsl::simple_sds::vector_size(this->sequences);
    result += this->kmers_present.simple_sds_size();
    return result;
}

void Haplotypes::TopLevelChain::simple_sds_serialize(std::ostream& out) const {
    sdsl::simple_sds::serialize_value<size_t>(this->offset, out);
    sdsl::simple_sds::serialize_value<size_t>(this->job_id, out);
    sdsl::simple_sds::serialize_string(this->contig_name, out);
    sdsl::simple_sds::serialize_value<size_t>(this->subchains.size(), out);
    for (auto& subchain : this->subchains) {
        subchain.simple_sds_serialize(out);
    }
}

void Haplotypes::TopLevelChain::simple_sds_load(std::istream& in) {
    this->offset = sdsl::simple_sds::load_value<size_t>(in);
    this->job_id = sdsl::simple_sds::load_value<size_t>(in);
    this->contig_name = sdsl::simple_sds::load_string(in);
    size_t subchain_count = sdsl::simple_sds::load_value<size_t>(in);
    this->subchains.resize(subchain_count);
    for (size_t i = 0; i < subchain_count; i++) {
        this->subchains[i].simple_sds_load(in);
    }
}

void Haplotypes::TopLevelChain::load_v1(std::istream& in) {
    this->offset = sdsl::simple_sds::load_value<size_t>(in);
    this->job_id = sdsl::simple_sds::load_value<size_t>(in);
    this->contig_name = "component_" + std::to_string(this->offset);
    size_t subchain_count = sdsl::simple_sds::load_value<size_t>(in);
    this->subchains.resize(subchain_count);
    for (size_t i = 0; i < subchain_count; i++) {
        this->subchains[i].load_v1(in);
    }
}

void Haplotypes::TopLevelChain::load_v2(std::istream& in) {
    this->offset = sdsl::simple_sds::load_value<size_t>(in);
    this->job_id = sdsl::simple_sds::load_value<size_t>(in);
    this->contig_name = sdsl::simple_sds::load_string(in);
    size_t subchain_count = sdsl::simple_sds::load_value<size_t>(in);
    this->subchains.resize(subchain_count);
    for (size_t i = 0; i < subchain_count; i++) {
        this->subchains[i].load_v1(in);
    }
}

size_t Haplotypes::TopLevelChain::simple_sds_size() const {
    size_t result = 3 * sdsl::simple_sds::value_size<size_t>();
    result += sdsl::simple_sds::string_size(this->contig_name);
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
    if (this->header.version < Header::MIN_VERSION || this->header.version > Header::VERSION) {
        std::string msg = "Haplotypes::simple_sds_load(): Expected version " + std::to_string(Header::MIN_VERSION)
            + " to " + std::to_string(Header::VERSION) + ", got version " + std::to_string(this->header.version);
        throw sdsl::simple_sds::InvalidData(msg);
    }

    this->jobs_for_cached_paths = sdsl::simple_sds::load_vector<size_t>(in);

    this->chains.resize(this->header.top_level_chains);
    for (auto& chain : this->chains) {
        if (this->header.version >= 3) {
            chain.simple_sds_load(in);
        } else if (this->header.version == 2) {
            chain.load_v2(in);
        } else {
            chain.load_v1(in);
        }
    }

    // Update to the current version.
    this->header.version = Header::VERSION;
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

    Haplotypes result;
    result.header.k = this->minimizer_index.k();

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
    this->distance_index.for_each_child(this->distance_index.get_root(), [&](const handlegraph::net_handle_t&) {
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
        result.chains[chain_id].contig_name = "component_" + std::to_string(chain_id); // Placeholder name.
    }

    // Assign chains to jobs and fill in contig names.
    auto chains_by_job = gbwtgraph::partition_chains(this->distance_index, this->gbz.graph, jobs);
    // We do not use a path filter, because a GBZ graph should not contain alt paths.
    auto contig_names = jobs.contig_names(this->gbz.graph);
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

    // Assign named and reference paths to jobs.
    result.jobs_for_cached_paths = std::vector<size_t>(this->gbz.graph.named_paths.size(), result.jobs());
    // Again, we do not use a path filter, because a GBZ graph should not contain alt paths.
    auto assignments = gbwtgraph::assign_paths(this->gbz.graph, jobs, nullptr, nullptr);
    for (size_t job = 0; job < assignments.size(); job++) {
        for (const path_handle_t& path : assignments[job]) {
            // We know that GBWTGraph path handles for reference/generic paths are offsets in named_paths.
            size_t path_id = handlegraph::as_integer(path);
            if (path_id >= result.jobs_for_cached_paths.size()) {
                std::string msg = "HaplotypePartitioner::partition_haplotypes(): path " + std::to_string(path_id) + " is not a reference/generic path";
                throw std::runtime_error(msg);
            }
            result.jobs_for_cached_paths[path_id] = job;
        }
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
                std::cerr << "Finished job " << job << " with " << chains.size() << " chains (";
                for (size_t i = 0; i < chains.size(); i++) {
                    if (i > 0) {
                        std::cerr << ", ";
                    }
                    std::cerr << result.chains[chains[i].offset].contig_name;
                }
                std::cerr << "), " << total_subchains << " subchains, and " << total_kmers << " kmers" << std::endl;
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
    net_handle_t curr = this->distance_index.get_bound(chain.chain, false, true);
    net_handle_t chain_end = this->distance_index.get_bound(chain.chain, true, false);
    while (curr != chain_end) {
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
    const HaplotypePartitioner& partitioner
) {
    // Determine the start for a prefix.
    bool single_fragment = (pos.first == gbwt::ENDMARKER || end == empty_gbwtgraph_handle());
    if (pos.first == gbwt::ENDMARKER) {
        pos = partitioner.gbz.index.start(sequence_id);
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
    gbwtgraph::view_type view = partitioner.gbz.graph.get_sequence_view(curr);
    size_t offset = (view.second > start_max ? view.second - start_max : 0);
    haplotype.emplace_back();
    haplotype.back().append(view.first + offset, view.second - offset);

    while (true) {
        pos = partitioner.gbz.index.LF(pos);
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
                pos = partitioner.gbz.index.start(sequence_id);
            }
            if (sequence_id == gbwt::invalid_sequence()) {
                break;
            }
            haplotype.emplace_back();
        }
        curr = gbwtgraph::GBWTGraph::node_to_handle(pos.first);
        view = partitioner.gbz.graph.get_sequence_view(curr);
        if (curr == end) {
            haplotype.back().append(view.first, std::min(view.second, end_max));
            break;
        } else {
            haplotype.back().append(view.first, view.second);
        }
    }

    return haplotype;
}

// Return the sorted set of kmers that are minimizers in the sequences and have a single
// occurrence in the graph.
std::vector<HaplotypePartitioner::kmer_type> take_unique_minimizers(
    const std::vector<std::string>& sequences, const HaplotypePartitioner::minimizer_index_type& minimizer_index
) {
    std::vector<HaplotypePartitioner::kmer_type> result;
    for (const std::string& sequence : sequences) {
        auto minimizers = minimizer_index.minimizers(sequence);
        for (auto& minimizer : minimizers) {
            if (minimizer_index.count(minimizer) == 1) {
                result.push_back(minimizer.key.get_key());
            }
        }
    }
    gbwt::removeDuplicates(result, false);
    return result;
}

std::vector<HaplotypePartitioner::kmer_type> HaplotypePartitioner::unique_minimizers(gbwt::size_type sequence_id) const {
    size_t limit = std::numeric_limits<size_t>::max();
    std::vector<std::string> haplotype = generate_haplotype(
        sequence_id,
        gbwt::edge_type(gbwt::ENDMARKER, 0), empty_gbwtgraph_handle(), limit, limit,
        *this
    );
    return take_unique_minimizers(haplotype, this->minimizer_index);
}

std::vector<HaplotypePartitioner::kmer_type> HaplotypePartitioner::unique_minimizers(sequence_type sequence, Subchain subchain, size_t& fragments) const {
    gbwt::edge_type pos(gbwt::ENDMARKER, 0);
    size_t start_max = std::numeric_limits<size_t>::max(), end_max = this->minimizer_index.k() - 1;
    if (subchain.has_start()) {
        pos = gbwt::edge_type(gbwtgraph::GBWTGraph::handle_to_node(subchain.start), sequence.second);
        start_max = this->minimizer_index.k() - 1;
    }
    std::vector<std::string> haplotype = generate_haplotype(
        sequence.first,
        pos, subchain.end, start_max, end_max,
        *this
    );
    fragments = haplotype.size();
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
    std::vector<HaplotypePartitioner::kmer_type>& all_kmers,
    sdsl::int_vector<0>& kmer_counts,
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
    kmer_counts = sdsl::int_vector<0>(present.size(), 0, sdsl::bits::length(sequences.size()));
    size_t offset = 0;
    for (auto iter = present.begin(); iter != present.end(); ++iter) {
        if (iter->second.second < sequences.size()) {
            all_kmers.push_back(iter->first);
            kmer_counts[offset] = iter->second.second;
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

    // Convert the subchains to actual subchains.
    size_t fragmented_subchains = 0, fragmented_sequences = 0, additional_fragments = 0;
    for (const Subchain& subchain : subchains) {
        std::vector<std::pair<Subchain, std::vector<sequence_type>>> to_process;
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
        for (auto iter = to_process.begin(); iter != to_process.end(); ++iter) {
            bool fragmented = false;
            output.subchains.push_back({
                iter->first.type,
                gbwtgraph::GBWTGraph::handle_to_node(iter->first.start), gbwtgraph::GBWTGraph::handle_to_node(iter->first.end),
                {}, sdsl::int_vector<0>(0, 0, 64), {}, sdsl::bit_vector()
            });
            Haplotypes::Subchain& subchain = output.subchains.back();
            std::vector<std::vector<kmer_type>> kmers_by_sequence;
            kmers_by_sequence.reserve(iter->second.size());
            for (sequence_type sequence : iter->second) {
                size_t fragments = 0;
                kmers_by_sequence.emplace_back(this->unique_minimizers(sequence, iter->first, fragments));
                if (fragments > 1) {
                    fragmented = true;
                    fragmented_sequences++;
                    additional_fragments += fragments - 1;
                }
            }
            present_kmers(kmers_by_sequence, subchain.kmers, subchain.kmer_counts, subchain.kmers_present);
            subchain.sequences = std::vector<Haplotypes::compact_sequence_type>(iter->second.size());
            for (size_t i = 0; i < iter->second.size(); i++) {
                subchain.sequences[i] = Haplotypes::compact_sequence_type(iter->second[i].first, iter->second[i].second);
            }
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
        output.subchains.push_back({
            Haplotypes::Subchain::full_haplotype,
            gbwt::ENDMARKER, gbwt::ENDMARKER,
            {}, sdsl::int_vector<0>(0, 0, 64), {}, sdsl::bit_vector()
        });
        Haplotypes::Subchain& subchain = output.subchains.back();
        gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(chain.handle);
        auto sequences = this->r_index.decompressDA(node);
        std::vector<std::vector<kmer_type>> kmers_by_sequence;
        kmers_by_sequence.reserve(sequences.size());
        for (auto seq_id : sequences) {
            kmers_by_sequence.emplace_back(this->unique_minimizers(seq_id));
        }
        present_kmers(kmers_by_sequence, subchain.kmers, subchain.kmer_counts, subchain.kmers_present);
        subchain.sequences.reserve(sequences.size());
        for (size_t i = 0; i < sequences.size(); i++) {
            subchain.sequences.push_back(Haplotypes::compact_sequence_type(sequences[i], 0));
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
        const std::string& contig_name, size_t id
    ) :
        recombinator(recombinator),
        builder(builder), metadata(metadata),
        contig_name(contig_name), id(id), fragment(0),
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
     */
    void extend(sequence_type sequence, const Haplotypes::Subchain& subchain);

    // Takes an existing haplotype from the GBWT index and inserts it into
    // the builder. This is intended for fragments that do not contain
    // subchains crossed by the original haplotypes. The call will fail if
    // `extend()` has been called.
    void take(gbwt::size_type sequence_id);

    // Finishes the current haplotype fragment and starts a new one.
    // If `with_suffix` is true, the current fragment will be extended
    // until the end of the corresponding path. In that case, the call will
    // fail if `extend()` has not been called for this fragment.
    void finish(bool with_suffix);

private:
    // Extends the haplotype over a unary path from a previous subchain.
    void connect(gbwt::node_type until);

    // Takes a prefix of the current sequence until `until`, assuming that
    // `sequence_id` has been set.
    void prefix(gbwt::node_type until);

    // Extends the haplotype from the previous subchain until the end.
    void suffix();

    // Inserts the current fragment into the builder.
    void insert();
};

void RecombinatorHaplotype::extend(sequence_type sequence, const Haplotypes::Subchain& subchain) {
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
        this->finish(true);
        return;
    }

    // This is a normal subchain.
    while (this->position.first != subchain.end) {
        this->position = index.LF(this->position);
        if (this->position.first == gbwt::ENDMARKER) {
            gbwt::size_type prev = this->sequence_id;
            this->finish(false);
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

void RecombinatorHaplotype::finish(bool with_suffix) {
    if (with_suffix) {
        if (this->position == gbwt::invalid_edge()) {
            throw std::runtime_error("Haplotype::finish(): there is no current position");
        }
        this->suffix();
    }
    this->insert();
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

void RecombinatorHaplotype::insert() {
    std::string sample_name = "recombination"; // TODO: Make this a static class variable.
    this->metadata.add_haplotype(sample_name, this->contig_name, this->id, this->fragment);
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
        while (pos.first == gbwt::ENDMARKER) {
            finish_fragment();
            sequence = fragment_map.oriented_next(sequence);
            if (sequence == gbwt::invalid_sequence()) {
                std::string msg = "RecombinatorFragment::generate(): no successor for GBWT sequence " + std::to_string(this->sequence_id);
                throw std::runtime_error(msg);
            }
            pos = index.start(sequence);
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
}

//------------------------------------------------------------------------------

Recombinator::Parameters::Parameters(preset_t preset) {
    if (preset == preset_haploid) {
        this->haploid_scoring = true;
        this->include_reference = true;
    } else if (preset == preset_diploid) {
        this->num_haplotypes = NUM_CANDIDATES;
        this->diploid_sampling = true;
        this->include_reference = true;
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
        out << "- include reference paths" << std::endl;
    }
}

//------------------------------------------------------------------------------

void add_path(const gbwt::GBWT& source, gbwt::size_type path_id, gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata) {
    // We know that sufficient metadata exists, because this is a cached path.
    gbwt::PathName path_name = source.metadata.path(path_id);
    std::string sample_name = source.metadata.sample(path_name.sample);
    std::string contig_name = source.metadata.contig(path_name.contig);
    if (sample_name == gbwtgraph::REFERENCE_PATH_SAMPLE_NAME) {
        metadata.add_generic_path(contig_name);
    } else {
        // Reference samples will be copied later.
        metadata.add_haplotype(sample_name, contig_name, path_name.phase, path_name.count);
    }

    gbwt::vector_type path = source.extract(gbwt::Path::encode(path_id, false));
    builder.insert(path, true);
}

//------------------------------------------------------------------------------

void recombinator_sanity_checks(const Recombinator::Parameters& parameters) {
    if (parameters.num_haplotypes == 0) {
        std::string msg = "recombinator_sanity_checks(): number of haplotypes cannot be 0";
        throw std::runtime_error(msg);
    }
    if (parameters.badness_threshold <= 0.0) {
        std::string msg = "recombinator_sanity_checks(): badness threshold must be positive";
        throw std::runtime_error(msg);
    }
    if (parameters.diploid_sampling && parameters.num_haplotypes < 2) {
        std::string msg = "recombinator_sanity_checks(): diploid sampling requires at least 2 haplotypes";
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

    // Use mode as the initial estimate for coverage.
    auto statistics = summary_statistics(count_to_frequency);
    double coverage = statistics.mode;
    bool reliable = true;
    if (verbosity >= Haplotypes::verbosity_detailed) {
        std::cerr << "Coverage: median " << statistics.median
            << ", mean " << statistics.mean
            << ", stdev " << statistics.stdev
            << ", mode " << statistics.mode;
    }

    // In the default (non-haploid) scoring model, if mode < median, we try
    // to find a secondary peak at ~2x mode and use it if it is good enough.
    if (statistics.mode < statistics.median && !parameters.haploid_scoring) {
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
        } else {
            reliable = false;
        }
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

    // Sanity checks (may throw).
    recombinator_sanity_checks(parameters);

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

    // Figure out GBWT path ids for reference paths in each job.
    std::vector<std::vector<gbwt::size_type>> reference_paths(this->haplotypes.jobs());
    if (parameters.include_reference) {
        for (size_t i = 0; i < this->gbz.graph.named_paths.size(); i++) {
            size_t job_id = this->haplotypes.jobs_for_cached_paths[i];
            if (job_id < this->haplotypes.jobs()) {
                reference_paths[job_id].push_back(this->gbz.graph.named_paths[i].id);
            }
        }
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
                Statistics chain_statistics = this->generate_haplotypes(
                    this->haplotypes.chains[chain_id], counts, builder, metadata, parameters, coverage
                );
                job_statistics.combine(chain_statistics);
            } catch (const std::runtime_error& e) {
                std::cerr << "error: [job " << job << "]: " << e.what() << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        builder.finish();
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
    if (parameters.include_reference) {
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

// Classify the kmers in the subchain according to the kmer counts and the
// coverage. Return a vector of kmer types and their initial scores. Update the
// statistics with the number of non-frequent kmers if necessary.
std::vector<std::pair<Recombinator::kmer_presence, double>> classify_kmers(
    const Haplotypes::Subchain& subchain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    double coverage,
    Recombinator::Statistics* statistics,
    const Recombinator::Parameters& parameters
) {
    // TODO: What are the proper thresholds?
    double absent_threshold = coverage * 0.1;
    double heterozygous_threshold = coverage / std::log(4.0);
    double homozygous_threshold = coverage * 2.5;

    // Determine the type of each kmer in the sample and the score for the kmer.
    // A haplotype with the kmer gets +1.0 * score, while a haplotype without it
    // gets -1.0 * score.
    std::vector<std::pair<Recombinator::kmer_presence, double>> kmer_types;
    size_t selected_kmers = 0;
    for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
        double count = kmer_counts.at(subchain.kmers[kmer_id]);
        if (count < absent_threshold) {
            kmer_types.push_back({ Recombinator::absent, -1.0 * parameters.absent_score });
            selected_kmers++;
        } else if (count < heterozygous_threshold && !parameters.haploid_scoring) {
            kmer_types.push_back({ Recombinator::heterozygous, 0.0 });
            selected_kmers++;
        } else if (count < homozygous_threshold) {
            kmer_types.push_back({ Recombinator::present, 1.0 });
            selected_kmers++;
        } else {
            kmer_types.push_back({ Recombinator::frequent, 0.0 });
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
    const std::vector<std::pair<Recombinator::kmer_presence, double>>& kmer_types,
    size_t original_selected,
    const Recombinator::Parameters& parameters
) {
    std::int64_t best_score = std::numeric_limits<std::int64_t>::min();
    size_t best_left = 0, best_right = 1;

    for (size_t left = 0; left < candidates.size(); left++) {
        size_t left_offset = candidates[left].first * subchain.kmers.size();
        for (size_t right = left + 1; right < candidates.size(); right++) {
            std::int64_t score = 0;
            size_t right_offset = candidates[right].first * subchain.kmers.size();
            for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                int64_t found = subchain.kmers_present[left_offset + kmer_id] + subchain.kmers_present[right_offset + kmer_id];
                switch (kmer_types[kmer_id].first) {
                case Recombinator::absent:
                    score += 1 - found; // +1 for 0, 0 for 1, -1 for 2
                    break;
                case Recombinator::heterozygous:
                    score += (found == 1 ? 1 : 0);
                    break;
                case Recombinator::present:
                    score += found - 1; // -1 for 0, 0 for 1, +1 for 2
                    break;
                default:
                    break;
                }
            }
            if (score > best_score) {
                best_score = score;
                best_left = left;
                best_right = right;
            }
        }
    }

    // If this is a bad subchain, we move the selected haplotypes to the front
    // and return the remaining non-duplicated haplotypes as extra fragments.
    // Otherwise we return only the selected haplotypes.
    if (parameters.extra_fragments) {
        double badness = subchain.badness(gbz);
        if (badness > parameters.badness_threshold) {
            std::vector<std::pair<size_t, double>> result;
            result.push_back(candidates[best_left]);
            result.push_back(candidates[best_right]);
            for (size_t i = 0; i < original_selected; i++) {
                if (i != best_left && i != best_right) {
                    result.push_back(candidates[i]);
                }
            }
            return result;
        }
    }

    return { candidates[best_left], candidates[best_right] };
}

// Returns (sequence offset in the subchain, score).
// Updates statistics with the number of kmers used if provided.
// Updates the local haplotypes with scores and ranks in each round of selection
// if provided.
std::vector<std::pair<size_t, double>> select_haplotypes(
    const gbwtgraph::GBZ& gbz,
    const Haplotypes::Subchain& subchain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    double coverage,
    Recombinator::Statistics* statistics,
    std::vector<Recombinator::LocalHaplotype>* local_haplotypes,
    const Recombinator::Parameters& parameters
) {
    // Classify the kmers.
    std::vector<std::pair<Recombinator::kmer_presence, double>> kmer_types = classify_kmers(subchain, kmer_counts, coverage, statistics, parameters);

    // Select the haplotypes greedily.
    std::vector<std::pair<size_t, double>> selected_haplotypes;
    std::vector<std::pair<size_t, double>> remaining_haplotypes;
    for (size_t seq_offset = 0; seq_offset < subchain.sequences.size(); seq_offset++) {
        remaining_haplotypes.push_back( { seq_offset, 0.0 });
    }
    while (selected_haplotypes.size() < parameters.num_haplotypes && !remaining_haplotypes.empty()) {
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

        // Select the highest-scoring haplotype.
        size_t selected = 0;
        for (size_t i = 1; i < remaining_haplotypes.size(); i++) {
            if (remaining_haplotypes[i].second > remaining_haplotypes[selected].second) {
                selected = i;
            }
        }
        selected_haplotypes.push_back(remaining_haplotypes[selected]);
        remaining_haplotypes.erase(remaining_haplotypes.begin() + selected);

        // Adjust kmer scores based on the selected haplotype.
        size_t offset = selected_haplotypes.back().first * subchain.kmers.size();
        for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
            switch (kmer_types[kmer_id].first) {
            case Recombinator::heterozygous:
                kmer_types[kmer_id].second += (subchain.kmers_present[offset + kmer_id] ? -1.0 : 1.0) * parameters.het_adjustment;
                break;
            case Recombinator::present:
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
    for (size_t i = original_selected; i < parameters.num_haplotypes; i++) {
        auto next = selected_haplotypes[i % original_selected];
        selected_haplotypes.push_back(next);
    }

    // Do diploid sampling. If this is a bad subchain, we also return the
    // extra fragments starting from the third haplotype. But we don't return
    // duplicated ones.
    if (parameters.diploid_sampling) {
        return select_diploid(gbz, subchain, selected_haplotypes, kmer_types, original_selected, parameters);
    } else {
        return selected_haplotypes;
    }
}

Recombinator::Statistics Recombinator::generate_haplotypes(const Haplotypes::TopLevelChain& chain,
    const hash_map<Haplotypes::Subchain::kmer_type, size_t>& kmer_counts,
    gbwt::GBWTBuilder& builder, gbwtgraph::MetadataBuilder& metadata,
    const Parameters& parameters,
    double coverage
) const {
    size_t final_haplotypes = (parameters.diploid_sampling ? 2 : parameters.num_haplotypes);
    std::vector<RecombinatorHaplotype> haplotypes;
    for (size_t i = 0; i < final_haplotypes; i++) {
        haplotypes.emplace_back(*this, builder, metadata, chain.contig_name, i + 1);
    }

    Statistics statistics;
    statistics.chains = 1; statistics.haplotypes = haplotypes.size();
    if (chain.subchains.size() == 1 && chain.subchains.front().type == Haplotypes::Subchain::full_haplotype) {
        // Full haplotypes should all be identical, because there are no snarls.
        // Therefore we do not need kmers.
        auto& subchain = chain.subchains.front();
        for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
            assert(!subchain.sequences.empty());
            size_t seq = haplotype % subchain.sequences.size();
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
            std::vector<std::pair<size_t, double>> selected_haplotypes = select_haplotypes(
                this->gbz, subchain, kmer_counts, coverage, &statistics, nullptr, parameters
            );
            if (parameters.diploid_sampling && selected_haplotypes.size() > 2) {
                for (size_t i = 2; i < selected_haplotypes.size(); i++) {
                    gbwt::size_type sequence_id = subchain.sequences[selected_haplotypes[i].first].first;
                    gbwt::edge_type start(subchain.start, subchain.sequences[selected_haplotypes[i].first].second);
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
                haplotypes[haplotype].extend(subchain.sequences[seq_offset], subchain);
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

} // namespace vg

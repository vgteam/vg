#include "recombinator.hpp"

#include <algorithm>
#include <unordered_set>

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t Recombinator::NUM_HAPLOTYPES;
constexpr size_t Recombinator::SUBCHAIN_LENGTH;
constexpr size_t Recombinator::APPROXIMATE_JOBS;

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

void Recombinator::Haplotype::extend(sequence_type sequence, Subchain subchain, const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
    if (subchain.type == Subchain::prefix) {
        if (!this->path.empty())  {
            std::cerr << "error: [Haplotype::extend()] got a prefix subchain after the start of a fragment" << std::endl;
            return;
        }
        this->prefix(recombinator.r_index.seqId(sequence.first), subchain.end, recombinator.gbz.index);
        return;
    }

    // Suffixes and normal subchains have a start node, so we must reach it first.
    if (!this->path.empty()) {
        this->connect(subchain.start, recombinator.gbz.graph);
    } else {
        this->prefix(recombinator.r_index.seqId(sequence.first), subchain.start, recombinator.gbz.index);
    }

    gbwt::edge_type curr(gbwtgraph::GBWTGraph::handle_to_node(subchain.start), sequence.second);
    if (!recombinator.gbz.index.contains(curr)) {
        std::cerr << "error: [Haplotype::extend()] the GBWT index does not contain position "; gbwt::operator<<(std::cerr, curr) << std::endl;
        return;
    }

    if (subchain.type == Subchain::suffix) {
        this->position = curr;
        this->finish(recombinator, builder);
        return;
    }

    // This is a normal subchain.
    gbwt::node_type end = gbwtgraph::GBWTGraph::handle_to_node(subchain.end);
    while (curr.first != end) {
        curr = recombinator.gbz.index.LF(curr);
        if (curr.first == gbwt::ENDMARKER) {
            std::cerr << "error: [Haplotype::extend()] the sequence did not reach " << to_string(subchain.end) << std::endl;
            return;
        }
        this->path.push_back(curr.first);
    }
    this->position = curr;
}

void Recombinator::Haplotype::take(gbwt::size_type sequence_id, const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
    if (!this->path.empty()) {
        std::cerr << "error: [Haplotype::take()] The current fragment is not empty" << std::endl;
        return;
    }
    if (sequence_id >= recombinator.gbz.index.sequences()) {
        std::cerr << "error: [Haplotype::take()] The GBWT index does not contain sequence " << sequence_id << std::endl;
        return;
    }
    this->path = recombinator.gbz.index.extract(sequence_id);
    this->insert(builder);
    this->fragment++;
    this->position = gbwt::invalid_edge();
    this->path.clear();
}

void Recombinator::Haplotype::finish(const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
    if (this->position == gbwt::invalid_edge()) {
        std::cerr << "error: [Haplotype::finish()] There is no current position" << std::endl;
        return;
    }
    this->suffix(recombinator.gbz.index);
    this->insert(builder);
    this->fragment++;
    this->position = gbwt::invalid_edge();
    this->path.clear();
}

void Recombinator::Haplotype::connect(handle_t until, const gbwtgraph::GBWTGraph& graph) {
    handle_t curr = gbwtgraph::GBWTGraph::node_to_handle(this->position.first);
    this->position = gbwt::invalid_edge();
    std::unordered_set<handle_t> visited;
    while (curr != until) {
        if (visited.find(curr) != visited.end()) {
            std::cerr << "error: [Haplotype::connect()] The path contains a cycle" << std::endl;
            return;
        }
        visited.insert(curr);
        handle_t successor = empty_handle();
        size_t successors = 0;
        graph.follow_edges(curr, false, [&](const handle_t& next) {
            successor = next;
            successors++;
        });
        if (successors != 1) {
            std::cerr << "error: [Haplotype::connect()] The path is not unary" << std::endl;
            return;
        }
        this->path.push_back(gbwtgraph::GBWTGraph::handle_to_node(successor));
        curr = successor;
    }
}

void Recombinator::Haplotype::prefix(gbwt::size_type sequence_id, handle_t until, const gbwt::GBWT& index) {
    this->position = gbwt::invalid_edge();
    if (sequence_id >= index.sequences()) {
        std::cerr << "error: [Haplotype::prefix()] Invalid GBWT sequence id " << sequence_id << std::endl;
        return;
    }
    gbwt::node_type end = gbwtgraph::GBWTGraph::handle_to_node(until);
    for (gbwt::edge_type curr = index.start(sequence_id); curr.first != gbwt::ENDMARKER; curr = index.LF(curr)) {
        this->path.push_back(curr.first);
        if (curr.first == end) {
            this->position = curr;
            return;
        }
    }
    std::cerr << "error: [Haplotype::prefix()] GBWT sequence " << sequence_id << " did not reach " << to_string(until) << std::endl;
    return;
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
    this->sequences += another.sequences;
    this->minimizers += another.minimizers;
}

std::ostream& Recombinator::Statistics::print(std::ostream& out) const {
    double seqs = static_cast<double>(this->sequences) / (this->subchains + this->full_haplotypes);
    double mins = static_cast<double>(this->minimizers) / this->sequences;
    out << (this->chains - this->full_haplotypes) << " chains with " << this->subchains << " subchains and " << this->fragments << " fragments; "
        << this->full_haplotypes << " chains with full haplotypes;"
        << " average " << seqs << " sequences with " << mins << " minimizers";
    return out;
}

//------------------------------------------------------------------------------

Recombinator::Recombinator(const gbwtgraph::GBZ& gbz,
    const gbwt::FastLocate& r_index,
    const SnarlDistanceIndex& distance_index,
    const gbwtgraph::DefaultMinimizerIndex& minimizer_index,
    Verbosity verbosity) :
    gbz(gbz), r_index(r_index), distance_index(distance_index), minimizer_index(minimizer_index),
    chains_by_job(),
    verbosity(verbosity)
{
    double start = gbwt::readTimer();
    if (this->verbosity >= verbosity_basic) {
        std::cerr << "Determining GBWT construction jobs" << std::endl;
    }
    size_t size_bound = this->gbz.graph.get_node_count() / APPROXIMATE_JOBS;
    gbwtgraph::ConstructionJobs jobs = gbwtgraph::gbwt_construction_jobs(this->gbz.graph, size_bound);
    this->chains_by_job = gbwtgraph::partition_chains(this->distance_index, this->gbz.graph, jobs);
    if (this->verbosity >= verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Partitioned " << jobs.components << " components into " << jobs.size() << " jobs in " << seconds << " seconds" << std::endl;
    }
}

//------------------------------------------------------------------------------

Recombinator::Statistics Recombinator::generate_haplotypes(const gbwtgraph::TopLevelChain& chain, gbwt::GBWTBuilder& builder) const {
    std::vector<Subchain> subchains = this->get_subchains(chain);
    std::vector<Haplotype> haplotypes;
    for (size_t i = 0; i < NUM_HAPLOTYPES; i++) {
        haplotypes.push_back({ chain.offset, i, 0, gbwt::invalid_edge(), {} });
    }

    bool have_haplotypes = false, generated_haplotypes = false;
    Statistics statistics; statistics.chains = 1;
    for (const Subchain& subchain : subchains) {
        std::vector<std::pair<Subchain, std::vector<sequence_type>>> to_process;
        auto sequences = this->get_sequences(subchain);
        if (sequences.empty()) {
            // There are no haplotypes crossing the subchain, so we break it into
            // a suffix and a prefix.
            to_process.push_back({ { Subchain::suffix, subchain.start, empty_handle() }, this->get_sequences(subchain.start) });
            to_process.push_back({ { Subchain::prefix, empty_handle(), subchain.end }, this->get_sequences(subchain.end) });
        } else {
            to_process.push_back({ subchain, std::move(sequences) });
        }
        for (auto iter = to_process.begin(); iter != to_process.end(); ++iter) {
            // TODO: It will eventually be faster to find all unique minimizers in the relevant subgraph,
            // align them to the haplotypes, and use the r-index to determine which haplotypes contain
            // that minimizer. On the other hand, that will cause some issues with haplotypes that visit
            // the same subchain multiple times.
            for (sequence_type sequence : iter->second) {
                statistics.sequences++;
                statistics.minimizers += this->unique_minimizers(sequence, iter->first);
            }
            /*
                FIXME Haplotype selection logic:
                * convert each sequence to a bitvector marking which selected minimizers occur in that sequence
                * use kmer counts from the reads to select sequences that are closest to the sample
                FIXME: What if there are multiple occurrences of the same minimizer?
            */
            for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
                // If this is a normal subchain, we have already checked that it's not empty.
                // For prefixes and suffixes, we list all GBWT path visits to the handle.
                // The list cannot be empty, because the handle exists in a GBWTGraph.
                assert(!iter->second.empty());
                haplotypes[haplotype].extend(iter->second[haplotype % iter->second.size()], iter->first, *this, builder);
            }
            have_haplotypes = iter->first.has_end();
            generated_haplotypes = true;
            statistics.subchains++;
        }
    }
    if (have_haplotypes) {
        for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
            haplotypes[haplotype].finish(*this, builder);
        }
        have_haplotypes = false;
    }

    // Take entire sequences if we could not generate any haplotypes.
    if (!generated_haplotypes) {
        gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(chain.handle);
        auto sequences = this->r_index.decompressDA(node);
        for (auto seq_id : sequences) {
            statistics.sequences++;
            statistics.minimizers += this->unique_minimizers(seq_id);
        }
        for (size_t haplotype = 0; haplotype < haplotypes.size(); haplotype++) {
            haplotypes[haplotype].take(sequences[haplotype % sequences.size()], *this, builder);
        }
        statistics.full_haplotypes = 1;
    } else {
        statistics.fragments = haplotypes.front().fragment;
    }

    return statistics;
}

//------------------------------------------------------------------------------

std::vector<Recombinator::Subchain>
Recombinator::get_subchains(const gbwtgraph::TopLevelChain& chain) const {
    std::vector<Subchain> result;

    // First pass: take all connected snarls as subchains.
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
                    snarls.push_back({ Subchain::prefix, empty_handle(), handle });
                } else {
                    size_t distance = this->get_distance(snarl_start, handle);
                    if (distance < std::numeric_limits<size_t>::max()) {
                        // Normal snarl with two boundary nodes.
                        snarls.push_back({ Subchain::normal, snarl_start, handle });
                    } else {
                        // The snarl is not connected, so we break it into two.
                        snarls.push_back({ Subchain::suffix, snarl_start, empty_handle() });
                        snarls.push_back({ Subchain::prefix, empty_handle(), handle });
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
            #pragma omp critical
            {
                std::cerr << "error: [Recombinator::get_subchains()] chain " << chain.offset << " has " << successors << " successors for a child" << std::endl;
            }
            return result;
        }
        curr = next;
    }
    if (was_snarl && has_start) {
        // If the chain ends with a snarl, we take it as a suffix.
        snarls.push_back({ Subchain::suffix, snarl_start, empty_handle() });
    }

    // Second pass: Combine snarls into subchains.
    size_t head = 0;
    while (head < snarls.size()) {
        if (snarls[head].type != Subchain::normal) {
            // Prefixes and suffixes should be rare, so we can simply use them directly.
            result.push_back(snarls[head]);
            head++;
            continue;
        }
        size_t tail = head;
        while (tail + 1 < snarls.size()) {
            if (snarls[tail + 1].type != Subchain::normal) {
                break;
            }
            size_t candidate = this->get_distance(snarls[head].start, snarls[tail + 1].end);
            if (candidate <= SUBCHAIN_LENGTH) {
                tail++;
            } else {
                break;
            }
        }
        result.push_back({ Subchain::normal, snarls[head].start, snarls[tail].end });
        head = tail + 1;
    }

    return result;
}

size_t Recombinator::get_distance(handle_t from, handle_t to) const {
    return this->distance_index.minimum_distance(
        this->gbz.graph.get_id(from), this->gbz.graph.get_is_reverse(from), this->gbz.graph.get_length(from) - 1,
        this->gbz.graph.get_id(to), this->gbz.graph.get_is_reverse(to), 0,
        false, &this->gbz.graph
    );
}

//------------------------------------------------------------------------------

std::vector<Recombinator::sequence_type> Recombinator::get_sequences(handle_t handle) const {
    std::vector<gbwt::size_type> sa = this->r_index.decompressSA(gbwtgraph::GBWTGraph::handle_to_node(handle));
    std::vector<sequence_type> result;
    result.reserve(sa.size());
    for (size_t i = 0; i < sa.size(); i++) {
        result.push_back({ sa[i], i });
    }
    std::sort(result.begin(), result.end(), [&](sequence_type a, sequence_type b) -> bool {
        return (r_index.seqId(a.first) < r_index.seqId(b.first));
    });
    return result;
}

std::vector<Recombinator::sequence_type> Recombinator::get_sequences(Subchain subchain) const {
    if (subchain.type == Subchain::prefix) {
        return this->get_sequences(subchain.end);
    }
    if (subchain.type == Subchain::suffix) {
        return this->get_sequences(subchain.start);
    }
    auto from = this->get_sequences(subchain.start);
    auto to = this->get_sequences(subchain.end);

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

    return result;
}

//------------------------------------------------------------------------------

void append_to(std::string& haplotype, handle_t handle, const gbwtgraph::GBWTGraph& graph) {
    gbwtgraph::view_type view = graph.get_sequence_view(handle);
    haplotype.append(view.first, view.second);
}

size_t count_unique_minimizers(const std::string& sequence, const gbwtgraph::DefaultMinimizerIndex& minimizer_index) {
    size_t result = 0;
    auto minimizers = minimizer_index.minimizers(sequence);
    // FIXME report duplicate occurrences on debug settings
    for (auto& minimizer : minimizers) {
        if (minimizer_index.count(minimizer) == 1) {
            result++;
        }
    }
    return result;
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

size_t Recombinator::unique_minimizers(gbwt::size_type sequence_id) const {
    gbwt::edge_type pos = this->gbz.index.start(sequence_id);
    std::string haplotype = generate_haplotype(pos, empty_handle(), this->gbz.graph);
    return count_unique_minimizers(haplotype, this->minimizer_index);
}

size_t Recombinator::unique_minimizers(sequence_type sequence, Subchain subchain) const {
    gbwt::edge_type pos;
    if (subchain.has_start()) {
        pos = gbwt::edge_type(gbwtgraph::GBWTGraph::handle_to_node(subchain.start), sequence.second);
    } else {
        pos = this->gbz.index.start(sequence.first);
    }
    std::string haplotype = generate_haplotype(pos, subchain.end, this->gbz.graph);
    return count_unique_minimizers(haplotype, this->minimizer_index);
}

//------------------------------------------------------------------------------

} // namespace vg

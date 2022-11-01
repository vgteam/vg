#include "recombinator.hpp"

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

// An empty subchain.
Recombinator::subchain_type empty_subchain() {
    return { empty_handle(), empty_handle() };
}

// Returns a GBWTGraph handle as a string (id, orientation).
std::string to_string(handle_t handle) {
    gbwt::node_type node = gbwtgraph::GBWTGraph::handle_to_node(handle);
    return std::string("(") + std::to_string(gbwt::Node::id(node)) + std::string(", ") + std::to_string(gbwt::Node::is_reverse(node)) + std::string(")");
}

//------------------------------------------------------------------------------

void Recombinator::Haplotype::extend(sequence_type sequence, subchain_type subchain, const Recombinator& recombinator) {
    if (this->position != gbwt::invalid_edge()) {
        this->connect(subchain.first, recombinator.gbz.graph);
    } else {
        this->prefix(recombinator.r_index.seqId(sequence.first), subchain.first, recombinator.gbz.index);
    }

    gbwt::edge_type curr(gbwtgraph::GBWTGraph::handle_to_node(subchain.first), sequence.second);
    if (!recombinator.gbz.index.contains(curr)) {
        std::cerr << "error: [Haplotype::extend()] The GBWT index does not contain position "; gbwt::operator<<(std::cerr, curr) << std::endl;
        return;
    }
    gbwt::node_type end = gbwtgraph::GBWTGraph::handle_to_node(subchain.second);
    while (curr.first != end) {
        curr = recombinator.gbz.index.LF(curr);
        if (curr.first == gbwt::ENDMARKER) {
            std::cerr << "error: [Haplotype::extend()] The sequence did not reach " << to_string(subchain.second) << std::endl;
            return;
        }
        this->path.push_back(curr.first);
    }
    this->position = curr;
}

void Recombinator::Haplotype::take(gbwt::size_type sequence, const Recombinator& recombinator, gbwt::GBWTBuilder& builder) {
    if (!this->path.empty()) {
        std::cerr << "error: [Haplotype::take()] The current fragment is not empty" << std::endl;
        return;
    }
    if (sequence >= recombinator.gbz.index.sequences()) {
        std::cerr << "error: [Haplotype::take()] The GBWT index does not contain sequence " << sequence << std::endl;
        return;
    }
    this->path = recombinator.gbz.index.extract(sequence);
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

void Recombinator::Haplotype::prefix(gbwt::size_type sequence, handle_t until, const gbwt::GBWT& index) {
    this->position = gbwt::invalid_edge();
    if (sequence >= index.sequences()) {
        std::cerr << "error: [Haplotype::prefix()] Invalid GBWT sequence id " << sequence << std::endl;
        return;
    }
    gbwt::node_type end = gbwtgraph::GBWTGraph::handle_to_node(until);
    for (gbwt::edge_type curr = index.start(sequence); curr.first != gbwt::ENDMARKER; curr = index.LF(curr)) {
        this->path.push_back(curr.first);
        if (curr.first == end) {
            return;
        }
    }
    std::cerr << "error: [Haplotype::prefix()] GBWT sequence " << sequence << " did not reach " << to_string(until) << std::endl;
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

Recombinator::Recombinator(const gbwtgraph::GBZ& gbz, const gbwt::FastLocate& r_index, const SnarlDistanceIndex& distance_index, bool progress) :
    gbz(gbz), r_index(r_index), distance_index(distance_index),
    chains_by_job(),
    progress(progress)
{
    double start = gbwt::readTimer();
    if (this->progress) {
        std::cerr << "Determining GBWT construction jobs" << std::endl;
    }
    size_t size_bound = this->gbz.graph.get_node_count() / APPROXIMATE_JOBS;
    gbwtgraph::ConstructionJobs jobs = gbwtgraph::gbwt_construction_jobs(this->gbz.graph, size_bound);
    this->chains_by_job = gbwtgraph::partition_chains(this->distance_index, this->gbz.graph, jobs);
    if (this->progress) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Partitioned " << jobs.components << " components into " << jobs.size() << " jobs in " << seconds << " seconds" << std::endl;
    }
}

//------------------------------------------------------------------------------

} // namespace vg

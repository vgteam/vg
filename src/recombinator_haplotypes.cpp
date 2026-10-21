#include "recombinator_haplotypes.hpp"

#include "kff.hpp"

#include <cmath>

//#define debug_selected_haplotypes

namespace vg {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_t Haplotypes::SUBCHAIN_LENGTH;
constexpr size_t Haplotypes::KFF_BLOCK_SIZE;

constexpr std::uint32_t Haplotypes::Header::MAGIC_NUMBER;
constexpr std::uint32_t Haplotypes::Header::VERSION;
constexpr std::uint32_t Haplotypes::Header::VERSION_WITH_TAGS;
constexpr std::uint32_t Haplotypes::Header::MIN_VERSION;
constexpr std::uint64_t Haplotypes::Header::DEFAULT_K;

//------------------------------------------------------------------------------

namespace {

// Returns the Hamming distance between the given rows of an uncompressed matrix.
// Stops as soon as the distance exceeds `limit` and returns a value greater than
// `limit` in that case.
size_t row_distance(
    const sdsl::bit_vector& matrix, size_t num_kmers,
    size_t row_a, size_t row_b, size_t limit
) {
    size_t start_a = row_a * num_kmers, start_b = row_b * num_kmers;
    size_t result = 0;
    for (size_t kmer_id = 0; kmer_id < num_kmers && result <= limit; kmer_id += 64) {
        size_t bits = std::min(static_cast<size_t>(64), num_kmers - kmer_id);
        std::uint64_t word_a = matrix.get_int(start_a + kmer_id, bits);
        std::uint64_t word_b = matrix.get_int(start_b + kmer_id, bits);
        result += sdsl::bits::cnt(word_a ^ word_b);
    }
    return result;
}

/*
    Returns the largest Hamming distance `d` for which it is worth storing a sequence
    as a child instead of a parent.

    An Elias-Fano encoded bitvector uses around `log2(gap) + 2.5` bits per set bit, and
    the average gap between the flipped bits in a child is around `num_kmers / d`. But
    because we concatenate the child bitvectors, the relevant gap measure is the
    average over all children in all clusters.

    TODO: In some situations, we could get better compression with a lower max distance
    and more clusters. Does that happen in practice and is the effect significant?
    Maybe we could expose a threshold multiplier as a parameter.

    NOTE: We could also improve compression a bit by selecting the final parents only
    after we have determined the clusters.
*/
size_t max_child_distance(size_t num_kmers) {
    if (num_kmers == 0) {
        return 0;
    }
    auto worth_it = [num_kmers](size_t d) -> bool {
        if (d == 0) {
            return true;
        }
        double cost = d * (std::log2(static_cast<double>(num_kmers) / d) + 2.5);
        return (cost <= static_cast<double>(num_kmers));
    };

    // Invariant: `worth_it(low)` is true and `worth_it(high)` is false.
    size_t low = 0, high = num_kmers + 1;
    while (high - low > 1) {
        size_t mid = low + (high - low) / 2;
        if (worth_it(mid)) { low = mid; } else { high = mid; }
    }
    return low;
}

} // anonymous namespace

KmerPresenceMatrix::KmerPresenceMatrix() : num_sequences(0), num_kmers(0) {
}

KmerPresenceMatrix::KmerPresenceMatrix(size_t num_sequences, size_t num_kmers, sdsl::bit_vector&& matrix, bool compress) :
    num_sequences(num_sequences), num_kmers(num_kmers),
    arrays(3 * num_sequences, 0, sdsl::bits::length((num_sequences == 0 ? 0 : num_sequences - 1)))
{
    if (matrix.size() != num_sequences * num_kmers) {
        throw std::runtime_error("KmerPresenceMatrix: inconsistent matrix size");
    }
    if (this->num_sequences == 0) {
        // Without sequences, there is no way of telling the number of kmers from the
        // serialized object. Normalize the dimensions to make the object round-trip.
        this->num_kmers = 0;
    }

    if (!compress) {
        for (size_t i = 0; i < this->num_sequences; i++) {
            this->arrays[i] = i;
            this->arrays[this->num_sequences + i] = i;
            this->arrays[2 * this->num_sequences + i] = i;
        }
        this->parents = std::move(matrix);
        this->children = sdsl::sd_vector<>();
        return;
    }

    // Greedy clustering: assign each sequence to the first parent within the threshold.
    size_t threshold = max_child_distance(this->num_kmers);
    std::vector<size_t> parent_of_cluster;
    std::vector<std::vector<size_t>> children_of_cluster;
    size_t total_flips = 0;
    for (size_t sequence = 0; sequence < this->num_sequences; sequence++) {
        size_t best_cluster = parent_of_cluster.size(), best_distance = threshold + 1;
        for (size_t cluster = 0; cluster < parent_of_cluster.size(); cluster++) {
            size_t distance = row_distance(
                matrix, this->num_kmers, sequence, parent_of_cluster[cluster], best_distance - 1
            );
            if (distance < best_distance) {
                best_cluster = cluster; best_distance = distance;
                if (distance == 0) { break; }
            }
        }
        if (best_cluster >= parent_of_cluster.size()) {
            parent_of_cluster.push_back(sequence);
            children_of_cluster.emplace_back();
        } else {
            children_of_cluster[best_cluster].push_back(sequence);
            total_flips += best_distance;
        }
    }
    size_t num_parents = parent_of_cluster.size();
    size_t num_children = this->num_sequences - num_parents;

    // Store the sequences in cluster order, with the parent of each cluster first.
    size_t rank = 0;
    for (size_t cluster = 0; cluster < num_parents; cluster++) {
        auto add_sequence = [&](size_t sequence) {
            this->arrays[sequence] = rank;
            this->arrays[this->num_sequences + rank] = sequence;
            this->arrays[2 * this->num_sequences + rank] = cluster;
            rank++;
        };
        add_sequence(parent_of_cluster[cluster]);
        for (size_t sequence : children_of_cluster[cluster]) {
            add_sequence(sequence);
        }
    }

    // Copy the parents.
    this->parents = sdsl::bit_vector(num_parents * this->num_kmers, 0);
    for (size_t cluster = 0; cluster < num_parents; cluster++) {
        size_t from = parent_of_cluster[cluster] * this->num_kmers;
        size_t to = cluster * this->num_kmers;
        for (size_t kmer_id = 0; kmer_id < this->num_kmers; kmer_id += 64) {
            size_t bits = std::min(static_cast<size_t>(64), this->num_kmers - kmer_id);
            this->parents.set_int(to + kmer_id, matrix.get_int(from + kmer_id, bits), bits);
        }
    }

    // Build the children by marking the bits that differ from the parent.
    if (num_children == 0 || this->num_kmers == 0) {
        this->children = sdsl::sd_vector<>();
        return;
    }
    sdsl::sd_vector_builder builder(num_children * this->num_kmers, total_flips);
    size_t child = 0;
    for (size_t cluster = 0; cluster < num_parents; cluster++) {
        size_t parent_start = parent_of_cluster[cluster] * this->num_kmers;
        for (size_t sequence : children_of_cluster[cluster]) {
            size_t child_start = child * this->num_kmers;
            size_t sequence_start = sequence * this->num_kmers;
            for (size_t kmer_id = 0; kmer_id < this->num_kmers; kmer_id += 64) {
                size_t bits = std::min(static_cast<size_t>(64), this->num_kmers - kmer_id);
                std::uint64_t word = matrix.get_int(sequence_start + kmer_id, bits)
                    ^ matrix.get_int(parent_start + kmer_id, bits);
                while (word != 0) {
                    size_t offset = sdsl::bits::lo(word);
                    builder.set_unsafe(child_start + kmer_id + offset);
                    word &= word - 1; // Clear the lowest set bit.
                }
            }
            child++;
        }
    }
    this->children = sdsl::sd_vector<>(builder);
}

bool KmerPresenceMatrix::operator==(const KmerPresenceMatrix& other) const {
    if (this->num_sequences != other.num_sequences || this->num_kmers != other.num_kmers) {
        return false;
    }
    if (this->arrays != other.arrays) {
        return false;
    }
    if (this->parents != other.parents) {
        return false;
    }
    if (this->children != other.children) {
        return false;
    }
    return true;
}

size_t KmerPresenceMatrix::get_num_present(size_t sequence) const {
    if (sequence >= this->num_sequences || this->num_kmers == 0) {
        return 0;
    }

    size_t rank = this->sequence_to_rank(sequence);
    size_t parent = this->parent_rank(rank);
    size_t parent_start = parent * this->num_kmers;
    size_t result = 0;

    // First count the number of present kmers in the parent.
    for (size_t kmer_id = 0; kmer_id < this->num_kmers; kmer_id += 64) {
        size_t bits = std::min(static_cast<size_t>(64), this->num_kmers - kmer_id);
        std::uint64_t word = this->parents.get_int(parent_start + kmer_id, bits);
        result += sdsl::bits::cnt(word);
    }

    // Then handle the bits that are different in the child.
    if (!this->is_parent(rank)) {
        ChildIterator child_iter(*this, rank);
        while (!child_iter.end()) {
            size_t kmer_id = child_iter.kmer_id();
            if (this->parents[parent_start + kmer_id]) {
                result -= 1;
            } else {
                result += 1;
            }
            ++child_iter;
        }
    }

    return result;
}

size_t KmerPresenceMatrix::get_total_present() const {
    size_t result = 0;
    size_t parent_start = 0;
    size_t parent_present = 0;

    // This duplicates the logic in `get_num_present` to avoid redundant work.
    for (size_t rank = 0; rank < this->num_sequences; rank++) {
        if (this->is_parent(rank)) {
            parent_start = this->parent_rank(rank) * this->num_kmers;
            parent_present = 0;
            for (size_t kmer_id = 0; kmer_id < this->num_kmers; kmer_id += 64) {
                size_t bits = std::min(static_cast<size_t>(64), this->num_kmers - kmer_id);
                std::uint64_t word = this->parents.get_int(parent_start + kmer_id, bits);
                parent_present += sdsl::bits::cnt(word);
            }
            result += parent_present;
        } else {
            result += parent_present;
            ChildIterator child_iter(*this, rank);
            while (!child_iter.end()) {
                size_t kmer_id = child_iter.kmer_id();
                if (this->parents[parent_start + kmer_id]) {
                    result -= 1;
                } else {
                    result += 1;
                }
                ++child_iter;
            }
        }
    }

    return result;
}

void KmerPresenceMatrix::for_each_kmer(size_t sequence, const std::function<void(size_t, bool)>& callback) const {
    Iterator iter(*this, sequence);
    while (!iter.end()) {
        callback(iter.kmer_id, iter.is_present);
        ++iter;
    }
}

void KmerPresenceMatrix::for_each_kmer(size_t first, size_t second, const std::function<void(size_t, size_t)>& callback) const {
    Iterator first_iter(*this, first);
    Iterator second_iter(*this, second);

    while (!first_iter.end()) {
        callback(first_iter.kmer_id, first_iter.is_present + second_iter.is_present);
        ++first_iter;
        ++second_iter;
    }
}

void KmerPresenceMatrix::score_sequences(
    const std::vector<std::pair<KmerPresenceMatrix::kmer_classification, double>>& kmers,
    const std::function<void(size_t, double)>& sequence_score
) const {
    size_t parent_start = 0;
    double parent_score = 0.0;

    for (size_t rank = 0; rank < this->num_sequences; rank++) {
        size_t sequence = this->rank_to_sequence(rank);
        if (this->is_parent(rank)) {
            parent_start = this->parent_rank(rank) * this->num_kmers;
            parent_score = 0.0;
            for (size_t kmer_id = 0; kmer_id < this->num_kmers; kmer_id++) {
                // +1 if the kmer is present and -1 if it is absent.
                double multiplier = -1.0 + 2.0 * this->parents[parent_start + kmer_id];
                parent_score += multiplier * kmers[kmer_id].second;
            }
            sequence_score(sequence, parent_score);
        } else {
            ChildIterator child_iter(*this, rank);
            double child_score = parent_score;
            while (!child_iter.end()) {
                size_t kmer_id = child_iter.kmer_id();
                auto was_present = this->parents[parent_start + kmer_id];
                // -2 if the kmer is present in the parent and +2 if it is absent.
                double multiplier = 2.0 - 4.0 * was_present;
                child_score += multiplier * kmers[kmer_id].second;
                ++child_iter;
            }
            sequence_score(sequence, child_score);
        }
    }
}

void KmerPresenceMatrix::simple_sds_serialize(std::ostream& out) const {
    // Since the dimensions are redundant, we do not store them.
    this->arrays.simple_sds_serialize(out);
    this->parents.simple_sds_serialize(out);
    this->children.simple_sds_serialize(out);
}

void KmerPresenceMatrix::simple_sds_load(std::istream& in) {
    this->arrays.simple_sds_load(in);
    this->parents.simple_sds_load(in);
    this->children.simple_sds_load(in);

    this->num_sequences = this->arrays.size() / 3;
    if (this->arrays.size() != this->num_sequences * 3) {
        throw sdsl::simple_sds::InvalidData("KmerPresenceMatrix: inconsistent array size");
    }

    if (this->num_sequences == 0) {
        this->num_kmers = 0;
    } else {
        this->num_kmers = (this->parents.size() + this->children.size()) / this->num_sequences;
    }
    if (this->parents.size() + this->children.size() != this->num_kmers * this->num_sequences) {
        throw sdsl::simple_sds::InvalidData("KmerPresenceMatrix: inconsistent matrix size");
    }

    // FIXME: is this too slow?
    // The queries use the arrays for indexing into `parents` and `children` without
    // further checks, so we validate the structure before trusting it.
    for (size_t rank = 0; rank < this->num_sequences; rank++) {
        size_t parent = this->parent_rank(rank);
        size_t previous = (rank == 0 ? 0 : this->parent_rank(rank - 1));
        if (rank == 0 ? parent != 0 : (parent != previous && parent != previous + 1)) {
            throw sdsl::simple_sds::InvalidData("KmerPresenceMatrix: invalid parent ranks");
        }
        size_t sequence = this->rank_to_sequence(rank);
        if (sequence >= this->num_sequences || this->sequence_to_rank(sequence) != rank) {
            throw sdsl::simple_sds::InvalidData("KmerPresenceMatrix: invalid sequence order");
        }
    }
    if (this->get_num_parents() * this->num_kmers != this->parents.size()) {
        throw sdsl::simple_sds::InvalidData("KmerPresenceMatrix: inconsistent number of parents");
    }
}

size_t KmerPresenceMatrix::simple_sds_size() const {
    return this->arrays.simple_sds_size() + this->parents.simple_sds_size() + this->children.simple_sds_size();
}

KmerPresenceMatrix::Iterator::Iterator(const KmerPresenceMatrix& matrix, size_t sequence) :
    matrix(matrix),
    parent_start(0), child_start(0), child_iter(matrix.children.one_end()),
    kmer_id(0), is_present(false)
{
    if (sequence >= this->matrix.num_sequences) {
        // Create an iterator that is already at the end.
        this->kmer_id = this->matrix.num_kmers;
        return;
    }

    size_t rank = this->matrix.sequence_to_rank(sequence);
    this->parent_start = this->matrix.parent_rank(rank) * this->matrix.num_kmers;
    if (!this->matrix.is_parent(rank)) {
        this->child_start = this->matrix.child_rank(rank) * this->matrix.num_kmers;
        this->child_iter = this->matrix.children.successor(this->child_start);
    }

    this->update_is_present();
}

void KmerPresenceMatrix::Iterator::update_is_present() {
    if (!this->end()) {
        this->is_present = this->matrix.parents[this->parent_start + this->kmer_id];
        if (this->child_iter != this->matrix.children.one_end() && this->child_iter->second == this->child_start + this->kmer_id) {
            this->is_present = !this->is_present;
            ++this->child_iter;
        }
    }
}

KmerPresenceMatrix::ChildIterator::ChildIterator(const KmerPresenceMatrix& matrix, size_t rank) :
    matrix(matrix)
{
    this->start = this->matrix.child_rank(rank) * this->matrix.num_kmers;
    this->iter = this->matrix.children.successor(this->start);
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
                std::vector<std::pair<ParallelKFFReader::kmer_type, size_t>> block = reader.read(KFF_BLOCK_SIZE);
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

// TODO: What is the right formula? Should subchain length be a parameter?
double Haplotypes::Subchain::badness(const gbwtgraph::GBZ& gbz) const {
    double result = 0.0;

    // Factor 1: Subchain length, ideally over a reference path.
    if (this->type == normal) {
        size_t selected = 0;
        for (size_t i = 0; i < this->sequences.size(); i++) {
            gbwt::size_type path_id = gbwt::Path::id(this->sequences[i].first);
            auto found = gbz.graph.id_to_path.find(path_id);
            if (found != gbz.graph.id_to_path.end()) {
                selected = i; break;
            }
        }
        size_t length = this->distance(gbz, selected);
        if (length > SUBCHAIN_LENGTH) {
            result += std::log(static_cast<double>(length) / static_cast<double>(SUBCHAIN_LENGTH));
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

// FIXME: also check individual dimensions
void load_subchain_kmers_present(Haplotypes::Subchain& subchain, std::istream& in) {
    subchain.kmers_present.simple_sds_load(in);
    if (subchain.kmers_present.size() != subchain.kmers.size() * subchain.sequences.size()) {
        throw sdsl::simple_sds::InvalidData("Invalid length for the kmer presence bitvector in subchain from " +
            std::to_string(subchain.start) + " to " + std::to_string(subchain.end));
    }
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
    this->tags.simple_sds_serialize(out);
    sdsl::simple_sds::serialize_vector(this->jobs_for_paths, out);
    for (auto& chain : this->chains) {
        chain.simple_sds_serialize(out);
    }
}

void Haplotypes::serialize_to(const std::string& filename) const {
    try {
        sdsl::simple_sds::serialize_to(*this, filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [Haplotypes] Serialization to " << filename << " failed: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
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

    // Tags are present from version 6 onwards.
    if (this->header.version >= Header::VERSION_WITH_TAGS) {
        this->tags.simple_sds_load(in);
    }

    this->jobs_for_paths = sdsl::simple_sds::load_vector<size_t>(in);

    this->chains.resize(this->header.top_level_chains);
    for (auto& chain : this->chains) {
        chain.simple_sds_load(in);
    }

    // Update to the current version.
    this->header.version = Header::VERSION;
}

void Haplotypes::load_from(const std::string& filename) {
    try {
        sdsl::simple_sds::load_from(*this, filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [Haplotypes] Loading from " << filename << " failed: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

size_t Haplotypes::simple_sds_size() const {
    size_t result = sdsl::simple_sds::value_size<Header>();
    result += sdsl::simple_sds::vector_size(this->jobs_for_paths);
    for (auto& chain : this->chains) {
        result += chain.simple_sds_size();
    }
    return result;
}

std::vector<size_t> Haplotypes::assign_reference_paths(const gbwtgraph::GBZ& gbz, Verbosity verbosity) const {
    if (verbosity >= verbosity_basic) {
        std::cerr << "Assigning reference paths to GBWT construction jobs" << std::endl;
    }
    double start = gbwt::readTimer();

    if (this->jobs_for_paths.size() != gbz.index.metadata.paths()) {
        std::string msg = "Haplotypes::assign_reference_paths(): Haplotype information was built for "
            + std::to_string(this->jobs_for_paths.size()) + " paths, but the graph contains "
            + std::to_string(gbz.index.metadata.paths()) + " paths";
        throw std::runtime_error(msg);
    }

    // All paths are initially unassigned.
    std::vector<size_t> result (gbz.named_paths(), this->jobs());
    size_t unassigned = 0;
    for (size_t i = 0; i < gbz.named_paths(); i++) {
        gbwt::size_type path_id = gbz.graph.named_paths[i].id;
        result[i] = this->jobs_for_paths[path_id];
        if (result[i] >= this->jobs()) {
            unassigned++;
        }
    }

    if (verbosity >= verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Assigned " << (result.size() - unassigned) << " reference paths (" << unassigned << " unassigned) in " << seconds << " seconds" << std::endl;
    }
    return result;
}

void Haplotypes::set_graph_name(const gbwtgraph::GBZ& gbz) {
    gbwtgraph::GraphName name = gbz.graph_name();
    name.set_tags(this->tags);
}

//------------------------------------------------------------------------------

} // namespace vg

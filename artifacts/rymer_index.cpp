#include "rymer_index.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <omp.h>

#include <gbwt/utils.h>
#include <gbwtgraph/gbwtgraph.h>

namespace vg {

using namespace std;

const string RYMER_INDEX_TAG = "rymer_index";
const string RYMER_INDEX_TAG_VALUE = "1";

gbwtgraph::Key64 rymer_key(gbwtgraph::Key64 dna_key, size_t k) {
    gbwtgraph::Key64::code_type dna = dna_key.get_key();
    gbwtgraph::Key64::code_type packed = 0;
    for (size_t i = 0; i < k; i++) {
        // The low bit of a base's two-bit code is 0 for a purine and 1 for a pyrimidine.
        size_t shift = (k - 1 - i) * gbwtgraph::KmerEncoding::PACK_WIDTH;
        packed = (packed << 1) | ((dna >> shift) & 1);
    }
    return gbwtgraph::Key64(packed);
}

gbwtgraph::DefaultMinimizerIndex::minimizer_type
to_rymer_minimizer(const gbwtgraph::DefaultMinimizerIndex::minimizer_type& dna_minimizer, size_t k) {
    gbwtgraph::DefaultMinimizerIndex::minimizer_type result = dna_minimizer;
    result.key = rymer_key(dna_minimizer.key, k);
    // The index holds no frequent k-mers, so the hash it uses for a key is the key's own hash.
    result.hash = result.key.hash();
    return result;
}

string rymer_string(const string& sequence) {
    string result;
    result.reserve(sequence.size());
    for (char c : sequence) {
        switch (c) {
        case 'A': case 'a': case 'G': case 'g':
            result.push_back('R');
            break;
        case 'C': case 'c': case 'T': case 't':
            result.push_back('Y');
            break;
        default:
            result.push_back('N');
        }
    }
    return result;
}

namespace {

/**
 * The loading dock of an RYmer index under construction.
 *
 * Each thread stacks its records here, and a full stack goes into the index at once, under a
 * lock. Threads then rarely wait for each other to insert.
 */
class RymerLoadingDock {
public:
    using minimizer_type = gbwtgraph::DefaultMinimizerIndex::minimizer_type;

    /// Makes a dock for the index, with a stack for each of `threads` OpenMP threads.
    RymerLoadingDock(gbwtgraph::DefaultMinimizerIndex& index, size_t threads) :
        index(index), stacks(threads) {
        // Nothing else to do
    }

    /// Stacks a DNA minimizer found at a graph position, to go into the index under its RYmer
    /// with its DNA key as payload. Exits the program if the index cannot store the position.
    void add(const minimizer_type& dna_minimizer, gbwtgraph::pos_t pos) {
        if (!gbwtgraph::Position::valid_offset(pos)) {
            #pragma omp critical (cerr)
            std::cerr << "error[index_rymers]: node offset " << gbwtgraph::offset(pos) << " is too large" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        size_t thread_id = omp_get_thread_num();
        stacks[thread_id].emplace_back(to_rymer_minimizer(dna_minimizer, index.k()), dna_minimizer.key.get_key(), pos);
        if (stacks[thread_id].size() >= STACK_SIZE) {
            load(thread_id);
        }
    }

    /// Puts every stacked record into the index. Call outside any parallel region.
    void load_all() {
        for (size_t thread_id = 0; thread_id < stacks.size(); thread_id++) {
            load(thread_id);
        }
    }

private:
    /// A minimizer to insert: its RYmer, its DNA key and its graph position.
    struct Record {
        minimizer_type minimizer;
        gbwtgraph::KmerEncoding::code_type dna_key;
        gbwtgraph::pos_t pos;
    };

    /// The number of records a thread stacks before it inserts them.
    static constexpr size_t STACK_SIZE = 1024;

    gbwtgraph::DefaultMinimizerIndex& index;
    vector<vector<Record>> stacks;

    /// Inserts the stack of one thread into the index and empties it.
    void load(size_t thread_id) {
        vector<Record>& stack = stacks[thread_id];
        #pragma omp critical (rymer_index)
        {
            for (Record& record : stack) {
                gbwtgraph::KmerEncoding::value_type value(gbwtgraph::Position(record.pos), &record.dna_key);
                index.insert(record.minimizer, value);
            }
        }
        stack.clear();
    }
};

}

void index_rymers(const gbwtgraph::GBZ& gbz, gbwtgraph::DefaultMinimizerIndex& index) {
    using minimizer_type = gbwtgraph::DefaultMinimizerIndex::minimizer_type;

    if (index.payload_size() != RYMER_PAYLOAD_SIZE) {
        throw std::invalid_argument("index_rymers(): an RYmer index needs a payload of "
                                    + std::to_string(RYMER_PAYLOAD_SIZE) + " word");
    }
    if (index.uses_syncmers()) {
        throw std::invalid_argument("index_rymers(): an RYmer index cannot use syncmers");
    }

    const gbwtgraph::GBWTGraph& graph = gbz.graph;
    index.set_graph_name(gbz.graph_name());
    index.set_tag(RYMER_INDEX_TAG, RYMER_INDEX_TAG_VALUE);
    int threads = omp_get_max_threads();
    RymerLoadingDock dock(index, threads);

    auto find_minimizers = [&graph, &index, &dock](const vector<gbwtgraph::handle_t>& traversal, const string& seq) {
        // The index has no frequent k-mers and no syncmers, so this chooses the same DNA
        // minimizers that a minimizer index with the same k and w would. They come in order of
        // offset, so we walk along the traversal once to find the node of each.
        auto iter = traversal.begin();
        size_t node_start = 0;
        for (const minimizer_type& minimizer : index.minimizers(seq)) {
            if (minimizer.empty()) {
                continue;
            }
            size_t node_length = graph.get_length(*iter);
            while (node_start + node_length <= minimizer.offset) {
                node_start += node_length;
                ++iter;
                node_length = graph.get_length(*iter);
            }
            gbwtgraph::pos_t pos = gbwtgraph::make_pos_t(graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start);
            dock.add(minimizer, minimizer.is_reverse ? gbwtgraph::reverse_base_pos(pos, node_length) : pos);
        }
    };

    gbwtgraph::for_each_haplotype_window(graph, index.window_bp(), find_minimizers, threads > 1);
    dock.load_all();
}

bool is_rymer_index(const gbwtgraph::DefaultMinimizerIndex& index) {
    return index.get_tag(RYMER_INDEX_TAG) == RYMER_INDEX_TAG_VALUE;
}

vector<gbwtgraph::Key64> graph_kmers_for_rymer(const gbwtgraph::DefaultMinimizerIndex& index,
                                               const gbwtgraph::DefaultMinimizerIndex::minimizer_type& rymer_minimizer) {
    vector<gbwtgraph::Key64> result;
    gbwtgraph::DefaultMinimizerIndex::multi_value_type hits = index.find(rymer_minimizer);
    result.reserve(hits.second);
    for (size_t i = 0; i < hits.second; i++) {
        auto [position, payload] = index.get_value(hits, i);
        result.emplace_back(payload[0]);
    }
    return result;
}

}

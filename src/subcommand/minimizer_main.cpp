/** \file minimizer_main.cpp
 *
 * Defines the "vg minimizer" subcommand, which builds the experimental
 * minimizer index.
 *
 * The index contains the lexicographically smallest kmer in a window of w
 * successive kmers. If the kmer contains characters other than A, C, G, and
 * T, it will not be indexed.
 *
 * By default, the index contains all minimizers in the graph. Option
 * --max-occs can be used to specify the maximum number of occurrences for
 * a kmer. Kmers more frequent than that will be removed from the index.
 *
 * The index contains either all minimizers or haplotype-consistent minimizers
 * (option --gbwt-name). Indexing all minimizers from complex graph regions
 * can take a long time (e.g. 65 hours vs 30 minutes for 1000GP), because many
 * windows have the same minimizer. As the total number of minimizers is
 * manageable (e.g. 2.1 billion vs. 1.4 billion for 1000GP), it should be
 * possible to develop a better algorithm for finding the minimizers.
 *
 * A quick idea:
 * - For each node v, extract the subgraph for the windows starting in v.
 * - Extract all k'-mers from the subgraph and use them to determine where the
 *   minimizers can start.
 */

#include "../gapless_extender.hpp"
#include "../gbwt_helper.hpp"
#include "../minimizer.hpp"
#include "subcommand.hpp"

#include "../stream/vpkg.hpp"
#include "../stream/stream.hpp"

#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>

#include <iostream>
#include <vector>

#include <getopt.h>
#include <omp.h>

using namespace vg;
using namespace vg::subcommand;


void help_minimizer(char** argv) {
    std::cerr << "usage: " << argv[0] << " minimizer [options] graph.xg" << std::endl;
    std::cerr << "Builds a minimizer index of the graph in the XG index." << std::endl;
    std::cerr << "    -k, --kmer-length N    length of the kmers in the index (default: " << MinimizerIndex::KMER_LENGTH << ")" << std::endl;
    std::cerr << "    -w, --window-length N  index the smallest kmer in a window of N kmers (default: " << MinimizerIndex::WINDOW_LENGTH << ")" << std::endl;
    std::cerr << "    -m, --max-occs N       do not index minimizers with more than N occurrences" << std::endl;
    std::cerr << "    -i, --index-name X     store the index to file X (required)" << std::endl;
    std::cerr << "    -l, --load-index X     load the index from file X and insert the new kmers into it" << std::endl;
    std::cerr << "                           (overrides --kmer-length, --window-length, and --max-occs)" << std::endl;
    std::cerr << "    -g, --gbwt-name X      index only haplotype-consistent kmers using the GBWT index in file X" << std::endl;
    std::cerr << "    -p, --progress         show progress information" << std::endl;
    std::cerr << "    -t, --threads N        use N threads for index construction (default: " << omp_get_max_threads() << ")" << std::endl;
    std::cerr << "benchmark options:" << std::endl;
    std::cerr << "    -b, --benchmark X      query performance benchmarks with the sequences in file X" << std::endl;
    std::cerr << "    -i, --index-name X     benchmark the minimizer index in file X (required)" << std::endl;
    std::cerr << "    -G, --gcsa-name X      also benchmark the GCSA2 index in file X" << std::endl;
    std::cerr << "    -L, --locate           locate the minimizer/MEM occurrences" << std::endl;
    std::cerr << "    -m, --max-occs N       do not locate minimizers with more than N occurrences" << std::endl;
    std::cerr << "    -e, --extend N         try gapless extension of unique hits with up to N mismatches (requires -g, -L)" << std::endl;
    std::cerr << "    -g, --gbwt-name X      extend using the GBWT index in file X" << std::endl;
    std::cerr << "    -M, --min-hits N       only extend reads with at least N minimizer hits" << std::endl;
}

int query_benchmarks(const std::unique_ptr<MinimizerIndex>& index, const std::unique_ptr<GBWTGraph>& gbwt_graph, const std::string& reads_name, const std::string& gcsa_name, bool locate, size_t max_occs, bool gapless_extend, size_t max_errors, size_t min_hits, bool progress);

int main_minimizer(int argc, char** argv) {

    if (argc == 2) {
        help_minimizer(argv);
        return 1;
    }

    // Command-line options.
    size_t kmer_length = MinimizerIndex::KMER_LENGTH;
    size_t window_length = MinimizerIndex::WINDOW_LENGTH;
    size_t max_occs = MinimizerIndex::MAX_OCCS;
    size_t max_errors = 0, min_hits = 1;
    std::string index_name, load_index, gbwt_name, xg_name, reads_name, gcsa_name;
    bool progress = false, locate = false, gapless_extend = false;
    int threads = omp_get_max_threads();

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            { "kmer-length", required_argument, 0, 'k' },
            { "window-length", required_argument, 0, 'w' },
            { "max-occs", required_argument, 0, 'm' },
            { "index-name", required_argument, 0, 'i' },
            { "load-index", required_argument, 0, 'l' },
            { "gbwt-name", required_argument, 0, 'g' },
            { "progress", no_argument, 0, 'p' },
            { "threads", required_argument, 0, 't' },
            { "benchmark", required_argument, 0, 'b' },
            { "gcsa-name", required_argument, 0, 'G' },
            { "locate", no_argument, 0, 'L' },
            { "extend", required_argument, 0, 'e' },
            { "min-hits", required_argument, 0, 'M' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "k:w:m:i:l:g:pt:hb:G:Le:M:", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'k':
            kmer_length = parse<size_t>(optarg);
            break;
        case 'w':
            window_length = parse<size_t>(optarg);
            break;
        case 'm':
            max_occs = parse<size_t>(optarg);
            break;
        case 'i':
            index_name = optarg;
            break;
        case 'l':
            load_index = optarg;
            break;
        case 'g':
            gbwt_name = optarg;
            break;
        case 'p':
            progress = true;
            break;
        case 't':
            threads = parse<int>(optarg);
            threads = std::min(threads, omp_get_max_threads());
            threads = std::max(threads, 1);
            omp_set_num_threads(threads);
            break;
        case 'b':
            reads_name = optarg;
            break;
        case 'G':
            gcsa_name = optarg;
            break;
        case 'L':
            locate = true;
            break;
        case 'e':
            max_errors = parse<size_t>(optarg);
            gapless_extend = true;
            break;
        case 'M':
            min_hits = std::max(static_cast<size_t>(1), parse<size_t>(optarg));
            break;

        case 'h':
        case '?':
            help_minimizer(argv);
            return 1;
        default:
            std::abort();
        }
    }
    if (index_name.empty()) {
        std::cerr << "[vg minimizer]: option --index-name is required" << std::endl;
        return 1;
    }
    if (gapless_extend && (!locate || gbwt_name.empty())) {
        std::cerr << "[vg minimizer]: option --extend requires --gbwt-name and --locate" << std::endl;
        return 1;
    }
    if (!reads_name.empty()) {
        load_index = index_name;
    }
    if (optind + 1 != argc) {
        help_minimizer(argv);
        return 1;
    }
    xg_name = argv[optind];

    double start = gbwt::readTimer();

    // Input graph.
    if (progress) {
        std::cerr << "Loading XG index " << xg_name << std::endl;
    }
    std::unique_ptr<xg::XG> xg_index;
    xg_index = stream::VPKG::load_one<xg::XG>(xg_name);

    // Minimizer index.
    std::unique_ptr<MinimizerIndex> index(new MinimizerIndex(kmer_length, window_length, max_occs));
    if (!load_index.empty()) {
        if (progress) {
            std::cerr << "Loading minimizer index " << load_index << std::endl;
        }
        index = stream::VPKG::load_one<MinimizerIndex>(load_index);
    }

    // GBWT-backed graph.
    std::unique_ptr<gbwt::GBWT> gbwt_index;
    std::unique_ptr<GBWTGraph> gbwt_graph;
    if (!gbwt_name.empty()) {
        if (progress) {
            std::cerr << "Loading GBWT index " << gbwt_name << std::endl;
        }
        gbwt_index = stream::VPKG::load_one<gbwt::GBWT>(gbwt_name);
        if (progress) {
            std::cerr << "Building GBWT-backed graph" << std::endl;
        }
        gbwt_graph.reset(new GBWTGraph(*gbwt_index, *xg_index));
        xg_index.reset(nullptr); // The XG index is no longer needed.
    }

    // Run the benchmarks and return.
    if (!reads_name.empty()) {
        return query_benchmarks(index, gbwt_graph, reads_name, gcsa_name, locate, max_occs, gapless_extend, max_errors, min_hits, progress);
    }

    // Minimizer caching.
    std::vector<std::vector<std::pair<MinimizerIndex::key_type, pos_t>>> cache(threads);
    constexpr size_t MINIMIZER_CACHE_SIZE = 1024;
    auto flush_cache = [&](int thread_id) {
        gbwt::removeDuplicates(cache[thread_id], false);
        #pragma omp critical (minimizer_index)
        {
            for (auto& minimizer : cache[thread_id]) {
                index->insert(minimizer.first, minimizer.second);
            }
        }
        cache[thread_id].clear();
    };

    // Build the index.
    if (progress) {
        std::cerr << "Building the index" << std::endl;
    }
    const HandleGraph& graph = *(gbwt_name.empty() ?
                                 static_cast<const HandleGraph*>(xg_index.get()) :
                                 static_cast<const HandleGraph*>(gbwt_graph.get()));
    auto lambda = [&](const std::vector<handle_t>& traversal, const std::string& seq) {
        std::vector<MinimizerIndex::minimizer_type> minimizers = index->minimizers(seq);
        auto iter = traversal.begin();
        size_t node_start = 0;
        int thread_id = omp_get_thread_num();
        for (MinimizerIndex::minimizer_type minimizer : minimizers) {
            if (minimizer.first == MinimizerIndex::NO_KEY) {
                continue;
            }
            // Find the node covering minimizer starting position.
            size_t node_length = graph.get_length(*iter);
            while (node_start + node_length <= minimizer.second) {
                node_start += node_length;
                ++iter;
                node_length = graph.get_length(*iter);
            }
            pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.second - node_start };
            cache[thread_id].emplace_back(minimizer.first, pos);
        }
        if (cache[thread_id].size() >= MINIMIZER_CACHE_SIZE) {
            flush_cache(thread_id);
        }
    };
    if (gbwt_name.empty()) {
        for_each_window(*xg_index, index->k() + index->w() - 1, lambda, (threads > 1));
    } else {
        for_each_haplotype_window(*gbwt_graph, index->k() + index->w() - 1, lambda, (threads > 1));
    }
    for (int thread_id = 0; thread_id < threads; thread_id++) {
        flush_cache(thread_id);
    }
    xg_index.reset(nullptr);
    gbwt_graph.reset(nullptr);
    gbwt_index.reset(nullptr);

    // Index statistics.
    if (progress) {
        std::cerr << index->size() << " keys (" << index->unique_keys() << " unique, " << index->frequent_keys() << " too frequent)" << std::endl;
        std::cerr << "Minimizer occurrences: " << index->values() << std::endl;
        std::cerr << "Load factor: " << index->load_factor() << std::endl;
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Construction so far: " << seconds << " seconds" << std::endl;
    }

    // Serialize the index.
    if (progress) {
        std::cerr << "Writing the index to " << index_name << std::endl;
    }
    stream::VPKG::save(*index, index_name);

    if (progress) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Time usage: " << seconds << " seconds" << std::endl;
        std::cerr << "Memory usage: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
    }

    return 0;
}

int query_benchmarks(const std::unique_ptr<MinimizerIndex>& index, const std::unique_ptr<GBWTGraph>& gbwt_graph, const std::string& reads_name, const std::string& gcsa_name, bool locate, size_t max_occs, bool gapless_extend, size_t max_errors, size_t min_hits, bool progress) {

    double start = gbwt::readTimer();

    // Load the GCSA index.
    bool benchmark_gcsa = !(gcsa_name.empty());
    std::unique_ptr<gcsa::GCSA> gcsa_index;
    std::unique_ptr<gcsa::LCPArray> lcp_index;
    if (benchmark_gcsa) {
        if (progress) {
            std::cerr << "Loading GCSA index " << gcsa_name << std::endl;
        }
        gcsa_index = stream::VPKG::load_one<gcsa::GCSA>(gcsa_name);
        std::string lcp_name = gcsa_name + gcsa::LCPArray::EXTENSION;
        lcp_index = stream::VPKG::load_one<gcsa::LCPArray>(lcp_name);
    }

    // Load the reads.
    std::vector<std::string> reads;
    if (progress) {
        std::cerr << "Loading the reads from " << reads_name << std::endl;
    }
    size_t total_size = gcsa::readRows(reads_name, reads, true);
    if (progress) {
        std::cerr << reads.size() << " reads of total length " << total_size << std::endl;
    }
    if (reads.empty()) {
        return 0;
    }

    size_t threads = omp_get_max_threads();
    if (progress) {
        std::cerr << "Query threads: " << threads << std::endl;
        std::cerr << std::endl;
    }

    // Minimizers.
    {
        double phase_start = gbwt::readTimer();

        GaplessExtender extender;
        if (gapless_extend) {
            extender = GaplessExtender(*gbwt_graph);
        }
        std::vector<size_t> min_counts(threads, 0);
        std::vector<size_t> occ_counts(threads, 0);
        std::vector<size_t> extend_counts(threads, 0);
        std::vector<size_t> seed_counts(threads, 0);
        std::vector<size_t> success_counts(threads, 0);
        std::vector<size_t> success_seed_counts(threads, 0);
        std::vector<size_t> partial_match_counts(threads, 0);
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < reads.size(); i++) {
            size_t thread = omp_get_thread_num();
            std::vector<MinimizerIndex::minimizer_type> result = index->minimizers(reads[i]);
            min_counts[thread] += result.size();
            if (locate) {
                std::vector<std::pair<size_t, pos_t>> hits; // (read offset, pos) for unique hits.
                for (auto minimizer : result) {
                    if (index->count(minimizer.first) <= max_occs) {
                        std::vector<pos_t> occs = index->find(minimizer.first);
                        if (occs.size() != 1 || !vg::is_empty(occs.front())) {
                            occ_counts[thread] += occs.size();
                            if (gapless_extend) {
                                for (pos_t occ : occs) {
                                    hits.emplace_back(minimizer.second, occ);
                                }
                            }
                        }
                    }
                }
                if (hits.size() >= min_hits) {
                    extend_counts[thread]++;
                    seed_counts[thread] += hits.size();
                    auto result = extender.extend_seeds(hits, reads[i], max_errors);
                    if (result.second <= max_errors) {
                        success_counts[thread]++;
                        success_seed_counts[thread] += hits.size();
                    } else {
                        auto partial_matches = extender.maximal_extensions(hits, reads[i]);
                        partial_match_counts[thread] += partial_matches.size();
                    }
                }
            } else {
                for (auto minimizer : result) {
                    occ_counts[thread] += index->count(minimizer.first);
                }
            }
        }
        size_t min_count = 0, occ_count = 0;
        size_t extend_count = 0, seed_count = 0, success_count = 0, success_seed_count = 0, partial_match_count = 0;
        for (size_t i = 0; i < threads; i++) {
            min_count += min_counts[i];
            occ_count += occ_counts[i];
            extend_count += extend_counts[i];
            seed_count += seed_counts[i];
            success_count += success_counts[i];
            success_seed_count += success_seed_counts[i];
            partial_match_count += partial_match_counts[i];
        }

        double phase_seconds = gbwt::readTimer() - phase_start;
        std::string query_type = "count";
        if (locate) { 
            query_type = "locate";
        }
        if (gapless_extend) {
            query_type = "extend";
        }
        std::cerr << "Minimizers (" << query_type << "): " << phase_seconds << " seconds (" << (reads.size() / phase_seconds) << " reads/second)" << std::endl;
        std::cerr << min_count << " minimizers with " << occ_count << " occurrences" << std::endl;
        if (gapless_extend) {
            std::cerr << extend_count << " reads with " << seed_count << " seeds, " << success_count << " extended with up to " << max_errors << " mismatches" << std::endl;
            std::cerr << "Extended " << (extend_count - success_count) << " reads with " << (seed_count - success_seed_count) << " seeds into " << partial_match_count << " partial matches" << std::endl;
        }
        std::cerr << std::endl;
    }

    // Count MEMs.
    if (benchmark_gcsa) {
        double phase_start = gbwt::readTimer();

        std::vector<size_t> mem_counts(threads, 0);
        std::vector<size_t> mem_lengths(threads, 0);
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < reads.size(); i++) {
            size_t thread = omp_get_thread_num();
            const std::string& read = reads[i];
            auto iter = read.rbegin();
            gcsa::range_type full(0, gcsa_index->size() - 1);
            gcsa::range_type curr = full;
            std::vector<gcsa::node_type> occurrences;
            while (iter != read.rend()) {
                gcsa::range_type prev = curr;
                curr = gcsa_index->LF(curr, gcsa_index->alpha.char2comp[*iter]);
                if (gcsa::Range::empty(curr)) {
                    if (prev != full) {
                        mem_counts[thread]++;
                        if (locate) {
                            gcsa_index->locate(prev, occurrences);
                            mem_lengths[thread] += occurrences.size();
                        } else {
                            mem_lengths[thread] += gcsa::Range::length(prev);
                        }
                        gcsa::STNode parent = lcp_index->parent(prev);
                        curr = parent.range();
                    } else {
                        curr = full;
                    }
                } else if (iter + 1 == read.rend()) {
                    if (prev != full) {
                        mem_counts[thread]++;
                        if (locate) {
                            gcsa_index->locate(prev, occurrences);
                            mem_lengths[thread] += occurrences.size();
                        } else {
                            mem_lengths[thread] += gcsa::Range::length(prev);
                        }
                    }
                }
                ++iter;
            }
        }
        size_t mem_count = 0, mem_length = 0;
        for (size_t i = 0; i < threads; i++) {
            mem_count += mem_counts[i];
            mem_length += mem_lengths[i];
        }

        double phase_seconds = gbwt::readTimer() - phase_start;
        std::string query_type = (locate ? "locate" : "count");
        std::cerr << "MEMs (" << query_type << "): " << phase_seconds << " seconds (" << (reads.size() / phase_seconds) << " reads/second)" << std::endl;
        std::cerr << mem_count << " MEMs with " << mem_length << " occurrences" << std::endl;
        std::cerr << std::endl;
    }


    // Locate minimizers
    // Locate MEMs

    double seconds = gbwt::readTimer() - start;
    std::cerr << "Benchmarks completed in " << seconds << " seconds" << std::endl;
    std::cerr << "Memory usage: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;

    return 0;
}

// FIXME change from DEVELOPMENT to TOOLKIT or PIPELINE later
// Register subcommand
static Subcommand vg_minimizer("minimizer", "build a minimizer index", DEVELOPMENT, main_minimizer);


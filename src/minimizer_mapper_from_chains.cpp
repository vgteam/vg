/**
 * \file minimizer_mapper_from_chains.cpp
 * Defines the code for the long-read code path for the
 * minimizer-and-GBWT-based mapper (long read Giraffe).
 */

#include "minimizer_mapper.hpp"

#include "annotation.hpp"
#include "crash.hpp"
#include "path_subgraph.hpp"
#include "multipath_alignment.hpp"
#include "split_strand_graph.hpp"
#include "subgraph.hpp"
#include "statistics.hpp"
#include "algorithms/alignment_path_offsets.hpp"
#include "algorithms/count_covered.hpp"
#include "algorithms/intersect_path_offsets.hpp"
#include "algorithms/extract_containing_graph.hpp"
#include "algorithms/extract_connecting_graph.hpp"
#include "algorithms/extract_extending_graph.hpp"
#include "algorithms/chain_items.hpp"
#include "algorithms/pad_band.hpp"

#include <bdsg/overlays/strand_split_overlay.hpp>
#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/cached_gbwtgraph.h>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cfloat>

// Turn on debugging prints
//#define debug
// Turn on printing of minimizer fact tables
//#define print_minimizer_table
// Dump the zip code forest
//#define debug_print_forest
// Dump local graphs that we align against 
//#define debug_dump_graph
// Dump fragment length distribution information
//#define debug_fragment_distr
//Do a brute force check that clusters are correct
//#define debug_validate_clusters
//#define debug_write_minimizers
// Debug generation of alignments from chains
#define debug_chain_alignment

namespace vg {

using namespace std;

/// Class for an error representing that chaining has backed us into some kind
/// of corner and we can't actually produce an alignment. We can throw this to
/// leave the read unmapped, complain, and try the next read.
class ChainAlignmentFailedError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

static void set_coverage_flags(std::vector<bool>& flags, size_t start, size_t end) {
    for (size_t i = start; i < end; i++) {
        flags[i] = true;
    }
}

static double get_fraction_covered(const std::vector<bool>& flags) {
    size_t covered_bases = 0;
    for (bool flag : flags) {
        if (flag) {
            covered_bases++;
        }
    }
    return (double) covered_bases / flags.size();
}

/// Get the forward-relative-to-the-read version of a seed's position. Will
/// have the correct orientation, but won't necessarily be to any particular
/// (i.e. first or last) base of the seed.
static pos_t forward_pos(const MinimizerMapper::Seed& seed, const VectorView<MinimizerMapper::Minimizer>& minimizers, const HandleGraph& graph) {
    pos_t position = seed.pos;
    if (minimizers[seed.source].value.is_reverse) {
        // Need to flip the position, for which we need to fetch the node length.
        position = reverse_base_pos(position, graph.get_length(graph.get_handle(id(position), is_rev(position))));
    }
    return position;
}

/// Figure out if the chains that start and end at the given seeds represent equivalent mappings
/// based on the range they cover in their top-level chain
static bool chain_ranges_are_equivalent(const MinimizerMapper::Seed& start_seed1, const MinimizerMapper::Seed& end_seed1,
                                        const MinimizerMapper::Seed& start_seed2, const MinimizerMapper::Seed& end_seed2) {
#ifdef debug
    assert(start_seed1.zipcode_decoder->get_distance_index_address(0) ==
           end_seed1.zipcode_decoder->get_distance_index_address(0));
    assert(start_seed2.zipcode_decoder->get_distance_index_address(0) ==
           end_seed2.zipcode_decoder->get_distance_index_address(0));
#endif
    if (start_seed1.zipcode_decoder->get_distance_index_address(0) !=
        start_seed2.zipcode_decoder->get_distance_index_address(0)) {
        //If the two ranges are on different connected components
        return false;
    }
    if (start_seed1.zipcode_decoder->get_code_type(0) == ZipCode::ROOT_SNARL) {
        //If this is in a root snarl
        if (start_seed1.zipcode_decoder->get_rank_in_snarl(1) !=
            start_seed2.zipcode_decoder->get_rank_in_snarl(1) 
            ||
            start_seed1.zipcode_decoder->get_rank_in_snarl(1) !=
            end_seed1.zipcode_decoder->get_rank_in_snarl(1) 
            ||
            start_seed2.zipcode_decoder->get_rank_in_snarl(1) !=
            end_seed2.zipcode_decoder->get_rank_in_snarl(1)) {
            //If the two ranges are on different children of the snarl
            return false;
        }
    }

    //Get the offset used for determining the range
    //On the top-level chain, node, or child of the top-level snarl 
    auto get_seed_offset = [&] (const MinimizerMapper::Seed& seed) {
        if (seed.zipcode_decoder->get_code_type(0) == ZipCode::ROOT_CHAIN) {
            return seed.zipcode_decoder->get_offset_in_chain(1);
        } else if (seed.zipcode_decoder->get_code_type(0) == ZipCode::ROOT_NODE) {
            return is_rev(seed.pos) ? seed.zipcode_decoder->get_length(0) - offset(seed.pos)
                                    : offset(seed.pos);
        } else {
            //Otherwise, this is a top-level snarl, and we've already made sure that it's on the 
            //same child chain/node
            if (seed.zipcode_decoder->get_code_type(1) == ZipCode::CHAIN) {
                //On a chain
                return seed.zipcode_decoder->get_offset_in_chain(2);
            } else {
                //On a node
                return is_rev(seed.pos) ? seed.zipcode_decoder->get_length(1) - offset(seed.pos)
                                        : offset(seed.pos);
            }
        }
    };
    size_t offset_start1 = get_seed_offset(start_seed1); 
    size_t offset_end1 = get_seed_offset(end_seed1);
    size_t offset_start2 = get_seed_offset(start_seed2); 
    size_t offset_end2 = get_seed_offset(end_seed2);

    if (offset_start1 > offset_end1) {
        size_t temp = offset_start1;
        offset_start1 = offset_end1;
        offset_end1 = temp;
    }
    if (offset_start2 > offset_end2) {
        size_t temp = offset_start2;
        offset_start2 = offset_end2;
        offset_end2 = temp;
    }

    if (offset_start1 > offset_end2 || offset_start2 > offset_end1 ){
        //If the ranges are disconnected
        return false;
    }if ( (offset_start1 <= offset_start2 && offset_end1 >= offset_end2) ||
         (offset_start2 <= offset_start1 && offset_end2 >= offset_end1)) {
        //If one range contains the other
        return true;
    } else {
        //Otherwise the two ranges must overlap on just one side

        if (offset_start1 > offset_start2) { 
            //Flip them so that range1 is first
            size_t tmp_start = offset_start1;
            size_t tmp_end = offset_end1;
            offset_start1 = offset_start2;
            offset_end1 = offset_end2;
            offset_start2 = tmp_start;
            offset_end2 = tmp_end;
        }

        size_t overlap_size = offset_end1 - offset_start2;
        //The two ranges count as equivalent if the length of the overlap is more than half the 
        //length of the shorter range
        return overlap_size > (std::min(offset_end1-offset_start1, offset_end2-offset_start2) / 2);

    }
}

void MinimizerMapper::dump_debug_dotplot(const std::string& name, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, const std::vector<std::pair<std::string, std::vector<std::vector<size_t>>>>& seed_sets, const PathPositionHandleGraph* path_graph) {
    if (!path_graph) {
        // We don't have a path positional graph for this
        return;
    }

    // Log the best bucket's seed positions in read and linear reference
    TSVExplainer exp(true, name + "-dotplot");

    // Determine the positions of all the involved seeds.
    std::unordered_map<size_t, algorithms::path_offset_collection_t> seed_positions;
    for (auto& kv : seed_sets) {
        for (const std::vector<size_t> included_seeds : kv.second) {
            for (auto& seed_num : included_seeds) {
                // For each seed in the run
                auto& seed = seeds.at(seed_num);

                auto found = seed_positions.find(seed_num);
                if (found == seed_positions.end()) {
                    // If we don't know the seed's positions yet, get them
                    seed_positions.emplace_hint(found, seed_num, algorithms::nearest_offsets_in_paths(path_graph, seed.pos, 100));
                }
            }
        }
    }
                
    for (auto& kv : seed_sets) {
        // For each named seed set
        const std::string& marker = kv.first;
        for (size_t run_number = 0; run_number < kv.second.size(); run_number++) {
            // For each run of seeds in it
            const std::vector<size_t>& included_seeds = kv.second[run_number];
            for (auto& seed_num : included_seeds) {
                // For each seed in the run
                auto& seed = seeds.at(seed_num);
                
                // Get its effective path positions
                auto& offsets = seed_positions.at(seed_num);

                for (auto& handle_and_positions : offsets) {
                    std::string path_name = path_graph->get_path_name(handle_and_positions.first);
                    for (auto& position : handle_and_positions.second) {
                        // For each position on a ref path that this seed is at, log a line
                        exp.line();
                        if (!marker.empty()) {
                            // Contig and a marker and a subscript
                            exp.field(path_name + "-" + marker + "-" + std::to_string(run_number));
                        } else {
                            // Contig alone
                            exp.field(path_name);
                        }
                        // Offset on contig
                        exp.field(position.first);
                        // Offset in read
                        exp.field(minimizers[seed.source].forward_offset());
                    }
                }
            }

        }
    }
}

void MinimizerMapper::dump_debug_graph(const HandleGraph& graph) {
    SubgraphExplainer exp(true);
    exp.subgraph(graph);
}

std::pair<double, double> MinimizerMapper::score_tree(const ZipCodeForest& zip_code_forest, size_t i, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t seq_length, Funnel& funnel) const {
    // Initialize the values.
    std::pair<double, double> to_return;
    auto& score = to_return.first;
    auto& coverage = to_return.second;
    
    // Start score at 0.
    score = 0;
    // Coverage gets set all at once.

    // Track if minimizers are present
    SmallBitset present(minimizers.size());
    // And if read bases are covered
    sdsl::bit_vector covered(seq_length, 0);

    vector<size_t> tree_seeds;
    for (ZipCodeTree::oriented_seed_t found : zip_code_forest.trees[i]) {
        if (this->track_provenance) {
            // Remember the seeds
            tree_seeds.push_back(found.seed);
        }
        // For each seed in the tree, find what minimizer it comes from
        if (found.seed >= seeds.size()) {
            throw std::out_of_range("Tree " + std::to_string(i) + " has seed " + std::to_string(found.seed) + " but we only have " + std::to_string(seeds.size()) + " seeds");
        }
        size_t source = seeds.at(found.seed).source;
        if (!present.contains(source)) {
            // If it's a new minimizer, count its score
            score += minimizers[source].score;

            // Mark its read bases covered.
            // The offset of a reverse minimizer is the endpoint of the kmer
            size_t start_offset = minimizers[source].forward_offset();
            size_t k = minimizers[source].length;

            // Set the k bits starting at start_offset.
            covered.set_int(start_offset, sdsl::bits::lo_set[k], k);

            // Mark it present
            present.insert(source);
        }
    }

    // Count up the covered positions and turn it into a fraction.
    coverage = sdsl::util::cnt_one_bits(covered) / static_cast<double>(seq_length);

    if (this->track_provenance) {
        // Record the tree in the funnel as a group of the size of the number of items.
        funnel.merge_group(tree_seeds.begin(), tree_seeds.end());
        funnel.score(funnel.latest(), score);

        // TODO: Should we tell the funnel we produced an output?

        if (show_work && track_correctness) {
            // We will have positions early, for all the seeds.
            auto tree_positions = funnel.get_positions(funnel.latest());
            #pragma omp critical (cerr)
            {
                std::cerr << log_name() << "Positions for tree " << i << " score " << score << " coverage " << coverage << ":" << std::endl;
                for (auto& handle_and_range : tree_positions) {
                    // Log each range on a path associated with the tree.
                    std::cerr << log_name() << "\t"
                        << this->path_graph->get_path_name(handle_and_range.first)
                        << ":" << handle_and_range.second.first
                        << "-" << handle_and_range.second.second << std::endl;
                }
                if (track_correctness && funnel.is_correct(funnel.latest())) {
                    cerr << log_name() << "\t\tCORRECT!" << endl;
                }
            }
        }
    }

    return to_return;
}

/**
 * Given a read interval for a gapless extension, the read positions of
 * mismatches, and the read positions of seeds, compute anchor intervals.
 *
 * Inputs and outputs are all sorted.
 *
 * Anchor intervals do not overlap.
 *
 * There will be at least one seed in each anchor interval.
 * 
 * Anchor intervals will begin and end at the bounds of the read interval, or
 * just outside mismatches.
 *
 * Anchor intervals will not go over logn runs of mismatches that give them
 * deceptively terrible scores.
 */
std::vector<std::pair<size_t, size_t>> find_anchor_intervals(
    const std::pair<size_t, size_t>& read_interval,
    const std::vector<size_t>& mismatch_positions,
    const std::vector<size_t>& seed_positions) {

    assert(!seed_positions.empty());

    std::vector<std::pair<size_t, size_t>> anchor_intervals;

    if (mismatch_positions.empty()) {
        // Everything will form one giant anchor and there will be no
        // mismatches to key on being after. So just handle it here.
        anchor_intervals.push_back(read_interval);
        return anchor_intervals;
    }


    // We are going to sweep line.
    auto mismatch_it = mismatch_positions.begin();
    auto seed_it = seed_positions.begin();

    // We need to track:
    // The previous seed.
    auto prev_seed = seed_positions.end();
    // The first mismatch we saw after the previous seed.
    auto mismatch_after_prev_seed = mismatch_positions.end();
    // The last mismatch we saw before the current seed.
    auto mismatch_before_current_seed = mismatch_positions.end();

    size_t interval_start = read_interval.first;

    auto visit_seed = [&]() {
#ifdef debug_anchor_intervals
        if (seed_it != seed_positions.end()) {
            std::cerr << "Visit seed at " << *seed_it << std::endl;
        } else {
            std::cerr << "Visit fake final seed" << std::endl;
        }
#endif
        
        // Process the seed at seed_it (which may be the end), which comes next.
        if (prev_seed == seed_positions.end()) {
            // This is the first seed, so we need to trim from the left end of the read.
#ifdef debug_anchor_intervals
            std::cerr << "This is the first seed" << std::endl;
#endif
            assert(seed_it != seed_positions.end());
            int score = 0;
            auto here = mismatch_before_current_seed;
            int max_score = score;
            auto max_cut = here;
            if (here != mismatch_positions.end()) {
                // There are mismatches to score 
                while (here != mismatch_positions.begin()) {
                    auto next = here;
                    --next;
                    // Score taking that mismatch and then going up to the next one
                    size_t matches = *here - *next - 1;
                    score += matches;
                    score -= 4; // TODO: use real scoring
                    if (score > max_score) {
                        max_score = score;
                        max_cut = next;
                    }
                    here = next;
                }
                // Now we're at the first mismatch, so score from there to the bound of the read interval.
                size_t matches = *here - read_interval.first;
                score += matches;
                score -= 4; // TODO: use real scoring
                if (score > max_score) {
                    max_score = score;
                    // Use end to represent going all the way to the read bound
                    max_cut = mismatch_positions.end();
                }
            }
            if (max_cut != mismatch_positions.end()) {
                // Trim the anchor interval start
                interval_start = *max_cut + 1;
            }
            // Otherwise leave the anchor interval start at the read interval start.
#ifdef debug_anchor_intervals
            std::cerr << "First seed interval should start at " << interval_start << std::endl;
#endif
        } else if (mismatch_after_prev_seed != mismatch_positions.end()) {
            // This is the first seed after some mismatches (or we did all the seeds and mismatches)
            assert(mismatch_before_current_seed != mismatch_positions.end());

#ifdef debug_anchor_intervals
            std::cerr << "Mismatch after previous seed was at " << *mismatch_after_prev_seed << std::endl;
            std::cerr << "Mismatch before current seed was at " << *mismatch_before_current_seed << std::endl;
#endif

            // So we have to finish off the last seed's interval.

            std::vector<size_t>::const_iterator split_mismatch;
            if (seed_it != seed_positions.end()) {
                // Pick a middle mismatch to divide the two intervals with initially.
                size_t separating_mismatches = mismatch_before_current_seed - mismatch_after_prev_seed  + 1;
                size_t middle_offset = separating_mismatches / 2;
                // TODO: Feed in information that would let us round in a
                // consistent direction even if we flip the read.
                split_mismatch = mismatch_after_prev_seed + middle_offset;
            } else {
                // Do the split at the past-end mismatch
                split_mismatch = mismatch_positions.end();
            }

            // Trim left for the old seed's interval.
            //
            // Starting at mismatch_after_prev_seed and going right to
            // split_mismatch, get the score we have taking up to just before
            // each mismatch, and the mismatch we cut at to get it.
            int score = 0;
            auto here = mismatch_after_prev_seed;
            int max_score = score;
            auto max_cut = here;
            while (here != split_mismatch) {
                auto next = here;
                ++next;
                // Score taking that mismatch and then going up to the next one
                size_t matches = (next == mismatch_positions.end() ? read_interval.second : *next) - *here - 1;
                score += matches;
                score -= 4; // TODO: use real scoring
                if (score > max_score) {
                    max_score = score;
                    max_cut = next;
                }
                here = next;
            }
            auto left_separating_mismatch = max_cut;
            size_t interval_end = (left_separating_mismatch == mismatch_positions.end() ? read_interval.second : *left_separating_mismatch);
#ifdef debug_anchor_intervals
            std::cerr << "Previous seed interval should end at " << interval_end << std::endl;
#endif
            // So that's where the old interval ends.
            anchor_intervals.emplace_back(interval_start, interval_end);
            
            if (seed_it != seed_positions.end()) {
                // Trim right for the new seed's interval.
                //
                // Starting at mismatch_before_current_seed and going left to
                // split_mismatch, get the score we have taking up to just before
                // each mismatch, and the mismatch we cut at to get it.
                score = 0;
                here = mismatch_before_current_seed;
                max_score = score;
                max_cut = here;
                while (here != split_mismatch) {
                    auto next = here;
                    --next;
                    // Score taking that mismatch and then going up to the next one
                    size_t matches = *here - *next - 1;
                    score += matches;
                    score -= 4; // TODO: use real scoring
                    if (score > max_score) {
                        max_score = score;
                        max_cut = next;
                    }
                    here = next;
                }
                auto right_separating_mismatch = max_cut;
                // And after it is where our interval starts.
                interval_start = *right_separating_mismatch + 1;
#ifdef debug_anchor_intervals
                std::cerr << "Current seed interval should start at " << interval_start << std::endl;
#endif
            }
        } else if (seed_it == seed_positions.end()) {
            // We ran out of seeds and there are no mismatches between the last seed and the itnerval end.
            // TODO: Combine with above case?
            size_t interval_end =read_interval.second;
#ifdef debug_anchor_intervals
            std::cerr << "Previous seed interval should end at end of extension at " << interval_end << std::endl;
#endif
            // So that's where the old interval ends.
            anchor_intervals.emplace_back(interval_start, interval_end);
        }

        // Now this seed is the previous seed.
        prev_seed = seed_it;
        // And no mismatch has been seen after it yet.
        mismatch_after_prev_seed = mismatch_positions.end();
    };

    auto visit_mismatch = [&]() {
        // Process the mismatch at mismatch_it (which is not the end), which comes next.
#ifdef debug_anchor_intervals
        std::cerr << "Visit mismatch at " << *mismatch_it << std::endl;
#endif

        if (prev_seed != seed_positions.end() && mismatch_after_prev_seed == mismatch_positions.end()) {
            // This is the first mismatch since we saw a seed, so save it.
            mismatch_after_prev_seed = mismatch_it;
        }
        // This is now the last mismatch we've seen.
        mismatch_before_current_seed = mismatch_it;
    };

    while (mismatch_it != mismatch_positions.end() && seed_it != seed_positions.end()) {
        if (*mismatch_it < *seed_it) {
            // Next is a mismatch
            visit_mismatch();
            ++mismatch_it;
        } else {
            // Next is a seed
            visit_seed();
            ++seed_it;
        }
    }
    while (mismatch_it != mismatch_positions.end()) {
        // Next is a mismatch
        visit_mismatch();
        ++mismatch_it;
    }
    while (seed_it != seed_positions.end()) {
        // Next is a seed
        visit_seed();
        ++seed_it;
    }
    // Visit the end seed to finish off the last interval
    visit_seed();

    assert(!anchor_intervals.empty());

    return anchor_intervals;
}

vector<Alignment> MinimizerMapper::map_from_chains(Alignment& aln) {

    //Do gapless extension if the read length is less than the limit
    bool do_gapless_extension = aln.sequence().size() <= gapless_extension_limit;

    
    if (show_work) {
        #pragma omp critical (cerr)
        dump_debug_query(aln);
    }
    
    // Make a new funnel instrumenter to watch us map this read.
    Funnel funnel;
    funnel.start(aln.name());
    
    // Prepare the RNG for shuffling ties, if needed
    LazyRNG rng([&]() {
        return aln.sequence();
    });


    // Minimizers sorted by position
    std::vector<Minimizer> minimizers_in_read = this->find_minimizers(aln.sequence(), funnel);
    // Indexes of minimizers, sorted into score order, best score first
    std::vector<size_t> minimizer_score_order = sort_minimizers_by_score(minimizers_in_read);
    // Minimizers sorted by best score first
    VectorView<Minimizer> minimizers{minimizers_in_read, minimizer_score_order};


    // Find the seeds and mark the minimizers that were located.
    vector<Seed> seeds = this->find_seeds(minimizers_in_read, minimizers, aln, funnel);

    if (seeds.empty()) {
        #pragma omp critical (cerr)
        std::cerr << log_name() << "warning[MinimizerMapper::map_from_chains]: No seeds found for " << aln.name() << "!" << std::endl;
    }
    
    if (this->track_provenance) {
        funnel.stage("tree");
    }

    // Make them into a zip code tree
    ZipCodeForest zip_code_forest;
    crash_unless(distance_index);
    zip_code_forest.fill_in_forest(seeds, minimizers, *distance_index, 
                                   max_lookback_bases, aln.sequence().size() * zipcode_tree_scale);

#ifdef debug_print_forest
    if (show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "Zip code forest:";
            zip_code_forest.print_self(&seeds, &minimizers);
        }
    }
#endif

    // Now score all the zip code trees in the forest by summing the scores of their involved minimizers.
    vector<double> tree_scores;
    double best_tree_score = 0;
    double second_best_tree_score = 0;
    tree_scores.reserve(zip_code_forest.trees.size());

    vector<double> tree_coverages;
    double best_tree_coverage = 0;
    double second_best_tree_coverage = 0;
    tree_coverages.reserve(zip_code_forest.trees.size());

    for (size_t i = 0; i < zip_code_forest.trees.size(); i++) {
        // For each zip code tree
        
        // Score it
        std::pair<double, double> metrics = this->score_tree(zip_code_forest, i, minimizers, seeds, aln.sequence().size(), funnel);
        auto& score = metrics.first;
        auto& coverage = metrics.second;

        tree_scores.push_back(score);
        tree_coverages.push_back(coverage);

        if (score > best_tree_score) {
            second_best_tree_score = best_tree_score;
            best_tree_score = score;
        } else if (score > second_best_tree_score) {
            second_best_tree_score = score;
        }

        if (coverage > best_tree_coverage) {
            second_best_tree_coverage = best_tree_coverage;
            best_tree_coverage = coverage;
        } else if (coverage > second_best_tree_coverage) {
            second_best_tree_coverage = coverage;
        }
    }

    // We will set a score cutoff based on the best, but move it down to the
    // second best if it does not include the second best and the second best
    // is within pad_zipcode_tree_score_threshold of where the cutoff would
    // otherwise be. This ensures that we won't throw away all but one 
    // based on score alone, unless it is really bad.
    double tree_score_cutoff = best_tree_score - zipcode_tree_score_threshold;
    if (tree_score_cutoff - pad_zipcode_tree_score_threshold < second_best_tree_score) {
        tree_score_cutoff = std::min(tree_score_cutoff, second_best_tree_score);
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "Found " << zip_code_forest.trees.size() << " zip code trees, scores " << best_tree_score << " best, " << second_best_tree_score << " second best, coverages " << best_tree_coverage << " best, " << second_best_tree_coverage << " second best" << std::endl;
        }
    }

    // Turn all the seeds into anchors. Either we'll fragment them directly or
    // use them to make gapless extension anchors over them.
    // TODO: Can we only use the seeds that are in trees we keep?
    vector<algorithms::Anchor> seed_anchors = this->to_anchors(aln, minimizers, seeds);

    // If we don't do gapless extension, we need one-item vectors for all the
    // seeds of their own numbers, to show what seed each anchor represents.
    // TODO: Can we only do this for the seeds that are in trees we keep?
    std::vector<std::vector<size_t>> seed_seed_sequences;
    if (!do_gapless_extension) {
        seed_seed_sequences.reserve(seed_anchors.size());
        for (size_t i = 0; i < seed_anchors.size(); ++i) {
            seed_seed_sequences.push_back({i});
        }
    }
    // Now we need to chain into fragments.
    // Each fragment needs to end up with a seeds array of seed numbers, and a
    // coverage float on the read, for downstream
    // processing.
    if (track_provenance) {
        funnel.stage("fragment");
        funnel.substage("fragment");
    }
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "=====Creating fragments=====" << endl;
        }
    }

    // Now compute fragments into these variables.
    // What seeds are visited in what order in the fragment?
    std::vector<std::vector<size_t>> fragments;
    // What score does each fragment have?
    std::vector<double> fragment_scores;
    // What are the fragments themselves as combined anchors, for chaining later?
    std::vector<algorithms::Anchor> fragment_anchors;
    // Which zip code tree did each fragment come from, so we know how to chain them?
    std::vector<size_t> fragment_source_tree;
    // How many of each minimizer ought to be considered explored by each fragment?
    // TODO: This is a lot of counts and a lot of allocations and should maybe be a 2D array if we really need it?
    std::vector<std::vector<size_t>> minimizer_kept_fragment_count;
    // For capping mapq, we want the multiplicity of each alignment. Start keeping track of this
    // here with the multiplicity of the trees for each fragment
    // For now, this just stores how many trees had equal or better score. After going through all
    // trees and counting how many are kept, each value will be divided by the number of trees kept
    std::vector<double> multiplicity_by_fragment;
    size_t kept_tree_count = 0;

    process_until_threshold_c<double>(zip_code_forest.trees.size(), [&](size_t i) -> double {
            return tree_coverages[i];
        }, [&](size_t a, size_t b) -> bool {
            return tree_coverages[a] > tree_coverages[b] || (tree_coverages[a] == tree_coverages[b] && tree_scores[a] > tree_scores[b]); 
        }, zipcode_tree_coverage_threshold, this->min_to_fragment, this->max_to_fragment, rng, [&](size_t item_num, size_t item_count) -> bool {
            // Handle sufficiently good fragmenting problems in descending score order
            
            if (track_provenance) {
                funnel.pass("zipcode-tree-coverage-threshold", item_num, tree_coverages[item_num]);
                funnel.pass("max-to-fragment", item_num);
            }

            // First check against the additional score filter
            if (zipcode_tree_score_threshold != 0 && tree_scores[item_num] < tree_score_cutoff 
                && kept_tree_count >= min_to_fragment) {
                // If the score isn't good enough and we already kept at least min_to_fragment trees,
                // ignore this tree
                if (track_provenance) {
                    funnel.fail("zipcode-tree-score-threshold", item_num, tree_scores[item_num]);
                }
                return false;
            }
            
            if (track_provenance) {
                funnel.pass("zipcode-tree-score-threshold", item_num, tree_scores[item_num]); 
            }

            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Making fragments for zip code tree " << item_num << " with score " << tree_scores[item_num] << " and coverage " << tree_coverages[item_num] << endl;
                }
            }
            
            kept_tree_count++;

            if (track_provenance) {
                // Say we're working on this 
                funnel.processing_input(item_num);
            }
          
            // Also make a list of all the seeds in the problem.
            // This lets us select the single-seed anchors to use.

            //Make sure that each seed gets added only once
            vector<bool> added_seed (seeds.size(), false);
            vector<size_t> selected_seeds;
            for (ZipCodeTree::oriented_seed_t found : zip_code_forest.trees[item_num]) {
                if (!added_seed[found.seed]) {
                    selected_seeds.push_back(found.seed);
                    added_seed[found.seed] = true;
                }
            }
            
            if (show_work) {
                dump_debug_seeds(minimizers, seeds, selected_seeds);
            }

            // If we do gapless extension, we will use these anchors to fragment instead of the seed ones.
            std::vector<algorithms::Anchor> extension_anchors;
            // And each of them (or of the seed anchors, if we use those) represents this run of seed numbers to put into the final chain.
            std::vector<std::vector<size_t>> extension_seed_sequences;
            // Extensions use a distinct list of included seeds vs. seeds we actually paste in, so we can glom up overlapping seeds.
            std::vector<std::vector<size_t>> extension_represented_seeds;
            // We need a list of all extension anchor indexes that we can sort.
            std::vector<size_t> extension_anchor_indexes;

            if (do_gapless_extension) {
                // Instead of fragmenting directly on the seeds, fragment on gapless extensions of the seeds.

                if (track_provenance) {
                    funnel.substage("gapless_extension");
                }

                // Extend the seeds and keep track of the seeds that went into each extension.
                // We'll use this to make anchors later.
                std::vector<std::vector<size_t>> seeds_for_extension;
                std::vector<GaplessExtension> tree_extensions = this->extend_seed_group(
                    selected_seeds,
                    item_num,
                    minimizers,
                    seeds,
                    aln.sequence(),
                    this->max_extension_mismatches,
                    nullptr,
                    nullptr,
                    &seeds_for_extension);
                // Note that we don't use the funnel here; we don't actually
                // track a gapless extension stage.
                
                // We can't actually handle the same seed being used as the
                // endpoint of multiple anchors in the chaining. So we need to
                // go through the gapless extensions in score order and make
                // them into anchors using the seeds not yet used by previous
                // ones.
                auto extension_score_order = sort_permutation(tree_extensions.begin(), tree_extensions.end(), [&](const GaplessExtension& a, const GaplessExtension& b) {
                    // Return true if the first gapless extension needs to be first.
                    // TODO: use real scores from the aligner.
                    int a_score = (a.read_interval.second - a.read_interval.first) - a.mismatch_positions.size() * 5;
                    int b_score = (b.read_interval.second - b.read_interval.first) - b.mismatch_positions.size() * 5;
                    // We want to sort descending so larger scores come first.
                    return a_score > b_score;
                });

                // This holds the seeds used to make previous anchors.
                std::unordered_set<size_t> used_seeds;

                for (auto& extension_index : extension_score_order) {
                    // For each extension
                    const GaplessExtension& extension = tree_extensions[extension_index];
                    // And the seeds that made it, sorted by stapled base
                    const std::vector<size_t>& extension_seeds = seeds_for_extension[extension_index];

                    // Make a list of all the seed positions still available
                    std::vector<size_t> seed_positions;
                    seed_positions.reserve(extension_seeds.size());
                    for (auto& seed_index : extension_seeds) {
                        if (!used_seeds.count(seed_index)) {
                            seed_positions.push_back(minimizers[seeds.at(seed_index).source].value.offset);
                        }
                    }

                    if (seed_positions.empty()) {
                        if (show_work) {
                            #pragma omp critical (cerr)
                            {
                                cerr << log_name() << "Extension on read " << extension.read_interval.first << "-" << extension.read_interval.second << " has no distinct seeds left to use for anchors" << endl;
                            }
                        }
                        continue;
                    }


                    // We want to break up the extension into read intervals
                    // and the seeds that go with them. Each of those will
                    // become an anchor.
                    std::vector<std::pair<size_t, size_t>> anchor_intervals = find_anchor_intervals(extension.read_interval, extension.mismatch_positions, seed_positions);

                    // Then convert those intervals into anchors.
                    auto mismatch_it = extension.mismatch_positions.begin();
                    auto seed_it = extension_seeds.begin();
                    for (auto& anchor_interval : anchor_intervals) {
                        // Find the relevant mismatch range
                        while (mismatch_it != extension.mismatch_positions.end() && *mismatch_it < anchor_interval.first) {
                            // Move mismatch iterator to inside or past the interval
                            ++mismatch_it;
                        }
                        auto internal_mismatch_begin = mismatch_it;
                        while (mismatch_it != extension.mismatch_positions.end() && *mismatch_it < anchor_interval.second) {
                            // Move mismatch iterator to past the interval
                            ++mismatch_it;
                        }
                        auto internal_mismatch_end = mismatch_it;

                        // Find the relevant seed range
                        std::vector<size_t> anchor_seeds;
                        while (seed_it != extension_seeds.end() && minimizers[seeds.at(*seed_it).source].value.offset < anchor_interval.first) {
                            // Move seed iterator to inside or past the interval (should really always be already inside).
                            ++seed_it;
                        }
                        while (seed_it != extension_seeds.end() && minimizers[seeds.at(*seed_it).source].value.offset < anchor_interval.second) {
                            // Take all the seeds into the vector of anchor seeds.
                            auto found = used_seeds.find(*seed_it);
                            if (found == used_seeds.end()) {
                                // As long as they haven't been used
                                anchor_seeds.push_back(*seed_it);
                                // And mark them used
                                used_seeds.insert(found, *seed_it);
                            }
                            ++seed_it;
                        }

                        if (anchor_seeds.empty()) {
                            // All the seeds we wanted for this piece specifically are already represented by pieces of previous extensions
                            if (show_work) {
                                #pragma omp critical (cerr)
                                {
                                    cerr << log_name() << "Extension on read " << extension.read_interval.first << "-" << extension.read_interval.second << " would produce anchor " << anchor_interval.first << "-" << anchor_interval.second << " but all seeds in the interval were used already" << endl;
                                }
                            }
                            // Go on to the next anchor interval
                        } else {
                            // We have seeds here and can make an anchor

                            // Note the index of the new anchor
                            extension_anchor_indexes.push_back(extension_anchors.size());
                            // Make the actual anchor out of this range of seeds and this read range.
                            extension_anchors.push_back(to_anchor(aln, anchor_interval.first, anchor_interval.second, anchor_seeds, seed_anchors, internal_mismatch_begin, internal_mismatch_end, gbwt_graph, this->get_regular_aligner()));
                            if (show_work) {
                                #pragma omp critical (cerr)
                                {
                                    cerr << log_name() << "Extension on read " << extension.read_interval.first << "-" << extension.read_interval.second << " produces anchor " << anchor_interval.first << "-" << anchor_interval.second << " with " << anchor_seeds.size() << " seeds involved and " << (internal_mismatch_end - internal_mismatch_begin) << " internal mismatches, score " << extension_anchors.back().score() << endl;
                                }
                            }

                            // And if we take that anchor, we'll grab these underlying
                            // seeds into the elaborating chain. Just use the bounding
                            // seeds and connect between them where it is easy.
                            extension_seed_sequences.push_back({anchor_seeds.front()});
                            if (seed_anchors.at(anchor_seeds.front()).read_end() <= seed_anchors.at(anchor_seeds.back()).read_start()) {
                                // There are multiple seeds in the extension and the last
                                // one doesn't overlap the first, so take the last one too.
                                extension_seed_sequences.back().push_back(anchor_seeds.back());
                            }

                            // Keep all the seeds that this anchor counts as using.
                            extension_represented_seeds.emplace_back(std::move(anchor_seeds));
                        }
                    }
                }
            }
            
            // Figure out what anchors we want to view.
            const std::vector<algorithms::Anchor>& anchors_to_fragment = do_gapless_extension ? extension_anchors : seed_anchors;
            // And what seeds each represents
            const std::vector<std::vector<size_t>>& anchor_seed_sequences = do_gapless_extension ? extension_seed_sequences : seed_seed_sequences;
            // And what subset/in what order
            std::vector<size_t>& anchor_indexes = do_gapless_extension ? extension_anchor_indexes : selected_seeds;
            // Sort anchors by read start of seeded region
            algorithms::sort_anchor_indexes(anchors_to_fragment, anchor_indexes);

            // And what seeds should count as explored when we take an anchor
            const std::vector<std::vector<size_t>>& anchor_represented_seeds = do_gapless_extension ? extension_represented_seeds : anchor_seed_sequences;
            
            

            if (track_provenance) {
                funnel.substage("fragment");
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Computing fragments over " << anchor_indexes.size() << " anchors" << endl;
                }
            }

#ifdef debug
            if (show_work) {
                // Log the chaining problem so we can try it again elsewhere.
                this->dump_chaining_problem(anchors_to_fragment, anchor_indexes, gbwt_graph);
            }
#endif
            
            // Compute lookback and indel limits based on read length.
            // Important since seed density goes down on longer reads.
            size_t lookback_limit = std::max(this->fragment_max_lookback_bases, (size_t)(this->fragment_max_lookback_bases_per_base * aln.sequence().size()));
            size_t indel_limit = std::max(this->fragment_max_indel_bases, (size_t)(this->fragment_max_indel_bases_per_base * aln.sequence().size()));

            // Find fragments over the seeds in the zip code tree
            algorithms::transition_iterator for_each_transition = algorithms::zip_tree_transition_iterator(
                seeds,
                zip_code_forest.trees[item_num],
                lookback_limit
            ); 
            // Make a view of the anchors we will fragment over
            VectorView<algorithms::Anchor> anchor_view {anchors_to_fragment, anchor_indexes}; 
            std::vector<std::pair<int, std::vector<size_t>>> results = algorithms::find_best_chains(
                anchor_view,
                *distance_index,
                gbwt_graph,
                get_regular_aligner()->gap_open,
                get_regular_aligner()->gap_extension,
                this->max_fragments,
                for_each_transition,
                this->item_bonus,
                this->item_scale,
                this->fragment_gap_scale,
                indel_limit,
                false
            );
            if (show_work) {
                #pragma omp critical (cerr)
                cerr << log_name() << "Found " << results.size() << " fragments in zip code tree " << item_num
                    << " running " << anchors_to_fragment[anchor_indexes.front()] << " to " << anchors_to_fragment[anchor_indexes.back()] << std::endl;
            }
            for (size_t result = 0; result < results.size(); result++) {
                // For each result
                auto& scored_fragment = results[result];
                if (show_work) {
#ifdef debug
                    if(true)
#else
                    if (result < MANY_LIMIT)
#endif
                    {
                        if (!scored_fragment.second.empty()) {
                            #pragma omp critical (cerr)
                            {
                                cerr << log_name() << "\tFragment with score " << scored_fragment.first
                                    << " and length " << scored_fragment.second.size()
                                    << " running " << anchor_view[scored_fragment.second.front()]
                                    << " to " << anchor_view[scored_fragment.second.back()] << std::endl;
#ifdef debug
                                
                                for (auto& anchor_number : scored_fragment.second) {
                                    std::cerr << log_name() << "\t\t" << anchor_view[anchor_number] << std::endl;
                                }
#endif

                            }
                        }
                    } else if (result == MANY_LIMIT) {
                        #pragma omp critical (cerr)
                        std::cerr << log_name() << "\t<" << (results.size() - result) << " more fragments>" << std::endl;
                    }
                }

                // Count how many of each minimizer is in each fragment produced
                minimizer_kept_fragment_count.emplace_back(minimizers.size(), 0);

                // Translate fragments into seed numbers and not local anchor numbers.
                fragments.emplace_back();
                fragments.back().reserve(scored_fragment.second.size() * 2);
                for (auto& selected_number : scored_fragment.second) {
                    // For each anchor in the chain, get its number in the whole group of anchors.
                    size_t anchor_number = anchor_indexes.at(selected_number);
                    for (auto& seed_number : anchor_seed_sequences.at(anchor_number)) {
                        // And get all the seeds it actually uses in sequence and put them in the fragment.
                        fragments.back().push_back(seed_number);
                    }
                    for (auto& seed_number : anchor_represented_seeds.at(anchor_number)) {
                        // And get all the seeds it represents exploring and mark their minimizers explored.
                        // TODO: Can we get the gapless extension logic to count this for us for that codepath?
                        minimizer_kept_fragment_count.back()[seeds[seed_number].source]++;
                    }
                }
                // Remember the score
                fragment_scores.push_back(scored_fragment.first);
                // And make an anchor of it right now, for chaining later.
                // Make sure to do it by combining the gapless extension anchors if applicable.
                fragment_anchors.push_back(algorithms::Anchor(anchors_to_fragment.at(anchor_indexes.at(scored_fragment.second.front())), anchors_to_fragment.at(anchor_indexes.at(scored_fragment.second.back())), 0, 0, fragment_scores.back()));
                // Remember how we got it
                fragment_source_tree.push_back(item_num);
                //Remember the number of better or equal-scoring trees
                multiplicity_by_fragment.emplace_back((float)item_count);

                if (track_provenance) {
                    // Tell the funnel
                    funnel.introduce();
                    funnel.score(funnel.latest(), scored_fragment.first);
                    // We come from all the seeds directly
                    // TODO: Include all the middle seeds when gapless extending!
                    funnel.also_merge_group(2, fragments.back().begin(), fragments.back().end());
                    // And are related to the problem
                    funnel.also_relevant(1, item_num);
                }

                if (track_position && result < MANY_LIMIT) {
                    // Add position annotations for the good-looking fragments.
                    // Should be much faster than full correctness tracking from every seed.
                    crash_unless(this->path_graph);
                    for (auto& boundary : {anchor_view[scored_fragment.second.front()].graph_start(), anchor_view[scored_fragment.second.back()].graph_end()}) {
                        // For each end of the fragment
                        auto offsets = algorithms::nearest_offsets_in_paths(this->path_graph, boundary, 100);
                        for (auto& handle_and_positions : offsets) {
                            for (auto& position : handle_and_positions.second) {
                                // Tell the funnel all the effective positions, ignoring orientation
                                funnel.position(funnel.latest(), handle_and_positions.first, position.first);
                            }
                        }

                    }
                }
                if (track_provenance && show_work && result < MANY_LIMIT) {
                    for (auto& handle_and_range : funnel.get_positions(funnel.latest())) {
                        // Log each range on a path associated with the fragment.
                        #pragma omp critical (cerr)
                        std::cerr << log_name() << "\t\tAt linear reference "
                            << this->path_graph->get_path_name(handle_and_range.first)
                            << ":" << handle_and_range.second.first
                            << "-" << handle_and_range.second.second << std::endl;
                    }
                    if (track_correctness && funnel.is_correct(funnel.latest())) {
                        #pragma omp critical (cerr)
                        cerr << log_name() << "\t\tCORRECT!" << endl;
                    }
                }
            }

            
            if (track_provenance) {
                // Say we're done with this 
                funnel.processed_input();
            }
            
            return true;
            
        }, [&](size_t item_num) -> void {
            // There are too many sufficiently good problems to do
            if (track_provenance) {
                funnel.pass("zipcode-tree-coverage-threshold", item_num, tree_coverages[item_num]);
                funnel.fail("max-to-fragment", item_num);
            }
            
        }, [&](size_t item_num) -> void {
            // This item is not sufficiently good.
            if (track_provenance) {
                funnel.fail("zipcode-tree-coverage-threshold", item_num, tree_coverages[item_num]);
            }
        });

    //Get the actual multiplicity from the counts
    for (size_t i = 0 ; i < multiplicity_by_fragment.size() ; i++) {
        multiplicity_by_fragment[i] = multiplicity_by_fragment[i] >= kept_tree_count
                                    ?  multiplicity_by_fragment[i] - (float)kept_tree_count
                                    : 0.0;
    }
    // Now glom the fragments together into chains 
    if (track_provenance) {
        funnel.stage("chain");
    }
    
    if (track_provenance) {
        funnel.substage("chain");
    }
    
    // For each chain, we need:
    // The chain itself, pointing into seeds
    std::vector<std::vector<size_t>> chains;
    // The zip code tree it came from
    std::vector<size_t> chain_source_tree;
    // An estimated alignment score
    std::vector<int> chain_score_estimates;
    // A count, for each minimizer, of how many hits of it could have been in the chain, or were considered when making the chain.
    std::vector<std::vector<size_t>> minimizer_kept_chain_count;
    // The multiplicity for each chain. For now, just the multiplicity of the tree it came from
    std::vector<double> multiplicity_by_chain;
    
    // Get all the fragment numbers for each zip code tree we actually used, so we can chain each independently again.
    // TODO: Stop reswizzling so much.
    std::unordered_map<size_t, std::vector<size_t>> tree_to_fragments;
    vector<double> multiplicity_by_tree(zip_code_forest.trees.size(), 0);
    for (size_t i = 0; i < fragment_source_tree.size(); i++) {
        tree_to_fragments[fragment_source_tree[i]].push_back(i);
#ifdef debug
        if (multiplicity_by_tree[fragment_source_tree[i]] != 0) {
            assert(multiplicity_by_tree[fragment_source_tree[i]] == multiplicity_by_fragment[i]);
        }
#endif
        multiplicity_by_tree[fragment_source_tree[i]] = multiplicity_by_fragment[i];
    }
    
    // Get the score of the top-scoring fragment in each collection.
    std::unordered_map<size_t, double> best_fragment_score_in;
    // And overall
    double best_fragment_score = 0;
    for (auto& kv : tree_to_fragments) {
        for (auto& fragment_num : kv.second) {
            // Max in the score of each fragment 
            best_fragment_score_in[kv.first] = std::max(best_fragment_score_in[kv.first], fragment_scores.at(fragment_num));
            best_fragment_score = std::max(best_fragment_score, best_fragment_score_in[kv.first]);
        }
    }
    
    // Decide on how good fragments have to be to keep.
    double fragment_score_threshold = std::min(best_fragment_score * fragment_score_fraction, fragment_max_min_score);
    double fragment_score_threshold_overall = std::max(fragment_score_threshold, fragment_min_score);

    // Filter down to just the good ones, sorted by read start
    std::unordered_map<size_t, std::vector<size_t>> good_fragments_in;
    for (auto& kv : tree_to_fragments) {
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Keeping, of the " << kv.second.size() << " fragments in " << kv.first << ", those with score of at least "  << fragment_score_threshold_overall << endl;
            }
        }
        
        size_t fragments_kept = 0;

        // Keep the fragments that have good scores.
        for (auto& fragment_num : kv.second) {
            // For each fragment
            auto fragment_score = fragment_scores.at(fragment_num);
            if (fragment_score >= fragment_score_threshold) {
                // If its score is high enough vs. the best
                if (track_provenance) {
                    // Tell the funnel
                    funnel.pass("fragment-score-fraction||fragment-max-min-score", fragment_num, best_fragment_score != 0 ? (fragment_score / best_fragment_score) : 0.0);
                }

                if (fragment_score >= fragment_min_score) {
                    // And its score is high enough overall

                    if (track_provenance) {
                        // Tell the funnel
                        funnel.pass("fragment-min-score", fragment_num, fragment_score);
                    }

                    // Keep it.
                    good_fragments_in[kv.first].push_back(fragment_num);
                    fragments_kept++;
                } else {
                    // If its score is not high enough overall
                    if (track_provenance) {
                        // Tell the funnel
                        funnel.fail("fragment-min-score", fragment_num, fragment_score);
                    }
                }
            } else {
                // If its score is not high enough vs. the best
                if (track_provenance) {
                    // Tell the funnel
                    funnel.fail("fragment-score-fraction||fragment-max-min-score", fragment_num, best_fragment_score != 0 ? (fragment_score / best_fragment_score) : 0.0);
                } 
            }
        }
        
        if (fragments_kept > 1) {
            // Only access the vector if we put stuff in it, to avoid making
            // empty vectors. And only sort if there are multiple fragments. 
            
            // Now sort anchors by read start. Don't bother with shadowing.
            algorithms::sort_anchor_indexes(fragment_anchors, good_fragments_in[kv.first]);
        }

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "\tKept " << fragments_kept << "/" << kv.second.size() << " fragments." << endl;
            }
        }
    }

    // Draft trees to chain all the fragments of based on how good their fragment sets look. 
    std::vector<size_t> trees_with_good_fragments;
    std::vector<double> fragment_set_scores;
    trees_with_good_fragments.reserve(good_fragments_in.size());
    fragment_set_scores.reserve(good_fragments_in.size());
    for (auto& kv : good_fragments_in) {
        // Make a vector of the numbers of all the still-eligible trees
        trees_with_good_fragments.push_back(kv.first);
        // And score each set of fragments
        double fragment_set_score = 0;
        for (auto& anchor_index : kv.second) {
            fragment_set_score += fragment_anchors.at(anchor_index).score();
        }
        fragment_set_scores.push_back(fragment_set_score);
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "=====Creating chains=====" << endl;
        }
    }

    process_until_threshold_b<double>(fragment_set_scores,
        fragment_set_score_threshold, min_chaining_problems, max_chaining_problems, rng, 
        [&](size_t processed_num, size_t item_count) -> bool {
            // This tree's fragment set is good enough.
            // Called in descending score order
            
            // TODO: How should this connect to multiplicity_by_tree? Given that we're dropping whole trees again?

            // Look up which tree this is
            size_t tree_num = trees_with_good_fragments.at(processed_num);
            auto& tree_fragments = good_fragments_in[tree_num];

            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Tree " << tree_num << " has a good enough fragment set (score=" << fragment_set_scores[processed_num] << ")" << endl;
                    if (track_correctness) {
                        for (auto& fragment_num : tree_fragments) {
                            if (funnel.was_correct(fragment_num)) {
                                cerr << log_name() << "\tCORRECT!" << endl;
                                break;
                            }
                        }
                    }
                }
            }
            if (track_provenance) {
                for (auto& fragment_num : tree_fragments) {
                    funnel.pass("fragment-set-score-threshold", fragment_num, fragment_set_scores[processed_num]);
                    funnel.pass("max-chaining-problems", fragment_num);
                }
            }

            // Get a view of all the good fragments.
            // TODO: Should we just not make a global fragment anchor list?
            VectorView<algorithms::Anchor> fragment_view {fragment_anchors, tree_fragments};
            
            // We should not be making empty entries
            crash_unless(!fragment_view.empty());
            
            if (show_work) {
                #pragma omp critical (cerr)
                std::cerr << log_name() << "Chaining fragments from zip code tree " << tree_num << std::endl;
            } 

            // Compute lookback and indel limits based on read length.
            // Important since seed density goes down on longer reads.
            size_t lookback_limit = std::max(this->max_lookback_bases, (size_t)(this->max_lookback_bases_per_base * aln.sequence().size()));
            size_t indel_limit = std::max(this->max_indel_bases, (size_t)(this->max_indel_bases_per_base * aln.sequence().size()));

            // Chain up the fragments
            algorithms::transition_iterator for_each_transition = algorithms::zip_tree_transition_iterator(
                seeds,
                zip_code_forest.trees[tree_num],
                lookback_limit
            );
            std::vector<std::pair<int, std::vector<size_t>>> chain_results = algorithms::find_best_chains(
                fragment_view,
                *distance_index,
                gbwt_graph,
                get_regular_aligner()->gap_open,
                get_regular_aligner()->gap_extension,
                this->max_alignments,
                for_each_transition,
                this->item_bonus,
                this->item_scale,
                this->gap_scale,
                indel_limit,
                false
            );
            
            for (size_t result = 0; result < chain_results.size(); result++) {
                auto& chain_result = chain_results[result];
                // Each chain of fragments becomes a chain of seeds
                chains.emplace_back();
                auto& chain = chains.back();
                // With a source
                chain_source_tree.push_back(tree_num);
                // With a score
                chain_score_estimates.emplace_back(0);
                int& score = chain_score_estimates.back();
                // And counts of each minimizer kept
                minimizer_kept_chain_count.emplace_back();
                auto& minimizer_kept = minimizer_kept_chain_count.back();
                //Remember the multiplicity from the fragments. For now, it is just based on
                //the trees so it doesn't matter which fragment this comes from
                multiplicity_by_chain.emplace_back(multiplicity_by_tree[tree_num]);
                
                // We record the fragments that merge into each chain for reporting.
                std::vector<size_t> chain_fragment_nums_overall;
                chain_fragment_nums_overall.reserve(chain_result.second.size());
                
                for (const size_t& local_fragment: chain_result.second) {
                    // For each fragment in the chain
                               
                    // Get its fragment number out of all fragments
                    size_t fragment_num_overall = tree_fragments.at(local_fragment);
                    
                    // Save it
                    chain_fragment_nums_overall.push_back(fragment_num_overall);
                    
                    // Go get that fragment
                    auto& fragment = fragments.at(fragment_num_overall);
                        
                    // And append all the seed numbers to the chain
                    std::copy(fragment.begin(), fragment.end(), std::back_inserter(chain));
                    
                    // And count the score
                    score += fragment_scores.at(fragment_num_overall);
                    
                    // And count the kept minimizers
                    auto& fragment_minimizer_kept = minimizer_kept_fragment_count.at(fragment_num_overall);
                    if (minimizer_kept.size() < fragment_minimizer_kept.size()) {
                        minimizer_kept.resize(fragment_minimizer_kept.size());
                    }
                    for (size_t i = 0; i < fragment_minimizer_kept.size(); i++) {
                        minimizer_kept[i] += fragment_minimizer_kept[i];
                    }
                }
                if (track_provenance) {
                    // Say all those fragments became a chain
                    funnel.merge_group(chain_fragment_nums_overall.begin(), chain_fragment_nums_overall.end());
                    // With the total score
                    funnel.score(funnel.latest(), score);
                }
                if (show_work) {
                    if (result < MANY_LIMIT) {
                        #pragma omp critical (cerr)
                        {
                            std::cerr << log_name() << "Chain " << (chains.size() - 1) << " with score " << score << " is composed from local fragments:";
                            for (auto& f : chain_result.second) {
                                std::cerr << " " << f;
                            } 
                            std::cerr << std::endl;
                            std::cerr << log_name() << "Chain " << (chains.size() - 1) << " with score " << score << " is composed from global fragments:";
                            for (auto& f : chain_fragment_nums_overall) {
                                std::cerr << " " << f;
                            } 
                            std::cerr << std::endl;
                            std::cerr << log_name() << "Chain " << (chains.size() - 1) << " with score " << score << " contains seeds:";
                            for (auto& s : chains.back()) {
                                std::cerr << " " << s;
                            } 
                            std::cerr << std::endl;
                        }
                        if (track_provenance) {
                            for (auto& handle_and_range : funnel.get_positions(funnel.latest())) {
                                // Log each range on a path associated with the chain.
                                #pragma omp critical (cerr)
                                std::cerr << log_name() << "\tAt linear reference "
                                    << this->path_graph->get_path_name(handle_and_range.first)
                                    << ":" << handle_and_range.second.first
                                    << "-" << handle_and_range.second.second << std::endl;
                            }
                        }
                        if (track_correctness && funnel.is_correct(funnel.latest())) {
                            #pragma omp critical (cerr)
                            cerr << log_name() << "\tCORRECT!" << endl;
                        }
                    } else if (result == MANY_LIMIT) {
                        #pragma omp critical (cerr)
                        std::cerr << log_name() << "<" << (chain_results.size() - result) << " more chains>" << std::endl;
                    }
                } 
            }

            return true;

        }, [&](size_t processed_num) -> void {
            // There are too many sufficiently good fragment sets.
            size_t tree_num = trees_with_good_fragments.at(processed_num);
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Tree " << tree_num << " skipped because too many trees have good enough fragment sets (score=" << fragment_set_scores[processed_num] << ")" << endl;
                    if (track_correctness) {
                        for (auto& fragment_num : good_fragments_in[tree_num]) {
                            if (funnel.was_correct(fragment_num)) {
                                cerr << log_name() << "\tCORRECT!" << endl;
                                break;
                            }
                        }
                    }
                }
            }
            if (track_provenance) {
                for (auto& fragment_num : good_fragments_in[tree_num]) {
                    funnel.pass("fragment-set-score-threshold", fragment_num, fragment_set_scores[processed_num]);
                    funnel.fail("max-chaining-problems", fragment_num);
                }
            }
        }, [&](size_t processed_num) -> void {
            // This fragment set is not sufficiently good.
            size_t tree_num = trees_with_good_fragments.at(processed_num);
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Tree " << tree_num << " skipped because its fragment set is not good enough (score=" << fragment_set_scores[processed_num] << ")" << endl;
                    if (track_correctness) {
                        for (auto& fragment_num : good_fragments_in[tree_num]) {
                            if (funnel.was_correct(fragment_num)) {
                                cerr << log_name() << "\tCORRECT!" << endl;
                                break;
                            }
                        }
                    }
                }
            }
            if (track_provenance) {
                for (auto& fragment_num : good_fragments_in[tree_num]) {
                    funnel.fail("fragment-set-score-threshold", fragment_num, fragment_set_scores[processed_num]);
                }
            }
        });
    
    // Find the best chain
    size_t best_chain = std::numeric_limits<size_t>::max();
    int best_chain_score = 0;
    for (size_t i = 0; i < chains.size(); i++) {
        if (best_chain == std::numeric_limits<size_t>::max() || chain_score_estimates.at(i) > best_chain_score) {
            // Friendship ended with old chain
            best_chain = i;
            best_chain_score = chain_score_estimates[i];
        }
    }
    bool best_chain_correct = false;
    if (track_correctness && best_chain != std::numeric_limits<size_t>::max()) {
        // We want to explicitly check if the best chain was correct, for looking at stats about it later.
        if (funnel.is_correct(best_chain)) {
            best_chain_correct = true;
        }
    }

    if (show_work && best_chain != std::numeric_limits<size_t>::max()) {
        // Dump the best chain

        auto& tree_num = chain_source_tree.at(best_chain);
        
        // Find all the seeds in its zip tree
        vector<size_t> involved_seeds;
        for (ZipCodeTree::oriented_seed_t found : zip_code_forest.trees.at(tree_num)) {
            involved_seeds.push_back(found.seed);   
        }

        // Start making a list of things to show. 
        std::vector<std::pair<std::string, std::vector<std::vector<size_t>>>> seed_sets;
        seed_sets.emplace_back("", std::vector<std::vector<size_t>>{std::move(involved_seeds)});
        seed_sets.emplace_back("chain", std::vector<std::vector<size_t>>{chains.at(best_chain)});

        // Find all the fragments we passed for this tree
        std::vector<std::vector<size_t>> relevant_fragments;
        auto& tree_fragments = good_fragments_in[tree_num];
        for (auto& fragment_num : tree_fragments) {
            // Get all the seeds in each fragment
            const std::vector<size_t>& fragment = fragments.at(fragment_num);
            relevant_fragments.push_back(fragment);
        }
        seed_sets.emplace_back("frag", std::move(relevant_fragments));

        // Sort everything in read order
        for (auto& seed_set : seed_sets) {
            for (auto& run : seed_set.second) {
                std::sort(run.begin(), run.end(), [&](const size_t& seed_index_a, const size_t& seed_index_b) {
                    auto& seed_a = seeds.at(seed_index_a);
                    auto& seed_b = seeds.at(seed_index_b);
                    
                    return minimizers[seed_a.source].forward_offset() < minimizers[seed_b.source].forward_offset();
    
                });
            }
        }


        dump_debug_dotplot("best-chain", minimizers, seeds, seed_sets, this->path_graph);

    }
    
    // Find its coverage
    double best_chain_coverage = 0;
    if (best_chain != std::numeric_limits<size_t>::max()) {
        best_chain_coverage = get_read_coverage(aln, std::vector<std::vector<size_t>> {chains.at(best_chain)}, seeds, minimizers);
    }
    
    // Find out how gappy it is. We can get the longest and the average distance maybe.
    size_t best_chain_longest_jump = 0;
    size_t best_chain_total_jump = 0;
    double best_chain_average_jump = 0;
    if (best_chain != std::numeric_limits<size_t>::max()) {
        for (size_t i = 1; i < chains.at(best_chain).size(); i++) {
            // Find the pair of anchors we go between
            auto& left_anchor = seed_anchors.at(chains.at(best_chain).at(i - 1));
            auto& right_anchor = seed_anchors.at(chains.at(best_chain).at(i));
            // And get the distance between them in the read
            size_t jump = right_anchor.read_start() - left_anchor.read_end();
            // Max and add it in
            best_chain_longest_jump = std::max(best_chain_longest_jump, jump);
            best_chain_total_jump += jump;
        }
        best_chain_average_jump = chains.at(best_chain).size() > 1 ? best_chain_total_jump / (chains.at(best_chain).size() - 1) : 0.0;
    }

    // Also count anchors in the chain
    size_t best_chain_anchors = 0;
    if (best_chain != std::numeric_limits<size_t>::max()) {
        best_chain_anchors = chains.at(best_chain).size();
    }

    // And total length of anchors in the chain
    size_t best_chain_anchor_length = 0;
    if (best_chain != std::numeric_limits<size_t>::max()) {
        for (auto& item : chains.at(best_chain)) {
            best_chain_anchor_length += seed_anchors.at(item).length();
        }
    }
    
    if (track_provenance) {
        funnel.stage("align");
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "=====Creating alignments=====" << endl;
        }
    }

#ifdef print_minimizer_table
    //How many of each minimizer ends up in a chain that actually gets turned into an alignment?
    vector<size_t> minimizer_kept_count(minimizers.size(), 0);
#endif
    
    // Now start the alignment step. Everything has to become an alignment.

    // We will fill this with all computed alignments in estimated score order.
    vector<Alignment> alignments;
    alignments.reserve(chain_score_estimates.size());
    // This maps from alignment index back to chain index, for
    // tracing back to minimizers for MAPQ. Can hold
    // numeric_limits<size_t>::max() for an unaligned alignment.
    vector<size_t> alignments_to_source;
    alignments_to_source.reserve(chain_score_estimates.size());
    //For finding the multiplicity of each alignment, first get the count
    // of equal scoring chains
    vector<size_t> chain_count_by_alignment (alignments.size(), 0);
    //The multiplicity for each alignment, projected from previous stages
    vector<double> multiplicity_by_alignment;

    // Create a new alignment object to get rid of old annotations.
    {
      Alignment temp;
      temp.set_sequence(aln.sequence());
      temp.set_name(aln.name());
      temp.set_quality(aln.quality());
      aln = std::move(temp);
    }

    // Annotate the read with metadata
    if (!sample_name.empty()) {
        aln.set_sample_name(sample_name);
    }
    if (!read_group.empty()) {
        aln.set_read_group(read_group);
    }
    
    // We need to be able to discard a chain because its score isn't good enough.
    // We have more components to the score filter than process_until_threshold_b supports.
    auto discard_chain_by_score = [&](size_t processed_num) -> void {
        // This chain is not good enough.
        if (track_provenance) {
            funnel.fail("min-chain-score-per-base||max-min-chain-score", processed_num, chain_score_estimates[processed_num]);
        }
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "chain " << processed_num << " failed because its score was not good enough (score=" << chain_score_estimates[processed_num] << ")" << endl;
                if (track_correctness && funnel.was_correct(processed_num)) {
                    cerr << log_name() << "\tCORRECT!" << endl;
                }
            }
        }
    };
    
    // Compute lower limit on chain score to actually investigate
    int chain_min_score = std::min((int) (min_chain_score_per_base * aln.sequence().size()), max_min_chain_score);

    // Track if minimizers were explored by alignments
    SmallBitset minimizer_explored(minimizers.size());

    // Track how many tree chains were used
    std::unordered_map<size_t, size_t> chains_per_tree;

    // Track what read offset, graph node pairs were used in previously generated alignments, so we can fish out alignments to different placements.
    std::unordered_set<std::pair<size_t, pos_t>> used_matchings;

    // Track statistics about how many bases were aligned by diffrent methods, and how much time was used.
    aligner_stats_t stats; 
    
    // Go through the chains in estimated-score order.
    process_until_threshold_b<int>(chain_score_estimates,
        chain_score_threshold, min_chains, max_alignments, rng, 
        [&](size_t processed_num, size_t item_count) -> bool {
            // This chain is good enough.
            // Called in descending score order.
        
            if (chain_score_estimates[processed_num] < chain_min_score) {
                // Actually discard by score
                discard_chain_by_score(processed_num);
                return false;
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Chain " << processed_num << " is good enough (score=" << chain_score_estimates[processed_num] << "/" << chain_min_score << ")" << endl;
                    if (track_correctness && funnel.was_correct(processed_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                    }
                }
            }
            if (track_provenance) {
                funnel.pass("min-chain-score-per-base||max-min-chain-score", processed_num, chain_score_estimates[processed_num]);
                funnel.pass("max-alignments", processed_num);
            }

            for (auto& seed_num : chains[processed_num]) {
                auto matching = std::make_pair(minimizers[seeds.at(seed_num).source].forward_offset(), seeds.at(seed_num).pos);
                if (used_matchings.count(matching)) {
                    if (track_provenance) {
                        funnel.fail("no-chain-overlap", processed_num);
                    }
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "Chain " << processed_num << " overlaps a previous alignment at read position " << matching.first << " and graph position " << matching.second << endl;
                        }
                    }
                    return false;
                }
            }
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Chain " << processed_num << " overlaps none of the " << used_matchings.size() << " read-node matchings used in previous alignments" << endl;
                }
            }
            if (track_provenance) {
                funnel.pass("no-chain-overlap", processed_num);
            }

            // Make sure we aren't doing too many chains from this one tree.
            auto& tree_count = chains_per_tree[chain_source_tree[processed_num]];
            if (tree_count >= max_chains_per_tree) {
                if (track_provenance) {
                    funnel.fail("max-chains-per-tree", processed_num, tree_count);
                }
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Chain " << processed_num << " is chain " << tree_count << " in its tree " << chain_source_tree[processed_num] << " and is rejected (score=" << chain_score_estimates[processed_num] << ")" << endl;
                    }
                }
                tree_count++;
                return false;
            } else {
                if (track_provenance) {
                    funnel.pass("max-chains-per-tree", processed_num, tree_count);
                }
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Chain " << processed_num << " is chain " << tree_count << " in its tree " << chain_source_tree[processed_num] << " and is kept" << endl;
                    }
                }
                tree_count++;
            }

            if (track_provenance) {
                funnel.processing_input(processed_num);
            }

            // Collect the top alignments. Make sure we have at least one always, starting with unaligned.
            vector<Alignment> best_alignments(1, aln);

            // Align from the chained-up seeds
            if (do_dp) {
                // We need to do base-level alignment.
            
                if (track_provenance) {
                    funnel.substage("align");
                }
                
                // We currently just have the one best score and chain per zip code tree
                vector<size_t>& chain = chains[processed_num];
                
                try {
                    // Do the DP between the items in the chain

                    // Collect stats into here
                    aligner_stats_t alignment_stats;
                    best_alignments[0] = find_chain_alignment(aln, seed_anchors, chain, &alignment_stats);
                    alignment_stats.add_annotations(best_alignments[0], "alignment");

                    // Remember the stats' usages
                    stats += alignment_stats;
                } catch (ChainAlignmentFailedError& e) {
                    // We can't actually make an alignment from this chain
                    #pragma omp critical (cerr)
                    cerr << log_name() << "Error creating alignment from chain for " << aln.name() << ": " << e.what() << endl;
                    // Leave the read unmapped.
                }

                if (track_provenance) {
                    funnel.substage_stop();
                }
                    
                // TODO: Come up with a good secondary somehow.
            } else {
                // We would do base-level alignment but it is disabled.
                // Leave best_alignment unaligned
            }
           
            // Have a function to process the best alignments we obtained
            auto observe_alignment = [&](Alignment& aln) {
                alignments.emplace_back(std::move(aln));
                alignments_to_source.push_back(processed_num);
                multiplicity_by_alignment.emplace_back(multiplicity_by_chain[processed_num]);
                chain_count_by_alignment.emplace_back(item_count);
                
                size_t read_pos = 0;
                for (auto& mapping : alignments.back().path().mapping()) {
                    // Mark all the read-node matches it visits used.
                    used_matchings.emplace(read_pos, make_pos_t(mapping.position()));
                    read_pos += mapping_to_length(mapping);
                }

                if (track_provenance) {
                    funnel.project(processed_num);
                    funnel.score(alignments.size() - 1, alignments.back().score());
                }
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Produced alignment from chain " << processed_num
                            << " with score " << alignments.back().score() << ": " << log_alignment(alignments.back()) << endl;
                    }
                }
            };
            
            if (!best_alignments.empty() && best_alignments[0].score() <= 0) {
                if (show_work) {
                    // Alignment won't be observed but log it anyway.
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Produced terrible best alignment from chain " << processed_num << ": " << log_alignment(best_alignments[0]) << endl;
                    }
                }
            }
            for(auto aln_it = best_alignments.begin() ; aln_it != best_alignments.end() && aln_it->score() != 0 && aln_it->score() >= best_alignments[0].score() * 0.8; ++aln_it) {
                //For each additional alignment with score at least 0.8 of the best score
                observe_alignment(*aln_it);
            }
           
            if (track_provenance) {
                // We're done with this input item
                funnel.processed_input();
            }

            if (track_provenance) {
                funnel.substage("minimizers_kept");
            }

            for (size_t i = 0 ; i < minimizer_kept_chain_count[processed_num].size() ; i++) {
#ifdef print_minimizer_table
                minimizer_kept_count[i] += minimizer_kept_chain_count[processed_num][i];
#endif
                if (use_explored_cap && minimizer_kept_chain_count[processed_num][i] > 0) {
                    // This minimizer is in a zip code tree that gave rise
                    // to at least one alignment, so it is explored.
                    minimizer_explored.insert(i);
                }
            }

            if (track_provenance) {
                funnel.substage_stop();
            }
            
            return true;
        }, [&](size_t processed_num) -> void {
            // There are too many sufficiently good chains
            if (track_provenance) {
                funnel.pass("min-chain-score-per-base||max-min-chain-score", processed_num, chain_score_estimates[processed_num]);
                funnel.fail("max-alignments", processed_num);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "chain " << processed_num << " failed because there were too many good chains (score=" << chain_score_estimates[processed_num] << ")" << endl;
                    if (track_correctness && funnel.was_correct(processed_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                    }
                }
            }
        }, discard_chain_by_score);
    
    // We want to be able to feed in an unaligned alignment on the normal
    // codepath, but we don't want it to really participate in the funnel
    // filters anymore. So we set this flag if the funnle is really empty of
    // items so we stop talking about filters.
    bool funnle_depleted = false;

    if (alignments.size() == 0) {
        // Produce an unaligned Alignment
        alignments.emplace_back(aln);
        alignments_to_source.push_back(numeric_limits<size_t>::max());
        multiplicity_by_alignment.emplace_back(0);
        // Stop telling the funnel about filters and items.
        funnle_depleted = true;
    } else {
        //chain_count_by_alignment is currently the number of better or equal chains that were used
        // We really want the number of chains not including the ones that represent the same mapping
        // TODO: This isn't very efficient
        for (size_t i = 0 ; i < chain_count_by_alignment.size() ; ++i) {
            size_t chain_i = alignments_to_source[i];
            for (size_t j = 0 ; j < chain_count_by_alignment.size() ; ++j) {
                size_t chain_j = alignments_to_source[j];
                if (i != j &&
                    chain_score_estimates[chain_i] >= chain_score_estimates[chain_j] &&
                    chain_ranges_are_equivalent(seeds[chains[chain_i].front()],
                                         seeds[chains[chain_i].back()],
                                         seeds[chains[chain_j].front()],
                                         seeds[chains[chain_j].back()])) {
                    --chain_count_by_alignment[i];
                }
            }
        }
        for (size_t i = 0 ; i < multiplicity_by_alignment.size() ; ++i) {
            multiplicity_by_alignment[i] += (chain_count_by_alignment[i] >= alignments.size()
                                          ? ((double)chain_count_by_alignment[i] - (double) alignments.size())
                                          : 0.0);
        }
    }
    
    if (track_provenance) {
        // Now say we are finding the winner(s)
        funnel.stage("winner");
    }
    
    // Fill this in with the alignments we will output as mappings
    vector<Alignment> mappings;
    mappings.reserve(min(alignments.size(), max_multimaps));
    
    // Grab all the scores in order for MAPQ computation.
    vector<double> scores;
    scores.reserve(alignments.size());
    
    // Go through the alignments in descending score order, with ties at the top end shuffled.
    process_until_threshold_a(alignments.size(), (std::function<double(size_t)>) [&](size_t i) -> double {
        return alignments.at(i).score();
    }, 0, 1, max_multimaps, rng, [&](size_t alignment_num, size_t item_count) {
        // This alignment makes it
        // Called in score order
        
        // Remember the score at its rank
        scores.emplace_back(alignments[alignment_num].score());
        
        // Remember the output alignment
        mappings.emplace_back(std::move(alignments[alignment_num]));
        
        if (track_provenance && !funnle_depleted) {
            // Tell the funnel
            funnel.pass("max-multimaps", alignment_num);
            funnel.project(alignment_num);
            funnel.score(funnel.latest(), scores.back());
        }
        
        return true;
    }, [&](size_t alignment_num) {
        // We already have enough alignments, although this one has a good score
        
        // Remember the score at its rank anyway
        scores.emplace_back(alignments[alignment_num].score());
        
        if (track_provenance && !funnle_depleted) {
            funnel.fail("max-multimaps", alignment_num);
        }
    }, [&](size_t alignment_num) {
        // This alignment does not have a sufficiently good score
        // Score threshold is 0; this should never happen
        crash_unless(false);
    });
    
    if (track_provenance) {
        funnel.substage("mapq");
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Picked best alignment " << log_alignment(mappings[0]) << endl;
            cerr << log_name() << "For scores";
            for (auto& score : scores) cerr << " " << score << ":" << endl;
        }
    }

    vector<double> scaled_scores;
    scaled_scores.reserve(scores.size());
    for (auto& score : scores) {
        scaled_scores.push_back(score * mapq_score_scale);
    }

    crash_unless(!mappings.empty());
    // Compute MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
    // Use exact mapping quality 
    double mapq = (mappings.front().path().mapping_size() == 0) ? 0 : 
        get_regular_aligner()->compute_max_mapping_quality(scaled_scores, false, &multiplicity_by_alignment) ;

#ifdef debug_write_minimizers
#pragma omp critical
    {
        std::ofstream out;
        out.open("minimizers.tsv", std::ios::app);
        out << aln.name() << "\t" << mapq << "\t" << aln.sequence().size();
        for (size_t i = 0 ; i < minimizers.size() ; i++) {
            out << "\t";
            out << minimizer_kept[i]
                << "," << passed_downsampling[minimizer_score_order[i]]
                << "," << minimizers[i].hits 
                << "," << minimizers[i].score
                << "," << minimizers[i].forward_offset()
                << "," << minimizers[i].length;
        }
        out << endl;
        out.close(); 
    }
#endif
    
#ifdef print_minimizer_table
    double uncapped_mapq = mapq;
#endif

    set_annotation(mappings.front(), "mapq_uncapped", mapq);

    if (use_explored_cap) {

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "uncapped MAPQ is " << mapq << endl;
            }
        }
    
        // TODO: give SmallBitset iterators so we can use it instead of an index vector.
        vector<size_t> explored_minimizers;
        for (size_t i = 0; i < minimizers.size(); i++) {
            if (minimizer_explored.contains(i)) {
                explored_minimizers.push_back(i);
            }
        }
        // Compute caps on MAPQ. TODO: avoid needing to pass as much stuff along.
        double escape_bonus = mapq < std::numeric_limits<int32_t>::max() ? 1.0 : 2.0;
        double mapq_explored_cap = escape_bonus * faster_cap(minimizers, explored_minimizers, aln.sequence(), aln.quality());

        set_annotation(mappings.front(), "mapq_explored_cap", mapq_explored_cap);

        // Apply the caps and transformations
        mapq = round(min(mapq_explored_cap, mapq));

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Explored cap is " << mapq_explored_cap << endl;
            }
        }
    }


    // Make sure to clamp 0-60.
    mapq = max(mapq, 0.0);
    mapq = min(mapq, 60.0);
    // And save the MAPQ
    mappings.front().set_mapping_quality(mapq);

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "MAPQ is " << mapq << endl;
        }
    }

    // Remember the scores
    set_compressed_annotation(mappings.front(),"secondary_scores", scores);

    if (track_provenance) {
        funnel.substage_stop();
    }
    
    for (size_t i = 0; i < mappings.size(); i++) {
        // For each output alignment in score order
        auto& out = mappings[i];
        
        // Assign primary and secondary status
        out.set_is_secondary(i > 0);
    }

    if (this->set_refpos) {
        if (track_provenance) {
            // Time how long setting reference positions takes
            funnel.substage("refpos");
        }

        crash_unless(path_graph != nullptr);
        for (auto& m : mappings) {
            // Annotate the reads with the positions of the nodes they are actually on (fast)
            vg::algorithms::annotate_with_node_path_positions(*path_graph, m, -1);
        }
    }
    
    // Stop this alignment
    funnel.stop();

    // Annotate with whatever's in the funnel
    funnel.annotate_mapped_alignment(mappings[0], track_correctness);
    
    if (track_provenance) {
        if (track_correctness) {
            annotate_with_minimizer_statistics(mappings[0], minimizers, seeds, seeds.size(), fragments.size(), funnel);
        }
    }
    
    // Special fragment and chain statistics
    set_compressed_annotation(mappings[0], "fragment_scores", fragment_scores);
    if (track_correctness) {
        set_annotation(mappings[0], "best_chain.correct", best_chain_correct);
    }
    set_annotation(mappings[0], "best_chain.coverage", best_chain_coverage);
    set_annotation(mappings[0], "best_chain.longest_jump", (double) best_chain_longest_jump);
    set_annotation(mappings[0], "best_chain.average_jump", best_chain_average_jump);
    set_annotation(mappings[0], "best_chain.anchors", (double) best_chain_anchors);
    set_annotation(mappings[0], "best_chain.anchor_length", (double) best_chain_anchor_length);

    stats.add_annotations(mappings[0], "read");
    
#ifdef print_minimizer_table
    cerr << aln.sequence() << "\t";
    for (char c : aln.quality()) {
        cerr << (char)(c+33);
    }
    cerr << "\t" << zip_code_forest.trees.size();
    for (size_t i = 0 ; i < minimizers.size() ; i++) {
        auto& minimizer = minimizers[i];
        cerr << "\t"
             << minimizer.value.key.decode(minimizer.length) << "\t"
             << minimizer.forward_offset() << "\t"
             << minimizer.agglomeration_start << "\t"
             << minimizer.agglomeration_length << "\t"
             << minimizer.hits << "\t"
             << minimizer_kept_count[i];
         if (minimizer_kept_count[i]>0) {
             assert(minimizer.hits<=hard_hit_cap) ;
         }
    }
    cerr << "\t" << uncapped_mapq << "\t" << mapq_explored_cap << "\t"  << mappings.front().mapping_quality() << "\t";
    cerr << "\t";
    for (auto& score : scores) {
        cerr << score << ",";
    }
    if (track_correctness) {
        cerr << "\t" << funnel.last_correct_stage() << endl;
    } else {
        cerr << "\t" << "?" << endl;
    }
#endif

    if (track_provenance) {
        if (show_work && aln.sequence().size() < LONG_LIMIT) {
            // Dump the funnel info graph to standard error
            #pragma omp critical (cerr)
            {
                funnel.to_dot(cerr);
            }
        }
        
        // Otherwise/also, if we are dumping explanations, dump it to a file
        DotDumpExplainer<Funnel> explainer(true, funnel);
    }

    return mappings;
}

double MinimizerMapper::get_read_coverage(
    const Alignment& aln,
    const VectorView<std::vector<size_t>>& seed_sets,
    const std::vector<Seed>& seeds,
    const VectorView<Minimizer>& minimizers) const {
    
    std::vector<bool> covered(aln.sequence().size(), false);
    
    for (auto& list : seed_sets) {
        // We will fill in the range it occupies in the read
        std::pair<size_t, size_t> read_range {std::numeric_limits<size_t>::max(), 0};
        
        for (auto& seed_index : list) {
            // Which means we look at the minimizer for each seed
            auto& seed = seeds.at(seed_index);
            crash_unless(seed.source < minimizers.size());
            auto& minimizer = minimizers[seed.source];
            
            if (minimizer.forward_offset() < read_range.first) {
                // Min all their starts to get the start
                read_range.first = minimizer.forward_offset();
            }
            
            if (minimizer.forward_offset() + minimizer.length > read_range.second) {
                // Max all their past-ends to get the past-end
                read_range.second = minimizer.forward_offset() + minimizer.length;
            }
        }
        
        // Then mark its coverage
        set_coverage_flags(covered, read_range.first, read_range.second);
    }
    
    // And return the fraction covered.
    return get_fraction_covered(covered);
}

Alignment MinimizerMapper::find_chain_alignment(
    const Alignment& aln,
    const VectorView<algorithms::Anchor>& to_chain,
    const std::vector<size_t>& chain,
    aligner_stats_t* stats
) const {
    
    if (chain.empty()) {
        throw ChainAlignmentFailedError("Cannot find an alignment for an empty chain!");
    }
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Align chain of";
            if (chain.size() < MANY_LIMIT) {
                cerr << ": ";
                for (auto item_number : chain) {
                    cerr << " " << item_number;
                }
            } else {
                cerr << " " << chain.size() << " items";
            }
            cerr << " in " << to_chain.size() << " items" << endl;
        }
    }
    
    // We need an Aligner for scoring.
    const Aligner& aligner = *get_regular_aligner();

    // We need an ErrorModel to limit what our WFAExtender is allowed to do.
    // The ErrorModel is in terms of mismatches, gaps, and gap extensions, but if you fill them all in then a problem is allowed to have that many of *all* of those.
    // So we set a limit just in mismatches, and if fewer mismatches than that are used some gaps will be allowed.
    WFAExtender::ErrorModel wfa_error_model {
        {wfa_max_mismatches_per_base, wfa_max_mismatches, wfa_max_max_mismatches},
        {0, 0, 0},
        {0, 0, 0},
        {wfa_distance_per_base, wfa_distance, wfa_max_distance}
    };
    
    // We need a WFAExtender to do tail and intervening alignments.
    // Note that the extender expects anchoring matches!!!
    WFAExtender extender(gbwt_graph, aligner, wfa_error_model); 
    
    // Keep a couple cursors in the chain: extension before and after the linking up we need to do.
    auto here_it = chain.begin();
    auto next_it = here_it;
    ++next_it;
    
    // Track the anchor we're at.
    // Note that, although it has a score, that's an anchor score; it isn't the
    // right score for the perfect-match alignment it represents.
    const algorithms::Anchor* here = &to_chain[*here_it];
    
#ifdef debug_chain_alignment
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "First item " << *here_it
                << " with overall index " << to_chain.backing_index(*here_it)
                << " aligns " << (*here).read_start() << "-" << (*here).read_end()
                << " with " << (*here).graph_start() << "-" << (*here).graph_end()
                << endl;
        }
    }
#endif

    // We time each alignment operation using this scratch.
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point stop_time;

    
    // We compose into a Path, since sometimes we may have to drop back to
    // aligners that aren't the WFAAligner and don't make WFAAlignments.
    Path composed_path;
    // We also track the total score of all the pieces.
    int composed_score = 0;

    // Do the left tail, if any.
    size_t left_tail_length = (*here).read_start();
    if (left_tail_length > 0) {
        // We need to do a left tail.
        // Anchor position will not be covered.
        
        string left_tail = aln.sequence().substr(0, left_tail_length);
        WFAAlignment left_alignment;
        pos_t right_anchor = (*here).graph_start();
        if (left_tail.size() <= max_tail_length) {
            // Tail is short so keep to the GBWT.
            // We align the left tail with prefix(), which creates a prefix of the alignment.
            if (stats) {
                start_time = std::chrono::high_resolution_clock::now();
            }
            left_alignment = extender.prefix(left_tail, right_anchor);
            if (stats) {
                stop_time = std::chrono::high_resolution_clock::now();
                stats->bases.wfa_tail += left_tail_length;
                stats->time.wfa_tail += std::chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count();
                stats->invocations.wfa_tail += 1;
            }
            if (left_alignment && left_alignment.seq_offset != 0) {
                // We didn't get all the way to the left end of the read without
                // running out of score.
                // Prepend a softclip.
                // TODO: Can we let the aligner know it can softclip for free?
                WFAAlignment prepend = WFAAlignment::make_unlocalized_insertion(0, left_alignment.seq_offset, 0);
                prepend.join(left_alignment);
                left_alignment = std::move(prepend);
            }
            if (left_alignment.length != (*here).read_start()) {
                // We didn't get the alignment we expected.
                stringstream ss;
                ss << "Aligning left tail " << left_tail << " from " << (*here).graph_start() << " produced wrong-length alignment ";
                left_alignment.print(ss);
                throw ChainAlignmentFailedError(ss.str());
            }
        }
        if (left_alignment) {
            // We got an alignment, so make it a path
            left_alignment.check_lengths(gbwt_graph);
            
#ifdef debug_chain_alignment
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Start with left tail of " << left_alignment.length << " with score of " << left_alignment.score << endl;
                }
            }
#endif
            
            composed_path = left_alignment.to_path(this->gbwt_graph, aln.sequence());
            composed_score = left_alignment.score;
        } else {
            // We need to fall back on alignment against the graph
            
            if (left_tail_length > MAX_DP_LENGTH) {
                // Left tail is too long to align.
                
#ifdef debug_chain_alignment
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Refusing to align " << left_tail_length << " bp left tail against " << right_anchor << " in " << aln.name() << " to avoid overflow" << endl;
                }
#endif
                
                // Make a softclip for it.
                left_alignment = WFAAlignment::make_unlocalized_insertion(0, left_tail.size(), 0);
                composed_path = left_alignment.to_path(this->gbwt_graph, aln.sequence());
                composed_score = left_alignment.score;
            } else {
            
#ifdef debug_chain_alignment
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Start with long left tail fallback alignment" << endl;
                    }
                }
#endif
                
                Alignment tail_aln;
                tail_aln.set_sequence(left_tail);
                if (!aln.quality().empty()) {
                    tail_aln.set_quality(aln.quality().substr(0, left_tail_length));
                }
                
                // Work out how far the tail can see
                size_t max_gap_length = std::min(this->max_tail_gap, longest_detectable_gap_in_range(aln, aln.sequence().begin(), aln.sequence().begin() + left_tail_length, this->get_regular_aligner()));
                size_t graph_horizon = left_tail_length + max_gap_length;

#ifdef warn_on_fallback
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << left_tail_length << " bp left tail against " << right_anchor << " allowing " << max_gap_length << " bp gap in " << aln.name() << endl;
                }
#endif

                // Align the left tail, anchoring the right end.
                if (stats) {
                    start_time = std::chrono::high_resolution_clock::now();
                }
                auto nodes_and_bases = align_sequence_between(empty_pos_t(), right_anchor, graph_horizon, max_gap_length, &this->gbwt_graph, this->get_regular_aligner(), tail_aln, &aln.name(), this->max_dp_cells, this->choose_band_padding);
                if (stats) {
                    stop_time = std::chrono::high_resolution_clock::now();
                    if (nodes_and_bases.first > 0) {
                        // Actually did the alignment
                        stats->bases.dozeu_tail += left_tail_length;
                        stats->time.dozeu_tail += std::chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count();
                        stats->invocations.dozeu_tail += 1;
                    }
                }

                
                if (show_work && max_tail_length > 0) {
                    #pragma omp critical (cerr)
                    {
                        cerr << "warning[MinimizerMapper::find_chain_alignment]: Fallback score: " << tail_aln.score() << endl;
                    }
                }

                // Since it's the left tail we can just clobber the path
                composed_path = tail_aln.path();
                composed_score = tail_aln.score();
            }
        }
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Aligned left tail length " << left_tail_length << std::endl;
            }
        }

    }
        
    size_t longest_attempted_connection = 0;
    while(next_it != chain.end()) {
        // Do each region between successive gapless extensions
        
        // We have to find the next item we can actually connect to
        const algorithms::Anchor* next;
        // And the actual connecting alignment to it
        WFAAlignment link_alignment;
        // Where did it come from?
        std::string link_alignment_source;
        
        while (next_it != chain.end()) {
            next = &to_chain[*next_it];
            // Try and find a next thing to connect to
            
            if (algorithms::get_read_distance(*here, *next) == std::numeric_limits<size_t>::max()) {
                // There's overlap between these items. Keep here and skip next.
#ifdef debug_chain_alignment
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Don't try and connect " << *here_it << " to " << *next_it << " because they overlap" << endl;
                    }
                }
#endif
            
                ++next_it;
            } else {
                // No overlap, so try it.
                break;
            }
        }
        
        if (next_it == chain.end()) {
            // We couldn't find anything to connect to
            break;
        }
            
#ifdef debug_chain_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Add current item " << *here_it << " of length " << (*here).length() << endl;
            }
        }
#endif
        
        // Make an alignment for the bases used in this item, and
        // concatenate it in.
        WFAAlignment here_alignment = this->to_wfa_alignment(*here, aln, &aligner);

#ifdef debug_chain_alignment
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "\tScore " << here_alignment.score << endl;
        }
    }
#endif

        append_path(composed_path, here_alignment.to_path(this->gbwt_graph, aln.sequence()));
        composed_score += here_alignment.score;
        
#ifdef debug_chain_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Next connectable item " << *next_it
                    << " with overall index " << to_chain.backing_index(*next_it)
                    << " aligns " << (*next).read_start() << "-" << (*next).read_end()
                    << " with " << (*next).graph_start() << "-" << (*next).graph_end()
                    << endl;
            }
        }
#endif

        // Pull out the intervening string to the next, if any.
        size_t link_start = (*here).read_end();
        size_t link_length = (*next).read_start() - link_start;
        string linking_bases = aln.sequence().substr(link_start, link_length);
        size_t graph_length = algorithms::get_graph_distance(*here, *next, *distance_index, gbwt_graph);
        
#ifdef debug_chain_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Need to align graph from " << (*here).graph_end() << " to " << (*next).graph_start()
                    << " separated by " << graph_length << " bp and sequence \"" << linking_bases << "\"" << endl;
            }
        }
#endif
        
        if (link_length == 0 && graph_length == 0) {
            // These items abut in the read and the graph, so we assume we can just connect them.
            // WFAExtender::connect() can't handle an empty read sequence, and
            // our fallback method to align just against the graph can't handle
            // an empty graph region.
            // TODO: We can be leaving the GBWT's space here!
            
#ifdef debug_chain_alignment
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Treat as empty link" << endl;
                }
            }
#endif
            
            link_alignment = WFAAlignment::make_empty();
            link_alignment_source = "empty";
        } else if (link_length > 0 && link_length <= max_chain_connection) {
            // If it's not empty and is a reasonable size, align it.
            // Make sure to walk back the left anchor so it is outside of the region to be aligned.
            pos_t left_anchor = (*here).graph_end();
            get_offset(left_anchor)--;
            
            if (stats) {
                start_time = std::chrono::high_resolution_clock::now();
            }
            link_alignment = extender.connect(linking_bases, left_anchor, (*next).graph_start());
            if (stats) {
                stop_time = std::chrono::high_resolution_clock::now();
                stats->bases.wfa_middle += link_length;
                stats->time.wfa_middle += std::chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count();
                stats->invocations.wfa_middle += 1;
                if (!link_alignment) {
                    // Note that we had to fall back from WFA
                    stats->fallbacks.wfa_middle += 1;
                } else {
                    stats->fallbacks.wfa_middle += 0;
                }
            }
            link_alignment_source = "WFAExtender";

            longest_attempted_connection = std::max(longest_attempted_connection, linking_bases.size());
            
            if (!link_alignment) {
                // We couldn't align.
                if (graph_length == 0) {
                    // We had read sequence but no graph sequence.
                    // Try falling back to a pure insertion.
                    // TODO: We can be leaving the GBWT's space here!
                    // TODO: What if this is forcing an insertion that could also be in the graph already?
#ifdef debug_chain_alignment
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "connect() failed; treat as insertion" << endl;
                        }
                    }
#endif
                    link_alignment = WFAAlignment::make_unlocalized_insertion((*here).read_end(), link_length, aligner.score_gap(link_length));
                    link_alignment_source = "unlocalized_insertion";
                }
            } else if (link_alignment.length != linking_bases.size()) {
                // We could align, but we didn't get the alignment we expected. This shouldn't happen for a middle piece that can't softclip.
                stringstream ss;
                ss << "Aligning anchored link " << linking_bases << " (" << linking_bases.size() << " bp) from " << left_anchor << " - " << (*next).graph_start() << " against graph distance " << graph_length << " produced wrong-length alignment ";
                link_alignment.print(ss);
                throw ChainAlignmentFailedError(ss.str());
            } else {
                // We got the right alignment.
                // Put the alignment back into full read space
                link_alignment.seq_offset += (*here).read_end();
            }
        }
        
        if (link_alignment) {
            // We found a link alignment
            
#ifdef debug_chain_alignment
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Add link of length " << link_alignment.length << " with score of " << link_alignment.score << endl;
                }
            }
#endif
        
            link_alignment.check_lengths(gbwt_graph);
            
            // Then the link (possibly empty)
            append_path(composed_path, link_alignment.to_path(this->gbwt_graph, aln.sequence()));
            composed_score += link_alignment.score;
        } else {
            // The sequence to the next thing is too long, or we couldn't reach it doing connect().
            // Fall back to another alignment method
            
            if (linking_bases.size() > MAX_DP_LENGTH) {
                // This would be too long for GSSW to handle and might overflow 16-bit scores in its matrix.
#ifdef debug_chain_alignment
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Refusing to align " << link_length << " bp connection between chain items " << to_chain.backing_index(*here_it) << " and " << to_chain.backing_index(*next_it) << " which are " << graph_length << " apart at " << (*here).graph_end() << " and " << (*next).graph_start() << " in " << aln.name() << " to avoid overflow" << endl;
                }
#endif
                // Just jump to right tail
                break;
            }
           
#ifdef warn_on_fallback
            // We can't actually do this alignment, we'd have to align too
            // long of a sequence to find a connecting path.
            #pragma omp critical (cerr)
            {
                cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << link_length << " bp connection between chain items " << to_chain.backing_index(*here_it) << " and " << to_chain.backing_index(*next_it) << " which are " << graph_length << " apart at " << (*here).graph_end() << " and " << (*next).graph_start() << " in " << aln.name() << endl;
            }
#endif
            
            Alignment link_aln;
            link_aln.set_sequence(linking_bases);
            if (!aln.quality().empty()) {
                link_aln.set_quality(aln.quality().substr(link_start, link_length));
            }
            // Guess how long of a graph path we ought to allow in the alignment.
            size_t max_gap_length = longest_detectable_gap_in_range(aln, aln.sequence().begin() + link_start, aln.sequence().begin() + link_start + link_length, this->get_regular_aligner());
            size_t path_length = std::max(graph_length, link_length);
            if (stats) {
                start_time = std::chrono::high_resolution_clock::now();
            }
            auto nodes_and_bases = MinimizerMapper::align_sequence_between((*here).graph_end(), (*next).graph_start(), path_length, max_gap_length, &this->gbwt_graph, this->get_regular_aligner(), link_aln, &aln.name(), this->max_dp_cells, this->choose_band_padding);
            if (stats) {
                stop_time = std::chrono::high_resolution_clock::now();
                if (nodes_and_bases.first > 0) {
                    // Actually did the alignment
                    stats->bases.bga_middle += link_length;
                    stats->time.bga_middle += std::chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count();
                    stats->invocations.bga_middle += 1;
                }
            }
            link_alignment_source = "align_sequence_between";
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Add link of length " << path_to_length(link_aln.path()) << " with score of " << link_aln.score() << endl;
                }
            }
            
            // Then tack that path and score on
            append_path(composed_path, link_aln.path());
            composed_score += link_aln.score();
        }

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Aligned and added link of " << link_length << " via " << link_alignment_source << std::endl;
            }
        }
        
        // Advance here to next and start considering the next after it
        here_it = next_it;
        ++next_it;
        here = next;
    }
    
#ifdef debug_chain_alignment
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Add last extension " << *here_it << " of length " << (*here).length() << endl;
        }
    }
#endif
    
    WFAAlignment here_alignment = this->to_wfa_alignment(*here, aln, &aligner);

#ifdef debug_chain_alignment
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "\tScore " << here_alignment.score << endl;
        }
    }
#endif
    
    here_alignment.check_lengths(gbwt_graph);
    
    // Do the final GaplessExtension itself (may be the first)
    append_path(composed_path, here_alignment.to_path(this->gbwt_graph, aln.sequence()));
    composed_score += here_alignment.score;
   
    // Do the right tail, if any. Do as much of it as we can afford to do.
    size_t right_tail_length = aln.sequence().size() - (*here).read_end();
    if (right_tail_length > 0) {
        // We need to do a right tail
        
        string right_tail = aln.sequence().substr((*here).read_end(), right_tail_length);
        WFAAlignment right_alignment;
        // Grab the past-end graph position from the last thing in the chain. It is included in the tail as a base to align against.
        pos_t left_anchor_included = (*here).graph_end();
        // Pull back a base to get the outside-the-alignment anchoring position.
        pos_t left_anchor_excluded = left_anchor_included;
        get_offset(left_anchor_excluded)--;
        if (right_tail_length <= max_tail_length) {
            // We align the right tail with suffix(), which creates a suffix of the alignment.
            // Make sure to use the anchor outside of the region to be aligned.
            if (stats) {
                start_time = std::chrono::high_resolution_clock::now();
            }
            right_alignment = extender.suffix(right_tail, left_anchor_excluded);
            if (stats) {
                stop_time = std::chrono::high_resolution_clock::now();
                stats->bases.wfa_tail += right_tail_length;
                stats->time.wfa_tail += std::chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count();
                stats->invocations.wfa_tail += 1;
            }
        }
        
        if (right_alignment) {
            // Right tail did align. Put the alignment back into full read space.
            right_alignment.seq_offset += (*here).read_end();
            if (right_alignment.seq_offset + right_alignment.length != aln.sequence().size()) {
                // We didn't get all the way to the right end of the read without
                // running out of score.
                // Append a softclip.
                // TODO: Can we let the aligner know it can softclip for free?
                size_t right_end = right_alignment.seq_offset + right_alignment.length;
                size_t remaining = aln.sequence().size() - right_end;
                right_alignment.join(WFAAlignment::make_unlocalized_insertion(right_end, remaining, 0));
            }
            if (right_alignment.length != right_tail_length) {
                // We didn't get the alignment we expected.
                stringstream ss;
                ss << "Aligning right tail " << right_tail << " from " << left_anchor_excluded << " produced wrong-length alignment ";
                right_alignment.print(ss);
                throw ChainAlignmentFailedError(ss.str());
            }
#ifdef debug_chain_alignment
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Add right tail of " << right_tail.size() << " with score of " << right_alignment.score << endl;
                }
            }
#endif
            
            right_alignment.check_lengths(gbwt_graph);
            
            append_path(composed_path, right_alignment.to_path(this->gbwt_graph, aln.sequence()));
            composed_score += right_alignment.score;
        } else {
            // We need to fall back on alignment against the graph
            
#ifdef debug_chain_alignment
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "End with long right tail fallback alignment" << endl;
                }
            }
#endif

            if (right_tail.size() > MAX_DP_LENGTH) {
                // Right tail is too long to align.
               
#ifdef debug_chain_alignment
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Refusing to align " << right_tail.size() << " bp right tail against " << left_anchor_included << " in " << aln.name() << " to avoid overflow" << endl;
                }
#endif
                
                // Make a softclip for it.
                right_alignment = WFAAlignment::make_unlocalized_insertion((*here).read_end(), aln.sequence().size() - (*here).read_end(), 0);
                append_path(composed_path, right_alignment.to_path(this->gbwt_graph, aln.sequence()));
                composed_score += right_alignment.score;
            } else {

                Alignment tail_aln;
                tail_aln.set_sequence(right_tail);
                if (!aln.quality().empty()) {
                    tail_aln.set_quality(aln.quality().substr((*here).read_end(), right_tail_length));
                }

                // Work out how far the tail can see
                size_t max_gap_length = std::min(this->max_tail_gap, longest_detectable_gap_in_range(aln, aln.sequence().begin() + (*here).read_end(), aln.sequence().end(), this->get_regular_aligner()));
                size_t graph_horizon = right_tail_length + max_gap_length;

#ifdef warn_on_fallback
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << right_tail_length << " bp right tail against " << left_anchor_included << " allowing " << max_gap_length << " bp gap in " << aln.name() << endl;
                }
#endif

                // Align the right tail, anchoring the left end.
                // We need to use the included-in-the-alignment left anchor position.
                // TODO: What if it is past a node end? Is it guaranteed to be handled right?
                if (stats) {
                    start_time = std::chrono::high_resolution_clock::now();
                }
                auto nodes_and_bases = align_sequence_between(left_anchor_included, empty_pos_t(), graph_horizon, max_gap_length, &this->gbwt_graph, this->get_regular_aligner(), tail_aln, &aln.name(), this->max_dp_cells, this->choose_band_padding);
                if (stats) {
                    stop_time = std::chrono::high_resolution_clock::now();
                    if (nodes_and_bases.first > 0) {
                        // Actually did the alignment
                        stats->bases.dozeu_tail += right_tail_length;
                        stats->time.dozeu_tail += std::chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count();
                        stats->invocations.dozeu_tail += 1;
                    }
                }

                if (show_work && max_tail_length > 0) {
                    #pragma omp critical (cerr)
                    {
                        cerr << "warning[MinimizerMapper::find_chain_alignment]: Fallback score: " << tail_aln.score() << endl;
                    }
                }

                // Since it's the right tail we have to add it on
                append_path(composed_path, tail_aln.path());
                composed_score += tail_aln.score();
            }
        } 
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Aligned right tail length " << right_tail_length << std::endl;
            }
        }

    }
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Composed alignment is length " << path_to_length(composed_path) << " with score of " << composed_score << endl;
            if (composed_path.mapping_size() > 0) {
                cerr << log_name() << "Composed alignment starts with: " << pb2json(composed_path.mapping(0)) << endl;
                cerr << log_name() << "Composed alignment ends with: " << pb2json(composed_path.mapping(composed_path.mapping_size() - 1)) << endl;
            }
        }
    }
    
    // Convert to a vg Alignment.
    Alignment result(aln);
    // Simplify the path but keep internal deletions; we want to assert the
    // read deleted relative to some graph, and avoid jumps along nonexistent
    // edges.
    *result.mutable_path() = std::move(simplify(composed_path, false));
    result.set_score(composed_score);
    if (!result.sequence().empty()) {
        result.set_identity(identity(result.path()));
    }
    
    set_annotation(result, "left_tail_length", (double) left_tail_length);
    set_annotation(result, "longest_attempted_connection", (double) longest_attempted_connection); 
    set_annotation(result, "right_tail_length", (double) right_tail_length);
    
    return result;
}

void MinimizerMapper::wfa_alignment_to_alignment(const WFAAlignment& wfa_alignment, Alignment& alignment) const {
    *(alignment.mutable_path()) = wfa_alignment.to_path(this->gbwt_graph, alignment.sequence());
    alignment.set_score(wfa_alignment.score);
    if (!alignment.sequence().empty()) {
        alignment.set_identity(identity(alignment.path()));
    }
}

void MinimizerMapper::with_dagified_local_graph(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, const HandleGraph& graph, const std::function<void(DeletableHandleGraph&, const std::function<std::pair<nid_t, bool>(const handle_t&)>&)>& callback) {
    
    if (is_empty(left_anchor) && is_empty(right_anchor)) {
        throw ChainAlignmentFailedError("Cannot align sequence between two unset positions");
    }
    
    // We need to get the graph to align to.
    bdsg::HashGraph local_graph;
    unordered_map<id_t, id_t> local_to_base;
    if (!is_empty(left_anchor) && !is_empty(right_anchor)) {
        // We want a graph actually between two positions
        local_to_base = algorithms::extract_connecting_graph(
            &graph,
            &local_graph,
            max_path_length,
            left_anchor, right_anchor
        );
    } else if (!is_empty(left_anchor)) {
        // We only have the left anchor
        local_to_base = algorithms::extract_extending_graph(
            &graph,
            &local_graph,
            max_path_length,
            left_anchor,
            false,
            false
        );
    } else {
        // We only have the right anchor
        local_to_base = algorithms::extract_extending_graph(
            &graph,
            &local_graph,
            max_path_length,
            right_anchor,
            true,
            false
        );
    }

#ifdef debug
    std::cerr << "Local graph:" << std::endl;
    dump_debug_graph(local_graph);
#endif
    
    // To find the anchoring nodes in the extracted graph, we need to scan local_to_base.
    nid_t local_left_anchor_id = 0;
    nid_t local_right_anchor_id = 0;
    for (auto& kv : local_to_base) {
        if (kv.second == id(left_anchor) && kv.second == id(right_anchor)) {
            // The left and right anchors are on the same node, and this is a copy of it.
            // It could be that the anchors face each other, and we extracted one intervening piece of node.
            // In which case we go through this section once.
            if (local_left_anchor_id == 0 && local_right_anchor_id == 0) {
                // First time through, say we probably cut out the middle piece of a node
                local_left_anchor_id = kv.first;
                local_right_anchor_id = kv.first;
            } else {
                // Or it could be that we have two pieces of the original
                // shared node represented as separate nodes, because the
                // connecting path has to come back to the other end of this
                // shared node.
                //
                // In that case, we assume that extract_connecting_graph
                // assigns IDs so the start copy has a lower ID than the end
                // copy.
                if (local_left_anchor_id != local_right_anchor_id) {
                    // We thought we already figured out the start and end
                    // nodes; there are too many copies of our shared node to
                    // work out which is which.
                    std::stringstream ss;
                    ss << "Extracted graph from " << left_anchor;
                    if (!is_empty(right_anchor)) {
                        ss << " to " << right_anchor;
                    }
                    ss << " with max path length of " << max_path_length;
                    ss << " but shared node appeared more than twice in the resulting translation";
                    local_graph.serialize("crashdump.vg");
                    throw ChainAlignmentFailedError(ss.str());
                }
                // Whichever copy has the lower ID is the left one and
                // whichever copy has the higher ID is the right one.
                local_left_anchor_id = std::min(local_left_anchor_id, kv.first);
                local_right_anchor_id = std::max(local_right_anchor_id, kv.second);
            }
        } else if (kv.second == id(left_anchor)) {
            local_left_anchor_id = kv.first;
        } else if (kv.second == id(right_anchor)) {
            local_right_anchor_id = kv.first;
        }
        // TODO: Stop early when we found them all.
    }

    if (!is_empty(left_anchor) && local_left_anchor_id == 0) {
        #pragma omp critical (cerr)
        {
            for (auto& kv : local_to_base) {
                std::cerr << "Local ID " << kv.first << " = base graph ID " << kv.second << std::endl;
            }
        }
        // Somehow the left anchor didn't come through. Complain.
        std::stringstream ss;
        ss << "Extracted graph from " << left_anchor;
        if (!is_empty(right_anchor)) {
            ss << " to " << right_anchor;
        }
        ss << " with max path length of " << max_path_length;
        ss << " but from node was not present in the resulting translation";
        local_graph.serialize("crashdump.vg");
        throw ChainAlignmentFailedError(ss.str());
    }

    if (!is_empty(right_anchor) && local_right_anchor_id == 0) {
        // Somehow the right anchor didn't come through. Complain.
        std::stringstream ss;
        ss << "Extracted graph";
        if (!is_empty(left_anchor)) {
            ss << " from " << left_anchor;
        }
        ss << " to " << right_anchor;
        ss << " with max path length of " << max_path_length;
        ss << " but to node was not present in the resulting translation";
        local_graph.serialize("crashdump.vg");
        throw ChainAlignmentFailedError(ss.str());
    }
    
    // And split by strand since we can only align to one strand
    StrandSplitGraph split_graph(&local_graph);

#ifdef debug
    std::cerr << "Split graph:" << std::endl;
    dump_debug_graph(split_graph);
#endif
    
    // And make sure it's a DAG of the stuff reachable from our anchors
    bdsg::HashGraph dagified_graph;
    // For which we need the handles that anchor the graph, facing inwards
    std::vector<handle_t> bounding_handles;
    if (!is_empty(left_anchor)) {
        // Dagify from the forward version of the left anchor
        
        // Grab the left anchor in the local graph
        if (!local_graph.has_node(local_left_anchor_id)) {
            std::stringstream ss;
            ss << "Extracted graph from " << left_anchor;
            if (!is_empty(right_anchor)) {
                ss << " to " << right_anchor;
            }
            ss << " with max path length of " << max_path_length;
            ss << " but from node local ID " << local_left_anchor_id << " was not present in the resulting graph";
            local_graph.serialize("crashdump.vg");
            throw ChainAlignmentFailedError(ss.str());
        }
        handle_t local_handle = local_graph.get_handle(local_left_anchor_id, is_rev(left_anchor));
        
        // And get the node that that orientation of it is in the strand-split graph
        handle_t overlay_handle = split_graph.get_overlay_handle(local_handle);
        
        // And use that
        bounding_handles.push_back(overlay_handle);
    }
    if (!is_empty(right_anchor)) {
        // Dagify from the reverse version of the node for the forward version of the right anchor
        
        // Grab the right anchor from the local graph
        if (!local_graph.has_node(local_right_anchor_id)) {
            std::stringstream ss;
            ss << "Extracted graph";
            if (!is_empty(left_anchor)) {
                ss << " from " << left_anchor;
            }
            ss << " to " << right_anchor;
            ss << " with max path length of " << max_path_length;
            ss << " but to node local ID " << local_right_anchor_id << " was not present in the resulting graph";
            local_graph.serialize("crashdump.vg");
            throw ChainAlignmentFailedError(ss.str());
        }
        handle_t local_handle = local_graph.get_handle(local_right_anchor_id, is_rev(right_anchor));
        
        // And get the node that that orientation of it is in the strand-split graph
        // But flip it because we want to dagify going inwards from the right
        handle_t overlay_handle = split_graph.flip(split_graph.get_overlay_handle(local_handle));
        
        // And use that
        bounding_handles.push_back(overlay_handle);
    }
    
    auto dagified_to_split = handlegraph::algorithms::dagify_from(&split_graph, bounding_handles, &dagified_graph, max_path_length);
    
#ifdef debug
    std::cerr << "Dagified from " << bounding_handles.size() << " bounding handles in " << split_graph.get_node_count() << " node strand-split graph to " << dagified_graph.get_node_count() << " node DAG" << std::endl;
#endif

    // Make an accessor for getting back to the base graph space
    auto dagified_handle_to_base = [&](const handle_t& h) -> pair<nid_t, bool> {
        nid_t dagified_id = dagified_graph.get_id(h);
        bool dagified_is_reverse = dagified_graph.get_is_reverse(h);
        auto found_in_split = dagified_to_split.find(dagified_id);
        if (found_in_split == dagified_to_split.end()) {
            throw ChainAlignmentFailedError("ID " + std::to_string(dagified_id) + " from dagified graph not found in strand-split graph");
        }
        nid_t split_id = found_in_split->second;
        handle_t split_handle = split_graph.get_handle(split_id, dagified_is_reverse);
        // We rely on get_underlying_handle understanding reversed handles in the split graph
        handle_t local_handle = split_graph.get_underlying_handle(split_handle);
        nid_t local_id = local_graph.get_id(local_handle);
        bool local_is_reverse = local_graph.get_is_reverse(local_handle);
        auto found_in_base = local_to_base.find(local_id);
        if (found_in_base == local_to_base.end()) {
            throw ChainAlignmentFailedError("ID " + std::to_string(local_id) + " from local graph not found in full base graph");
        }
        nid_t base_id = found_in_base->second;
        return std::make_pair(base_id, local_is_reverse);
    };
    
    // Show the graph we made and the translation function
    callback(dagified_graph, dagified_handle_to_base);
}

size_t MinimizerMapper::longest_detectable_gap_in_range(const Alignment& aln, const std::string::const_iterator& sequence_begin, const std::string::const_iterator& sequence_end, const GSSWAligner* aligner) {
    
    // TODO: Should we take numbers and not iterators? This API could convert
    // better to quality adjustment later though.

    // If the range covers the middle, the longest detectable gap is the one from the middle.
    // TODO: Won't always be true anymore if we add quality adjustment
    size_t middle_index = aln.sequence().size() / 2;
    size_t begin_index = sequence_begin - aln.sequence().begin();
    size_t end_index = sequence_end - aln.sequence().begin();
    if (end_index > middle_index && begin_index <= middle_index) {
        return aligner->longest_detectable_gap(aln, aln.sequence().begin() + middle_index);
    }
    
    // Otherwise it is the length from the boundary nearest to the middle.
    // And we know the while range is on one side or the other of the middle.
    if (begin_index > middle_index) {
        // Beginning is on the inside
        return aligner->longest_detectable_gap(aln, sequence_begin);
    }

    // Otherwise the end is on the inside
    return aligner->longest_detectable_gap(aln, sequence_end);
}

std::pair<size_t, size_t> MinimizerMapper::align_sequence_between(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, size_t max_gap_length, const HandleGraph* graph, const GSSWAligner* aligner, Alignment& alignment, const std::string* alignment_name, size_t max_dp_cells, const std::function<size_t(const Alignment&, const HandleGraph&)>& choose_band_padding) {
    
    std::pair<size_t, size_t> to_return;

    // Get the dagified local graph, and the back translation
    MinimizerMapper::with_dagified_local_graph(left_anchor, right_anchor, max_path_length, *graph,
        [&](DeletableHandleGraph& dagified_graph, const std::function<std::pair<nid_t, bool>(const handle_t&)>& dagified_handle_to_base) {

//#ifdef debug
        dump_debug_graph(dagified_graph);
//#endif
        
        // Then trim off the tips that are either in the wrong orientation relative
        // to whether we want them to be a source or a sink, or extraneous
        
        std::vector<handle_t> tip_handles = handlegraph::algorithms::find_tips(&dagified_graph);
        bool trimmed;
        size_t trim_count = 0;
        do {
            trimmed = false;
            // We need to make sure to remove only one orientation of each handle
            // we remove.
            std::unordered_set<nid_t> to_remove_ids;
            std::vector<handle_t> to_remove_handles;
            for (auto& h : tip_handles) {
                auto base_coords = dagified_handle_to_base(h);
                if (!dagified_graph.get_is_reverse(h) && (is_empty(left_anchor) || base_coords.first == id(left_anchor))) {
                    // Tip is inward forward, so it's a source.
                    // This is a head in the subgraph, and either matches a left
                    // anchoring node or we don't have any, so keep it.
#ifdef debug
                    std::cerr << "Dagified graph node " << dagified_graph.get_id(h) << " " << dagified_graph.get_is_reverse(h) << " is an acceptable source (" << base_coords.first << " " << base_coords.second << ")" << std::endl;
#endif
                } else if (dagified_graph.get_is_reverse(h) && (is_empty(right_anchor) || base_coords.first == id(right_anchor))) {
                    // Tip is inward reverse, so it's a sink.
                    // This is a tail in the subgraph, and either matches a right
                    // anchoring node or we don't have any, so keep it.
#ifdef debug
                    std::cerr << "Dagified graph node " << dagified_graph.get_id(h) << " " << dagified_graph.get_is_reverse(h) << " is an acceptable sink (" << base_coords.first << " " << base_coords.second << ")" << std::endl;
#endif
                } else {
                    // This is a wrong orientation of an anchoring node, or some other tip.
                    // We don't want to keep this handle
#ifdef debug
                    std::cerr << "Dagified graph node " << dagified_graph.get_id(h) << " " << dagified_graph.get_is_reverse(h) << " is an unacceptable tip (" << base_coords.first << " " << base_coords.second << ")" << std::endl;
#endif
                    nid_t dagified_id = dagified_graph.get_id(h);
                    if (!to_remove_ids.count(dagified_id)) {
                        to_remove_ids.insert(dagified_id);
                        to_remove_handles.push_back(h);
                    }
                }
            }
            for (auto& h : to_remove_handles) {
                dagified_graph.destroy_handle(h);
                trimmed = true;
            }
            if (trimmed) {
                // TODO: This is going to be O(slow) if we actually have to
                // prune back a dangling run. We should look at what is
                // connected to the tip and the tip only, and make that the new
                // tip. Or keep some kind of online tip info. Or use an
                // algorithm function that we make actually good.
                tip_handles = handlegraph::algorithms::find_tips(&dagified_graph);
                trim_count++;

#ifdef debug
                std::cerr << "Dagified graph trim " << trim_count << ":" << std::endl;
                dump_debug_graph(dagified_graph);
#endif
            }
        } while (trimmed);
        if (trim_count > 0) {
            #pragma omp critical (cerr)
            {
                std::cerr << "warning[MinimizerMapper::align_sequence_between]: Trimmed back tips " << trim_count << " times on graph between " << left_anchor << " and " << right_anchor << " leaving " <<  dagified_graph.get_node_count() << " nodes and " << tip_handles.size() << " tips";
                if (alignment_name) {
                    std::cerr << " for read " << *alignment_name;
                }
                std::cerr << std::endl;
            }
        }
        
        if (!is_empty(left_anchor) && !is_empty(right_anchor)) {
            // Then align the linking bases, with global alignment so they have
            // to go from a source to a sink. Banded alignment means we can
            // safely do big problems.
            //
            // We need to pick band padding based on what we are aligning, and
            // we want to use permissive banding.
            size_t band_padding = choose_band_padding(alignment, dagified_graph);
#ifdef debug
            std::cerr << "Aligning with band padding: " << band_padding << " for alignment length " << alignment.sequence().size() << std::endl;
#endif
            aligner->align_global_banded(alignment, dagified_graph, band_padding, true);
            to_return.first = dagified_graph.get_node_count();
            to_return.second = dagified_graph.get_total_length();
        } else {
            // Do pinned alignment off the anchor we actually have.
            // Work out how big it will be.
            size_t cell_count = dagified_graph.get_total_length() * alignment.sequence().size();
            if (cell_count > max_dp_cells) {
                #pragma omp critical (cerr)
                {
                    std::cerr << "warning[MinimizerMapper::align_sequence_between]: Refusing to fill " << cell_count << " DP cells in tail with Xdrop";
                    if (alignment_name) {
                        std::cerr << " for read " << *alignment_name;
                    }
                    std::cerr << std::endl;
                }
                // Fake a softclip right in input graph space
                alignment.clear_path();
                Mapping* m = alignment.mutable_path()->add_mapping();
                // TODO: Is this fake position OK regardless of anchoring side?
                m->mutable_position()->set_node_id(is_empty(left_anchor) ? id(right_anchor) : id(left_anchor));
                m->mutable_position()->set_is_reverse(is_empty(left_anchor) ? is_rev(right_anchor) : is_rev(left_anchor));
                m->mutable_position()->set_offset(is_empty(left_anchor) ? offset(right_anchor) : offset(left_anchor));
                Edit* e = m->add_edit();
                e->set_to_length(alignment.sequence().size());
                e->set_sequence(alignment.sequence());
                to_return.first = 0;
                to_return.second = 0;
                return;
            } else {
#ifdef debug
                #pragma omp critical (cerr)
                std::cerr << "debug[MinimizerMapper::align_sequence_between]: Fill " << cell_count << " DP cells in tail with Xdrop" << std::endl;
#endif
                aligner->align_pinned(alignment, dagified_graph, !is_empty(left_anchor), true, max_gap_length);
                to_return.first = dagified_graph.get_node_count();
                to_return.second = dagified_graph.get_total_length();
            }
        }
        
        // And translate back into original graph space
        for (size_t i = 0; i < alignment.path().mapping_size(); i++) {
            // Translate each mapping's ID and orientation down to the base graph
            Mapping* m = alignment.mutable_path()->mutable_mapping(i);
            
            handle_t dagified_handle = dagified_graph.get_handle(m->position().node_id(), m->position().is_reverse());
            auto base_coords = dagified_handle_to_base(dagified_handle); 
            
            m->mutable_position()->set_node_id(base_coords.first);
            m->mutable_position()->set_is_reverse(base_coords.second);
        }
        if (!is_empty(left_anchor) && alignment.path().mapping_size() > 0 && offset(left_anchor) != 0 && offset(left_anchor) < graph->get_length(graph->get_handle(id(left_anchor)))) {
            // There is some of the left anchor's node actually in the
            // extracted graph. The left anchor isn't past the end of its node.
            
            // Get the positions of the leftmost mapping
            Position* left_pos = alignment.mutable_path()->mutable_mapping(0)->mutable_position();

            // The alignment must actually start on the anchor node.
            assert(left_pos->node_id() == id(left_anchor));

            // Add on the offset for the missing piece of the left anchor node
            left_pos->set_offset(left_pos->offset() + offset(left_anchor));
        }
        if (alignment.path().mapping_size() > 0) {
            // Make sure we don't have an empty mapping on the end
            auto* last_mapping = alignment.mutable_path()->mutable_mapping(alignment.path().mapping_size() - 1);
            if (last_mapping->edit_size() > 0) {
                // Make sure we don't have an empty edit on the end
                auto& last_edit = last_mapping->edit(last_mapping->edit_size() - 1);
                if (last_edit.from_length() == 0 && last_edit.to_length() == 0 && last_edit.sequence().empty()) {
                    // Last edit is empty so drop from the mapping
                    last_mapping->mutable_edit()->RemoveLast();
                }
            }
            if (last_mapping->edit_size() == 0) {
                // Last mapping is empty, so drop it.
                alignment.mutable_path()->mutable_mapping()->RemoveLast();
            }
        }
    
        // Now the alignment is filled in!
    });

    return to_return;
}

std::vector<algorithms::Anchor> MinimizerMapper::to_anchors(const Alignment& aln, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds) const {
    std::vector<algorithms::Anchor> to_return;
    to_return.reserve(seeds.size());
    for (size_t i = 0; i < seeds.size(); i++) {
        to_return.push_back(MinimizerMapper::to_anchor(aln, minimizers, seeds, i, gbwt_graph, get_regular_aligner()));
    }
    return to_return;
}

algorithms::Anchor MinimizerMapper::to_anchor(const Alignment& aln, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, size_t seed_number, const HandleGraph& graph, const Aligner* aligner) {
    // Turn each seed into the part of its match on the node where the
    // anchoring end (start for forward-strand minimizers, end for
    // reverse-strand minimizers) falls.
    auto& seed = seeds[seed_number];
    auto& source = minimizers[seed.source];
    size_t length;
    pos_t graph_start;
    size_t read_start;
    size_t hint_start;
    size_t margin_left;
    size_t margin_right;
    if (source.value.is_reverse) {
        // Seed stores the final base of the match in the graph.
        // So get the past-end position.
        pos_t graph_end = make_pos_t(id(seed.pos), is_rev(seed.pos), offset(seed.pos) + 1);
        
        // Work out how much of the node it could use before there.
        length = std::min((size_t) source.length, offset(graph_end));
        // And how much we cut off the start
        margin_left = (size_t)source.length - length;
        // We cut nothing off the end
        margin_right = 0;
        // And derive the graph start
        graph_start = make_pos_t(id(graph_end), is_rev(graph_end), offset(graph_end) - length);
        // And the read start
        read_start = source.value.offset + 1 - length;
        // The seed is actually the last 1bp interval
        hint_start = length - 1;
    } else {
        // Seed stores the first base of the match in the graph
        graph_start = seed.pos;
        
        // Get the handle to the node it's on.
        handle_t start_handle = graph.get_handle(id(graph_start), is_rev(graph_start));
        // Work out how much of the node it could use before there.
        length = std::min((size_t) source.length, graph.get_length(start_handle) - offset(graph_start));
        // We cut nothing off the start
        margin_left = 0;
        // How much do we cut off the end?
        margin_right = (size_t)source.length - length;
        // And we store the read start position already in the item
        read_start = source.value.offset;
        // The seed is actually at the start
        hint_start = 0;
    }

#ifdef debug
    std::cerr << "Minimizer at read " << source.forward_offset() << " length " << source.length
              << " orientation " << source.value.is_reverse << " pinned at " << source.value.offset
              << " is anchor of length " << length << " matching graph " << graph_start << " and read " << read_start
              << " forward, with hint " << hint_start << " bases later on the read" << std::endl;
#endif

    // Work out how many points the anchor is.
    // TODO: Always make sequence and quality available for scoring!
    // We're going to score the anchor as the full minimizer, and rely on the margins to stop us from taking overlapping anchors.
    int score = aligner->score_exact_match(aln, read_start - margin_left, length + margin_right);
    return algorithms::Anchor(read_start, graph_start, length, margin_left, margin_right, score, seed_number, seed.zipcode_decoder.get(), hint_start); 
}

algorithms::Anchor MinimizerMapper::to_anchor(const Alignment& aln, size_t read_start, size_t read_end, const std::vector<size_t>& sorted_seeds, const std::vector<algorithms::Anchor>& seed_anchors, const std::vector<size_t>::const_iterator& mismatch_begin, const std::vector<size_t>::const_iterator& mismatch_end, const HandleGraph& graph, const Aligner* aligner) {
    if (sorted_seeds.empty()) {
        // This should never happen
        throw std::runtime_error("Can't make an anchor from no seeds");
    }

    // Score all the matches and mismatches.
    int score = 0;
    size_t scored_until = read_start;
    auto mismatch_it = mismatch_begin;
    while(mismatch_it != mismatch_end) {
        // Score the perfect match up to mismatch_it, and the mismatch at mismatch_it.
        score += aligner->score_exact_match(aln, scored_until, *mismatch_it - scored_until);
        score += aligner->score_mismatch(aln.sequence().begin() + *mismatch_it,
                                         aln.sequence().begin() + *mismatch_it + 1,
                                         aln.quality().begin() + *mismatch_it); 
        scored_until = *mismatch_it + 1;
        ++mismatch_it;
    }
    // Score the perfect match from where we are to the end.
    score += aligner->score_exact_match(aln, scored_until, read_end - scored_until);
    
    // Get the anchors we are going to weld together. These may be the same one.
    const algorithms::Anchor& left_anchor = seed_anchors.at(sorted_seeds.front());
    const algorithms::Anchor& right_anchor = seed_anchors.at(sorted_seeds.back());

    // Work out the additional left and right margin we need to block out other
    // overlapping extensions and justify our score. The range can extend
    // beyond even the outermost minimizers.
    size_t extra_left_margin = left_anchor.read_exclusion_start() - read_start;
    size_t extra_right_margin = read_end - right_anchor.read_exclusion_end();

    // Now make an anchor with the score of the range, with the anchors of
    // the first and last seeds, and enough margin to cover the distance out
    // from the outer seeds that we managed to extend.
    algorithms::Anchor result(left_anchor, right_anchor, extra_left_margin, extra_right_margin, score);

    assert(result.read_exclusion_start() == read_start);
    assert(result.read_exclusion_end() == read_end);

    return result;
}

WFAAlignment MinimizerMapper::to_wfa_alignment(const algorithms::Anchor& anchor, const Alignment& aln, const Aligner* aligner) const {
    // Get the score without full length bonuses
    auto score = aligner->score_exact_match(aln, anchor.read_start(), anchor.length());
    if (anchor.read_start() == 0) {
        // Apply full elngth bonus on the left if we abut the left end of the read.
        score += aligner->score_full_length_bonus(true, aln);
    }
    if (anchor.read_end() == aln.sequence().length()) {
        // Apply full lenght bonus on the right if we abut the riht end of the read.
        score += aligner->score_full_length_bonus(false, aln);
    }

    return {
        {gbwt_graph.get_handle(id(anchor.graph_start()), is_rev(anchor.graph_start()))},
        {{WFAAlignment::match, (uint32_t)anchor.length()}},
        (uint32_t)offset(anchor.graph_start()),
        (uint32_t)anchor.read_start(),
        (uint32_t)anchor.length(),
        score,
        true
    };
}

}

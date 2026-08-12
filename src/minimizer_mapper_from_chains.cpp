/**
 * \file minimizer_mapper_from_chains.cpp
 * Defines the code for the long-read code path for the
 * minimizer-and-GBWT-based mapper (long read Giraffe).
 */

#include "minimizer_mapper.hpp"

#include "logged_gap_alignment_scorer.hpp"

#include "annotation.hpp"
#include "banded_global_aligner.hpp"
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
#include <unordered_set>
#include <bitset>

// Turn on debugging prints
//#define debug
// Turn on recombintion debugging prints
//#define debug_rec
// Turn on printing of minimizer fact tables
//#define print_minimizer_table
// Dump the zip code forest
//#define debug_print_forest
// Dump local graphs during tip trimming
//#define debug_dump_tips
// Dump local graphs that we align against 
//#define debug_dump_graph
// Dump fragment length distribution information
//#define debug_fragment_distr
//Do a brute force check that clusters are correct
//#define debug_validate_clusters
//#define debug_write_minimizers
// Debug generation of alignments from chains
//#define debug_base_level_alignment

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

void MinimizerMapper::dump_debug_chains(const ZipCodeForest& zip_code_forest,
                                        const std::vector<Seed>& seeds,
                                        const VectorView<Minimizer>& minimizers,
                                        const vector<algorithms::Anchor>& seed_anchors,
                                        const std::vector<algorithms::SubchainGroup>& subchain_groups,
                                        const std::vector<size_t>& subchain_source_tree,
                                        const PathPositionHandleGraph* path_graph,
                                        bool haplotype_positions) {
    if (!path_graph) {
        // We don't have a path positional graph for this
        return;
    }

    // Loop through all trees' chaining results
    size_t overall_index = 0;
    for (size_t group_num = 0; group_num < subchain_groups.size(); group_num++) {
        for (size_t chain_num = 0; chain_num < subchain_groups.at(group_num).subchains.size(); chain_num++) {
            // For each chain, create a separate TSV file
            std::string cur_name = std::to_string(group_num) + "chain" + std::to_string(chain_num);

            auto& tree_num = subchain_source_tree.at(overall_index);
            overall_index++;

            // Find all the seeds in its zip tree
            vector<size_t> involved_seeds;
            for (ZipCodeTree::oriented_seed_t found : zip_code_forest.trees.at(tree_num).get_all_seeds()) {
                involved_seeds.push_back(found.seed);
            }

            // Start making a list of things to show.
            std::vector<std::pair<std::string, std::vector<std::vector<size_t>>>> seed_sets;
            seed_sets.emplace_back("", std::vector<std::vector<size_t>>{std::move(involved_seeds)});
            seed_sets.emplace_back("chain", std::vector<std::vector<size_t>>{subchain_groups[group_num].subchains.at(chain_num).anchors});

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

            // Create a name for this chain's file
            std::string chain_file_name = cur_name + "-dotplot";

            // Now we need to actually dump the data to a TSV
            // The read context should already be set by the caller
            TSVExplainer exp(true, chain_file_name);

            // We make another TSV that's more parseable, with all the seeds.
            TSVExplainer seedpos(true, cur_name + "-seeds");

            // Determine the positions of all the involved seeds.
            std::unordered_map<size_t, algorithms::path_offset_collection_t> seed_positions;
            std::unordered_set<PathSense> wanted_senses {PathSense::REFERENCE, PathSense::GENERIC};
            if (haplotype_positions) {
                wanted_senses.insert(PathSense::HAPLOTYPE);
            }
            for (auto& kv : seed_sets) {
                for (const std::vector<size_t> included_seeds : kv.second) {
                    for (auto& seed_num : included_seeds) {
                        // For each seed in the run
                        auto& seed = seeds.at(seed_num);

                        auto found = seed_positions.find(seed_num);
                        if (found == seed_positions.end()) {
                            // If we don't know the seed's positions yet, get them.
                            // We are working with the *pin point* (seed pos in the
                            // graph and minimizer pin_offset() in the read), not
                            // anything to do with the anchor.
                            
                            // Find that in the graph, on paths.
                            found = seed_positions.emplace_hint(found, seed_num, algorithms::nearest_offsets_in_paths(path_graph, seed.pos, 100, wanted_senses));
                            for (auto& handle_and_positions : found->second) {
                                std::string path_name = path_graph->get_path_name(handle_and_positions.first);
                                for (auto& position : handle_and_positions.second) {
                                    // Dump all the seed positions so we can select seeds we want to know about.
                                    // These are used with scripts/make-chain-viz.py to make interactive chaining problem visualizations.
                                    seedpos.line();
                                    seedpos.field(minimizers[seed.source].pin_offset());
                                    seedpos.field(path_name);
                                    seedpos.field(position.first);
                                    seedpos.field(position.second ? "-" : "+");
                                    seedpos.field(seed_num);
                                    std::stringstream ss;
                                    ss << seed_anchors.at(seed_num);
                                    seedpos.field(ss.str());
                                }
                            }
                            if (found->second.empty()) {
                                // The seed doesn't have any linear positions, but might still participate in the winning chain traceback.
                                // Report it.
                                seedpos.line();
                                seedpos.field(minimizers[seed.source].pin_offset());
                                seedpos.field("");
                                seedpos.field("");
                                seedpos.field("");
                                seedpos.field(seed_num);
                                std::stringstream ss;
                                ss << seed_anchors.at(seed_num);
                                seedpos.field(ss.str());
                            }
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
                    for (size_t idx = 0; idx < included_seeds.size(); ++idx) {
                        // For each seed in the run (index-based so we can consult chain flags)
                        size_t seed_num = included_seeds[idx];
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
                                // Offset on contig of the pin point
                                exp.field(position.first);
                                // Offset in read *of the pin point* (not of the forward-strand start of the minimizer)
                                exp.field(minimizers[seed.source].pin_offset());
                            }
                        }
                        if (offsets.empty()) {
                            // Note that we don't actually have a position
                            exp.line();
                            if (!marker.empty()) {
                                // Sentinel and a marker and a subscript
                                exp.field("NO_PATH-" + marker + "-" + std::to_string(run_number));
                            } else {
                                // Sentinel alone
                                exp.field("NO_PATH");
                            }
                            // Put it at 0 on no path
                            exp.field(0);
                            // Offset in read *of the pin point* (not of the forward-strand start of the minimizer)
                            exp.field(minimizers[seed.source].pin_offset());
                        }
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
    for (ZipCodeTree::oriented_seed_t found : zip_code_forest.trees[i].get_all_seeds()) {
        if (this->track_provenance) {
            // Remember the seeds
            tree_seeds.push_back(found.seed);
        }
        // For each seed in the tree, find what minimizer it comes from
        if (found.seed >= seeds.size()) {
            throw std::out_of_range("Tree " + std::to_string(i) + " has seed " + std::to_string(found.seed)
                                    + " but we only have " + std::to_string(seeds.size()) + " seeds");
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
            if (!tree_positions.empty()) {
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
    
    Explainer::set_context(aln.name());

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

    // Create a new alignment object to get rid of old annotations.
    {
      Alignment temp;
      temp.set_sequence(aln.sequence());
      temp.set_name(aln.name());
      temp.set_quality(aln.quality());
      if (has_annotation(aln, "tags")) {
        // Preserve any BAM tags, which might really have come in as FASTQ
        // comments and which we want to keep.
        // TODO: What if these came from a previous Giraffe run though???
        set_annotation(temp, "tags", get_annotation<string>(aln, "tags"));
      }
      aln = std::move(temp);
    }

    // Annotate the read with metadata
    if (!sample_name.empty()) {
        aln.set_sample_name(sample_name);
    }
    if (!read_group.empty()) {
        aln.set_read_group(read_group);
    }

    // Minimizers sorted by position
    std::vector<Minimizer> minimizers_in_read = this->find_minimizers(aln.sequence(), funnel);
    // Flag minimizers as being in repetitive regions of the read or not
    this->flag_repetitive_minimizers(minimizers_in_read);
    // Indexes of minimizers, sorted into score order, best score first
    std::vector<size_t> minimizer_score_order = sort_minimizers_by_score(minimizers_in_read, rng);
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
    zip_code_forest.fill_in_forest(seeds, *distance_index, aln.sequence().size() * zipcode_tree_scale);

#ifdef debug_print_forest
    if (show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "Zip code forest:" << std::endl;
            zip_code_forest.print_self(&seeds);
        }
    }
#endif

    // Turn all the seeds into anchors. Either we'll chain them directly or
    // use them to make gapless extension anchors over them.
    // TODO: Can we only use the seeds that are in trees we keep?
    vector<algorithms::Anchor> seed_anchors = this->to_anchors(aln, minimizers, seeds);
    
    // If we do gapless extension, then it is possible to find full-length gapless extensions at this stage
    // If we have at least one good gapless extension, then we will turn them directly into alignments
    // and skip the later stages. Store alignments from gapless extensions here

    // We will fill this with all computed alignments in estimated score order
    std::vector<Alignment> alignments;
    //The multiplicity for each alignment, projected from previous stages
    vector<double> multiplicity_by_alignment;
    // Track if minimizers were explored by alignments
    SmallBitset minimizer_explored(minimizers.size());

    // For each initial chaining DP, we need:
    // The subchains and connections between them
    std::vector<algorithms::SubchainGroup> subchain_groups;
    // The zip code tree it came from
    std::vector<size_t> subchain_source_tree;
    // A count, for each minimizer, of how many hits of it could have been in the chain, or were considered when making the chain.
    std::vector<std::vector<size_t>> minimizer_kept_chain_count;
    // The multiplicity for each chain. For now, just the multiplicity of the tree it came from
    std::vector<double> multiplicity_by_chain;

    do_chaining_on_trees(aln, zip_code_forest, seeds, minimizers, seed_anchors,
                         subchain_groups, subchain_source_tree,
                         minimizer_kept_chain_count, multiplicity_by_chain,
                         alignments, minimizer_explored, multiplicity_by_alignment,
                         rng, funnel);

    //Fill in chain stats for annotating the final alignment
    bool best_chain_correct = false;
    double best_chain_coverage = 0;
    size_t best_chain_longest_jump = 0;
    double best_chain_average_jump = 0;
    size_t best_chain_anchors = 0;
    size_t best_chain_anchor_length = 0;

    // Dump all chains if requested (do this before alignments, while chains still exist)
    if (show_work && !subchain_groups.empty() && this->path_graph != nullptr) {
        dump_debug_chains(zip_code_forest, seeds, minimizers, seed_anchors, subchain_groups, subchain_source_tree, this->path_graph, this->haplotype_positions);
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "=====Creating alignments=====" << endl;
        }
    }

    
    // Now start the alignment step. Everything has to become an alignment.

    // Track statistics about how many bases were aligned by diffrent methods, and how much time was used.
    aligner_stats_t stats; 

    // This maps from alignment index back to chain index, for
    // tracing back to minimizers for MAPQ. Can hold
    // numeric_limits<size_t>::max() for an unaligned alignment.
    vector<size_t> alignments_to_source;
    alignments_to_source.reserve(subchain_groups.size());

    if (alignments.size() == 0) {
        do_alignment_on_chains(aln, seeds, minimizers, seed_anchors, subchain_groups,
                               subchain_source_tree, multiplicity_by_chain,
                               minimizer_kept_chain_count, alignments, multiplicity_by_alignment, 
                               alignments_to_source, minimizer_explored, stats, rng, funnel);
    }
    
    if (track_provenance) {
        // Now say we are finding the winner(s)
        funnel.stage("winner");
    }
    
    // Fill this in with the alignments we will output as mappings
    vector<Alignment> mappings;
    mappings.reserve(min(alignments.size(), max_multimaps));
    //The scores of the mappings
    vector<double> scores;
    //The multiplicities of mappings
    vector<double> multiplicity_by_mapping;
   
    // Collect the chosen mappings, or create the unmapped-read mapping.
    pick_mappings_from_alignments(aln, alignments, multiplicity_by_alignment, alignments_to_source,
                                  mappings, scores, multiplicity_by_mapping, rng, funnel);
    
    if (track_provenance) {
        funnel.substage("mapq");
    }

    // Note that it is possible for the top base-level alignment score *not* to be the winning alignment!

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Picked best alignment " << log_alignment(mappings[0]) << endl;
            cerr << log_name() << "For scores:";
            for (size_t i = 0; i < scores.size(); i++) {
                cerr << " " << scores[i];
                if (i + 1 < scores.size()) {
                    cerr << ",";
                }
            }
            cerr << endl;
        }
    }

    vector<double> scaled_scores;
    scaled_scores.reserve(scores.size());
    for (auto& score : scores) {
        double scaled_score = score;
        if (mapq_score_window > 0) {
            // Rescale to the size of the score window
            scaled_score = scaled_score * mapq_score_window / aln.sequence().size();
        }
        // Rescale by a constant factor
        scaled_score *= mapq_score_scale;
        scaled_scores.push_back(scaled_score);
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Scaled scores:";
            for (size_t i = 0; i < scaled_scores.size(); i++) {
                cerr << " " << scaled_scores[i];
                if (i + 1 < scaled_scores.size()) {
                    cerr << ",";
                }
            }
            cerr << endl;
        }
    }

    crash_unless(!mappings.empty());
    // Compute MAPQ if not unmapped. Otherwise use 0 instead of the 50% this would give us.
    // Use exact mapping quality.
    // Because the winning alignment won't necessarily *always* have the
    // maximum score, we need to use compute_first_mapping_quality and not
    // compute_max_mapping_quality.
    double mapq = (mappings.front().path().mapping_size() == 0) ? 0 : 
        get_regular_aligner()->mapq_calc->compute_first_mapping_quality(scaled_scores, false, &multiplicity_by_alignment) ;

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

    if (track_provenance) {
        funnel.stage("demapping");
    }

    if (mapq == 0 && !scores.empty() && scores.front() < min_mapq0_score) {
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Failing MAPQ 0 alignment for having top score "
                     << scores.front() << " which is below " << min_mapq0_score << endl;
            }
        }
        if (track_provenance) {
            // Fail all remaining mappings
            for (size_t i = 0; i < mappings.size(); i++) {
                funnel.fail("mapq0-score", i, scores.front());
            }
        }

        // Reset scores / mappings
        scores.clear();
        mappings.clear();

        scores.emplace_back(0);
        mappings.emplace_back(aln);
    } else if (track_provenance) {
        // Pass all remaining mappings
        for (size_t i = 0; i < mappings.size(); i++) {
            funnel.pass("mapq0-score", i);
            funnel.project(i);
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
    
    if (track_provenance && track_correctness) {
        annotate_with_minimizer_statistics(mappings[0], minimizers, seeds, seeds.size(), subchain_groups.size(), funnel);
    }

    // Special chain statistics
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

    Explainer::clear_context();

    return mappings;
}

void MinimizerMapper::do_chaining_on_trees(const Alignment& aln, const ZipCodeForest& zip_code_forest,
    const std::vector<Seed>& seeds, const VectorView<MinimizerMapper::Minimizer>& minimizers,
    const vector<algorithms::Anchor>& seed_anchors,
    std::vector<algorithms::SubchainGroup>& subchain_groups, std::vector<size_t>& subchain_source_tree,
    std::vector<std::vector<size_t>>& minimizer_kept_chain_count,
    std::vector<double>& multiplicity_by_chain,
    std::vector<Alignment>& alignments, SmallBitset& minimizer_explored, vector<double>& multiplicity_by_alignment,
    LazyRNG& rng, Funnel& funnel) const {

    // Keep track of which chain each alignment comes from for the funnel
    std::vector<size_t> alignment_source_chain;

    // After going through all trees and counting how many are kept,
    // each value will be divided by the number of trees kept
    size_t kept_tree_count = 0;

    //Do gapless extension if the read length is less than the limit
    bool do_gapless_extension = aln.sequence().size() <= gapless_extension_limit;

    // First score all the zip code trees in the forest by summing the scores of their involved minimizers.
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
            std::cerr << log_name() << "Found "
                      << zip_code_forest.trees.size() << " zip code trees, scores "
                      << best_tree_score << " best, "
                      << second_best_tree_score << " second best, coverages "
                      << best_tree_coverage << " best, "
                      << second_best_tree_coverage << " second best" << std::endl;
        }
    }



    if (track_provenance) {
        funnel.stage("chain");
        funnel.substage("chain");
    }

    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "=====Creating chains=====" << endl;
        }
    }

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

    process_until_threshold_c<double>(zip_code_forest.trees.size(), [&](size_t i) -> double {
            return tree_coverages[i];
        }, [&](size_t a, size_t b) -> bool {
            auto equalish = [&] (const double x, const double y) {
                if (x == y) {
                    return true;
                } else if (x > y) {
                    return x - y <= std::numeric_limits<double>::round_error();
                } else {
                    return y - x <= std::numeric_limits<double>::round_error();
                }
            };
            auto greater_than = [&] (const double x, const double y) {
                if (equalish(x, y)) {
                    return false;
                } else {
                    return x > y;
                }
            };

            return greater_than(tree_coverages[a], tree_coverages[b])
                || (equalish(tree_coverages[a], tree_coverages[b]) && greater_than(tree_scores[a], tree_scores[b]));

        }, this->zipcode_tree_coverage_threshold, this->min_chaining_problems, this->max_chaining_problems, rng, [&](size_t item_num, size_t item_count) -> bool {
            // Handle sufficiently good chaining problems in descending score order

            if (track_provenance) {
                funnel.pass("zipcode-tree-coverage-threshold", item_num, tree_coverages[item_num]);
                funnel.pass("max-chaining-problems", item_num);
            }

            // First check against the additional score filter
            if (zipcode_tree_score_threshold != 0 && tree_scores[item_num] < tree_score_cutoff
                && kept_tree_count >= min_chaining_problems) {
                // If the score isn't good enough and we already kept at least min_chaining_problems trees,
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
                    cerr << log_name() << "Making chains for zip code tree " << item_num 
                         << " with score " << tree_scores[item_num]
                         << " and coverage " << tree_coverages[item_num] << endl;
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
            for (ZipCodeTree::oriented_seed_t found : zip_code_forest.trees[item_num].get_all_seeds()) {
                selected_seeds.push_back(found.seed);
                added_seed[found.seed] = true;
            }

            if (show_work) {
                dump_debug_seeds(minimizers, seeds, selected_seeds);
            }

            // If we do gapless extension, we will use these anchors to chain instead of the seed ones.
            std::vector<algorithms::Anchor> extension_anchors;
            // And each of them (or of the seed anchors, if we use those) represents this run of seed numbers to put into the final chain.
            std::vector<std::vector<size_t>> extension_seed_sequences;
            // Extensions use a distinct list of included seeds vs. seeds we actually paste in, so we can glom up overlapping seeds.
            std::vector<std::vector<size_t>> extension_represented_seeds;
            // We need a list of all extension anchor indexes that we can sort.
            std::vector<size_t> extension_anchor_indexes;

            if (do_gapless_extension) {
                // Instead of chaining directly on the seeds, chain on gapless extensions of the seeds.

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

                //If there are full-length extensions that are good enough, then just turn them into alignments.
                if (GaplessExtender::full_length_extensions(tree_extensions)) {
                    for (size_t extension_i = 0 ; extension_i < tree_extensions.size() ; extension_i++) {
                        if (tree_extensions[extension_i].full() &&
                            tree_extensions[extension_i].mismatches() <= this->default_max_extension_mismatches) {

                            // For all good-scoring full-length extensions, make them into alignments
                            // TODO When we pair:
                            // We want them all to go on to the pairing stage so we don't miss a possible pairing in a tandem repeat.

                            alignments.emplace_back(aln);
                            alignments.back().clear_refpos();
                            alignments.back().clear_path();
                            alignments.back().set_score(0);
                            alignments.back().set_identity(0);
                            alignments.back().set_mapping_quality(0);
                            this->extension_to_alignment(tree_extensions[extension_i], alignments.back());

                            if (track_provenance) {
                                //We want to know which "chain" this came from
                                alignment_source_chain.emplace_back(subchain_groups.size());
                            }

                            multiplicity_by_alignment.emplace_back(item_count);
                            for (size_t seed_i : seeds_for_extension[extension_i]) {
                                minimizer_explored.insert(seeds.at(seed_i).source);
                            }

                            if (show_work) {
                                #pragma omp critical (cerr)
                                {
                                    cerr << log_name() << "Produced additional alignment "
                                         << "directly from full length gapless extension " << extension_i << endl;
                                }
                            }
                        }
                    }
                }
                // If we got at least one full-length extension as an alignment,
                // Then skip chaining for this tree
                if (alignments.size() >= 1) {
                    if (track_provenance) {
                        //We might have already done some chaining so the funnel might already have started on that stage
                        //So to get the funnel to track the gapless extensions properly, we need to make a fake chaining
                        //stage for these too
                        // Tell the funnel
                        //TODO: idk what score to give it funnel.score(funnel.latest(), scored_chain.first);!

                        funnel.project(item_num);

                        funnel.processed_input();

                        //Add an entry to the list of chains so we know which chain num to give the alignments
                        //This is just so the funnel can track everything
                        subchain_groups.emplace_back();

                    }
                    return true;
                }


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
                            seed_positions.push_back(minimizers[seeds.at(seed_index).source].pin_offset());
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
                        while (seed_it != extension_seeds.end() && minimizers[seeds.at(*seed_it).source].pin_offset() < anchor_interval.first) {
                            // Move seed iterator to inside or past the interval (should really always be already inside).
                            ++seed_it;
                        }
                        while (seed_it != extension_seeds.end() && minimizers[seeds.at(*seed_it).source].pin_offset() < anchor_interval.second) {
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
                            // seeds into the elaborating chain (only non-overlapping)
                            extension_seed_sequences.push_back({anchor_seeds.front()});
                            for (const auto& cur_seed : anchor_seeds) {
                                if (seed_anchors.at(extension_seed_sequences.back().back()).read_end() <= seed_anchors.at(cur_seed).read_start()) {
                                    // This seed doesn't overlap the previous one
                                    extension_seed_sequences.back().push_back(cur_seed);
                                }
                            }

                            // Keep all the seeds that this anchor counts as using.
                            extension_represented_seeds.emplace_back(std::move(anchor_seeds));
                        }
                    }
                }
            }

            // Figure out what anchors we want to view.
            const std::vector<algorithms::Anchor>& anchors_to_chain = do_gapless_extension ? extension_anchors : seed_anchors;
            // And what seeds each represents
            const std::vector<std::vector<size_t>>& anchor_seed_sequences = do_gapless_extension ? extension_seed_sequences : seed_seed_sequences;
            // And what subset/in what order
            std::vector<size_t>& anchor_indexes = do_gapless_extension ? extension_anchor_indexes : selected_seeds;
            // Sort anchors by read start of seeded region
            algorithms::sort_anchor_indexes(anchors_to_chain, anchor_indexes);

            // And what seeds should count as explored when we take an anchor
            const std::vector<std::vector<size_t>>& anchor_represented_seeds = do_gapless_extension ? extension_represented_seeds : anchor_seed_sequences;

            if (track_provenance) {
                funnel.substage("chain");
            }

            // Make a view of the anchors we will chain over
            VectorView<algorithms::Anchor> anchor_view {anchors_to_chain, anchor_indexes};

            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Computing chains over " << anchor_view.size() << " anchors" << endl;
                }
            }
#ifdef debug
            if (show_work) {
                // Log the chaining problem so we can try it again elsewhere.
                this->dump_chaining_problem(anchors_to_chain, anchor_indexes, gbwt_graph);
            }
#endif

            // Compute lookback and indel limits based on read length.
            // Important since seed density goes down on longer reads.
            size_t graph_lookback_limit = std::max(this->max_graph_lookback_bases, (size_t)(this->max_graph_lookback_bases_per_base * aln.sequence().size()));
            size_t read_lookback_limit = std::max(this->max_read_lookback_bases, (size_t)(this->max_read_lookback_bases_per_base * aln.sequence().size()));
            size_t indel_limit = std::max(this->max_indel_bases, (size_t)(this->max_indel_bases_per_base * aln.sequence().size()));

            // Find chains over the seeds in the zip code tree
            algorithms::transition_iterator for_each_transition = algorithms::zip_tree_transition_iterator(
                seeds,
                zip_code_forest.trees[item_num],
                graph_lookback_limit,
                read_lookback_limit
            );
            // TODO: Should we just inherit from ChainScoringScheme? Or should
            // we set one up as a member?
            algorithms::ChainScoringScheme scheme {
                this->item_bonus,
                this->gap_scale,
                this->rec_penalty,
                // TODO: Do this once at setup?
                this->rec_consistency_bonus == -1 ? this->rec_penalty : this->rec_consistency_bonus,
            };
            vector<algorithms::SubchainGroup> new_groups = algorithms::find_best_chains(
                anchor_view,
                *distance_index,
                gbwt_graph,
                aln.sequence().size(),
                scheme,
                this->max_chains_per_tree,
                for_each_transition,
                indel_limit,
                this->max_alt_lookback_score,
                this->extra_tail_grace_window,
                show_work);

            for (auto& group : new_groups) {
                if (show_work) {
                    #pragma omp critical (cerr)
                    cerr << log_name() << "\t[" << aln.name() << "] Found "
                         << group.subchains.size() << " subchains with "
                         << group.connections.size() << " inter-subchain connections in zip code tree " << item_num
                         << " running " << anchors_to_chain[anchor_indexes.front()] << " to " << anchors_to_chain[anchor_indexes.back()] << std::endl;
                }

                for (size_t subchain_i = 0; subchain_i < group.subchains.size(); subchain_i++) {
                    // For each subchain
                    vector<size_t>& subchain = group.subchains[subchain_i].anchors;
    #ifdef debug_rec
                    if (true)
    #else
                    if (show_work)
    #endif
                    {
    #ifdef debug
                        if(true)
    #else
                        if (subchain_i < MANY_LIMIT)
    #endif
                        {
                            #pragma omp critical (cerr)
                            {
                                cerr << log_name() << "\t[" << aln.name() << "] Subchain " << subchain_i
                                     << " and length " << subchain.size()
                                     << " with " << group.subchains[subchain_i].rec_count << " recombinations "
                                     << " running " << anchor_view[subchain.front()]
                                     << " to " << anchor_view[subchain.back()] << std::endl;
    #ifdef debug_rec
                                algorithms::path_flags_t current_paths = 0;
                                bool first = true;
                                for (auto& selected_number : subchain) {
                                    auto& anchor = anchor_view[selected_number];
                                    auto new_paths = anchor.anchor_paths();
                                    if (first) {
                                        current_paths = new_paths.second;
                                        first = false;
                                    } else {
                                        if (new_paths.first == new_paths.second) {
                                            if ((current_paths & new_paths.first) == 0) {
                                                current_paths = new_paths.first;
                                            } else {
                                                current_paths &= new_paths.first;
                                            }
                                        } else {
                                            current_paths = new_paths.second;
                                        }
                                    }
                                    
                                    std::cerr << log_name() << "\t\t" << anchor 
                                              << " anchor_paths: " << std::bitset<64>(new_paths.first).count() << " " << std::bitset<64>(new_paths.first) 
                                              << " chain_paths: " << std::bitset<64>(current_paths).count() << " " << std::bitset<64>(current_paths) << std::endl;
                                }
    #endif
                            }
                        } else if (subchain_i == MANY_LIMIT) {
                            #pragma omp critical (cerr)
                            std::cerr << log_name() << "\t[" << aln.name() << "] <" << (group.subchains.size() - subchain_i) << " more chains>" << std::endl;
                        }
                    }

                    // Count how many of each minimizer is in each chain produced
                    minimizer_kept_chain_count.emplace_back(minimizers.size(), 0);

                    // Translate subchains into seed numbers and not local anchor numbers.
                    vector<size_t> seed_nums;
                    seed_nums.reserve(subchain.size() * 2);

                    for (auto& selected_number : subchain) {
                        // For each anchor in the chain, get its number in the whole group of anchors.
                        size_t anchor_number = anchor_indexes.at(selected_number);
                        for (auto& seed_number : anchor_seed_sequences.at(anchor_number)) {
                            // And get all the seeds it actually uses in sequence and put them in the chain.
                            seed_nums.push_back(seed_number);
                        }
                        for (auto& seed_number : anchor_represented_seeds.at(anchor_number)) {
                            // And get all the seeds it represents exploring and mark their minimizers explored.
                            // TODO: Can we get the gapless extension logic to count this for us for that codepath?
                            minimizer_kept_chain_count.back()[seeds[seed_number].source]++;
                        }
                    }
                    subchain = seed_nums;

                    // Remember how we got it
                    subchain_source_tree.push_back(item_num);
                    // Remember the number of better or equal-scoring things
                    multiplicity_by_chain.emplace_back((float)item_count);

                    if (track_provenance) {
                        // Tell the funnel
                        funnel.introduce();
                        /// TODO: no score provided because these are intentionally just pieces
                        // We come from all the seeds directly
                        // TODO: Include all the middle seeds when gapless extending!
                        funnel.also_merge_group(2, subchain.begin(), subchain.end());
                        // And are related to the problem
                        funnel.also_relevant(1, item_num);
                    }

                    if (track_position && subchain_i < MANY_LIMIT) {
                        // Add position annotations for some chains.
                        // Should be much faster than full correctness tracking from every seed.
                        crash_unless(this->path_graph);
                        std::unordered_set<PathSense> wanted_senses {PathSense::REFERENCE, PathSense::GENERIC};
                        if (haplotype_positions) {
                            wanted_senses.insert(PathSense::HAPLOTYPE);
                        }
                        for (auto& boundary : {anchor_view[subchain.front()].graph_start(), anchor_view[subchain.back()].graph_end()}) {
                            // For each end of the chain
                            auto offsets = algorithms::nearest_offsets_in_paths(this->path_graph, boundary, 100, wanted_senses);
                            for (auto& handle_and_positions : offsets) {
                                for (auto& position : handle_and_positions.second) {
                                    // Tell the funnel all the effective positions, ignoring orientation
                                    funnel.position(funnel.latest(), handle_and_positions.first, position.first);
                                }
                            }

                        }
                    }

                    if (track_provenance && show_work && subchain_i < MANY_LIMIT) {
                        for (auto& handle_and_range : funnel.get_positions(funnel.latest())) {
                            // Log each range on a path associated with the chain.
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

                // Save our work
                subchain_groups.push_back(group);
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
                funnel.fail("max-chaining-problems", item_num);
            }

        }, [&](size_t item_num) -> void {
            // This item is not sufficiently good.
            if (track_provenance) {
                funnel.fail("zipcode-tree-coverage-threshold", item_num, tree_coverages[item_num]);
            }
        });

    if (alignments.size() >= 1) {
        //If we did get alignments from chaining, boot them through the funnel all at once
        funnel.stage("extension_to_alignment");
        for (size_t chain_num : alignment_source_chain) {
            funnel.project(chain_num);
        }
        //Get the actual multiplicity from the counts
        for (size_t i = 0 ; i < multiplicity_by_alignment.size() ; i++) {
            multiplicity_by_alignment[i] = multiplicity_by_alignment[i] >= kept_tree_count
                                        ?  multiplicity_by_alignment[i] - (float)kept_tree_count
                                        : 0.0;
        }

    } else {

        //Get the actual multiplicity from the counts
        for (size_t i = 0 ; i < multiplicity_by_chain.size() ; i++) {
            multiplicity_by_chain[i] = multiplicity_by_chain[i] >= kept_tree_count
                                        ?  multiplicity_by_chain[i] - (float)kept_tree_count
                                        : 0.0;
        }
    }

}

void MinimizerMapper::do_alignment_on_chains(const Alignment& aln, const std::vector<Seed>& seeds, 
                                             const VectorView<MinimizerMapper::Minimizer>& minimizers,
                                             const vector<algorithms::Anchor>& seed_anchors,
                                             const std::vector<algorithms::SubchainGroup>& subchain_groups, 
                                             const std::vector<size_t>& subchain_source_tree,
                                             const std::vector<double>& multiplicity_by_chain,
                                             const std::vector<std::vector<size_t>>& minimizer_kept_chain_count,
                                             vector<Alignment>& alignments, vector<double>& multiplicity_by_alignment,
                                             vector<size_t>& alignments_to_source,
                                             SmallBitset& minimizer_explored, aligner_stats_t& stats,
                                             LazyRNG& rng, Funnel& funnel) const {
  
    if (track_provenance) {
        funnel.stage("align");
    }
    //For finding the multiplicity of each alignment, first get the count
    // of equal scoring chains
    vector<size_t> chain_count_by_alignment (alignments.size(), 0);

#ifdef print_minimizer_table
    //How many of each minimizer ends up in a chain that actually gets turned into an alignment?
    vector<size_t> minimizer_kept_count(minimizers.size(), 0);
#endif

    // Compute lower limit on chain score to actually investigate
    int chain_min_score = (int) (min_chain_score_per_base * aln.sequence().size());
    // Apply the max min chain score limit
    chain_min_score = std::min(chain_min_score, max_min_chain_score);
    vector<int> max_sparse_chain_scores;
    for (const auto& group : subchain_groups) {
        max_sparse_chain_scores.emplace_back(group.max_sparse_chain_score);
    }

    // Remember: we also have chain_score_threshold, which counts down from best chain score
    
    // We need to be able to discard a SubchainGroup because its score isn't good enough.
    // We have more components to the score filter than process_until_threshold_b supports.
    auto discard_chain_by_score = [&](size_t processed_num) -> void {
        // This SubchainGroup is not good enough.
        if (track_provenance) {
            funnel.fail("min-chain-score-per-base||max-min-chain-score", processed_num, max_sparse_chain_scores[processed_num]);
        }
        
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "subchain group " << processed_num
                     << " failed because its score was not good enough (max score="
                     << max_sparse_chain_scores[processed_num]
                     << ", min=" << chain_min_score
                     << ", threshold " << chain_score_threshold << " off best)" << endl;
                if (track_correctness && funnel.was_correct(processed_num)) {
                    cerr << log_name() << "\tCORRECT!" << endl;
                }
            }
        }
    };
    
    // Track how many tree chains were used
    std::unordered_map<size_t, size_t> chains_per_tree;
    // Track total count of alignments made (we will make at most max_alignments)
    size_t alns_made = 0;

    // Track what node ID, orientation, read-minus-node offset tuples were used
    // in previously generated alignments, so we can fish out alignments to
    // different placements.
    // Use pairs since we can't hash tuples.
    std::unordered_set<std::pair<std::pair<nid_t, bool>, int64_t>> used_matchings;

    
    // Go through the chains in estimated-score order.
    process_until_threshold_b<int>(max_sparse_chain_scores,
        chain_score_threshold, min_chains, max_alignments, rng, 
        [&](size_t processed_num, size_t item_count) -> bool {
            // This subchain group is good enough.
            // Called in descending score order.

            if (alns_made >= max_alignments) {
                // Earlier SubchainGroups made enough alignments already
                return false;
            }
        
            if (max_sparse_chain_scores[processed_num] < chain_min_score && alns_made >= min_chains) {
                // Actually discard by score
                discard_chain_by_score(processed_num);
                return false;
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "subchain group " << processed_num
                         << " is good enough (max score=" << max_sparse_chain_scores[processed_num]
                         << ", min=" << chain_min_score
                         << ", threshold " << chain_score_threshold << " off best)" << endl;
                    if (track_correctness && funnel.was_correct(processed_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                    }
                }
            }
            if (track_provenance) {
                funnel.pass("min-chain-score-per-base||max-min-chain-score", processed_num, max_sparse_chain_scores[processed_num]);
                funnel.pass("max-alignments", processed_num);
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
                
                try {
                    // Do the DP between the items in the chain

                    // Collect stats into here
                    aligner_stats_t alignment_stats;
                    best_alignments = do_base_level_alignment(aln, seed_anchors, subchain_groups.at(processed_num), max_alignments, funnel, &alignment_stats);
                    /// TODO: need to rethink how this works now that I return multiple alignments
                    //alignment_stats.add_annotations(best_alignments[0], "alignment");

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
                    pos_t graph_pos = make_pos_t(mapping.position());

                    nid_t node_id = id(graph_pos);
                    bool orientation = is_rev(graph_pos);
                    size_t graph_offset = offset(graph_pos);

                    for (auto& edit : mapping.edit()) {
                        if (edit.sequence().empty() && edit.from_length() == edit.to_length()) {
                            // It's an actual match so make a matching
                            int64_t read_minus_node_offset = (int64_t)read_pos - (int64_t)graph_offset;
                            auto matching = std::make_pair(std::make_pair(node_id, orientation), read_minus_node_offset);

#ifdef debug
                            if (show_work) {
                                #pragma omp critical (cerr)
                                {
                                    cerr << log_name() << "Create matching " << matching.first.first << ", " << matching.first.second << ", " << matching.second << endl;
                                }
                            }
#endif

                            used_matchings.emplace(std::move(matching));
                        }
                        read_pos += edit.to_length();
                        graph_offset += edit.from_length();
                    }
                    
                }

                if (track_provenance) {
                    funnel.project(processed_num);
                    funnel.score(alignments.size() - 1, alignments.back().score());
                }
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Produced alignment from chain group " << processed_num
                             << " with score " << alignments.back().score() << ": " << log_alignment(alignments.back()) << endl;
                    }
                }
            };
            
            for(auto aln_it = best_alignments.begin() ; 
                aln_it != best_alignments.end() && aln_it->score() != 0 
                    && (aln_it->score() >= best_alignments[0].score() * 0.8 || aln_it->score() == best_alignments[0].score()) ;
                ++aln_it) {
                //For each additional alignment with score at least 0.8 of the best score
                //Guarantee that all alignments with top score (even if negative) are used
                observe_alignment(*aln_it);
                alns_made++;
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
                funnel.pass("min-chain-score-per-base||max-min-chain-score", processed_num, max_sparse_chain_scores[processed_num]);
                funnel.fail("max-alignments", processed_num);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "subchain group " << processed_num 
                         << " failed because there were too many good groups (max score="
                         << max_sparse_chain_scores[processed_num] << ")" << endl;
                    if (track_correctness && funnel.was_correct(processed_num)) {
                        cerr << log_name() << "\tCORRECT!" << endl;
                    }
                }
            }
        }, discard_chain_by_score);
    
    for (size_t i = 0 ; i < multiplicity_by_alignment.size() ; ++i) {
        multiplicity_by_alignment[i] += (chain_count_by_alignment[i] >= alignments.size()
                                      ? ((double)chain_count_by_alignment[i] - (double) alignments.size())
                                      : 0.0);
    }
}

void MinimizerMapper::pick_mappings_from_alignments(const Alignment& aln, const std::vector<Alignment>& alignments, 
                                                    const std::vector<double>& multiplicity_by_alignment,
                                                    const std::vector<size_t>& alignments_to_source,
                                                    std::vector<Alignment>& mappings,
                                                    std::vector<double>& scores,
                                                    std::vector<double>& multiplicity_by_mapping,
                                                    LazyRNG& rng,
                                                    Funnel& funnel) const {
    // Grab all the scores in order for MAPQ computation.
    scores.reserve(alignments.size());
    
    // Go through the alignments in descending score order, with ties at the top end shuffled.
    process_until_threshold_a(alignments.size(), (std::function<double(size_t)>) [&](size_t i) -> double {
        // Tiebreak by identity (which is always 0 to 1)
        return alignments.at(i).score() + identity(alignments.at(i).path());
    }, 0, 1, max_multimaps, rng, [&](size_t alignment_num, size_t item_count) {
        // This alignment makes it
        // Called in score order
        
        // Filter to alignments with strictly positive scores
        if (alignments[alignment_num].score() <= 0) {
            if (track_provenance) {
                funnel.fail("nonzero-score", alignment_num);
            }
            return false;
        } else {
            if (track_provenance) {
                funnel.pass("nonzero-score", alignment_num);
            }
        }

        if (track_provenance) {
            // Tell the funnel
            funnel.pass("max-multimaps", alignment_num);
        }

        // Remember the score at its rank
        scores.emplace_back(alignments[alignment_num].score());
        
        // Remember the output alignment
        mappings.emplace_back(std::move(alignments[alignment_num]));

        // Remember the multiplicity
        multiplicity_by_mapping.emplace_back(multiplicity_by_alignment[alignment_num]);
        
        if (track_provenance) {
            // Tell the funnel
            funnel.project(alignment_num);
            funnel.score(funnel.latest(), scores.back());
        }
        
        return true;
    }, [&](size_t alignment_num) {
        // We already have enough alignments, although this one has a good score

        // TODO: We end up having to duplicate a bunch of filters here so the
        // filters are always in order.
        
        // Go back and do the nonzero score filter first.
        // Filter to alignments with strictly positive scores
        if (alignments[alignment_num].score() <= 0) {
            if (track_provenance) {
                funnel.fail("nonzero-score", alignment_num);
            }
            // If we fail the nonzero score filter, we won't count as a secondary for MAPQ
            return;
        } else {
            if (track_provenance) {
                funnel.pass("nonzero-score", alignment_num);
            }
        }

        // Remember the score at its rank even if it won't be output as a multimapping
        scores.emplace_back(alignments[alignment_num].score());
        multiplicity_by_mapping.emplace_back(multiplicity_by_alignment[alignment_num]);
        
        if (track_provenance) {
            funnel.fail("max-multimaps", alignment_num);
        }
    }, [&](size_t alignment_num) {
        // This alignment does not have a sufficiently good score.
        // It may have been penalized into negative score.
        
        if (track_provenance) {
            // Call this a fail of the nonzero score filter, even though that
            // really filters on alignment score and not sorting score, which
            // can be different.
            //
            // TODO: Remove feature to sort by chain score?
            funnel.fail("nonzero-score", alignment_num);
        }
    });

    if (mappings.empty()) {
        // We didn't find any mappings, so make an unmapped one.

        scores.emplace_back(0);
        mappings.emplace_back(aln);
        multiplicity_by_mapping.emplace_back(0);
        
        if (track_provenance) {
            // Tell the funnel
            funnel.introduce();
            funnel.score(funnel.latest(), scores.back());
        }
    }
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

MinimizerMapper::ScoredPath MinimizerMapper::find_tail_alignment(
    const Alignment& aln, 
    const algorithms::Anchor& tail_anchor, 
    const WFAExtender& wfa_extender, 
    bool is_left_tail,
    aligner_stats_t* stats
) const {
    // Set up alignment parameters
    string tail_side = is_left_tail ? "left" : "right";
    size_t tail_length = is_left_tail ? tail_anchor.read_start()
                                      : aln.sequence().size() - tail_anchor.read_end();
    string tail_seq = is_left_tail ? aln.sequence().substr(0, tail_length)
                                   : aln.sequence().substr(tail_anchor.read_end(), tail_length);
    pos_t anchor_pos = is_left_tail ? tail_anchor.graph_start()
                                    : tail_anchor.graph_end();
    if (!is_left_tail) {
        // Pull back a base to get the outside-the-alignment anchoring position.
        get_offset(anchor_pos)--;
    }

    ScoredPath output;
    // Handle empty input
    if (tail_length == 0) {
        return output;
    }

    // Scratch to calculate time differences
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point stop_time;

    WFAAlignment tail_wfa_aln;

    // ---- Initial GBWT extension ----

    if (tail_length <= this->max_tail_length) {
        if (stats) {
            start_time = std::chrono::high_resolution_clock::now();
        }

        // Generate a prefix/suffix alignment for the tail
        tail_wfa_aln = is_left_tail ? wfa_extender.prefix(tail_seq, anchor_pos)
                                    : wfa_extender.suffix(tail_seq, anchor_pos);

        if (stats) {
            stop_time = std::chrono::high_resolution_clock::now();
            stats->bases.wfa_tail += tail_length;
            stats->time.wfa_tail += std::chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count();
            stats->invocations.wfa_tail += 1;
        }

        if (tail_wfa_aln) {
            if (!is_left_tail) {
                // Shift the alignment over to where the tail actually is
                tail_wfa_aln.seq_offset += tail_anchor.read_end();
            }
            // Did we get all the way to the end of the read?
            // If not, add a softclip.
            // TODO: Can we let the aligner know it can softclip for free?
            size_t softclipped_bases = 0;
            if (is_left_tail && tail_wfa_aln.seq_offset != 0) {
                // Prepend softclip
                softclipped_bases = tail_wfa_aln.seq_offset;
                WFAAlignment prepend = WFAAlignment::make_unlocalized_insertion(0, softclipped_bases, 0);
                prepend.join(tail_wfa_aln);
                tail_wfa_aln = std::move(prepend);
            } else if (!is_left_tail && tail_wfa_aln.seq_offset + tail_wfa_aln.length != aln.sequence().size()) {
                // Append softclip
                size_t right_end = tail_wfa_aln.seq_offset + tail_wfa_aln.length;
                softclipped_bases = aln.sequence().size() - right_end;
                tail_wfa_aln.join(WFAAlignment::make_unlocalized_insertion(right_end, softclipped_bases, 0));
            }

            if (this->softclip_penalty != 0.0 && softclipped_bases > 0) {
                double penalty = this->softclip_penalty * softclipped_bases;
                if (show_work) {
                    #pragma omp critical (cerr)
                    cerr << log_name() << "Applied softclip penalty of " << penalty 
                         << " for " << softclipped_bases << " total softclipped bases" << endl;
                }
                tail_wfa_aln.score -= penalty;
            }
        }

        if (tail_wfa_aln.length != tail_length) {
            // We didn't get the alignment we expected.
            stringstream ss;
            ss << "Aligning " << tail_side << " tail " << tail_seq 
               << " anchored at " << anchor_pos << " produced wrong-length alignment ";
            tail_wfa_aln.print(ss);
            throw ChainAlignmentFailedError(ss.str());
        }
    }

    // ---- Check for early return ----

    // Do we already have an alignment (from GBWT step?)
    if (tail_wfa_aln) {
        tail_wfa_aln.check_lengths(this->gbwt_graph);
            
#ifdef debug_base_level_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Add " << tail_side << " tail of "
                     << tail_wfa_aln.length << "bp with score of " << tail_wfa_aln.score << endl;
            }
        }
#endif
        output.path = tail_wfa_aln.to_path(this->gbwt_graph, aln.sequence());
        output.score = tail_wfa_aln.score;
        return output;
    }

    // Is the tail too long to align?
    if (tail_length > this->max_tail_dp_length) {
#ifdef debug_base_level_alignment
        #pragma omp critical (cerr)
        {
            cerr << "warning[MinimizerMapper::find_tail_alignment]: Refusing to align "
                 << tail_length << " bp " << tail_side << " tail against "
                 << anchor_pos << " in " << aln.name() << endl;
        }
#endif
                
        // Make a softclip for it.
        tail_wfa_aln = WFAAlignment::make_unlocalized_insertion(is_left_tail ? 0 : tail_anchor.read_end(), tail_length, 0);
        output.path = tail_wfa_aln.to_path(this->gbwt_graph, aln.sequence());
        output.score = tail_wfa_aln.score;
        return output;
    }

    // ---- Fall back on alignment against graph ----
            
#ifdef debug_base_level_alignment
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Long " << tail_side << " tail fallback alignment" << endl;
        }
    }
#endif
                
    Alignment tail_fallback_aln;
    tail_fallback_aln.set_sequence(tail_seq);
    if (!aln.quality().empty()) {
        tail_fallback_aln.set_quality(aln.quality().substr(0, tail_length));
    }
                
    // Work out how far the tail can see
    auto seq_start = is_left_tail ? aln.sequence().begin()
                                  : aln.sequence().begin() + tail_anchor.read_end();
    auto seq_end = is_left_tail ? aln.sequence().begin() + tail_length
                                : aln.sequence().end();
    size_t max_gap_length = longest_detectable_gap_in_range(aln, seq_start, seq_end, this->get_regular_aligner());
    max_gap_length = std::min(this->max_tail_gap, max_gap_length);
    size_t graph_horizon = tail_length + max_gap_length;

    // Align the tail, anchoring one end.
    if (stats) {
        start_time = std::chrono::high_resolution_clock::now();
    }
    bool did_aln = align_sequence_between_consistently(is_left_tail ? empty_pos_t() : tail_anchor.graph_end(),
                                                       is_left_tail ? tail_anchor.graph_start() : empty_pos_t(),
                                                       graph_horizon,
                                                       max_gap_length,
                                                       &this->gbwt_graph,
                                                       this->get_regular_aligner(),
                                                       tail_fallback_aln,
                                                       &aln.name(),
                                                       this->max_dp_cells,
                                                       this->choose_band_padding);
    if (stats) {
        stop_time = std::chrono::high_resolution_clock::now();
        if (did_aln) {
            // Actually did the alignment
            stats->bases.dozeu_tail += tail_length;
            stats->time.dozeu_tail += std::chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count();
            stats->invocations.dozeu_tail += 1;
        }
    }
                
    if (show_work && max_tail_length > 0) {
        #pragma omp critical (cerr)
        {
            cerr << "warning[MinimizerMapper::find_tail_alignment]: Fallback score: " << tail_fallback_aln.score() << endl;
        }
    }

    output.path = tail_fallback_aln.path();
    output.score = tail_fallback_aln.score();
        
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Aligned " << tail_side << " tail length " << tail_length << std::endl;
        }
    }

    return output;
}

void MinimizerMapper::find_next_non_overlapping(
    const VectorView<algorithms::Anchor>& to_chain,
    const std::vector<size_t>& chain, 
    const algorithms::Anchor*& here, 
    vector<size_t>::const_iterator& next_it) const {
    // Don't move past the end of the chain
    while (next_it != chain.end()) {
        // Try and find a next thing to connect to
        if (algorithms::get_read_distance(*here, *(&to_chain[*next_it])) == std::numeric_limits<size_t>::max()) {
            // There's overlap between these items. Keep here and skip next.
#ifdef debug_base_level_alignment
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Don't try and connect to " << *next_it << " because they overlap" << endl;
                }
            }
#endif
        
            ++next_it;
        } else {
            // No overlap, so try it.
            break;
        }
    }
}

void MinimizerMapper::find_next_to_skip_to(
    const VectorView<algorithms::Anchor>& to_chain,
    const std::vector<size_t>& chain,
    const algorithms::Anchor*& here,
    vector<size_t>::const_iterator& next_it) const {
    // Keep track of the total distance from the previous seed to the next one we choose in the graph
    const algorithms::Anchor* next = &to_chain[*next_it];
    size_t total_graph_distance = algorithms::get_graph_distance(*here, *next, *this->distance_index, this->gbwt_graph);
    size_t prev_read_distance = algorithms::get_read_distance(*here, *next);

    // The sum of the differences between read and graph lengths
    size_t gap_lengths = (std::max(total_graph_distance, prev_read_distance) 
                        - std::min(total_graph_distance, prev_read_distance));
    // Where we will be moving next_it to (maybe)
    auto skip_to_it = next_it;

    while (skip_to_it != chain.end()) {
        const algorithms::Anchor* skip_to = &to_chain[*skip_to_it];
        // Try and find a next thing to connect to
        
        //TODO: Getting the graph distance is probably slow, might want to save it from chaining
        size_t cur_graph_distance;
        if (skip_to->is_skippable() && skip_to_it+1 != chain.end()) {
            // Distance from end of this anchor to start of next anchor
            cur_graph_distance = algorithms::get_graph_distance(*skip_to, to_chain[*(skip_to_it+1)], 
                                                                *this->distance_index, this->gbwt_graph);
            // Also add in distance from start of this anchor to end of this anchor
            cur_graph_distance += skip_to->length();
            // Combined those are graph start -> start dist we are trying to skip past
        }

        if (skip_to->is_skippable() && skip_to_it+1 != chain.end() && 
            total_graph_distance + cur_graph_distance < this->max_skipped_bases) {
            // This anchor is repetitive and the next one is close enough to connect
#ifdef debug_base_level_alignment
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Try to avoid connecting to " << *skip_to_it << " because it is repetitive" << endl;
                }
            }
#endif
            // Read start -> start distance for what we're skipping
            size_t cur_read_distance = algorithms::get_read_distance(*skip_to, to_chain[*(skip_to_it+1)]);
            cur_read_distance += skip_to->length();
            // Total gap so far
            gap_lengths += (std::max(cur_read_distance, cur_graph_distance) 
                            - std::min(cur_read_distance, cur_graph_distance));
            total_graph_distance += cur_graph_distance;
            
            ++skip_to_it;
        } else {
            // skip_to is either not skippable or too far away so stop
            if (gap_lengths > this->min_indel_avoid_bases) {
#ifdef debug_base_level_alignment
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Skipping to " << *skip_to_it 
                             << " to avoid gap of " << gap_lengths << " in repetitive region" << endl;
                    }
                }
#endif
                // If there was a big gap
                next_it = skip_to_it;
            } else {
#ifdef debug_base_level_alignment
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "Not bothering to skip to " << *skip_to_it 
                             << " because total gaps are only " << gap_lengths << endl;
                    }
                }
#endif
            }
            // If there wasn't a gap then don't skip anything
            break;
        }
    }
}

MinimizerMapper::ScoredPath MinimizerMapper::find_link_alignment(
    const VectorView<algorithms::Anchor>& to_chain,
    const Alignment& aln,
    const vector<size_t>::const_iterator& here_it,
    const vector<size_t>::const_iterator& next_it,
    const WFAExtender& wfa_extender,
    const Aligner& aligner,
    aligner_stats_t* stats) const {

    // Where are we?
    const algorithms::Anchor* here = &to_chain[*here_it];
    const algorithms::Anchor* next = &to_chain[*next_it];

#ifdef debug_base_level_alignment
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Next connectable item " << *next
                 << " with overall index " << to_chain.backing_index(*next_it)
                 << " aligns " << (*next).read_start() << "-" << (*next).read_end()
                 << " with " << (*next).graph_start() << "-" << (*next).graph_end()
                 << endl;
        }
    }
#endif
    // Where we save output
    ScoredPath output;
    // Pull out the intervening string to the next, if any.
    size_t link_start = (*here).read_end();
    size_t link_length = (*next).read_start() - link_start;
    string linking_bases = aln.sequence().substr(link_start, link_length);
    size_t graph_dist = algorithms::get_graph_distance(*here, *next, *this->distance_index, this->gbwt_graph);
    // Storage for this alignment once we find it
    WFAAlignment link_alignment;
    // Where did it come from?
    std::string link_alignment_source;

    // Scratch to calculate time differences
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point stop_time;
        
#ifdef debug_base_level_alignment
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Need to align graph from " << (*here).graph_end() << " to " << (*next).graph_start()
                 << " separated by ~" << graph_dist << " bp";
            if (linking_bases.size() < 200) {
                cerr << " and sequence \"" << linking_bases << "\"";
            }
            cerr << endl;
        }
    }
#endif
        
    if (link_length == 0 && graph_dist == 0) {
        // These items abut in the read and the graph, so we assume we can just connect them.
        // WFAExtender::connect() can't handle an empty read sequence, and
        // our fallback method to align just against the graph can't handle
        // an empty graph region.
        // TODO: We can be leaving the GBWT's space here!
        
#ifdef debug_base_level_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Treat as empty link" << endl;
            }
        }
#endif
        
        link_alignment = WFAAlignment::make_empty();
        link_alignment_source = "empty";
    } else if (link_length > 0 && link_length <= this->max_chain_connection && graph_dist * link_length <= max_dp_cells) {
        // If it's not empty and is a reasonable size, align it.
        // Make sure to walk back the left anchor so it is outside of the region to be aligned.
        pos_t left_anchor = (*here).graph_end();
        get_offset(left_anchor)--;
        
        if (stats) {
            start_time = std::chrono::high_resolution_clock::now();
        }
        link_alignment = connect_consistently(linking_bases, left_anchor, (*next).graph_start(), wfa_extender);
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
        
        if (!link_alignment) {
            // We couldn't align.
            if (graph_dist == 0) {
                // We had read sequence but no graph sequence.
                // Try falling back to a pure insertion.
                // TODO: We can be leaving the GBWT's space here!
                // TODO: What if this is forcing an insertion that could also be in the graph already?
#ifdef debug_base_level_alignment
                if (show_work) {
                    #pragma omp critical (cerr)
                    {
                        cerr << log_name() << "connect() failed; treat as insertion" << endl;
                    }
                }
#endif
                link_alignment = WFAAlignment::make_unlocalized_insertion((*here).read_end(), link_length, aligner.scorer->score_gap(link_length));
                link_alignment_source = "unlocalized_insertion";
            }
        } else if (link_alignment.length != linking_bases.size()) {
            // We could align, but we didn't get the alignment we expected. This shouldn't happen for a middle piece that can't softclip.
            stringstream ss;
            ss << "Aligning anchored link " << linking_bases << " (" << linking_bases.size()
               << " bp) from " << left_anchor << " - " << (*next).graph_start()
               << " against graph distance " << graph_dist << " produced wrong-length alignment ";
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
        
#ifdef debug_base_level_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Add link of length " << link_alignment.length
                     << " with score of " << link_alignment.score << endl;
            }
        }
#endif
    
        link_alignment.check_lengths(gbwt_graph);
        
        // Then the link (possibly empty)
        output.path = link_alignment.to_path(this->gbwt_graph, aln.sequence());
        output.score = link_alignment.score;
#ifdef debug_base_level_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "\t" << pb2json(output.path) << endl;
            }
        }
#endif
    } else {
        // The sequence to the next thing is too long, or we couldn't reach it doing connect().
        // Fall back to another alignment method
        
        if (linking_bases.size() > max_middle_dp_length || graph_dist * link_length > max_dp_cells) {
            // This would be too long for the middle aligner(s) to handle and might overflow a score somewhere.
            #pragma omp critical (cerr)
            {
                cerr << "warning[MinimizerMapper::find_link_alignment]: Refusing to align "
                     << link_length << " bp connection between chain items " 
                     << to_chain.backing_index(*here_it) << " and " << to_chain.backing_index(*next_it) 
                     << " which are " << graph_dist << " apart at " 
                     << (*here).graph_end() << " and " << (*next).graph_start() 
                     << " in " << aln.name() << " due to DP size limits, creating " 
                     << (aln.sequence().size() - (*here).read_end()) << " bp right tail" << endl;
            }
            // Give up
            output.score = -std::numeric_limits<int32_t>::max();
            return output;
        }
        
        Alignment link_aln;
        link_aln.set_sequence(linking_bases);
        if (!aln.quality().empty()) {
            link_aln.set_quality(aln.quality().substr(link_start, link_length));
        }
        // Guess how long of a graph path we ought to allow in the alignment.
        size_t max_gap_length = longest_detectable_gap_in_range(aln, aln.sequence().begin() + link_start, aln.sequence().begin() + link_start + link_length, this->get_regular_aligner());
        max_gap_length = std::min(this->max_middle_gap, max_gap_length);
        size_t path_length = std::max(graph_dist, link_length);
        if (stats) {
            start_time = std::chrono::high_resolution_clock::now();
        }
        bool did_aln = MinimizerMapper::align_sequence_between_consistently(
            (*here).graph_end(), 
            (*next).graph_start(), 
            path_length + max_gap_length,
            max_gap_length,
            &this->gbwt_graph,
            this->get_regular_aligner(),
            link_aln, &aln.name(),
            this->max_dp_cells,
            this->choose_band_padding);
        if (stats) {
            stop_time = std::chrono::high_resolution_clock::now();
            if (did_aln) {
                // Actually did the alignment
                stats->bases.bga_middle += link_length;
                stats->time.bga_middle += std::chrono::duration_cast<chrono::duration<double>>(stop_time - start_time).count();
                stats->invocations.bga_middle += 1;
            }
        }
        
        if (linking_bases.size() > 0 && link_aln.path().mapping_size() == 0) {
            // Connecting alignment bailed out. Assume that this is due to size.
            // TODO: Should we let the exceptions propagate up to here instead?
            #pragma omp critical (cerr)
            {
                cerr << "warning[MinimizerMapper::find_link_alignment]: BGA alignment too big for "
                     << link_length << " bp connection between chain items " << to_chain.backing_index(*here_it) 
                     << " and " << to_chain.backing_index(*next_it) 
                     << " which are " << graph_dist << " apart at " 
                     << (*here).graph_end() << " and " << (*next).graph_start() << " in " << aln.name() << endl;
            }

            // Give up
            output.score = -std::numeric_limits<int32_t>::max();
            return output;
        }
        
        // Otherwise we actually have a link alignment result.
        link_alignment_source = "align_sequence_between";
        
        // Then tack that path and score on
        output.path = link_aln.path();
        output.score = link_aln.score();
    }

#ifdef debug_base_level_alignment
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Aligned link from " << *here << " to " << *next << " of " 
                 << link_length << " bp read and " << path_from_length(output.path) << " bp graph via "
                 << link_alignment_source << " with score of " << output.score << std::endl;
        }
    }
#endif
    return output;
}

pair<MinimizerMapper::ScoredPath, size_t> MinimizerMapper::find_all_inner_chain_links(
    const VectorView<algorithms::Anchor>& to_chain,
    const Alignment& aln,
    const vector<size_t>& chain,
    const WFAExtender& wfa_extender,
    const Aligner& aligner,
    aligner_stats_t* stats) const {
    // Will be gradually built up by alignment
    ScoredPath output;
    // Keep a couple cursors in the chain: extension before and after the linking up we need to do.
    auto here_it = chain.begin();
    auto next_it = here_it;
    ++next_it;
    
    // Track the anchor we're at.
    // Note that, although it has a score, that's an anchor score; it isn't the
    // right score for the perfect-match alignment it represents.
    const algorithms::Anchor* here = &to_chain[*here_it];

#ifdef debug_base_level_alignment
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

    while (next_it != chain.end()) {
        // Do each region between successive gapless extensions
        
        // We have to find the next item we can actually connect to
        find_next_non_overlapping(to_chain, chain, here, next_it);
        // Next, we want to skip seeds that are in repetitive regions of the read
        // Since skipping all repetitive seeds would leave too many gaps in the chain,
        // only skip seeds if they are involved in gaps,
        // i.e. the distances in the read and graph are different
        find_next_to_skip_to(to_chain, chain, here, next_it);

        if (next_it == chain.end()) {
            // We couldn't find anything to connect to
            break;
        }

        // We have something to connect to! Make an alignment
        const algorithms::Anchor* next = &to_chain[*next_it];
            
#ifdef debug_base_level_alignment
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

#ifdef debug_base_level_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "\tScore " << here_alignment.score << endl;
            }
        }
#endif

        append_path(output.path, here_alignment.to_path(this->gbwt_graph, aln.sequence()));
        output.score += here_alignment.score;
        
#ifdef debug_base_level_alignment
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

        ScoredPath link_aln = find_link_alignment(to_chain, aln, here_it, next_it, wfa_extender, aligner, stats);

        if (link_aln.score == -std::numeric_limits<int32_t>::max()) {
            // We gave up. Jump to right tail.
            break;
        }

        append_path(output.path, std::move(link_aln.path));
        output.score += link_aln.score;
        
        // Advance here to next and start considering the next after it
        here_it = next_it;
        ++next_it;
        here = next;
    }

    if (next_it == chain.end()) {
        // We didn't bail out to treat a too-long connection as a tail. We still need to add the final extension anchor.
    
#ifdef debug_base_level_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Add last extension " << *here_it << " of length " << (*here).length() << endl;
            }
        }
#endif
    
        WFAAlignment here_alignment = this->to_wfa_alignment(*here, aln, &aligner);

#ifdef debug_base_level_alignment
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "\tScore " << here_alignment.score << endl;
            }
        }
#endif

        here_alignment.check_lengths(gbwt_graph);
    
        // Do the final GaplessExtension itself (may be the first)
        append_path(output.path, here_alignment.to_path(this->gbwt_graph, aln.sequence()));
        output.score += here_alignment.score;
    }

    return make_pair(output, *here_it);
}

/// TODO: split this function up it's getting really big
vector<Alignment> MinimizerMapper::do_base_level_alignment(
    const Alignment& aln,
    const VectorView<algorithms::Anchor>& to_chain,
    const algorithms::SubchainGroup& subchain_group,
    const size_t& max_alignments,
    Funnel& funnel,
    aligner_stats_t* stats
) const {
    
    if (subchain_group.subchains.empty()) {
        throw ChainAlignmentFailedError("Cannot find an alignment for an empty chain!");
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
    WFAExtender wfa_extender(gbwt_graph, aligner, wfa_error_model); 

    
    // Used to figure out which subchains should get tail alignments
    vector<bool> seen_as_source(to_chain.size(), false);
    vector<bool> seen_as_sink(to_chain.size(), false);
    for (const auto& extra_edge : subchain_group.connections) {
        seen_as_source[extra_edge.first] = true;
        seen_as_sink[extra_edge.second] = true;
    }

    multipath_alignment_t mp_aln;
    // We will store Paths for each node/edge in the MP alignment
    // The edges (connection_t) can't store Paths anyhow
    // and the nodes (subpath_t) require annoying conversion
    // We will piece the alignment back together with this memory
    size_t n_subchains = subchain_group.subchains.size();
    vector<Path> node_paths(n_subchains);
    unordered_map<pair<size_t, size_t>, Path> edge_paths;
        // We want to annotate alignments with their tail lengths
    // so for any subchains with a tail, save their length
    vector<double> left_tail_len(subchain_group.subchains.size());
    vector<double> right_tail_len(subchain_group.subchains.size());
    // Subchains where we bailed out of link alignments
    // We will have to ignore any connections which start from them
    unordered_set<size_t> early_bail_subchains;

    // Set up pseudo-subpaths
    for (size_t i = 0; i < subchain_group.subchains.size(); i++) {
        // Keep track of total base-level results for this subchain
        Path composed_path;
        int composed_score = 0;

        // Should this subchain get a left tail?
        if (!seen_as_sink[i]) {
#ifdef debug_base_level_alignment
            cerr << "Doing left tail alignment for subchain " << i << endl;
#endif
            ScoredPath left_tail = find_tail_alignment(
                aln, to_chain[subchain_group.subchains[i].anchors.front()], wfa_extender, true, stats);
            composed_path = left_tail.path;
            // Rescore if not just a softclip
            composed_score = left_tail.score;
            left_tail_len[i] = left_tail.path.length();
        }

        // Add in everything internal to the subchain
        ScoredPath inner_links;
        size_t last_anchor;
#ifdef debug_base_level_alignment
        cerr << "Doing inner link alignment for subchain " << i << endl;
#endif
        vg::tie(inner_links, last_anchor) = find_all_inner_chain_links(
            to_chain, aln, subchain_group.subchains[i].anchors, wfa_extender, aligner, stats);
        append_path(composed_path, inner_links.path);
        composed_score += inner_links.score;

        if (last_anchor != subchain_group.subchains[i].anchors.back()) {
#ifdef debug_base_level_alignment
            cerr << "Bailed out of subchain " << i << endl;
#endif
            // Oh no, we bailed out of a too-long chain connection
            early_bail_subchains.emplace(i);
        }

        // Should this subchain get a right tail?
        if (!seen_as_source[i] || last_anchor != subchain_group.subchains[i].anchors.back()) {
#ifdef debug_base_level_alignment
            cerr << "Doing right tail alignment for subchain " << i << endl;
#endif
            ScoredPath right_tail = find_tail_alignment(
                aln, to_chain[last_anchor], wfa_extender, false, stats);
            append_path(composed_path, right_tail.path);
            // Rescore if not just a softclip
            composed_score += right_tail.score;
            right_tail_len[i] = right_tail.path.length();
        }

        subpath_t* subpath = mp_aln.add_subpath();
        // Remember the path & score
        node_paths[i] = composed_path;
        subpath->set_score(composed_score);
        // Add a fake position to this subpath to store the subchain ID
        position_t* position = subpath->mutable_path()->add_mapping()->mutable_position();
        position->set_node_id(i);
        position->set_is_reverse(false);
        position->set_offset(0);
    }

    // Set up connections between subpaths
    for (const auto& extra_edge : subchain_group.connections) {
        // Only use edge if we didn't bail out of its source
        if (!early_bail_subchains.count(extra_edge.first)) {
#ifdef debug_base_level_alignment
            cerr << "Extra edge " << extra_edge.first << " -> " << extra_edge.second << endl;
#endif
            // Calculate base-level alignment for this connection
            vector<size_t> edge = {subchain_group.subchains[extra_edge.first].anchors.back(),
                                   subchain_group.subchains[extra_edge.second].anchors.front()};
            ScoredPath link_aln = find_link_alignment(to_chain, aln, edge.begin(), edge.begin() + 1, wfa_extender, aligner, stats);

            if (link_aln.score == -std::numeric_limits<int32_t>::max()) {
                // We gave up on this link. It isn't usable, but other,
                // unrelated connections in this group may still be fine.
                continue;
            }

            // Remember the path
            edge_paths[extra_edge] = link_aln.path;

            // Create & save edge
            connection_t* connection = mp_aln.mutable_subpath(extra_edge.first)->add_connection();
            connection->set_next(extra_edge.second);
            connection->set_score(link_aln.score);
        }
    }

    // Do DP, finding all possible paths through the graph
    vector<Alignment> tracebacks = optimal_alignments(mp_aln, std::numeric_limits<int32_t>::max());

    // Convert back to real alignments
    vector<Alignment> output;

    // Look for duplicate alignments by using this collection of node IDs and orientations
    std::unordered_set<std::pair<nid_t, bool>> used_nodes;
    // Compute the fraction of an alignment that is unique
    auto get_fraction_unique = [&](const Path& new_path) {
        // Work out how much of this alignment is from nodes not claimed by previous alignments
        size_t from_length_from_used = 0;
        size_t from_length_total = 0;
        for (size_t i = 0; i < new_path.mapping_size(); i++) {
            // For every mapping
            auto& mapping = new_path.mapping(i);
            auto& position = mapping.position();
            size_t from_length = mapping_from_length(mapping);
            std::pair<nid_t, bool> key{position.node_id(), position.is_reverse()};
            if (used_nodes.count(key)) {
                // Count the from_length on already-used nodes
                from_length_from_used += from_length;
            }
            // And the overall from length
            from_length_total += from_length;
        }
        double unique_node_fraction = from_length_total > 0 ? ((double)(from_length_total - from_length_from_used) / from_length_total) : 1.0;
        return unique_node_fraction;
    };

    // Mark the nodes visited by an alignment as used for uniqueness.
    auto mark_nodes_used = [&](const Path& new_path) {
        for (size_t i = 0; i < new_path.mapping_size(); i++) {
            // For every mapping
            auto& position = new_path.mapping(i).position();
            std::pair<nid_t, bool> key{position.node_id(), position.is_reverse()};
            // Make sure we know we used the oriented node.
            used_nodes.insert(key);
        }
    };

    for (const auto& trace : tracebacks) {
        // Figure out which subchains were used
        vector<size_t> subchains_used;
        subchains_used.reserve(trace.path().mapping_size());
#ifdef debug_base_level_alignment
        cerr << "Possible alignment " << output.size() << " uses subchains: ";
#endif
        for (const auto& mapping : trace.path().mapping()) {
            subchains_used.push_back(mapping.position().node_id());
#ifdef debug_base_level_alignment
            cerr << subchains_used.back() << " ";
#endif
        }
#ifdef debug_base_level_alignment
        cerr << endl;
#endif

        // Build up total base-level results for this alignment
        Path composed_path = node_paths[subchains_used.front()];
        for (size_t i = 1; i < subchains_used.size(); i++) {
            // Add path for edge from previous thing
            append_path(composed_path, edge_paths[make_pair(subchains_used[i-1], subchains_used[i])]);
            // Add path for current thing
            append_path(composed_path, node_paths[subchains_used[i]]);
        }

        // Do the unique node fraction filter
        double unique_node_fraction = get_fraction_unique(composed_path);
        if (unique_node_fraction < min_unique_node_fraction) {
            // If not enough of the alignment is from unique nodes, drop it.
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Possible alignment " << output.size()
                         << " rejected because only " << unique_node_fraction << " of it is from nodes not already used" << endl;
                }
            }
            // Don't bother saving this alignment; it's a duplicate
            continue;
        } else {
            mark_nodes_used(composed_path);
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Possible alignment " << output.size()
                         << " accepted because " << unique_node_fraction << " of it is from nodes not already used" << endl;
                }
            }
        }

        if (track_provenance) {
            // Tell the funnel
            funnel.introduce();
            funnel.score(funnel.latest(), trace.score());
            // We come from all the subchains directly
            funnel.also_merge_group(1, subchains_used.begin(), subchains_used.end());
        }

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Composed alignment is length " << path_to_length(composed_path) << " with score of " << trace.score() << endl;
                if (composed_path.mapping_size() > 0) {
                    cerr << log_name() << "Composed alignment starts with: " << pb2json(composed_path.mapping(0)) << endl;
                    cerr << log_name() << "Composed alignment ends with: " << pb2json(composed_path.mapping(composed_path.mapping_size() - 1)) << endl;
                }
            }
        }

        // Stick into alignment
        output.emplace_back(aln);
        *(output.back()).mutable_path() = std::move(simplify(composed_path, false));
        // Rescore with log-gap penalties
        LoggedGapAlignmentScorer scheme(output.back());
        // Score the alignment
        int32_t logged_gaps_score = scheme.score_alignment(output.back());

        // Penalize score according to the number of recombinations their chains required.
        // This allows alignments that required fewer recombinations in their chains to win.
        // TODO: We'd also eventaully like to count recombinations that we don't know are needed until base-level DP.
        size_t total_rec_count = algorithms::count_total_recombinations(subchain_group, subchains_used);
        if (rec_penalty != 0) {
            logged_gaps_score -= (rec_penalty_aln == -1 ? rec_penalty : rec_penalty_aln) * total_rec_count;
        }

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << " Original score: " << trace.score()
                                   << " Matches: " << scheme.matches
                                   << " Mismatches: " << scheme.mismatches
                                   << " Gap opens: " << scheme.gap_lengths.size()
                                   << " Recombinations: " << total_rec_count
                                   << " New score: " << logged_gaps_score << endl;
            }
        }
        output.back().set_score(logged_gaps_score);
        if (!output.back().sequence().empty()) {
            output.back().set_identity(identity(output.back().path()));
        }
        // Annotate with tail lengths
        set_annotation(output.back(), "left_tail_length", left_tail_len[subchains_used.front()]); 
        set_annotation(output.back(), "right_tail_length", right_tail_len[subchains_used.back()]);
        set_annotation(output.back(), "chain.rec_count", total_rec_count);
    }

    // Sort by new score
    std::sort(output.begin(), output.end(), 
    [](const Alignment& a, const Alignment& b) {
        // Return true if a has the higher score and belongs first
        return a.score() > b.score();
    });
    
    return output;
}

void MinimizerMapper::wfa_alignment_to_alignment(const WFAAlignment& wfa_alignment, Alignment& alignment) const {
    *(alignment.mutable_path()) = wfa_alignment.to_path(this->gbwt_graph, alignment.sequence());
    alignment.set_score(wfa_alignment.score);
    if (!alignment.sequence().empty()) {
        alignment.set_identity(identity(alignment.path()));
    }
}

void MinimizerMapper::with_dagified_local_graph(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, const HandleGraph& graph, const std::function<void(DeletableHandleGraph&, const handle_t&, const handle_t&, const std::function<std::pair<nid_t, bool>(const handle_t&)>&)>& callback) {
    
    if (is_empty(left_anchor) && is_empty(right_anchor)) {
        throw ChainAlignmentFailedError("Cannot align sequence between two unset positions");
    }
    
    // We need to get the graph to align to.
    bdsg::HashGraph local_graph;
    unordered_map<id_t, id_t> local_to_base;
    if (!is_empty(left_anchor) && !is_empty(right_anchor)) {
        // We want a graph actually between two positions.
        // Enforce strict max length to avoid extra tips.
        local_to_base = algorithms::extract_connecting_graph(
            &graph,
            &local_graph,
            max_path_length,
            left_anchor, right_anchor,
            true
        );

        if (local_to_base.empty()) {
            // A possible result is that the one anchor is not reachable from
            // the other within the length limit (at least without doubling
            // back through one of the anchor nodes). In that case, we get an
            // empty graph and an empty translation.

            // Explain the problem but bail out on the chain connection safely,
            // because we we expect this sometimes.
            std::stringstream ss;
            ss << "Cannot find an acceptable path from " << left_anchor;
            ss << " to " << right_anchor;
            ss << " with max path length of " << max_path_length;
            throw ChainAlignmentFailedError(ss.str());
        }
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
        auto& local_id = kv.first;
        auto& base_id = kv.second;
        if (base_id == id(left_anchor) && base_id == id(right_anchor)) {
            // The left and right anchors are on the same node, and this is a copy of it.
            // It could be that the anchors face each other, and we extracted one intervening piece of node.
            // In which case we go through this section once.
            if (local_left_anchor_id == 0 && local_right_anchor_id == 0) {
                // First time through, say we probably cut out the middle piece of a node
                local_left_anchor_id = local_id;
                local_right_anchor_id = local_id;
#ifdef debug
                std::cerr << "Assume left and right anchors are both node " << local_id << " representing " << base_id << std::endl;
#endif
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
                    ss << "Extracted graph of " << local_graph.get_node_count() << " nodes";
                    ss << " from " << left_anchor << " to " << right_anchor;
                    ss << " with max path length of " << max_path_length;
                    ss << " but shared node appeared more than twice in the resulting translation.";
                    ss << " Graph dumped as crashdump.vg.";
                    local_graph.serialize("crashdump.vg");
                    throw std::runtime_error(ss.str());
                }
                // Whichever copy has the lower ID is the left one and
                // whichever copy has the higher ID is the right one.
                local_left_anchor_id = std::min(local_left_anchor_id, local_id);
                local_right_anchor_id = std::max(local_right_anchor_id, local_id);
#ifdef debug
                std::cerr << "Second shared anchor copy as " << local_id << " representing " << base_id << "; left is now " << local_left_anchor_id << " and right is " << local_right_anchor_id << std::endl;
#endif
            }
        } else if (base_id == id(left_anchor)) {
            if (local_left_anchor_id != 0) {
                // We thought we already figured out the start node; there are
                // too many copies of our start node to find it. 
                std::stringstream ss;
                ss << "Extracted graph of " << local_graph.get_node_count() << " nodes";
                ss << " from " << left_anchor << " to " << right_anchor;
                ss << " with max path length of " << max_path_length;
                ss << " but start node appeared twice in the resulting translation.";
                ss << " Graph dumped as crashdump.vg.";
                local_graph.serialize("crashdump.vg");
                throw std::runtime_error(ss.str());
            }
            local_left_anchor_id = local_id;
#ifdef debug
            std::cerr << "Left anchor is " << local_left_anchor_id << std::endl;
#endif
        } else if (base_id == id(right_anchor)) {
            if (local_right_anchor_id != 0) {
                // We thought we already figured out the end node; there are
                // too many copies of our end node to find it. 
                std::stringstream ss;
                ss << "Extracted graph of " << local_graph.get_node_count() << " nodes";
                ss << " from " << left_anchor << " to " << right_anchor;
                ss << " with max path length of " << max_path_length;
                ss << " but end node appeared twice in the resulting translation.";
                ss << " Graph dumped as crashdump.vg.";
                local_graph.serialize("crashdump.vg");
                throw std::runtime_error(ss.str());
            }
            local_right_anchor_id = local_id;
#ifdef debug
            std::cerr << "Right anchor is " << local_right_anchor_id << std::endl;
#endif
        }
        // TODO: Stop early when we found them all.
    }

    if (!is_empty(left_anchor) && local_left_anchor_id == 0) {
        // Somehow the left anchor didn't come through. Complain.
        std::stringstream ss;
        ss << "Extracted graph of " << local_graph.get_node_count() << " nodes";
        ss << " from " << left_anchor << " to " << right_anchor;
        ss << " with max path length of " << max_path_length;
        ss << " but from node was not present in the resulting translation.";
        ss << " Graph dumped as crashdump.vg.";
        local_graph.serialize("crashdump.vg");
        throw std::runtime_error(ss.str());
    }

    if (!is_empty(right_anchor) && local_right_anchor_id == 0) {
        // Somehow the right anchor didn't come through. Complain.
        std::stringstream ss;
        ss << "Extracted graph of " << local_graph.get_node_count() << " nodes";
        ss << " from " << left_anchor << " to " << right_anchor;
        ss << " with max path length of " << max_path_length;
        ss << " but to node was not present in the resulting translation.";
        ss << " Graph dumped as crashdump.vg.";
        local_graph.serialize("crashdump.vg");
        throw std::runtime_error(ss.str());
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
            ss << " but from node local ID " << local_left_anchor_id << " from translation was not present in the resulting graph.";
            ss << " Graph dumped as crashdump.vg.";
            local_graph.serialize("crashdump.vg");
            throw std::runtime_error(ss.str());
        }
        handle_t local_handle = local_graph.get_handle(local_left_anchor_id, is_rev(left_anchor));
        
        // And get the node that that orientation of it is in the strand-split graph
        handle_t overlay_handle = split_graph.get_overlay_handle(local_handle);

#ifdef debug
        std::cerr << "Left anchor " << local_graph.get_id(local_handle) << (local_graph.get_is_reverse(local_handle) ? "-" : "+") << " in local graph is " << split_graph.get_id(overlay_handle) << (split_graph.get_is_reverse(overlay_handle) ? "-" : "+") << " in strand-split graph." << std::endl;
#endif
        
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
            ss << " but to node local ID " << local_right_anchor_id << " from translation was not present in the resulting graph.";
            ss << " Graph dumped as crashdump.vg.";
            local_graph.serialize("crashdump.vg");
            throw std::runtime_error(ss.str());
        }
        handle_t local_handle = local_graph.get_handle(local_right_anchor_id, is_rev(right_anchor));
        
        // And get the node that that orientation of it is in the strand-split graph
        // But flip it because we want to dagify going inwards from the right
        handle_t overlay_handle = split_graph.flip(split_graph.get_overlay_handle(local_handle));

#ifdef debug
        std::cerr << "Right anchor " << local_graph.get_id(local_handle) << (local_graph.get_is_reverse(local_handle) ? "-" : "+") << " in local graph is " << split_graph.get_id(overlay_handle) << (split_graph.get_is_reverse(overlay_handle) ? "-" : "+") << " in strand-split graph." << std::endl;
#endif
        
        // And use that
        bounding_handles.push_back(overlay_handle);
    }
    
    // Do the dagification from those input handles.
    // TODO: Note that this can add tips! We should come up with a dagification method that is guaranteed not to!
    auto dagification_result = handlegraph::algorithms::dagify_from(&split_graph, bounding_handles, &dagified_graph, max_path_length);
    auto& dagified_to_split = dagification_result.first;
    auto& anchor_handles = dagification_result.second;
    
    // Figure out which anchoring handle is which anchor.
    handle_t left_anchor_handle = is_empty(left_anchor) ? handle_t() : anchor_handles.front();
    // The right handle will be in facing-into-the-graph orientation, but we
    // want to produce it in facing outwards orientation like the right anchor
    // position was. 
    handle_t right_anchor_handle = is_empty(right_anchor) ? handle_t() : dagified_graph.flip(anchor_handles.back());

    // Note that in addition to cut nodes at these anchor handles, we can have
    // strand-split version of the cut nodes' other strands.

    
#ifdef debug
    std::cerr << "Dagified from " << bounding_handles.size() << " bounding handles in " << split_graph.get_node_count() << " node strand-split graph to " << dagified_graph.get_node_count() << " node DAG" << std::endl;
#endif

    // Make an accessor for getting back to the base graph space
    auto dagified_handle_to_base = [&](const handle_t& h) -> pair<nid_t, bool> {
        nid_t dagified_id = dagified_graph.get_id(h);
        bool dagified_is_reverse = dagified_graph.get_is_reverse(h);
        auto found_in_split = dagified_to_split.find(dagified_id);
        if (found_in_split == dagified_to_split.end()) {
            throw std::runtime_error("ID " + std::to_string(dagified_id) + " from dagified graph not found in strand-split graph");
        }
        nid_t split_id = found_in_split->second;
        handle_t split_handle = split_graph.get_handle(split_id, dagified_is_reverse);
        // We rely on get_underlying_handle understanding reversed handles in the split graph
        handle_t local_handle = split_graph.get_underlying_handle(split_handle);
        nid_t local_id = local_graph.get_id(local_handle);
        bool local_is_reverse = local_graph.get_is_reverse(local_handle);
        auto found_in_base = local_to_base.find(local_id);
        if (found_in_base == local_to_base.end()) {
            throw std::runtime_error("ID " + std::to_string(local_id) + " from local graph not found in full base graph");
        }
        nid_t base_id = found_in_base->second;
        return std::make_pair(base_id, local_is_reverse);
    };
    
    // Show the graph we made and the translation function
    callback(dagified_graph, left_anchor_handle, right_anchor_handle, dagified_handle_to_base);
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
        return aligner->scorer->longest_detectable_gap(aln, aln.sequence().begin() + middle_index);
    }
    
    // Otherwise it is the length from the boundary nearest to the middle.
    // And we know the while range is on one side or the other of the middle.
    if (begin_index > middle_index) {
        // Beginning is on the inside
        return aligner->scorer->longest_detectable_gap(aln, sequence_begin);
    }

    // Otherwise the end is on the inside
    return aligner->scorer->longest_detectable_gap(aln, sequence_end);
}

bool MinimizerMapper::align_sequence_between(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, size_t max_gap_length, const HandleGraph* graph, const GSSWAligner* aligner, Alignment& alignment, const std::string* alignment_name, size_t max_dp_cells, const std::function<size_t(const Alignment&, const HandleGraph&)>& choose_band_padding) {

    bool did_aln = true;
    // Get the dagified local graph, and the back translation
    MinimizerMapper::with_dagified_local_graph(left_anchor, right_anchor, max_path_length, *graph,
        [&](DeletableHandleGraph& dagified_graph,
            const handle_t& left_anchor_handle,
            const handle_t& right_anchor_handle,
            const std::function<std::pair<nid_t, bool>(const handle_t&)>& dagified_handle_to_base) {

#ifdef debug
        dump_debug_graph(dagified_graph);
#endif
        
        // Then trim off the tips that are either in the wrong orientation relative
        // to whether we want them to be a source or a sink, or extraneous.
        // Also find the unique tip handles for the potentially cut-down anchor nodes.
        // TODO: We re-find them on every trim.
        
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
                if (!dagified_graph.get_is_reverse(h) && (is_empty(left_anchor) || h == left_anchor_handle)) {
                    // Tip is inward forward, so it's a source.
                    // This is a head in the subgraph, and either matches the
                    // left anchor or we don't have any, so keep it.
#ifdef debug
                    std::cerr << "Dagified graph node " << dagified_graph.get_id(h) << " " << dagified_graph.get_is_reverse(h) << " is an acceptable source (" << base_coords.first << " " << base_coords.second << ") length " << dagified_graph.get_length(h) << std::endl;
#endif
                } else if (dagified_graph.get_is_reverse(h) && (is_empty(right_anchor) || h == dagified_graph.flip(right_anchor_handle))) {
                    // Tip is inward reverse, so it's a sink.
                    // This is a tail in the subgraph, and either matches the right
                    // anchor (which faces outwards) or we don't have any, so keep it.
#ifdef debug
                    std::cerr << "Dagified graph node " << dagified_graph.get_id(h) << " " << dagified_graph.get_is_reverse(h) << " is an acceptable sink (" << base_coords.first << " " << base_coords.second << ") length " << dagified_graph.get_length(h) << std::endl;
#endif
                } else {
                    // This is a wrong orientation or other copy of an anchoring node, or some other tip.
                    // We don't want to keep this handle
#ifdef debug
                    std::cerr << "Dagified graph node " << dagified_graph.get_id(h) << " " << dagified_graph.get_is_reverse(h) << " is an unacceptable tip (" << base_coords.first << " " << base_coords.second << ") length " << dagified_graph.get_length(h) << std::endl;
#endif
                    nid_t dagified_id = dagified_graph.get_id(h);
                    if (!to_remove_ids.count(dagified_id)) {
                        to_remove_ids.insert(dagified_id);
                        to_remove_handles.push_back(h);
                    }
                }
            }
#ifdef debug_dump_tips
            if (!to_remove_handles.empty() && trim_count == 0) {
                // We're going to trim, so dump the graph.
                std::cerr << "warning[MinimizerMapper::align_sequence_between]: Going to trim, so dumping pre-trim.vg" << std::endl;
                dynamic_cast<SerializableHandleGraph*>(&dagified_graph)->serialize("pre-trim.vg");
            }
#endif
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
                std::cerr << "warning[MinimizerMapper::align_sequence_between]: Trimmed back tips " << trim_count << " times on graph between " << left_anchor;
                if (!is_empty(left_anchor)) {
                    std::cerr << " as " << dagified_graph.get_id(left_anchor_handle) << (dagified_graph.get_is_reverse(left_anchor_handle) ? "-" : "+");
                }
                std::cerr << " and " << right_anchor;
                if (!is_empty(right_anchor)) {
                    std::cerr << " as " << dagified_graph.get_id(right_anchor_handle) << (dagified_graph.get_is_reverse(right_anchor_handle) ? "-" : "+");
                }
                std::cerr << " leaving " <<  dagified_graph.get_node_count() << " nodes and " << tip_handles.size() << " tips";
                if (alignment_name) {
                    std::cerr << " for read " << *alignment_name;
                }
#ifdef debug_dump_tips
                std::cerr << " to make post-trim.vg";
                dynamic_cast<SerializableHandleGraph*>(&dagified_graph)->serialize("post-trim.vg");
#endif
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
            try {
                aligner->align_global_banded(alignment, dagified_graph, band_padding, true, max_dp_cells);
            } catch (BandMatricesTooBigException& e) {
                // We would use too many DP cells.
                #pragma omp critical (cerr)
                {
                    std::cerr << "warning[MinimizerMapper::align_sequence_between]: " << e.what() << std::endl;
                }
                // Clear out the alignment path to indicate that we didn't actually compute an alignment.
                alignment.mutable_path()->clear_mapping();
            }
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
                did_aln = false;
                return;
            } else {
#ifdef debug
                #pragma omp critical (cerr)
                std::cerr << "debug[MinimizerMapper::align_sequence_between]: Fill " << cell_count << " DP cells in tail with Xdrop" << std::endl;
#endif
                aligner->align_pinned(alignment, dagified_graph, !is_empty(left_anchor), true, max_gap_length);
            }
        }

        // And translate back into original graph ID and orientation space.
        for (size_t i = 0; i < alignment.path().mapping_size(); i++) {
            // Translate each mapping's ID and orientation down to the base graph
            Mapping* m = alignment.mutable_path()->mutable_mapping(i);
            
            handle_t dagified_handle = dagified_graph.get_handle(m->position().node_id(), m->position().is_reverse());
            auto base_coords = dagified_handle_to_base(dagified_handle);
           
            if (i == 0) {
                if (!is_empty(left_anchor) && base_coords.first == id(left_anchor) && base_coords.second == is_rev(left_anchor)) {
                    // The alignment starts on the left anchor's node. It may be on a cut copy.
                    if (dagified_graph.get_length(dagified_handle) == dagified_graph.get_length(left_anchor_handle)) {
                        // We're on a cut copy (the anchor itself or the strand-split version), or it was not cut.

                        // Add on the offset for the cut-off outside piece of the left anchor node
                        m->mutable_position()->set_offset(m->position().offset() + offset(left_anchor));
                    }
                } else if (!is_empty(right_anchor) && base_coords.first == id(right_anchor) && base_coords.second != is_rev(right_anchor)) {
                    // The alignment starts on the right anchor's node, reading into the graph. It may be on a cut copy.
                    if (dagified_graph.get_length(dagified_handle) == dagified_graph.get_length(right_anchor_handle)) {
                        // We're on a cut copy (the anchor itself or the strand-split version), or it was not cut.

                        // Add on the offset for the cut-off outside piece of the right anchor node
                        m->mutable_position()->set_offset(m->position().offset() + graph->get_length(graph->get_handle(id(right_anchor))) - offset(right_anchor));
                    }
                }
            }

            m->mutable_position()->set_node_id(base_coords.first);
            m->mutable_position()->set_is_reverse(base_coords.second);
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

    return did_aln;
}

bool MinimizerMapper::align_sequence_between_consistently(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, size_t max_gap_length, const HandleGraph* graph, const GSSWAligner* aligner, Alignment& alignment, const std::string* alignment_name, size_t max_dp_cells, const std::function<size_t(const Alignment&, const HandleGraph&)>& choose_band_padding) {
    if (left_anchor < right_anchor) {
        // Left anchor is unambiguously first, so align as-is
        return align_sequence_between(left_anchor, right_anchor, max_path_length, max_gap_length, graph, aligner, alignment, alignment_name, max_dp_cells, choose_band_padding);
    }

    // Otherwise left anchor is equal or greater.

    // Make a node length getter for flipping alignments
    auto get_node_length = [&](id_t node_id) -> int64_t {
        return graph->get_length(graph->get_handle(node_id));
    };
    

    // Compute the reverse-complement sequence, which we either need or might break the tie.
    Alignment flipped_query = reverse_complement_alignment(alignment, get_node_length);

    if (left_anchor == right_anchor && flipped_query.sequence() >= alignment.sequence()) {
        // The anchors are tied and the sequence doesn't demand a switch. Align as-is.
        //
        // TODO: For palindromic sequences aligned between identical endpoints,
        // we still might get inconsistencies by read strand in the final
        // output, since the read around it might be in either orientation
        // relative to the flow of the reference.
        return align_sequence_between(left_anchor, right_anchor, max_path_length, max_gap_length, graph, aligner, alignment, alignment_name, max_dp_cells, choose_band_padding);
    }

    // Now we know a swap is required.
    
    // The anchors face left to right so we need to flip their orientations in addition to swapping them.
    // align_sequence_between uses between-base positions for anchoring
    pos_t flipped_left_anchor = is_empty(right_anchor) ? empty_pos_t() : reverse(right_anchor, get_node_length(id(right_anchor)));
    pos_t flipped_right_anchor = is_empty(left_anchor) ? empty_pos_t() : reverse(left_anchor, get_node_length(id(left_anchor)));

    // Do the alignment
    auto result = align_sequence_between(flipped_left_anchor, flipped_right_anchor, max_path_length, max_gap_length, graph, aligner, flipped_query, alignment_name, max_dp_cells, choose_band_padding);

    // Flip and send the answer
    reverse_complement_alignment_in_place(&flipped_query, get_node_length);
    alignment = std::move(flipped_query);

    // We shouldn't use any nonzero offsets after the first node. We're not
    // meant to be making a split alignment.
    for (size_t i = 0; i < alignment.path().mapping_size(); i++) {
        if (i > 0) {
            crash_unless(alignment.path().mapping(i).position().offset() == 0);
        }
    }

    // Return the metadata we track
    return result;
}

WFAAlignment MinimizerMapper::connect_consistently(const std::string& sequence, const pos_t& left_anchor, const pos_t& right_anchor, const WFAExtender& wfa_extender) {

    // TODO: Deduplicate swap logic with align_sequence_between_consistently

    if (left_anchor < right_anchor) {
        // Left anchor is unambiguously first, so align as-is
        return wfa_extender.connect(sequence, left_anchor, right_anchor);
    }

    // Otherwise left anchor is equal or greater.
    // Compute the reverse-complement sequence, which we either need or might break the tie.
    std::string flipped_sequence = reverse_complement(sequence);

    if (left_anchor == right_anchor && flipped_sequence >= sequence) {
        // The anchors are tied and the sequence doesn't demand a switch. Align as-is.
        //
        // TODO: For palindromic sequences aligned between identical endpoints,
        // we still might get inconsistencies by read strand in the final
        // output, since the read around it might be in either orientation
        // relative to the flow of the reference.
        return wfa_extender.connect(sequence, left_anchor, right_anchor);
    }

    // Now we know a swap is required.
    
    // TODO: We probably don't *really* need to track orientation here
    handle_t left_handle = wfa_extender.graph->get_handle(id(left_anchor), is_rev(left_anchor));
    handle_t right_handle = wfa_extender.graph->get_handle(id(right_anchor), is_rev(right_anchor));
    
    // The anchors face left to right so we need to flip their orientations in addition to swapping them.
    // Also note that WFAExtender works with base positions and not intervening positions.
    pos_t flipped_left_anchor = reverse_base_pos(right_anchor, wfa_extender.graph->get_length(right_handle));
    pos_t flipped_right_anchor = reverse_base_pos(left_anchor, wfa_extender.graph->get_length(left_handle));
   
    // Make the reverse alignment
    WFAAlignment result = wfa_extender.connect(flipped_sequence, flipped_left_anchor, flipped_right_anchor);
    
    // Put the alignment back, which needs the final alignment's sequence (see WFAExtender's prefix() implementation)
    result.flip(*wfa_extender.graph, sequence);
    
    // And ship it
    return result;
}

std::vector<algorithms::Anchor> MinimizerMapper::to_anchors(const Alignment& aln, const VectorView<Minimizer>& minimizers, std::vector<Seed>& seeds) const {
    std::vector<algorithms::Anchor> to_return;
    to_return.reserve(seeds.size());
    for (size_t i = 0; i < seeds.size(); i++) {
        to_return.push_back(MinimizerMapper::to_anchor(aln, minimizers, seeds, i, gbwt_graph, get_regular_aligner()));
    }
    return to_return;
}

algorithms::Anchor MinimizerMapper::to_anchor(const Alignment& aln, const VectorView<Minimizer>& minimizers, std::vector<Seed>& seeds, size_t seed_number, const HandleGraph& graph, const Aligner* aligner) {
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
    auto paths = seed.paths;
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
        read_start = source.pin_offset() + 1 - length;
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
        read_start = source.pin_offset();
        // The seed is actually at the start
        hint_start = 0;
    }

#ifdef debug
    std::cerr << "Minimizer at read " << source.forward_offset() << " length " << source.length
              << " orientation " << source.value.is_reverse << " pinned at " << source.pin_offset()
              << " is anchor of length " << length << " matching graph " << graph_start << " and read " << read_start
              << " forward, with hint " << hint_start << " bases later on the read" << std::endl;
#endif

    // Work out how many points the anchor is.
    // TODO: Always make sequence and quality available for scoring!
    // We're going to score the anchor as the full minimizer, and rely on the margins to stop us from taking overlapping anchors.
    int score = aligner->scorer->score_exact_match(aln, read_start - margin_left, margin_left + length + margin_right);
    return algorithms::Anchor(read_start, graph_start, length, margin_left, margin_right, score, seed_number, &(seed.zipcode), hint_start, source.is_repetitive, paths); 
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
        score += aligner->scorer->score_exact_match(aln, scored_until, *mismatch_it - scored_until);
        score += aligner->scorer->score_mismatch(aln.sequence().begin() + *mismatch_it,
                                         aln.sequence().begin() + *mismatch_it + 1,
                                         aln.quality().begin() + *mismatch_it); 
        scored_until = *mismatch_it + 1;
        ++mismatch_it;
    }
    // Score the perfect match from where we are to the end.
    score += aligner->scorer->score_exact_match(aln, scored_until, read_end - scored_until);
    
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
    auto score = aligner->scorer->score_exact_match(aln, anchor.read_start(), anchor.length());
    if (anchor.read_start() == 0) {
        // Apply full elngth bonus on the left if we abut the left end of the read.
        score += aligner->scorer->score_full_length_bonus(true, aln);
    }
    if (anchor.read_end() == aln.sequence().length()) {
        // Apply full lenght bonus on the right if we abut the riht end of the read.
        score += aligner->scorer->score_full_length_bonus(false, aln);
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

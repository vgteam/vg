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

void MinimizerMapper::dump_debug_dotplot(const std::string& name, const std::string& marker, const VectorView<Minimizer>& minimizers, const std::vector<Seed>& seeds, const std::vector<size_t>& included_seeds, const std::vector<size_t>& highlighted_seeds, const PathPositionHandleGraph* path_graph) {
    if (!path_graph) {
        // We don't have a path positional graph for this
        return;
    }

    // Log the best bucket's seed positions in read and linear reference
    TSVExplainer exp(true, name + "-dotplot");

    // We need to know which seeds to highlight
    std::unordered_set<size_t> highlight_set;
    for (auto& seed_num : highlighted_seeds) {
        highlight_set.insert(seed_num);
    }

    for (auto& seed_num : included_seeds) {
        // For each seed in the best bucket
        auto& seed = seeds.at(seed_num);
        
        // Get its effective path positions again
        auto offsets = algorithms::nearest_offsets_in_paths(path_graph, seed.pos, 100);

        for (auto& handle_and_positions : offsets) {
            std::string path_name = path_graph->get_path_name(handle_and_positions.first);
            for (auto& position : handle_and_positions.second) {
                // For each position on a ref path that this seed is at, log a line
                exp.line();
                if (highlight_set.count(seed_num)) {
                    // Contig and a marker
                    exp.field(path_name + "-" + marker);
                } else {
                    // Contig
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

void MinimizerMapper::dump_debug_graph(const HandleGraph& graph) {
    graph.for_each_handle([&](const handle_t& h) {
        std::cerr << "Node " << graph.get_id(h) << ": " << graph.get_sequence(h) << std::endl;
    });
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
                std::cerr << log_name() << "Positions for tree " << i << ":" << std::endl;
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

vector<Alignment> MinimizerMapper::map_from_chains(Alignment& aln) {
    
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

    vector<ZipCodeDecoder> decoders;
    
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

    if (show_work) {
        #pragma omp critical (cerr)
        {
            std::cerr << log_name() << "Found " << zip_code_forest.trees.size() << " zip code trees, scores " << best_tree_score << " best, " << second_best_tree_score << " second best, coverages " << best_tree_coverage << " best, " << second_best_tree_coverage << " second best" << std::endl;
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

    // Convert the seeds into chainable anchors in the same order
    vector<algorithms::Anchor> seed_anchors = this->to_anchors(aln, minimizers, seeds);

    // Now compute fragments into these variables.
    // What seeds are visited in what order in the fragment?
    std::vector<std::vector<size_t>> fragments;
    // What score does each fragment have?
    std::vector<double> fragment_scores;
    // Which zip code tree did each fragment come from, so we know how to chain them?
    std::vector<size_t> fragment_source_tree;
    // How many of each minimizer ought to be considered explored by each fragment?
    // TODO: This is a lot of counts and a lot of allocations and should maybe be a 2D array if we really need it?
    std::vector<std::vector<size_t>> minimizer_kept_fragment_count;

    process_until_threshold_c<double>(zip_code_forest.trees.size(), [&](size_t i) -> double {
            return tree_scores[i];
        }, [&](size_t a, size_t b) -> bool {
            return tree_scores[a] > tree_scores[b] || (tree_scores[a] == tree_scores[b] && tree_coverages[a] > tree_coverages[b]); 
        }, zipcode_tree_score_threshold, this->min_to_fragment, this->max_to_fragment, rng, [&](size_t item_num) -> bool {
            // Handle sufficiently good fragmenting problems in descending score order
            
            if (track_provenance) {
                funnel.pass("fragmenting-score", item_num, tree_scores[item_num]);
                funnel.pass("max-to-fragment", item_num);
                funnel.pass("fragmenting-coverage", item_num, tree_coverages[item_num]);
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Making fragments for zip code tree " << item_num << " with score " << tree_scores[item_num] << " and coverage " << tree_coverages[item_num] << endl;
                }
            }
            
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
           
            // Sort seeds by read start of seeded region
            algorithms::sort_anchor_indexes(seed_anchors, selected_seeds);
            
            if (track_provenance) {
                funnel.substage("find_fragment");
            }
            
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Computing fragments over " << selected_seeds.size() << " seeds" << endl;
                }
            }

#ifdef debug
            if (show_work) {
                // Log the chaining problem so we can try it again elsewhere.
                this->dump_chaining_problem(seed_anchors, selected_seeds, gbwt_graph);
            }
#endif
            
            // Find fragments over the seeds in the zip code tree
            algorithms::transition_iterator for_each_transition = algorithms::zip_tree_transition_iterator(
                seeds,
                zip_code_forest.trees[item_num],
                this->fragment_max_lookback_bases
            ); 
            VectorView<algorithms::Anchor> anchor_view {seed_anchors, selected_seeds};
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
                this->fragment_max_indel_bases,
                false // Don't show work for fragmenting, there are too many seeds.
            );
            if (show_work) {
                #pragma omp critical (cerr)
                cerr << log_name() << "Found " << results.size() << " fragments in zip code tree " << item_num
                    << " running " << seed_anchors[selected_seeds.front()] << " to " << seed_anchors[selected_seeds.back()] << std::endl;
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
                fragments.back().reserve(scored_fragment.second.size());
                for (auto& selected_number : scored_fragment.second) {
                    // Translate from selected seed/anchor space to global seed space.
                    fragments.back().push_back(selected_seeds[selected_number]);
                    // And count the minimizer as being in the fragment
                    minimizer_kept_fragment_count.back()[seeds[fragments.back().back()].source]++;
                }
                // Remember the score
                fragment_scores.push_back(scored_fragment.first);
                // Remember how we got it
                fragment_source_tree.push_back(item_num);

                if (track_provenance) {
                    // Tell the funnel
                    funnel.introduce();
                    funnel.score(funnel.latest(), scored_fragment.first);
                    // We come from all the seeds directly
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
                funnel.pass("fragmenting-score", item_num, tree_scores[item_num]);
                funnel.fail("max-to-fragment", item_num);
            }
            
        }, [&](size_t item_num) -> void {
            // This item is not sufficiently good.
            if (track_provenance) {
                funnel.fail("fragmenting-score", item_num, tree_scores[item_num]);
            }
        });

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
    
    // Make a list of anchors where we have each fragment as itself an anchor
    std::vector<algorithms::Anchor> fragment_anchors;
    fragment_anchors.reserve(fragments.size());
    for (size_t i = 0; i < fragments.size(); i++) {
        auto& fragment = fragments.at(i);
        auto& score = fragment_scores.at(i);
        fragment_anchors.push_back(algorithms::Anchor(seed_anchors.at(fragment.front()), seed_anchors.at(fragment.back()), score));
    }
    
    // Get all the fragment numbers for each zip code tree we actually used, so we can chain each independently again.
    // TODO: Stop reswizzling so much.
    std::unordered_map<size_t, std::vector<size_t>> tree_to_fragments;
    for (size_t i = 0; i < fragment_source_tree.size(); i++) {
        tree_to_fragments[fragment_source_tree[i]].push_back(i);
    }
    
    // Get the score of the top-scoring fragment in each collection.
    std::unordered_map<size_t, double> best_fragment_score_in;
    for (auto& kv : tree_to_fragments) {
        for (auto& fragment_num : kv.second) {
            // Max in the score of each fragment 
            best_fragment_score_in[kv.first] = std::max(best_fragment_score_in[kv.first], fragment_scores.at(fragment_num));
        }
    }
    
    // Filter down to just the good ones, sorted by read start
    // TODO: Should we drop short fragments in one place because of long fragments in a *different* place?
    // TODO: If not, can we just immediately chain the results of each fragmenting run?
    std::unordered_map<size_t, std::vector<size_t>> good_fragments_in;
    for (auto& kv : tree_to_fragments) {
        // Decide on how good fragments have to be to keep.
        double fragment_score_threshold = best_fragment_score_in.at(kv.first) * fragment_score_fraction;
    
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Keeping, of the " << kv.second.size() << " fragments in " << kv.first << ", those with score of at least "  << fragment_score_threshold << endl;
            }
        }
    
        // Keep the fragments that have good scores.
        for (auto& fragment_num : kv.second) {
            // For each fragment
            if (fragment_scores.at(fragment_num) >= fragment_score_threshold) {
                // If its score is high enough, keep it.
                // TODO: Tell the funnel.
                good_fragments_in[kv.first].push_back(fragment_num);
            }
        }
        
        // Now sort anchors by read start. Don't bother with shadowing.
        algorithms::sort_anchor_indexes(fragment_anchors, good_fragments_in[kv.first]);

        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "\tKept " << good_fragments_in[kv.first].size() << "/" << kv.second.size() << " fragments." << endl;
            }
        }
    }
    
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "=====Creating chains=====" << endl;
        }
    }

    for (auto& kv : good_fragments_in) {
        auto& tree_num = kv.first;
        // Get a view of all the good fragments.
        // TODO: Should we just not make a global fragment anchor list?
        VectorView<algorithms::Anchor> fragment_view {fragment_anchors, kv.second};

        if (fragment_view.empty()) {
            // Nothing to chain!
            if (show_work) {
                #pragma omp critical (cerr)
                std::cerr << log_name() << "Zip code tree " << tree_num << " has no good fragments to chain!" << std::endl;
            } 
            continue;
        }
        
        if (show_work) {
            #pragma omp critical (cerr)
            std::cerr << log_name() << "Chaining fragments from zip code tree " << tree_num << std::endl;
        } 

        // Chain up the fragments
        algorithms::transition_iterator for_each_transition = algorithms::zip_tree_transition_iterator(
            seeds,
            zip_code_forest.trees[tree_num],
            this->max_lookback_bases
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
            this->max_indel_bases,
            this->show_work
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
            
            // We record the fragments that merge into each chain for reporting.
            std::vector<size_t> chain_fragment_nums_overall;
            chain_fragment_nums_overall.reserve(chain_result.second.size());
            
            for (const size_t& local_fragment: chain_result.second) {
                // For each fragment in the chain
                           
                // Get its fragment number out of all fragments
                size_t fragment_num_overall = kv.second.at(local_fragment);
                
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
                        std::cerr << log_name() << "Chain " << (chains.size() - 1) << " with score " << score << " is composed from fragments:";
                        for (auto& f : chain_fragment_nums_overall) {
                            std::cerr << " " << f;
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
    }
    
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
        vector<size_t> involved_seeds;
        for (ZipCodeTree::oriented_seed_t found : zip_code_forest.trees.at(chain_source_tree.at(best_chain))) {
            involved_seeds.push_back(found.seed);   
        }
        dump_debug_dotplot("best-chain", "chain", minimizers, seeds, involved_seeds, chains.at(best_chain), this->path_graph);
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
            funnel.fail("chain-score", processed_num, chain_score_estimates[processed_num]);
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
    
    // Go through the chains in estimated-score order.
    process_until_threshold_b<int>(chain_score_estimates,
        chain_score_threshold, min_chains, max_alignments, rng, [&](size_t processed_num) -> bool {
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
                funnel.pass("chain-score", processed_num, chain_score_estimates[processed_num]);
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
                
                // We currently just have the one best score and chain per zip code tree
                vector<size_t>& chain = chains[processed_num];
                
                try {
                    // Do the DP between the items in the chain. 
                    best_alignments[0] = find_chain_alignment(aln, seed_anchors, chain);
                } catch (ChainAlignmentFailedError& e) {
                    // We can't actually make an alignment from this chain
                    #pragma omp critical (cerr)
                    cerr << log_name() << "Error creating alignment from chain for " << aln.name() << ": " << e.what() << endl;
                    // Leave the read unmapped.
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
            
            return true;
        }, [&](size_t processed_num) -> void {
            // There are too many sufficiently good chains
            if (track_provenance) {
                funnel.pass("chain-score", processed_num, chain_score_estimates[processed_num]);
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
    
    if (alignments.size() == 0) {
        // Produce an unaligned Alignment
        alignments.emplace_back(aln);
        alignments_to_source.push_back(numeric_limits<size_t>::max());
        
        if (track_provenance) {
            // Say it came from nowhere
            funnel.introduce();
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
    
    process_until_threshold_a(alignments.size(), (std::function<double(size_t)>) [&](size_t i) -> double {
        return alignments.at(i).score();
    }, 0, 1, max_multimaps, rng, [&](size_t alignment_num) {
        // This alignment makes it
        // Called in score order
        
        // Remember the score at its rank
        scores.emplace_back(alignments[alignment_num].score());
        
        // Remember the output alignment
        mappings.emplace_back(std::move(alignments[alignment_num]));
        
        if (track_provenance) {
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
        
        if (track_provenance) {
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
        get_regular_aligner()->compute_max_mapping_quality(scaled_scores, false) ;
    
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
    set_annotation(mappings.front(),"secondary_scores", scores);

    if (track_provenance) {
        funnel.substage_stop();
    }
    
    for (size_t i = 0; i < mappings.size(); i++) {
        // For each output alignment in score order
        auto& out = mappings[i];
        
        // Assign primary and secondary status
        out.set_is_secondary(i > 0);
    }
    
    // Stop this alignment
    funnel.stop();
    
    // Annotate with whatever's in the funnel
    funnel.annotate_mapped_alignment(mappings[0], track_correctness);
    
    if (track_provenance) {
        if (track_correctness) {
            annotate_with_minimizer_statistics(mappings[0], minimizers, seeds, seeds.size(), fragments.size(), funnel);
        }
        // Annotate with parameters used for the filters and algorithms.
        
        set_annotation(mappings[0], "param_hit-cap", (double) hit_cap);
        set_annotation(mappings[0], "param_hard-hit-cap", (double) hard_hit_cap);
        set_annotation(mappings[0], "param_score-fraction", (double) minimizer_score_fraction);
        set_annotation(mappings[0], "param_max-unique-min", (double) max_unique_min);
        set_annotation(mappings[0], "param_num-bp-per-min", (double) num_bp_per_min);
        set_annotation(mappings[0], "param_exclude-overlapping-min", exclude_overlapping_min);
        set_annotation(mappings[0], "param_align-from-chains", align_from_chains);
        set_annotation(mappings[0], "param_zipcode-tree-score-threshold", (double) zipcode_tree_score_threshold);
        set_annotation(mappings[0], "param_min-to-fragment", (double) min_to_fragment);
        set_annotation(mappings[0], "param_max-to-fragment", (double) max_to_fragment);
        
        // Chaining algorithm parameters
        set_annotation(mappings[0], "param_max-lookback-bases", (double) max_lookback_bases);
        set_annotation(mappings[0], "param_item-bonus", (double) item_bonus);
        set_annotation(mappings[0], "param_item-scale", (double) item_scale);
        set_annotation(mappings[0], "param_max-indel-bases", (double) max_indel_bases);
        
        set_annotation(mappings[0], "param_max-chain-connection", (double) max_chain_connection);
        set_annotation(mappings[0], "param_max-tail-length", (double) max_tail_length);
        set_annotation(mappings[0], "param_max-alignments", (double) max_alignments);
        set_annotation(mappings[0], "param_chain-score", (double) chain_score_threshold);
        set_annotation(mappings[0], "param_min-chain-score-per-base", min_chain_score_per_base);
        set_annotation(mappings[0], "param_max-min-chain-score", (double) max_min_chain_score);
        set_annotation(mappings[0], "param_min-chains", (double) min_chains);
        
    }
    
    // Special fragment and chain statistics
    set_annotation(mappings[0], "fragment_scores", fragment_scores);
    if (track_correctness) {
        set_annotation(mappings[0], "best_chain_correct", best_chain_correct);
    }
    set_annotation(mappings[0], "best_chain_coverage", best_chain_coverage);
    set_annotation(mappings[0], "best_chain_longest_jump", (double) best_chain_longest_jump);
    set_annotation(mappings[0], "best_chain_average_jump", best_chain_average_jump);
    set_annotation(mappings[0], "best_chain_anchors", (double) best_chain_anchors);
    set_annotation(mappings[0], "best_chain_anchor_length", (double) best_chain_anchor_length);
    
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
    const std::vector<size_t>& chain) const {
    
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
    
    // We need a WFAExtender to do tail and intervening alignments.
    // Note that the extender expects anchoring matches!!!
    WFAExtender extender(gbwt_graph, aligner); 
    
    // Keep a couple cursors in the chain: extension before and after the linking up we need to do.
    auto here_it = chain.begin();
    auto next_it = here_it;
    ++next_it;
    
    const algorithms::Anchor* here = &to_chain[*here_it];
    
#ifdef debug_chaining
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
            left_alignment = extender.prefix(left_tail, right_anchor);
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
            
#ifdef debug_chaining
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
                
#ifdef debug
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
            
#ifdef debug_chaining
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
                size_t max_gap_length = this->get_regular_aligner()->longest_detectable_gap(aln, aln.sequence().begin() + left_tail_length);
                size_t graph_horizon = left_tail_length + max_gap_length;

#ifdef warn_on_fallback
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << left_tail_length << " bp left tail against " << right_anchor << " allowing " << max_gap_length << " bp gap in " << aln.name() << endl;
                }
#endif

                // Align the left tail, anchoring the right end.
                align_sequence_between(empty_pos_t(), right_anchor, graph_horizon, max_gap_length, &this->gbwt_graph, this->get_regular_aligner(), tail_aln, &aln.name(), this->max_dp_cells, this->choose_band_padding);
                
                if (show_work) {
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
    }
        
    size_t longest_attempted_connection = 0;
    while(next_it != chain.end()) {
        // Do each region between successive gapless extensions
        
        // We have to find the next item we can actually connect to
        const algorithms::Anchor* next;
        // And the actual connecting alignment to it
        WFAAlignment link_alignment;
        
        while (next_it != chain.end()) {
            next = &to_chain[*next_it];
            // Try and find a next thing to connect to
            
            if (algorithms::get_read_distance(*here, *next) == std::numeric_limits<size_t>::max()) {
                // There's overlap between these items. Keep here and skip next.
#ifdef debug_chaining
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
            
#ifdef debug_chaining
        if (show_work) {
            #pragma omp critical (cerr)
            {
                cerr << log_name() << "Add current item " << *here_it << " of length " << (*here).length() << " with score of " << (*here).score() << endl;
            }
        }
#endif
        
        // Make an alignment for the bases used in this item, and
        // concatenate it in.
        WFAAlignment here_alignment = this->to_wfa_alignment(*here);
        append_path(composed_path, here_alignment.to_path(this->gbwt_graph, aln.sequence()));
        composed_score += here_alignment.score;
        
#ifdef debug_chaining
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
        
#ifdef debug_chaining
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
            
#ifdef debug_chaining
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "Treat as empty link" << endl;
                }
            }
#endif
            
            link_alignment = WFAAlignment::make_empty();
        } else if (link_length > 0 && link_length <= max_chain_connection) {
            // If it's not empty and is a reasonable size, align it.
            // Make sure to walk back the left anchor so it is outside of the region to be aligned.
            pos_t left_anchor = (*here).graph_end();
            get_offset(left_anchor)--;
            
            link_alignment = extender.connect(linking_bases, left_anchor, (*next).graph_start());
            
            longest_attempted_connection = std::max(longest_attempted_connection, linking_bases.size());
            
            if (!link_alignment) {
                // We couldn't align.
                if (graph_length == 0) {
                    // We had read sequence but no graph sequence.
                    // Try falling back to a pure insertion.
                    // TODO: We can be leaving the GBWT's space here!
                    // TODO: What if this is forcing an insertion that could also be in the graph already?
#ifdef debug_chaining
                    if (show_work) {
                        #pragma omp critical (cerr)
                        {
                            cerr << log_name() << "connect() failed; treat as insertion" << endl;
                        }
                    }
#endif
                    link_alignment = WFAAlignment::make_unlocalized_insertion((*here).read_end(), link_length, aligner.score_gap(link_length));
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
            
#ifdef debug_chaining
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
#ifdef debug
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
            size_t max_gap_length = this->get_regular_aligner()->longest_detectable_gap(aln, aln.sequence().begin() + link_start);
            size_t path_length = std::max(graph_length, link_length);
            MinimizerMapper::align_sequence_between((*here).graph_end(), (*next).graph_start(), path_length, max_gap_length, &this->gbwt_graph, this->get_regular_aligner(), link_aln, &aln.name(), this->max_dp_cells, this->choose_band_padding);
            
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
        
        // Advance here to next and start considering the next after it
        here_it = next_it;
        ++next_it;
        here = next;
    }
    
#ifdef debug_chaining
    if (show_work) {
        #pragma omp critical (cerr)
        {
            cerr << log_name() << "Add last extension " << *here_it << " of length " << (*here).length() << " with score of " << (*here).score() << endl;
        }
    }
#endif
    
    WFAAlignment here_alignment = this->to_wfa_alignment(*here);
    
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
        pos_t left_anchor = (*here).graph_end();
        get_offset(left_anchor)--;
        if (right_tail_length <= max_tail_length) {
            // We align the right tail with suffix(), which creates a suffix of the alignment.
            // Make sure to walk back the anchor so it is outside of the region to be aligned.
            right_alignment = extender.suffix(right_tail, left_anchor);
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
                ss << "Aligning right tail " << right_tail << " from " << left_anchor << " produced wrong-length alignment ";
                right_alignment.print(ss);
                throw ChainAlignmentFailedError(ss.str());
            }
#ifdef debug_chaining
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
            
#ifdef debug_chaining
            if (show_work) {
                #pragma omp critical (cerr)
                {
                    cerr << log_name() << "End with long right tail fallback alignment" << endl;
                }
            }
#endif

            if (right_tail.size() > MAX_DP_LENGTH) {
                // Right tail is too long to align.
               
#ifdef debug
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Refusing to align " << right_tail.size() << " bp right tail against " << left_anchor << " in " << aln.name() << " to avoid overflow" << endl;
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
                size_t max_gap_length = this->get_regular_aligner()->longest_detectable_gap(aln, aln.sequence().begin() + (*here).read_end());
                size_t graph_horizon = right_tail_length + max_gap_length;

#ifdef warn_on_fallback
                #pragma omp critical (cerr)
                {
                    cerr << "warning[MinimizerMapper::find_chain_alignment]: Falling back to non-GBWT alignment of " << right_tail_length << " bp right tail against " << left_anchor << " allowing " << max_gap_length << " bp gap in " << aln.name() << endl;
                }
#endif

                // Align the right tail, anchoring the left end.
                align_sequence_between(left_anchor, empty_pos_t(), graph_horizon, max_gap_length, &this->gbwt_graph, this->get_regular_aligner(), tail_aln, &aln.name(), this->max_dp_cells, this->choose_band_padding);
                
                if (show_work) {
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
    *result.mutable_path() = std::move(simplify(composed_path));
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

void MinimizerMapper::align_sequence_between(const pos_t& left_anchor, const pos_t& right_anchor, size_t max_path_length, size_t max_gap_length, const HandleGraph* graph, const GSSWAligner* aligner, Alignment& alignment, const std::string* alignment_name, size_t max_dp_cells, const std::function<size_t(const Alignment&, const HandleGraph&)>& choose_band_padding) {
    
    // Get the dagified local graph, and the back translation
    MinimizerMapper::with_dagified_local_graph(left_anchor, right_anchor, max_path_length, *graph,
        [&](DeletableHandleGraph& dagified_graph, const std::function<std::pair<nid_t, bool>(const handle_t&)>& dagified_handle_to_base) {

#ifdef debug
        std::cerr << "Dagified graph:" << std::endl;
        dump_debug_graph(dagified_graph);
#endif
    
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
                return;
            } else {
#ifdef debug
                #pragma omp critical (cerr)
                std::cerr << "debug[MinimizerMapper::align_sequence_between]: Fill " << cell_count << " DP cells in tail with Xdrop" << std::endl;
#endif
                aligner->align_pinned(alignment, dagified_graph, !is_empty(left_anchor), true, max_gap_length);
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
        if (!is_empty(left_anchor) && alignment.path().mapping_size() > 0 && offset(left_anchor) != 0) {
            // Get the positions of the leftmost mapping
            Position* left_pos = alignment.mutable_path()->mutable_mapping(0)->mutable_position();
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
    if (source.value.is_reverse) {
        // Seed stores the final base of the match in the graph.
        // So get the past-end position.
        pos_t graph_end = make_pos_t(id(seed.pos), is_rev(seed.pos), offset(seed.pos) + 1);
        
        // Work out how much of the node it could use before there.
        length = std::min((size_t) source.length, offset(graph_end));
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

    // Work out how many points the anchor is
    // TODO: Always make sequence and quality available for scoring!
    int score = aligner->score_exact_match(aln, read_start, length);
    return algorithms::Anchor(read_start, graph_start, length, score, seed_number, seed.zipcode_decoder.get(), hint_start); 
}

WFAAlignment MinimizerMapper::to_wfa_alignment(const algorithms::Anchor& anchor) const {
    return {
        {gbwt_graph.get_handle(id(anchor.graph_start()), is_rev(anchor.graph_start()))},
        {{WFAAlignment::match, (uint32_t)anchor.length()}},
        (uint32_t)offset(anchor.graph_start()),
        (uint32_t)anchor.read_start(),
        (uint32_t)anchor.length(),
        anchor.score(),
        true
    };
}

}

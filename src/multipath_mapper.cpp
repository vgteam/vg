/**
 * \file multipath_mapper.cpp
 *
 * Implements the MultipathMapper class
 */

//#define debug_multipath_mapper
//#define debug_multipath_mapper_alignment
//#define debug_validate_multipath_alignments
//#define debug_report_startup_training
//#define debug_pretty_print_alignments
//#define debug_time_phases
//#define debug_log_splice_align_stats
//#define debug_output_distance_correlation
//#define debug_check_adapters

#ifdef debug_time_phases
#include <ctime>
#endif

#ifdef debug_check_adapters
#include "ssw_aligner.hpp"
#endif

#include "multipath_mapper.hpp"

#include "multipath_alignment_graph.hpp"
#include "kmp.hpp"
#include "hash_map.hpp"
#include "position.hpp"
#include "nodeside.hpp"
#include "path.hpp"
#include "utility.hpp"
#include "annotation.hpp"
#include "statistics.hpp"

#include "identity_overlay.hpp"
#include "reverse_graph.hpp"
#include "split_strand_graph.hpp"
#include "dagified_graph.hpp"

#include "algorithms/count_covered.hpp"
#include "algorithms/extract_containing_graph.hpp"
#include "algorithms/locally_expand_graph.hpp"
#include "algorithms/jump_along_path.hpp"
#include "algorithms/ref_path_distance.hpp"
#include "algorithms/component.hpp"

namespace vg {
    
    //size_t MultipathMapper::PRUNE_COUNTER = 0;
    //size_t MultipathMapper::SUBGRAPH_TOTAL = 0;
    //size_t MultipathMapper::SECONDARY_RESCUE_COUNT = 0;
    //size_t MultipathMapper::SECONDARY_RESCUE_ATTEMPT = 0;
    //size_t MultipathMapper::SECONDARY_RESCUE_TOTAL = 0;
    
    const size_t MultipathMapper::RESCUED = numeric_limits<size_t>::max();

    MultipathMapper::MultipathMapper(PathPositionHandleGraph* graph, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_array,
                                     haplo::ScoreProvider* haplo_score_provider, SnarlManager* snarl_manager,
                                     SnarlDistanceIndex* distance_index) :
        BaseMapper(graph, gcsa_index, lcp_array, haplo_score_provider),
        snarl_manager(snarl_manager),
        distance_index(distance_index),
        path_component_index(distance_index ? nullptr : new PathComponentIndex(graph)),
        splice_stats(*get_regular_aligner())
    {
        set_max_merge_supression_length();
    }

    MultipathMapper::~MultipathMapper() {
        
    }
    
    void MultipathMapper::multipath_map(const Alignment& alignment,
                                        vector<multipath_alignment_t>& multipath_alns_out) {
        
#ifdef debug_multipath_mapper
        cerr << "multipath mapping read " << pb2json(alignment) << endl;
        cerr << "querying MEMs..." << endl;
#endif
        
        vector<deque<pair<string::const_iterator, char>>> mem_fanouts;
        auto mems = find_mems(alignment, &mem_fanouts);
        unique_ptr<match_fanouts_t> fanouts(mem_fanouts.empty() ? nullptr :
                                            new match_fanouts_t(record_fanouts(mems, mem_fanouts)));
        
#ifdef debug_multipath_mapper
        cerr << "obtained MEMs:" << endl;
        for (MaximalExactMatch mem : mems) {
            cerr << "\t" << mem << " (" << mem.nodes.size() << " hits)" << endl;
            if (fanouts.get() && fanouts->count(&mem)) {
                cerr << "\t\tfan-outs:" << endl;
                for (auto fanout : fanouts->at(&mem)) {
                    cerr << "\t\t\t" << (fanout.first - mem.begin) << ": " << *fanout.first << " -> " << fanout.second << endl;
                }
            }
        }
        cerr << "clustering MEMs..." << endl;
#endif
        
        // TODO: use the automatic expected MEM length algorithm to restrict the MEMs used for clustering?
        
        // cluster the MEMs
        MemoizingGraph memoizing_graph(xindex);
        unique_ptr<OrientedDistanceMeasurer> distance_measurer = get_distance_measurer(memoizing_graph);
        
        vector<memcluster_t> clusters = get_clusters(alignment, mems, &(*distance_measurer), fanouts.get());
        
#ifdef debug_multipath_mapper
        cerr << "obtained clusters:" << endl;
        for (int i = 0; i < clusters.size(); i++) {
            cerr << "\tcluster " << i << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : clusters[i].first) {
                cerr << "\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
                if (fanouts.get() && fanouts->count(hit.first)) {
                    cerr << "\t\t\tfan-outs:" << endl;
                    for (auto fanout : fanouts->at(hit.first)) {
                        cerr << "\t\t\t\t" << (fanout.first - hit.first->begin) << ": " << *fanout.first << " -> " << fanout.second << endl;
                    }
                }
            }
        }
        cerr << "extracting subgraphs..." << endl;
#endif
        
        // extract graphs around the clusters
        auto cluster_graphs = query_cluster_graphs(alignment, mems, clusters);
        
        // actually perform the alignments and post-process to meet multipath_alignment_t invariants
        // TODO: do i still need cluster_idx? i think it might have only been used for capping
        vector<double> multiplicities;
        vector<size_t> cluster_idxs;
        align_to_cluster_graphs(alignment, cluster_graphs, multipath_alns_out, multiplicities,
                                num_mapping_attempts, fanouts.get(), &cluster_idxs);
        
        if (multipath_alns_out.empty()) {
            // add a null alignment so we know it wasn't mapped
            multipath_alns_out.emplace_back();
            cluster_graphs.emplace_back();
            multiplicities.push_back(1.0);
            to_multipath_alignment(alignment, multipath_alns_out.back());
            
            // in case we're realigning GAMs that have paths already
            multipath_alns_out.back().clear_subpath();
            multipath_alns_out.back().clear_start();
        }
        
        if (do_spliced_alignment && !likely_mismapping(multipath_alns_out.front())) {
            find_spliced_alignments(alignment, multipath_alns_out, multiplicities, cluster_idxs,
                                    mems, cluster_graphs, fanouts.get());
        }
        
        if (agglomerate_multipath_alns) {
            // we want the disconnected alignments combined into one
            agglomerate_alignments(multipath_alns_out, &multiplicities);
        }
        
        // if we computed extra alignments to get a mapping quality, remove them
        if (multipath_alns_out.size() > max_alt_mappings) {
            multipath_alns_out.resize(max_alt_mappings);
        }
        
        // mark unmapped reads and get rid of noise alignments
        purge_unmapped_alignments(multipath_alns_out);
        
        for (size_t i = 1; i < multipath_alns_out.size(); ++i) {
            multipath_alns_out[i].set_annotation("secondary", true);
        }
        
        if (strip_bonuses) {
            for (multipath_alignment_t& multipath_aln : multipath_alns_out) {
                strip_full_length_bonuses(multipath_aln);
            }
        }
        
#ifdef debug_pretty_print_alignments
        cerr << "final alignments being returned:" << endl;
        for (const multipath_alignment_t& multipath_aln : multipath_alns_out) {
            view_multipath_alignment(cerr, multipath_aln, *xindex);
        }
#endif
        
#ifdef mpmap_instrument_mem_statistics
        size_t num_mems = mems.size();
        size_t min_mem_length = numeric_limits<size_t>::max();
        size_t max_mem_length = 0;
        double avg_mem_length = 0.0;
        for (const auto& mem : mems) {
            min_mem_length = min<size_t>(min_mem_length, mem.length());
            max_mem_length = max<size_t>(max_mem_length, mem.length());
            avg_mem_length += mem.length();
        }
        avg_mem_length /= mems.size();
        double avg_mem_overlap = 0.0;
        for (size_t i = 1; i < mems.size(); ++i) {
            avg_mem_overlap += max<int64_t>(mems[i - 1].end - mems[i].begin, 0);
        }
        avg_mem_overlap /= (mems.size() - 1);
        
        vector<size_t> hit_lengths;
        for (const auto& mem : mems) {
            for (const auto& n : mem.nodes) {
                hit_lengths.push_back(mem.length());
            }
        }
        sort(hit_lengths.begin(), hit_lengths.end(), std::greater<size_t>());
        
        size_t num_clusters = clusters.size();
        vector<size_t> winning_lengths;
        size_t winning_cluster_num_mems = clusters.empty() ? 0 : clusters[cluster_idxs.front()].size();
        size_t winning_cluster_total_bases = 0;
        size_t winning_cluster_min_mem_length = numeric_limits<size_t>::max();
        size_t winning_cluster_max_mem_length = 0;
        size_t winning_cluster_tail_bases = 0;
        double winning_cluster_avg_intermem_gap = 0.0;
        vector<size_t> order;

        if (!clusters.empty()) {
            for (const auto& hit : clusters[cluster_idxs.front()]) {
                winning_cluster_min_mem_length = min<size_t>(winning_cluster_min_mem_length, hit.first->length());
                winning_cluster_max_mem_length = max<size_t>(winning_cluster_max_mem_length, hit.first->length());
                winning_cluster_total_bases += hit.first->length();
            }
            for (size_t i = 0; i < clusters[cluster_idxs.front()].size(); ++i) {
                order.push_back(i);
                winning_lengths.push_back(clusters[cluster_idxs.front()][i].first->length());
            }
            sort(order.begin(), order.end(), [&](size_t i, size_t j) {
                return clusters[cluster_idxs.front()][i].first->begin < clusters[cluster_idxs.front()][j].first->begin;
            });
            sort(winning_lengths.begin(), winning_lengths.end(), std::greater<size_t>());
            
            winning_cluster_tail_bases = ((clusters[cluster_idxs.front()][order.front()].first->begin - alignment.sequence().begin())
                                          + (alignment.sequence().end() - clusters[cluster_idxs.front()][order.back()].first->end));
            if (clusters[cluster_idxs.front()].size() == 0) {
                winning_cluster_avg_intermem_gap = numeric_limits<double>::quiet_NaN();
            }
            else {
                for (size_t i = 1; i < order.size(); ++i) {
                    winning_cluster_avg_intermem_gap += (clusters[cluster_idxs.front()][order[i]].first->begin
                                                         - clusters[cluster_idxs.front()][order[i - 1]].first->end);
                }
                winning_cluster_avg_intermem_gap /= order.size() - 1;
            }
        }
        
        vector<size_t> secondary_lengths;
        if (cluster_idxs.size() > 1 && clusters.size() > 1) {
            for (const auto& hit : clusters.at(cluster_idxs[1])) {
                secondary_lengths.push_back(hit.first->length());
            }
        }
        sort(secondary_lengths.begin(), secondary_lengths.end(), greater<size_t>());
        
        vector<size_t> secondary_lengths;
        if (cluster_idxs.size() > 1 && clusters.size() > 1) {
            for (const auto& hit : clusters[cluster_idxs[1]]) {
                secondary_lengths.push_back(hit.first->length());
            }
        }
        sort(secondary_lengths.begin(), secondary_lengths.end(), greater<size_t>());
        
        int64_t max_non_winning_mem_length = 0;
        for (size_t i = 0; i < mems.size(); ++i) {
            bool found = false;
            if (!clusters.empty()) {
                for (const auto hit : clusters[cluster_idxs.front()]) {
                    if (hit.first == &mems[i]) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                max_non_winning_mem_length = max<int64_t>(max_non_winning_mem_length, mems[i].length());
            }
        }
        
#pragma omp critical
        {
            if (!_wrote_mem_stats_header) {
                _mem_stats << "name\tread_len\tnum_mems\tmin_mem_length\tmax_mem_length\tavg_mem_length\tavg_mem_overlap\tnum_clusters\twinning_cluster_num_mems\twinning_cluster_min_mem_length\twinning_cluster_max_mem_length\twinning_cluster_total_bases\twinning_cluster_tail_bases\twinning_cluster_avg_intermem_gap\tmax_non_winning_mem_length\tmapping_quality\thit_lengths\twinning_lengths\tsecondary_lengths" << endl;
                _wrote_mem_stats_header = true;
            }
            _mem_stats << alignment.name() << "\t" << alignment.sequence().size() << "\t" << num_mems << "\t" << min_mem_length << "\t" << max_mem_length << "\t" << avg_mem_length << "\t" << avg_mem_overlap << "\t" << num_clusters << "\t" << winning_cluster_num_mems << "\t" << winning_cluster_min_mem_length << "\t" << winning_cluster_max_mem_length << "\t" << winning_cluster_total_bases << "\t" << winning_cluster_tail_bases << "\t" << winning_cluster_avg_intermem_gap << "\t" << max_non_winning_mem_length << "\t" << multipath_alns_out.front().mapping_quality();
            _mem_stats << "\t";
            for (size_t i = 0; i < hit_lengths.size(); ++i) {
                if (i > 0) {
                    _mem_stats << ",";
                }
                _mem_stats << hit_lengths[i];
            }
            _mem_stats << "\t";
            for (size_t i = 0; i < winning_lengths.size(); ++i) {
                if (i > 0) {
                    _mem_stats << ",";
                }
                _mem_stats << winning_lengths[i];
            }
            _mem_stats << "\t";
            if (secondary_lengths.empty()) {
                _mem_stats << "NA";
            }
            else {
                for (size_t i = 0; i < secondary_lengths.size(); ++i) {
                    if (i > 0) {
                        _mem_stats << ",";
                    }
                    _mem_stats << secondary_lengths[i];
                }
            }
            _mem_stats << endl;
        }
#endif
    }
    
    vector<MultipathMapper::memcluster_t> MultipathMapper::get_clusters(const Alignment& alignment, const vector<MaximalExactMatch>& mems,
                                                                        OrientedDistanceMeasurer* distance_measurer,
                                                                        const match_fanouts_t* fanouts) const {
        
        // note: we don't want to generate the distance measurer in this function because we want
        // to be able to re-use its memoization if we cluster pairs later
        
        // choose a clusterer (ordered by expected most likely use for better branch prediction)
        unique_ptr<MEMClusterer> clusterer;
        if (!no_clustering && use_min_dist_clusterer && component_min_dist) {
            clusterer = unique_ptr<MEMClusterer>(new ComponentMinDistanceClusterer(distance_index));
        }
        else if (!no_clustering && !use_min_dist_clusterer && !use_tvs_clusterer) {
            clusterer = unique_ptr<MEMClusterer>(new OrientedDistanceClusterer(*distance_measurer,
                                                                               max_expected_dist_approx_error));
        }
        else if (no_clustering) {
            clusterer = unique_ptr<MEMClusterer>(new NullClusterer());
        }
        else if (use_min_dist_clusterer && !greedy_min_dist) {
            clusterer = unique_ptr<MEMClusterer>(new MinDistanceClusterer(distance_index));
        }
        else if (use_min_dist_clusterer && greedy_min_dist) {
            clusterer = unique_ptr<MEMClusterer>(new GreedyMinDistanceClusterer(distance_index));
        }
        else {
            clusterer = unique_ptr<MEMClusterer>(new TVSClusterer(xindex, distance_index));
        }
        clusterer->max_gap = max_alignment_gap;
        
        // generate clusters
        return clusterer->clusters(alignment, mems, get_aligner(!alignment.quality().empty()),
                                   min_clustering_mem_length, max_mapping_quality, log_likelihood_approx_factor,
                                   min_median_mem_coverage_for_split, 0.75, truncation_multiplicity_mq_limit,
                                   fanouts);
    }
    
    vector<MaximalExactMatch> MultipathMapper::find_mems(const Alignment& alignment,
                                                         vector<deque<pair<string::const_iterator, char>>>* mem_fanout_breaks) {
        if (!use_stripped_match_alg &&
            (!use_fanout_match_alg || (use_fanout_match_alg && alignment.quality().empty()))) {
            double dummy1, dummy2;
            return find_mems_deep(alignment.sequence().begin(), alignment.sequence().end(), dummy1, dummy2,
                                  0, min_mem_length, mem_reseed_length, false, true, true, false);
        }
        else if (use_fanout_match_alg) {
            return find_fanout_mems(alignment.sequence().begin(), alignment.sequence().end(),
                                    alignment.quality().begin(), max_fans_out, max_fanout_base_quality,
                                    mem_fanout_breaks);
        }
        else {
            return find_stripped_matches(alignment.sequence().begin(), alignment.sequence().end(),
                                         stripped_match_alg_strip_length, stripped_match_alg_max_length,
                                         stripped_match_alg_target_count);
        }
    }

    vector<pair<pair<size_t, size_t>, int64_t>> MultipathMapper::get_cluster_pairs(const Alignment& alignment1,
                                                                                   const Alignment& alignment2,
                                                                                   vector<clustergraph_t>& cluster_graphs1,
                                                                                   vector<clustergraph_t>& cluster_graphs2,
                                                                                   OrientedDistanceMeasurer* distance_measurer) {
        // make vectors of cluster pointers to shim into the cluster pairing function
        vector<memcluster_t*> cluster_mems_1(cluster_graphs1.size()), cluster_mems_2(cluster_graphs2.size());
        for (size_t i = 0; i < cluster_mems_1.size(); i++) {
            cluster_mems_1[i] = &(get<1>(cluster_graphs1[i]));
        }
        for (size_t i = 0; i < cluster_mems_2.size(); i++) {
            cluster_mems_2[i] = &(get<1>(cluster_graphs2[i]));
        }
        
        // Find the clusters that have a tie for the longest MEM, and create alternate anchor points for those clusters
        vector<pair<size_t, size_t>> alt_anchors_1, alt_anchors_2;
        for (size_t i = 0; i < cluster_mems_1.size(); i++) {
            auto& mem_cluster = cluster_mems_1[i]->first;
            for (size_t j = 1; j < mem_cluster.size(); j++) {
                if (mem_cluster[j].first->length() + alt_anchor_max_length_diff >= mem_cluster.front().first->length()) {
                    alt_anchors_1.emplace_back(i, j);
                }
                else {
                    break;
                }
            }
        }
        for (size_t i = 0; i < cluster_mems_2.size(); i++) {
            auto& mem_cluster = cluster_mems_2[i]->first;
            for (size_t j = 1; j < mem_cluster.size(); j++) {
                if (mem_cluster[j].first->length() + alt_anchor_max_length_diff >= mem_cluster.front().first->length()) {
                    alt_anchors_2.emplace_back(i, j);
                }
                else {
                    break;
                }
            }
        }
        
        // Compute the pairs of cluster graphs and their approximate distances from each other
        // (ordered by expected most likely use case for better branch prediction)
        unique_ptr<MEMClusterer> clusterer;
        if (use_min_dist_clusterer && !no_clustering) {
            // greedy and non-greedy algorithms are the same, so don't bother distinguishing
            clusterer = unique_ptr<MEMClusterer>(new MinDistanceClusterer(distance_index));
        }
        else if (no_clustering) {
            clusterer = unique_ptr<MEMClusterer>(new NullClusterer());
        }
        else if (!use_tvs_clusterer) {
            clusterer = unique_ptr<MEMClusterer>(new OrientedDistanceClusterer(*distance_measurer));
        }
        else {
            clusterer = unique_ptr<MEMClusterer>(new TVSClusterer(xindex, distance_index));
        }
        
        return clusterer->pair_clusters(alignment1, alignment2, cluster_mems_1, cluster_mems_2,
                                        alt_anchors_1, alt_anchors_2,
                                        fragment_length_distr.mean(),
                                        ceil(10.0 * fragment_length_distr.std_dev()));
    }
    
    void MultipathMapper::align_to_cluster_graphs(const Alignment& alignment,
                                                  vector<clustergraph_t>& cluster_graphs,
                                                  vector<multipath_alignment_t>& multipath_alns_out,
                                                  vector<double>& multiplicities_out,
                                                  size_t num_mapping_attempts,
                                                  const match_fanouts_t* fanouts,
                                                  vector<size_t>* cluster_idxs) {
        
        
#ifdef debug_multipath_mapper
        cerr << "aligning to (up to) " << cluster_graphs.size() << " subgraphs..." << endl;
#endif
      
        // we may need to compute an extra mapping above the one we'll report if we're computing mapping quality
        size_t num_mappings_to_compute = max(num_mapping_attempts, (size_t) 2);
        
        multipath_alns_out.clear();
        
//#pragma omp atomic
//        SUBGRAPH_TOTAL += cluster_graphs.size();
        
        // align to each cluster subgraph
        size_t num_mappings = 0;
        for (auto& cluster_graph : cluster_graphs) {
            // if we have a cluster graph with small enough MEM coverage compared to the best one or we've made
            // the maximum number of alignments we stop producing alternate alignments
            if (get<2>(cluster_graph) < mem_coverage_min_ratio * get<2>(cluster_graphs.front())
                || num_mappings >= num_mappings_to_compute) {
#ifdef debug_multipath_mapper
                cerr << "halting further alignments, either because MEM coverage of " << get<2>(cluster_graph) << " is too far below optimum of " << get<2>(cluster_graphs.front()) << " or because already made " << num_mappings << " of " << num_mappings_to_compute << " mappings" << endl;
#endif
                
//#pragma omp atomic
//                PRUNE_COUNTER += cluster_graphs.size() - num_mappings;
                break;
            }
            
#ifdef debug_multipath_mapper_alignment
            cerr << "performing alignment to subgraph with coverage " << get<2>(cluster_graph) << " and multiplicity " << get<1>(cluster_graph).second << endl;
#endif
            double multiplicity = cluster_multiplicity(get<1>(cluster_graph));
            size_t cluster_size = get<1>(cluster_graph).first.size();
            multipath_alns_out.emplace_back();
            multipath_align(alignment, cluster_graph, multipath_alns_out.back(), fanouts);
            multiplicities_out.emplace_back(multiplicity);
            num_mappings++;
        }
        
        if (!multipath_alns_out.empty()) {
            // find clusters whose likelihoods are approximately the same as the low end of the clusters we aligned
            auto aligner = get_aligner(!alignment.quality().empty());
            int64_t score_diff = round(aligner->mapping_quality_score_diff(truncation_multiplicity_mq_limit));
            int64_t max_tail_idx = multipath_alns_out.size();
            while (max_tail_idx < cluster_graphs.size()
                   && get<2>(cluster_graphs[max_tail_idx]) >= get<2>(cluster_graphs[multipath_alns_out.size() - 1]) - score_diff) {
                ++max_tail_idx;
            }
            if (max_tail_idx > multipath_alns_out.size()) {
                // there are some (nearly) identical cluster that we ignored, so we'll account for them in the multiplicity
                
                // find the clusters that are approximately the same
                int64_t min_tail_idx = multipath_alns_out.size() - 1;
                while (min_tail_idx > 0 &&
                       get<2>(cluster_graphs[min_tail_idx - 1]) <= get<2>(cluster_graphs[multipath_alns_out.size() - 1]) + score_diff) {
                    --min_tail_idx;
                }
                
                // multiply their multiplicity by the inverse of the fraction aligned
                double trunc_multiplicity = double(max_tail_idx - min_tail_idx) / double(multipath_alns_out.size() - min_tail_idx);
                for (size_t i = min_tail_idx; i < multipath_alns_out.size(); ++i) {
                    multiplicities_out[i] *= trunc_multiplicity;
                }
            }
        }
        
        if (cluster_idxs) {
            // initially all of the multipath alignments are in the order of the clusters
            *cluster_idxs = range_vector(multipath_alns_out.size());
        }
        
#ifdef debug_multipath_mapper
        cerr << "splitting multicomponent alignments..." << endl;
#endif
        
        if (!suppress_multicomponent_splitting) {
            // split up any alignments that ended up being disconnected
            split_multicomponent_alignments(multipath_alns_out, &alignment, &cluster_graphs, cluster_idxs, &multiplicities_out);
        }
        
        // remove low-scoring bits of complicated multipath alignments
        for (auto& multipath_aln : multipath_alns_out) {
            topologically_order_subpaths(multipath_aln);
            simplify_complicated_multipath_alignment(multipath_aln);
        }
        
#ifdef debug_multipath_mapper
        cerr << "computing mapping quality and sorting mappings" << endl;
#endif
        sort_and_compute_mapping_quality(multipath_alns_out, cluster_idxs, &multiplicities_out);
        
        if (!multipath_alns_out.empty() && likely_mismapping(multipath_alns_out.front())) {
            multipath_alns_out.front().set_mapping_quality(0);
        }
        
        // for debugging: an expensive check for invariant validity that can be turned on
        // with a preprocessor flag
#ifdef debug_validate_multipath_alignments
        for (multipath_alignment_t& multipath_aln : multipath_alns_out) {
#ifdef debug_multipath_mapper
            cerr << "validating multipath alignment:" << endl;
            cerr << debug_string(multipath_aln) << endl;
#endif
            if (!validate_multipath_alignment(multipath_aln, *xindex)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln.sequence() << " failed to validate" << endl;
            }
        }
#endif
        
    }
    
    bool MultipathMapper::attempt_unpaired_multipath_map_of_pair(const Alignment& alignment1, const Alignment& alignment2,
                                                                 vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                                 vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer) {
        
        // compute single ended mappings, and make sure we also compute mapping qualities to assess
        // mapping ambiguity
        vector<multipath_alignment_t> multipath_alns_1, multipath_alns_2;
        multipath_map(alignment1, multipath_alns_1);
        multipath_map(alignment2, multipath_alns_2);
        
        bool is_ambiguous = true;
        
        if (!multipath_alns_1.empty() && !multipath_alns_2.empty()) {
            multipath_alignment_t& multipath_aln_1 = multipath_alns_1.front();
            multipath_alignment_t& multipath_aln_2 = multipath_alns_2.front();
            
            auto aligner = get_aligner(!alignment1.quality().empty() && !alignment2.quality().empty());
            
            // score possible of a perfect match (at full base quality)
            int32_t max_score_1 = multipath_aln_1.sequence().size() * aligner->match + 2 * aligner->full_length_bonus * !strip_bonuses;
            int32_t max_score_2 = multipath_aln_2.sequence().size() * aligner->match + 2 * aligner->full_length_bonus * !strip_bonuses;
            
#ifdef debug_multipath_mapper
            cerr << "single ended mappings achieves scores " << optimal_alignment_score(multipath_aln_1) << " and " << optimal_alignment_score(multipath_aln_2) << ", looking for scores " << .8 * max_score_1 << " and " << .8 * max_score_2 << endl;
            cerr << "single ended mappings achieves mapping qualities " << multipath_aln_1.mapping_quality() << " and " << multipath_aln_2.mapping_quality() << ", looking for mapq " << min(max_mapping_quality, 45) << endl;
#endif
            
            // are these reads unambiguously mapped and well-aligned?
            // TODO: i don't like having constants floating around in here
            if (num_connected_components(multipath_aln_1) == 1
                && num_connected_components(multipath_aln_2) == 1
                && multipath_aln_1.mapping_quality() >= min(max_mapping_quality, 45)
                && multipath_aln_2.mapping_quality() >= min(max_mapping_quality, 45)
                && optimal_alignment_score(multipath_aln_1) >= .8 * max_score_1
                && optimal_alignment_score(multipath_aln_2) >= .8 * max_score_2) {
                
                int64_t fragment_length = distance_between(multipath_aln_1, multipath_aln_2, true);
                
#ifdef debug_multipath_mapper
                cerr << "fragment length between mappings measured at " << fragment_length << endl;
#endif
                
                // can we obtain a distance between these positions?
                if (fragment_length != numeric_limits<int64_t>::max()) {
                    
                    // record the unambiguous mappings and the fragment length
                    
                    
#ifdef debug_multipath_mapper
                    cerr << "registering measurement, now have " << fragment_length_distr.curr_sample_size() << " of " << fragment_length_distr.max_sample_size() << endl;
#endif
                    
                    multipath_aln_pairs_out.emplace_back(move(multipath_aln_1), move(multipath_aln_2));
                    fragment_length_distr.register_fragment_length(fragment_length);
                    
                    is_ambiguous = false;
                }
            }
        }
        
        if (is_ambiguous) {
            // we didn't find an unambiguous pairing in single-ended mode, buffer these for once
            // the paired mode is finalized
#ifdef debug_multipath_mapper
            cerr << "couldn't find unambiguous mapping, adding pair to ambiguous buffer" << endl;
#endif
            
            ambiguous_pair_buffer.emplace_back(alignment1, alignment2);
            
            if (ambiguous_pair_buffer.size() + fragment_length_distr.curr_sample_size()
                == fragment_length_distr.max_sample_size() * fragment_length_warning_factor) {
                cerr << "warning:[vg mpmap] Mapped " << ambiguous_pair_buffer.size() + fragment_length_distr.curr_sample_size() << " read pairs as unpaired reads to learn fragment length distribution, but only obtained " << fragment_length_distr.curr_sample_size() << " unambiguous, consistently mapped pairs. Often this indicates data issues, such as reads that are pre-sorted with unmappable reads at the front, reads that are not actually paired, or mismatched indexes." << endl;
            }
        }
        
        
        // for debugging:
        // we must have just finalized the distribution or else we wouldn't have entered this function
#ifdef debug_report_startup_training
        if (fragment_length_distr.is_finalized()) {
            cerr << "finalized read distribution with " << fragment_length_distr.max_sample_size() << " measurements on read pair " << alignment1.name() << ", " << alignment2.name() << endl;
            cerr << "mean: " << fragment_length_distr.mean() << endl;
            cerr << "std dev: " << fragment_length_distr.std_dev() << endl;
            cerr << "ambiguous buffer contains pairs:" << endl;
            for (pair<Alignment,Alignment>& aln_pair : ambiguous_pair_buffer) {
                cerr << "\t" << aln_pair.first.name() << ", " << aln_pair.second.name() << endl;
            }
            cerr << "distance measurements:" << endl;
            auto iter = fragment_length_distr.measurements_begin();
            if (iter != fragment_length_distr.measurements_end()) {
                cerr << *iter;
                iter++;
            }
            for (; iter != fragment_length_distr.measurements_end(); iter++) {
                cerr << ", " << *iter;
            }
            cerr << endl;
        }
#endif
        
        return !is_ambiguous;
    }
    
    bool MultipathMapper::attempt_rescue(const multipath_alignment_t& multipath_aln, const Alignment& other_aln,
                                         bool rescue_forward, multipath_alignment_t& rescue_multipath_aln) {
        
#ifdef debug_multipath_mapper
        cerr << "attemping pair rescue in " << (rescue_forward ? "forward" : "backward") << " direction from " << debug_string(multipath_aln) << endl;
#endif
        bool succeeded = do_rescue_alignment(multipath_aln, other_aln, rescue_forward, rescue_multipath_aln,
                                             fragment_length_distr.mean(), rescue_graph_std_devs);
        
        if (!succeeded) {
            return false;
        }
        
        auto aligner = get_aligner(!multipath_aln.quality().empty() && !other_aln.quality().empty());
        vector<double> score(1, optimal_alignment_score(rescue_multipath_aln));
        int32_t solo_mapq = mapq_scaling_factor * aligner->compute_max_mapping_quality(score, mapping_quality_method == Approx);
        int32_t adjusted_mapq = min<int32_t>(solo_mapq, min(max_mapping_quality, multipath_aln.mapping_quality()));
        rescue_multipath_aln.set_mapping_quality(adjusted_mapq);
        
#ifdef debug_multipath_mapper
        cerr << "converted multipath alignment is" << endl;
        cerr << debug_string(rescue_multipath_aln) << endl;
        cerr << "rescued alignment has effective match length " << pseudo_length(rescue_multipath_aln) << ", which gives p-value " << random_match_p_value(pseudo_length(rescue_multipath_aln), rescue_multipath_aln.sequence().size()) << endl;
#endif

        // TODO: magic number
        if (solo_mapq < min(25, max_mapping_quality)) {
#ifdef debug_multipath_mapper
            cerr << "rescue fails because raw_mapq " << solo_mapq << " < " << min(25, max_mapping_quality) << endl;
#endif
            return false;
        }
        
        if (likely_misrescue(rescue_multipath_aln)) {
#ifdef debug_multipath_mapper
            cerr << "rescue fails with p value above " << max_rescue_p_value << endl;
#endif
            return false;
        }
        
        return true;
    }

    bool MultipathMapper::do_rescue_alignment(const multipath_alignment_t& multipath_aln, const Alignment& other_aln,
                                              bool rescue_forward, multipath_alignment_t& rescue_multipath_aln,
                                              double rescue_mean_length, double num_std_devs) const {
        bdsg::HashGraph rescue_graph;
        extract_rescue_graph(multipath_aln, other_aln, rescue_forward, &rescue_graph,
                             rescue_mean_length, rescue_graph_std_devs);
        
        if (rescue_graph.get_node_count() == 0) {
            return false;
        }
        
#ifdef debug_multipath_mapper_alignment
        cerr << "got rescue graph" << endl;
        rescue_graph.for_each_handle([&](const handle_t& h) {
            cerr << rescue_graph.get_id(h) << " " << rescue_graph.get_sequence(h) << endl;
            rescue_graph.follow_edges(h, true, [&](const handle_t& p) {
                cerr << "\t" << rescue_graph.get_id(p) << (rescue_graph.get_is_reverse(p) ? "-" : "+") << " <-" << endl;
            });
            rescue_graph.follow_edges(h, false, [&](const handle_t& n) {
                cerr << "\t-> " << rescue_graph.get_id(n) << (rescue_graph.get_is_reverse(n) ? "-" : "+") << endl;
            });
        });
#endif
        
        // TODO: repetitive code with multipath_align
        
        // the longest path we could possibly align to (full gap and a full sequence)
        auto aligner = get_aligner(!multipath_aln.quality().empty() && !other_aln.quality().empty());
        size_t target_length = other_aln.sequence().size() + min(aligner->longest_detectable_gap(other_aln), max_alignment_gap);
        
        // convert from bidirected to directed
        StrandSplitGraph align_digraph(&rescue_graph);
        
        // if necessary, convert from cyclic to acylic (this is expensive, so only do it if we need to)
        IdentityOverlay undagified(&align_digraph);
        unique_ptr<DagifiedGraph> dagified;
        
        ExpandingOverlayGraph* align_dag = nullptr;
        if (handlealgs::is_directed_acyclic(&align_digraph)) {
            align_dag = &undagified;
        }
        else {
#ifdef debug_multipath_mapper_alignment
            cerr << "graph contains directed cycles, performing dagification" << endl;
#endif
            dagified = unique_ptr<DagifiedGraph>(new DagifiedGraph(&align_digraph, target_length, max_dagify_duplications));
            align_dag = dagified.get();
        }
        
        // put local alignment here
        Alignment aln = other_aln;
        // in case we're realigning a GAM, get rid of the path and score
        aln.clear_path();
        aln.set_score(0);
        
        aligner->align(aln, *align_dag, true);
        
        // get the IDs back into the space of the reference graph
        function<pair<id_t, bool>(id_t)> translator = [&](const id_t& node_id) {
            handle_t original = align_digraph.get_underlying_handle(align_dag->get_underlying_handle(align_dag->get_handle(node_id)));
            return make_pair(rescue_graph.get_id(original), rescue_graph.get_is_reverse(original));
        };
        translate_oriented_node_ids(*aln.mutable_path(), translator);
        
#ifdef debug_multipath_mapper
        cerr << "rescued direct alignment is" << endl;
        cerr << pb2json(aln) << endl;
#endif
        
        if (num_alt_alns > 1 && (snarl_manager != nullptr || distance_index != nullptr)) {
            // make an interesting multipath alignment by realigning the single path alignment inside snarls
            make_nontrivial_multipath_alignment(aln, *align_dag, translator, rescue_multipath_aln);
            
        }
        else {
            // just convert the single alignment into a trivial multipath alignment
            to_multipath_alignment(aln, rescue_multipath_aln);
        }
        
        identify_start_subpaths(rescue_multipath_aln);
        
        return true;
    }

    void MultipathMapper::extract_rescue_graph(const multipath_alignment_t& multipath_aln, const Alignment& other_aln,
                                               bool rescue_forward, MutableHandleGraph* rescue_graph,
                                               double rescue_mean_length, double num_std_devs) const {
        
        // get the position to jump from and the distance to jump
        Alignment opt_anchoring_aln;
        optimal_alignment(multipath_aln, opt_anchoring_aln);
        
        if (opt_anchoring_aln.path().mapping_size() == 0) {
            return;
        }
                
        if (get_rescue_graph_from_paths || !distance_index) {
            // we're either not using the distance index or we don't have one
            pos_t pos_from = rescue_forward ? initial_position(opt_anchoring_aln.path()) : final_position(opt_anchoring_aln.path());
            int64_t jump_dist = rescue_forward ? rescue_mean_length : -rescue_mean_length;
            
            // get the seed position(s) for the rescue by jumping along paths
            vector<pos_t> jump_positions = algorithms::jump_along_closest_path(xindex, pos_from, jump_dist, 250);
            
#ifdef debug_multipath_mapper
            cerr << "found jump positions:" << endl;
            for (pos_t& pos : jump_positions) {
                cerr << "\t" << pos << endl;
            }
#endif
            if (jump_positions.empty()) {
                return;
            }
            
            size_t search_dist_bwd, search_dist_fwd;
            if (rescue_forward) {
                search_dist_bwd = size_t(round(num_std_devs * fragment_length_distr.std_dev())) + other_aln.sequence().size();
                search_dist_fwd = num_std_devs * fragment_length_distr.std_dev();
            }
            else {
                search_dist_bwd = num_std_devs * fragment_length_distr.std_dev();
                search_dist_fwd = size_t(round(num_std_devs * fragment_length_distr.std_dev())) + other_aln.sequence().size();
            }
            
            vector<size_t> backward_dist(jump_positions.size(), search_dist_bwd);
            vector<size_t> forward_dist(jump_positions.size(), search_dist_fwd);
            algorithms::extract_containing_graph(xindex, rescue_graph, jump_positions, backward_dist, forward_dist,
                                                 num_alt_alns > 1 ? reversing_walk_length : 0);
            
        }
        else {
            // we have a distance index and we want to use it
            
            // get the set of nodes that we want to extrat
            unordered_set<id_t> subgraph_nodes_to_add;
            int64_t min_distance = max(0.0, rescue_mean_length - other_aln.sequence().size()
                                       - rescue_graph_std_devs * fragment_length_distr.std_dev());
            int64_t max_distance = rescue_mean_length + rescue_graph_std_devs * fragment_length_distr.std_dev();
            subgraph_in_distance_range(*distance_index, opt_anchoring_aln.path(), xindex, min_distance, max_distance,
                                              subgraph_nodes_to_add, rescue_forward);
            
            // this algorithm is better matched to the GBWTGraph, we need to extract the subgraph manually now.
            // we'll use an algorithm that tries to follow edges to find nodes. this way we minimize calls to XG's
            // get_handle, and we have to follow all edges anyway to add them
            
            while (!subgraph_nodes_to_add.empty()) {
                // there's at least one node that we haven't added yet
                
                // initialize a search out from an arbitrary unadded node
                id_t node_id = *subgraph_nodes_to_add.begin();
                subgraph_nodes_to_add.erase(node_id);
                handle_t start_handle = xindex->get_handle(node_id);
                rescue_graph->create_handle(xindex->get_sequence(start_handle),
                                            xindex->get_id(start_handle));
                vector<handle_t> stack(1, start_handle);
                while (!stack.empty()) {
                    handle_t super_handle = stack.back();
                    stack.pop_back();
                    for (bool go_left : {true, false}) {
                        xindex->follow_edges(super_handle, go_left, [&](const handle_t& neighbor) {
                            
                            if (subgraph_nodes_to_add.count(xindex->get_id(neighbor))) {
                                // we've found a new node that we haven't added yet, add it to the graph
                                // and the queue, and erase from the nodes left to add
                                subgraph_nodes_to_add.erase(xindex->get_id(neighbor));
                                rescue_graph->create_handle(xindex->get_sequence(xindex->forward(neighbor)),
                                                            xindex->get_id(neighbor));
                                stack.push_back(neighbor);
                            }
                            
                            if (rescue_graph->has_node(xindex->get_id(neighbor))) {
                                // we always know that the handle we're coming from is in the subgraph, it
                                // seems that this neighbor is as well, so add the edge
                                
                                handle_t sub_handle = rescue_graph->get_handle(xindex->get_id(super_handle),
                                                                               xindex->get_is_reverse(super_handle));
                                handle_t sub_neighbor = rescue_graph->get_handle(xindex->get_id(neighbor),
                                                                                 xindex->get_is_reverse(neighbor));
                                
                                // edges automatically deduplicate, so don't worry about checking whether
                                // it exists
                                if (go_left) {
                                    rescue_graph->create_edge(sub_neighbor, sub_handle);
                                }
                                else {
                                    rescue_graph->create_edge(sub_handle, sub_neighbor);
                                }
                            }
                        });
                    }
                }
            }
        }
    }
    
    void MultipathMapper::init_band_padding_memo() {
        band_padding_memo.clear();
        band_padding_memo.resize(band_padding_memo_size);
        
        for (size_t i = 0; i < band_padding_memo.size(); i++) {
            band_padding_memo[i] = size_t(band_padding_multiplier * sqrt(i)) + 1;
        }
    }

    void MultipathMapper::set_alignment_scores(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend,
                                               int8_t full_length_bonus) {
        AlignerClient::set_alignment_scores(match, mismatch, gap_open, gap_extend, full_length_bonus);
        splice_stats.update_scoring(*get_regular_aligner());
        set_min_softclip_length_for_splice(min_softclip_length_for_splice);
        set_log_odds_against_splice(no_splice_natural_log_odds);
        set_max_merge_supression_length();
    }

    void MultipathMapper::set_alignment_scores(std::istream& matrix_stream, int8_t gap_open, int8_t gap_extend,
                                               int8_t full_length_bonus) {
        AlignerClient::set_alignment_scores(matrix_stream, gap_open, gap_extend, full_length_bonus);
        splice_stats.update_scoring(*get_regular_aligner());
        set_min_softclip_length_for_splice(min_softclip_length_for_splice);
        set_log_odds_against_splice(no_splice_natural_log_odds);
        set_max_merge_supression_length();
    }

    void MultipathMapper::set_alignment_scores(const int8_t* score_matrix, int8_t gap_open, int8_t gap_extend,
                                               int8_t full_length_bonus) {
        AlignerClient::set_alignment_scores(score_matrix, gap_open, gap_extend, full_length_bonus);
        splice_stats.update_scoring(*get_regular_aligner());
        set_min_softclip_length_for_splice(min_softclip_length_for_splice);
        set_log_odds_against_splice(no_splice_natural_log_odds);
        set_max_merge_supression_length();
    }

    
    bool MultipathMapper::likely_mismapping(const multipath_alignment_t& multipath_aln) {
        if (!suppress_mismapping_detection) {
            
            auto p_val = random_match_p_value(pseudo_length(multipath_aln), multipath_aln.sequence().size());
            
#ifdef debug_multipath_mapper
            cerr << "effective match length of read " << multipath_aln.sequence() << " is " << pseudo_length(multipath_aln) << " in read length " << multipath_aln.sequence().size() << ", yielding p-value " << p_val << endl;
#endif
            
            return p_val > max_mapping_p_value;
        }
        else {
            return false;
        }
    }

    bool MultipathMapper::likely_misrescue(const multipath_alignment_t& multipath_aln) {
        auto p_val = random_match_p_value(pseudo_length(multipath_aln), multipath_aln.sequence().size());
        
#ifdef debug_multipath_mapper
        cerr << "effective match length of rescued read " << multipath_aln.sequence() << " is " << pseudo_length(multipath_aln) << " in read length " << multipath_aln.sequence().size() << ", yielding p-value " << p_val << endl;
#endif
        
        return p_val > max_rescue_p_value;
    }

    
    int64_t MultipathMapper::pseudo_length(const multipath_alignment_t& multipath_aln) const {
        return optimal_alignment_score(multipath_aln);
    }
    
    // make the memo live in this .o file
    thread_local unordered_map<pair<int64_t, size_t>, double> MultipathMapper::p_value_memo;
    
    double MultipathMapper::random_match_p_value(int64_t match_length, size_t read_length) {
        // memoized to avoid transcendental functions (at least in cases where read lengths don't vary too much)
        auto iter = p_value_memo.find(make_pair(match_length, read_length));
        if (iter != p_value_memo.end()) {
            return iter->second;
        }
        else {
            double rate = max_exponential_rate_intercept + max_exponential_rate_slope * read_length;
            double shape = exp(max_exponential_shape_intercept + max_exponential_shape_slope * read_length);
            double p_value = 1.0 - max_exponential_cdf(match_length, rate, shape);
            if (p_value_memo.size() < max_p_value_memo_size && !suppress_p_value_memoization) {
                p_value_memo[make_pair(match_length, read_length)] = p_value;
            }
            return p_value;
        }
    }
    
    void MultipathMapper::calibrate_mismapping_detection(size_t num_simulations, const vector<size_t>& simulated_read_lengths) {
        
        // we are calibrating the parameters, so we don't want to memoize any p-values using the default values
        suppress_p_value_memoization = true;
        
        // we don't want to do base quality adjusted alignments for this stage since we are just simulating random sequences
        // with no base qualities
        bool reset_quality_adjustments = adjust_alignments_for_base_quality;
        adjust_alignments_for_base_quality = false;
        
        // these p-values will eventually be used internally where the scores still have the bonus applied, so we
        // want to make sure we don't strip them off here
        bool reset_strip_bonuses = strip_bonuses;
        strip_bonuses = false;
        
        // we want to avoid returning empty alignments that indicated unmapped reads
        bool reset_suppress_mismapping_detection = suppress_mismapping_detection;
        suppress_mismapping_detection = true;
        
        // and we expect small MEMs, so don't filter them out
        int reset_min_mem_length = min_mem_length;
        size_t reset_min_clustering_mem_length = min_clustering_mem_length;
        min_mem_length = 1;
        min_clustering_mem_length = 1;
        
        // and these reads are slow to map, but we only need primaries
        size_t reset_max_alt_mappings = max_alt_mappings;
        max_alt_mappings = 1;
        
        // reset the memo of p-values (which we are calibrating) for any updates using the default parameter during the null mappings
        p_value_memo.clear();
        
        // the logarithms of the MLE estimators at each read length
        vector<double> mle_max_exponential_rates;
        vector<double> log_mle_max_exponential_shapes;
        
        for (const size_t simulated_read_length : simulated_read_lengths) {
            // compute the pseudo length of a bunch of randomly generated sequences
            vector<double> pseudo_lengths(num_simulations, 0.0);
#pragma omp parallel for
            for (size_t i = 0; i < num_simulations; i++) {
                
                Alignment alignment;
                alignment.set_sequence(pseudo_random_sequence(simulated_read_length, i * 716293332 + simulated_read_length));
                vector<multipath_alignment_t> multipath_alns;
                multipath_map(alignment, multipath_alns);
                
                if (!multipath_alns.empty()) {
                    pseudo_lengths[i] = pseudo_length(multipath_alns.front());
                }
            }
            
            auto max_exp_params = fit_max_exponential(pseudo_lengths);
            mle_max_exponential_rates.push_back(max_exp_params.first);
            log_mle_max_exponential_shapes.push_back(log(max_exp_params.second));
                        
#ifdef debug_report_startup_training
            unordered_map<int64_t, size_t> length_counts;
            for (auto length : pseudo_lengths) {
                length_counts[round(length)]++;
            }
            vector<pair<size_t, size_t>> sorted_length_counts(length_counts.begin(), length_counts.end());
            sort(sorted_length_counts.begin(), sorted_length_counts.end());
            cerr << "data for length " << simulated_read_length << endl;
            for (auto length_count : sorted_length_counts) {
                cerr << "\t" << length_count.first << ": " << length_count.second << endl;
            }
            cerr << "trained parameters for length " << simulated_read_length << ": " << endl;
            cerr << "\tmax exp rate: " << max_exp_params.first << endl;
            cerr << "\tmax exp shape: " << max_exp_params.second << endl;
#endif
        }
        
        // make a design matrix for a log regression and a linear regression
        vector<vector<double>> X(simulated_read_lengths.size());
        for (size_t i = 0; i < X.size(); ++i) {
            X[i].resize(2, 1.0);
            X[i][1] = simulated_read_lengths[i];
        }
        
        auto max_exp_rate_coefs = regress(X, mle_max_exponential_rates);
        auto max_exp_shape_coefs = regress(X, log_mle_max_exponential_shapes);
        
        max_exponential_rate_intercept = max_exp_rate_coefs[0];
        max_exponential_rate_slope = max_exp_rate_coefs[1];
        max_exponential_shape_intercept = max_exp_shape_coefs[0];
        max_exponential_shape_slope = max_exp_shape_coefs[1];
        
#ifdef debug_report_startup_training
        cerr << "final regression parameters:" << endl;
        cerr << "\tmax exp rate = " << max_exponential_rate_intercept << " + " << max_exponential_rate_slope << " * L" << endl;
        cerr << "\tmax exp shape = exp(" << max_exponential_shape_intercept << " + " << max_exponential_shape_slope << " * L)" << endl;
#endif
        
        // reset mapping parameters to their original values
        adjust_alignments_for_base_quality = reset_quality_adjustments;
        strip_bonuses = reset_strip_bonuses;
        min_clustering_mem_length = reset_min_clustering_mem_length;
        min_mem_length = reset_min_mem_length;
        max_alt_mappings = reset_max_alt_mappings;
        suppress_p_value_memoization = false;
        suppress_mismapping_detection = reset_suppress_mismapping_detection;
    }

    void MultipathMapper::determine_distance_correlation() {
        
        // FIXME: very experimental, not sure if i actually want to do this
        
        if (!distance_index || ref_path_handles.empty()) {
            return;
        }
        
        vector<path_handle_t> refs(ref_path_handles.begin(), ref_path_handles.end());
        sort(refs.begin(), refs.end(), [&](path_handle_t a, path_handle_t b) {
            return xindex->get_path_name(a) < xindex->get_path_name(b);
        });
        vector<size_t> ref_weights(refs.size());
        for (size_t i = 0; i < refs.size(); ++i) {
            ref_weights[i] = xindex->get_path_length(refs[i]);
        }
        
        SnarlOrientedDistanceMeasurer measurer(distance_index);
        
        discrete_distribution<size_t> ref_distr(ref_weights.begin(), ref_weights.end());
        mt19937 gen;
        gen.seed(749753582ul);
        
        vector<tuple<path_handle_t, size_t, size_t, pos_t, pos_t>> positions;
        vector<vector<double>> distances;
        
        int num_measurements = 100 * 1000;
        int radius = 50 * 1000;
        for (int i = 0; i < num_measurements; ++i) {
            
            path_handle_t ref = refs[ref_distr(gen)];
            if (xindex->get_path_length(ref) < 2) {
                continue;
            }
            
            // don't include the past-the-last, which makes things tricky
            uniform_int_distribution<size_t> off_distr_1(1, xindex->get_path_length(ref) - 1);
            size_t path_off_1 = off_distr_1(gen);
            
            uniform_int_distribution<size_t> off_distr_2(max<int64_t>(path_off_1 - radius, 1),
                                                         min<int64_t>(path_off_1 + radius, xindex->get_path_length(ref) - 1));
            size_t path_off_2 = off_distr_2(gen);
            
            auto step_1 = xindex->get_step_at_position(ref, path_off_1);
            auto step_2 = xindex->get_step_at_position(ref, path_off_2);
            
            pos_t pos_1(xindex->get_id(xindex->get_handle_of_step(step_1)),
                       xindex->get_is_reverse(xindex->get_handle_of_step(step_1)),
                       path_off_1 - xindex->get_position_of_step(step_1));
            pos_t pos_2(xindex->get_id(xindex->get_handle_of_step(step_2)),
                       xindex->get_is_reverse(xindex->get_handle_of_step(step_2)),
                       path_off_2 - xindex->get_position_of_step(step_2));
            
            positions.emplace_back();
            get<0>(positions.back()) = ref;
            get<1>(positions.back()) = path_off_1;
            get<2>(positions.back()) = path_off_2;
            get<3>(positions.back()) = pos_1;
            get<4>(positions.back()) = pos_2;
            
            pos_t rev_pos_1(id(pos_1), !is_rev(pos_1),
                            xindex->get_length(xindex->get_handle_of_step(step_1)) - offset(pos_1));
            
            vector<int64_t> row;
            row.push_back(algorithms::ref_path_distance(xindex, pos_1, pos_2, ref_path_handles,
                                                        max_splice_ref_search_length));
            row.push_back(measurer.oriented_distance(pos_1, pos_2));
            row.push_back(algorithms::ref_path_distance(xindex, rev_pos_1, pos_2, ref_path_handles,
                                                        max_splice_ref_search_length));
            row.push_back(measurer.oriented_distance(rev_pos_1, pos_2));
            
            distances.emplace_back();
            for (auto d : row) {
                distances.back().push_back(d == numeric_limits<int64_t>::max() ? numeric_limits<double>::quiet_NaN() : d);
            }
            
            
        }
#ifdef debug_output_distance_correlation
        for (size_t j = 0; j < distances.size(); ++j) {
            auto& pos_row = positions[j];
            cerr << xindex->get_path_name(get<0>(pos_row));
            cerr << "\t" << get<1>(pos_row);
            cerr << "\t" << get<2>(pos_row);
            cerr << "\t" << get<3>(pos_row);
            cerr << "\t" << get<4>(pos_row);
            auto& row = distances[j];
            for (size_t i = 0; i < row.size(); ++i) {
                cerr << "\t";
                cerr << row[i];
            }
            cerr << endl;
        }
#endif
    }

    unique_ptr<OrientedDistanceMeasurer> MultipathMapper::get_distance_measurer(MemoizingGraph& memoizing_graph) const {
        
        unique_ptr<OrientedDistanceMeasurer> distance_measurer;
        if (distance_index) {
#ifdef debug_multipath_mapper
            cerr << "using a snarl-based distance measurer (if doing oriented distance clustering)" << endl;
#endif
            distance_measurer = unique_ptr<OrientedDistanceMeasurer>(new SnarlOrientedDistanceMeasurer(distance_index));
        }
        else {
#ifdef debug_multipath_mapper
            cerr << "using a path-based distance measurer (if doing oriented distance clustering)" << endl;
#endif
            distance_measurer = unique_ptr<OrientedDistanceMeasurer>(new PathOrientedDistanceMeasurer(&memoizing_graph,
                                                                                                      path_component_index.get()));
        }
        return distance_measurer;
    }
    
    int64_t MultipathMapper::distance_between(const multipath_alignment_t& multipath_aln_1,
                                              const multipath_alignment_t& multipath_aln_2,
                                              bool full_fragment, bool forward_strand) const {
                                              
        if (multipath_aln_1.subpath_size() == 0 || multipath_aln_2.subpath_size() == 0) {
            // Something is an unmapped alignment
            return numeric_limits<int64_t>::max();
        }
        
        Alignment aln_1;
        optimal_alignment(multipath_aln_1, aln_1);
        pos_t pos_1;
        if (aln_1.path().mapping_size() != 0) {
            pos_1 = initial_position(aln_1.path());
        }
        else {
            // the whole thing is negative scoring, take an arbitray postition just to not break
            pos_1 = make_pos_t(multipath_aln_1.subpath().front().path().mapping().front().position());
        }
        
        Alignment aln_2;
        optimal_alignment(multipath_aln_2, aln_2);
        pos_t pos_2;
        if (aln_2.path().mapping_size() != 0) {
            pos_2 = full_fragment ? final_position(aln_2.path()) : initial_position(aln_2.path());
        }
        else {
            // the whole thing is negative scoring, take an arbitray postition just to not break
            pos_2 = make_pos_t(multipath_aln_2.subpath().front().path().mapping().front().position());
        }
#ifdef debug_multipath_mapper
        cerr << "measuring left-to-" << (full_fragment ? "right" : "left") << " end distance between " << pos_1 << " and " << pos_2 << endl;
#endif
        
        int64_t dist;
        if (use_min_dist_clusterer || use_tvs_clusterer) {
            assert(!forward_strand);
            // measure the distance in both directions and choose the minimum (or the only) absolute distance
            size_t forward_dist = minimum_distance(*distance_index,pos_1, pos_2, false, xindex);
            size_t reverse_dist = minimum_distance(*distance_index,pos_2, pos_1, false, xindex);
            if (forward_dist == std::numeric_limits<size_t>::max() && reverse_dist == std::numeric_limits<size_t>::max()) {
                // unreachable both ways, convert to the sentinel that the client code expects
                dist = numeric_limits<int64_t>::max();
            }
            else if (forward_dist == std::numeric_limits<size_t>::max() || 
                     (reverse_dist < forward_dist && reverse_dist != std::numeric_limits<size_t>::max())) {
                dist = - (int64_t)reverse_dist;
            }
            else {
                dist = forward_dist;
            }
        }
        else {
            PathOrientedDistanceMeasurer measurer(xindex);
            dist = measurer.oriented_distance(pos_1, pos_2);
        }
        return dist;
    }

    int64_t MultipathMapper::distance(const pos_t& pos_1, const pos_t& pos_2) const {
        if (distance_index) {
            size_t distance = minimum_distance(*distance_index, pos_1, pos_2, false, xindex);
            return distance == std::numeric_limits<size_t>::max() ? -1 : (int64_t)distance;
        }
        else {
            return PathOrientedDistanceMeasurer(xindex).oriented_distance(pos_1, pos_2);
        }
    }
    
    bool MultipathMapper::is_consistent(int64_t distance) const {
        return (distance < fragment_length_distr.mean() + 10.0 * fragment_length_distr.std_dev()
                && distance > fragment_length_distr.mean() - 10.0 * fragment_length_distr.std_dev());
    }
    
    bool MultipathMapper::are_consistent(const multipath_alignment_t& multipath_aln_1,
                                         const multipath_alignment_t& multipath_aln_2) const {
        
        return is_consistent(distance_between(multipath_aln_1, multipath_aln_2, true));
    }
    
    bool MultipathMapper::share_terminal_positions(const multipath_alignment_t& multipath_aln_1,
                                                   const multipath_alignment_t& multipath_aln_2) const {
        
        unordered_set<pos_t> terminal_positions;
        
        // first look for matching starts
        for (size_t i = 0; i < multipath_aln_1.start_size(); i++) {
            terminal_positions.insert(make_pos_t(multipath_aln_1.subpath(multipath_aln_1.start(i)).path().mapping(0).position()));
        }
        
        for (size_t i = 0; i < multipath_aln_2.start_size(); i++) {
            if (terminal_positions.count(make_pos_t(multipath_aln_2.subpath(multipath_aln_2.start(i)).path().mapping(0).position()))) {
                return true;
            }
        }
        
        // remove the starts
        terminal_positions.clear();
        
        // now look for matching ends
        for (size_t i = 0; i < multipath_aln_1.subpath_size(); i++) {
            const subpath_t& subpath = multipath_aln_1.subpath(i);
            if (subpath.next_size() == 0) {
                terminal_positions.insert(final_position(subpath.path()));
            }
        }
        
        for (size_t i = 0; i < multipath_aln_2.subpath_size(); i++) {
            const subpath_t& subpath = multipath_aln_2.subpath(i);
            if (subpath.next_size() == 0) {
                if (terminal_positions.count(final_position(subpath.path()))) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    bool MultipathMapper::align_to_cluster_graphs_with_rescue(const Alignment& alignment1, const Alignment& alignment2,
                                                              vector<clustergraph_t>& cluster_graphs1,
                                                              vector<clustergraph_t>& cluster_graphs2,
                                                              vector<MaximalExactMatch>& mems1,
                                                              vector<MaximalExactMatch>& mems2,
                                                              vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                              vector<pair<pair<size_t, size_t>, int64_t>>& pair_distances,
                                                              vector<double>& pair_multiplicities,
                                                              const match_fanouts_t* fanouts1, const match_fanouts_t* fanouts2) {
        
        // align the two ends independently
        vector<multipath_alignment_t> multipath_alns_1, multipath_alns_2;
        vector<size_t> cluster_idxs_1, cluster_idxs_2;
        vector<double> multiplicities_1, multiplicities_2;
        align_to_cluster_graphs(alignment1,
                                cluster_graphs1, multipath_alns_1, multiplicities_1, max_single_end_mappings_for_rescue,
                                fanouts1, &cluster_idxs_1);
        align_to_cluster_graphs(alignment2,
                                cluster_graphs2, multipath_alns_2, multiplicities_2, max_single_end_mappings_for_rescue,
                                fanouts2, &cluster_idxs_2);
        
        if (!multipath_alns_1.empty() &&
            !multipath_alns_2.empty() &&
            multipath_alns_1.front().mapping_quality() >= min(60, max_mapping_quality) &&
            multipath_alns_2.front().mapping_quality() >= min(60, max_mapping_quality) &&
            are_consistent(multipath_alns_1.front(), multipath_alns_2.front())) {
            
            // we are able to obtain confident matches that satisfy the pairing constraints
#ifdef debug_multipath_mapper
            cerr << "found consistent, confident pair mapping from independent end mapping" << endl;
#endif
            multipath_aln_pairs_out.emplace_back(move(multipath_alns_1.front()), move(multipath_alns_2.front()));
            pair_distances.emplace_back(make_pair(cluster_idxs_1.front(), cluster_idxs_2.front()),
                                        distance_between(multipath_aln_pairs_out.back().first, multipath_aln_pairs_out.back().second, true));
            pair_multiplicities.emplace_back(min(multiplicities_1.front(), multiplicities_2.front()));
            return true;
        }
        
        // figure out how many rescues we will do and could do from each side
        
        int32_t max_score_diff = get_aligner(!alignment1.quality().empty() &&
                                             !alignment2.quality().empty())->mapping_quality_score_diff(max_mapping_quality);
        
        int32_t top_score_1 = multipath_alns_1.empty() ? 0 : optimal_alignment_score(multipath_alns_1.front());
        int32_t top_score_2 = multipath_alns_2.empty() ? 0 : optimal_alignment_score(multipath_alns_2.front());
        
        size_t num_rescuable_alns_1 = multipath_alns_1.size();
        size_t num_rescuable_alns_2 = multipath_alns_2.size();
        for (size_t i = 0; i < num_rescuable_alns_1; i++){
            if (likely_mismapping(multipath_alns_1[i]) ||
                (i > 0 ? optimal_alignment_score(multipath_alns_1[i]) < top_score_1 - max_score_diff : false)) {
                num_rescuable_alns_1 = i;
                break;
            }
        }
        for (size_t i = 0; i < num_rescuable_alns_2; i++){
            if (likely_mismapping(multipath_alns_2[i]) ||
                (i > 0 ? optimal_alignment_score(multipath_alns_2[i]) < top_score_2 - max_score_diff : false)) {
                num_rescuable_alns_2 = i;
                break;
            }
        }
        size_t num_to_rescue_1 = min(num_rescuable_alns_1, max_rescue_attempts);
        size_t num_to_rescue_2 = min(num_rescuable_alns_2, max_rescue_attempts);
        
#ifdef debug_multipath_mapper
        cerr << "rescuing from " << num_to_rescue_1 << " read1's and " << num_to_rescue_2 << " read2's" << endl;
#endif
        
        // calculate the estimated multiplicity of a pair found from each of the two ends
        double estimated_multiplicity_from_1 = num_to_rescue_1 > 0 ? double(num_rescuable_alns_1) / num_to_rescue_1 : 1.0;
        double estimated_multiplicity_from_2 = num_to_rescue_2 > 0 ? double(num_rescuable_alns_2) / num_to_rescue_2 : 1.0;
        
        // actually doe the rescues and record which ones succeeded
        vector<multipath_alignment_t> rescue_multipath_alns_1(num_to_rescue_2), rescue_multipath_alns_2(num_to_rescue_1);
        unordered_set<size_t> rescued_from_1, rescued_from_2;
        
        for (size_t i = 0; i < num_to_rescue_1; i++) {
            multipath_alignment_t rescue_multipath_aln;
            if (attempt_rescue(multipath_alns_1[i], alignment2, true, rescue_multipath_aln)) {
                rescued_from_1.insert(i);
                rescue_multipath_alns_2[i] = move(rescue_multipath_aln);
            }
        }
        
        for (size_t i = 0; i < num_to_rescue_2; i++) {
            multipath_alignment_t rescue_multipath_aln;
            if (attempt_rescue(multipath_alns_2[i], alignment1, false, rescue_multipath_aln)) {
                rescued_from_2.insert(i);
                rescue_multipath_alns_1[i] = move(rescue_multipath_aln);
            }
        }
        
        // follow some complicated logic to check if any of the rescued alignments are duplicates
        // of the original alignments
        
        bool found_consistent = false;
        if (!rescued_from_1.empty() && !rescued_from_2.empty()) {
#ifdef debug_multipath_mapper
            cerr << "successfully rescued from both read ends" << endl;
#endif
            
            unordered_set<size_t> found_duplicate;
            
            // for each rescue attempt from a read 1
            for (size_t i  : rescued_from_1) {
                bool duplicate = false;
                for (size_t j : rescued_from_2) {
                    if (found_duplicate.count(j)) {
                        continue;
                    }
                    
#ifdef debug_multipath_mapper
                    cerr << "checking duplication between mapped read1 " << i << " and rescued read1 " << j << endl;
#endif
                    if (share_terminal_positions(multipath_alns_1[i], rescue_multipath_alns_1[j])) {
#ifdef debug_multipath_mapper
                        cerr << "found duplicate, now checking rescued read2 " << i << " and mapped read2 " << j << endl;
#endif
                        if (share_terminal_positions(rescue_multipath_alns_2[i], multipath_alns_2[j])) {
#ifdef debug_multipath_mapper
                            cerr << "found duplicate, marking entire pair as duplicate" << endl;
#endif
                            // these two alignments found each other with their rescue, we don't want to add duplicate mappings
                            duplicate = true;
                            found_duplicate.insert(j);
                            
                            // move the original mappings
                            int64_t dist = distance_between(multipath_alns_1[i], multipath_alns_2[j], true);
                            if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
                                multipath_aln_pairs_out.emplace_back(move(multipath_alns_1[i]), move(multipath_alns_2[j]));
                                pair_distances.emplace_back(make_pair(cluster_idxs_1[i], cluster_idxs_2[j]), dist);
                                pair_multiplicities.emplace_back(min(estimated_multiplicity_from_1 * multiplicities_1[i],
                                                                     estimated_multiplicity_from_2 * multiplicities_2[j]));
                                found_consistent = true;
                            }
                            
                            break;
                        }
                    }
                }
                
                // if we haven't already moved the pair and marked it as a duplicate, move the rescued pair into the output vector
                if (!duplicate) {
                    int64_t dist = distance_between(multipath_alns_1[i], rescue_multipath_alns_2[i], true);
                    if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
#ifdef debug_multipath_mapper
                        cerr << "adding read1 and rescued read2 " << i << " to output vector" << endl;
#endif
                        multipath_aln_pairs_out.emplace_back(move(multipath_alns_1[i]), move(rescue_multipath_alns_2[i]));
                        pair_distances.emplace_back(make_pair(cluster_idxs_1[i], RESCUED), dist);
                        pair_multiplicities.emplace_back(estimated_multiplicity_from_1 * multiplicities_1[i]);
                        found_consistent = true;
                    }
                }
            }
            
            // for each rescue attempt from a read 2
            for (size_t j : rescued_from_2) {
                if (found_duplicate.count(j)) {
                    // we already moved it as part of a duplicate pair
                    continue;
                }
                int64_t dist = distance_between(rescue_multipath_alns_1[j], multipath_alns_2[j], true);
                if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
#ifdef debug_multipath_mapper
                    cerr << "adding rescued read1 and read2 " << j << " to output vector" << endl;
#endif
                    multipath_aln_pairs_out.emplace_back(move(rescue_multipath_alns_1[j]), move(multipath_alns_2[j]));
                    pair_distances.emplace_back(make_pair(RESCUED, cluster_idxs_2[j]), dist);
                    pair_multiplicities.emplace_back(estimated_multiplicity_from_2 * multiplicities_2[j]);
                    found_consistent = true;
                }
            }
        }
        else if (!rescued_from_1.empty()) {
#ifdef debug_multipath_mapper
            cerr << "successfully rescued from only read 1" << endl;
#endif
            for (size_t i: rescued_from_1) {
                int64_t dist = distance_between(multipath_alns_1[i], rescue_multipath_alns_2[i], true);
                if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
                    multipath_aln_pairs_out.emplace_back(move(multipath_alns_1[i]), move(rescue_multipath_alns_2[i]));
                    pair_distances.emplace_back(make_pair(cluster_idxs_1[i], RESCUED), dist);
                    pair_multiplicities.emplace_back(estimated_multiplicity_from_1 * multiplicities_1[i]);
                    found_consistent = true;
                }
            }
        }
        else if (!rescued_from_2.empty()) {
#ifdef debug_multipath_mapper
            cerr << "successfully rescued from only read 2" << endl;
#endif
            for (size_t i : rescued_from_2) {
                int64_t dist = distance_between(rescue_multipath_alns_1[i], multipath_alns_2[i], true);
                if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
                    multipath_aln_pairs_out.emplace_back(move(rescue_multipath_alns_1[i]), move(multipath_alns_2[i]));
                    pair_distances.emplace_back(make_pair(RESCUED, cluster_idxs_2[i]), dist);
                    pair_multiplicities.emplace_back(estimated_multiplicity_from_2 * multiplicities_2[i]);
                    found_consistent = true;
                }
            }
        }
        
        if (!found_consistent && do_spliced_alignment) {
#ifdef debug_multipath_mapper
            cerr << "rescue failed, doing independent spliced alignment and then re-attempting pairing" << endl;
#endif
            
#ifdef debug_time_phases
            clock_t start = clock();
#endif
            
            multipath_alignment_t* rescue_anchor_1 = nullptr;
            multipath_alignment_t* rescue_anchor_2 = nullptr;
            double rescue_multiplicity_1 = 1.0, rescue_multiplicity_2 = 1.0;
            if (!multipath_alns_1.empty()) {
                rescue_anchor_1 = &multipath_alns_1.front();
                rescue_multiplicity_1 = multiplicities_1.front();
            }
            if (!multipath_alns_2.empty()) {
                rescue_anchor_2 = &multipath_alns_2.front();
                rescue_multiplicity_2 = multiplicities_2.front();
            }
            
            // find splices independently, also use the mate to rescue missing splice segments
            bool did_splice_1 = false, did_splice_2 = false;
            if (!multipath_alns_1.empty() && !likely_mismapping(multipath_alns_1.front())) {
                did_splice_1 = find_spliced_alignments(alignment1, multipath_alns_1, multiplicities_1, cluster_idxs_1,
                                                       mems1, cluster_graphs1, fanouts1,
                                                       rescue_anchor_2, true, rescue_multiplicity_2);
            }
            if (!multipath_alns_2.empty() && !likely_mismapping(multipath_alns_2.front())) {
                did_splice_2 = find_spliced_alignments(alignment2, multipath_alns_2, multiplicities_2, cluster_idxs_2,
                                                       mems2, cluster_graphs2, fanouts2,
                                                       rescue_anchor_1, false, rescue_multiplicity_1);
            }
            
            if (did_splice_1 || did_splice_2) {
                // it may now be possible to identify some pairs as properly paired using the spliced alignment
                found_consistent = retry_pairing_spliced_alignments(alignment1, alignment2, multipath_alns_1, multipath_alns_2,
                                                                    cluster_idxs_1, cluster_idxs_2, multiplicities_1,
                                                                    multiplicities_2, multipath_aln_pairs_out,
                                                                    pair_distances, pair_multiplicities);
            }
            
#ifdef debug_time_phases
            cerr << "\trescue failed, spliced alignment used: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
        }
        
        if (!found_consistent) {
            
#ifdef debug_multipath_mapper
            cerr << "failed to successfully rescue from either read end, reporting independent mappings" << endl;
#endif
            
            // agglomerate them them independently if necessary
            if (agglomerate_multipath_alns) {
                agglomerate_alignments(multipath_alns_1, &multiplicities_1);
                agglomerate_alignments(multipath_alns_2, &multiplicities_2);
            }
            
            // rescue failed, so we just report these as independent mappings
            size_t num_pairs_to_report = min(max_alt_mappings, max(multipath_alns_1.size(), multipath_alns_2.size()));
            
            // move the multipath alignments to the return vector
            multipath_aln_pairs_out.reserve(num_pairs_to_report);
            for (size_t i = 0; i < num_pairs_to_report; i++) {
                if (i < multipath_alns_1.size() && i < multipath_alns_2.size()) {
                    multipath_aln_pairs_out.emplace_back(move(multipath_alns_1[i]), move(multipath_alns_2[i]));
                    
                }
                else if (i < multipath_alns_1.size()) {
                    multipath_aln_pairs_out.emplace_back(move(multipath_alns_1[i]), multipath_alignment_t());
                    to_multipath_alignment(alignment2, multipath_aln_pairs_out.back().second);
                    multipath_aln_pairs_out.back().second.clear_subpath();
                    multipath_aln_pairs_out.back().second.clear_start();
                }
                else {
                    multipath_aln_pairs_out.emplace_back(multipath_alignment_t(), move(multipath_alns_2[i]));
                    to_multipath_alignment(alignment1, multipath_aln_pairs_out.back().first);
                    multipath_aln_pairs_out.back().first.clear_subpath();
                    multipath_aln_pairs_out.back().first.clear_start();
                }
            }
        }
        
        if (found_consistent) {
            // compute the paired mapping quality
            sort_and_compute_mapping_quality(multipath_aln_pairs_out, pair_distances, nullptr, &pair_multiplicities);
        }
        
#ifdef debug_validate_multipath_alignments
        for (pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : multipath_aln_pairs_out) {
#ifdef debug_multipath_mapper
            cerr << "validating multipath alignments:" << endl;
            cerr << debug_string(multipath_aln_pair.first) << endl;
            cerr << debug_string(multipath_aln_pair.second) << endl;
#endif
            if (!validate_multipath_alignment(multipath_aln_pair.first, *xindex)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.first.sequence() << " failed to validate" << endl;
            }
            if (!validate_multipath_alignment(multipath_aln_pair.second, *xindex)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.second.sequence() << " failed to validate" << endl;
            }
        }
#endif
        
        return found_consistent;
    }
    
    void MultipathMapper::attempt_rescue_for_secondaries(const Alignment& alignment1, const Alignment& alignment2,
                                                         vector<clustergraph_t>& cluster_graphs1,
                                                         vector<clustergraph_t>& cluster_graphs2,
                                                         vector<pair<size_t, size_t>>& duplicate_pairs,
                                                         vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                         vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                         vector<double>& pair_multiplicities,
                                                         const match_fanouts_t* fanouts1, const match_fanouts_t* fanouts2) {
        
#ifdef debug_multipath_mapper
        cerr << "using rescue to find secondary mappings" << endl;
#endif

        unordered_set<size_t> paired_clusters_1, paired_clusters_2;
        
        for (size_t i = 0; i < multipath_aln_pairs_out.size(); i++) {
            // keep track of which clusters already have consistent pairs
            paired_clusters_1.insert(cluster_pairs[i].first.first);
            paired_clusters_2.insert(cluster_pairs[i].first.second);
        }
        
        for (size_t i = 0; i < duplicate_pairs.size(); i++) {
            // also mark which pairs have already been identified as duplicates
            paired_clusters_1.insert(duplicate_pairs[i].first);
            paired_clusters_2.insert(duplicate_pairs[i].second);
        }
        
        auto aligner = get_aligner(!alignment1.quality().empty() && !alignment2.quality().empty());
        int32_t cluster_score_1 = aligner->match * get<2>(cluster_graphs1[cluster_pairs.front().first.first]);
        int32_t cluster_score_2 = aligner->match * get<2>(cluster_graphs2[cluster_pairs.front().first.second]);
        int32_t max_score_diff = secondary_rescue_score_diff * aligner->mapping_quality_score_diff(max_mapping_quality);
        
        vector<pair<multipath_alignment_t, multipath_alignment_t>> rescued_secondaries;
        vector<pair<pair<size_t, size_t>, int64_t>> rescued_distances;
        vector<double> rescued_multiplicities;
        
        auto align_and_rescue = [&](const Alignment& anchor_aln, const Alignment& rescue_aln,
                                    vector<clustergraph_t>& cluster_graphs, unordered_set<size_t>& paired_clusters,
                                    int32_t max_score, bool anchor_is_read_1, const match_fanouts_t* anchor_fanouts) {
            
#ifdef debug_multipath_mapper
            cerr << "checking for rescues from read " << (anchor_is_read_1 ? 1 : 2) << endl;
#endif
            // remember how many pairs are already in the vector
            size_t num_preexisting_pairs = rescued_secondaries.size();
            
            size_t num_rescuable = 0;
            size_t num_rescues = 0;
            for (size_t i = 0; i < cluster_graphs.size(); ++i) {
                if (paired_clusters.count(i)) {
                    // we already have a consistent pair from this cluster
#ifdef debug_multipath_mapper
                    cerr << "cluster " << i << " is already in a pair" << endl;
#endif

                    continue;
                }
                
#ifdef debug_multipath_mapper
                cerr << "cluster " << i << "'s approximate score is " << get<2>(cluster_graphs[i]) * aligner->match << ", looking for " << max_score - max_score_diff << endl;
#endif
                
                if (get<2>(cluster_graphs[i]) * aligner->match < max_score - max_score_diff) {
#ifdef debug_multipath_mapper
                    cerr << "the approximate score of the remaining is too low to consider" << endl;
#endif
                    // the approximate score of the remaining is too low to consider
                    break;
                }
                
                // count this one as potentially rescuable from
                ++num_rescuable;
                
                if (num_rescues >= secondary_rescue_attempts) {
                    // we have used up our budget of rescue attempts
                    continue;
                }
                
                ++num_rescues;
                
                // TODO: repetitive with align_to_cluster_graphs
                
                // make the alignment
                vector<multipath_alignment_t> cluster_multipath_alns;
                cluster_multipath_alns.emplace_back();
                multipath_align(anchor_aln, cluster_graphs[i], cluster_multipath_alns.back(), anchor_fanouts);
                
                if (!suppress_multicomponent_splitting) {
                    // split it up if it turns out to be multiple components
                    split_multicomponent_alignments(cluster_multipath_alns);
                }
                
                // order the subpaths
                for (multipath_alignment_t& multipath_aln : cluster_multipath_alns) {
                    topologically_order_subpaths(multipath_aln);
                }
                
                // if we split it up, move the best one to the front
                if (cluster_multipath_alns.size() > 1) {
                    sort_and_compute_mapping_quality(cluster_multipath_alns);
                }
                
                // rescue from the alignment
                multipath_alignment_t rescue_multipath_aln;
                if (!likely_mismapping(cluster_multipath_alns.front())) {
                    bool rescued = attempt_rescue(cluster_multipath_alns.front(), rescue_aln, anchor_is_read_1, rescue_multipath_aln);
#ifdef debug_multipath_mapper
                    cerr << "rescued alignment is " << debug_string(rescue_multipath_aln) << endl;
#endif
                    if (rescued) {
#ifdef debug_multipath_mapper
                        cerr << "rescue succeeded, adding to rescue pair vector" << endl;
#endif
                        if (anchor_is_read_1) {
                            int64_t dist = distance_between(cluster_multipath_alns.front(), rescue_multipath_aln, true);
                            if (dist >= 0 && dist != numeric_limits<int64_t>::max()) {
                                simplify_complicated_multipath_alignment(cluster_multipath_alns.front());
                                rescued_secondaries.emplace_back(move(cluster_multipath_alns.front()), move(rescue_multipath_aln));
                                rescued_distances.emplace_back(make_pair(i, RESCUED), dist);
                                
                            }
                        }
                        else {
                            int64_t dist = distance_between(rescue_multipath_aln, cluster_multipath_alns.front(), true);
                            if (dist >= 0 && dist != numeric_limits<int64_t>::max()) {
                                simplify_complicated_multipath_alignment(cluster_multipath_alns.front());
                                rescued_secondaries.emplace_back(move(rescue_multipath_aln), move(cluster_multipath_alns.front()));
                                rescued_distances.emplace_back(make_pair(RESCUED, i), dist);
                                
                            }
                        }
                    }
#ifdef debug_multipath_mapper
                    else {
                        cerr << "rescue failed" << endl;
                    }
#endif
                }
#ifdef debug_multipath_mapper
                else {
                    cerr << "alignment we're rescuing from is likely a mismapping" << endl;
                }
#endif
            }
            
            // estimate how many of these alignments there probably are in total
            double rescue_multiplicity = double(num_rescuable) / double(num_rescues);
            
            // fill out the multiplicity with estimated multiplicity based on rescue and cluster
            for (size_t i = num_preexisting_pairs; i < rescued_secondaries.size(); ++i) {
                const auto& rescued_cluster_pair = rescued_distances[i];
                double clust_multiplicity;
                if (rescued_cluster_pair.first.first == RESCUED) {
                    // the read 1 mapping is from a rescue, get the cluster multiplicity for read 2
                    clust_multiplicity = cluster_multiplicity(get<1>(cluster_graphs2[rescued_cluster_pair.first.second]));
                }
                else {
                    // the read 2 mapping is from a rescue, get the cluster multiplicity for read 1
                    clust_multiplicity = cluster_multiplicity(get<1>(cluster_graphs1[rescued_cluster_pair.first.first]));
                }
                rescued_multiplicities.push_back(rescue_multiplicity * clust_multiplicity);
            }
            
//#pragma omp atomic
//            SECONDARY_RESCUE_TOTAL++;
//#pragma omp atomic
//            SECONDARY_RESCUE_ATTEMPT += num_rescues;
//            if (num_rescues > 0) {
//#pragma omp atomic
//                SECONDARY_RESCUE_COUNT++;
//            }
        };
        
        // perform routine for both read ends
        align_and_rescue(alignment1, alignment2, cluster_graphs1, paired_clusters_1,
                         cluster_score_1, true, fanouts1);
        align_and_rescue(alignment2, alignment1, cluster_graphs2, paired_clusters_2,
                         cluster_score_2, false, fanouts2);
        
#ifdef debug_validate_multipath_alignments
        for (pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : rescued_secondaries) {
#ifdef debug_multipath_mapper
            cerr << "validating rescued secondary multipath alignments:" << endl;
            cerr << debug_string(multipath_aln_pair.first) << endl;
            cerr << debug_string(multipath_aln_pair.second) << endl;
#endif
            if (!validate_multipath_alignment(multipath_aln_pair.first, *xindex)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.first.sequence() << " failed to validate" << endl;
            }
            if (!validate_multipath_alignment(multipath_aln_pair.second, *xindex)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.second.sequence() << " failed to validate" << endl;
            }
        }
#endif
        
        if (!rescued_secondaries.empty()) {
            // we found mappings that could be rescues of each other
            
#ifdef debug_multipath_mapper
            cerr << "some rescues succeeded, deduplicating rescued pairs" << endl;
#endif
            
            // find any rescued pairs that are duplicates of each other
            vector<bool> duplicate(rescued_secondaries.size(), false);
            for (size_t i = 1; i < rescued_secondaries.size(); i++) {
                for (size_t j = 0; j < i; j++) {
                    if (share_terminal_positions(rescued_secondaries[i].first, rescued_secondaries[j].first)) {
                        if (share_terminal_positions(rescued_secondaries[i].second, rescued_secondaries[j].second)) {
                            duplicate[i] = true;
                            duplicate[j] = true;
                        }
                    }
                }
            }
            
            // move the duplicates to the end of the vector
            size_t end = rescued_secondaries.size();
            for (size_t i = 0; i < end; ) {
                if (duplicate[i]) {
                    
                    std::swap(rescued_secondaries[i], rescued_secondaries[end - 1]);
                    std::swap(rescued_distances[i], rescued_distances[end - 1]);
                    std::swap(rescued_multiplicities[i], rescued_multiplicities[end - 1]);
                    std::swap(duplicate[i], duplicate[end - 1]);
                    
                    end--;
                }
                else {
                    i++;
                }
            }
            
            // remove duplicates
            if (end < rescued_secondaries.size()) {
                rescued_secondaries.resize(end);
                rescued_distances.resize(end);
                rescued_multiplicities.resize(end);
            }
            
            // merge the rescued secondaries into the return vector
            merge_rescued_mappings(multipath_aln_pairs_out, cluster_pairs, pair_multiplicities,
                                   rescued_secondaries, rescued_distances, rescued_multiplicities);
        } else {
#ifdef debug_multipath_mapper
            cerr << "no rescues succeeded" << endl;
#endif
        }
    }
    
    bool MultipathMapper::multipath_map_paired(const Alignment& alignment1, const Alignment& alignment2,
                                               vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                               vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer) {

        //cerr << (to_string(omp_get_thread_num()) + " " + alignment1.name() + "\n");
#ifdef debug_multipath_mapper
        cerr << "multipath mapping paired reads " << pb2json(alignment1) << " and " << pb2json(alignment2) << endl;
#endif
        
#ifdef debug_time_phases
        cerr << "starting time measurement for pair " << alignment1.name() << endl;
        clock_t start = clock();
#endif
        
        // empty the output vector (just for safety)
        multipath_aln_pairs_out.clear();
        
        if (!fragment_length_distr.is_finalized()) {
            // we have not estimated a fragment length distribution yet, so we revert to single ended mode and look
            // for unambiguous pairings
            
#ifdef debug_multipath_mapper
            cerr << "no fragment length distribution yet, looking for unambiguous single ended pairs" << endl;
#endif
            
            return attempt_unpaired_multipath_map_of_pair(alignment1, alignment2, multipath_aln_pairs_out, ambiguous_pair_buffer);
        }
        
        // the fragment length distribution has been estimated, so we can do full-fledged paired mode
        vector<deque<pair<string::const_iterator, char>>> mem_fanouts1, mem_fanouts2;
        auto mems1 = find_mems(alignment1, &mem_fanouts1);
        auto mems2 = find_mems(alignment2, &mem_fanouts2);
        unique_ptr<match_fanouts_t> fanouts1(mem_fanouts1.empty() ? nullptr
                                             : new match_fanouts_t(record_fanouts(mems1, mem_fanouts1)));
        unique_ptr<match_fanouts_t> fanouts2(mem_fanouts2.empty() ? nullptr
                                             : new match_fanouts_t(record_fanouts(mems2, mem_fanouts2)));
#ifdef debug_time_phases
        cerr << "collected MEMs, time elapsed: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
        
#ifdef debug_multipath_mapper
        cerr << "obtained read1 MEMs:" << endl;
        for (MaximalExactMatch mem : mems1) {
            cerr << "\t" << mem << " (" << mem.nodes.size() << " hits filled out of " << mem.match_count << ")" << endl;
        }
        cerr << "obtained read2 MEMs:" << endl;
        for (MaximalExactMatch mem : mems2) {
            cerr << "\t" << mem << " (" << mem.nodes.size() << " hits filled out of " << mem.match_count << ")" << endl;
        }
#endif
        
        MemoizingGraph memoizing_graph(xindex);
        unique_ptr<OrientedDistanceMeasurer> distance_measurer = get_distance_measurer(memoizing_graph);
        
#ifdef debug_multipath_mapper
        cerr << "clustering MEMs on both read ends..." << endl;
#endif
        
        // try to rescue high count runs of order-length MEMs for both reads before clustering
        rescue_high_count_order_length_mems(mems1, order_length_repeat_hit_max);
        rescue_high_count_order_length_mems(mems2, order_length_repeat_hit_max);
        
        // do the clustering
        vector<memcluster_t> clusters1 = get_clusters(alignment1, mems1, &(*distance_measurer), fanouts1.get());
        vector<memcluster_t> clusters2 = get_clusters(alignment2, mems2, &(*distance_measurer), fanouts2.get());
        
#ifdef debug_time_phases
        cerr << "clustered MEMs, time elapsed: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
        
        // extract graphs around the clusters and get the assignments of MEMs to these graphs
        vector<clustergraph_t> cluster_graphs1 = query_cluster_graphs(alignment1, mems1, clusters1);
        vector<clustergraph_t> cluster_graphs2 = query_cluster_graphs(alignment2, mems2, clusters2);
        
#ifdef debug_time_phases
        cerr << "extracted subgraphs, time elapsed: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
        
#ifdef debug_multipath_mapper
        cerr << "obtained independent clusters:" << endl;
        cerr << "read 1" << endl;
        for (int i = 0; i < cluster_graphs1.size(); i++) {
            cerr << "\tcluster " << i << " (multiplicity " << get<1>(cluster_graphs1[i]).second << ")" << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : get<1>(cluster_graphs1[i]).first) {
                cerr << "\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
            }
        }
        cerr << "read 2" << endl;
        for (int i = 0; i < cluster_graphs2.size(); i++) {
            cerr << "\tcluster " << i << " (multiplicity " << get<1>(cluster_graphs2[i]).second << ")" << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : get<1>(cluster_graphs2[i]).first) {
                cerr << "\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
            }
        }
#endif
        
        // we haven't already obtained a paired mapping by rescuing into a repeat, so we should try to get one
        // by cluster pairing
        
        vector<pair<pair<size_t, size_t>, int64_t>> cluster_pairs = get_cluster_pairs(alignment1, alignment2,
                                                                                      cluster_graphs1, cluster_graphs2,
                                                                                      &(*distance_measurer));
#ifdef debug_time_phases
        cerr << "paired clusters, time elapsed: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
        
#ifdef debug_multipath_mapper
        cerr << "obtained cluster pairs:" << endl;
        for (int i = 0; i < cluster_pairs.size(); i++) {
            cerr << "\tpair "  << i << " at distance " << cluster_pairs[i].second << endl;
            cerr << "\t\t read 1 (cluster " << cluster_pairs[i].first.first <<  ")" << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : get<1>(cluster_graphs1[cluster_pairs[i].first.first]).first) {
                cerr << "\t\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
            }
            cerr << "\t\t read 2 (cluster " << cluster_pairs[i].first.second << ")" << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : get<1>(cluster_graphs2[cluster_pairs[i].first.second]).first) {
                cerr << "\t\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
            }
        }
#endif
        
        // initialize some pair variables
        vector<double> pair_multiplicities;
        vector<pair<size_t, size_t>> duplicate_pairs;

        bool proper_paired = true;
        // do we find any pairs that satisfy the distance requirements?
        if (!cluster_pairs.empty()) {
            // We got some pairs that satisfy the distance requirements.
            
            // only perform the mappings that satisfy the expectations on distance
            
            align_to_cluster_graph_pairs(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                         multipath_aln_pairs_out, cluster_pairs, pair_multiplicities,
                                         duplicate_pairs, fanouts1.get(), fanouts2.get());
            
#ifdef debug_time_phases
            cerr << "made alignments, time elapsed: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
            
            // do we produce at least one good looking pair alignments from the clustered clusters?
            if (multipath_aln_pairs_out.empty()
                || likely_mismapping(multipath_aln_pairs_out.front().first)
                || likely_mismapping(multipath_aln_pairs_out.front().second)) {
                
#ifdef debug_multipath_mapper
                cerr << "pair may be mismapped, attempting individual end mappings" << endl;
#endif
                // we're not happy with the pairs we got, try to get a good pair by rescuing from single ended alignments
                
                vector<pair<multipath_alignment_t, multipath_alignment_t>> rescue_aln_pairs;
                vector<pair<pair<size_t, size_t>, int64_t>> rescue_distances;
                vector<double> rescue_multiplicities;
                bool rescued = align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2, mems1,
                                                                   mems2, rescue_aln_pairs, rescue_distances, rescue_multiplicities,
                                                                   fanouts1.get(), fanouts2.get());
#ifdef debug_time_phases
                cerr << "attempted rescue, time elapsed: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
                
                if (rescued) {
                    // we found consistent pairs by rescue, merge the two lists
                    
#ifdef debug_multipath_mapper
                    cerr << "found some rescue pairs, merging into current list of consistent mappings" << endl;
#endif
                    
                    merge_rescued_mappings(multipath_aln_pairs_out, cluster_pairs, pair_multiplicities,
                                           rescue_aln_pairs, rescue_distances, rescue_multiplicities);

                }
                else if (multipath_aln_pairs_out.empty() ||
                         (!(!likely_mismapping(multipath_aln_pairs_out.front().first) &&
                            !likely_misrescue(multipath_aln_pairs_out.front().second)) ||
                          !(!likely_misrescue(multipath_aln_pairs_out.front().first) &&
                            !likely_mismapping(multipath_aln_pairs_out.front().second)))) {
                    
                    // rescue didn't find any consistent mappings and we didn't have any pairings
                    // that we would have accepted from rescue beforehand. just take the single ended
                    // mappings that were computed for the sake of rescue
                    
                    proper_paired = false;
                    std::swap(multipath_aln_pairs_out, rescue_aln_pairs);
                    
                    // Don't sort and compute mapping quality; preserve the single-ended MAPQs
                }
            }
            else {
                
                // We don't think any of our hits are likely to be mismapped
                
                if (multipath_aln_pairs_out.front().first.mapping_quality() >= max_mapping_quality - secondary_rescue_subopt_diff &&
                    multipath_aln_pairs_out.front().second.mapping_quality() >= max_mapping_quality - secondary_rescue_subopt_diff) {
                    
                    // we're very confident about this pair, but it might be because we over-pruned at the clustering stage
                    // or because of problems with the seeds. we use this routine to use rescue on other very good looking
                    // independent end clusters
                    
                    attempt_rescue_for_secondaries(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                   duplicate_pairs, multipath_aln_pairs_out, cluster_pairs,
                                                   pair_multiplicities, fanouts1.get(), fanouts2.get());
                    
#ifdef debug_time_phases
                    cerr << "attempted secondary rescue, time elapsed: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
                }
            }
        }
        else {
            // We got no pairs that satisfy the distance requirements
            
            // revert to independent single ended mappings, but skip any rescues that we already tried
            
#ifdef debug_multipath_mapper
            cerr << "could not find a consistent pair, reverting to single ended mapping" << endl;
#endif
            
            vector<double> rescue_multiplicities;
            proper_paired = align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2, mems1,
                                                                mems2, multipath_aln_pairs_out, cluster_pairs, rescue_multiplicities,
                                                                fanouts1.get(), fanouts2.get());
            
#ifdef debug_time_phases
            cerr << "no pairs, did alignments with rescue, time elapsed: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
            
            if (proper_paired) {
                // we'll want to remember the multiplicities
                pair_multiplicities = move(rescue_multiplicities);
            }
        }
        
        if (multipath_aln_pairs_out.empty()) {
            // we tried all of our tricks and still didn't find a mapping
            
            // add a null alignment so we know it wasn't mapped
            multipath_aln_pairs_out.emplace_back();
            to_multipath_alignment(alignment1, multipath_aln_pairs_out.back().first);
            to_multipath_alignment(alignment2, multipath_aln_pairs_out.back().second);
            pair_multiplicities.emplace_back();
            cluster_pairs.emplace_back();
            
            // in case we're realigning GAMs that have paths already
            multipath_aln_pairs_out.back().first.clear_subpath();
            multipath_aln_pairs_out.back().first.clear_start();
            multipath_aln_pairs_out.back().second.clear_subpath();
            multipath_aln_pairs_out.back().second.clear_start();
        }
        
        // do paired spliced alignment only if we have real pairs
        if (proper_paired && do_spliced_alignment) {
            find_spliced_alignments(alignment1, alignment2, multipath_aln_pairs_out, cluster_pairs, pair_multiplicities,
                                    mems1, mems2, cluster_graphs1, cluster_graphs2);
#ifdef debug_time_phases
            cerr << "formed spliced alignments, time elapsed: " << double(clock() - start) / CLOCKS_PER_SEC << " secs" << endl;
#endif
        }
        
        // only agglomerate if the pairs are true pairs, otherwise it gets too complicated
        // to estimate mapping qualities
        if (proper_paired && agglomerate_multipath_alns) {
            agglomerate_alignment_pairs(multipath_aln_pairs_out, cluster_pairs, pair_multiplicities);
        }
        
        // if we computed extra alignments to get a mapping quality or investigate ambiguous clusters, remove them
        if (multipath_aln_pairs_out.size() > max_alt_mappings) {
            multipath_aln_pairs_out.resize(max_alt_mappings);
        }
        
        // mark if any of the alignments are just noise
        purge_unmapped_alignments(multipath_aln_pairs_out, proper_paired);
        
        for (size_t i = 0; i < multipath_aln_pairs_out.size(); ++i) {
            multipath_aln_pairs_out[i].first.set_annotation("proper_pair", proper_paired);
            multipath_aln_pairs_out[i].second.set_annotation("proper_pair", proper_paired);
            if (i != 0) {
                multipath_aln_pairs_out[i].first.set_annotation("secondary", true);
                multipath_aln_pairs_out[i].second.set_annotation("secondary", true);
            }
        }
        
        // remove the full length bonus if we don't want it in the final score
        if (strip_bonuses) {
            for (pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : multipath_aln_pairs_out) {
                strip_full_length_bonuses(multipath_aln_pair.first);
                strip_full_length_bonuses(multipath_aln_pair.second);
            }
        }
        
        // Compute the fragment length distribution.
        string distribution = "-I " + to_string(fragment_length_distr.mean()) + " -D " + to_string(fragment_length_distr.std_dev());
        
        for (pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : multipath_aln_pairs_out) {
            // Annotate with paired end distribution
            multipath_aln_pair.first.set_annotation("fragment_length_distribution", distribution);
            multipath_aln_pair.second.set_annotation("fragment_length_distribution", distribution);
        }
        
#ifdef debug_pretty_print_alignments
        cerr << "final alignments being returned:" << endl;
        for (const pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : multipath_aln_pairs_out) {
            cerr << "read 1: " << endl;
            view_multipath_alignment(cerr, multipath_aln_pair.first, *xindex);
            cerr << "read 2: " << endl;
            view_multipath_alignment(cerr, multipath_aln_pair.second, *xindex);
        }
#endif
        return proper_paired;
    }
    
    void MultipathMapper::reduce_to_single_path(const multipath_alignment_t& multipath_aln, vector<Alignment>& alns_out,
                                                size_t max_number) const {
    
#ifdef debug_multipath_mapper
        cerr << "linearizing multipath alignment to assess positional diversity" << endl;
#endif
        // Compute a few optimal alignments using disjoint sets of subpaths.
        // This hopefully gives us a feel for the positional diversity of the MultipathMapping.
        // But we still may have duplicates or overlaps in vg node space.
        auto alns = optimal_alignments_with_disjoint_subpaths(multipath_aln, max_number + 1);
        
        if (alns.empty()) {
            // This happens only if the read is totally unmapped
            assert(multipath_aln.subpath_size() == 0);
            
            // Output an unmapped alignment.
            alns_out.emplace_back();
            Alignment& aln = alns_out.back();
            
            // Transfer read information over to alignment
            transfer_read_metadata(multipath_aln, aln);
            
            // Score and MAPQ and path and stuff will all be 0.
            return;
        }
        
        // Otherwise we know there is at least one non-unmapped mapping
        
        // Make a list of all the scores
        vector<double> scores(1, alns[0].score());
        // Emit the alignment
        alns_out.push_back(alns[0]);
        
#ifdef debug_multipath_mapper
        cerr << "found optimal mapping with score " << alns_out[0].score() << endl;
        cerr << "\t" << pb2json(alns_out[0]) << endl;
#endif
        
        // Find all the nodes touched by the best alignment
        unordered_set<id_t> in_best;
        for (auto& m : alns[0].path().mapping()) {
            // Put each node in the set
            in_best.insert(m.position().node_id());
        }
        
        for (size_t i = 1; i < alns.size(); i++) {
            // For each other alignment, decide if it overlaps the best one
            size_t overlapped = 0;
            for (auto& m : alns[i].path().mapping()) {
                if (in_best.count(m.position().node_id())) {
                    overlapped++;
                }
            }
            
#ifdef debug_multipath_mapper
            cerr << "found suboptimal mapping overlapping " << overlapped << "/" << alns[i].path().mapping_size() << " with score "
                << alns[i].score() << endl;
            cerr << "\t" << pb2json(alns[i]) << endl;
#endif
            
            if (overlapped == 0) {
                // This is a nonoverlapping alignment so we want to emit it
                // Save its score
                scores.push_back(alns[i].score());
                // Emit the alignment
                alns_out.push_back(alns[i]);
                
                // Don't overlap with it either.
                for (auto& m : alns[i].path().mapping()) {
                    // Put each node in the set
                    in_best.insert(m.position().node_id());
                }
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "overall found optimal mapping with score " << alns_out[0].score() << " plus " << (alns_out.size() - 1)
            << " of " << max_number << " alternate linearizations";
        if (alns_out.size() >= 2) {
            cerr << " with best score " << alns_out[1].score();
        }
        cerr << endl;
#endif
       
        // Now compute the MAPQ for the best alignment
        auto placement_mapqs = compute_raw_mapping_qualities_from_scores(scores,
                                                                         !multipath_aln.quality().empty());
        // And min it in with what;s there already.
        alns_out[0].set_mapping_quality(min(alns_out[0].mapping_quality(), placement_mapqs.front()));
        for (size_t i = 1; i < alns_out.size(); i++) {
            // And zero all the others
            alns_out[i].set_mapping_quality(0);
        }
    }

    vector<pair<int64_t, int64_t>> MultipathMapper::covered_intervals(const Alignment& alignment,
                                                                      const clustergraph_t& cluster) const {
        
        // convert MEM subsequences to integer intervals
        vector<pair<int64_t, int64_t>> mem_intervals;
        mem_intervals.reserve(get<1>(cluster).first.size());
        for (const auto& hit : get<1>(cluster).first) {
            mem_intervals.emplace_back(hit.first->begin - alignment.sequence().begin(),
                                       hit.first->end - alignment.sequence().begin());
        }
        
        // put them in order
        sort(mem_intervals.begin(), mem_intervals.end());
        
        // do a sweep line
        vector<pair<int64_t, int64_t>> interval_union;
        int64_t begin, end;
        tie(begin, end) = mem_intervals.front();
        for (size_t i = 1; i < mem_intervals.size(); ++i) {
            if (mem_intervals[i].first > end) {
                interval_union.emplace_back(begin, end);
                tie(begin, end) = mem_intervals[i];
            }
            else {
                end = max(end, mem_intervals[i].second);
            }
        }
        interval_union.emplace_back(begin, end);
        return interval_union;
    }

    bool MultipathMapper::test_splice_candidates(const Alignment& alignment, bool searching_left,
                                                 multipath_alignment_t& anchor_mp_aln, double* anchor_multiplicity_out,
                                                 SpliceStrand& strand, int64_t num_candidates,
                                                 const function<const multipath_alignment_t&(int64_t)>& get_candidate,
                                                 const function<double(int64_t)>& get_multiplicity,
                                                 const function<multipath_alignment_t&&(int64_t)>& consume_candidate) const {
        
        
        
        /*
         * The region around a candidate's end, which could contain a splice junction
         */
        struct PrejoinSide {
            int64_t candidate_idx;
            SpliceRegion* splice_region;
            pos_t search_pos;
            int64_t clip_length;
            int32_t untrimmed_score;
        };
        
        /*
         * Two consistent pairs of identified splice sites with a joining alignment
         */
        struct PutativeJoin {
            PutativeJoin(const PathPositionHandleGraph& graph,
                         const SpliceStats& splice_stats, const Alignment& opt,
                         const GSSWAligner& aligner,
                         const PrejoinSide& left, const PrejoinSide& right,
                         const tuple<handle_t, size_t, int64_t>& left_location,
                         const tuple<handle_t, size_t, int64_t>& right_location,
                         size_t motif_idx)
                : joined_graph(graph, left.splice_region->get_subgraph(),
                               get<0>(left_location), get<1>(left_location),
                               right.splice_region->get_subgraph(),
                               get<0>(right_location), get<1>(right_location)),
                  left_search_dist(get<2>(left_location)),
                  right_search_dist(get<2>(right_location)),
                  left_clip_length(left.clip_length),
                  right_clip_length(right.clip_length),
                  left_candidate_idx(left.candidate_idx),
                  right_candidate_idx(right.candidate_idx),
                  estimated_intron_length(-1),
                  intron_score(0),
                  motif_idx(motif_idx),
                  untrimmed_score(left.untrimmed_score + right.untrimmed_score)
            {
                // memoize the best score
                max_score = pre_align_max_score(aligner, splice_stats, opt);
            }
            
            JoinedSpliceGraph joined_graph;
            int64_t left_search_dist;
            int64_t right_search_dist;
            int64_t left_clip_length;
            int64_t right_clip_length;
            size_t left_candidate_idx;
            size_t right_candidate_idx;
            int32_t max_score;
            int32_t untrimmed_score;
            size_t motif_idx;
            // intron stats start uninitialized until measuring length
            int32_t intron_score;
            int64_t estimated_intron_length;
            // these two filled out after doing alignment
            Alignment connecting_aln;
            size_t splice_idx;
            
            int32_t fixed_score_components(const SpliceStats& splice_stats,
                                           const Alignment& opt) {
                return splice_stats.motif_score(motif_idx) + untrimmed_score - opt.score();
            }
            
            void set_intron_length(int64_t estimated_intron_length,
                                   const SpliceStats& splice_stats) {
                
                estimated_intron_length = estimated_intron_length;
                intron_score = splice_stats.intron_length_score(estimated_intron_length);
                // memoize the max score again
                max_score += intron_score;
            }
            
            int32_t pre_align_max_score(const GSSWAligner& aligner,
                                        const SpliceStats& splice_stats,
                                        const Alignment& opt) {
                // compute a bound on the best possible score this join could get
                int64_t min_dist = joined_graph.min_link_length();
                int64_t max_dist = joined_graph.max_link_length();
                int32_t min_gap_penalty = 0;
                int64_t link_length = left_clip_length + right_clip_length - opt.sequence().size();
                if (link_length < min_dist) {
                    min_gap_penalty = aligner.score_gap(min_dist - link_length);
                }
                else if (link_length > max_dist) {
                    min_gap_penalty = aligner.score_gap(link_length - max_dist);
                }
                return (min_gap_penalty
                        + aligner.score_exact_match(opt, opt.sequence().size() - left_clip_length,
                                                    link_length)
                        + fixed_score_components(splice_stats, opt));
            }
            
            int32_t post_align_net_score(const SpliceStats& splice_stats,
                                         const Alignment& opt) {
                return fixed_score_components(splice_stats, opt) + connecting_aln.score() + intron_score;
            }
        };
        
        /*
         * Socket for iteration that either enumerates motifs exhaustively or by selecting the
         * the most promising pairs if there are too many
         */
        class MotifPairIterable {
        public:
            MotifPairIterable(size_t max_num_pairs,
                              const PrejoinSide& left_side,
                              const PrejoinSide& right_side,
                              size_t motif_num, size_t seq_len) :
                motif_num(motif_num),
                left_side(left_side),
                right_side(right_side),
                seq_len(seq_len)
            {
                const auto& left_sites = left_side.splice_region->candidate_splice_sites(motif_num);
                const auto& right_sites = right_side.splice_region->candidate_splice_sites(motif_num);
                
                if (left_sites.size() * right_sites.size() > max_num_pairs) {
#ifdef debug_multipath_mapper
                    cerr << "number of pairs " << left_sites.size() * right_sites.size() << " for motif " << motif_num << " is above maximum " << max_num_pairs << ", sampling downward" << endl;
#endif
                    
                    // there are too many pairs to just iterate over all of them, we have to select
                    // the most promising subset
                    
                    vector<size_t> left_idxs = range_vector(left_sites.size());
                    vector<size_t> right_idxs = range_vector(right_sites.size());
                    
                    // sort the indexes by min distance in opposite ordering
                    stable_sort(left_idxs.begin(), left_idxs.end(), [&](size_t i, size_t j) {
                        return get<2>(left_sites[i]) < get<2>(left_sites[j]);
                    });
                    stable_sort(right_idxs.begin(), right_idxs.end(), [&](size_t i, size_t j) {
                        return get<2>(right_sites[i]) > get<2>(right_sites[j]);
                    });
                    
                    // the difference from the ideal distance to the approximate graph distance, takes indexes
                    // in the left_idxs and right_idxs vectors
                    int64_t target_len = 2 * seq_len - left_side.clip_length - right_side.clip_length;
                    auto distance_diff = [&](size_t l, size_t r) {
                        return abs<int64_t>(get<2>(left_sites[left_idxs[l]]) + get<2>(right_sites[right_idxs[r]]) - target_len);
                    };
                    
                    // sweep to identify pairs that most nearly align
                    // records of (left idx, right idx, searching left?)
                    vector<tuple<size_t, size_t, bool>> nearest_idx;
                    for (size_t l = 0, r = 0; l < left_sites.size(); ++l) {
                        while (r + 1 < right_idxs.size() && distance_diff(l, r + 1) < distance_diff(l, r)) {
                            ++r;
                        }
                        nearest_idx.emplace_back(l, r, true);
                        if (r + 1 < right_idxs.size()) {
                            nearest_idx.emplace_back(l, r + 1, false);
                        }
                    }
                    
                    // opposite order to get smallest differences first
                    auto cmp = [&](const tuple<size_t, size_t, bool>& a, const tuple<size_t, size_t, bool>& b) {
                        return distance_diff(get<0>(a), get<1>(a)) > distance_diff(get<0>(b), get<1>(b));
                    };
                    
                    // preprocess to find the minimum diff quickly
                    make_heap(nearest_idx.begin(), nearest_idx.end(), cmp);
                    
                    while (idx_pairs.size() < max_num_pairs && !nearest_idx.empty()) {
                        auto front = nearest_idx.front();
                        idx_pairs.emplace_back(left_idxs[get<0>(front)], right_idxs[get<1>(front)]);
                        pop_heap(nearest_idx.begin(), nearest_idx.end(), cmp);
                        if (get<2>(front) && get<1>(front) > 0) {
                            // continue to the left
                            --get<1>(nearest_idx.back());
                            push_heap(nearest_idx.begin(), nearest_idx.end(), cmp);
                        }
                        else if (!get<2>(front) && get<1>(front) + 1 < right_idxs.size()) {
                            // continue to the right
                            ++get<1>(nearest_idx.back());
                            push_heap(nearest_idx.begin(), nearest_idx.end(), cmp);
                        }
                        else {
                            // we can't keep looking in this direction
                            nearest_idx.pop_back();
                        }
                    }
                }
            }
            
            /*
             * simple iterator interface
             */
            class iterator {
            public:
                iterator(const MotifPairIterable& iteratee, bool begin) : iteratee(iteratee) {
                    if (!begin) {
                        if (iteratee.idx_pairs.empty()) {
                            i = (iteratee.left_side.splice_region->candidate_splice_sites(iteratee.motif_num).size()
                                 * iteratee.right_side.splice_region->candidate_splice_sites(iteratee.motif_num).size());
                        }
                        else {
                            i = iteratee.idx_pairs.size();
                        }
                    }
                }
                
                iterator& operator++() {
                    ++i;
                    return *this;
                }
                
                bool operator==(const iterator& other) const {
                    // this isn't valid if they are looking at different motifs or have a different maximum,
                    // but that should never happen within the limited scope of this function
                    return (&iteratee == &other.iteratee && i == other.i);
                }
                
                bool operator!=(const iterator& other) const {
                    return !(*this == other);
                }
                
                pair<tuple<handle_t, size_t, int64_t>, tuple<handle_t, size_t, int64_t>> operator*() const {
                    if (iteratee.idx_pairs.empty()) {
                        size_t size = iteratee.left_side.splice_region->candidate_splice_sites(iteratee.motif_num).size();
                        size_t j = i % size;
                        size_t k = i / size;
                        return make_pair(iteratee.left_side.splice_region->candidate_splice_sites(iteratee.motif_num)[j],
                                         iteratee.right_side.splice_region->candidate_splice_sites(iteratee.motif_num)[k]);
                    }
                    else {
                        return make_pair(iteratee.left_side.splice_region->candidate_splice_sites(iteratee.motif_num)[iteratee.idx_pairs[i].first],
                                         iteratee.right_side.splice_region->candidate_splice_sites(iteratee.motif_num)[iteratee.idx_pairs[i].second]);
                    }
                }
                
            private:
                size_t i = 0;
                const MotifPairIterable& iteratee;
            };
            
            
            // iterable methods
            iterator begin() const {
                return iterator(*this, true);
            }
            iterator end() const {
                return iterator(*this, false);
            }
            
        private:
            
            friend class iterator;
            
            size_t motif_num;
            size_t seq_len;
            const PrejoinSide& left_side;
            const PrejoinSide& right_side;
            vector<pair<size_t, size_t>> idx_pairs;
        };
        
        if (num_candidates == 0) {
#ifdef debug_multipath_mapper
            cerr << "no splice candidate to attempt join with" << endl;
            return false;
#endif
        }
        
        // we'll memoize the relatively expensive reference distance computations, since often
        // there are multiple splice motifs on the same node
        unordered_map<tuple<nid_t, bool, nid_t, bool>, int64_t> ref_length_memo;
        auto get_reference_dist = [&](const pos_t& pos_1, const pos_t& pos_2) -> int64_t {
            
            tuple<nid_t, bool, nid_t, bool> key(id(pos_1), is_rev(pos_1), id(pos_2), is_rev(pos_2));
            
            auto it = ref_length_memo.find(key);
            if (it != ref_length_memo.end()) {
                // the reference distance of these nodes is already memoized
                return it->second - offset(pos_1) + offset(pos_2);
            }
            else {
                int64_t dist = numeric_limits<int64_t>::max();
                if (xindex->get_path_count() != 0) {
                    // estimate the distance using the reference path
                    dist = algorithms::ref_path_distance(xindex, pos_1, pos_2, ref_path_handles,
                                                         max_splice_ref_search_length);
                    
                }
                
//                if (distance_index && dist == numeric_limits<int64_t>::max()) {
//                    // FIXME: this will still sometimes produce finite distances for reads that
//                    // can't reach each other along a surjection path in cyclic graphs
//                    // FIXME: it can also find finite distances to different strands of a path,
//                    // but this might be okay sometimes?
//                    
//                    // they're probably still reachable if they got this far, get a worse estimate of the
//                    // distance from the distance index
//                    int64_t min_dist = distance_index->min_distance(pos_1, pos_2);
//                    if (min_dist >= 0) {
//                        dist = min_dist;
//                    }
//                }
                
                if (dist != numeric_limits<int64_t>::max()) {
                    // not memoizing unreachable distances, since distance index should
                    // filter out most of those anyway, and they actually might change on
                    // different positions on the node
                    ref_length_memo[key] = dist + offset(pos_1) - offset(pos_2);
                }
                return dist;
            }
        };
        
        
#ifdef debug_multipath_mapper
        cerr << "testing splice candidates for mp aln: " << endl;
        cerr << debug_string(anchor_mp_aln) << endl;
#endif
        
        vector<unique_ptr<SpliceRegion>> splice_regions;
        vector<PrejoinSide> left_prejoin_sides, right_prejoin_sides;
        
        vector<PrejoinSide>& anchor_prejoin_sides = searching_left ? right_prejoin_sides : left_prejoin_sides;
        vector<PrejoinSide>& candidate_prejoin_sides = searching_left ? left_prejoin_sides : right_prejoin_sides;
        
        // examine the region along the possible splice region of the anchor
        Alignment opt;
        optimal_alignment(anchor_mp_aln, opt);
        if (opt.path().mapping_size() == 0) {
            return false;
        }
        auto anchor_pos = trimmed_end(opt, max_splice_overhang, !searching_left, *xindex,
                                      *get_aligner(!opt.quality().empty()));
        
        splice_regions.emplace_back(new SpliceRegion(get<0>(anchor_pos), searching_left, 2 * max_splice_overhang,
                                                     *xindex, dinuc_machine, splice_stats));
        
        anchor_prejoin_sides.emplace_back();
        anchor_prejoin_sides.front().candidate_idx = -1;
        anchor_prejoin_sides.front().splice_region = splice_regions.front().get();
        anchor_prejoin_sides.front().search_pos = get<0>(anchor_pos);
        anchor_prejoin_sides.front().clip_length = get<1>(anchor_pos);
        anchor_prejoin_sides.front().untrimmed_score = opt.score() - get<2>(anchor_pos);
        
#ifdef debug_multipath_mapper
        cerr << "anchor stats:" << endl;
        cerr << "\tsearch pos " << anchor_prejoin_sides.front().search_pos << endl;
        cerr << "\tclip length " << anchor_prejoin_sides.front().clip_length << endl;
        cerr << "\topt score " << opt.score() << endl;
        cerr << "\tuntrimmed score " << anchor_prejoin_sides.front().untrimmed_score << endl;
#endif
                
        // examine the possible splice regions for the candidates
        bool found_splice_aln = false;
        for (int64_t i = 0; i < num_candidates; ++i) {
            
            auto& candidate = get_candidate(i);
            
            if (candidate.subpath().empty()) {
#ifdef debug_multipath_mapper
                cerr << "skipping empty candidate " << i << endl;
#endif
                continue;
            }
            
#ifdef debug_multipath_mapper
            cerr << "extracting splice region for candidate " << i << ":" << endl;
            cerr << debug_string(candidate) << endl;
#endif
            
            Alignment candidate_opt;
            optimal_alignment(candidate, candidate_opt);
            if (candidate_opt.path().mapping_size() == 0) {
#ifdef debug_multipath_mapper
                cerr << "skipping negative scoring candidate " << i << endl;
#endif
                continue;
            }
            
            auto candidate_pos = trimmed_end(candidate_opt, max_splice_overhang, searching_left, *xindex,
                                             *get_aligner(!opt.quality().empty()));
            
            splice_regions.emplace_back(new SpliceRegion(get<0>(candidate_pos), !searching_left, 2 * max_splice_overhang,
                                                         *xindex, dinuc_machine, splice_stats));
            
            candidate_prejoin_sides.emplace_back();
            auto& candidate_side = candidate_prejoin_sides.back();
            candidate_side.candidate_idx = i;
            candidate_side.splice_region = splice_regions.back().get();
            candidate_side.search_pos = get<0>(candidate_pos);
            candidate_side.clip_length = get<1>(candidate_pos);
            candidate_side.untrimmed_score = candidate_opt.score() - get<2>(candidate_pos);
            
#ifdef debug_multipath_mapper
            cerr << "candidate stats:" << endl;
            cerr << "\tsearch pos " << candidate_side.search_pos << endl;
            cerr << "\tclip length " << candidate_side.clip_length << endl;
            cerr << "\topt score " << candidate_opt.score() << endl;
            cerr << "\tuntrimmed score " << candidate_side.untrimmed_score << endl;
#endif
        }
        
        // identify the possible joins down to a base level, including the intron length
        vector<PutativeJoin> putative_joins;
        for (auto& left_prejoin_side : left_prejoin_sides) {
            
            auto& left_region = *left_prejoin_side.splice_region;
            auto left_seed = left_region.get_seed_pos();
            auto left_seed_under = left_region.get_subgraph().get_underlying_handle(left_seed.first);
            pos_t left_seed_pos(xindex->get_id(left_seed_under), xindex->get_is_reverse(left_seed_under), left_seed.second);
            
            for (auto& right_prejoin_side : right_prejoin_sides) {
                
#ifdef debug_multipath_mapper
                cerr << "resolving joins for left candidate " << left_prejoin_side.candidate_idx << ", right candidate " << right_prejoin_side.candidate_idx << endl;
#endif
                auto& right_region = *right_prejoin_side.splice_region;
                auto right_seed = right_region.get_seed_pos();
                auto right_seed_under = right_region.get_subgraph().get_underlying_handle(right_seed.first);
                pos_t right_seed_pos(xindex->get_id(right_seed_under), xindex->get_is_reverse(right_seed_under), right_seed.second);
                
                if (distance_index) {
                    // check if these regions can reach each other
                    size_t test_dist = minimum_distance(*distance_index, left_seed_pos, right_seed_pos, false, xindex);
#ifdef debug_multipath_mapper
                    cerr << "got distance index test distance " << test_dist << " between seed positions " << left_seed_pos << " and " << right_seed_pos << endl;
#endif
                    if (test_dist == numeric_limits<size_t>::max() || test_dist > max_intron_length) {
#ifdef debug_multipath_mapper
                        cerr << "test distance shows that this pair of candidates cannot reach each other" << endl;
#endif
                        continue;
                    }
                }
                
                // figure out how many pairs there are of each motif and across all motifs
                size_t total_num_pairs = 0;
                vector<size_t> num_candidate_pairs(splice_stats.motif_size(), 0);
                for (size_t j = 0; j < splice_stats.motif_size(); ++j) {
                    if (strand != Undetermined && splice_stats.motif_is_reverse(j) != (strand == Reverse)) {
                        // we can only find splicing at junctions that have a consistent strand
                        continue;
                    }
                    
                    num_candidate_pairs[j] = (left_region.candidate_splice_sites(j).size()
                                              * right_region.candidate_splice_sites(j).size());
                    total_num_pairs += num_candidate_pairs[j];
                    
#ifdef debug_multipath_mapper
                    cerr << "motif " << j << " will have num candidate pairs: " << left_region.candidate_splice_sites(j).size() << " * " << right_region.candidate_splice_sites(j).size() << " = " << num_candidate_pairs[j] << endl;
#endif
                }
                
                // apportion the effort we'll spend across motifs
                vector<size_t> motif_max_num_pairs;
                if (total_num_pairs < max_motif_pairs) {
                    // we can afford to do all of the candidates
                    motif_max_num_pairs = move(num_candidate_pairs);
                }
                else {
                    // we have to budget out the number of candidates, so we have this procedure
                    // to apportion them among motifs based the motif frequency
                    // note: each round has to either clear out the entire budget or clear out all
                    // instances of at least one motif pair, which guarantees eventual termination
                    // in `splice_stats.motif_size()` rounds
                    // (the rounding up can actually cause us to break the maximum, but it
                    // shouldn't be by very much, and rounding up is necessary to guarantee the
                    // clearing property)
                    
                    motif_max_num_pairs.resize(splice_stats.motif_size(), 0);
                    
                    int64_t budget_remaining = max_motif_pairs;
                    while (budget_remaining > 0) {
#ifdef debug_multipath_mapper
                        cerr << "new round of apportioning motif pair budget:" << endl;
#endif
                        // compute the normalizing factor for the frequency among the motifs that
                        // still have motifs to assign
                        double total_freq = 0.0;
                        for (size_t j = 0; j < splice_stats.motif_size(); ++j) {
                            if (motif_max_num_pairs[j] < num_candidate_pairs[j]) {
                                total_freq += splice_stats.motif_frequency(j);
                            }
                        }
                        int64_t next_budget_remaining = budget_remaining;
                        for (size_t j = 0; j < splice_stats.motif_size(); ++j) {
                            int64_t motif_remaining = num_candidate_pairs[j] - motif_max_num_pairs[j];
                            if (motif_remaining > 0) {
                                // figure out how many pairs out of the budget to give to this motif on this
                                // round
                                double round_fraction = splice_stats.motif_frequency(j) / total_freq;
                                int64_t round_pairs = min<int64_t>(ceil(round_fraction * budget_remaining),
                                                                   motif_remaining);
                                // give it out from the budget
                                next_budget_remaining -= round_pairs;
                                motif_max_num_pairs[j] += round_pairs;
#ifdef debug_multipath_mapper
                                cerr << "\tadd " << round_pairs << " to motif " << j << ": " << motif_max_num_pairs[j] << endl;
#endif
                            }
                            
                        }
                        budget_remaining = next_budget_remaining;
                    }
#ifdef debug_multipath_mapper
                    cerr << "final apportionment:" << endl;
                    for (size_t j = 0; j < splice_stats.motif_size(); ++j) {
                        
                        size_t left = left_region.candidate_splice_sites(j).size();
                        size_t right = right_region.candidate_splice_sites(j).size();
                        cerr << "motif " << j << ": " << motif_max_num_pairs[j]  << " pairs out of " << left * right << endl;
                    }
#endif
                }
                
                
                for (size_t j = 0; j < splice_stats.motif_size(); ++j) {
                    if (motif_max_num_pairs[j] == 0) {
                        // skip this one
                        continue;
                    }
                    
                    // let the iterable decide which motif pairs we should look at
                    MotifPairIterable motif_pairs(motif_max_num_pairs[j], left_prejoin_side,
                                                  right_prejoin_side, j, alignment.sequence().size());
                    
                    for (auto it = motif_pairs.begin(), end = motif_pairs.end(); it != end; ++it) {
                        
                        tuple<handle_t, size_t, int64_t> left_location, right_location;
                        tie(left_location, right_location) = *it;
                        
                        auto l_under = left_region.get_subgraph().get_underlying_handle(get<0>(left_location));
                        auto r_under = right_region.get_subgraph().get_underlying_handle(get<0>(right_location));
                        
                        pos_t l_pos(xindex->get_id(l_under), xindex->get_is_reverse(l_under), get<1>(left_location));
                        pos_t r_pos(xindex->get_id(r_under), xindex->get_is_reverse(r_under), get<1>(right_location));
                        
#ifdef debug_multipath_mapper
                        cerr << "\tchecking shared motif " << j << " with has positions " << l_pos << ", and " << r_pos << endl;
#endif
                        
                        putative_joins.emplace_back(*xindex, splice_stats, opt,
                                                    *get_aligner(!alignment.quality().empty()),
                                                    left_prejoin_side, right_prejoin_side,
                                                    left_location, right_location, j);
                        
                        if (putative_joins.back().max_score < no_splice_log_odds) {
#ifdef debug_multipath_mapper
                            cerr << "\tscore bound of " << putative_joins.back().max_score << " ensures insigificant spliced alignment against prior log odds " << no_splice_log_odds << " before measuring intron length" << endl;
#endif
                            
                            // this has no chance of becoming significant, let's skip it
                            putative_joins.pop_back();
                            continue;
                        }
                        
                        // measure the intron length
                        int64_t dist = get_reference_dist(l_pos, r_pos);
                        if (dist <= 0 || dist > max_intron_length || dist == numeric_limits<int64_t>::max()) {
#ifdef debug_multipath_mapper
                            cerr << "\tinconsistent intron length " << dist << ", skipping putative join" << endl;
#endif
                            putative_joins.pop_back();
                            continue;
                        }
                        
                        putative_joins.back().set_intron_length(dist, splice_stats);
                        
                        // TODO: enforce pairing constraints?
                        
#ifdef debug_multipath_mapper
                        cerr << "\tshared motif has a spliceable path of length " << dist << " (intron score: " << putative_joins.back().intron_score << "), adding as a putative join with score bound " << putative_joins.back().max_score << endl;
#endif
                        if (putative_joins.back().max_score < no_splice_log_odds) {
#ifdef debug_multipath_mapper
                            cerr << "\tscore bound of " << putative_joins.back().max_score << " ensures insigificant spliced alignment against prior log odds " << no_splice_log_odds << " before doing alignment" << endl;
#endif
                            
                            // this has no chance of becoming significant, let's skip it
                            putative_joins.pop_back();
                        }
                    }
                }
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "realigning across " << putative_joins.size() << " putative join regions" << endl;
#endif
        
        // TODO: allow multiple splices in a multipath alignment
        int32_t best_net_score = -1;
        int64_t best_intron_length = numeric_limits<int64_t>::max();
        unique_ptr<PutativeJoin> best_join;
        
        auto score_bound_comp = [](const PutativeJoin& join_1, const PutativeJoin& join_2) {
            return join_1.max_score < join_2.max_score;
        };
        
        // order the joins by the highest upper bound on net score
        make_heap(putative_joins.begin(), putative_joins.end(), score_bound_comp);
        
        while (!putative_joins.empty() && putative_joins.front().max_score >= best_net_score) {
            
            auto& join = putative_joins.front();
            
            
            size_t connect_begin = alignment.sequence().size() - join.left_clip_length;
            size_t connect_len = join.left_clip_length + join.right_clip_length - alignment.sequence().size();
            
            join.connecting_aln.set_sequence(alignment.sequence().substr(connect_begin, connect_len));
            if (!alignment.quality().empty()) {
                join.connecting_aln.set_quality(alignment.quality().substr(connect_begin, connect_len));
            }
            
            // TODO: multi alignment?
            auto alnr = get_aligner(!alignment.quality().empty());
            alnr->align_global_banded(join.connecting_aln, join.joined_graph, 1);
            
            // the total score of extending the anchor by the candidate
            int32_t net_score = join.post_align_net_score(splice_stats, opt);
            
#ifdef debug_multipath_mapper
            cerr << "next candidate spliced alignment with score bound " << join.max_score << " has net score " << net_score << " after realigning read interval " << connect_begin << ":" << (connect_begin + connect_len) << ", must get " << no_splice_log_odds << " for significance"  << endl;
#endif
            
            if (net_score > no_splice_log_odds) {
                // this is a statistically significant spliced alignment
                
                // find which mapping is immediately after the splice
                auto path = join.connecting_aln.mutable_path();
                auto splice_id = join.joined_graph.get_id(join.joined_graph.right_splice_node());
                join.splice_idx = 1;
                while (path->mapping(join.splice_idx).position().node_id() != splice_id) {
                    ++join.splice_idx;
                }
                
                // and translate into the original ID space
                join.joined_graph.translate_node_ids(*path);
                
                if (net_score > best_net_score ||
                    (net_score == best_net_score && join.estimated_intron_length < best_intron_length)) {
#ifdef debug_multipath_mapper
                    cerr << "this score from motif " << join.motif_idx << " and splice site " << path->mapping(join.splice_idx - 1).position().node_id() << (path->mapping(join.splice_idx - 1).position().is_reverse() ? "-" : "+") << " -> " << path->mapping(join.splice_idx).position().node_id() << (path->mapping(join.splice_idx).position().is_reverse() ? "-" : "+") <<  " is the best so far, beating previous best " << best_net_score << endl;
#endif
                    best_intron_length = join.estimated_intron_length;
                    best_net_score = net_score;
                    best_join = unique_ptr<PutativeJoin>(new PutativeJoin(move(join)));
                    
                }
            }
            
            // queue up the next join in the heap front
            pop_heap(putative_joins.begin(), putative_joins.end(), score_bound_comp);
            putative_joins.pop_back();
        }
        
#ifdef debug_multipath_mapper
        cerr << "pruned " << putative_joins.size() << " putative joins for having low score bounds" << endl;
#endif
        
        if (best_join.get() == nullptr) {
#ifdef debug_multipath_mapper
            cerr << "no splice candidates were statistically significant" << endl;
#endif
            return false;
        }
        
        // greedily fix the strand
        // TODO: ideally we'd probably try fixing it each way and see which is better
        strand = (splice_stats.motif_is_reverse(best_join->motif_idx) ? Reverse : Forward);

        *anchor_multiplicity_out = min<double>(get_multiplicity(best_join->left_candidate_idx),
                                               get_multiplicity(best_join->right_candidate_idx));
        anchor_mp_aln = fuse_spliced_alignments(alignment,
                                                consume_candidate(best_join->left_candidate_idx),
                                                consume_candidate(best_join->right_candidate_idx),
                                                alignment.sequence().size() - best_join->left_clip_length,
                                                best_join->connecting_aln, best_join->splice_idx,
                                                splice_stats.motif_score(best_join->motif_idx) + best_join->intron_score - no_splice_log_odds,
                                                *get_aligner(!alignment.quality().empty()), *xindex);
        
#ifdef debug_multipath_mapper
        cerr << "found significant splice join, fused mp aln:" << endl;
        cerr << debug_string(anchor_mp_aln) << endl;
#endif
        
#ifdef debug_validate_multipath_alignments
#ifdef debug_multipath_mapper
        cerr << "validating spliced alignment:" << endl;
        cerr << debug_string(anchor_mp_aln) << endl;
#endif
        if (!validate_multipath_alignment(anchor_mp_aln, *xindex)) {
            cerr << "### WARNING ###" << endl;
            cerr << "multipath alignment of read " << anchor_mp_aln.sequence() << " failed to validate" << endl;
        }
#endif
        
        return true;
    }

    void MultipathMapper::align_to_splice_candidates(const Alignment& alignment,
                                                     vector<clustergraph_t>& cluster_graphs,
                                                     const vector<MaximalExactMatch>& mems,
                                                     const vector<size_t>& cluster_candidates,
                                                     const vector<pair<const MaximalExactMatch*, pos_t>>& hit_candidates,
                                                     const pair<int64_t, int64_t>& primary_interval,
                                                     bool searching_left,
                                                     bool is_read_1,
                                                     unordered_map<candidate_id_t, pair<multipath_alignment_t, double>>& unaligned_candidate_bank,
                                                     vector<candidate_id_t>& candidates_out,
                                                     const match_fanouts_t* mem_fanouts) const {
        
        
        // make an alignment to the a candidate's cluster
        auto align_to_candidate = [&](clustergraph_t& cluster_graph,
                                      const candidate_id_t& candidate_id) {
            
            auto& banked_candidate = unaligned_candidate_bank[candidate_id];
            
            multipath_align(alignment, cluster_graph, banked_candidate.first, mem_fanouts);
            topologically_order_subpaths(banked_candidate.first);
            
            banked_candidate.second = cluster_multiplicity(get<1>(cluster_graph));
        };
        
        // get the candidate alignment associated with the ID, making the alignment if necessary
        auto get_candidate = [&](const candidate_id_t& candidate_id) -> multipath_alignment_t& {
            auto it = unaligned_candidate_bank.find(candidate_id);
            if (it != unaligned_candidate_bank.end()) {
                // the candidate is already in the bank
                return it->second.first;
            }
            else {
                if (get<1>(candidate_id) >= 0) {
                    // this is an unaligned cluster
                    align_to_candidate(cluster_graphs[get<1>(candidate_id)], candidate_id);
                }
                else {
                    // this is an unclustered hit
                    vector<memcluster_t> dummy_clusters(1);
                    dummy_clusters.front().first.emplace_back(make_pair(get<2>(candidate_id),
                                                                        get<3>(candidate_id)));
                    dummy_clusters.front().second = 1.0;
                    auto dummy_cluster_graphs = query_cluster_graphs(alignment, mems, dummy_clusters);
                    align_to_candidate(dummy_cluster_graphs.front(), candidate_id);
                }
                return unaligned_candidate_bank[candidate_id].first;
            }
        };
        
        // decide if a candidate should be evaluated as a potential splice
        auto consider_candidate = [&](const candidate_id_t& candidate_id) {
            
            const multipath_alignment_t& candidate = get_candidate(candidate_id);
            
            // TODO: repetitive with identify
            // check if the fully realized alignment still looks approx disjoint with the primary
            auto interval = aligned_interval(candidate);
            if (searching_left) {
                if (interval.second >= primary_interval.first + 2 * max_softclip_overlap ||
                    min<int64_t>(interval.second, primary_interval.first) - interval.first < min_softclip_length_for_splice) {
#ifdef debug_multipath_mapper
                    cerr << "rejecting candidate because of overlap" << endl;
                    cerr << "\tprimary interval: " << primary_interval.first << " " << primary_interval.second << endl;
                    cerr << "\tcandidate interval: " << interval.first << " " << interval.second << endl;
#endif
                    return;
                }
            }
            else {
                if (interval.first < primary_interval.second - 2 * max_softclip_overlap ||
                    interval.second - max<int64_t>(interval.first, primary_interval.second) < min_softclip_length_for_splice) {
#ifdef debug_multipath_mapper
                    cerr << "rejecting candidate because of overlap" << endl;
                    cerr << "\tprimary interval: " << primary_interval.first << " " << primary_interval.second << endl;
                    cerr << "\tcandidate interval: " << interval.first << " " << interval.second << endl;
#endif
                    return;
                }
            }
            // did not reject the candidate
            candidates_out.push_back(candidate_id);
            return;
        };
        
        // unaligned clusters
        for (auto i : cluster_candidates) {
#ifdef debug_multipath_mapper
            size_t num_before_align = candidates_out.size();
#endif
            
            // make the associated ID and evaluate it
            candidate_id_t candidate_id(is_read_1, i, nullptr, pos_t());
            consider_candidate(candidate_id);
            
#ifdef debug_multipath_mapper
            if (candidates_out.size() > num_before_align) {
                cerr << "added cluster splice candidate " << i << ":" << endl;
                cerr << debug_string(get_candidate(candidates_out.back())) << endl;
            }
#endif
        }
        
        // unclustered hits
        for (const auto& hit : hit_candidates) {
            
#ifdef debug_multipath_mapper
            size_t num_before_align = candidates_out.size();
#endif
            // make the associated ID and evaluate it
            candidate_id_t candidate_id(is_read_1, -1, hit.first, hit.second);
            consider_candidate(candidate_id);
            
#ifdef debug_multipath_mapper
            if (candidates_out.size() > num_before_align) {
                cerr << "made alignment to a seed splice candidate " << hit.first->sequence() << " " << hit.second << endl;
                cerr << debug_string(get_candidate(candidates_out.back())) << endl;
            }
#endif
        }
        
#ifdef debug_validate_multipath_alignments
        for (size_t i = 0; i < candidates_out.size(); ++i) {
            const auto& multipath_aln = get_candidate(candidates_out[i]);
#ifdef debug_multipath_mapper
            cerr << "validating " << i << "-th splice candidate alignment:" << endl;
            cerr << debug_string(multipath_aln) << endl;
#endif
            if (!validate_multipath_alignment(multipath_aln, *xindex)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln.sequence() << " failed to validate" << endl;
            }
        }
#endif
    }

    bool MultipathMapper::attempt_rescue_for_splice_segment(const Alignment& alignment, const pair<int64_t, int64_t>& primary_interval,
                                                            const multipath_alignment_t& rescue_anchor,
                                                            bool rescue_left, multipath_alignment_t& rescued_out) const {
        
        // extract the portion of the read that we want to form a spliced alignment with
        string::const_iterator begin, end;
        if (rescue_left) {
            begin = alignment.sequence().begin() + max<int64_t>(0, primary_interval.second - max_softclip_overlap);
            end = alignment.sequence().end();
        }
        else {
            begin = alignment.sequence().begin();
            end = alignment.sequence().begin() + min<int64_t>(primary_interval.first + max_softclip_overlap,
                                                              alignment.sequence().size());
        }
        Alignment splice_aln;
        splice_aln.set_sequence(string(begin, end));
        if (!alignment.quality().empty()) {
            splice_aln.set_quality(string(alignment.quality().begin() + (begin - alignment.sequence().begin()),
                                          alignment.quality().begin() + (end - alignment.sequence().begin())));
        }
        
#ifdef debug_multipath_mapper
        cerr << "attempting to rescue a spliced alignment segment to " << (rescue_left ? "left" : "right") << ":" << endl;
        cerr << pb2json(splice_aln) << endl;
#endif
        
        // adjust the mean length to account for the part we're not realigning
        double rescue_mean_length = max<double>(fragment_length_distr.mean()
                                                - (alignment.sequence().size() - splice_aln.sequence().size()), 0.0);
        
        // try to align
        bool succeeded = do_rescue_alignment(rescue_anchor, splice_aln, !rescue_left, rescued_out,
                                             rescue_mean_length, splice_rescue_graph_std_devs);
        
        if (!succeeded) {
            // we couldn't do the alignment
            return false;
        }
        
        // check if we got enough new, disjoint matches to be worth looking at (we can be
        // pretty permissive here because we will test this candidate for signficance in next
        // step of spliced alignment algorithm
        auto rescued_interval = aligned_interval(rescued_out);
        auto aligned_length = rescued_interval.second - rescued_interval.first;
        auto threshold = log(2.0 * splice_rescue_graph_std_devs * fragment_length_distr.std_dev()) / log(4.0);
        if (aligned_length < threshold) {
#ifdef debug_multipath_mapper
            cerr << "rescue candidate with aligned length " << aligned_length << " does not reach significance threshold of " << threshold << ", not checking for spliced alignments" << endl;
#endif
            return false;
        }
        if (rescue_left) {
            // we need to adjust offsets to make it match the full read
            rescued_interval.first += (begin - alignment.sequence().begin());
            rescued_interval.second += (begin - alignment.sequence().begin());
            
            succeeded = (rescued_interval.first >= primary_interval.second - max_softclip_overlap
                         && rescued_interval.second - max(rescued_interval.first, primary_interval.second) >= min_splice_rescue_matches
                         && rescued_interval.second > primary_interval.first);
        }
        else {
            succeeded = (rescued_interval.second <= primary_interval.first + max_softclip_overlap
                         && min(rescued_interval.second, primary_interval.first) >= min_splice_rescue_matches
                         && rescued_interval.first < primary_interval.first);
        }
#ifdef debug_multipath_mapper
        cerr << "rescue candidate covers read interval " << rescued_interval.first << ":" << rescued_interval.second << " compared to primary interval " << primary_interval.first << ":" << primary_interval.second << ", considered successful? " << succeeded << endl;
#endif
        
        if (succeeded) {
#ifdef debug_multipath_mapper
            cerr << "re-introducing trimmed sequence to rescued multipath alignment:" << endl;
            cerr << debug_string(rescued_out) << endl;
#endif
            
            // add in the soft-clips for the part of the read we trimmed off
            if (rescue_left) {
                // the softclip should come at the beginning of the rescued alignment
                rescued_out.mutable_sequence()->insert(rescued_out.mutable_sequence()->begin(),
                                                       alignment.sequence().begin(), begin);
                if (!alignment.quality().empty()) {
                    rescued_out.mutable_quality()->insert(rescued_out.mutable_quality()->begin(),
                                                          alignment.quality().begin(),
                                                          alignment.quality().begin() + (begin - alignment.sequence().begin()));
                }
                for (auto i : rescued_out.start()) {
                    auto subpath = rescued_out.mutable_subpath(i);
                    auto mapping = subpath->mutable_path()->mutable_mapping(0);
                    if (mapping->edit().front().from_length() == 0) {
                        // there's already a softclip, just expand it
                        auto edit = mapping->mutable_edit(0);
                        edit->mutable_sequence()->insert(edit->mutable_sequence()->begin(),
                                                         alignment.sequence().begin(), begin);
                        edit->set_to_length(edit->sequence().size());
                    }
                    else {
                        // add a new edit for the softclip
                        edit_t softclip;
                        softclip.set_sequence(string(alignment.sequence().begin(), begin));
                        softclip.set_to_length(softclip.sequence().size());
                        softclip.set_from_length(0);
                        mapping->mutable_edit()->insert(mapping->edit().begin(), softclip);
                        // we need to take away the full length bonus
                        subpath->set_score(subpath->score()
                                           - get_aligner(alignment.quality().empty())->score_full_length_bonus(true, splice_aln));
                    }
                }
            }
            else {
                // the softclip should come at the end of the rescued alignment
                *rescued_out.mutable_sequence() += string(end, alignment.sequence().end());
                if (!alignment.quality().empty()) {
                    *rescued_out.mutable_quality() += string(alignment.quality().begin() + (end - alignment.sequence().begin()),
                                                             alignment.quality().end());
                }
                for (size_t i = 0; i < rescued_out.subpath_size(); ++i) {
                    auto subpath = rescued_out.mutable_subpath(i);
                    if (subpath->next().empty()) {
                        auto& mapping = subpath->mutable_path()->mutable_mapping()->back();
                        if (mapping.edit().back().from_length() == 0) {
                            // expand existing softclip
                            auto& edit = mapping.mutable_edit()->back();
                            *edit.mutable_sequence() += string(end, alignment.sequence().end());
                            edit.set_to_length(edit.sequence().size());
                        }
                        else {
                            // add new edit for softclip
                            auto softclip = mapping.add_edit();
                            softclip->set_sequence(string(end, alignment.sequence().end()));
                            softclip->set_to_length(softclip->sequence().size());
                            softclip->set_from_length(0);
                            // we need to take away the full length bonus
                            subpath->set_score(subpath->score()
                                               - get_aligner(alignment.quality().empty())->score_full_length_bonus(false, splice_aln));
                        }
                    }
                }
            }
#ifdef debug_multipath_mapper
            cerr << "final rescued multipath alignment segment:" << endl;
            cerr << debug_string(rescued_out) << endl;
#endif
        }
        return succeeded;
    }

    void MultipathMapper::identify_aligned_splice_candidates(const Alignment& alignment, bool search_left,
                                                             const pair<int64_t, int64_t>& primary_interval,
                                                             const vector<multipath_alignment_t>& multipath_alns,
                                                             const vector<size_t>& cluster_idxs,
                                                             const vector<int64_t>& current_index, int64_t anchor,
                                                             unordered_set<size_t>& clusters_used_out,
                                                             vector<size_t>& mp_aln_candidates_out) const {
        
        // just to make sure we don't bother looking at earlier clusters (which have already had a chance to splice)
        for (size_t idx = 0; idx <= anchor; ++idx) {
            int64_t i = current_index[idx];
            if (i >= 0) {
                clusters_used_out.insert(cluster_idxs[i]);
            }
        }
        
        for (size_t idx = anchor + 1; idx < current_index.size(); ++idx) {
            
            int64_t i = current_index[idx];
            if (i < 0) {
                continue;
            }
            
            // check that the alignment is mostly disjoint of the primary and that
            // that the independent aligned portion is significant to call this a potential splice alignment
            auto interval = aligned_interval(multipath_alns[i]);
            
#ifdef debug_multipath_mapper
            cerr << "aligned candidate " << i << " has interval " << interval.first << " " << interval.second << endl;
#endif
            
            
            if (search_left) {
                if (interval.second < primary_interval.first + max_softclip_overlap &&
                    min<int64_t>(interval.second, primary_interval.first) - interval.first >= min_softclip_length_for_splice) {
                    
                    mp_aln_candidates_out.push_back(i);
                    clusters_used_out.insert(cluster_idxs[i]);
                }
            }
            else {
                if (interval.first >= primary_interval.second - max_softclip_overlap &&
                    interval.second - max<int64_t>(interval.first, primary_interval.second) >= min_softclip_length_for_splice) {
                    
                    mp_aln_candidates_out.push_back(i);
                    clusters_used_out.insert(cluster_idxs[i]);
                }
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "found fully aligned splice candidates:" << endl;
        for (auto i : mp_aln_candidates_out) {
            cerr << "\t" << i << endl;
        }
#endif
    }

    void MultipathMapper::identify_aligned_splice_candidates(const Alignment& alignment, bool read_1, bool search_left,
                                                             const pair<int64_t, int64_t>& primary_interval,
                                                             const vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs,
                                                             const vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                             const vector<int64_t>& current_index, int64_t anchor,
                                                             unordered_set<size_t>& clusters_used_out,
                                                             vector<size_t>& mp_aln_candidates_out) const {
        
        
        // just to make sure we don't bother looking at earlier clusters (which have already had a chance to splice)
        for (size_t idx = 0; idx <= anchor; ++idx) {
            int64_t i = current_index[idx];
            if (i >= 0) {
                clusters_used_out.insert(read_1 ? cluster_pairs[i].first.first : cluster_pairs[i].first.second);
            }
        }
        
        for (size_t idx = anchor + 1; idx < current_index.size(); ++idx) {
            
            int64_t i = current_index[idx];
            if (i < 0) {
                continue;
            }
            
            const multipath_alignment_t& mp_aln = read_1 ? multipath_aln_pairs[i].first : multipath_aln_pairs[i].second;
            
            // check that the alignment is mostly disjoint of the primary and that
            // that the independent aligned portion is significant to call this a potential splice alignment
            auto interval = aligned_interval(mp_aln);
            
#ifdef debug_multipath_mapper
            cerr << "aligned candidate " << i << " has interval " << interval.first << " " << interval.second << endl;
#endif
            if (search_left) {
                if (interval.second < primary_interval.first + max_softclip_overlap &&
                    min<int64_t>(interval.second, primary_interval.first) - interval.first >= min_softclip_length_for_splice) {
                    
                    size_t cluster = read_1 ? cluster_pairs[i].first.first : cluster_pairs[i].first.second;
                    if (!clusters_used_out.count(cluster)) {
                        mp_aln_candidates_out.push_back(i);
                        clusters_used_out.insert(cluster);
                    }
                }
            }
            else {
                if (interval.first >= primary_interval.second - max_softclip_overlap &&
                    interval.second - max<int64_t>(interval.first, primary_interval.second) >= min_softclip_length_for_splice) {
                    
                    size_t cluster = read_1 ? cluster_pairs[i].first.first : cluster_pairs[i].first.second;
                    if (!clusters_used_out.count(cluster)) {
                        mp_aln_candidates_out.push_back(i);
                        clusters_used_out.insert(cluster);
                    }
                }
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "found fully aligned splice candidates:" << endl;
        for (auto i : mp_aln_candidates_out) {
            cerr << "\t" << i << endl;
        }
#endif
        
    }

    void MultipathMapper::identify_unaligned_splice_candidates(const Alignment& alignment, bool search_left,
                                                               const pair<int64_t, int64_t>& primary_interval,
                                                               const vector<MaximalExactMatch>& mems,
                                                               const vector<clustergraph_t>& cluster_graphs,
                                                               const unordered_set<size_t>& clusters_already_used,
                                                               vector<size_t>& cluster_candidates_out,
                                                               vector<pair<const MaximalExactMatch*, pos_t>>& hit_candidates_out) const {
        
#ifdef debug_multipath_mapper
        cerr << "looking for unaligned splice candidates" << endl;
#endif
        
        for (size_t i = 0; i < cluster_graphs.size(); ++i) {
            
            if (clusters_already_used.count(i) || get<1>(cluster_graphs[i]).first.empty()) {
                continue;
            }
            
            auto intervals = covered_intervals(alignment, cluster_graphs[i]);
            
#ifdef debug_multipath_mapper
            cerr << "cluster candidate " << i << " has intervals" << endl;
            for (auto interval : intervals) {
                cerr << "\t" << interval.first << " " << interval.second << endl;
            }
#endif
            
            // check to make sure cluster doesn't overlap too much with the alignment
            // and also covers a sufficient part of the alignment's softclip
            
            if (search_left) {
                int64_t overlap = intervals.back().second - primary_interval.first;
                int64_t ind_cov = 0;
                for (size_t i = 0; i < intervals.size(); ++i) {
                    auto& interval = intervals[i];
                    if (interval.first >= primary_interval.first) {
                        break;
                    }
                    else if (interval.second >= primary_interval.first) {
                        ind_cov += primary_interval.first - interval.first;
                        break;
                    }
                    else {
                        ind_cov += interval.second - interval.first;
                    }
                }
                
                if (overlap < max_softclip_overlap && ind_cov >= min_softclip_length_for_splice) {
                    cluster_candidates_out.emplace_back(i);
                }
            }
            else {
                int64_t overlap = primary_interval.second - intervals.front().first;
                int64_t ind_cov = 0;
                for (int64_t i = intervals.size() - 1; i >= 0; --i) {
                    auto& interval = intervals[i];
                    if (interval.second <= primary_interval.second) {
                        break;
                    }
                    else if (interval.first <= primary_interval.second) {
                        ind_cov += interval.second - primary_interval.second;
                        break;
                    }
                    else {
                        ind_cov += interval.second - interval.first;
                    }
                }
                
                if (overlap < max_softclip_overlap && ind_cov >= min_softclip_length_for_splice) {
                    cluster_candidates_out.emplace_back(i);
                }
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "found unaligned cluster candidates:" << endl;
        for (auto i : cluster_candidates_out) {
            cerr << "\t" << i << endl;
        }
#endif

        // if we already processed a hit's cluster, we don't need to look at the hit
        unordered_set<pair<const MaximalExactMatch*, pos_t>> hits_already_used;
        for (const auto& i : cluster_candidates_out) {
            // make sure it isn't a phony cluster number from rescue
            if (i < cluster_graphs.size()) {
                for (const auto& hit : get<1>(cluster_graphs[i]).first) {
                    hits_already_used.insert(hit);
                }
            }
        }
        for (const auto& i : clusters_already_used) {
            // make sure it isn't a phony cluster number from rescue
            if (i < cluster_graphs.size()) {
                for (const auto& hit : get<1>(cluster_graphs[i]).first) {
                    hits_already_used.insert(hit);
                }
            }
        }
        
        // TODO: tie in mem fanouts?
        
        for (size_t i = 0; i < mems.size(); ++i) {
            
            const auto& mem = mems[i];
            
            if (mem.length() < min_softclip_length_for_splice) {
                continue;
            }
            if (search_left) {
                int64_t overlap = (mem.end - alignment.sequence().begin()) - primary_interval.first;
                int64_t ind_cov = (min<int64_t>(mem.end - alignment.sequence().begin(), primary_interval.first)
                                   - (mem.begin - alignment.sequence().begin()));
                
                if (overlap < max_softclip_overlap && ind_cov >= min_softclip_length_for_splice) {
                    for (auto gcsa_node : mem.nodes) {
                        pos_t pos = make_pos_t(gcsa_node);
                        if (!hits_already_used.count(make_pair(&mem, pos))) {
                            hit_candidates_out.emplace_back(&mem, pos);
                        }
                    }
                }
            }
            else {
                int64_t overlap = primary_interval.second - (mem.begin - alignment.sequence().begin());
                int64_t ind_cov = ((mem.end - alignment.sequence().begin())
                                   - max<int64_t>(mem.begin - alignment.sequence().begin(), primary_interval.second));
                
                if (overlap < max_softclip_overlap && ind_cov >= min_softclip_length_for_splice) {
                    for (auto gcsa_node : mem.nodes) {
                        pos_t pos = make_pos_t(gcsa_node);
                        if (!hits_already_used.count(make_pair(&mem, pos))) {
                            hit_candidates_out.emplace_back(&mem, pos);
                        }
                    }
                }
            }
        }
#ifdef debug_multipath_mapper
        cerr << "found unclustered hit candidates:" << endl;
        for (auto hit : hit_candidates_out) {
            cerr << "\t" << hit.first->sequence() << " " << hit.second << endl;
        }
#endif
    }

    bool MultipathMapper::find_spliced_alignments(const Alignment& alignment,
                                                  vector<multipath_alignment_t>& multipath_alns_out,
                                                  vector<double>& multiplicities,
                                                  vector<size_t>& cluster_idxs,
                                                  const vector<MaximalExactMatch>& mems,
                                                  vector<clustergraph_t>& cluster_graphs,
                                                  const match_fanouts_t* fanouts,
                                                  const multipath_alignment_t* rescue_anchor,
                                                  bool rescue_left,
                                                  double rescue_multiplicity) {
        
        if (multipath_alns_out.empty()) {
            return false;
        }
        
        // TODO: it would be better to use likelihoods rather than scores (esp. in paired)
        
        // choose the score cutoff based on the original unspliced mappings
        int32_t min_score_to_attempt = (optimal_alignment_score(multipath_alns_out.front())
                                        - get_aligner()->mapping_quality_score_diff(max_mapping_quality));
        
        vector<int64_t> current_index(multipath_alns_out.size(), 0);
        for (int64_t i = 1; i < multipath_alns_out.size(); ++i) {
            current_index[i] = i;
        }
        vector<int64_t> original_index = current_index;
        
        // a shared bank of candidates that we can re-use across alignments
        unordered_map<candidate_id_t, pair<multipath_alignment_t, double>> unaligned_candidate_bank;
        
        // we'll keep track of whether any spliced alignments succeeded
        bool any_splices = false;
        // so far we haven't restricted to splice motifs on any particular strand
        SpliceStrand strand = Undetermined;
        
        for (size_t j = 0; j < current_index.size(); ) {
            if (current_index[j] < 0) {
                // this alignment has been consumed as a splice candidate
                ++j;
                continue;
            }
            
#ifdef debug_multipath_mapper
            cerr << "deciding whether to look for spliced alignments on mp aln " << current_index[j] << endl;
#endif
            
            if (optimal_alignment_score(multipath_alns_out[current_index[j]]) < min_score_to_attempt) {
                // the rest of the alignments are too low-scoring to look at
                break;
            }
            
            auto interval = aligned_interval(multipath_alns_out[current_index[j]]);
            if (interval.first == interval.second) {
                // this anchor is unmapped
                ++j;
                continue;
            }
            // TODO: repetitive with paired version
            auto alnr = get_aligner(!alignment.quality().empty());
            int64_t left_max_score = (alnr->score_exact_match(alignment, 0, interval.first)
                                      + (interval.first == 0 ? 0 : alnr->score_full_length_bonus(true, alignment)));
            int64_t right_max_score = (alnr->score_exact_match(alignment, interval.second,
                                                               alignment.sequence().size() - interval.second)
                                       + (interval.second == alignment.sequence().size() ? 0 : alnr->score_full_length_bonus(false, alignment)));
            bool search_left = left_max_score >= min_softclipped_score_for_splice;
            bool search_right = right_max_score >= min_softclipped_score_for_splice;
            
            if (!(search_left || search_right)) {
#ifdef debug_multipath_mapper
                cerr << "soft clips are not sufficiently large to look for spliced alignment on interval " << interval.first << ":" << interval.second << " with max tail scores " << left_max_score << " and " << right_max_score << ", max score required: " << min_softclipped_score_for_splice << endl;
#endif
                ++j;
                continue;
            }
            
#ifdef debug_multipath_mapper
            cerr << "looking for spliced alignments, to left? " << search_left << ", to right? " << search_right << ", interval " << interval.first << ":" << interval.second << ", with max tail scores " << left_max_score << " and " << right_max_score << ", max score required: " << min_softclipped_score_for_splice << endl;
#endif
            
            bool found_splice_for_anchor = false;
            
            for (bool do_left : {true, false}) {
                if ((do_left && !search_left) || (!do_left && !search_right)) {
                    continue;
                }
                // try to figure out the softclip that we're seeing looks like it could be due to readthrough
                // of a adapter across the paired reads
                // TODO: magic number
                static const int64_t adapter_overlap_slosh = 2;
                if (!do_left && !read_1_adapter.empty()) {
                    size_t begin = max<int64_t>(0, interval.second - adapter_overlap_slosh);
                    size_t end = min<size_t>(begin + read_1_adapter.size() + 2 * adapter_overlap_slosh,
                                             alignment.sequence().size());
                    size_t pos = kmp_search(alignment.sequence().c_str() + begin, end - begin,
                                            read_1_adapter.c_str(), read_1_adapter.size(),
                                            read_1_adapter_lps);
                    if (pos != string::npos) {
                        // this softclip is bracketed by a known adapter sequence, it is much more likely
                        // that it should be adapter trimmed rather than meriting a spliced alignment
                        continue;
                    }
                }
                
                // move the anchor out of the vector to protect it from any shuffling that goes on
                multipath_alignment_t splice_anchor = move(multipath_alns_out[current_index[j]]);
                
                // identify which alignments, clusters, and hits could be part of a spliced alignment
                vector<size_t> mp_aln_candidates;
                unordered_set<size_t> clusters_used;
                vector<size_t> cluster_candidates;
                vector<pair<const MaximalExactMatch*, pos_t>> hit_candidates;
                identify_aligned_splice_candidates(alignment, do_left, interval, multipath_alns_out, cluster_idxs,
                                                   current_index, j, clusters_used, mp_aln_candidates);
                identify_unaligned_splice_candidates(alignment, do_left, interval, mems,
                                                     cluster_graphs, clusters_used, cluster_candidates,
                                                     hit_candidates);
                
                // make alignments for any unaligned candidates
                vector<candidate_id_t> unaligned_candidates;
                align_to_splice_candidates(alignment, cluster_graphs, mems, cluster_candidates, hit_candidates,
                                           interval, do_left, true, unaligned_candidate_bank, unaligned_candidates,
                                           fanouts);
                
                // set up functions for the spliced alignment decision process
                
                function<const multipath_alignment_t&(int64_t)> get_candidate = [&](int64_t i) -> const multipath_alignment_t& {
                    if (i < mp_aln_candidates.size()) {
                        return multipath_alns_out[mp_aln_candidates[i]];
                    }
                    else {
                        return unaligned_candidate_bank.at(unaligned_candidates[i - mp_aln_candidates.size()]).first;
                    }
                };
                
                function<double(int64_t)> get_multiplicity = [&](int64_t i) {
                    if (i < 0) {
                        return multiplicities[current_index[j]];
                    }
                    else if (i < mp_aln_candidates.size()) {
                        return multiplicities[mp_aln_candidates[i]];
                    }
                    else {
                        return unaligned_candidate_bank.at(unaligned_candidates[i - mp_aln_candidates.size()]).second;
                    }
                };
                
                // TODO: this solution is kinda ugly
                multipath_alignment_t tmp;
                function<multipath_alignment_t&&(int64_t)> consume_candidate = [&](int64_t i) -> multipath_alignment_t&& {
                    if (i < 0) {
                        // consume the anchor
                        return move(splice_anchor);
                    }
                    else if (i < mp_aln_candidates.size()) {
                        
#ifdef debug_multipath_mapper
                        cerr << "consuming alignment at index " << mp_aln_candidates[i] << endl;
#endif
                        
                        // pull the alignment out
                        tmp = move(multipath_alns_out[mp_aln_candidates[i]]);
                        
                        // replace it in the vectors and clear the final position
                        multipath_alns_out[mp_aln_candidates[i]] = move(multipath_alns_out.back());
                        multiplicities[mp_aln_candidates[i]] = multiplicities.back();
                        cluster_idxs[mp_aln_candidates[i]] = cluster_idxs.back();
                        multipath_alns_out.pop_back();
                        multiplicities.pop_back();
                        cluster_idxs.pop_back();
                        
                        // do the bookkeeping to track the original indexes
                        int64_t consuming_original_index = original_index[mp_aln_candidates[i]];
                        int64_t moving_original_index = original_index[multipath_alns_out.size()];
                        
                        current_index[moving_original_index] = mp_aln_candidates[i];
                        current_index[consuming_original_index] = -1;
                        
                        original_index[mp_aln_candidates[i]] = moving_original_index;
                        original_index[multipath_alns_out.size()] = -1;
                        
                        // TODO: also do bookkeeping on the clusters to claim the hits
                        
#ifdef debug_multipath_mapper
                        cerr << "current indexes are now:" << endl;
                        for (size_t i = 0; i < current_index.size(); ++i) {
                            cerr << "\t" << i << " " << current_index[i] << endl;
                        }
                        cerr << "original indexes are now:" << endl;
                        for (size_t i = 0; i < original_index.size(); ++i) {
                            cerr << "\t" << i << " " << original_index[i] << endl;
                        }
#endif
                        
                        return move(tmp);
                    }
                    else {
                        // remove the candidate from the bank
                        auto candidate_id = unaligned_candidates[i - mp_aln_candidates.size()];
                        tmp = move(unaligned_candidate_bank.at(candidate_id).first);
                        unaligned_candidate_bank.erase(candidate_id);
                        // and return it
                        return move(tmp);
                    }
                };
                
                double anchor_multiplicity = multiplicities[current_index[j]];
                bool did_splice = test_splice_candidates(alignment, do_left, splice_anchor, &anchor_multiplicity,
                                                         strand, mp_aln_candidates.size() + unaligned_candidates.size(),
                                                         get_candidate, get_multiplicity, consume_candidate);
                
                if (!did_splice && rescue_anchor && do_left != rescue_left) {
                    // we didn't find any splice junctions, but we might be able to rescue the spliced portion off
                    // of the other read
                    did_splice = find_rescuable_spliced_alignments(alignment, splice_anchor, anchor_multiplicity,
                                                                   strand, *rescue_anchor, rescue_multiplicity,
                                                                   rescue_left, interval);
                }
                
#ifdef debug_log_splice_align_stats
                string line = alignment.name() + '\t' + to_string(did_splice) + '\t' + to_string(interval.first) + '\t' + to_string(interval.second) + '\t' + to_string(do_left) + '\t' + to_string(mp_aln_candidates.size()) + '\t' + to_string(cluster_candidates.size()) + '\t' + to_string(hit_candidates.size()) + '\n';
#pragma omp critical
                cerr << line;
#endif
                
                if (did_splice) {
                    // we may need to update the multiplicity based on the splicing
                    multiplicities[current_index[j]] = anchor_multiplicity;
                }
                
                any_splices = any_splices || did_splice;
                found_splice_for_anchor = found_splice_for_anchor || did_splice;
                
                // move the alignment back now that all of the shuffling is finished
                multipath_alns_out[current_index[j]] = move(splice_anchor);
            }
            
            if (!found_splice_for_anchor) {
                // there's no more splicing to be found for this anchor alignment
                ++j;
            }
        }
        
        if (any_splices) {
            // we'll need to re-score and re-sort
            sort_and_compute_mapping_quality(multipath_alns_out,
                                             &cluster_idxs, &multiplicities);
        }
        
        return any_splices;
    }

    bool MultipathMapper::find_spliced_alignments(const Alignment& alignment1, const Alignment& alignment2,
                                                  vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                  vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                  vector<double>& pair_multiplicities,
                                                  const vector<MaximalExactMatch>& mems1, const vector<MaximalExactMatch>& mems2,
                                                  vector<clustergraph_t>& cluster_graphs1, vector<clustergraph_t>& cluster_graphs2,
                                                  const match_fanouts_t* fanouts) {
        
        if (multipath_aln_pairs_out.empty()) {
            return false;
        }
        
#ifdef debug_check_adapters
        bool attempt_1_left = false, attempt_1_right = false, attempt_2_left = false, attempt_2_right = false;
        bool success_1_left = false, success_1_right = false, success_2_left = false, success_2_right = false;
        bool adapter_1_left = false, adapter_1_right = false, adapter_2_left = false, adapter_2_right = false;
#endif
        
        // TODO: it would be better to use likelihoods rather than scores (esp. in paired)
        
        // choose the score cutoff based on the original unspliced mappings
        int32_t min_score_to_attempt_1 = (optimal_alignment_score(multipath_aln_pairs_out.front().first)
                                          - get_aligner()->mapping_quality_score_diff(max_mapping_quality));
        int32_t min_score_to_attempt_2 = (optimal_alignment_score(multipath_aln_pairs_out.front().second)
                                          - get_aligner()->mapping_quality_score_diff(max_mapping_quality));
        
        vector<int64_t> current_index(multipath_aln_pairs_out.size(), 0);
        for (int64_t i = 1; i < multipath_aln_pairs_out.size(); ++i) {
            current_index[i] = i;
        }
        vector<int64_t> original_index = current_index;
        
        // a shared bank of candidates that we can re-use across alignments
        unordered_map<candidate_id_t, pair<multipath_alignment_t, double>> unaligned_candidate_bank;
        
        // we'll keep track of whether any spliced alignments succeeded
        bool any_splices = false;
        
        // so far we haven't restricted to splice motifs on any particular strand
        SpliceStrand strand = Undetermined;
        
        for (size_t j = 0; j < current_index.size(); ++j) {
            if (current_index[j] < 0) {
                // this alignment has been consumed as a splice candidate
                continue;
            }
            
            if (optimal_alignment_score(multipath_aln_pairs_out[current_index[j]].first) < min_score_to_attempt_1
                && optimal_alignment_score(multipath_aln_pairs_out[current_index[j]].second) < min_score_to_attempt_2) {
                // the rest of the alignments are too low scoring to consider
                break;
            }
            
#ifdef debug_multipath_mapper
            cerr << "determining whether to make spliced alignment for pair at index " << current_index[j] << endl;
#endif
            
#ifdef debug_check_adapters
            bool attempt_1_left = false, attempt_1_right = false, attempt_2_left = false, attempt_2_right = false;
            bool success_1_left = false, success_1_right = false, success_2_left = false, success_2_right = false;
#endif
            
            for (int read_num = 0; read_num < 2; ) {
                
                bool do_read_1 = (read_num == 0);
                
                // decide if this read looks like it could benefit from a spliced alignment
                auto interval = aligned_interval(do_read_1 ? multipath_aln_pairs_out[current_index[j]].first
                                                           : multipath_aln_pairs_out[current_index[j]].second);
                if (interval.first == interval.second) {
                    // this anchor is unmapped
                    read_num++;
                    continue;
                }
                
                auto alnr = get_aligner(!alignment1.quality().empty() && !alignment2.quality().empty());
                
                const Alignment& aln = do_read_1 ? alignment1 : alignment2;
                
                int64_t left_max_score = (alnr->score_exact_match(aln, 0, interval.first)
                                          + (interval.first == 0 ? 0 : alnr->score_full_length_bonus(true, aln)));
                int64_t right_max_score = (alnr->score_exact_match(aln, interval.second,
                                                                   aln.sequence().size() - interval.second)
                                           + (interval.second == aln.sequence().size() ? 0 : alnr->score_full_length_bonus(false, aln)));
                bool search_left = left_max_score >= min_softclipped_score_for_splice;
                bool search_right = right_max_score >= min_softclipped_score_for_splice;
                
#ifdef debug_check_adapters
                if (do_read_1) {
                    attempt_1_left |= search_left;
                    attempt_1_right |= search_right;
                }
                else {
                    attempt_2_left |= search_left;
                    attempt_2_right |= search_right;
                }
#endif
#ifdef debug_multipath_mapper
                cerr << "on read " << (do_read_1 ? 1 : 2) << " looking for spliced alignments, to left? " << search_left << ", to right? " << search_right << ", interval " << interval.first << ":" << interval.second << " with max tail scores " << left_max_score << " " << right_max_score << ", and required max score " << min_softclipped_score_for_splice << endl;
#endif
                
                bool found_splice_for_anchor = false;
                
                for (bool do_left : {true, false}) {
                    if ((do_left && !search_left) || (!do_left && !search_right)) {
                        continue;
                    }
                    // try to figure out the softclip that we're seeing looks like it could be due to readthrough
                    // of a adapter across the paired reads
                    // TODO: magic number
                    static const int64_t adapter_overlap_slosh = 2;
                    if (do_read_1 && !do_left && !read_1_adapter.empty()) {
                        size_t begin = max<int64_t>(0, interval.second - adapter_overlap_slosh);
                        size_t end = min<size_t>(begin + read_1_adapter.size() + 2 * adapter_overlap_slosh,
                                                 aln.sequence().size());
                        size_t pos = kmp_search(aln.sequence().c_str() + begin, end - begin,
                                                read_1_adapter.c_str(), read_1_adapter.size(),
                                                read_1_adapter_lps);
                        if (pos != string::npos) {
                            // this softclip is bracketed by a known adapter sequence, it is much more likely
                            // that it should be adapter trimmed rather than meriting a spliced alignment
                            continue;
                        }
                    }
                    if (!do_read_1 && do_left && !read_2_adapter.empty()) {
                        size_t end = min<size_t>(interval.first + adapter_overlap_slosh, aln.sequence().size());
                        size_t begin = max<int64_t>(0, end - read_2_adapter.size() - 2 * adapter_overlap_slosh);
                        size_t pos = kmp_search(aln.sequence().c_str() + begin, end - begin,
                                                read_2_adapter.c_str(), read_2_adapter.size(),
                                                read_2_adapter_lps);
                        if (pos != string::npos) {
                            // this softclip is bracketed by a known adapter sequence, it is much more likely
                            // that it should be adapter trimmed rather than meriting a spliced alignment
                            continue;
                        }
                    }
                    
                    // select the candidate sources for the corresponding read
                    // note: we move the whole anchor to protect it from any shuffling within the vector
                    // so we will have to move it back at the end of this iteration
                    multipath_alignment_t anchor_mp_aln;
                    const vector<MaximalExactMatch>* mems;
                    vector<clustergraph_t>* cluster_graphs;
                    if (do_read_1) {
                        anchor_mp_aln = move(multipath_aln_pairs_out[current_index[j]].first);
                        mems = &mems1;
                        cluster_graphs = &cluster_graphs1;
                    }
                    else {
                        anchor_mp_aln = move(multipath_aln_pairs_out[current_index[j]].second);
                        mems = &mems2;
                        cluster_graphs = &cluster_graphs2;
                    }
                    
                    // identify the splice candidate
                    vector<size_t> mp_aln_candidates;
                    unordered_set<size_t> clusters_used;
                    vector<size_t> cluster_candidates;
                    vector<pair<const MaximalExactMatch*, pos_t>> hit_candidates;
                    identify_aligned_splice_candidates(aln, do_read_1, do_left, interval, multipath_aln_pairs_out, cluster_pairs,
                                                       current_index, j, clusters_used, mp_aln_candidates);
                    identify_unaligned_splice_candidates(aln, do_left, interval, *mems,
                                                         *cluster_graphs, clusters_used, cluster_candidates,
                                                         hit_candidates);
                    
                    // align splice candidates that haven't been aligned yet
                    vector<candidate_id_t> unaligned_candidates;
                    align_to_splice_candidates(aln, *cluster_graphs, *mems, cluster_candidates, hit_candidates,
                                               interval, do_left, do_read_1, unaligned_candidate_bank, unaligned_candidates,
                                               fanouts);
                    
                    
                    function<const multipath_alignment_t&(int64_t)> get_candidate = [&](int64_t i) -> const multipath_alignment_t& {
                        if (i < mp_aln_candidates.size()) {
                            if (do_read_1) {
                                return multipath_aln_pairs_out[mp_aln_candidates[i]].first;
                            }
                            else {
                                return multipath_aln_pairs_out[mp_aln_candidates[i]].second;
                            }
                        }
                        else {
                            return unaligned_candidate_bank.at(unaligned_candidates[i - mp_aln_candidates.size()]).first;
                        }
                    };
                    
                    function<double (int64_t)> get_multiplicity = [&](int64_t i) -> double {
                        if (i < 0) {
                            return pair_multiplicities[current_index[j]];
                        }
                        else if (i < mp_aln_candidates.size()) {
                            return pair_multiplicities[i];
                        }
                        else {
                            return unaligned_candidate_bank.at(unaligned_candidates[i - mp_aln_candidates.size()]).second;
                        }
                    };
                    
                    // TODO: this solution is kinda ugly
                    multipath_alignment_t tmp;
                    function<multipath_alignment_t&&(int64_t)> consume_candidate = [&](int64_t i) -> multipath_alignment_t&& {
                        if (i < 0) {
#ifdef debug_multipath_mapper
                            cerr << "consuming anchor for read " << (do_read_1 ? 1 : 2) << " in pair at current index " << current_index[j] << " and original index " << original_index[current_index[j]] << endl;
#endif
                            // consume the anchor
                            return move(anchor_mp_aln);
                        }
                        else if (i < mp_aln_candidates.size()) {
#ifdef debug_multipath_mapper
                            cerr << "consuming read " << (do_read_1 ? 1 : 2) << " of pair " << mp_aln_candidates[i] << endl;
#endif
                            
                            // look to see if the opposite side of this pair exists in multiple pairs
                            size_t opposite_cluster = do_read_1 ? cluster_pairs[mp_aln_candidates[i]].first.second
                                                                : cluster_pairs[mp_aln_candidates[i]].first.first;
                            bool opposite_duplicated = false;
                            if (opposite_cluster != RESCUED) {
                                // we don't want to count rescued alignments as the same cluster as each other
                                for (size_t j = 0; j < cluster_pairs.size() && !opposite_duplicated; ++j) {
                                    opposite_duplicated = (j != mp_aln_candidates[i] &&
                                                           ((do_read_1 && cluster_pairs[j].first.second == opposite_cluster)
                                                            || (!do_read_1 && cluster_pairs[j].first.first == opposite_cluster)));
                                }
                            }
                            
                            if (opposite_duplicated) {
                                // the other side will continue to live in its other pair, so we can scavenge
                                // this multipath alignment
#ifdef debug_multipath_mapper
                                cerr << "the opposite side is duplicated, removing pair" << endl;
                                cerr << "pair is index " << i << " among candidates, which is current index " << mp_aln_candidates[i] << " and original index " << original_index[mp_aln_candidates[i]] << endl;
                                cerr << "will swap with current index " << multipath_aln_pairs_out.size() - 1 << ", which has original index " << original_index[multipath_aln_pairs_out.size() - 1] << endl;
#endif
                                
                                tmp = do_read_1 ? move(multipath_aln_pairs_out[mp_aln_candidates[i]].first)
                                                : move(multipath_aln_pairs_out[mp_aln_candidates[i]].second);
                                
                                // replace it in the vectors and clear the final position
                                multipath_aln_pairs_out[mp_aln_candidates[i]] = move(multipath_aln_pairs_out.back());
                                pair_multiplicities[mp_aln_candidates[i]] = pair_multiplicities.back();
                                cluster_pairs[mp_aln_candidates[i]] = cluster_pairs.back();
                                multipath_aln_pairs_out.pop_back();
                                pair_multiplicities.pop_back();
                                cluster_pairs.pop_back();
                                
                                // do the bookkeeping to track the original indexes
                                int64_t consuming_original_index = original_index[mp_aln_candidates[i]];
                                int64_t moving_original_index = original_index[multipath_aln_pairs_out.size()];
                                
                                current_index[moving_original_index] = mp_aln_candidates[i];
                                current_index[consuming_original_index] = -1;
                                
                                original_index[mp_aln_candidates[i]] = moving_original_index;
                                original_index[multipath_aln_pairs_out.size()] = -1;
                                // TODO: also do bookkeeping on the clusters to claim the hits
                            }
                            else {
                                // we don't want to mess up pairs, so just copy it out
                                tmp = do_read_1 ? multipath_aln_pairs_out[mp_aln_candidates[i]].first
                                                : multipath_aln_pairs_out[mp_aln_candidates[i]].second;
                            }
                            
#ifdef debug_multipath_mapper
                            cerr << "pair current indexes are now:" << endl;
                            for (size_t i = 0; i < current_index.size(); ++i) {
                                cerr << "\t" << i << " " << current_index[i] << endl;
                            }
                            cerr << "pair original indexes are now:" << endl;
                            for (size_t i = 0; i < original_index.size(); ++i) {
                                cerr << "\t" << i << " " << original_index[i] << endl;
                            }
                            cerr << "pair order:" << endl;
                            for (auto& mp_aln_pair : multipath_aln_pairs_out) {
                                pos_t p1, p2;
                                Alignment aln1, aln2;
                                optimal_alignment(mp_aln_pair.first, aln1);
                                optimal_alignment(mp_aln_pair.second, aln2);
                                if (aln1.path().mapping_size() != 0) {
                                    p1 = make_pos_t(aln1.path().mapping(0).position());
                                }
                                if (aln2.path().mapping_size() != 0) {
                                    p2 = make_pos_t(aln2.path().mapping(0).position());
                                }
                                cerr << "\t" << p1 << " " << p2 << endl;
                            }
#endif
                            
                            return move(tmp);
                        }
                        else {
                            // remove the candidate from the bank
                            auto candidate_id = unaligned_candidates[i - mp_aln_candidates.size()];
                            tmp = move(unaligned_candidate_bank.at(candidate_id).first);
                            unaligned_candidate_bank.erase(candidate_id);
                            // and return it
                            return move(tmp);
                        }
                    };
                    
                    // see if we can actually make spliced alignments
                    double anchor_multiplicity = pair_multiplicities[current_index[j]];
                    bool spliced_side = test_splice_candidates(aln, do_left, anchor_mp_aln, &anchor_multiplicity,
                                                               strand, mp_aln_candidates.size() + unaligned_candidates.size(),
                                                               get_candidate, get_multiplicity, consume_candidate);
                    
                    if (!spliced_side && do_read_1 != do_left) {
                        // we might be able to rescue a spliced alignment segment
                        const auto& rescue_anchor = do_read_1 ? multipath_aln_pairs_out[current_index[j]].second
                                                              : multipath_aln_pairs_out[current_index[j]].first;
                        spliced_side = find_rescuable_spliced_alignments(aln, anchor_mp_aln, anchor_multiplicity,
                                                                         strand, rescue_anchor, anchor_multiplicity,
                                                                         do_read_1, interval);
                    }
                    
                    if (spliced_side) {
                        // we may need to update the pair's multiplicity
                        pair_multiplicities[current_index[j]] = anchor_multiplicity;
                    }
                    
#ifdef debug_check_adapters
                    if (do_read_1) {
                        if (do_left) {
                            success_1_left |= spliced_side;
                        }
                        else {
                            success_1_right |= spliced_side;
                        }
                    }
                    else {
                        if (do_left) {
                            success_2_left |= spliced_side;
                        }
                        else {
                            success_2_right |= spliced_side;
                        }
                    }
#endif
                    
                    any_splices = any_splices || spliced_side;
                    found_splice_for_anchor = found_splice_for_anchor || spliced_side;
                    
                    // move the anchor back now that we've done all of the shuffling we were going to do
                    if (do_read_1) {
                        multipath_aln_pairs_out[current_index[j]].first = move(anchor_mp_aln);
                    }
                    else {
                        multipath_aln_pairs_out[current_index[j]].second = move(anchor_mp_aln);
                    }
                }
                
                if (!found_splice_for_anchor) {
                    // there are no more splices to be found for this read end
                    ++read_num;
                }
            }
            
#ifdef debug_check_adapters
            vector<int32_t> adapter_aln_scores;
            vector<string> adapters{"AGATCGGAAGAG"}; // the illumina universal adapter
            adapters.push_back(reverse_complement(adapters.back()));
            for (auto mp_aln : {multipath_aln_pairs_out[current_index[j]].first, multipath_aln_pairs_out[current_index[j]].second}) {
                for (auto adapter : adapters) {
                    adapter_aln_scores.push_back(SSWAligner().align(mp_aln.sequence(), adapter).score());
                }
            }
            string line = (alignment1.name() +
                           "\t" + to_string(attempt_1_left) + "\t" + to_string(attempt_1_right) +
                           "\t" + to_string(attempt_2_left) + "\t" + to_string(attempt_2_right) +
                           "\t" + to_string(success_1_left) + "\t" + to_string(success_1_right) +
                           "\t" + to_string(success_2_left) + "\t" + to_string(success_2_right));
            for (auto s : adapter_aln_scores) {
                line += (string("\t") + to_string(s));
            }
            line += "\n";
            cerr << line;
#endif
        }
        
        if (any_splices) {
            sort_and_compute_mapping_quality(multipath_aln_pairs_out, cluster_pairs, nullptr, &pair_multiplicities);
        }
        
        return any_splices;
    }

    bool MultipathMapper::find_rescuable_spliced_alignments(const Alignment& alignment,
                                                            multipath_alignment_t& splice_anchor,
                                                            double& anchor_multiplicity,
                                                            SpliceStrand& strand,
                                                            const multipath_alignment_t& rescue_anchor,
                                                            double rescue_multiplicity,
                                                            bool rescue_left,
                                                            const pair<int64_t, int64_t>& primary_interval) const {
        
        auto rescue_interval = aligned_interval(rescue_anchor);
        
        // TODO: it's possible that the original estimated mapping qualities aren't reflective
        // of the score ordering anymore now that we may have added spliced alignments, but this
        // should be a rare problem and we want to be sparing with rescues
        if (rescue_anchor.mapping_quality() >= min<int>(30, max_mapping_quality)
            && ((rescue_left && rescue_interval.first == 0) ||
                (!rescue_left && rescue_interval.second == rescue_anchor.sequence().size()))) {
            
            multipath_alignment_t rescued;
            bool succeeded = attempt_rescue_for_splice_segment(alignment, primary_interval, rescue_anchor,
                                                               rescue_left, rescued);
            
            if (succeeded) {
                // set up simple functions to provide the rescued alignment as a candidate
                function<const multipath_alignment_t&(int64_t)> get_rescued_candidate = [&](int64_t i) -> const multipath_alignment_t& {
                    if (i < 0) {
                        // TODO: is this branch actually necessary?
                        return splice_anchor;
                    }
                    else {
                        return rescued;
                    }
                };
                function<double(int64_t)> get_rescued_multiplicity = [&](int64_t i) {
                    if (i < 0) {
                        return anchor_multiplicity;
                    }
                    else {
                        return rescue_multiplicity;
                    }
                };
                function<multipath_alignment_t&&(int64_t)> consume_rescued = [&](int64_t i) -> multipath_alignment_t&& {
                    if (i < 0) {
                        // consume the anchor
                        return move(splice_anchor);
                    }
                    else {
                        return move(rescued);
                    }
                };
#ifdef debug_multipath_mapper
                cerr << "testing rescued spliced alignment segment as a candidate" << endl;
#endif
                
                // do a test again with only the rescued candidate
                return test_splice_candidates(alignment, !rescue_left, splice_anchor, &anchor_multiplicity, strand, 1,
                                              get_rescued_candidate, get_rescued_multiplicity, consume_rescued);
                
            }
        }
        return false;
    }

    bool MultipathMapper::retry_pairing_spliced_alignments(const Alignment& alignment1, const Alignment& alignment2,
                                                           vector<multipath_alignment_t>& multipath_alns_1,
                                                           vector<multipath_alignment_t>& multipath_alns_2,
                                                           const vector<size_t>& cluster_idxs_1,
                                                           const vector<size_t>& cluster_idxs_2,
                                                           const vector<double>& multiplicities_1,
                                                           const vector<double>& multiplicities_2,
                                                           vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                           vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs_out,
                                                           vector<double>& pair_multiplicities_out) const {
        
#ifdef debug_multipath_mapper
        cerr << "trying to re-pair spliced alignments" << endl;
#endif
        
        MemoizingGraph memoizing_graph(xindex);
        unique_ptr<OrientedDistanceMeasurer> distance_measurer = get_distance_measurer(memoizing_graph);
        
        // we want to restrict our attention to alignments that have been spliced
        vector<bool> is_spliced_1(multipath_alns_1.size()), is_spliced_2(multipath_alns_2.size());
        for (size_t i = 0; i < multipath_alns_1.size(); ++i) {
            is_spliced_1[i] = contains_connection(multipath_alns_1[i]);
        }
        for (size_t i = 0; i < multipath_alns_2.size(); ++i) {
            is_spliced_2[i] = contains_connection(multipath_alns_2[i]);
        }
        
        // TODO: this could all be made more efficient by not doing these computations
        // multiple times
        bool found_consistent = false;
        auto attempt_to_pair = [&](size_t i, size_t j) {

            if (multipath_alns_1[i].subpath_size() == 0 || multipath_alns_2[j].subpath_size() == 0) {
                return;
            }
            
            Alignment opt_1, opt_2;
            optimal_alignment(multipath_alns_1[i], opt_1);
            optimal_alignment(multipath_alns_2[j], opt_2);
            if (opt_1.path().mapping_size() == 0 || opt_2.path().mapping_size() == 0) {
                return;
            }
            
            pos_t inner_pos_1 = final_position(opt_1.path());
            int64_t aligned_length_1 = path_from_length(opt_1.path());
            
            pos_t inner_pos_2 = initial_position(opt_2.path());
            int64_t aligned_length_2 = path_from_length(opt_2.path());

#ifdef debug_multipath_mapper
            cerr << "trying to re-pair alns " << i << " and " << j << " with inner positions " << inner_pos_1 << " and " << inner_pos_2 << ", and aligned lengths " << aligned_length_1 << " and " << aligned_length_2 << endl;
#endif
            if (aligned_length_1 == 0 || aligned_length_2 == 0) {
                return;
            }
            
            int64_t dist = distance_measurer->oriented_distance(inner_pos_1, inner_pos_2);
            if (dist != numeric_limits<int64_t>::max()) {
                int64_t total_dist = dist + aligned_length_1 + aligned_length_2;
#ifdef debug_multipath_mapper
                cerr << "re-estimated disatnce: " << total_dist << endl;
#endif
                if (is_consistent(total_dist)) {
                    // note: we're kind of abusing cluster pairs here by temporarily making it
                    // point to alignments instead of clusters
                    cluster_pairs_out.emplace_back(make_pair(i, j), total_dist);
                    found_consistent = true;
#ifdef debug_multipath_mapper
                    cerr << "re-pair succeeded" << endl;
#endif
                }
            }
        };
        
        for (size_t i = 0; i < is_spliced_1.size(); ++i) {
            if (is_spliced_1[i]) {
                for (size_t j = 0; j < multipath_alns_2.size(); ++j) {
                    attempt_to_pair(i, j);
                }
            }
        }
        for (size_t j = 0; j < multipath_alns_2.size(); ++j) {
            if (is_spliced_2[j]) {
                for (size_t i = 0; i < multipath_alns_1.size(); ++i) {
                    if (!is_spliced_1[i]) {
                        // we haven't tried this pairing in the opposite direction
                        attempt_to_pair(i, j);
                    }
                }
            }
        }
        
        if (found_consistent) {
            // count up how many rescued pairs use each original alignment so we can
            // know when it's safe to consume one
            unordered_map<size_t, int> left_count, right_count;
            for (const auto& cluster_pair : cluster_pairs_out) {
                ++left_count[cluster_pair.first.first];
                ++right_count[cluster_pair.first.second];
            }
            
            for (auto& cluster_pair : cluster_pairs_out) {
                // either copy or move the individual end mappings
                multipath_aln_pairs_out.emplace_back();
                if (--left_count[cluster_pair.first.first] == 0) {
                    multipath_aln_pairs_out.back().first = move(multipath_alns_1[cluster_pair.first.first]);
                }
                else {
                    multipath_aln_pairs_out.back().first = multipath_alns_1[cluster_pair.first.first];
                }
                if (--right_count[cluster_pair.first.second] == 0) {
                    multipath_aln_pairs_out.back().second = move(multipath_alns_2[cluster_pair.first.second]);
                }
                else {
                    multipath_aln_pairs_out.back().second = multipath_alns_2[cluster_pair.first.second];
                }
                
                // pair is at least as unique as either of its two ends
                pair_multiplicities_out.emplace_back(min(multiplicities_1[cluster_pair.first.first],
                                                         multiplicities_2[cluster_pair.first.second]));
                
                // fix cluster pair to point at clusters instead of mappings
                cluster_pair.first.first = cluster_idxs_1[cluster_pair.first.first];
                cluster_pair.first.second = cluster_idxs_2[cluster_pair.first.second];
            }
        }
        return found_consistent;
    }

    void MultipathMapper::agglomerate(size_t idx, multipath_alignment_t& agglomerating, const multipath_alignment_t& multipath_aln,
                                      vector<size_t>& agglomerated_group, unordered_set<pos_t>& agg_start_positions,
                                      unordered_set<pos_t>& agg_end_positions) const {
        // does it look like we've already agglomerated a mapping at this position
        vector<pos_t> start_positions, end_positions;
        for (auto j : multipath_aln.start()) {
            start_positions.emplace_back(initial_position(multipath_aln.subpath(j).path()));
        }
        for (size_t j = 0; j < multipath_aln.subpath_size(); ++j) {
            if (multipath_aln.subpath(j).next_size() == 0) {
                end_positions.emplace_back(final_position(multipath_aln.subpath(j).path()));
            }
        }
        
        bool do_agglomerate = true;
        for (size_t j = 0; j < start_positions.size() && do_agglomerate; ++j) {
            do_agglomerate = !agg_start_positions.count(start_positions[j]);
        }
        for (size_t j = 0; j < end_positions.size() && do_agglomerate; ++j) {
            do_agglomerate = !agg_end_positions.count(end_positions[j]);
        }
        
        if (do_agglomerate) {
            // we want to merge this one in with the primary mapping
            agglomerated_group.push_back(idx);
            // skip the merge step if it's the first mapping
            if (idx > 0) {
                append_multipath_alignment(agglomerating, multipath_aln);
            }
            // record that the aggregated mapping now has these start and end positions
            agg_start_positions.insert(start_positions.begin(), start_positions.end());
            agg_end_positions.insert(end_positions.begin(), end_positions.end());
        }
    }

    void MultipathMapper::agglomerate_alignments(vector<multipath_alignment_t>& multipath_alns_out,
                                                 vector<double>* multiplicities) {
        
        if (multipath_alns_out.empty()) {
            return;
        }
        
        // the likelihoods of each alignment, which we assume to be sorted
        vector<double> scores = mapping_likelihoods(multipath_alns_out);
        auto alnr = get_aligner(!multipath_alns_out.front().quality().empty());
        double min_score = scores.front() - alnr->mapping_quality_score_diff(max_mapping_quality);
        
        size_t i;
        vector<size_t> agglomerated_group;
        unordered_set<pos_t> agg_start_positions, agg_end_positions;
        for (i = 0; i < multipath_alns_out.size(); ++i) {
            // is the score good enough to agglomerate?
            if (scores[i] < min_score || likely_mismapping(multipath_alns_out[i])) {
                // none of the following will be either
                break;
            }
            
            // apply the agglomeration procedure
            agglomerate(i, multipath_alns_out.front(),
                        multipath_alns_out[i], agglomerated_group,
                        agg_start_positions, agg_end_positions);
        }
        
        if (i > 1) {
            
            // figure out the mapping quality for the whole aggregated alignment
            double raw_mapq = alnr->compute_group_mapping_quality(scores, agglomerated_group,
                                                                  multiplicities);
            int32_t mapq = min<int32_t>(max_mapping_quality, int32_t(mapq_scaling_factor * raw_mapq));
            multipath_alns_out.front().set_mapping_quality(mapq);
            
            multipath_alns_out.front().set_annotation("disconnected", true);
            
            // move the remaining alignments up in the return vector and resize the remnants away
            for (size_t j = i, k = 1; j < multipath_alns_out.size(); ++j, ++k) {
                multipath_alns_out[k] = move(multipath_alns_out[j]);
            }
            multipath_alns_out.resize(multipath_alns_out.size() - i + 1);
        }
    }

    void MultipathMapper::agglomerate_alignment_pairs(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                      vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                      vector<double>& multiplicities) {

        if (multipath_aln_pairs_out.empty()) {
            return;
        }
        
        // the likelihoods of each alignment, which we assume to be sorted
        vector<double> scores = pair_mapping_likelihoods(multipath_aln_pairs_out, cluster_pairs);
        auto alnr = get_aligner(!multipath_aln_pairs_out.front().first.quality().empty()
                                && !multipath_aln_pairs_out.front().second.quality().empty());
        double min_score = scores.front() - alnr->mapping_quality_score_diff(max_mapping_quality);
        
        size_t i;
        vector<size_t> agglomerated_group_1, agglomerated_group_2;
        unordered_set<pos_t> agg_start_positions_1, agg_end_positions_1, agg_start_positions_2, agg_end_positions_2;
        for (i = 0; i < multipath_aln_pairs_out.size(); ++i) {
            // is the score good enough to agglomerate?
            if (scores[i] < min_score ||
                (likely_mismapping(multipath_aln_pairs_out[i].first) && likely_mismapping(multipath_aln_pairs_out[i].second))) {
                // none of the following will be either
                break;
            }
            
            auto& multipath_aln_pair = multipath_aln_pairs_out[i];
            
            // repeat agglomeration procedure for both ends of the pair
            agglomerate(i, multipath_aln_pairs_out.front().first,
                        multipath_aln_pair.first, agglomerated_group_1,
                        agg_start_positions_1, agg_end_positions_1);
            agglomerate(i, multipath_aln_pairs_out.front().second,
                        multipath_aln_pair.second, agglomerated_group_2,
                        agg_start_positions_2, agg_end_positions_2);
        }
        
        if (i > 1) {
            
            // figure out the mapping quality for the whole aggregated alignment
            double raw_mapq_1 = alnr->compute_group_mapping_quality(scores, agglomerated_group_1,
                                                                    &multiplicities);
            double raw_mapq_2 = alnr->compute_group_mapping_quality(scores, agglomerated_group_2,
                                                                    &multiplicities);
            int32_t mapq_1 = min<int32_t>(max_mapping_quality, int32_t(mapq_scaling_factor * raw_mapq_1));
            int32_t mapq_2 = min<int32_t>(max_mapping_quality, int32_t(mapq_scaling_factor * raw_mapq_2));
            multipath_aln_pairs_out.front().first.set_mapping_quality(mapq_1);
            multipath_aln_pairs_out.front().second.set_mapping_quality(mapq_2);
            multipath_aln_pairs_out.front().first.set_annotation("disconnected", true);
            multipath_aln_pairs_out.front().second.set_annotation("disconnected", true);
            
            // move the remaining alignments up in the return vector and resize the remnants away
            for (size_t j = i, k = 1; j < multipath_aln_pairs_out.size(); ++j, ++k) {
                multipath_aln_pairs_out[k] = move(multipath_aln_pairs_out[j]);
            }
            multipath_aln_pairs_out.resize(multipath_aln_pairs_out.size() - i + 1);
        }
    }

    void MultipathMapper::purge_unmapped_alignments(vector<multipath_alignment_t>& multipath_alns_out) {
        for (size_t i = 0; i < multipath_alns_out.size(); ++i) {
            // TODO: could do this more efficiently by bisect search
            if (likely_mismapping(multipath_alns_out[i])) {
                // we can't distinguish this alignment from the longest MEM of a random sequence
                // so we don't report this mapping
                if (i == 0) {
#ifdef debug_multipath_mapper
                    cerr << "mapping is not distinguishable from a random sequence, reporting as unmapped" << endl;
#endif
                    // leave an unmapped placeholder
                    multipath_alns_out.resize(1);
                    clear_alignment(multipath_alns_out.front());
                }
                else {
                    // truncate the output from this point on
                    multipath_alns_out.resize(i);
                }
                break;
            }
        }
    }

    void MultipathMapper::purge_unmapped_alignments(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out, bool proper_paired) {
        
        // decide if the read is unmapped
        if (proper_paired) {
            for (size_t i = 0; i < multipath_aln_pairs_out.size(); ++i) {
                // if they're part of a proper pair, we count even pretty bad alignments
                // as mapped
                if (likely_mismapping(multipath_aln_pairs_out[i].first) &&
                    likely_mismapping(multipath_aln_pairs_out[i].second)) {
                    // this pair is actually unmapped
                    // TODO: maybe we should accept it though, in a large genome they are very
                    // unlikely to be paired by chance...
                    if (i == 0) {
                        // the read is completely unmapped, get rid of multimappings
                        multipath_aln_pairs_out.resize(1);
                        clear_alignment(multipath_aln_pairs_out.front().first);
                        clear_alignment(multipath_aln_pairs_out.front().second);
                    }
                    else {
                        // truncate the list here
                        multipath_aln_pairs_out.resize(i);
                    }
                    break;
                }
            }
        }
        else {
            // we have to do some complicated logic to look for unmapped reads in both primary
            // and secondary because of vg's silly interleaving logic...
            
            // find where unmapped reads begin in each list
            size_t i, j;
            for (i = 0; i < multipath_aln_pairs_out.size(); ++i) {
                if (multipath_aln_pairs_out[i].first.subpath().empty()
                    || likely_mismapping(multipath_aln_pairs_out[i].first)) {
                    break;
                }
            }
            for (j = 0; j < multipath_aln_pairs_out.size(); ++j) {
                if (multipath_aln_pairs_out[j].second.subpath().empty()
                    || likely_mismapping(multipath_aln_pairs_out[j].second)) {
                    break;
                }
            }
            
            if (i != multipath_aln_pairs_out.size() || j != multipath_aln_pairs_out.size()) {
                // we need to clear out unmapped reads in at least one of the read lists
                
                // truncate
                multipath_aln_pairs_out.resize(max<size_t>(1, max(i, j)));
                // add placeholders for unmapped reads if necessary
                for (; i < multipath_aln_pairs_out.size(); ++i) {
                    clear_alignment(multipath_aln_pairs_out[i].first);
                }
                for (; j < multipath_aln_pairs_out.size(); ++j) {
                    clear_alignment(multipath_aln_pairs_out[j].second);
                }
            }
        }
    }

    void MultipathMapper::simplify_complicated_multipath_alignment(multipath_alignment_t& multipath_aln) const {
        
        if (multipath_aln.subpath_size() > prune_subpaths_multiplier * multipath_aln.sequence().size()) {
            // this is a very complicated multipath alignment relative to the length of the sequence,
            // so we'll see if we can maybe get rid of parts of it for having low score
            
            auto aligner = get_aligner(!multipath_aln.quality().empty());
            int32_t max_diff = ceil(aligner->mapping_quality_score_diff(max_mapping_quality) / mapq_scaling_factor);
            
#ifdef debug_multipath_mapper
            cerr << "multipath alignment has " << multipath_aln.subpath_size() << " subpaths, which is large relative to sequence length of " << multipath_aln.sequence().size() << ", attempting to simplify by pruning low-scoring sections, diff = " << max_diff << endl;
#endif
            
            remove_low_scoring_sections(multipath_aln, max_diff);
            
            // TODO: it would be nice to merge non branching subpaths here, but we don't know what
            // the prohibited merges were...
#ifdef debug_multipath_mapper
            cerr << "after pruning, alignment has " << multipath_aln.subpath_size() << " subpaths" << endl;
#endif
        }
    }
    
    void MultipathMapper::split_multicomponent_alignments(vector<multipath_alignment_t>& multipath_alns_out,
                                                          const Alignment* alignment,
                                                          vector<clustergraph_t>* cluster_graphs,
                                                          vector<size_t>* cluster_idxs,
                                                          vector<double>* multiplicities) const {
        
        size_t num_original_alns = multipath_alns_out.size();
        vector<size_t> split_idxs;
        for (size_t i = 0; i < num_original_alns; i++) {
            
            vector<vector<int64_t>> comps = connected_components(multipath_alns_out[i]);
            
            if (comps.size() > 1) {
#ifdef debug_multipath_mapper
                cerr << "splitting multicomponent alignment " << debug_string(multipath_alns_out[i]) << endl;
#endif
                // split this multipath alignment into its connected components
                for (size_t j = 1; j < comps.size(); j++) {
                    split_idxs.push_back(multipath_alns_out.size());
                    multipath_alns_out.emplace_back();
                    extract_sub_multipath_alignment(multipath_alns_out[i], comps[j],
                                                    multipath_alns_out.back());
                    // also label the split alignment with its cluster of origin, if we're keeping track of that
                    if (cluster_idxs) {
                        cluster_idxs->emplace_back(cluster_idxs->at(i));
                    }
                    if (multiplicities) {
                        multiplicities->emplace_back(multiplicities->at(i));
                    }
                }
                // put the first component into the original location
                multipath_alignment_t last_component;
                extract_sub_multipath_alignment(multipath_alns_out[i], comps[0], last_component);
                multipath_alns_out[i] = move(last_component);
                split_idxs.push_back(i);
            }
        }
        
        if (alignment && cluster_graphs && cluster_idxs && do_spliced_alignment && !split_idxs.empty()) {
            // we only do this in spliced alignment because we want to clustering to
            // unclaim certain hits so they can be seen as spliced alignment candidates
            
            vector<const multipath_alignment_t*> split_mp_alns(split_idxs.size());
            vector<size_t*> cluster_assignments(split_idxs.size());
            for (size_t i = 0; i < split_idxs.size(); ++i) {
                auto& mp_aln = multipath_alns_out[split_idxs[i]];
                // TODO: we need to have these be ordered to find the MEMs, but this will be wastefully repeated later
                topologically_order_subpaths(mp_aln);
                split_mp_alns[i] = &mp_aln;
                cluster_assignments[i] = &(*cluster_idxs)[split_idxs[i]];
            }
            
            vector<size_t*> all_cluster_assignments(cluster_idxs->size());
            for (size_t i = 0; i < all_cluster_assignments.size(); ++i) {
                all_cluster_assignments[i] = &(*cluster_idxs)[i];
            }
            
            reassign_split_clusters(*alignment, *cluster_graphs, split_mp_alns, cluster_assignments,
                                    all_cluster_assignments);
        }
    }

    void MultipathMapper::reassign_split_clusters(const Alignment& alignment,
                                                  vector<clustergraph_t>& cluster_graphs,
                                                  const vector<const multipath_alignment_t*>& split_mp_alns,
                                                  const vector<size_t*>& cluster_assignments,
                                                  const vector<size_t*>& all_cluster_assignments) const {
#ifdef debug_multipath_mapper
        cerr << "reassigning split clusters for mp alns" << endl;
        for (auto mp_aln : split_mp_alns) {
            cerr << debug_string(*mp_aln) << endl;
        }
#endif
        
        // it's often possible to divvy up the hits into smaller clusters now
        
        // reorganize the split alignments by their original cluster
        unordered_map<size_t, vector<size_t>> original_cluster;
        for (size_t i = 0; i < cluster_assignments.size(); ++i) {
            original_cluster[*cluster_assignments[i]].push_back(i);
#ifdef debug_multipath_mapper
            cerr << "mp aln " << i << " associated with cluster " << *cluster_assignments[i] << endl;
#endif
        }
        
        bool any_new_clusters = false;
        for (const auto& record : original_cluster) {
            
#ifdef debug_multipath_mapper
            cerr << "reassigning clusters for mp alns originally from cluster " << record.first << ":" << endl;
            for (auto i : record.second) {
                cerr << "\t" << i << endl;
            }
#endif
            
            bool all_hits_found = true;
            vector<vector<size_t>> hits_found_in_aln(record.second.size());
            
            // note: we use the ugly "get" function every time because we're also going to be modifying
            // the cluster graphs vector and we don't want there to be trouble with a local reference
            
            for (size_t i = 0; i < get<1>(cluster_graphs[record.first]).first.size() && all_hits_found; ++i) {
                all_hits_found = false;
                auto& hit = get<1>(cluster_graphs[record.first]).first[i];
                for (size_t j = 0; j < record.second.size(); ++j) {
#ifdef debug_multipath_mapper
                    cerr << "checking for " << i << "-th hit in this cluster the in " << j << "-th mp aln" << endl;
#endif
                    
                    // check if the i-th MEM is contained in the j-th split up component
                    if (contains_match(*split_mp_alns[record.second[j]], hit.second,
                                       hit.first->begin - alignment.sequence().begin(), hit.first->length())) {
                        
#ifdef debug_multipath_mapper
                        cerr << "hit found" << endl;
#endif
                        
                        all_hits_found = true;
                        hits_found_in_aln[j].push_back(i);
                    }
                }
            }
            
            if (!all_hits_found) {
                // our partition is incomplete, which makes the process of divvying up the
                // hits too complicated (some will get lost), so just skip it
                // TODO: is it really a problem if we lose some hits?
#ifdef debug_multipath_mapper
                cerr << "not all hits are still found, skipping reassignment" << endl;
#endif
                continue;
            }
            
            map<vector<size_t>, size_t> new_cluster;
            for (size_t j = 0; j < record.second.size(); ++j) {
#ifdef debug_multipath_mapper
                cerr << "reassigning cluster of " << j << "-th mp aln, " << record.second[j] << endl;
#endif
                if (hits_found_in_aln[j].size() == get<1>(cluster_graphs[record.first]).first.size()) {
                    // this alignment still contains the whole cluster, so we don't need
                    // to point it at anything else
#ifdef debug_multipath_mapper
                    cerr << "still in original cluster " << record.first << endl;
#endif
                    continue;
                }
                
                auto it = new_cluster.find(hits_found_in_aln[j]);
                if (it == new_cluster.end()) {
#ifdef debug_multipath_mapper
                    cerr << "making a new cluster at index " << cluster_graphs.size() << " with hits:" << endl;
                    for (auto hit_idx : hits_found_in_aln[j]) {
                        cerr << "\t" << hit_idx << endl;
                    }
#endif
                    
                    // this is the first time we're encountering this combination of hits, so we need to make
                    // a corresponding cluster
                    it = new_cluster.insert(make_pair(move(hits_found_in_aln[j]), cluster_graphs.size())).first;
                    
                    // and now make the cluster graph itself
                    cluster_graphs.emplace_back();
                    auto& cluster_graph = cluster_graphs.back();
                    get<0>(cluster_graph) = unique_ptr<bdsg::HashGraph>(new bdsg::HashGraph());
                    handlealgs::copy_handle_graph(get<0>(cluster_graphs[record.first]).get(),
                                                  get<0>(cluster_graph).get());
                    
                    for (auto i : it->first) {
                        get<1>(cluster_graph).first.emplace_back(get<1>(cluster_graphs[record.first]).first[i]);
                    }
                    get<1>(cluster_graph).second = get<1>(cluster_graphs[record.first]).second;
                    
#ifdef debug_multipath_mapper
                    cerr << "\tmultiplicity " << get<1>(cluster_graph).second << endl;
#endif
                    
                    set_read_coverage(cluster_graph);
                    
#ifdef debug_multipath_mapper
                    cerr << "assign total coverage " << get<2>(cluster_graphs.back()) << endl;
#endif
                    any_new_clusters = true;
                }
#ifdef debug_multipath_mapper
                cerr << "assigned to cluster " << it->second << endl;
#endif
                // reassign the cluster idx to the newer cluster
                *cluster_assignments[record.second[j]] = it->second;
            }
        }
        
        if (any_new_clusters) {
            // we reassigned some clusters, now let's check if we can clear any of the old clusters
            
            // clear the original clusters out of the map if at least one alignment is
            // still using them
            for (auto cluster_ptr : all_cluster_assignments) {
                auto it = original_cluster.find(*cluster_ptr);
                if (it != original_cluster.end()) {
                    original_cluster.erase(it);
                }
            }
            
            // any remaining clusters aren't being used by any alignment and they can
            // have their hits released
            for (const auto& record : original_cluster) {
#ifdef debug_multipath_mapper
                cerr << "cluster " << record.first << " has no more assignees, unclaiming all of its hits" << endl;
#endif
                get<1>(cluster_graphs[record.first]).first.clear();
                get<2>(cluster_graphs[record.first]) = 0;
            }
            
            
            // reorder the clusters based on the new coverage
            vector<size_t> order(cluster_graphs.size(), 0);
            for (size_t i = 1; i < cluster_graphs.size(); ++i) {
                order[i] = i;
            }
            
            stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
                return get<2>(cluster_graphs[i]) > get<2>(cluster_graphs[j]);
            });
            
            vector<size_t> index(order.size());
            for (size_t i = 0; i < order.size(); ++i) {
                index[order[i]] = i;
            }
            
#ifdef debug_multipath_mapper
            cerr << "reordering clusters based on coverage" << endl;
            for (size_t i = 0; i < index.size(); ++i) {
                cerr << "\t" << i << " -> " << index[i] << endl;
            }
#endif
            
            for (size_t* cluster_assignment : all_cluster_assignments) {
                if (*cluster_assignment != RESCUED) {
                    *cluster_assignment = index[*cluster_assignment];
                }
            }
            
            for (size_t i = 0; i < index.size(); ++i) {
                while (index[i] != i) {
                    std::swap(cluster_graphs[index[i]], cluster_graphs[i]);
                    std::swap(index[index[i]], index[i]);
                    
                }
            }
            
#ifdef debug_multipath_mapper
            cerr << "new cluster ordering" << endl;
            for (int i = 0; i < cluster_graphs.size(); i++) {
                cerr << "cluster " << i << ", coverage " << get<2>(cluster_graphs[i]) << ", multiplicity " << get<1>(cluster_graphs[i]).second << endl;
                for (pair<const MaximalExactMatch*, pos_t>  hit : get<1>(cluster_graphs[i]).first) {
                    cerr << "\t" << hit.second << " " <<  hit.first->sequence() << endl;
                }
            }
#endif
        }
    }

    void MultipathMapper::merge_rescued_mappings(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                 vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                 vector<double>& pair_multiplicities,
                                                 vector<pair<multipath_alignment_t, multipath_alignment_t>>& rescued_multipath_aln_pairs,
                                                 vector<pair<pair<size_t, size_t>, int64_t>>& rescued_cluster_pairs,
                                                 vector<double>& rescued_multiplicities) const {
                
        size_t num_unrescued_pairs = multipath_aln_pairs_out.size();
        
        for (size_t j = 0; j < rescued_multipath_aln_pairs.size(); j++) {
            
            // make sure this pair isn't a duplicate with any of the original pairs
            bool duplicate = false;
            for (size_t i = 0; i < num_unrescued_pairs; i++) {
#ifdef debug_multipath_mapper
                cerr << "checking if rescue pair " << j << " is duplicate of original pair " << i << endl;
#endif
                if (share_terminal_positions(multipath_aln_pairs_out[i].first, rescued_multipath_aln_pairs[j].first)) {
                    if (share_terminal_positions(multipath_aln_pairs_out[i].second, rescued_multipath_aln_pairs[j].second)) {
#ifdef debug_multipath_mapper
                        cerr << "found a duplicate" << endl;
#endif
                        duplicate = true;
                        break;
                    }
                }
            }
            
            if (!duplicate) {
#ifdef debug_multipath_mapper
                cerr << "no duplicate, adding to return vector if distance is finite and positive" << endl;
#endif
                multipath_aln_pairs_out.emplace_back(move(rescued_multipath_aln_pairs[j]));
                cluster_pairs.emplace_back(rescued_cluster_pairs[j]);
                pair_multiplicities.emplace_back(rescued_multiplicities[j]);
            }
        }
        
        sort_and_compute_mapping_quality(multipath_aln_pairs_out, cluster_pairs, nullptr, &pair_multiplicities);
    }

    double MultipathMapper::estimate_missed_rescue_multiplicity(size_t which_pair,
                                                                const vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                                const vector<clustergraph_t>& cluster_graphs1,
                                                                const vector<clustergraph_t>& cluster_graphs2,
                                                                bool from_secondary_rescue) const {
#ifdef debug_multipath_mapper
        cerr << "checking whether we should enter rescue multiplicity routine" << endl;
#endif
        
        
        double multiplicity = 1.0;
        
        // did we use an out-of-bounds cluster index to flag either end as coming from a rescue?
        bool opt_aln_1_is_rescued = cluster_pairs[which_pair].first.first == RESCUED;
        bool opt_aln_2_is_rescued = cluster_pairs[which_pair].first.second == RESCUED;
        
        // was the optimal cluster pair obtained by rescue?
        if (opt_aln_1_is_rescued || opt_aln_2_is_rescued) {
            // let's figure out if we should reduce its mapping quality to reflect the fact that we may not have selected the
            // correct cluster as a rescue candidate
            
#ifdef debug_multipath_mapper
            cerr << "the optimal alignment is a rescue, checking if we need to cap the mapping quality" << endl;
#endif
            
            const vector<clustergraph_t>& anchor_clusters = opt_aln_1_is_rescued ? cluster_graphs2 : cluster_graphs1;
            size_t anchor_idx = opt_aln_1_is_rescued ? cluster_pairs[which_pair].first.second : cluster_pairs[which_pair].first.first;
            
            // find the range of clusters that could plausibly be about as good as the one that rescue succeeded from
            size_t plausible_clusters_end_idx = anchor_idx;
            for (; plausible_clusters_end_idx < anchor_clusters.size(); plausible_clusters_end_idx++) {
                if (get<2>(anchor_clusters[plausible_clusters_end_idx]) < get<2>(anchor_clusters[anchor_idx]) - plausible_rescue_cluster_coverage_diff) {
                    break;
                }
            }
            
            // TODO: it's a bit ugly/fragile to have the secondary rescue logic recapitulated here
            
            // figure out which index corresponds to the end of the range we would have rescued
            size_t max_rescues_attempted_idx;
            if (from_secondary_rescue) {
                // find the indexes that were added by pair clustering
                unordered_set<size_t> paired_idxs;
                for (auto& cluster_pair : cluster_pairs) {
                    paired_idxs.insert(opt_aln_1_is_rescued ? cluster_pair.first.second : cluster_pair.first.first);
                }
                
                // the "budget" of rescues we could have performed
                size_t rescues_left = secondary_rescue_attempts;
                size_t i = 0;
                for (; i < anchor_clusters.size(); i++) {
                    // did we skip thi index because it was already in a pair?
                    if (!paired_idxs.count(i)) {
                        if (rescues_left) {
                            // we would have tried a secondary rescue here
                            rescues_left--;
                        }
                        else {
                            // we would have run out of secondary rescues here
                            break;
                        }
                    }
                }
                // this is the first index we didn't rescue
                max_rescues_attempted_idx = i;
            }
            else {
                // simpler without secondary rescue, we would have just rescued up to the maximum allowable
                max_rescues_attempted_idx = max_rescue_attempts;
            }
            
#ifdef debug_multipath_mapper
            cerr << "performed up to " << max_rescues_attempted_idx << " out of " << plausible_clusters_end_idx << " plausible rescues" << endl;
#endif
            
            multiplicity = max(1.0, (double)(plausible_clusters_end_idx) / (double)(max_rescues_attempted_idx));
        }
        
        return multiplicity;
    }

    double MultipathMapper::cluster_multiplicity(const memcluster_t& cluster) const {
        if (cluster.first.empty()) {
            return cluster.second;
        }
        double max_fraction_sampled = 0.0;
        for (const pair<const MaximalExactMatch*, pos_t>& hit : cluster.first) {
            const MaximalExactMatch& mem = *hit.first;
            max_fraction_sampled = max(max_fraction_sampled, (double)(mem.queried_count) / (double)(mem.match_count));
        }
        return cluster.second / max_fraction_sampled;
    }

    double MultipathMapper::pair_cluster_multiplicity(const memcluster_t& cluster_1, const memcluster_t& cluster_2) const {
        return min(cluster_multiplicity(cluster_1), cluster_multiplicity(cluster_2));
    }

    MultipathMapper::match_fanouts_t MultipathMapper::record_fanouts(const vector<MaximalExactMatch>& mems,
                                                                     vector<deque<pair<string::const_iterator, char>>>& fanouts) const {
        
        match_fanouts_t match_fanouts;
        if (!fanouts.empty()) {
            assert(fanouts.size() == mems.size());
            for (size_t i = 0; i < mems.size(); ++i) {
                if (!fanouts[i].empty()) {
                    match_fanouts[&mems[i]] = move(fanouts[i]);
                }
            }
        }
        return match_fanouts;
    }
    
    void MultipathMapper::split_multicomponent_alignments(const Alignment& alignment1, const Alignment& alignment2,
                                                          vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                          vector<clustergraph_t>& cluster_graphs1,
                                                          vector<clustergraph_t>& cluster_graphs2,
                                                          vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                          vector<double>& pair_multiplicities) const {
        
        size_t original_num_pairs = multipath_aln_pairs_out.size();
        vector<size_t> split_idxs_1, split_idxs_2;
        for (size_t i = 0; i < original_num_pairs; i++) {
            vector<vector<int64_t>> connected_components_1 = connected_components(multipath_aln_pairs_out[i].first);
            vector<vector<int64_t>> connected_components_2 = connected_components(multipath_aln_pairs_out[i].second);
            
#ifdef debug_multipath_mapper
            cerr << "finding connected components for mapping:" << endl;
#endif
#ifdef debug_multipath_mapper_alignment
            view_multipath_alignment_as_dot(cerr, multipath_aln_pairs_out[i].first);
            view_multipath_alignment_as_dot(cerr, multipath_aln_pairs_out[i].second);
#endif
#ifdef debug_multipath_mapper
            cerr << "read 1 connected components:" << endl;
            for (vector<int64_t>& comp : connected_components_1) {
                cerr << "\t";
                for (int64_t j : comp) {
                    cerr << j << " ";
                }
                cerr << endl;
            }
            cerr << "read 2 connected components:" << endl;
            for (vector<int64_t>& comp : connected_components_2) {
                cerr << "\t";
                for (int64_t j : comp) {
                    cerr << j << " ";
                }
                cerr << endl;
            }
#endif
            // we will put pairs of split up components in here
            vector<pair<multipath_alignment_t, multipath_alignment_t>> split_multipath_alns;
            
            if (connected_components_1.size() > 1 && connected_components_2.size() > 1) {
#ifdef debug_multipath_mapper
                cerr << "splitting both multicomponent alignments" << endl;
#endif
                // need to split both ends
                for (size_t j = 0; j < connected_components_1.size(); j++) {
                    for (size_t k = 0; k < connected_components_2.size(); k++) {
                        split_multipath_alns.emplace_back();
                        extract_sub_multipath_alignment(multipath_aln_pairs_out[i].first, connected_components_1[j],
                                                        split_multipath_alns.back().first);
                        extract_sub_multipath_alignment(multipath_aln_pairs_out[i].second, connected_components_2[k],
                                                        split_multipath_alns.back().second);
                    }
                }
            }
            else if (connected_components_1.size() > 1) {
#ifdef debug_multipath_mapper
                cerr << "splitting read 1 multicomponent alignments" << endl;
#endif
                // only need to split first end
                for (size_t j = 0; j < connected_components_1.size(); j++) {
                    split_multipath_alns.emplace_back(multipath_alignment_t(), multipath_aln_pairs_out[i].second);
                    extract_sub_multipath_alignment(multipath_aln_pairs_out[i].first, connected_components_1[j],
                                                    split_multipath_alns.back().first);
                }
            }
            else if (connected_components_2.size() > 1) {
#ifdef debug_multipath_mapper
                cerr << "splitting read 2 multicomponent alignments" << endl;
#endif
                // only need to split second end
                for (size_t j = 0; j < connected_components_2.size(); j++) {
                    split_multipath_alns.emplace_back(multipath_aln_pairs_out[i].first, multipath_alignment_t());
                    extract_sub_multipath_alignment(multipath_aln_pairs_out[i].second, connected_components_2[j],
                                                    split_multipath_alns.back().second);
                }
            }
            
            // are there split up pairs to add to the output vector?
            if (!split_multipath_alns.empty()) {
                
                bool replaced_original = false;
                for (pair<multipath_alignment_t, multipath_alignment_t>& split_multipath_aln_pair : split_multipath_alns) {
                    // we also need to measure the disance for scoring
                    int64_t dist = distance_between(split_multipath_aln_pair.first, split_multipath_aln_pair.second, true);
                    
                    // if we can't measure a distance, then don't add the pair
                    if (dist != numeric_limits<int64_t>::max()) {
                        
#ifdef debug_multipath_mapper
                        cerr << "adding component pair at distance " << dist << ":" << endl;
                        cerr  << debug_string(split_multipath_aln_pair.first) << endl;
                        cerr  << debug_string(split_multipath_aln_pair.second) << endl;
#endif
                        
                        if (!replaced_original) {
                            // put the first one back into the original position in the output vector
                            multipath_aln_pairs_out[i] = move(split_multipath_aln_pair);
                            cluster_pairs[i].second = dist;
                            replaced_original = true;
                            if (connected_components_1.size() > 1) {
                                split_idxs_1.push_back(i);
                            }
                            if (connected_components_2.size() > 1) {
                                split_idxs_2.push_back(i);
                            }
                        }
                        else {
                            // append the rest of them to the end
                            if (connected_components_1.size() > 1) {
                                split_idxs_1.push_back(multipath_aln_pairs_out.size());
                            }
                            if (connected_components_2.size() > 1) {
                                split_idxs_2.push_back(multipath_aln_pairs_out.size());
                            }
                            multipath_aln_pairs_out.emplace_back(move(split_multipath_aln_pair));
                            cluster_pairs.emplace_back(cluster_pairs[i].first, dist);
                            pair_multiplicities.emplace_back(pair_multiplicities[i]);
                        }
                    }
                }
            }
        }
        
        if (do_spliced_alignment) {
            // we only do this in spliced alignment because we want to clustering to
            // unclaim certain hits so they can be seen as spliced alignment candidates
            
            if (!split_idxs_1.empty()) {
                
                vector<const multipath_alignment_t*> split_mp_alns_1(split_idxs_1.size());
                vector<size_t*> cluster_assignments_1(split_idxs_1.size());
                for (size_t i = 0; i < split_idxs_1.size(); ++i) {
                    auto& mp_aln = multipath_aln_pairs_out[split_idxs_1[i]].first;
                    // TODO: we need to have these be ordered to find the MEMs, but this will be wastefully repeated later
                    topologically_order_subpaths(mp_aln);
                    split_mp_alns_1[i] = &mp_aln;
                    cluster_assignments_1[i] = &cluster_pairs[split_idxs_1[i]].first.first;
                }
                
                vector<size_t*> all_cluster_assignments_1(cluster_pairs.size());
                for (size_t i = 0; i < all_cluster_assignments_1.size(); ++i) {
                    all_cluster_assignments_1[i] = &cluster_pairs[i].first.first;
                }
                
                reassign_split_clusters(alignment1, cluster_graphs1, split_mp_alns_1, cluster_assignments_1,
                                        all_cluster_assignments_1);
            }
            
            if (!split_idxs_2.empty()) {
                vector<const multipath_alignment_t*> split_mp_alns_2(split_idxs_2.size());
                vector<size_t*> cluster_assignments_2(split_idxs_2.size());
                for (size_t i = 0; i < split_idxs_2.size(); ++i) {
                    auto& mp_aln = multipath_aln_pairs_out[split_idxs_2[i]].second;
                    // TODO: we need to have these be ordered to find the MEMs, but this will be wastefully repeated later
                    topologically_order_subpaths(mp_aln);
                    split_mp_alns_2[i] = &mp_aln;
                    cluster_assignments_2[i] = &cluster_pairs[split_idxs_2[i]].first.second;
                }
                
                vector<size_t*> all_cluster_assignments_2(cluster_pairs.size());
                for (size_t i = 0; i < all_cluster_assignments_2.size(); ++i) {
                    all_cluster_assignments_2[i] = &cluster_pairs[i].first.second;
                }
                
                reassign_split_clusters(alignment2, cluster_graphs2, split_mp_alns_2, cluster_assignments_2,
                                        all_cluster_assignments_2);
            }
        }
    }
    
    void MultipathMapper::align_to_cluster_graph_pairs(const Alignment& alignment1, const Alignment& alignment2,
                                                       vector<clustergraph_t>& cluster_graphs1,
                                                       vector<clustergraph_t>& cluster_graphs2,
                                                       vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                       vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                       vector<double>& pair_multiplicities,
                                                       vector<pair<size_t, size_t>>& duplicate_pairs_out,
                                                       const match_fanouts_t* fanouts1, const match_fanouts_t* fanouts2) {
        
        assert(multipath_aln_pairs_out.empty());
        
        auto aligner = get_aligner(!alignment1.quality().empty() && !alignment2.quality().empty());
        auto get_pair_approx_likelihood = [&](const pair<pair<size_t, size_t>, int64_t>& cluster_pair) {
            return ((get<2>(cluster_graphs1[cluster_pair.first.first])
                     + get<2>(cluster_graphs2[cluster_pair.first.second])) * aligner->match
                    + fragment_length_log_likelihood(cluster_pair.second) / aligner->log_base);
        };
        
        // sort the pairs descending by approximate likelihood
        stable_sort(cluster_pairs.begin(), cluster_pairs.end(),
                    [&](const pair<pair<size_t, size_t>, int64_t>& a, const pair<pair<size_t, size_t>, int64_t>& b) {
            // compute approximate likelihood in similar way to how the mapping quality routine will
            double likelihood_1 = get_pair_approx_likelihood(a);
            double likelihood_2 = get_pair_approx_likelihood(b);
            size_t hash_1 = wang_hash<pair<pair<size_t, size_t>, int64_t>>()(a);
            size_t hash_2 = wang_hash<pair<pair<size_t, size_t>, int64_t>>()(b);
            return (likelihood_1 > likelihood_2 || (likelihood_1 == likelihood_2 && hash_1 < hash_2));
        });
        
#ifdef debug_multipath_mapper
        cerr << "sorting cluster pairs by approximate likelihood:" << endl;
        for (size_t i = 0; i < cluster_pairs.size(); i++) {
            cerr << i << "-th cluster: " << cluster_pairs[i].first.first << " " << cluster_pairs[i].first.second << ", likelihood " << get_pair_approx_likelihood(cluster_pairs[i]) << endl;
        }
        
        cerr << "aligning to cluster pairs..." << endl;
#endif
        
        // we may need to compute an extra mapping above the one we'll report if we're computing mapping quality
        size_t num_mappings_to_compute = max(num_mapping_attempts, (size_t) 2);
        
        // TODO: some cluster pairs will produce redundant subgraph pairs.
        // We'll end up with redundant pairs being output.
        
        // the same index may occur in multiple pairs, if so we can copy it rather than needing to recompute it
        // we keep track of where the original occurrence was here
        unordered_map<size_t, size_t> previous_multipath_alns_1, previous_multipath_alns_2;
        
        // align to each cluster pair
        multipath_aln_pairs_out.reserve(min(num_mappings_to_compute, cluster_pairs.size()));
        size_t num_mappings = 0;
        for (size_t i = 0; i < cluster_pairs.size(); ++i) {
            // For each cluster pair
            const pair<pair<size_t, size_t>, int64_t>& cluster_pair = cluster_pairs[i];
            
            // TODO: using a multiplier here instead of a difference is pretty ugly, really. it also has
            // weird effects, like not producing any alignments if the log likelihood is negative (which
            // shouldn't matter). but in practice that only happens on very small clusters with bad fragment
            // lengths.
            
            // if we have a cluster graph pair with small enough MEM coverage
            // compared to the best one or we've made the maximum number of
            // alignments we stop producing alternate alignments
            if (get_pair_approx_likelihood(cluster_pair) < mem_coverage_min_ratio * get_pair_approx_likelihood(cluster_pairs.front())
                || num_mappings >= num_mappings_to_compute) {
                
                // remove the rest of the cluster pairs to establish the invariant that there are the
                // same number of cluster pairs as alternate mappings
                cluster_pairs.resize(i);
                
                break;
            }
            
#ifdef debug_multipath_mapper
            cerr << "doing pair " << cluster_pair.first.first << " " << cluster_pair.first.second << endl;
#endif
            // create multipath alignments to fill
            multipath_aln_pairs_out.emplace_back();
            pair_multiplicities.push_back(pair_cluster_multiplicity(get<1>(cluster_graphs1[cluster_pair.first.first]),
                                                                    get<1>(cluster_graphs2[cluster_pair.first.second])));
                        
            auto prev_1 = previous_multipath_alns_1.find(cluster_pair.first.first);
            if (prev_1 == previous_multipath_alns_1.end()) {
                // we haven't done this alignment yet, so we have to complete it for the first time
                
#ifdef debug_multipath_mapper
                cerr << "performing alignment of read 1 to subgraph" << endl;
#endif
                
                multipath_align(alignment1, cluster_graphs1[cluster_pair.first.first],
                                multipath_aln_pairs_out.back().first,
                                fanouts1);
                
                // keep track of the fact that we have completed this multipath alignment
                previous_multipath_alns_1[cluster_pair.first.first] = i;
            }
            else {
#ifdef debug_multipath_mapper
                cerr << "copying alignment from read1 at index " << prev_1->second << endl;
#endif
                // we've already completed this multipath alignment, so we can copy it
                multipath_aln_pairs_out.back().first = multipath_aln_pairs_out[prev_1->second].first;
            }
            
            // repeat this routine for the second read end
            // TODO: repetitive code
            auto prev_2 = previous_multipath_alns_2.find(cluster_pair.first.second);
            if (prev_2 == previous_multipath_alns_2.end()) {
                // we haven't done this alignment yet, so we have to complete it for the first time
                
#ifdef debug_multipath_mapper
                cerr << "performing alignment of read 2 to subgraph" << endl;
#endif
                
                multipath_align(alignment2, cluster_graphs2[cluster_pair.first.second],
                                multipath_aln_pairs_out.back().second, fanouts2);
                
                // keep track of the fact that we have completed this multipath alignment
                previous_multipath_alns_2[cluster_pair.first.second] = i;
            }
            else {
#ifdef debug_multipath_mapper
                cerr << "copying alignment from read2 at index " << prev_2->second << endl;
#endif
                // we've already completed this multipath alignment, so we can copy it
                multipath_aln_pairs_out.back().second = multipath_aln_pairs_out[prev_2->second].second;
            }
            
            num_mappings++;
        }
        
        if (!multipath_aln_pairs_out.empty()) {
            
            double likelihood_diff = aligner->mapping_quality_score_diff(truncation_multiplicity_mq_limit);
            double tail_likelihood = get_pair_approx_likelihood(cluster_pairs[multipath_aln_pairs_out.size() - 1]);
            
            // find clusters whose likelihoods are approximately the same as the low end of the clusters we aligned
            int64_t max_tail_idx = multipath_aln_pairs_out.size();
            while (max_tail_idx < cluster_pairs.size()
                   && get_pair_approx_likelihood(cluster_pairs[max_tail_idx]) >= tail_likelihood - likelihood_diff) {
                ++max_tail_idx;
            }
            
            if (max_tail_idx > multipath_aln_pairs_out.size()) {
                // there are some (nearly) identical cluster pairs that we ignored, so we'll account for them in the multiplicity
                
                // find the pairs that are approximately the same as the last
                int64_t min_tail_idx = multipath_aln_pairs_out.size() - 1;
                while (min_tail_idx > 0 &&
                       get_pair_approx_likelihood(cluster_pairs[min_tail_idx]) <= tail_likelihood + likelihood_diff) {
                    --min_tail_idx;
                }
                
                // multiply their multiplicity by the inverse of the fraction aligned
                double trunc_multiplicity = double(max_tail_idx - min_tail_idx) / double(multipath_aln_pairs_out.size() - min_tail_idx);
                for (size_t i = min_tail_idx; i < multipath_aln_pairs_out.size(); ++i) {
                    pair_multiplicities[i] *= trunc_multiplicity;
                }
            }
        }
        
        if (!suppress_multicomponent_splitting) {
            // split up any multi-component multipath alignments
            split_multicomponent_alignments(alignment1, alignment2,
                                            multipath_aln_pairs_out, cluster_graphs1, cluster_graphs2,
                                            cluster_pairs, pair_multiplicities);
            
            // it's possible to get supraquadratic growth in the number of alignments from doing this, so
            // we add a bit of restraint
            size_t max_split_pairs = num_mappings_to_compute * num_mappings_to_compute;
            if (multipath_aln_pairs_out.size() > max_split_pairs) {
#ifdef debug_multipath_mapper
                cerr << "too many pairs after splitting multicomponent alignments, truncating from " << multipath_aln_pairs_out.size() << " to " << max_split_pairs << " pairs" << endl;
#endif
                
                // TODO: repetitive with the previous truncation routine
                // TODO: these likelihood were just computed in the sort routine, it would be nice to do
                // be able to re-use the results
                
                // figure out a truncation multiplier for this pruning step
                auto scores = pair_mapping_likelihoods(multipath_aln_pairs_out, cluster_pairs);
                double score_diff = aligner->mapping_quality_score_diff(truncation_multiplicity_mq_limit);
                size_t last_idx = max_split_pairs - 1;
                size_t min_tail_idx = last_idx;
                while (min_tail_idx > 0 && scores[min_tail_idx - 1] <= scores[last_idx] + score_diff) {
                    --min_tail_idx;
                }
                size_t max_tail_idx = max_split_pairs;
                while (max_tail_idx < scores.size() && scores[max_tail_idx] >= scores[last_idx] - score_diff) {
                    ++max_tail_idx;
                }
                double trunc_multiplicity = double(max_tail_idx - min_tail_idx) / double(max_split_pairs - min_tail_idx);
                
                // discard the low-likelihood pairs
                multipath_aln_pairs_out.resize(max_split_pairs);
                cluster_pairs.resize(max_split_pairs);
                pair_multiplicities.resize(max_split_pairs);
                
                // increase the multiplicity for pairs near the cutoff point
                for (size_t i = min_tail_idx; i < pair_multiplicities.size(); ++i) {
                    pair_multiplicities[i] *= trunc_multiplicity;
                }
            }
        }
        
        // downstream algorithms assume multipath alignments are topologically sorted (including the scoring
        // algorithm in the next step)
        for (auto& multipath_aln_pair : multipath_aln_pairs_out) {
            topologically_order_subpaths(multipath_aln_pair.first);
            topologically_order_subpaths(multipath_aln_pair.second);
            simplify_complicated_multipath_alignment(multipath_aln_pair.first);
            simplify_complicated_multipath_alignment(multipath_aln_pair.second);
        }
        
        // put pairs in score sorted order and compute mapping quality of best pair using the score
        sort_and_compute_mapping_quality(multipath_aln_pairs_out, cluster_pairs,
                                         &duplicate_pairs_out, &pair_multiplicities);
        
#ifdef debug_validate_multipath_alignments
        for (pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : multipath_aln_pairs_out) {
#ifdef debug_multipath_mapper
            cerr << "validating multipath alignments:" << endl;
            cerr << debug_string(multipath_aln_pair.first) << endl;
            cerr << debug_string(multipath_aln_pair.second) << endl;
#endif
            if (!validate_multipath_alignment(multipath_aln_pair.first, *xindex)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.first.sequence() << " failed to validate" << endl;
            }
            if (!validate_multipath_alignment(multipath_aln_pair.second, *xindex)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.second.sequence() << " failed to validate" << endl;
            }
        }
#endif
        
    }

    pair<unique_ptr<bdsg::HashGraph>, bool> MultipathMapper::extract_maximal_graph(const Alignment& alignment,
                                                                                   const memcluster_t& mem_cluster) const {
        
        // Figure out the aligner to use
        auto aligner = get_aligner(!alignment.quality().empty());
        // get the seed hits
        const auto& cluster = mem_cluster.first;
        
        vector<pos_t> positions;
        vector<size_t> forward_max_dist;
        vector<size_t> backward_max_dist;
        
        positions.reserve(cluster.size());
        forward_max_dist.reserve(cluster.size());
        backward_max_dist.reserve(cluster.size());
        
        for (auto& mem_hit : cluster) {
            // get the start position of the MEM
            positions.push_back(mem_hit.second);
            // search far enough away to get any hit detectable without soft clipping
            forward_max_dist.push_back(min(aligner->longest_detectable_gap(alignment, mem_hit.first->end), max_alignment_gap)
                                       + (alignment.sequence().end() - mem_hit.first->begin));
            backward_max_dist.push_back(min(aligner->longest_detectable_gap(alignment, mem_hit.first->begin), max_alignment_gap)
                                        + (mem_hit.first->begin - alignment.sequence().begin()));
        }
        
        // TODO: a progressive expansion of the subgraph if the MEM hit is already contained in
        // a cluster graph somewhere?
        
        // extract the subgraph within the search distance
        
        unique_ptr<bdsg::HashGraph> cluster_graph(new bdsg::HashGraph());
        
        algorithms::extract_containing_graph(xindex, cluster_graph.get(), positions, forward_max_dist, backward_max_dist,
                                             num_alt_alns > 1 ? reversing_walk_length : 0);
        
        return move(make_pair(move(cluster_graph), cluster.size() == 1));
    }

    // TODO: entirely duplicative with MultipathAlignmentGraph...
    const size_t MultipathMapper::gap_memo_max_size = 1000;
    thread_local unordered_map<double, vector<int64_t>> MultipathMapper::pessimistic_gap_memo;
    int64_t MultipathMapper::pessimistic_gap(int64_t length, double multiplier) const {
        int64_t gap_length;
        if (length >= gap_memo_max_size) {
            gap_length = multiplier * sqrt(length);
        }
        else {
            vector<int64_t>& memo = pessimistic_gap_memo[multiplier];
            while (memo.size() <= length) {
                memo.emplace_back(multiplier * sqrt(memo.size()));
            }
            gap_length = memo[length];
        }
        return gap_length;
    }

    pair<unique_ptr<bdsg::HashGraph>, bool> MultipathMapper::extract_restrained_graph(const Alignment& alignment,
                                                                                      const memcluster_t& mem_cluster) const {
        
        // Figure out the aligner to use
        auto aligner = get_aligner(!alignment.quality().empty());
        // get the seed hits
        const auto& cluster = mem_cluster.first;
        
        // the MEMs are size sorted, we want to know the read order so we can
        // use the inter-MEM distance to figure out how much to extract
        vector<size_t> order(cluster.size(), 0);
        for (size_t i = 1; i < order.size(); ++i) {
            order[i] = i;
        }
        stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) {
            return cluster[i].first->begin < cluster[j].first->begin;
        });
        
        // and we'll also want to
        vector<size_t> index(order.size());
        for (size_t i = 0; i < index.size(); ++i) {
            index[order[i]] = i;
        }
        
        vector<pos_t> positions(cluster.size());
        
        // determine an initial restrained set of distances to extract from
        vector<size_t> forward_dist(cluster.size()), backward_dist(cluster.size());
        for (size_t i = 0; i < cluster.size(); ++i) {
            size_t idx = index[i];
            if (idx == 0) {
                // this is the left tail
                if (use_pessimistic_tail_alignment) {
                    int64_t tail_length = cluster[i].first->begin - alignment.sequence().begin();
                    backward_dist[i] = tail_length + pessimistic_gap(tail_length, pessimistic_gap_multiplier);
                }
                else {
                    backward_dist[i] = aligner->longest_detectable_gap(alignment, cluster[i].first->begin);
                }
            }
            else {
                // there is another MEM leftward
                int64_t between_length = max<int64_t>(0, cluster[i].first->begin - cluster[order[idx - 1]].first->end);
                backward_dist[i] = between_length + pessimistic_gap(between_length, pessimistic_gap_multiplier);
            }
            
            if (idx + 1 == cluster.size()) {
                // this is the right tail
                if (use_pessimistic_tail_alignment) {
                    int64_t tail_length = alignment.sequence().end() - cluster[i].first->end;
                    forward_dist[i] = tail_length + pessimistic_gap(tail_length, pessimistic_gap_multiplier) + cluster[i].first->length();
                }
                else {
                    forward_dist[i] = aligner->longest_detectable_gap(alignment, cluster[i].first->end) + cluster[i].first->length();
                }
            }
            else {
                // there is another MEM rightward
                int64_t between_length = max<int64_t>(0, cluster[order[idx + 1]].first->begin - cluster[i].first->end);
                forward_dist[i] = between_length + pessimistic_gap(between_length, pessimistic_gap_multiplier) + cluster[i].first->length();
            }
            
            positions[i] = cluster[i].second;
        }
        
        // expand the restrained search distances until we extract a connected graph or
        // expand the distances up to the maximum detectable length
        
        unique_ptr<bdsg::HashGraph> cluster_graph;
        bool do_extract = true;
        bool connected = false;
        while (do_extract) {
            
            // get rid of the old graph (if there is one)
            cluster_graph = unique_ptr<bdsg::HashGraph>(new bdsg::HashGraph());
            
            // extract according to the current search distances
            algorithms::extract_containing_graph(xindex, cluster_graph.get(), positions, forward_dist, backward_dist,
                                                 num_alt_alns > 1 ? reversing_walk_length : 0);
            
            // we can avoid a costly algorithm when the cluster was extracted from one position (and therefore
            // must be connected)
            if (cluster.size() == 1 || handlealgs::is_weakly_connected(cluster_graph.get())) {
                // we consider enough of the graph extracted once it is connected
                // stop doing further exttraction
                do_extract = false;
                connected = true;
            }
            else {
                // double the search distances, up to the maximum detectable gap
                bool any_dists_changed = false;
                for (size_t i = 0; i < cluster.size(); ++i) {
                    size_t bwd_dist = min<size_t>(backward_dist[i] * 2,
                                                  aligner->longest_detectable_gap(alignment, cluster[i].first->begin));
                    size_t fwd_dist = min<size_t>(forward_dist[i] * 2,
                                                  aligner->longest_detectable_gap(alignment, cluster[i].first->end) + cluster[i].first->length());
                    if (bwd_dist > backward_dist[i]) {
                        backward_dist[i] = bwd_dist;
                        any_dists_changed = true;
                    }
                    if (fwd_dist > forward_dist[i]) {
                        forward_dist[i] = fwd_dist;
                        any_dists_changed = true;
                    }
                }
                // do another extraction as long as we increased at least one search distance
                do_extract = any_dists_changed;
            }
        }
        return move(make_pair(move(cluster_graph), connected));
    }

    pair<unique_ptr<bdsg::HashGraph>, bool> MultipathMapper::extract_cluster_graph(const Alignment& alignment,
                                                                                   const memcluster_t& mem_cluster) const {
        if (restrained_graph_extraction) {
            return extract_restrained_graph(alignment, mem_cluster);
        }
        else {
            return extract_maximal_graph(alignment, mem_cluster);
        }
    }
    
    vector<MultipathMapper::clustergraph_t> MultipathMapper::query_cluster_graphs(const Alignment& alignment,
                                                                                  const vector<MaximalExactMatch>& mems,
                                                                                  const vector<memcluster_t>& clusters) const {
        
        // some settings want us to not merge clusters that have overlapping nodes, and
        // we can also save some bookkeeping work if we neglect to do cluster merging
        // when there's only one cluster anyway
        bool do_merge_suppression = suppress_cluster_merging || clusters.size() <= 1;
        
        // We populate this with all the cluster graphs.
        vector<clustergraph_t> cluster_graphs_out;
        
        // unless suppressing cluster merging, we will ensure that nodes are in only one
        // cluster and we use this to record which one
        unordered_map<id_t, size_t> node_id_to_cluster;
        
        // to hold the clusters as they are (possibly) merged, bools indicate
        // whether we've verified that the graph is connected
        // doubles are the cluster graph's multiplicity
        unordered_map<size_t, tuple<unique_ptr<bdsg::HashGraph>, bool, double>> cluster_graphs;
        
        // to keep track of which clusters have been merged
        UnionFind union_find(clusters.size(), false);
        
        // (for the suppressed merge code path)
        // maps the hits that make up a cluster to the index of the cluster
        unordered_map<pair<const MaximalExactMatch*, pos_t>, size_t> hit_to_cluster;
        
        for (size_t i = 0; i < clusters.size(); i++) {
            
#ifdef debug_multipath_mapper
            cerr << "extracting subgraph for cluster " << i << endl;
#endif
            
            // gather the parameters for subgraph extraction from the MEM hits
            auto& cluster = clusters[i];
            auto extracted = extract_cluster_graph(alignment, cluster);
            tuple<unique_ptr<bdsg::HashGraph>, bool, double> cluster_graph(move(extracted.first), extracted.second, cluster.second);
            
            // check if this subgraph overlaps with any previous subgraph (indicates a probable clustering failure where
            // one cluster was split into multiple clusters)
            unordered_set<size_t> overlapping_graphs;
            
            if (!do_merge_suppression) {
                get<0>(cluster_graph)->for_each_handle([&](const handle_t& handle) {
                    id_t node_id = get<0>(cluster_graph)->get_id(handle);
                    if (node_id_to_cluster.count(node_id)) {
                        overlapping_graphs.insert(node_id_to_cluster[node_id]);
                    }
                    else {
                        node_id_to_cluster[node_id] = i;
                    }
                    return true;
                });
            }
            else {
                // assign the hits to clusters
                for (auto& mem_hit : cluster.first) {
                    hit_to_cluster[mem_hit] = i;
                }
            }
            
            if (overlapping_graphs.empty()) {
                // there is no overlap with any other graph, suggesting a new unique hit
                
#ifdef debug_multipath_mapper
                cerr << "cluster graph does not overlap with any other cluster graphs, adding as cluster " << i << endl;
#endif
                cluster_graphs[i] = move(cluster_graph);
            }
            else {
                // this graph overlaps at least one other graph, so we merge them into one
                
#ifdef debug_multipath_mapper
                cerr << "cluster graph overlaps with the following graphs:" << endl;
                for (auto idx : overlapping_graphs) {
                    cerr << "\t" << idx << endl;
                }
#endif
                
                // merge the groups and decide one graph to stay in the record
                for (size_t j : overlapping_graphs) {
                    union_find.union_groups(i, j);
                }
                size_t remaining_idx = union_find.find_group(i);
                
#ifdef debug_multipath_mapper
                cerr << "merging as cluster " << remaining_idx << endl;
#endif
                
                bdsg::HashGraph* merging_graph;
                bool all_connected;
                double multiplicity;
                if (remaining_idx == i) {
                    // the new graph was chosen to remain, so add it to the record
                    cluster_graphs[i] = move(cluster_graph);
                    merging_graph = get<0>(cluster_graph).get();
                    all_connected = get<1>(cluster_graph);
                    multiplicity = get<2>(cluster_graph);
                }
                else {
                    // the new graph will be merged into an existing graph
                    merging_graph = get<0>(cluster_graphs[remaining_idx]).get();
                    
                    // add in the new graph
                    handlealgs::extend(get<0>(cluster_graph).get(), merging_graph);
                    all_connected = get<1>(cluster_graphs[remaining_idx]) && get<1>(cluster_graph);
                    multiplicity = min(get<2>(cluster_graphs[remaining_idx]), get<2>(cluster_graph));
                }
                
                // merge any other chained graphs into the remaining graph
                for (size_t j : overlapping_graphs) {
                    if (j != remaining_idx) {
                        auto removing_graph = move(cluster_graphs[j]);
                        handlealgs::extend(get<0>(removing_graph).get(), merging_graph);
                        all_connected = all_connected && get<1>(removing_graph);
                        multiplicity = min(multiplicity, get<2>(removing_graph));
                        cluster_graphs.erase(j);
                    }
                }
                get<1>(cluster_graphs[remaining_idx]) = all_connected;
                get<2>(cluster_graphs[remaining_idx]) = multiplicity;
                
                // update the node-to-cluster mapping
                merging_graph->for_each_handle([&](const handle_t& handle) {
                    node_id_to_cluster[merging_graph->get_id(handle)] = remaining_idx;
                    return true;
                });
            }
        }
        
        // check if any cluster graphs pulled disconnected components (as a result of a clustering failure)
        // and if so split them up
        
        // keeps track of the connected components of any multicomponent graph
        vector<pair<size_t, vector<unordered_set<id_t>>>> multicomponent_graphs;
        unordered_map<size_t, vector<size_t>> multicomponent_splits;
        
        size_t max_graph_idx = 0;
        for (const auto& cluster_graph : cluster_graphs) {
            if (!get<1>(cluster_graph.second)) {
                vector<unordered_set<id_t>> connected_components = handlealgs::weakly_connected_components(get<0>(cluster_graph.second).get());
                if (connected_components.size() > 1) {
                    multicomponent_graphs.emplace_back(cluster_graph.first, std::move(connected_components));
                }
            }
            max_graph_idx = max(cluster_graph.first, max_graph_idx);
        }
            
        // did we find any graphs that consist of disconnected components?
        for (pair<size_t, vector<unordered_set<id_t>>>& multicomponent_graph : multicomponent_graphs) {
            max_graph_idx++;
            // make a new graph for each of the components
#ifdef debug_multipath_mapper
            cerr << "cluster graph " << multicomponent_graph.first << " has multiple connected components, splitting now" << endl;
            for (size_t i = 0; i < multicomponent_graph.second.size(); i++) {
                cerr << "component " << max_graph_idx + i << ":" << endl;
                
                for (auto node_id : multicomponent_graph.second[i]) {
                    cerr << "\t" << node_id << endl;
                }
            }
#endif
            
            for (size_t i = 0; i < multicomponent_graph.second.size(); i++) {
                cluster_graphs[max_graph_idx + i] = make_tuple(unique_ptr<bdsg::HashGraph>(new bdsg::HashGraph()), true,
                                                               get<2>(cluster_graphs[multicomponent_graph.first]));
            }
            
            // divvy up the nodes
            auto joined_graph = get<0>(cluster_graphs[multicomponent_graph.first]).get();
            joined_graph->for_each_handle([&](const handle_t& handle) {
                for (size_t j = 0; j < multicomponent_graph.second.size(); j++) {
                    if (multicomponent_graph.second[j].count(joined_graph->get_id(handle))) {
                        get<0>(cluster_graphs[max_graph_idx + j])->create_handle(joined_graph->get_sequence(handle),
                                                                                 joined_graph->get_id(handle));
                        // if we're suppressing cluster merging, we don't maintain this index
                        if (!do_merge_suppression) {
                            node_id_to_cluster[joined_graph->get_id(handle)] = max_graph_idx + j;
                        }
                        break;
                    }
                }
                return true;
            });
            
            // divvy up the edges
            joined_graph->for_each_edge([&](const edge_t& edge) {
                for (size_t j = 0; j < multicomponent_graph.second.size(); j++) {
                    if (multicomponent_graph.second[j].count(joined_graph->get_id(edge.first))) {
                        auto comp_graph = get<0>(cluster_graphs[max_graph_idx + j]).get();
                        comp_graph->create_edge(comp_graph->get_handle(joined_graph->get_id(edge.first),
                                                                       joined_graph->get_is_reverse(edge.first)),
                                                comp_graph->get_handle(joined_graph->get_id(edge.second),
                                                                       joined_graph->get_is_reverse(edge.second)));
                        break;
                    }
                }
                return true;
            });
            
            // remove the old graph
            cluster_graphs.erase(multicomponent_graph.first);
            
            if (do_merge_suppression) {
                // we need to re-assign the hits to the new cluster graphs
                for (auto& mem_hit : clusters[multicomponent_graph.first].first) {
                    for (size_t i = 0; i < multicomponent_graph.second.size(); i++) {
                        if (multicomponent_graph.second[i].count(id(mem_hit.second))) {
                            hit_to_cluster[mem_hit] = max_graph_idx + i;
                            break;
                        }
                    }
                }
            }
            
            max_graph_idx += multicomponent_graph.second.size();
        }
        
        // move the pointers to the return vector and figure out which graph in the return
        // vector each MEM cluster ended up in
        cluster_graphs_out.reserve(cluster_graphs.size());
        unordered_map<size_t, size_t> cluster_to_idx;
        for (pair<const size_t, tuple<unique_ptr<bdsg::HashGraph>, bool, double>>& cluster_graph : cluster_graphs) {
            cluster_to_idx[cluster_graph.first] = cluster_graphs_out.size();
            cluster_graphs_out.emplace_back();
            get<0>(cluster_graphs_out.back()) = move(get<0>(cluster_graph.second));
            get<1>(cluster_graphs_out.back()).second = get<2>(cluster_graph.second);
#ifdef debug_multipath_mapper
            cerr << "adding cluster graph " << cluster_graph.first << " to return vector at index " << cluster_graphs_out.size() << " with multiplicity " << get<1>(cluster_graphs_out.back()).second << endl;
#endif
        }

        
#ifdef debug_multipath_mapper
        cerr << "computing MEM assignments to cluster graphs" << endl;
#endif
        if (!do_merge_suppression) {
            // which MEMs are in play for which cluster?
            for (const MaximalExactMatch& mem : mems) {
                for (gcsa::node_type hit : mem.nodes) {
                    id_t node_id = gcsa::Node::id(hit);
                    if (node_id_to_cluster.count(node_id)) {
                        size_t cluster_idx = cluster_to_idx[node_id_to_cluster[node_id]];
                        get<1>(cluster_graphs_out[cluster_idx]).first.push_back(make_pair(&mem, make_pos_t(hit)));
#ifdef debug_multipath_mapper
                        cerr << "\tMEM " << mem.sequence() << " at " << make_pos_t(hit) << " found in cluster " << node_id_to_cluster[node_id] << " at index " << cluster_idx << endl;
#endif
                    }
                }
            }
        }
        else {
            // we haven't been maintaining the node ID to cluster index (since are not enforcing
            // that each node is in on cluster), so we do something analogous here
            
            // TODO: kinda dumb redundant code
            
#ifdef debug_multipath_mapper
            cerr << "suppressed merging path, creating index to identify nodes with cluster graphs" << endl;
#endif
            
            // identify all of the clusters that contain each node
            unordered_map<id_t, vector<size_t>> node_id_to_cluster_idxs;
            for (size_t i = 0; i < cluster_graphs_out.size(); i++) {
                auto cluster_graph = get<0>(cluster_graphs_out[i]).get();
                cluster_graph->for_each_handle([&](const handle_t& handle){
                    node_id_to_cluster_idxs[cluster_graph->get_id(handle)].push_back(i);
                    return true;
                });
            }
            
            for (const MaximalExactMatch& mem : mems) {
                for (gcsa::node_type hit : mem.nodes) {
                    auto mem_hit = make_pair(&mem, make_pos_t(hit));
                    // force the hits that generated a cluster to be assigned to it
                    auto iter = hit_to_cluster.find(mem_hit);
                    if (iter != hit_to_cluster.end()) {
                        get<1>(cluster_graphs_out[cluster_to_idx[iter->second]]).first.push_back(mem_hit);
#ifdef debug_multipath_mapper
                        cerr << "\tMEM " << mem.sequence() << " at " << mem_hit.second << " assigned as seed to cluster at index " << cluster_to_idx[iter->second] << endl;
#endif
                    }
                    else {
                        // also try to find other, unassigned MEMs in the clusters
                        auto id_iter = node_id_to_cluster_idxs.find(id(mem_hit.second));
                        if (id_iter != node_id_to_cluster_idxs.end()) {
                            for (size_t cluster_idx : id_iter->second) {
                                get<1>(cluster_graphs_out[cluster_idx]).first.push_back(mem_hit);
#ifdef debug_multipath_mapper
                                cerr << "\tMEM " << mem.sequence() << " at " << mem_hit.second << " found in cluster at index " << cluster_idx << endl;
#endif
                            }
                        }
                    }
                }
            }
        }
        
        // compute the read coverage of each cluster graph and sort the assigned MEMs by length
        // and then lexicographically by read index
        for (size_t i = 0; i < cluster_graphs_out.size(); i++) {
            set_read_coverage(cluster_graphs_out[i]);
#ifdef debug_multipath_mapper
            cerr << "compute read coverage of cluster at index " << i << " to be " << get<2>(cluster_graphs_out[i]) << endl;
#endif
        }
            
        // find the node ID range for the cluster graphs to help set up a stable, system-independent ordering
        // note: technically this is not quite a total ordering, but it should be close to one
        unordered_map<bdsg::HashGraph*, uint64_t> graph_hash;
        graph_hash.reserve(cluster_graphs_out.size());
        for (const auto& cluster_graph : cluster_graphs_out) {
            graph_hash[get<0>(cluster_graph).get()] = wang_hash<pair<nid_t, nid_t>>()(make_pair(get<0>(cluster_graph)->min_node_id(),
                                                                                                get<0>(cluster_graph)->max_node_id()));
        }
        
        // sort the cluster graphs descending by unique sequence coverage, breaking ties by scrambling according to a hash
        sort(cluster_graphs_out.begin(), cluster_graphs_out.end(), [&](const clustergraph_t& cluster_graph_1,
                                                                       const clustergraph_t& cluster_graph_2) {
            return (get<2>(cluster_graph_1) > get<2>(cluster_graph_2) ||
                    (get<2>(cluster_graph_1) == get<2>(cluster_graph_2) &&
                     graph_hash[get<0>(cluster_graph_1).get()] < graph_hash[get<0>(cluster_graph_2).get()]));
        });
        
        return cluster_graphs_out;
    }

    bool MultipathMapper::expand_for_softclips(clustergraph_t& cluster_graph,
                                               const multipath_alignment_t& multipath_aln) const {
        
        // TODO: this post-hoc check might make it possible to extract less graph initially...
        // but that could get fiddly when there is a mismatch at the graph boundary (makes soft clip)
        // TODO: this should work on non-branching path tips, not node tips
        
        unordered_map<handle_t, int64_t> expansion_points;
        auto& graph = *get<0>(cluster_graph);
        
        // find expansion points on source subpaths
        for (auto i : multipath_aln.start()) {
            const auto& mapping = multipath_aln.subpath(i).path().mapping().front();
            const auto& edit = mapping.edit().front();
            if (edit.from_length() == 0 && edit.to_length() != 0) {
                // this is soft-clipped
                size_t max_len = get_aligner()->longest_detectable_gap(multipath_aln.sequence().size(), edit.to_length());
                if (mapping.position().offset() < max_len) {
                    // the clipped sequence's alignment might not fit on this node
                    handle_t handle = graph.get_handle(mapping.position().node_id(), mapping.position().is_reverse());
                    bool is_source = graph.follow_edges(handle, true, [](const handle_t& h) {return false;});
                    if (is_source) {
                        // this soft clip might arise because we didn't extract enough graph
                        auto flipped = graph.flip(handle);
                        expansion_points[flipped] =  max<int64_t>(expansion_points[flipped],
                                                                  max_len - mapping.position().offset());
                    }
                }
            }
        }
        
        // find expansion points on sink subpaths
        for (const auto& subpath : multipath_aln.subpath()) {
            if (subpath.next().empty()) {
                const auto& mapping = subpath.path().mapping().back();
                const auto& edit = mapping.edit().back();
                if (edit.from_length() == 0 && edit.to_length() != 0) {
                    // this is soft-clipped
                    handle_t handle = graph.get_handle(mapping.position().node_id(), mapping.position().is_reverse());
                    int64_t final_offset = mapping.position().offset() + mapping_from_length(mapping);
                    size_t max_len = get_aligner()->longest_detectable_gap(multipath_aln.sequence().size(), edit.to_length());
                    if (graph.get_length(handle) - final_offset < max_len) {
                        // the clipped sequence's alignment might not fit on this node
                        bool is_sink = graph.follow_edges(handle, false, [](const handle_t& h) {return false;});
                        if (is_sink) {
                            // this soft clip might arise because we didn't extract enough graph
                            expansion_points[handle] =  max<int64_t>(expansion_points[handle],
                                                                     max_len - (graph.get_length(handle) - final_offset));
                        }
                    }
                }
            }
        }
        
        if (!expansion_points.empty()) {
            int64_t num_edges_before_expansion = graph.get_edge_count();
            for (const pair<handle_t, int64_t>& expansion_point : expansion_points) {
                // we might need to extend the subgraph here to get a full alignment
                algorithms::locally_expand_graph(*xindex, graph, expansion_point.first, expansion_point.second);
            }
            return num_edges_before_expansion != graph.get_edge_count();
        }
        else {
            return false;
        }
    }
    
    void MultipathMapper::multipath_align(const Alignment& alignment, clustergraph_t& cluster_graph,
                                          multipath_alignment_t& multipath_aln_out,
                                          const match_fanouts_t* fanouts) const {

        // we put this in a loop so that we can check to make sure we didn't miss part of the
        // alignment because there wasn't enough graph extracted
        bool new_graph_material = true;
        while (new_graph_material) {
            // there are parts of this graph we haven't tried to align to yet
            
            multipath_aln_out.mutable_subpath()->clear();
            multipath_aln_out.mutable_start()->clear();
            
            auto graph = get<0>(cluster_graph).get();
            auto& graph_mems = get<1>(cluster_graph);
            
#ifdef debug_multipath_mapper_alignment
            cerr << "constructing alignment graph for cluster of " << get<1>(cluster_graph).first.size() << " hits" << endl;
#endif
            
            if (graph_mems.first.empty()) {
#ifdef debug_multipath_mapper_alignment
                cerr << "cluster is empty, aborting" << endl;
#endif
                transfer_read_metadata(alignment, multipath_aln_out);
                return;
            }
            
            // the longest path we could possibly align to (full gap and a full sequence)
            auto aligner = get_aligner(!alignment.quality().empty());
            size_t target_length = alignment.sequence().size() + min(aligner->longest_detectable_gap(alignment), max_alignment_gap);
            
            // check if we can get away with using only one strand of the graph
            bool use_single_stranded = handlealgs::is_single_stranded(graph);
            bool mem_strand = false;
            if (use_single_stranded) {
                mem_strand = is_rev(graph_mems.first[0].second);
                for (size_t i = 1; i < graph_mems.first.size(); i++) {
                    if (is_rev(graph_mems.first[i].second) != mem_strand) {
                        use_single_stranded = false;
                        break;
                    }
                }
            }
            
            // make the graph we need to align to
#ifdef debug_multipath_mapper_alignment
            cerr << "use_single_stranded: " << use_single_stranded << " mem_strand: " << mem_strand << endl;
#endif
            
#ifdef debug_multipath_mapper_alignment
            cerr << "initial alignment graph with " << graph->get_node_count() << " nodes and " << graph->get_edge_count() << ":" << endl;
            graph->for_each_handle([&](const handle_t& h) {
                cerr << graph->get_id(h) << " " << graph->get_sequence(h) << endl;
                graph->follow_edges(h, false, [&](const handle_t& n) {
                    cerr << "\t-> " << graph->get_id(n) << " " << (graph->get_is_reverse(n) ? "-" : "+") << endl;
                });
                graph->follow_edges(h, true, [&](const handle_t& n) {
                    cerr << "\t " << graph->get_id(n) << " " << (graph->get_is_reverse(n) ? "-" : "+") << " <-" << endl;
                });
            });
#endif
            
            // make our options for single stranded graphs
            IdentityOverlay fwd_graph(graph);
            ReverseGraph rev_graph(graph, true);
            StrandSplitGraph split_graph(graph);
            
            // choose which one we want
            ExpandingOverlayGraph* align_digraph = nullptr;
            if (!use_single_stranded) {
                align_digraph = &split_graph;
            }
            else if (mem_strand) {
                align_digraph = &rev_graph;
            }
            else {
                align_digraph = &fwd_graph;
            }
            
            // if necessary, convert from cyclic to acylic (can be expensive to construct, only do it
            // if we need to)
            
            IdentityOverlay undagified(align_digraph);
            unique_ptr<DagifiedGraph> dagified;
            
            ExpandingOverlayGraph* align_dag = nullptr;
            if (handlealgs::is_directed_acyclic(align_digraph)) {
                align_dag = &undagified;
            }
            else {
#ifdef debug_multipath_mapper_alignment
                cerr << "graph contains directed cycles, performing dagification" << endl;
#endif
                dagified = unique_ptr<DagifiedGraph>(new DagifiedGraph(align_digraph, target_length, max_dagify_duplications));
                align_dag = dagified.get();
            }
            
            // a function to translate from the transformed graphs ID space to the original graph's
            function<pair<id_t, bool>(id_t)> translator = [&](const id_t node_id) {
                handle_t original = align_digraph->get_underlying_handle(align_dag->get_underlying_handle(align_dag->get_handle(node_id)));
                return make_pair(graph->get_id(original), graph->get_is_reverse(original));
            };
            
#ifdef debug_multipath_mapper_alignment
            cerr << "final alignment graph of size " << align_dag->get_node_count() << " nodes and " << align_dag->get_edge_count() << ":" << endl;
            align_dag->for_each_handle([&](const handle_t& h) {
                auto tr = translator(align_dag->get_id(h));
                cerr << align_dag->get_id(h) << " (" << tr.first << (tr.second ? "-" : "+") << ") " << align_dag->get_sequence(h) << endl;
                align_dag->follow_edges(h, false, [&](const handle_t& n) {
                    cerr << "\t-> " << align_dag->get_id(n) << " " << (align_dag->get_is_reverse(n) ? "-" : "+") << endl;
                });
                align_dag->follow_edges(h, true, [&](const handle_t& n) {
                    cerr << "\t " << align_dag->get_id(n) << " " << (align_dag->get_is_reverse(n) ? "-" : "+") << " <-" << endl;
                });
            });
#endif
            
            // construct a graph that summarizes reachability between MEMs
            
#ifdef debug_multipath_mapper_alignment
            cerr << "making multipath alignment MEM graph" << endl;
#endif
            
            vector<size_t> hit_provenance;
            MultipathAlignmentGraph multi_aln_graph(*align_dag, graph_mems, translator, hit_provenance,
                                                    max_branch_trim_length, gcsa, fanouts);
            
            {
                // Compute a topological order over the graph
                vector<size_t> topological_order;
                multi_aln_graph.topological_sort(topological_order);
                
                // it's sometimes possible for transitive edges to survive the original construction algorithm, so remove them
                multi_aln_graph.remove_transitive_edges(topological_order);
                
                // prune this graph down the paths that have reasonably high likelihood
                size_t size_before_prune = multi_aln_graph.size();
                multi_aln_graph.prune_to_high_scoring_paths(alignment, aligner, max_suboptimal_path_score_ratio,
                                                            topological_order, hit_provenance);
                
                if (multi_aln_graph.size() != size_before_prune && do_spliced_alignment) {
                    // we pruned away some path nodes, so let's check if we pruned away any entire hits
                    // and, if so, un-claim them from this cluster
                    
                    vector<bool> found(graph_mems.first.size(), false);
                    for (auto i : hit_provenance) {
                        found[i] = true;
                    }
                    
                    // if necessary remove the ones we didn't keep
                    size_t removed = 0;
                    for (size_t i = 0; i < graph_mems.first.size(); ++i) {
                        if (!found[i]) {
#ifdef debug_multipath_alignment
                            cerr << "completely pruned hit from cluster: " << graph_mems.first[i].second << " " << graph_mems.first[i].first->sequence() << endl;
#endif
                            ++removed;
                        }
                        else if (removed) {
                            graph_mems.first[i - removed] = graph_mems.first[i];
                        }
                    }
                    graph_mems.first.resize(graph_mems.first.size() - removed);
                    
                    // we may need to recompute the coverage of the cluster because some MEMs were pruned out of it
                    set_read_coverage(cluster_graph);
                }
            }
            
            if (snarl_manager || distance_index) {
                // We want to do snarl cutting
                
                if (!suppress_tail_anchors) {
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "Synthesizing tail anchors for snarl cutting" << endl;
#endif
                    
                    // Make fake anchor paths to cut the snarls out of in the tails
                    multi_aln_graph.synthesize_tail_anchors(alignment, *align_dag, aligner, min_tail_anchor_length, num_alt_alns,
                                                            false, max_alignment_gap,
                                                            use_pessimistic_tail_alignment ? pessimistic_gap_multiplier : 0.0);
                    
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << "MultipathAlignmentGraph going into snarl cutting:" << endl;
                multi_aln_graph.to_dot(cerr, &alignment);
#endif
                
                // Do the snarl cutting, which modifies the nodes in the multipath alignment graph
                if (max_snarl_cut_size) {
                    multi_aln_graph.resect_snarls_from_paths(snarl_manager, distance_index, translator, max_snarl_cut_size);
                }
            }
            
            
#ifdef debug_multipath_mapper_alignment
            cerr << "MultipathAlignmentGraph going into alignment:" << endl;
            multi_aln_graph.to_dot(cerr, &alignment);
            
            for (auto& ids : multi_aln_graph.get_connected_components()) {
                cerr << "Component: ";
                for (auto& id : ids) {
                    cerr << id << " ";
                }
                cerr << endl;
            }
#endif
            
            function<size_t(const Alignment&, const HandleGraph&)> choose_band_padding = [&](const Alignment& seq, const HandleGraph& graph) {
                size_t read_length = seq.sequence().size();
                return read_length < band_padding_memo.size() ? band_padding_memo.at(read_length)
                                                              : size_t(band_padding_multiplier * sqrt(read_length)) + 1;
            };
            
            // do the connecting alignments and fill out the multipath_alignment_t object
            multi_aln_graph.align(alignment, *align_dag, aligner, true, num_alt_alns, dynamic_max_alt_alns, max_alignment_gap,
                                  use_pessimistic_tail_alignment ? pessimistic_gap_multiplier : 0.0, simplify_topologies,
                                  max_tail_merge_supress_length, choose_band_padding, multipath_aln_out, snarl_manager,
                                  distance_index, &translator);
            
            // Note that we do NOT topologically order the multipath_alignment_t. The
            // caller has to do that, after it is finished breaking it up into
            // connected components or whatever.
            
#ifdef debug_multipath_mapper_alignment
            cerr << "multipath alignment before translation: " << debug_string(multipath_aln_out) << endl;
#endif
            for (size_t j = 0; j < multipath_aln_out.subpath_size(); j++) {
                translate_oriented_node_ids(*multipath_aln_out.mutable_subpath(j)->mutable_path(), translator);
            }
            
#ifdef debug_multipath_mapper_alignment
            cerr << "completed multipath alignment: " << debug_string(multipath_aln_out) << endl;
#endif
            
            // check if we need to expand the graph to get a full alignment
            new_graph_material = expand_for_softclips(cluster_graph, multipath_aln_out);
#ifdef debug_multipath_mapper_alignment
            if (new_graph_material) {
                cerr << "found soft clip near graph boundary, expanding subraph and re-aligning" << endl;
            }
#endif
        }
    }
            
    void MultipathMapper::make_nontrivial_multipath_alignment(const Alignment& alignment, const HandleGraph& subgraph,
                                                              const function<pair<id_t, bool>(id_t)>& translator,
                                                              multipath_alignment_t& multipath_aln_out) const {
        
#ifdef debug_multipath_mapper_alignment
        cerr << "attempting to make nontrivial alignment for " << alignment.name() << endl;
#endif
        
        auto aligner = get_aligner(!alignment.quality().empty());
        
        // create an alignment graph with the internals of snarls removed
        MultipathAlignmentGraph multi_aln_graph(subgraph, alignment, snarl_manager, distance_index, max_snarl_cut_size, translator);
        
        // remove any transitive edges that may have found their way in there
        // TODO: is this necessary? all edges should be across snarls, how could they be transitive? from trimmed indels maybe?
        // in any case, the optimization on the transitive edge algorithm will make this linear if there actually aren't any
        vector<size_t> topological_order;
        multi_aln_graph.topological_sort(topological_order);
        multi_aln_graph.remove_transitive_edges(topological_order);
        
        function<size_t(const Alignment&, const HandleGraph&)> choose_band_padding = [&](const Alignment& seq, const HandleGraph& graph) {
            size_t read_length = seq.sequence().end() - seq.sequence().begin();
            return read_length < band_padding_memo.size() ? band_padding_memo.at(read_length)
                                                          : size_t(band_padding_multiplier * sqrt(read_length)) + 1;
        };
        
        // do the connecting alignments and fill out the multipath_alignment_t object
        multi_aln_graph.align(alignment, subgraph, aligner, false, num_alt_alns, dynamic_max_alt_alns, max_alignment_gap,
                              use_pessimistic_tail_alignment ? pessimistic_gap_multiplier : 0.0, simplify_topologies,
                              max_tail_merge_supress_length, choose_band_padding, multipath_aln_out);
        
        for (size_t j = 0; j < multipath_aln_out.subpath_size(); j++) {
            translate_oriented_node_ids(*multipath_aln_out.mutable_subpath(j)->mutable_path(), translator);
        }
        
        topologically_order_subpaths(multipath_aln_out);
        
#ifdef debug_multipath_mapper_alignment
        cerr << "completed multipath alignment: " << debug_string(multipath_aln_out) << endl;
#endif
    }
    
    void MultipathMapper::set_read_coverage(clustergraph_t& cluster_graph) {
        
        auto& mem_hits = get<1>(cluster_graph);
        if (mem_hits.first.empty()) {
            get<2>(cluster_graph) = 0;
            return;
        }
        
        // lexicographic comparison for the sort this algorithm needs
        auto lex_cmp = [](const pair<const MaximalExactMatch*, pos_t>& hit_1,
                         const pair<const MaximalExactMatch*, pos_t>& hit_2) {
            return (hit_1.first->begin < hit_2.first->begin ||
                    (hit_1.first->begin == hit_2.first->begin && hit_1.first->end < hit_2.first->end));
        };
        
        // length first comparison for later steps
        auto length_first_cmp = [](const pair<const MaximalExactMatch*, pos_t>& hit_1,
                                   const pair<const MaximalExactMatch*, pos_t>& hit_2) {
            return (hit_1.first->length() > hit_2.first->length() ||
                    (hit_1.first->length() == hit_2.first->length() &&
                     (hit_1.first->begin < hit_2.first->begin ||
                      (hit_1.first->begin == hit_2.first->begin && hit_1.first->end < hit_2.first->end))));
        };
        
        
        if (!is_sorted(get<1>(cluster_graph).first.begin(), get<1>(cluster_graph).first.end(), lex_cmp)) {
            stable_sort(get<1>(cluster_graph).first.begin(), get<1>(cluster_graph).first.end(), lex_cmp);
        }
        
        vector<pair<string::const_iterator, string::const_iterator>> mem_read_segments;
        mem_read_segments.reserve(mem_hits.first.size());
        for (auto& mem_hit : mem_hits.first) {
            mem_read_segments.emplace_back(mem_hit.first->begin, mem_hit.first->end);
        }
        get<2>(cluster_graph) = algorithms::count_covered(mem_read_segments);
        
        stable_sort(get<1>(cluster_graph).first.begin(), get<1>(cluster_graph).first.end(), length_first_cmp);
    }
    
    void MultipathMapper::strip_full_length_bonuses(multipath_alignment_t& multipath_aln) const {
        
        // TODO: this could technically be wrong if only one read in a pair has qualities
        const auto& aligner = *get_aligner(!multipath_aln.quality().empty());
        // strip bonus from source paths
        if (multipath_aln.start_size()) {
            // use the precomputed list of sources if we have it
            for (size_t i = 0; i < multipath_aln.start_size(); i++) {
                subpath_t* source_subpath = multipath_aln.mutable_subpath(multipath_aln.start(i));
                const edit_t& edit = source_subpath->path().mapping(0).edit(0);
                if (edit.to_length() != 0 && edit.from_length() != 0) {
                    source_subpath->set_score(source_subpath->score()
                                              - aligner.score_full_length_bonus(true, multipath_aln.sequence().begin(),
                                                                                multipath_aln.sequence().end(),
                                                                                multipath_aln.quality().begin()));
                }
            }
        }
        else {
            // find sources
            vector<bool> is_source(multipath_aln.subpath_size(), true);
            for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
                const subpath_t& subpath = multipath_aln.subpath(i);
                for (size_t j = 0; j < subpath.next_size(); j++) {
                    is_source[subpath.next(j)] = false;
                }
            }
            // strip the bonus from the sources
            for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
                if (!is_source[i]) {
                    continue;
                }
                subpath_t* source_subpath = multipath_aln.mutable_subpath(i);
                const edit_t& edit = source_subpath->path().mapping(0).edit(0);
                if (edit.to_length() != 0 && edit.from_length() != 0) {
                    source_subpath->set_score(source_subpath->score()
                                              - aligner.score_full_length_bonus(true, multipath_aln.sequence().begin(),
                                                                                multipath_aln.sequence().end(),
                                                                                multipath_aln.quality().begin()));
                }
            }
        }
        // strip bonus from sink paths
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            subpath_t* subpath = multipath_aln.mutable_subpath(i);
            if (subpath->next_size() == 0) {
                const path_mapping_t& final_mapping = subpath->path().mapping(subpath->path().mapping_size() - 1);
                const edit_t& edit = final_mapping.edit(final_mapping.edit_size() - 1);
                if (edit.to_length() != 0 && edit.from_length() != 0) {
                    subpath->set_score(subpath->score()
                                       - aligner.score_full_length_bonus(false, multipath_aln.sequence().begin(),
                                                                         multipath_aln.sequence().end(),
                                                                         multipath_aln.quality().begin()));
                }
            }
        }
    }


    vector<double> MultipathMapper::mapping_likelihoods(vector<multipath_alignment_t>& multipath_alns) const {
        
        // Only do the population MAPQ if it might disambiguate two paths
        // (since it's not as cheap as just using the score), or if we set the
        // setting to always do it.
        bool include_population_component = (use_population_mapqs && (multipath_alns.size() > 1 || always_check_population));
        // Records whether, for each multipath alignment, at least one of the
        // paths enumerated followed only edges in the index. We count totally
        // unmapped reads as pop consistent, as all 0 edges they cross are pop
        // consistent.
        bool all_multipaths_pop_consistent = true;
        
        double log_base = get_aligner(!multipath_alns.front().quality().empty())->log_base;
        
        // The score of the optimal Alignment for each multipath_alignment_t, not adjusted for population
        vector<double> scores(multipath_alns.size(), 0.0);
        // The scores of the best Alignment for each multipath_alignment_t, adjusted for population.
        // These can be negative but will be bumped up to all be positive later.
        vector<double> pop_adjusted_scores;
        if (include_population_component) {
            pop_adjusted_scores.resize(multipath_alns.size());
        }
        
        // We need to track the score adjustments so we can compensate for
        // negative values, turning the largest penalty into a 0 bonus.
        double min_adjustment = numeric_limits<double>::max();
        
        
        for (size_t i = 0; i < multipath_alns.size(); i++) {
            // Score all the multipath alignment candidates, optionally using
            // population adjustment
            
            // We will query the population database for this alignment if it
            // is turned on and it succeeded for the others.
            bool query_population = include_population_component && all_multipaths_pop_consistent;
            
            /// Get all the linearizations we are going to work with, possibly with duplicates.
            /// The first alignment will be optimal.
            vector<Alignment> alignments;
            int32_t aln_score = -1;
            
            if (query_population) {
                // We want to do population scoring
                if (!top_tracebacks && haplo_score_provider->has_incremental_search()) {
                    // We can use incremental haplotype search to find all the linearizations consistent with haplotypes
                    // Make sure to also always include the optimal alignment first, even if inconsistent.
                    // And also include up to population_max_paths non-consistent but hopefully scorable paths
                    alignments = haplotype_consistent_alignments(multipath_alns[i], *haplo_score_provider, population_max_paths,
                                                                 population_paths_hard_cap, true);
                } else {
                    // We will just find the top n best-alignment-scoring linearizations and hope some match haplotypes
                    alignments = optimal_alignments(multipath_alns[i], population_max_paths);
                }
                
                aln_score = alignments.front().score();
            } else {
                // Just compute a single optimal alignment
                aln_score = optimal_alignment_score(multipath_alns[i]);
            }
            
#ifdef debug_multipath_mapper
            cerr << "Got " << alignments.size() << " tracebacks for multipath " << i << endl;
#endif
#ifdef debug_multipath_mapper_alignment
            cerr << debug_string(multipath_alns[i]) << endl;
#endif
            
            // Now, we may have been fed a multipath_alignment_t where the best
            // single path alignment is to leave it unmapped alltogether. Maybe
            // we cut out a really terrible bit of the alignment graph somehow.
            // We used to fail an assert in that case, but we can handle it as
            // just an unmapped read with score 0.
            
            // Collect the score of the optimal alignment, to use if population
            // scoring fails for a multipath alignment. Put it in the optimal
            // base score.
            scores[i] = max<int32_t>(aln_score, 0);
            
            if (query_population) {
                
                // Work out the population size. Use the override, then try the score provider, and then fall back to the xg.
                auto haplotype_count = force_haplotype_count;
                
                if (haplotype_count == 0) {
                    haplotype_count = haplo_score_provider->get_haplotype_count();
                }
                
                if (haplotype_count == 0 || haplotype_count == -1) {
                    // The score provider doesn't ahve a haplotype count. Fall back to the count in the XG.
                    // No longer available!
                    //haplotype_count = xindex->get_haplotype_count();
                }
                
                if (haplotype_count == 0 || haplotype_count == -1) {
                    // We really should have a haplotype count
                    throw runtime_error("Cannot score any haplotypes with a 0 or -1 haplotype count; are haplotypes available?");
                }
                
                // Make sure to grab the memo
                auto& memo = get_rr_memo(recombination_penalty, haplotype_count);
                
                // Now we need to score the linearizations. The pop-adjusted
                // score of the pop-scorable linearization with the best
                // pop-adjusted score lives in pop_adjusted_scores[i].
                double& best_linearization_total_score = pop_adjusted_scores[i];
                // The population score for that alignment lives here.
                double best_linearization_pop_score = 0;
                // We set this to true if we found a best linearization.
                bool have_best_linearization = false;
                
                for (size_t j = 0; j < alignments.size(); j++) {
                    // Score each alignment if possible
                    auto pop_score = haplo_score_provider->score(alignments[j].path(), memo);
                    
#ifdef debug_multipath_mapper
                    cerr << "Got pop score " << pop_score.first << ", " << pop_score.second << " for alignment " << j
                    << " score " << alignments[j].score() << " of multipath " << i << endl;
#endif
#ifdef debug_multipath_mapper_alignment
                    cerr << pb2json(alignments[j]) << endl;
#endif
                    
                    if (std::isnan(pop_score.first) && pop_score.second) {
                        // This shouldn't happen. Bail out on haplotype adjustment for this read and warn.
                        cerr << "warning:[vg::MultipathMapper]: NAN population score obtained for read "
                        << alignments[j].name() << " with ostensibly successful query. Changing to failure." << endl;
                        pop_score.second = false;
                    }
                    
                    if (std::isnan(alignments[j].score())) {
                        // This is even worse! The alignment score itself is somehow NAN.
                        cerr << "warning:[vg::MultipathMapper]: NAN alignment score obtained in alignment being considered for read "
                        << alignments[j].name() << ". This should never happen! Changing to failure." << endl;
                        pop_score.second = false;
                    }
                    
                    if (pop_score.second) {
                        // If the alignment was pop-score-able, mix it in as a candidate for the best linearization
                        
                        // Compute its pop-adjusted score.
                        // Make sure to account for the aligner's log base to have consistent point values.
                        double total_score = alignments[j].score() + pop_score.first / log_base;
                        
                        if (!have_best_linearization || total_score > best_linearization_total_score) {
                            // This is the new best linearization
                            
                            best_linearization_total_score = total_score;
                            best_linearization_pop_score = pop_score.first / log_base;
                            have_best_linearization = true;
                        }
                        
                    }
                    
                    // Otherwise, skip it
                }
                
                
                if (!have_best_linearization) {
                    // If we have no best linear pop-scored Alignment, bail out on population score correction for this read entirely.
                    // We probably have a placement in a region not covered by the haplotype index at all.
                    all_multipaths_pop_consistent = false;
                    continue;
                }
                
                // Otherwise, we have population scores.
                
#ifdef debug_multipath_mapper
                cerr << "Best population-adjusted linearization score is " << best_linearization_total_score << endl;
#endif
                
                // Save the population score from the best total score Alignment.
                // TODO: This is not the pop score of the linearization that the multipath_alignment_t wants to give us by default.
                multipath_alns[i].set_annotation("haplotype_score", best_linearization_pop_score);
                
                // The multipath's base score is the base score of the
                // best-base-score linear alignment. This is the "adjustment"
                // we apply to the multipath's score to make it match the
                // pop-adjusted score of the best-pop-adjusted-score linear
                // alignment.
                double adjustment = best_linearization_total_score - scores[i];
                
                // See if we have a new minimum adjustment value, for the adjustment applicable to the chosen traceback.
                min_adjustment = min(min_adjustment, adjustment);
            }
        }
        
        if (include_population_component && all_multipaths_pop_consistent) {
            // We will go ahead with pop scoring for this read
            
#ifdef debug_multipath_mapper
            cerr << "Haplotype consistency score is being used." << endl;
#endif
            
            for (auto& score : pop_adjusted_scores) {
                // Adjust the adjusted scores up/down by the minimum adjustment to ensure no scores are negative
                score -= min_adjustment;
            }
            
            for (auto& mpaln : multipath_alns) {
                // Remember that we did use population scoring on all these multipath_alignment_ts
                mpaln.set_annotation("haplotype_score_used", true);
            }
        } else {
            // Clean up pop score annotations and remove scores on all the reads.
            
#ifdef debug_multipath_mapper
            cerr << "Haplotype consistency score is not being used." << endl;
#endif
            
            for (auto& mpaln : multipath_alns) {
                mpaln.clear_annotation("haplotype_score_used");
                mpaln.clear_annotation("haplotype_score");
            }
        }
        
        // Select whether to use base or adjusted scores depending on whether
        // we did population-aware alignment and succeeded for all the
        // multipath alignments.
        if (include_population_component && all_multipaths_pop_consistent) {
            scores = move(pop_adjusted_scores);
        }
        return scores;
    }

    vector<double> MultipathMapper::pair_mapping_likelihoods(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs,
                                                             const vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs) const {
        
#ifdef debug_multipath_mapper
        cerr << "computing paired read likelihoods" << endl;
#endif
        
        // Only do the population MAPQ if it might disambiguate two paths (since it's not
        // as cheap as just using the score), or if we set the setting to always do it.
        bool include_population_component = (use_population_mapqs &&
                                             (multipath_aln_pairs.size() > 1 || always_check_population));
        // Records whether, for each multipath alignment pair, at least one of
        // the paths enumerated for each end followed only edges in the index.
        // We count totally unmapped reads as pop consistent, as all 0 edges
        // they cross are pop consistent.
        bool all_multipaths_pop_consistent = true;
        
        double log_base = get_aligner(!multipath_aln_pairs.front().first.quality().empty() &&
                                      !multipath_aln_pairs.front().second.quality().empty())->log_base;
        
        // the scores of the optimal alignments and fragments, ignoring population
        vector<double> scores(multipath_aln_pairs.size(), 0.0);
        
        // the scores of the optimal alignments and fragments, accounting for population
        vector<double> pop_adjusted_scores;
        if (include_population_component) {
            pop_adjusted_scores.resize(multipath_aln_pairs.size());
        }
        // population + fragment score, for when population adjustment is used, to make scores nonnegative
        double min_extra_score = numeric_limits<double>::max();
        // just fragment score, for running without population adjustment, to make scores nonnegative
        double min_frag_score = numeric_limits<double>::max();
        
        
        for (size_t i = 0; i < multipath_aln_pairs.size(); i++) {
            // For each pair of read placements
            pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair = multipath_aln_pairs[i];
            
            // We will query the population database for this alignment pair if it
            // is turned on and it succeeded for the others.
            bool query_population = include_population_component && all_multipaths_pop_consistent;
            
            // Generate the top alignments on each side, or the top
            // population_max_paths alignments if we are doing multiple
            // alignments for population scoring.
            vector<vector<Alignment>> alignments;
            int32_t aln_score_1 = -1, aln_score_2 = -1;
            
            if (query_population) {
                // We want to do population scoring
                alignments.resize(2);
                if (!top_tracebacks && haplo_score_provider->has_incremental_search()) {
                    // We can use incremental haplotype search to find all the linearizations consistent with haplotypes
                    // Make sure to also always include the optimal alignment first, even if inconsistent.
                    // Also pad out with population_max_paths inconsistent or unscorable paths
                    alignments[0] = haplotype_consistent_alignments(multipath_aln_pair.first, *haplo_score_provider, population_max_paths,
                                                                    population_paths_hard_cap, true);
                    alignments[1] = haplotype_consistent_alignments(multipath_aln_pair.second, *haplo_score_provider, population_max_paths,
                                                                    population_paths_hard_cap, true);
                } else {
                    // We will just find the top n best-alignment-scoring linearizations and hope some match haplotypes
                    alignments[0] = optimal_alignments(multipath_aln_pair.first, population_max_paths);
                    alignments[1] = optimal_alignments(multipath_aln_pair.second, population_max_paths);
                }
                
                if (!alignments[0].empty()) {
                    aln_score_1 = alignments[0].front().score();
                }
                if (!alignments[1].empty()) {
                    aln_score_2 = alignments[1].front().score();
                }
                
#ifdef debug_multipath_mapper
                
                cerr << "Got " << alignments[0].size() << " and " << alignments[1].size() << " linearizations on each end" << endl;
#endif
            } else {
                // Just compute a single optimal alignment
                aln_score_1 = optimal_alignment_score(multipath_aln_pair.first);
                aln_score_2 = optimal_alignment_score(multipath_aln_pair.second);
            }
            
            
            // We used to fail an assert if either list of optimal alignments
            // was empty, but now we handle it as if that side is an unmapped
            // read with score 0.
            
            // Compute the optimal alignment score ignoring population
            int32_t alignment_score = max<int32_t>(aln_score_1, 0) + max<int32_t>(aln_score_2, 0);
            
            // This is the contribution to the alignment's score from the fragment length distribution
            double frag_score;
            if (aln_score_1 == -1 || aln_score_2 == -1) {
                // Actually there should be no fragment score, because one or both ends are unmapped
                frag_score = 0;
            } else {
                // compute the fragment distribution's contribution to the score
                frag_score = fragment_length_log_likelihood(cluster_pairs[i].second) / log_base;
            }
            min_frag_score = min(frag_score, min_frag_score);
            
            // Record the base score, including fragment contribution
            scores[i] = alignment_score + frag_score;
            
            if (query_population) {
                // We also want to select the optimal population-scored alignment on each side and compute a pop-adjusted score.
                
                // Work out the population size. Use the override, then try the score provider, and then fall back to the xg.
                auto haplotype_count = force_haplotype_count;
                
                if (haplotype_count == 0) {
                    haplotype_count = haplo_score_provider->get_haplotype_count();
                }
                
                if (haplotype_count == 0 || haplotype_count == -1) {
                    // The score provider doesn't ahve a haplotype count. Fall back to the count in the XG.
                    //haplotype_count = xindex->get_haplotype_count();
                }
                
                if (haplotype_count == 0 || haplotype_count == -1) {
                    // We really should have a haplotype count
                    throw runtime_error("Cannot score any haplotypes with a 0 or -1 haplotype count; are haplotypes available?");
                }
                
                // Make sure to grab the memo
                auto& memo = get_rr_memo(recombination_penalty, haplotype_count);
                
                // Now we need to score the linearizations.
                
                // This is the best pop-adjusted linearization score for each end.
                double best_total_score[2] = {0, 0};
                // This is the pop score that goes with it
                double best_pop_score[2] = {0, 0};
                // We set this to true if we find a best linearization for each end.
                bool have_best_linearization[2] = {false, false};
                // Note that for unmapped reads, the total and pop scores will stay 0.
                
                for (int end : {0, 1}) {
                    // For each read in the pair
                    
                    for (size_t j = 0; j < alignments[end].size(); j++) {
                        // For each alignment of the read in this location
                        
                        // Pop score the alignment
                        auto pop_score = haplo_score_provider->score(alignments[end][j].path(), memo);
                        
                        if (std::isnan(pop_score.first) && pop_score.second) {
                            // This shouldn't happen. Bail out on haplotype adjustment for this read and warn.
                            cerr << "warning:[vg::MultipathMapper]: NAN population adjusted score obtained for paired read "
                            << alignments[end][j].name() << " with ostensibly successful query. Changing to failure." << endl;
                            pop_score.second = false;
                        }
                        
                        if (std::isnan(alignments[end][j].score())) {
                            // This is even worse! The alignment score itself is somehow NAN.
                            cerr << "warning:[vg::MultipathMapper]: NAN alignment score obtained in alignment being considered for paired read "
                            << alignments[end][j].name() << ". This should never happen! Changing to failure." << endl;
                            pop_score.second = false;
                        }
                        
#ifdef debug_multipath_mapper
                        cerr << "Linearization " << j << " on end " << end << " gets pop score " << pop_score.first
                        << " and alignment score " << alignments[end][j].score() << endl;
#endif
                        
                        if (pop_score.second) {
                            // If the alignment was pop-score-able, mix it in as a candidate for the best linearization
                            
                            // Compute its pop-adjusted score.
                            // Make sure to account for the aligner's log base to have consistent point values.
                            double total_score = alignments[end][j].score() + pop_score.first / log_base;
                            
                            if (!have_best_linearization[end] || total_score > best_total_score[end]) {
                                // This is the new best linearization
                                
                                best_total_score[end] = total_score;
                                best_pop_score[end] = pop_score.first / log_base;
                                have_best_linearization[end] = true;
                            }
                            
                        }
                    }
                }
                
                if ((!alignments[0].empty() && !have_best_linearization[0]) ||
                    (!alignments[1].empty() && !have_best_linearization[1])) {
                    // If we couldn't find a linearization for each mapped end that we could score, bail on pop scoring.
                    all_multipaths_pop_consistent = false;
                    continue;
                }
                
                // Compute the total pop adjusted score for this multipath_alignment_t
                pop_adjusted_scores[i] = best_total_score[0] + best_total_score[1] + frag_score;
                
                // Save the pop scores without the base scores to the multipath alignments.
                // TODO: Should we be annotating unmapped reads with 0 pop scores when the other read in the pair is mapped?
                multipath_aln_pair.first.set_annotation("haplotype_score", best_pop_score[0]);
                multipath_aln_pair.second.set_annotation("haplotype_score", best_pop_score[1]);
                
                assert(!std::isnan(best_total_score[0]));
                assert(!std::isnan(best_total_score[1]));
                assert(!std::isnan(frag_score));
                assert(!std::isnan(pop_adjusted_scores[i]));
                
                // How much was extra over the score of the top-base-score alignment on each side?
                // This might be negative if e.g. that alignment looks terrible population-wise but we take it anyway.
                auto extra = pop_adjusted_scores[i] - alignment_score;
                
                // Record our extra score if it was a new minimum
                min_extra_score = min(extra, min_extra_score);
            }
        }
        
        // Decide which scores to use depending on whether we have pop adjusted scores we want to use
        if (include_population_component && all_multipaths_pop_consistent) {
            scores = move(pop_adjusted_scores);
        }
        
        // Pull the min frag or extra score out of the score so it will be nonnegative
        double zero_point = (include_population_component && all_multipaths_pop_consistent) ? min_extra_score : min_frag_score;
        for (auto& score : scores) {
            score -= zero_point;
        }
        
        if (include_population_component && all_multipaths_pop_consistent) {
            // Record that we used the population score
#ifdef debug_multipath_mapper
            cerr << "Haplotype consistency score is being used." << endl;
#endif
            for (auto& multipath_aln_pair : multipath_aln_pairs) {
                // We have to do it on each read in each pair.
                // TODO: Come up with a simpler way to dump annotations in based on what happens during mapping.
                multipath_aln_pair.first.set_annotation("haplotype_score_used", true);
                multipath_aln_pair.second.set_annotation("haplotype_score_used", true);
            }
        } else {
            // Clean up pop score annotations if present and remove scores from all the reads
#ifdef debug_multipath_mapper
            cerr << "Haplotype consistency score is not being used." << endl;
#endif
            for (auto& multipath_aln_pair : multipath_aln_pairs) {
                // We have to do it on each read in each pair.
                multipath_aln_pair.first.clear_annotation("haplotype_score_used");
                multipath_aln_pair.first.clear_annotation("haplotype_score");
                multipath_aln_pair.second.clear_annotation("haplotype_score_used");
                multipath_aln_pair.second.clear_annotation("haplotype_score");
            }
        }
        return scores;
    }
    
    vector<int32_t> MultipathMapper::compute_raw_mapping_qualities_from_scores(const vector<double>& scores, bool have_qualities,
                                                                               const vector<double>* multiplicities) const {
   
        auto aligner = get_aligner(have_qualities);
        
        vector<int32_t> raw_mapqs;
        
        if (scores.size() > 1 && max_alt_mappings > 1) {
            // we want MAPQs for all of the multi-mapped reads, so we need the exact algorithm
            raw_mapqs = aligner->compute_all_mapping_qualities(scores, multiplicities);
        }
        else {
            // we only need a MAPQ for the primary
            
            bool use_exact = (mapping_quality_method == Exact);
            if (!use_exact && scores.size() >= 2
                && (scores[1] > scores[0] ||
                    (mapping_quality_method == Adaptive && scores[1] < scores[0] - get_aligner()->mapping_quality_score_diff(max_mapping_quality)))) {
                use_exact = true;
            }
            
            raw_mapqs.push_back(aligner->compute_first_mapping_quality(scores, !use_exact, multiplicities));
        }
        
        // arbitrary scaling, seems to help performance
        for (auto& raw_mapq : raw_mapqs) {
            raw_mapq *= mapq_scaling_factor;
            
#ifdef debug_multipath_mapper
            cerr << "scores yield a raw MAPQ of " << raw_mapq << endl;
#endif
        }

        return raw_mapqs;

    }
    
    void MultipathMapper::sort_and_compute_mapping_quality(vector<multipath_alignment_t>& multipath_alns,
                                                           vector<size_t>* cluster_idxs,
                                                           vector<double>* multiplicities) const {
        if (cluster_idxs) {
            assert(cluster_idxs->size() == multipath_alns.size());
        }
        if (multiplicities) {
            assert(multiplicities->size() == multipath_alns.size());
        }
        
        if (multipath_alns.empty()) {
            return;
        }
        
        // get the log-likelihoods of each mapping
        vector<double> scores = mapping_likelihoods(multipath_alns);
        
        // find the order of the scores
        vector<size_t> order(multipath_alns.size(), 0);
        for (size_t i = 1; i < multipath_alns.size(); i++) {
            order[i] = i;
        }
        // Sort, shuffling based on the aligned sequence to break ties.
        LazyRNG rng([&]() {
            return make_shuffle_seed(multipath_alns.front());
        });
        sort_shuffling_ties(order.begin(), order.end(),
            [&](const size_t i, const size_t j) { return scores[i] > scores[j]; },
            rng);
        
        // translate the order to an index
        vector<size_t> index(multipath_alns.size());
        for (size_t i = 0; i < multipath_alns.size(); i++) {
            index[order[i]] = i;
        }
        
        // put the scores, multiplicities, clusters-of-origin, and alignments in order
        for (size_t i = 0; i < multipath_alns.size(); i++) {
            while (index[i] != i) {
                std::swap(scores[index[i]], scores[i]);
                std::swap(multipath_alns[index[i]], multipath_alns[i]);
                if (cluster_idxs) {
                    std::swap((*cluster_idxs)[index[i]], (*cluster_idxs)[i]);
                }
                if (multiplicities) {
                    std::swap((*multiplicities)[index[i]], (*multiplicities)[i]);
                }
                std::swap(index[index[i]], index[i]);
            }
        }
        
        // identify and remove duplicates
        size_t removed_so_far = 0;
        for (size_t i = 1; i < multipath_alns.size(); i++) {
            if (share_terminal_positions(multipath_alns.front(), multipath_alns[i])) {
                removed_so_far++;
            }
            else if (removed_so_far) {
                multipath_alns[i - removed_so_far] = move(multipath_alns[i]);
                scores[i - removed_so_far] = scores[i];
                if (cluster_idxs) {
                    (*cluster_idxs)[i - removed_so_far] = (*cluster_idxs)[i];
                }
                if (multiplicities) {
                    (*multiplicities)[i - removed_so_far] = (*multiplicities)[i];
                }
            }
        }
        if (removed_so_far) {
            multipath_alns.resize(multipath_alns.size() - removed_so_far);
            scores.resize(scores.size() - removed_so_far);
            if (cluster_idxs) {
                cluster_idxs->resize(cluster_idxs->size() - removed_so_far);
            }
            if (multiplicities) {
                multiplicities->resize(multiplicities->size() - removed_so_far);
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "scores obtained of multi-mappings:" << endl;
        for (size_t i = 0; i < scores.size(); i++) {
            Alignment aln;
            optimal_alignment(multipath_alns[i], aln);
            cerr << "\t" << scores[i] << " " << (aln.path().mapping_size() ? make_pos_t(aln.path().mapping(0).position()) : pos_t());
            if (multiplicities) {
                cerr << ", multiplicity " << multiplicities->at(i);
            }
            cerr << endl;
        }
#endif
        
        // Compute and set the mapping quality
        vector<int32_t> uncapped_mapqs = compute_raw_mapping_qualities_from_scores(scores, !multipath_alns.front().quality().empty(),
                                                                                   multiplicities);
        for (size_t i = 0; i < uncapped_mapqs.size(); ++i) {
            multipath_alns[i].set_mapping_quality(min<int32_t>(uncapped_mapqs[i], max_mapping_quality));
        }
        
        if (report_allelic_mapq) {
            // figure out what the mapping quality would be for the lowest-scoring combination of
            // alleles
            int32_t allelic_diff = optimal_alignment_score(multipath_alns.front()) - worst_alignment_score(multipath_alns.front());
            if (allelic_diff != 0) {
                scores[0] -= allelic_diff;
                vector<int32_t> uncapped_allelic_mapqs = compute_raw_mapping_qualities_from_scores(scores,
                                                                                                   !multipath_alns.front().quality().empty(),
                                                                                                   multiplicities);
                auto uncapped_allelic_mapq = uncapped_allelic_mapqs.front();
                int32_t allelic_mapq = min<int32_t>(uncapped_allelic_mapqs.front(), max_mapping_quality);
                if (allelic_mapq != multipath_alns.front().mapping_quality()) {
                    // other alleles do not place this read as confidently
                    multipath_alns.front().set_annotation("allelic_mapq", (double) allelic_mapq);
                }
                scores[0] += allelic_diff;
            }
        }
        
        if (report_group_mapq) {
            // TODO: this can include alignments that are later removed as being insigificant, but they shouldn't
            // affect the sum much at all
            size_t num_reporting = min(multipath_alns.size(), max_alt_mappings);
            vector<size_t> reporting_idxs(num_reporting, 0);
            for (size_t i = 1; i < num_reporting; ++i) {
                reporting_idxs[i] = i;
            }
            double raw_mapq = get_aligner(!multipath_alns.front().quality().empty())->compute_group_mapping_quality(scores, reporting_idxs,
                                                                                                                    multiplicities);
            // TODO: for some reason set_annotation will accept a double but not an int
            double group_mapq = min<double>(max_mapping_quality, mapq_scaling_factor * raw_mapq);
            
            for (size_t i = 0; i < num_reporting; ++i) {
                multipath_alns[i].set_annotation("group_mapq", group_mapq);
            }
        }
    }
    
    // TODO: pretty duplicative with the unpaired version
    void MultipathMapper::sort_and_compute_mapping_quality(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs,
                                                           vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                           vector<pair<size_t, size_t>>* duplicate_pairs_out,
                                                           vector<double>* multiplicities) const {
        
#ifdef debug_multipath_mapper
        cerr << "Sorting and computing mapping qualities for paired reads" << endl;
#endif
        
        assert(multipath_aln_pairs.size() == cluster_pairs.size());
        if (multiplicities) {
            assert(multipath_aln_pairs.size() == multiplicities->size());
        }
        
        if (multipath_aln_pairs.empty()) {
            return;
        }
        
        // get the log-likelihoods of each mapping
        vector<double> scores = pair_mapping_likelihoods(multipath_aln_pairs, cluster_pairs);
        
        // find the order of the scores
        vector<size_t> order(multipath_aln_pairs.size(), 0);
        for (size_t i = 1; i < multipath_aln_pairs.size(); i++) {
            order[i] = i;
        }
        // Sort, shuffling based on the aligned sequence to break ties.
        LazyRNG rng([&]() {
            return make_shuffle_seed(multipath_aln_pairs.front());
        });
        sort_shuffling_ties(order.begin(), order.end(),
            [&](const size_t i, const size_t j) {
                return scores[i] > scores[j];
            },
            rng);
        
        // translate the order to an index
        vector<size_t> index(multipath_aln_pairs.size());
        for (size_t i = 0; i < multipath_aln_pairs.size(); i++) {
            index[order[i]] = i;
        }
        
        // put the scores, distances, and alignments in order
        for (size_t i = 0; i < multipath_aln_pairs.size(); i++) {
            while (index[i] != i) {
                std::swap(scores[index[i]], scores[i]);
                std::swap(cluster_pairs[index[i]], cluster_pairs[i]);
                std::swap(multipath_aln_pairs[index[i]], multipath_aln_pairs[i]);
                if (multiplicities) {
                    std::swap((*multiplicities)[index[i]], (*multiplicities)[i]);
                }
                std::swap(index[index[i]], index[i]);
                
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "scores and distances obtained of multi-mappings:" << endl;
        for (int i = 0; i < multipath_aln_pairs.size(); i++) {
            Alignment aln1, aln2;
            optimal_alignment(multipath_aln_pairs[i].first, aln1);
            optimal_alignment(multipath_aln_pairs[i].second, aln2);
            int64_t id1 = 0, id2 = 0;
            bool rev1 = false, rev2 = false;
            if (aln1.path().mapping_size()) {
                const auto& pos = aln1.path().mapping(0).position();
                id1 = pos.node_id();
                rev1 = pos.is_reverse();
            }
            else {
                const auto& pos = multipath_aln_pairs[i].first.subpath(0).path().mapping(0).position();
                id1 = pos.node_id();
                rev1 = pos.is_reverse();
            }
            if (aln2.path().mapping_size()) {
                const auto& pos = aln2.path().mapping(0).position();
                id2 = pos.node_id();
                rev2 = pos.is_reverse();
            }
            else {
                const auto& pos = multipath_aln_pairs[i].second.subpath(0).path().mapping(0).position();
                id2 = pos.node_id();
                rev2 = pos.is_reverse();
            }
        
            cerr << "\tpos:" << id1 << (rev1 ? "-" : "+") << "(" << aln1.score() << ") - " << id2 << (rev2 ? "-" : "+") << "(" << aln2.score() << ")"
                << " align:" << optimal_alignment_score(multipath_aln_pairs[i].first) + optimal_alignment_score(multipath_aln_pairs[i].second)
                << ", length: " << cluster_pairs[i].second;
            cerr << ", combined: " << scores[i];
            if (multiplicities) {
                cerr << ", multiplicity: " << multiplicities->at(i);
            }
            cerr << endl;
        }
#endif
        
        // Compute the raw mapping quality
        vector<int32_t> uncapped_mapqs = compute_raw_mapping_qualities_from_scores(scores,
                                                                                   !multipath_aln_pairs.front().first.quality().empty() &&
                                                                                   !multipath_aln_pairs.front().second.quality().empty(),
                                                                                   multiplicities);
        // Limit it to the max.
        for (size_t i = 0; i < uncapped_mapqs.size(); ++i) {
            int32_t mapq = min<int32_t>(uncapped_mapqs[i], max_mapping_quality);
            multipath_aln_pairs[i].first.set_mapping_quality(mapq);
            multipath_aln_pairs[i].second.set_mapping_quality(mapq);
        }
        
        
        int32_t allelic_diff_1 = 0, allelic_diff_2 = 0;
        if (report_allelic_mapq) {
            // figure out what the mapping quality would be for the lowest-scoring combination of  alleles
            allelic_diff_1 = (optimal_alignment_score(multipath_aln_pairs.front().first)
                              - worst_alignment_score(multipath_aln_pairs.front().first));
            allelic_diff_2 = (optimal_alignment_score(multipath_aln_pairs.front().second)
                              - worst_alignment_score(multipath_aln_pairs.front().second));
            if (allelic_diff_1 != 0 || allelic_diff_2 != 0) {
                scores[0] -= allelic_diff_1 + allelic_diff_2;
                vector<int32_t> uncapped_allelic_mapqs = compute_raw_mapping_qualities_from_scores(scores,
                                                                                                   !multipath_aln_pairs.front().first.quality().empty() &&
                                                                                                   !multipath_aln_pairs.front().second.quality().empty(),
                                                                                                   multiplicities);
                auto uncapped_allelic_mapq = uncapped_allelic_mapqs.front();
                int32_t allelic_mapq = min<int32_t>(uncapped_allelic_mapq, max_mapping_quality);
                if (allelic_mapq != multipath_aln_pairs.front().first.mapping_quality() ||
                    allelic_mapq != multipath_aln_pairs.front().second.mapping_quality()) {
                    // other alleles might not place this read as confidently
                    multipath_aln_pairs.front().first.set_annotation("allelic_mapq", (double) allelic_mapq);
                    multipath_aln_pairs.front().second.set_annotation("allelic_mapq", (double) allelic_mapq);
                }
                scores[0] += allelic_diff_1 + allelic_diff_2;
            }
        }
        
        
        if (multipath_aln_pairs.size() > 1) {
            // TODO: it would be nice to also look for duplicates with other pairs, but i don't love
            // the quadratic work that this would require...
            
            // find the duplicates of the optimal pair (initially mark with only the pair itself)
            vector<size_t> duplicates_1(1, 0);
            vector<size_t> duplicates_2(1, 0);
            vector<size_t> to_remove;
            for (size_t i = 1; i < multipath_aln_pairs.size(); i++) {
                bool duplicate_1 = share_terminal_positions(multipath_aln_pairs[0].first, multipath_aln_pairs[i].first);
                bool duplicate_2 = share_terminal_positions(multipath_aln_pairs[0].second, multipath_aln_pairs[i].second);
                if (duplicate_1 && duplicate_2) {
#ifdef debug_multipath_mapper
                    cerr << "found double end duplication at index " << i << endl;
#endif
                    // this pair is a complete duplication (not just one end) we want it gone
                    to_remove.push_back(i);
                    if (duplicate_pairs_out) {
                        duplicate_pairs_out->push_back(cluster_pairs[i].first);
                    }
                }
                else if (duplicate_1) {
#ifdef debug_multipath_mapper
                    cerr << "found left end duplication at index " << i << endl;
#endif
                    duplicates_1.push_back(i);
                }
                else if (duplicate_2) {
#ifdef debug_multipath_mapper
                    cerr << "found right end duplication at index " << i << endl;
#endif
                    duplicates_2.push_back(i);
                }
            }
            
            if (!to_remove.empty()) {
                
                // remove the full duplicates from all relevant vectors
                for (size_t i = 1, removed_so_far = 0; i < multipath_aln_pairs.size(); i++) {
                    if (removed_so_far < to_remove.size() ? i == to_remove[removed_so_far] : false) {
                        removed_so_far++;
                    }
                    else if (removed_so_far > 0) {
                        // move these items into their new position
                        multipath_aln_pairs[i - removed_so_far] = move(multipath_aln_pairs[i]);
                        scores[i - removed_so_far] = scores[i];
                        cluster_pairs[i - removed_so_far] = move(cluster_pairs[i]);
                        if (multiplicities) {
                            (*multiplicities)[i - removed_so_far] = (*multiplicities)[i];
                        }
                    }
                }
                
                // remove the end positions that are now empty
                multipath_aln_pairs.resize(multipath_aln_pairs.size() - to_remove.size());
                scores.resize(scores.size() - to_remove.size());
                cluster_pairs.resize(cluster_pairs.size() - to_remove.size());
                if (multiplicities) {
                    multiplicities->resize(multiplicities->size() - to_remove.size());
                }
                
                // update the indexes of the marked single-end duplicates
                for (size_t i = 0, removed_so_far = 0; i < duplicates_1.size(); i++) {
                    while (removed_so_far < to_remove.size() ? to_remove[removed_so_far] < duplicates_1[i] : false) {
                        removed_so_far++;
                    }
                    duplicates_1[i] -= removed_so_far;
                }
                
                for (size_t i = 0, removed_so_far = 0; i < duplicates_2.size(); i++) {
                    while (removed_so_far < to_remove.size() ? to_remove[removed_so_far] < duplicates_2[i] : false) {
                        removed_so_far++;
                    }
                    duplicates_2[i] -= removed_so_far;
                }
            }
            
            // did we find any duplicates with the optimal pair?
            if (duplicates_1.size() > 1 || duplicates_2.size() > 1 || !to_remove.empty()) {
                // compute the mapping quality of the whole group of duplicates for each end
                auto aligner = get_aligner(!multipath_aln_pairs.front().first.quality().empty() &&
                                           !multipath_aln_pairs.front().second.quality().empty());
                
                int32_t raw_mapq_1 = aligner->compute_group_mapping_quality(scores, duplicates_1, multiplicities);
                int32_t raw_mapq_2 = aligner->compute_group_mapping_quality(scores, duplicates_2, multiplicities);
                
#ifdef debug_multipath_mapper
                cerr << "deduplicated raw MAPQs are " << raw_mapq_1 << " and " << raw_mapq_2 << endl;
#endif
                
                // arbitrary scaling, seems to help performance
                int32_t mapq_1 = min<int32_t>(raw_mapq_1 * mapq_scaling_factor, max_mapping_quality);
                int32_t mapq_2 = min<int32_t>(raw_mapq_2 * mapq_scaling_factor, max_mapping_quality);
                
#ifdef debug_multipath_mapper
                cerr << "processed MAPQs are " << mapq_1 << " and " << mapq_2 << endl;
#endif
                
                multipath_aln_pairs.front().first.set_mapping_quality(mapq_1);
                multipath_aln_pairs.front().second.set_mapping_quality(mapq_2);
                
                if (report_allelic_mapq && (allelic_diff_1 != 0 || allelic_diff_2)) {
                    for (auto i : duplicates_1) {
                        scores[i] -= allelic_diff_1;
                    }
                    for (auto i : duplicates_2) {
                        scores[i] -= allelic_diff_2;
                    }
                    
                    int32_t raw_allelic_mapq_1 = aligner->compute_group_mapping_quality(scores, duplicates_1, multiplicities);
                    int32_t raw_allelic_mapq_2 = aligner->compute_group_mapping_quality(scores, duplicates_2, multiplicities);
                    
                    int32_t allelic_mapq_1 = min<int32_t>(raw_mapq_1 * mapq_scaling_factor, max_mapping_quality);
                    int32_t allelic_mapq_2 = min<int32_t>(raw_mapq_2 * mapq_scaling_factor, max_mapping_quality);
                    
                    if (allelic_mapq_1 != mapq_1) {
                        // other alleles might not place this read as confidently
                        multipath_aln_pairs.front().first.set_annotation("allelic_mapq", (double) allelic_mapq_1);
                    }
                    if (allelic_mapq_2 != mapq_2) {
                        // other alleles might not place this read as confidently
                        multipath_aln_pairs.front().second.set_annotation("allelic_mapq", (double) allelic_mapq_2);
                    }
                    
                    for (auto i : duplicates_1) {
                        scores[i] += allelic_diff_1;
                    }
                    for (auto i : duplicates_2) {
                        scores[i] += allelic_diff_2;
                    }
                }
            }
        }
        
        if (report_group_mapq) {
            // TODO: this can include alignments that are later removed as being insigificant, but they shouldn't
            // affect the sum much at all
            size_t num_reporting = min(multipath_aln_pairs.size(), max_alt_mappings);
            vector<size_t> reporting_idxs(num_reporting, 0);
            for (size_t i = 1; i < num_reporting; ++i) {
                reporting_idxs[i] = i;
            }
            auto aligner = get_aligner(!multipath_aln_pairs.front().first.quality().empty() &&
                                       !multipath_aln_pairs.front().second.quality().empty());
            double raw_mapq = aligner->compute_group_mapping_quality(scores, reporting_idxs,
                                                                     multiplicities);
            
            // TODO: for some reason set_annotation will accept a double but not an int
            double group_mapq = min<double>(max_mapping_quality, mapq_scaling_factor * raw_mapq);
            for (size_t i = 0; i < num_reporting; ++i) {
                multipath_aln_pairs[i].first.set_annotation("group_mapq", group_mapq);
                multipath_aln_pairs[i].second.set_annotation("group_mapq", group_mapq);
            }
        }
    }
            
    double MultipathMapper::fragment_length_log_likelihood(int64_t length) const {
        double dev = length - fragment_length_distr.mean();
        return -dev * dev / (2.0 * fragment_length_distr.std_dev() * fragment_length_distr.std_dev());
    }
    
    void MultipathMapper::set_automatic_min_clustering_length(double random_mem_probability) {
        min_clustering_mem_length = max<int>(log(1.0 - pow(random_mem_probability, 1.0 / total_seq_length)) / log(0.25), 1);
    }

    void MultipathMapper::set_min_softclip_length_for_splice(size_t length) {
        min_softclip_length_for_splice = length;
        
        // find the lowest score that could correspond to a high quality match of this length
        // TODO: kinda ugly, but whatever
        string dummy_a(length, 'A');
        string dummy_c(length, 'C');
        string dummy_g(length, 'G');
        string dummy_t(length, 'T');
        int32_t score_a = get_regular_aligner()->score_exact_match(dummy_a);
        int32_t score_c = get_regular_aligner()->score_exact_match(dummy_c);
        int32_t score_g = get_regular_aligner()->score_exact_match(dummy_g);
        int32_t score_t = get_regular_aligner()->score_exact_match(dummy_t);
        int32_t lowest_score = min(score_a, min(score_c, min(score_g, score_t)));
        
        // add in the full length bonus, this is the criterion we will actually check against
        string dummy_qual(length, char(40));
        min_softclipped_score_for_splice = lowest_score + get_regular_aligner()->score_full_length_bonus(false, dummy_a.begin(),
                                                                                                         dummy_a.end(),
                                                                                                         dummy_qual.begin());
    }

    void MultipathMapper::set_log_odds_against_splice(double log_odds) {
        no_splice_natural_log_odds = log_odds;
        no_splice_log_odds = round(log_odds / get_regular_aligner()->log_base);
    }

    void MultipathMapper::set_intron_length_distribution(const vector<double>& intron_mixture_weights,
                                                         const vector<pair<double, double>>& intron_component_params) {
        splice_stats.update_intron_length_distribution(intron_mixture_weights, intron_component_params, *get_aligner());
    }

    void MultipathMapper::set_max_merge_supression_length() {
        max_tail_merge_supress_length = ceil(double(get_regular_aligner()->match) / double(get_regular_aligner()->mismatch));
    }

    void MultipathMapper::set_read_1_adapter(const string& adapter) {
        // TODO: magic number
        string trimmed_adapter = adapter;
        if (trimmed_adapter.size() > 12) {
            trimmed_adapter.resize(12);
        }
        read_1_adapter = trimmed_adapter;
        read_1_adapter_lps = make_prefix_suffix_table(read_1_adapter.c_str(), read_1_adapter.size());
    }

    void MultipathMapper::set_read_2_adapter(const string& adapter) {
        // TODO: magic number
        string trimmed_adapter = adapter;
        if (trimmed_adapter.size() > 12) {
            trimmed_adapter.resize(12);
        }
        read_2_adapter = reverse_complement(trimmed_adapter);
        read_2_adapter_lps = make_prefix_suffix_table(read_2_adapter.c_str(), read_2_adapter.size());
    }

    // make the memos live in this .o file
    thread_local unordered_map<pair<double, size_t>, haploMath::RRMemo> MultipathMapper::rr_memos;
    
    haploMath::RRMemo& MultipathMapper::get_rr_memo(double recombination_penalty, size_t population_size) const {
        auto iter = rr_memos.find(make_pair(recombination_penalty, population_size));
        if (iter != rr_memos.end()) {
            return iter->second;
        }
        else {
            rr_memos.insert(make_pair(make_pair(recombination_penalty, population_size),
                                      haploMath::RRMemo(recombination_penalty, population_size)));
            return rr_memos.at(make_pair(recombination_penalty, population_size));
        }
    }
}





//
//  multipath_mapper.cpp
//  
//
//

//#define debug_multipath_mapper
//#define debug_multipath_mapper_alignment
//#define debug_validate_multipath_alignments
//#define debug_report_startup_training
//#define debug_pretty_print_alignments

#include "multipath_mapper.hpp"

// include this here to avoid a circular dependency
#include "multipath_alignment_graph.hpp"

namespace vg {
    
    //size_t MultipathMapper::PRUNE_COUNTER = 0;
    //size_t MultipathMapper::SUBGRAPH_TOTAL = 0;
    //size_t MultipathMapper::SECONDARY_RESCUE_COUNT = 0;
    //size_t MultipathMapper::SECONDARY_RESCUE_ATTEMPT = 0;
    //size_t MultipathMapper::SECONDARY_RESCUE_TOTAL = 0;
    
    MultipathMapper::MultipathMapper(PathPositionHandleGraph* graph, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_array,
                                     haplo::ScoreProvider* haplo_score_provider, SnarlManager* snarl_manager,
                                     MinimumDistanceIndex* distance_index) :
        BaseMapper(graph, gcsa_index, lcp_array, haplo_score_provider),
        snarl_manager(snarl_manager),
        distance_index(distance_index),
        path_component_index(graph)
    {
        // nothing to do
    }
    
    MultipathMapper::~MultipathMapper() {
        
    }
    
    void MultipathMapper::multipath_map(const Alignment& alignment,
                                        vector<multipath_alignment_t>& multipath_alns_out) {
        multipath_map_internal(alignment, mapping_quality_method, multipath_alns_out);
    }
    
    void MultipathMapper::multipath_map_internal(const Alignment& alignment,
                                                 MappingQualityMethod mapq_method,
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
                                                                                                      &path_component_index));
        }
        
        vector<memcluster_t> clusters = get_clusters(alignment, mems, &(*distance_measurer), fanouts.get());
        
#ifdef debug_multipath_mapper
        cerr << "obtained clusters:" << endl;
        for (int i = 0; i < clusters.size(); i++) {
            cerr << "\tcluster " << i << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : clusters[i]) {
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
        vector<size_t> cluster_idxs = range_vector(cluster_graphs.size());
        align_to_cluster_graphs(alignment, mapq_method, cluster_graphs, multipath_alns_out, num_mapping_attempts,
                                fanouts.get(), &cluster_idxs);
        
        if (multipath_alns_out.empty()) {
            // add a null alignment so we know it wasn't mapped
            multipath_alns_out.emplace_back();
            cluster_graphs.emplace_back();
            cluster_idxs.push_back(0);
            to_multipath_alignment(alignment, multipath_alns_out.back());
            
            // in case we're realigning GAMs that have paths already
            multipath_alns_out.back().clear_subpath();
            multipath_alns_out.back().clear_start();
        }
        
        if (likely_mismapping(multipath_alns_out.front())) {
            // we can't distinguish this alignment from the longest MEM of a random sequence
#ifdef debug_multipath_mapper
            cerr << "mapping is not distinguishable from a random sequence, snapping MAPQ to 0" << endl;
#endif
            
            multipath_alns_out.front().set_mapping_quality(0);
        }
        else {
            // account for the possiblity that we missed the correct cluster due to sub-sampling MEM hits
            // within the cluster
            cap_mapping_quality_by_hit_sampling_probability(multipath_alns_out, cluster_idxs, cluster_graphs);
        }
        
        // if we computed extra alignments to get a mapping quality, remove them
        if (multipath_alns_out.size() > max_alt_mappings) {
            multipath_alns_out.resize(max_alt_mappings);
        }
        
        if (simplify_topologies) {
            for (multipath_alignment_t& multipath_aln : multipath_alns_out) {
                merge_non_branching_subpaths(multipath_aln);
            }
        }
        
        if (strip_bonuses) {
            for (multipath_alignment_t& multipath_aln : multipath_alns_out) {
                strip_full_length_bonuses(multipath_aln);
            }
        }
        
        // clean up the cluster graphs
        for (auto cluster_graph : cluster_graphs) {
            delete get<0>(cluster_graph);
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
                                   min_median_mem_coverage_for_split, 0.75, fanouts);;
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
            auto& mem_cluster = *cluster_mems_1[i];
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
            auto& mem_cluster = *cluster_mems_2[i];
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
                                                  MappingQualityMethod mapq_method,
                                                  vector<clustergraph_t>& cluster_graphs,
                                                  vector<multipath_alignment_t>& multipath_alns_out,
                                                  size_t num_mapping_attempts,
                                                  const match_fanouts_t* fanouts,
                                                  vector<size_t>* cluster_idxs) {
        
        
#ifdef debug_multipath_mapper
        cerr << "aligning to (up to) " << cluster_graphs.size() << " subgraphs..." << endl;
#endif
      
        // we may need to compute an extra mapping above the one we'll report if we're computing mapping quality
        size_t num_mappings_to_compute = mapq_method != None ? max(num_mapping_attempts, (size_t) 2) : num_mapping_attempts;
        
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
            cerr << "performing alignment to subgraph" << endl;
#endif
            
            multipath_alns_out.emplace_back();
            multipath_align(alignment, get<0>(cluster_graph), get<1>(cluster_graph), multipath_alns_out.back(),
                            fanouts);
            
            num_mappings++;
        }
        
        // if we didn't end up performing all of the cluster alignments, re-sync the list
        if (cluster_idxs) {
            if (cluster_idxs->size() > multipath_alns_out.size()) {
                cluster_idxs->resize(multipath_alns_out.size());
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "splitting multicomponent alignments..." << endl;
#endif
        
        // split up any alignments that ended up being disconnected
        split_multicomponent_alignments(multipath_alns_out, cluster_idxs);
        
#ifdef debug_multipath_mapper
        cerr << "topologically ordering " << multipath_alns_out.size() << " multipath alignments" << endl;
#endif
        for (multipath_alignment_t& multipath_aln : multipath_alns_out) {
            topologically_order_subpaths(multipath_aln);
        }
        
#ifdef debug_multipath_mapper
        cerr << "computing mapping quality and sorting mappings" << endl;
#endif
        sort_and_compute_mapping_quality(multipath_alns_out, mapq_method, cluster_idxs);
        
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
    
    void MultipathMapper::attempt_unpaired_multipath_map_of_pair(const Alignment& alignment1, const Alignment& alignment2,
                                                                 vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                                 vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer) {
        
        // compute single ended mappings, and make sure we also compute mapping qualities to assess
        // mapping ambiguity
        vector<multipath_alignment_t> multipath_alns_1, multipath_alns_2;
        multipath_map_internal(alignment1, mapping_quality_method == None ? Approx : mapping_quality_method,
                               multipath_alns_1);
        multipath_map_internal(alignment2, mapping_quality_method == None ? Approx : mapping_quality_method,
                               multipath_alns_2);
        
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
            if (multipath_aln_1.mapping_quality() >= min(max_mapping_quality, 45)
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
                cerr << "warning:[vg mpmap] Mapped " << ambiguous_pair_buffer.size() + fragment_length_distr.curr_sample_size() << " read pairs as unpaired reads to learn fragment length distribution, but only obtained " << fragment_length_distr.curr_sample_size() << " unambiguous, consistently mapped pairs. Often this indicates data issues, such as reads that are pre-sorted with unmappable reads at the front, or reads that are not actually paired." << endl;
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
    }
    
    bool MultipathMapper::attempt_rescue(const multipath_alignment_t& multipath_aln, const Alignment& other_aln,
                                         bool rescue_forward, multipath_alignment_t& rescue_multipath_aln) {
        
#ifdef debug_multipath_mapper
        cerr << "attemping pair rescue in " << (rescue_forward ? "forward" : "backward") << " direction from " << debug_string(multipath_aln) << endl;
#endif
        bdsg::HashGraph rescue_graph;
        extract_rescue_graph(multipath_aln, other_aln, rescue_forward, &rescue_graph);
        
        if (rescue_graph.get_node_count() == 0) {
            return false;
        }
        
#ifdef debug_multipath_mapper
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
        DagifiedGraph* dagified = nullptr;
        
        ExpandingOverlayGraph* align_dag = nullptr;
        if (algorithms::is_directed_acyclic(&align_digraph)) {
            align_dag = &undagified;
        }
        else {
#ifdef debug_multipath_mapper_alignment
            cerr << "graph contains directed cycles, performing dagification" << endl;
#endif
            dagified = new DagifiedGraph(&align_digraph, target_length);
            align_dag = dagified;
        }
        
        // put local alignment here
        Alignment aln = other_aln;
        // in case we're realigning a GAM, get rid of the path
        aln.clear_path();
        
        aligner->align(aln, *align_dag, true);
        
        // get the IDs back into the space of the reference graph
        function<pair<id_t, bool>(id_t)> translator = [&](const id_t& node_id) {
            handle_t original = align_digraph.get_underlying_handle(align_dag->get_underlying_handle(align_dag->get_handle(node_id)));
            return make_pair(rescue_graph.get_id(original), rescue_graph.get_is_reverse(original));
        };
        translate_oriented_node_ids(*aln.mutable_path(), translator);
        
#ifdef debug_multipath_mapper
        cerr << "resecued direct alignment is" << endl;
        cerr << pb2json(aln) << endl;
#endif
        
        if (num_alt_alns > 1 && snarl_manager != nullptr) {
            // make an interesting multipath alignment by realigning the single path alignment inside snarls
            make_nontrivial_multipath_alignment(aln, *align_dag, translator, rescue_multipath_aln);
            
        }
        else {
            // just convert the single alignment into a trivial multipath alignment
            to_multipath_alignment(aln, rescue_multipath_aln);
        }
        
        identify_start_subpaths(rescue_multipath_aln);
        
        vector<double> score(1, aln.score());
        int32_t raw_mapq = aligner->compute_mapping_quality(score, mapping_quality_method == None || mapping_quality_method == Approx);
        int32_t adjusted_mapq = min(raw_mapq, min(max_mapping_quality, multipath_aln.mapping_quality()));
        rescue_multipath_aln.set_mapping_quality(adjusted_mapq);
        
#ifdef debug_multipath_mapper
        cerr << "converted multipath alignment is" << endl;
        cerr << debug_string(rescue_multipath_aln) << endl;
        cerr << "rescued alignment has effective match length " << pseudo_length(rescue_multipath_aln) << ", which gives p-value " << random_match_p_value(pseudo_length(rescue_multipath_aln), rescue_multipath_aln.sequence().size()) << endl;
#endif

        
        if (raw_mapq < min(25, max_mapping_quality)) {
#ifdef debug_multipath_mapper
            cerr << "rescue fails because raw_mapq " << raw_mapq << " < " << min(25, max_mapping_quality) << endl;
#endif
            return false;
        }
        
        auto p_val = random_match_p_value(pseudo_length(rescue_multipath_aln), rescue_multipath_aln.sequence().size());
        
        if (p_val >= max_mapping_p_value * 0.1) {
#ifdef debug_multipath_mapper
            cerr << "rescue fails because p value " << p_val << " >= " <<  (max_mapping_p_value * 0.1) << endl;
#endif
            return false;
        }
        
        // clean up (safe on nullptr)
        delete dagified;
        
        return true;
    }

    void MultipathMapper::extract_rescue_graph(const multipath_alignment_t& multipath_aln, const Alignment& other_aln,
                                               bool rescue_forward, MutableHandleGraph* rescue_graph) const {
        
        // get the position to jump from and the distance to jump
        Alignment opt_anchoring_aln;
        optimal_alignment(multipath_aln, opt_anchoring_aln);
                
        if (get_rescue_graph_from_paths || !distance_index) {
            // we're either not using the distance index or we don't have one
            
            pos_t pos_from = rescue_forward ? initial_position(opt_anchoring_aln.path()) : final_position(opt_anchoring_aln.path());
            int64_t jump_dist = rescue_forward ? fragment_length_distr.mean() : -fragment_length_distr.mean();
            
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
                search_dist_bwd = size_t(round(rescue_graph_std_devs * fragment_length_distr.std_dev())) + other_aln.sequence().size();
                search_dist_fwd = rescue_graph_std_devs * fragment_length_distr.std_dev();
            }
            else {
                search_dist_bwd = rescue_graph_std_devs * fragment_length_distr.std_dev();
                search_dist_fwd = size_t(round(rescue_graph_std_devs * fragment_length_distr.std_dev())) + other_aln.sequence().size();
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
            int64_t min_distance = max(0.0, fragment_length_distr.mean() - other_aln.sequence().size()
                                       - rescue_graph_std_devs * fragment_length_distr.std_dev());
            int64_t max_distance = fragment_length_distr.mean() + rescue_graph_std_devs * fragment_length_distr.std_dev();
            distance_index->subgraph_in_range(opt_anchoring_aln.path(), xindex, min_distance, max_distance,
                                              subgraph_nodes_to_add, rescue_forward);
            
            
#ifdef debug_multipath_mapper
            cerr << "got distance based rescue graph:" << endl;
            for (auto nid : subgraph_nodes_to_add) {
                cerr << "\t" << nid << endl;
            }
#endif
            
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
#ifdef debug_multipath_mapper
                cerr << "initializing extraction search at node " << node_id << endl;
#endif
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
#ifdef debug_multipath_mapper
                                cerr << "reached " << xindex->get_id(neighbor) << " from searching" << endl;
#endif
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
#ifdef debug_multipath_mapper
                                    cerr << "create edge " << rescue_graph->get_id(sub_neighbor) << (rescue_graph->get_is_reverse(sub_neighbor) ? "-" : "+") << " -> " << rescue_graph->get_id(sub_handle) << (rescue_graph->get_is_reverse(sub_handle) ? "-" : "+")  << endl;
#endif
                                    rescue_graph->create_edge(sub_neighbor, sub_handle);
                                }
                                else {
#ifdef debug_multipath_mapper
                                    cerr << "create edge " << rescue_graph->get_id(sub_handle) << (rescue_graph->get_is_reverse(sub_handle) ? "-" : "+") << " -> " << rescue_graph->get_id(sub_neighbor) << (rescue_graph->get_is_reverse(sub_neighbor) ? "-" : "+")  << endl;
#endif
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
    
    bool MultipathMapper::likely_mismapping(const multipath_alignment_t& multipath_aln) {
        if (!suppress_mismapping_detection) {
            // empirically, we get better results by scaling the pseudo-length down, I have no good explanation for this probabilistically
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
    
    size_t MultipathMapper::pseudo_length(const multipath_alignment_t& multipath_aln) const {
        Alignment alignment;
        optimal_alignment(multipath_aln, alignment);
        const Path& path = alignment.path();
        
        int64_t net_matches = 0;
        for (size_t i = 0; i < path.mapping_size(); i++) {
            const Mapping& mapping = path.mapping(i);
            for (size_t j = 0; j < mapping.edit_size(); j++) {
                const Edit& edit = mapping.edit(j);
                
                // skip soft clips
                if (((i == 0 && j == 0) || (i == path.mapping_size() - 1 && j == mapping.edit_size() - 1))
                    && edit.from_length() == 0 && edit.to_length() > 0) {
                    continue;
                }
                
                // add matches and subtract mismatches/indels
                if (edit.from_length() == edit.to_length() && edit.sequence().empty()) {
                    net_matches += edit.from_length();
                }
                else {
                    net_matches -= max(edit.from_length(), edit.to_length());
                }
            }
        }
        
        return max<int64_t>(0, net_matches);
    }
    
    // make the memo live in this .o file
    thread_local unordered_map<pair<size_t, size_t>, double> MultipathMapper::p_value_memo;
    
    double MultipathMapper::random_match_p_value(size_t match_length, size_t read_length) {
        // memoized to avoid transcendental functions (at least in cases where read lengths don't vary too much)
        auto iter = p_value_memo.find(make_pair(match_length, read_length));
        if (iter != p_value_memo.end()) {
            return iter->second;
        }
        else {
            double p_value;
            if (use_weibull_calibration) {
                double scale = exp(weibull_scale_intercept + weibull_scale_slope * log(read_length));
                double shape = exp(weibull_shape_intercept + weibull_shape_slope * log(read_length));
                double offset = exp(weibull_offset_intercept + weibull_offset_slope * log(read_length));
                p_value = 1.0 - weibull_cdf(match_length, scale, shape, offset);
            }
            else {
                double rate = max_exponential_rate_intercept + max_exponential_rate_slope * read_length;
                double shape = exp(max_exponential_shape_intercept + max_exponential_shape_slope * read_length);
                p_value = 1.0 - max_exponential_cdf(match_length, rate, shape);
            }
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
        vector<double> log_mle_weibull_scales;
        vector<double> log_mle_weibull_shapes;
        vector<double> log_mle_weibull_offsets;
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
            
            // alternatively model lengths with a weibull distribution that has an offset
            auto weibull_params = fit_offset_weibull(pseudo_lengths);
            log_mle_weibull_scales.push_back(log(get<0>(weibull_params)));
            log_mle_weibull_shapes.push_back(log(get<1>(weibull_params)));
            log_mle_weibull_offsets.push_back(log(get<2>(weibull_params)));
            
            auto max_exp_params = fit_max_exponential(pseudo_lengths);
            mle_max_exponential_rates.push_back(max_exp_params.first);
            log_mle_max_exponential_shapes.push_back(log(max_exp_params.second));
                        
#ifdef debug_report_startup_training
            cerr << "trained parameters for length " << simulated_read_length << ": " << endl;
            cerr << "\tweibull scale: " << get<0>(weibull_params) << endl;
            cerr << "\tweibull shape: " << get<1>(weibull_params) << endl;
            cerr << "\tweibull offset: " << get<2>(weibull_params) << endl;
            cerr << "\tmax exp rate: " << max_exp_params.first << endl;
            cerr << "\tmax exp shape: " << max_exp_params.second << endl;
#endif
        }
        
        // make a design matrix for a log regression and a linear regression
        vector<vector<double>> X(simulated_read_lengths.size());
        vector<vector<double>> X_log(simulated_read_lengths.size());
        for (size_t i = 0; i < X.size(); ++i) {
            X[i].resize(2, 1.0);
            X_log[i].resize(2, 1.0);
            X[i][1] = simulated_read_lengths[i];
            X_log[i][1] = log(simulated_read_lengths[i]);
        }
        
        auto weibull_scale_coefs = regress(X, log_mle_weibull_scales);
        auto weibull_shape_coefs = regress(X, log_mle_weibull_shapes);
        auto weibull_offset_coefs = regress(X, log_mle_weibull_offsets);
        
        weibull_scale_intercept = weibull_scale_coefs[0];
        weibull_scale_slope = weibull_scale_coefs[1];
        weibull_shape_intercept = weibull_shape_coefs[0];
        weibull_shape_slope = weibull_shape_coefs[1];
        weibull_offset_intercept = weibull_offset_coefs[0];
        weibull_offset_slope = weibull_offset_coefs[1];
        
        auto max_exp_rate_coefs = regress(X, mle_max_exponential_rates);
        auto max_exp_shape_coefs = regress(X, log_mle_max_exponential_shapes);
        
        max_exponential_rate_intercept = max_exp_rate_coefs[0];
        max_exponential_rate_slope = max_exp_rate_coefs[1];
        max_exponential_shape_intercept = max_exp_shape_coefs[0];
        max_exponential_shape_slope = max_exp_shape_coefs[1];
        
#ifdef debug_report_startup_training
        cerr << "final regression parameters:" << endl;
        cerr << "\tweibull scale = exp(" << weibull_scale_intercept << " + " << weibull_scale_slope << " * log L)" << endl;
        cerr << "\tweibull shape = exp(" << weibull_shape_intercept << " + " << weibull_shape_slope << " * log L)" << endl;
        cerr << "\tweibull offset = exp(" << weibull_offset_intercept << " + " << weibull_offset_slope << " * log L)" << endl;
        cerr << "\tmax exp rate = " << max_exponential_rate_intercept << " + " << max_exponential_rate_slope << " * L" << endl;
        cerr << "\tmax exp shape = exp(" << max_exponential_shape_intercept << " + " << max_exponential_shape_slope << " * L)" << endl;
#endif
        
        // reset mapping parameters to their original values
        adjust_alignments_for_base_quality = reset_quality_adjustments;
        min_clustering_mem_length = reset_min_clustering_mem_length;
        min_mem_length = reset_min_mem_length;
        max_alt_mappings = reset_max_alt_mappings;
        suppress_p_value_memoization = false;
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
        // We already threw out unmapped things
        assert(aln_1.path().mapping_size() != 0);
        pos_t pos_1 = initial_position(aln_1.path());
        assert(id(pos_1) != 0);
        
        Alignment aln_2;
        optimal_alignment(multipath_aln_2, aln_2);
        assert(aln_2.path().mapping_size() != 0);
        pos_t pos_2 = full_fragment ? final_position(aln_2.path()) : initial_position(aln_2.path());
        assert(id(pos_2) != 0);
#ifdef debug_multipath_mapper
        cerr << "measuring left-to-" << (full_fragment ? "right" : "left") << " end distance between " << pos_1 << " and " << pos_2 << endl;
#endif
        
        int64_t dist;
        if (use_min_dist_clusterer || use_tvs_clusterer) {
            assert(!forward_strand);
            // measure the distance in both directions and choose the minimum (or the only) absolute distance
            int64_t forward_dist = distance_index->min_distance(pos_1, pos_2);
            int64_t reverse_dist = distance_index->min_distance(pos_2, pos_1);
            if (forward_dist == -1 && reverse_dist == -1) {
                // unreachable both ways, convert to the sentinel that the client code expects
                dist = numeric_limits<int64_t>::max();
            }
            else if (forward_dist == -1 || (reverse_dist < forward_dist && reverse_dist != -1)) {
                dist = -reverse_dist;
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
                                                              bool block_rescue_from_1, bool block_rescue_from_2,
                                                              vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                              vector<pair<pair<size_t, size_t>, int64_t>>& pair_distances,
                                                              vector<double>& pair_multiplicities,
                                                              const match_fanouts_t* fanouts1, const match_fanouts_t* fanouts2) {
        
        // align the two ends independently
        vector<multipath_alignment_t> multipath_alns_1, multipath_alns_2;
        vector<size_t> cluster_idxs_1, cluster_idxs_2;
        if (!block_rescue_from_1) {
            cluster_idxs_1 = range_vector(cluster_graphs1.size());
            align_to_cluster_graphs(alignment1, mapping_quality_method == None ? Approx : mapping_quality_method,
                                    cluster_graphs1, multipath_alns_1, max_single_end_mappings_for_rescue,
                                    fanouts1, &cluster_idxs_1);
        }
        if (!block_rescue_from_2) {
            cluster_idxs_2 = range_vector(cluster_graphs2.size());
            align_to_cluster_graphs(alignment2, mapping_quality_method == None ? Approx : mapping_quality_method,
                                    cluster_graphs2, multipath_alns_2, max_single_end_mappings_for_rescue,
                                    fanouts2, &cluster_idxs_2);
        }
        
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
            pair_multiplicities.emplace_back(1.0);
            return true;
        }
        
        // figure out how many rescues we will do and could do from each side
        
        int32_t max_score_diff = get_aligner(!alignment1.quality().empty() &&
                                             !alignment2.quality().empty())->mapping_quality_score_diff(max_mapping_quality);
        
        int32_t top_score_1 = multipath_alns_1.empty() ? 0 : optimal_alignment_score(multipath_alns_1.front());
        int32_t top_score_2 = multipath_alns_2.empty() ? 0 : optimal_alignment_score(multipath_alns_2.front());
        
        size_t num_rescuable_alns_1 = block_rescue_from_1 ? 0 : multipath_alns_1.size();
        size_t num_rescuable_alns_2 = block_rescue_from_2 ? 0 : multipath_alns_2.size();
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
        
        // label all current pair multiplicities (usually none) with multiplicity 1
        pair_multiplicities.clear();
        pair_multiplicities.resize(multipath_aln_pairs_out.size(), 1.0);
        
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
                    if (abs(distance_between(multipath_alns_1[i], rescue_multipath_alns_1[j])) < 20) {
#ifdef debug_multipath_mapper
                        cerr << "found duplicate, now checking rescued read2 " << i << " and mapped read2 " << j << endl;
#endif
                        if (abs(distance_between(rescue_multipath_alns_2[i], multipath_alns_2[j])) < 20) {
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
                                pair_multiplicities.emplace_back(min(estimated_multiplicity_from_1, estimated_multiplicity_from_2));
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
                        pair_distances.emplace_back(make_pair(cluster_idxs_1[i], cluster_graphs2.size()), dist);
                        pair_multiplicities.emplace_back(estimated_multiplicity_from_1);
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
                    pair_distances.emplace_back(make_pair(cluster_graphs1.size(), cluster_idxs_2[j]), dist);
                    pair_multiplicities.emplace_back(estimated_multiplicity_from_2);
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
                    pair_distances.emplace_back(make_pair(cluster_idxs_1[i], cluster_graphs2.size()), dist);
                    pair_multiplicities.emplace_back(estimated_multiplicity_from_1);
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
                    pair_distances.emplace_back(make_pair(cluster_graphs1.size(), cluster_idxs_2[i]), dist);
                    pair_multiplicities.emplace_back(estimated_multiplicity_from_2);
                    found_consistent = true;
                }
            }
        }
        
        if (found_consistent) {
            // compute the paired mapping quality
            sort_and_compute_mapping_quality(multipath_aln_pairs_out, pair_distances, nullptr, &pair_multiplicities);
        }
        else {
#ifdef debug_multipath_mapper
            cerr << "failed to successfully rescue from either read end, reporting independent mappings" << endl;
#endif
            
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
        
        if (mapping_quality_method == None) {
            for (pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : multipath_aln_pairs_out) {
                multipath_aln_pair.first.set_mapping_quality(0);
                multipath_aln_pair.second.set_mapping_quality(0);
            }
        }
        
        return found_consistent;
    }
    
    void MultipathMapper::attempt_rescue_for_secondaries(const Alignment& alignment1, const Alignment& alignment2,
                                                         vector<clustergraph_t>& cluster_graphs1,
                                                         vector<clustergraph_t>& cluster_graphs2,
                                                         vector<pair<size_t, size_t>>& duplicate_pairs,
                                                         vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                         vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
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
                multipath_align(anchor_aln, get<0>(cluster_graphs[i]), get<1>(cluster_graphs[i]),
                                cluster_multipath_alns.back(), anchor_fanouts);
                
                // split it up if it turns out to be multiple components
                split_multicomponent_alignments(cluster_multipath_alns);
                
                // order the subpaths
                for (multipath_alignment_t& multipath_aln : cluster_multipath_alns) {
                    topologically_order_subpaths(multipath_aln);
                }
                
                // if we split it up, move the best one to the front
                if (cluster_multipath_alns.size() > 1) {
                    sort_and_compute_mapping_quality(cluster_multipath_alns, None);
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
                                rescued_secondaries.emplace_back(move(cluster_multipath_alns.front()), move(rescue_multipath_aln));
                                rescued_distances.emplace_back(make_pair(i, cluster_graphs2.size()), dist);
                                
                            }
                        }
                        else {
                            int64_t dist = distance_between(rescue_multipath_aln, cluster_multipath_alns.front(), true);
                            if (dist >= 0 && dist != numeric_limits<int64_t>::max()) {
                                rescued_secondaries.emplace_back(move(rescue_multipath_aln), move(cluster_multipath_alns.front()));
                                rescued_distances.emplace_back(make_pair(cluster_graphs1.size(), i), dist);
                                
                            }
                        }
                    } else {
#ifdef debug_multipath_mapper
                        cerr << "rescue failed" << endl;
#endif 
                    }
                } else {
#ifdef debug_multipath_mapper
                    cerr << "alignment we're rescuing from is likely a mismapping" << endl;
#endif
                }
            }
            
            // estimate how many of these alignments there probably are in total
            double multiplicity = double(num_rescuable) / double(num_rescues);
            
            // fill out the multiplicity vector to match
            for (size_t i = num_preexisting_pairs; i < rescued_secondaries.size(); ++i) {
                rescued_multiplicities.push_back(multiplicity);
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
                    if (abs(distance_between(rescued_secondaries[i].first, rescued_secondaries[j].first)) < 20) {
                        if (abs(distance_between(rescued_secondaries[i].second, rescued_secondaries[j].second)) < 20) {
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
            merge_rescued_mappings(multipath_aln_pairs_out, cluster_pairs, rescued_secondaries, rescued_distances,
                                   rescued_multiplicities);
        } else {
#ifdef debug_multipath_mapper
            cerr << "no rescues succeeded" << endl;
#endif
        }
    }
    
    void MultipathMapper::multipath_map_paired(const Alignment& alignment1, const Alignment& alignment2,
                                               vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                               vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer) {
                
#ifdef debug_multipath_mapper
        cerr << "multipath mapping paired reads " << pb2json(alignment1) << " and " << pb2json(alignment2) << endl;
#endif
        
        // empty the output vector (just for safety)
        multipath_aln_pairs_out.clear();
        
        if (!fragment_length_distr.is_finalized()) {
            // we have not estimated a fragment length distribution yet, so we revert to single ended mode and look
            // for unambiguous pairings
            
#ifdef debug_multipath_mapper
            cerr << "no fragment length distribution yet, looking for unambiguous single ended pairs" << endl;
#endif
            
            attempt_unpaired_multipath_map_of_pair(alignment1, alignment2, multipath_aln_pairs_out, ambiguous_pair_buffer);
            
            return;
        }
        
        // the fragment length distribution has been estimated, so we can do full-fledged paired mode
        vector<deque<pair<string::const_iterator, char>>> mem_fanouts1, mem_fanouts2;
        auto mems1 = find_mems(alignment1, &mem_fanouts1);
        auto mems2 = find_mems(alignment2, &mem_fanouts2);
        unique_ptr<match_fanouts_t> fanouts1(mem_fanouts1.empty() ? nullptr
                                             : new match_fanouts_t(record_fanouts(mems1, mem_fanouts1)));
        unique_ptr<match_fanouts_t> fanouts2(mem_fanouts2.empty() ? nullptr
                                             : new match_fanouts_t(record_fanouts(mems2, mem_fanouts2)));
                
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
        
        // TODO: this method of identifying a unique mapping is fragile to having a low
        // count noise MEM. it might be worth fixing this if we ever want to reenable it
        
        // find the count of the most unique match among the MEMs to assess how repetitive the sequence is
        size_t min_match_count_1 = numeric_limits<int64_t>::max();
        size_t min_match_count_2 = numeric_limits<int64_t>::max();
        for (const MaximalExactMatch& mem : mems1) {
            min_match_count_1 = min(min_match_count_1, mem.match_count);
        }
        for (const MaximalExactMatch& mem : mems2) {
            min_match_count_2 = min(min_match_count_2, mem.match_count);
        }
        
        // initialize cluster variables
        vector<memcluster_t> clusters1, clusters2;
        vector<clustergraph_t> cluster_graphs1, cluster_graphs2;
        vector<pair<pair<size_t, size_t>, int64_t>> cluster_pairs;
        vector<pair<size_t, size_t>> duplicate_pairs;
        
        MemoizingGraph memoizing_graph(xindex);
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
                                                                                                      &path_component_index));
        }
        
        // do we want to try to only cluster one read end and rescue the other?
        bool do_repeat_rescue_from_1 = min_match_count_2 > rescue_only_min && min_match_count_1 <= rescue_only_anchor_max;
        bool do_repeat_rescue_from_2 = min_match_count_1 > rescue_only_min && min_match_count_2 <= rescue_only_anchor_max;
        
        bool rescued_order_length_runs_1 = false, rescued_order_length_runs_2 = false;
        
#ifdef debug_multipath_mapper
        cerr << "min hit count on read 1: " << min_match_count_1 << ", on read 2: " << min_match_count_2 << ", doing rescue from read 1? " << (do_repeat_rescue_from_1 ? "yes" : "no") << ", from read 2? " << (do_repeat_rescue_from_2 ? "yes" : "no") << endl;
#endif
        
        if (do_repeat_rescue_from_1 || do_repeat_rescue_from_2) {
            
            // one side appears to be repetitive and the other non-repetitive, so try to only align the non-repetitive side
            // and get the other side from rescue
            
            // try to rescue high count runs of order-length MEMs for any read we're going to perform clustering on
            if (order_length_repeat_hit_max && do_repeat_rescue_from_1 && !rescued_order_length_runs_1) {
                rescue_high_count_order_length_mems(mems1, order_length_repeat_hit_max);
                rescued_order_length_runs_1 = true;
            }
            if (order_length_repeat_hit_max && do_repeat_rescue_from_2 && !rescued_order_length_runs_2) {
                rescue_high_count_order_length_mems(mems2, order_length_repeat_hit_max);
                rescued_order_length_runs_2 = true;
            }
            
            attempt_rescue_of_repeat_from_non_repeat(alignment1, alignment2, mems1, mems2, do_repeat_rescue_from_1, do_repeat_rescue_from_2,
                                                     clusters1, clusters2, cluster_graphs1, cluster_graphs2, multipath_aln_pairs_out,
                                                     cluster_pairs, *distance_measurer, fanouts1.get(), fanouts2.get());
            
            if (multipath_aln_pairs_out.empty() && do_repeat_rescue_from_1 && !do_repeat_rescue_from_2) {
                // we've clustered and extracted read 1, but rescue failed, so do the same for read 2 to prepare for the
                // normal pair clustering routine
                
#ifdef debug_multipath_mapper
                cerr << "repeat rescue failed from read 1, extracting clusters for read 2 and transitioning to standard clustering approach" << endl;
#endif
                
                // rescue high count runs of order-length MEMs now that we're going to cluster here
                if (order_length_repeat_hit_max && !rescued_order_length_runs_2) {
                    rescue_high_count_order_length_mems(mems2, order_length_repeat_hit_max);
                    rescued_order_length_runs_2 = true;
                }
                
                // do the clustering
                clusters2 = get_clusters(alignment2, mems2, &(*distance_measurer), fanouts2.get());
                cluster_graphs2 = query_cluster_graphs(alignment2, mems2, clusters2);
            }
            
            if (multipath_aln_pairs_out.empty() && do_repeat_rescue_from_2 && !do_repeat_rescue_from_1) {
                // we've clustered and extracted read 2, but rescue failed, so do the same for read 1 to prepare for the
                // normal pair clustering routine
                
#ifdef debug_multipath_mapper
                cerr << "repeat rescue failed from read 2, extracting clusters for read 1 and transitioning to standard clustering approach" << endl;
#endif
                
                // rescue high count runs of order-length MEMs now that we're going to cluster here
                if (order_length_repeat_hit_max && !rescued_order_length_runs_1) {
                    rescue_high_count_order_length_mems(mems1, order_length_repeat_hit_max);
                    rescued_order_length_runs_1 = true;
                }
                
                // do the clustering
                clusters1 = get_clusters(alignment1, mems1, &(*distance_measurer), fanouts1.get());
                cluster_graphs1 = query_cluster_graphs(alignment1, mems1, clusters1);
            }
        }
        else {
            // we have reasonably unique hits on both reads, so cluster them both and extract the graphs (i.e. don't try
            // to rely on rescue for either end yet)
            
#ifdef debug_multipath_mapper
            cerr << "clustering MEMs on both read ends..." << endl;
#endif
            
            // try to rescue high count runs of order-length MEMs for both reads before clustering
            if (order_length_repeat_hit_max && !rescued_order_length_runs_1) {
                rescue_high_count_order_length_mems(mems1, order_length_repeat_hit_max);
                rescued_order_length_runs_1 = true;
            }
            if (order_length_repeat_hit_max && !rescued_order_length_runs_2) {
                rescue_high_count_order_length_mems(mems2, order_length_repeat_hit_max);
                rescued_order_length_runs_2 = true;
            }
            
            // do the clustering
            clusters1 = get_clusters(alignment1, mems1, &(*distance_measurer), fanouts1.get());
            clusters2 = get_clusters(alignment2, mems2, &(*distance_measurer), fanouts2.get());
            
            // extract graphs around the clusters and get the assignments of MEMs to these graphs
            cluster_graphs1 = query_cluster_graphs(alignment1, mems1, clusters1);
            cluster_graphs2 = query_cluster_graphs(alignment2, mems2, clusters2);
        }
        
#ifdef debug_multipath_mapper
        cerr << "obtained independent clusters:" << endl;
        cerr << "read 1" << endl;
        for (int i = 0; i < cluster_graphs1.size(); i++) {
            cerr << "\tcluster " << i << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : get<1>(cluster_graphs1[i])) {
                cerr << "\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
            }
        }
        cerr << "read 2" << endl;
        for (int i = 0; i < cluster_graphs2.size(); i++) {
            cerr << "\tcluster " << i << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : get<1>(cluster_graphs2[i])) {
                cerr << "\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
            }
        }
#endif
        
        if (multipath_aln_pairs_out.empty()) {
            // we haven't already obtained a paired mapping by rescuing into a repeat, so we should try to get one
            // by cluster pairing
            
            cluster_pairs = get_cluster_pairs(alignment1, alignment2, cluster_graphs1, cluster_graphs2, &(*distance_measurer));
            
#ifdef debug_multipath_mapper
            cerr << "obtained cluster pairs:" << endl;
            for (int i = 0; i < cluster_pairs.size(); i++) {
                cerr << "\tpair "  << i << " at distance " << cluster_pairs[i].second << endl;
                cerr << "\t\t read 1 (cluster " << cluster_pairs[i].first.first <<  ")" << endl;
                for (pair<const MaximalExactMatch*, pos_t>  hit : get<1>(cluster_graphs1[cluster_pairs[i].first.first])) {
                    cerr << "\t\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
                }
                cerr << "\t\t read 2 (cluster " << cluster_pairs[i].first.second << ")" << endl;
                for (pair<const MaximalExactMatch*, pos_t>  hit : get<1>(cluster_graphs2[cluster_pairs[i].first.second])) {
                    cerr << "\t\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
                }
            }
#endif
            
            // do we find any pairs that satisfy the distance requirements?
            if (!cluster_pairs.empty()) {
                // We got some pairs that satisfy the distance requirements.
                
                // only perform the mappings that satisfy the expectations on distance
                
                align_to_cluster_graph_pairs(alignment1, alignment2, cluster_graphs1, cluster_graphs2, cluster_pairs,
                                             multipath_aln_pairs_out, duplicate_pairs, fanouts1.get(), fanouts2.get());
                
                // do we produce at least one good looking pair alignments from the clustered clusters?
                if (multipath_aln_pairs_out.empty() ? true : (likely_mismapping(multipath_aln_pairs_out.front().first) ||
                                                              likely_mismapping(multipath_aln_pairs_out.front().second))) {
                    
#ifdef debug_multipath_mapper
                    cerr << "one end of the pair may be mismapped, attempting individual end mappings" << endl;
#endif
                    // we're not happy with the pairs we got, try to get a good pair by rescuing from single ended alignments
                    // but block rescue from any sides that we already tried rescue from in the repeat rescue routine
                    
                    vector<pair<multipath_alignment_t, multipath_alignment_t>> rescue_aln_pairs;
                    vector<pair<pair<size_t, size_t>, int64_t>> rescue_distances;
                    vector<double> rescue_multiplicities;
                    bool rescued = align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                                       do_repeat_rescue_from_1, do_repeat_rescue_from_2,
                                                                       rescue_aln_pairs, rescue_distances, rescue_multiplicities,
                                                                       fanouts1.get(), fanouts2.get());
                    
                    // if we find consistent pairs by rescue, merge the two lists
                    if (rescued) {
#ifdef debug_multipath_mapper
                        cerr << "found some rescue pairs, merging into current list of consistent mappings" << endl;
#endif

                        merge_rescued_mappings(multipath_aln_pairs_out, cluster_pairs, rescue_aln_pairs, rescue_distances,
                                               rescue_multiplicities);
                        
                        // if we still haven't found mappings that are distinguishable from matches to random sequences,
                        // don't let them have any mapping quality
                        if (likely_mismapping(multipath_aln_pairs_out.front().first) ||
                            likely_mismapping(multipath_aln_pairs_out.front().second)) {
                            multipath_aln_pairs_out.front().first.set_mapping_quality(0);
                            multipath_aln_pairs_out.front().second.set_mapping_quality(0);
                        }
                        else {
                            // TODO: is this still necessary with the multiplicity code?
                            // also account for the possiblity that we selected the wrong ends to rescue with
                            cap_mapping_quality_by_rescue_probability(multipath_aln_pairs_out, cluster_pairs,
                                                                      cluster_graphs1, cluster_graphs2, false);
                            
                            // and for the possibility that we missed the correct cluster because of hit sub-sampling
                            // within the MEMs of the cluster
                            cap_mapping_quality_by_hit_sampling_probability(multipath_aln_pairs_out, cluster_pairs,
                                                                            cluster_graphs1, cluster_graphs2, false);
                        }
                    }
                    else {
                        // rescue didn't find any consistent mappings, revert to the single ended mappings
                        std::swap(multipath_aln_pairs_out, rescue_aln_pairs);
                        
                        // Don't sort and compute mapping quality; preserve the single-ended MAPQs
                    }
                }
                else {
                
                    // We don't think any of our hits are likely to be mismapped
                    
                    // does it look like we might be overconfident about this pair because of our clustering strategy
                    bool do_secondary_rescue = (multipath_aln_pairs_out.front().first.mapping_quality() >= max_mapping_quality - secondary_rescue_subopt_diff &&
                                                multipath_aln_pairs_out.front().second.mapping_quality() >= max_mapping_quality - secondary_rescue_subopt_diff);
                    
                    if (do_secondary_rescue) {
                        // we're very confident about this pair, but it might be because we over-pruned at the clustering stage
                        // so we use this routine to use rescue on other very good looking independent end clusters
                        attempt_rescue_for_secondaries(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                       duplicate_pairs, multipath_aln_pairs_out, cluster_pairs,
                                                       fanouts1.get(), fanouts2.get());
                        
                        // TODO: is this still necessary with the multiplicity code?
                        // account for the possiblity that we selected the wrong ends to rescue with
                        cap_mapping_quality_by_rescue_probability(multipath_aln_pairs_out, cluster_pairs,
                                                                  cluster_graphs1, cluster_graphs2, true);
                    }
                    
                    // account for the possibility that we missed the correct cluster because of hit sub-sampling
                    // within the MEMs of the cluster
                    cap_mapping_quality_by_hit_sampling_probability(multipath_aln_pairs_out, cluster_pairs,
                                                                    cluster_graphs1, cluster_graphs2, do_secondary_rescue);
                    
                }
            }
            else {
                // We got no pairs that satisfy the distance requirements
                
                // revert to independent single ended mappings, but skip any rescues that we already tried
                
#ifdef debug_multipath_mapper
                cerr << "could not find a consistent pair, reverting to single ended mapping" << endl;
#endif
                
                // have it record the multiplicities, even though we don't need them in thiss code path
                vector<double> multiplicities;
                
                bool rescued = align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                                   do_repeat_rescue_from_1, do_repeat_rescue_from_2,
                                                                   multipath_aln_pairs_out, cluster_pairs, multiplicities,
                                                                   fanouts1.get(), fanouts2.get());
                
                if (rescued) {
                    // We found valid pairs from rescue
                
                    // TODO: is this still necessary with the multiplicity code?
                    // account for the possiblity that we selected the wrong ends to rescue with
                    cap_mapping_quality_by_rescue_probability(multipath_aln_pairs_out, cluster_pairs,
                                                              cluster_graphs1, cluster_graphs2, false);
                }
                
                // Otherwise, we are just using single-end alignments, so don't
                // sort and compute mapping qualities. Leave the single-ended
                // MAPQs.
            }
        }
        
        if (multipath_aln_pairs_out.empty()) {
            // we tried all of our tricks and still didn't find a mapping
            
            // add a null alignment so we know it wasn't mapped
            multipath_aln_pairs_out.emplace_back();
            to_multipath_alignment(alignment1, multipath_aln_pairs_out.back().first);
            to_multipath_alignment(alignment2, multipath_aln_pairs_out.back().second);
            
            // in case we're realigning GAMs that have paths already
            multipath_aln_pairs_out.back().first.clear_subpath();
            multipath_aln_pairs_out.back().first.clear_start();
            multipath_aln_pairs_out.back().second.clear_subpath();
            multipath_aln_pairs_out.back().second.clear_start();
        }
        
        // if we computed extra alignments to get a mapping quality or investigate ambiguous clusters, remove them
        if (multipath_aln_pairs_out.size() > max_alt_mappings) {
            multipath_aln_pairs_out.resize(max_alt_mappings);
        }
        
        if (simplify_topologies) {
            for (pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : multipath_aln_pairs_out) {
                merge_non_branching_subpaths(multipath_aln_pair.first);
                merge_non_branching_subpaths(multipath_aln_pair.second);
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
        // TODO: make this machine-readable instead of a copy-able string.
        string distribution = "-I " + to_string(fragment_length_distr.mean()) + " -D " + to_string(fragment_length_distr.std_dev());
        
        for (pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : multipath_aln_pairs_out) {
            // Annotate with paired end distribution
            multipath_aln_pair.first.set_annotation("fragment_length_distribution", distribution);
            multipath_aln_pair.second.set_annotation("fragment_length_distribution", distribution);
        }
        
        // clean up the graph objects on the heap
        for (auto cluster_graph : cluster_graphs1) {
            delete get<0>(cluster_graph);
        }
        for (auto cluster_graph : cluster_graphs2) {
            delete get<0>(cluster_graph);
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
       
        if (mapping_quality_method != None) {
            // Now compute the MAPQ for the best alignment
            auto placement_mapq = compute_raw_mapping_quality_from_scores(scores, mapping_quality_method,
                                                                          !multipath_aln.quality().empty());
            // And min it in with what;s there already.
            alns_out[0].set_mapping_quality(min(alns_out[0].mapping_quality(), placement_mapq));
            for (size_t i = 1; i < alns_out.size(); i++) {
                // And zero all the others
                alns_out[i].set_mapping_quality(0);
            }
        }
    }
    
    void MultipathMapper::split_multicomponent_alignments(vector<multipath_alignment_t>& multipath_alns_out,
                                                          vector<size_t>* cluster_idxs) const {
        
        size_t num_original_alns = multipath_alns_out.size();
        for (size_t i = 0; i < num_original_alns; i++) {
            
            vector<vector<int64_t>> comps = connected_components(multipath_alns_out[i]);
            
            if (comps.size() > 1) {
#ifdef debug_multipath_mapper
                cerr << "splitting multicomponent alignment " << debug_string(multipath_alns_out[i]) << endl;
#endif
                // split this multipath alignment into its connected components
                for (size_t j = 1; j < comps.size(); j++) {
                    multipath_alns_out.emplace_back();
                    extract_sub_multipath_alignment(multipath_alns_out[i], comps[j],
                                                    multipath_alns_out.back());
                    // also label the split alignment with its cluster of origin, if we're keeping track of that
                    if (cluster_idxs) {
                        cluster_idxs->emplace_back(cluster_idxs->at(i));
                    }
                }
                // put the first component into the original location
                multipath_alignment_t last_component;
                extract_sub_multipath_alignment(multipath_alns_out[i], comps[0], last_component);
                multipath_alns_out[i] = move(last_component);
            }
        }
    }
    
    void MultipathMapper::attempt_rescue_of_repeat_from_non_repeat(const Alignment& alignment1, const Alignment& alignment2,
                                                                   const vector<MaximalExactMatch>& mems1, const vector<MaximalExactMatch>& mems2,
                                                                   bool do_repeat_rescue_from_1, bool do_repeat_rescue_from_2,
                                                                   vector<memcluster_t>& clusters1, vector<memcluster_t>& clusters2,
                                                                   vector<clustergraph_t>& cluster_graphs1, vector<clustergraph_t>& cluster_graphs2,
                                                                   vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                                   vector<pair<pair<size_t, size_t>, int64_t>>& pair_distances,
                                                                   OrientedDistanceMeasurer& distance_measurer,
                                                                   const match_fanouts_t* fanouts1, const match_fanouts_t* fanouts2) {
        
        bool rescue_succeeded_from_1 = false, rescue_succeeded_from_2 = false;
        
        vector<double> all_multiplicities;
        
        if (do_repeat_rescue_from_1) {
            
#ifdef debug_multipath_mapper
            cerr << "attempting repeat rescue from read 1" << endl;
#endif
            
            // get the clusters for the non repeat
            clusters1 = get_clusters(alignment1, mems1, &distance_measurer, fanouts1);
            
            // extract the graphs around the clusters
            cluster_graphs1 = query_cluster_graphs(alignment1, mems1, clusters1);
            
            // attempt rescue from these graphs
            vector<pair<multipath_alignment_t, multipath_alignment_t>> rescued_pairs;
            vector<pair<pair<size_t, size_t>, int64_t>> rescued_distances;
            vector<double> rescued_multiplicities;
            rescue_succeeded_from_1 = align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                                          false, true, rescued_pairs, pair_distances, rescued_multiplicities,
                                                                          fanouts1, fanouts2);
            
            // move the rescued pairs to the output vectors
            if (rescue_succeeded_from_1) {
                
#ifdef debug_multipath_mapper
                cerr << "repeat rescue succeeded" << endl;
#endif
                
                for (auto& multipath_aln_pair : rescued_pairs) {
                    multipath_aln_pairs_out.emplace_back(move(multipath_aln_pair));
                }
                for (auto& pair_distance : rescued_distances) {
                    pair_distances.emplace_back(pair_distance);
                }
                for (auto& multiplicity : rescued_multiplicities) {
                    all_multiplicities.emplace_back(multiplicity);
                }
            }
        }
        
        if (do_repeat_rescue_from_2) {
            // TODO: duplicative code
            
#ifdef debug_multipath_mapper
            cerr << "attempting repeat rescue from read 2" << endl;
#endif
            
            // get the clusters for the non repeat
            clusters2 = get_clusters(alignment2, mems2, &distance_measurer, fanouts2);
            
            // extract the graphs around the clusters
            cluster_graphs2 = query_cluster_graphs(alignment2, mems2, clusters2);
            
            // attempt rescue from these graphs
            vector<pair<multipath_alignment_t, multipath_alignment_t>> rescued_pairs;
            vector<pair<pair<size_t, size_t>, int64_t>> rescued_distances;
            vector<double> rescued_multiplicities;
            rescue_succeeded_from_2 = align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                                          true, false, rescued_pairs, pair_distances, rescued_multiplicities,
                                                                          fanouts1, fanouts2);
            
            // move the rescued pairs to the output vectors
            if (rescue_succeeded_from_2) {
                
#ifdef debug_multipath_mapper
                cerr << "repeat rescue succeeded" << endl;
#endif
                for (auto& multipath_aln_pair : rescued_pairs) {
                    multipath_aln_pairs_out.emplace_back(move(multipath_aln_pair));
                }
                for (auto& pair_distance : rescued_distances) {
                    pair_distances.emplace_back(pair_distance);
                }
                for (auto& multiplicity : rescued_multiplicities) {
                    all_multiplicities.emplace_back(multiplicity);
                }
            }
        }
        
        // re-sort the rescued alignments if we actually did it from both sides
        if (rescue_succeeded_from_1 && rescue_succeeded_from_2) {
            // TODO: add multiplicity logic if I reactivate this code path
            sort_and_compute_mapping_quality(multipath_aln_pairs_out, pair_distances, nullptr, &all_multiplicities);
        }
        
        // consider whether we should cap the mapping quality based on the chance that we rescued from the wrong clusters
        // TODO: is this still necessary with the multiplicity code?
        if (rescue_succeeded_from_1 || rescue_succeeded_from_2) {
            cap_mapping_quality_by_rescue_probability(multipath_aln_pairs_out, pair_distances, cluster_graphs1, cluster_graphs2, false);
        }
    }
    
    void MultipathMapper::merge_rescued_mappings(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                 vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                 vector<pair<multipath_alignment_t, multipath_alignment_t>>& rescued_multipath_aln_pairs,
                                                 vector<pair<pair<size_t, size_t>, int64_t>>& rescued_cluster_pairs,
                                                 vector<double>& rescued_multiplicities) const {
                
        size_t num_unrescued_pairs = multipath_aln_pairs_out.size();
        
        // give all the original mappings a multiplicity of 1
        vector<double> merged_multiplicities(num_unrescued_pairs, 1.0);
        
        for (size_t j = 0; j < rescued_multipath_aln_pairs.size(); j++) {
            
            // make sure this pair isn't a duplicate with any of the original pairs
            bool duplicate = false;
            for (size_t i = 0; i < num_unrescued_pairs; i++) {
#ifdef debug_multipath_mapper
                cerr << "checking if rescue pair " << j << " is duplicate of original pair " << i << endl;
#endif
                if (abs(distance_between(multipath_aln_pairs_out[i].first, rescued_multipath_aln_pairs[j].first)) < 20) {
                    if (abs(distance_between(multipath_aln_pairs_out[i].second, rescued_multipath_aln_pairs[j].second)) < 20) {
#ifdef debug_multipath_mapper
                        cerr << "found a duplicate" << endl;
#endif
                        duplicate = true;
                        // TODO: should i take the max of the multiplicities, or let it stay at 1?
                        break;
                    }
                }
            }
            
            // TODO: does it make more sense to give a duplicate a multiplicity of 1 or the estimated rescue multiplicity?
            // currently using 1
            
            if (!duplicate) {
#ifdef debug_multipath_mapper
                cerr << "no duplicate, adding to return vector if distance is finite and positive" << endl;
#endif
                cluster_pairs.emplace_back(rescued_cluster_pairs[j]);
                multipath_aln_pairs_out.emplace_back(move(rescued_multipath_aln_pairs[j]));
                merged_multiplicities.emplace_back(rescued_multiplicities[j]);
            }
        }
        
        sort_and_compute_mapping_quality(multipath_aln_pairs_out, cluster_pairs, nullptr, &merged_multiplicities);
    }
    
    void MultipathMapper::cap_mapping_quality_by_rescue_probability(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                                    vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                                    vector<clustergraph_t>& cluster_graphs1,
                                                                    vector<clustergraph_t>& cluster_graphs2,
                                                                    bool from_secondary_rescue) const {
        
#ifdef debug_multipath_mapper
        cerr << "checking whether we should enter rescue mapping quality capping routine" << endl;
#endif
        
        // did we use an out-of-bounds cluster index to flag either end as coming from a rescue?
        bool opt_aln_1_is_rescued = cluster_pairs.front().first.first >= cluster_graphs1.size();
        bool opt_aln_2_is_rescued = cluster_pairs.front().first.second >= cluster_graphs2.size();
                
        // was the optimal cluster pair obtained by rescue?
        if (opt_aln_1_is_rescued || opt_aln_2_is_rescued) {
            // let's figure out if we should reduce its mapping quality to reflect the fact that we may not have selected the
            // correct cluster as a rescue candidate
            
#ifdef debug_multipath_mapper
            cerr << "the optimal alignment is a rescue, checking if we need to cap the mapping quality" << endl;
#endif
            
            vector<clustergraph_t>& anchor_clusters = opt_aln_1_is_rescued ? cluster_graphs2 : cluster_graphs1;
            size_t anchor_idx = opt_aln_1_is_rescued ? cluster_pairs.front().first.second : cluster_pairs.front().first.first;
            
            // find the range of clusters that could plausibly be about as good as the one that rescue succeeded from
            size_t plausible_clusters_end_idx = anchor_idx;
            for (; plausible_clusters_end_idx < anchor_clusters.size(); plausible_clusters_end_idx++) {
                if (get<2>(anchor_clusters[plausible_clusters_end_idx]) < get<2>(anchor_clusters[anchor_idx]) - plausible_rescue_cluster_coverage_diff) {
                    break;
                }
            }
            
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
            
            if (max_rescues_attempted_idx < plausible_clusters_end_idx) {
                // we didn't attempt rescue from all of the plausible clusters, so we have to account for the possibility that
                // we missed the mapping just because of the subset of rescues we attempted
                
                double fraction_considered = double(max_rescues_attempted_idx) / double(plausible_clusters_end_idx);
                int32_t rescue_mapq = round(prob_to_phred(1.0 - fraction_considered));
                
#ifdef debug_multipath_mapper
                cerr << "capping mapping quality at " << rescue_mapq << endl;
#endif
                
                // cap the mapping quality at this value
                multipath_aln_pairs_out.front().first.set_mapping_quality(min(multipath_aln_pairs_out.front().first.mapping_quality(),
                                                                              rescue_mapq));
                multipath_aln_pairs_out.front().second.set_mapping_quality(min(multipath_aln_pairs_out.front().second.mapping_quality(),
                                                                               rescue_mapq));
            }
        }
    }
    
    double MultipathMapper::prob_equivalent_clusters_hits_missed(const memcluster_t& cluster) const {
        // we will approximate the probability of missing a cluster with the same MEMs as this one due to hit sampling
        // by the probability of missing the correct hit for all of the MEMs in this one
        double prob_of_missing_all = 1.0;
        for (const pair<const MaximalExactMatch*, pos_t>& hit : cluster) {
            const MaximalExactMatch& mem = *hit.first;
            prob_of_missing_all *= (mem.queried_count >= mem.match_count ? 0.0 : 1.0 - double(mem.queried_count) / double(mem.match_count));
        }
        return prob_of_missing_all;
    }
    
    void MultipathMapper::cap_mapping_quality_by_hit_sampling_probability(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                                          vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                                          vector<clustergraph_t>& cluster_graphs1,
                                                                          vector<clustergraph_t>& cluster_graphs2,
                                                                          bool did_secondary_rescue) const {
        
#ifdef debug_multipath_mapper
        cerr << "checking whether to cap mapping quality based on hit sub-sampling" << endl;
#endif
        
        // did we use an out-of-bounds cluster index to flag either end as coming from a rescue?
        bool opt_aln_1_is_rescued = cluster_pairs.front().first.first >= cluster_graphs1.size();
        bool opt_aln_2_is_rescued = cluster_pairs.front().first.second >= cluster_graphs2.size();
        
        // what is the chance that we would have missed a cluster with the same MEMs because of hit sub-sampling
        double prob_missing_equiv_cluster_1 = 0.0, prob_missing_equiv_cluster_2 = 0.0;
        if (!opt_aln_1_is_rescued) {
            prob_missing_equiv_cluster_1 = prob_equivalent_clusters_hits_missed(get<1>(cluster_graphs1[cluster_pairs.front().first.first]));
        }
        if (!opt_aln_2_is_rescued) {
            prob_missing_equiv_cluster_2 = prob_equivalent_clusters_hits_missed(get<1>(cluster_graphs2[cluster_pairs.front().first.second]));
        }
        
#ifdef debug_multipath_mapper
        cerr << "estimate probability of missing correct cluster from hit sampling for read 1 as " << prob_missing_equiv_cluster_1 << " and read 2 as " << prob_missing_equiv_cluster_2 << endl;
#endif
        
        if (prob_missing_equiv_cluster_1 == 0.0 && prob_missing_equiv_cluster_2 == 0.0) {
            // we can bail out now if we don't think hit sub-sampling was a factor
            return;
        }
        
        // what is the chance that we would have missed a pair?
        double prob_missing_pair = 0.0;
        if (opt_aln_1_is_rescued) {
#ifdef debug_multipath_mapper
            cerr << "optimal mapping is a rescue on read 1, sampling error only applies to read 2" << endl;
#endif
            
            prob_missing_pair = prob_missing_equiv_cluster_2;
        }
        else if (opt_aln_2_is_rescued) {
#ifdef debug_multipath_mapper
            cerr << "optimal mapping is a rescue on read 2, sampling error only applies to read 1" << endl;
#endif
            
            prob_missing_pair = prob_missing_equiv_cluster_1;
        }
        else if (!did_secondary_rescue) {
            // in this case we are assured that there was no rescue (except perhaps a failed rescue
            // that was trying to improve over a likely mismapping--a rare case where we won't mind
            // being pessimistic about the mapping quality)
            
            prob_missing_pair = 1.0 - (1.0 - prob_missing_equiv_cluster_1) * (1.0 - prob_missing_equiv_cluster_2);
        }
        else {
            // the complicated case:
            
            // we did secondary rescue, so now we must account for the fact that we *could* have recovered a cluster that we missed
            // because of MEM sub-sampling through the secondary rescue code path.
            
            // we will first cap the estimate of the probability of missing an equivalent cluster by the probability of rescuing
            // from the other read's cluster (assuming the other cluster is correct)
            
#ifdef debug_multipath_mapper
            cerr << "we did secondary rescue, so we must also account for the ability to rescue a sub-sampled cluster" << endl;
#endif
            
            // checks if one cluster is a subset of another
            auto cluster_contained_in = [](const memcluster_t& sub, const memcluster_t& super) {
                size_t i = 0;
                if (sub.size() <= super.size()) {
                    // query_cluster_graphs sorts the cluster hits so that they are ordered by length
                    // and then lexicographically by read interval
                    for (size_t j = 0; i < sub.size() && j < super.size(); j++) {
                        const MaximalExactMatch* mem_sub = sub[i].first;
                        const MaximalExactMatch* mem_super = super[j].first;
                        if (mem_sub == mem_super) {
                            i++;
                        }
                        else if (mem_sub->length() > mem_super->length() || (mem_sub->length() == mem_super->length() &&
                                                                             make_pair(mem_sub->begin, mem_sub->end) > make_pair(mem_super->begin, mem_super->end))) {
                            // we've passed the place in the vector where this hit should have occurred
                            break;
                        }
                                 
                    }
                }
                // did we find all of the MEMs in the superset?
                return i == sub.size();;
            };
            
            // TODO: figuring out which clusters we rescued from requires me to duplicate the secondary rescue code...
            
            // keep track of which clusters already have consistent pairs (so we wouldn't have rescued from them)
            unordered_set<size_t> paired_clusters_1, paired_clusters_2;
            for (size_t i = 0; i < multipath_aln_pairs_out.size(); i++) {
                paired_clusters_1.insert(cluster_pairs[i].first.first);
                paired_clusters_2.insert(cluster_pairs[i].first.second);
            }
            
            // the score parameters from the secondary rescue code that we use to decide which clusters to rescue from
            auto aligner = get_aligner(multipath_aln_pairs_out.empty() ? false :
                                       !multipath_aln_pairs_out.front().first.quality().empty() && !multipath_aln_pairs_out.front().second.quality().empty());
            int32_t cluster_score_1 = aligner->match * get<2>(cluster_graphs1[cluster_pairs.front().first.first]);
            int32_t cluster_score_2 = aligner->match * get<2>(cluster_graphs2[cluster_pairs.front().first.second]);
            int32_t max_score_diff = secondary_rescue_score_diff * aligner->mapping_quality_score_diff(max_mapping_quality);
            
            // we'll count up how many of the equivalent clusters we would have done secondary rescue from
            size_t num_equivalent_1 = 0;
            size_t num_equivalent_2 = 0;
            size_t num_rescued_equivalent_1 = 0;
            size_t num_rescued_equivalent_2 = 0;
            
            size_t num_rescues_1 = 0;
            for (size_t i = 0; i < cluster_graphs1.size(); i++) {
                // is this an equivalent cluster?
                bool equiv = cluster_contained_in(get<1>(cluster_graphs1[cluster_pairs.front().first.first]), get<1>(cluster_graphs1[i]));
                // and would we have done used it for a secondary rescue?
                bool did_rescue = (!paired_clusters_1.count(i) && num_rescues_1 < secondary_rescue_attempts
                                   && get<2>(cluster_graphs1[i]) * aligner->match >= cluster_score_1 - max_score_diff);
                // keep track of the counts
                num_rescues_1 += did_rescue;
                num_equivalent_1 += equiv;
                num_rescued_equivalent_1 += did_rescue && equiv;
            }
            
            size_t num_rescues_2 = 0;
            for (size_t i = 0; i < cluster_graphs2.size(); i++) {
                // is this an equivalent cluster?
                bool equiv = cluster_contained_in(get<1>(cluster_graphs2[cluster_pairs.front().first.second]), get<1>(cluster_graphs2[i]));
                // and would we have done used it for a secondary rescue?
                bool did_rescue = (!paired_clusters_2.count(i) && num_rescues_2 < secondary_rescue_attempts
                                   && get<2>(cluster_graphs2[i]) * aligner->match >= cluster_score_2 - max_score_diff);
                // keep track of the counts
                num_rescues_2 += did_rescue;
                num_equivalent_2 += equiv;
                num_rescued_equivalent_2 += did_rescue && equiv;
            }
#ifdef debug_multipath_mapper
            cerr << "did secondary rescue from " << num_rescued_equivalent_1 << " of " << num_equivalent_1 << " equivalent clusters on read 1 and " << num_rescued_equivalent_2 << " of " << num_equivalent_2 << " equivalent clusters on read 2" << endl;
#endif
            
            // compute the probability that we failed to rescue each of the clusters from the other
            double prob_missed_rescue_of_cluster_1 = 1.0 - (1.0 - prob_missing_equiv_cluster_2) * double(num_rescued_equivalent_2) / double(num_equivalent_2);
            double prob_missed_rescue_of_cluster_2 = 1.0 - (1.0 - prob_missing_equiv_cluster_1) * double(num_rescued_equivalent_1) / double(num_equivalent_1);
            
#ifdef debug_multipath_mapper
            cerr << "estimate the probability we missed rescue of cluster 1 at " << prob_missed_rescue_of_cluster_1 << " and cluster 2 at " << prob_missed_rescue_of_cluster_2 << endl;
#endif
            
            // possibly lower our estimate of the probabilty of missing a cluster if there was a high probabilty
            // that we would have found it during secondary rescue
            prob_missing_equiv_cluster_1 = min(prob_missed_rescue_of_cluster_1, prob_missing_equiv_cluster_1);
            prob_missing_equiv_cluster_2 = min(prob_missed_rescue_of_cluster_2, prob_missing_equiv_cluster_2);
            
            prob_missing_pair = 1.0 - (1.0 - prob_missing_equiv_cluster_1) * (1.0 - prob_missing_equiv_cluster_2);
        }
        
#ifdef debug_multipath_mapper
        cerr << "estimate probability that hit sampling caused us to miss correct cluster pair at " << prob_missing_pair << endl;
#endif
        
        if (prob_missing_pair > 0.0) {
            int32_t hit_sampling_mapq = round(prob_to_phred(prob_missing_pair));
            
#ifdef debug_multipath_mapper
            cerr << "capping mapping quality to " << hit_sampling_mapq << endl;
#endif
            
            // cap the mapping quality at this value
            multipath_aln_pairs_out.front().first.set_mapping_quality(min(multipath_aln_pairs_out.front().first.mapping_quality(),
                                                                          hit_sampling_mapq));
            multipath_aln_pairs_out.front().second.set_mapping_quality(min(multipath_aln_pairs_out.front().second.mapping_quality(),
                                                                           hit_sampling_mapq));
        }
    }
    
    void MultipathMapper::cap_mapping_quality_by_hit_sampling_probability(vector<multipath_alignment_t>& multipath_alns_out,
                                                                          vector<size_t>& cluster_idxs,
                                                                          vector<clustergraph_t>& cluster_graphs) const {
        
        clustergraph_t& opt_cluster = cluster_graphs[cluster_idxs.front()];
        
        // what is the chance that we would have missed a cluster with the same MEMs because of hit sub-sampling
        double prob_missing_equiv_cluster = prob_equivalent_clusters_hits_missed(get<1>(opt_cluster));
        
        
#ifdef debug_multipath_mapper
        cerr << "estimate probability that hit sampling caused us to miss correct cluster at " << prob_missing_equiv_cluster << endl;
#endif
        
        if (prob_missing_equiv_cluster > 0.0) {
            int32_t hit_sampling_mapq = round(prob_to_phred(prob_missing_equiv_cluster));
            
#ifdef debug_multipath_mapper
            cerr << "capping mapping quality to " << hit_sampling_mapq << endl;
#endif
            
            // cap the mapping quality at this value
            multipath_alns_out.front().set_mapping_quality(min(multipath_alns_out.front().mapping_quality(),
                                                               hit_sampling_mapq));
        }
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
    
    void MultipathMapper::split_multicomponent_alignments(vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                          vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs) const {
        
        size_t original_num_pairs = multipath_aln_pairs_out.size();
        for (size_t i = 0; i < original_num_pairs; i++) {
            vector<vector<int64_t>> connected_components_1 = connected_components(multipath_aln_pairs_out[i].first);
            vector<vector<int64_t>> connected_components_2 = connected_components(multipath_aln_pairs_out[i].second);
            
#ifdef debug_multipath_mapper
            cerr << "finding connected components for mapping:" << endl;
            view_multipath_alignment_as_dot(cerr, multipath_aln_pairs_out[i].first);
            view_multipath_alignment_as_dot(cerr, multipath_aln_pairs_out[i].second);
            //cerr  << pb2json(multipath_aln_pairs_out[i].first) << endl;
            //cerr  << pb2json(multipath_aln_pairs_out[i].second) << endl;
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
                    int64_t dist = distance_between(split_multipath_aln_pair.first, split_multipath_aln_pair.second,
                                                    true);
                    
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
                        }
                        else {
                            // append the rest of them to the end
                            multipath_aln_pairs_out.emplace_back(move(split_multipath_aln_pair));
                            cluster_pairs.emplace_back(cluster_pairs[i].first, dist);
                        }
                    }
                }
            }
        }
    }
    
    void MultipathMapper::align_to_cluster_graph_pairs(const Alignment& alignment1, const Alignment& alignment2,
                                                       vector<clustergraph_t>& cluster_graphs1,
                                                       vector<clustergraph_t>& cluster_graphs2,
                                                       vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                       vector<pair<multipath_alignment_t, multipath_alignment_t>>& multipath_aln_pairs_out,
                                                       vector<pair<size_t, size_t>>& duplicate_pairs_out,
                                                       const match_fanouts_t* fanouts1, const match_fanouts_t* fanouts2) {
        
        assert(multipath_aln_pairs_out.empty());
        
        auto aligner = get_aligner(!alignment1.quality().empty() && !alignment2.quality().empty());
        auto get_pair_approx_likelihood = [&](const pair<pair<size_t, size_t>, int64_t>& cluster_pair) {
            return ((get<2>(cluster_graphs1[cluster_pair.first.first])
                     + get<2>(cluster_graphs2[cluster_pair.first.second])) * aligner->match * aligner->log_base
                    + fragment_length_log_likelihood(cluster_pair.second));
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
        size_t num_mappings_to_compute = mapping_quality_method != None ? max(num_mapping_attempts, (size_t) 2) : num_mapping_attempts;
        
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
            
            // if we have a cluster graph pair with small enough MEM coverage
            // compared to the best one or we've made the maximum number of
            // alignments we stop producing alternate alignments
            if (get_pair_approx_likelihood(cluster_pair) < mem_coverage_min_ratio * get_pair_approx_likelihood(cluster_pairs.front())
                || num_mappings >= num_mappings_to_compute) {
                
                // remove the rest of the cluster pairs so we maintain the invariant that there are the
                // same number of cluster pairs as alternate mappings
                cluster_pairs.resize(i);
                
                break;
            }
#ifdef debug_multipath_mapper
            cerr << "doing pair " << cluster_pair.first.first << " " << cluster_pair.first.second << endl;
#endif
            // create multipath alignments to fill
            multipath_aln_pairs_out.emplace_back();
            
            auto prev_1 = previous_multipath_alns_1.find(cluster_pair.first.first);
            if (prev_1 == previous_multipath_alns_1.end()) {
                // we haven't done this alignment yet, so we have to complete it for the first time
                bdsg::HashGraph* graph1 = get<0>(cluster_graphs1[cluster_pair.first.first]);
                memcluster_t& graph_mems1 = get<1>(cluster_graphs1[cluster_pair.first.first]);
                
#ifdef debug_multipath_mapper
                cerr << "performing alignment of read 1 to subgraph" << endl;
#endif
                
                multipath_align(alignment1, graph1, graph_mems1, multipath_aln_pairs_out.back().first,
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
                auto graph2 = get<0>(cluster_graphs2[cluster_pair.first.second]);
                memcluster_t& graph_mems2 = get<1>(cluster_graphs2[cluster_pair.first.second]);
                
#ifdef debug_multipath_mapper
                cerr << "performing alignment of read 2 to subgraph" << endl;
#endif
                
                multipath_align(alignment2, graph2, graph_mems2, multipath_aln_pairs_out.back().second,
                                fanouts2);
                
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
        
        // split up any multi-component multipath alignments
        split_multicomponent_alignments(multipath_aln_pairs_out, cluster_pairs);
        
        // downstream algorithms assume multipath alignments are topologically sorted (including the scoring
        // algorithm in the next step)
        for (pair<multipath_alignment_t, multipath_alignment_t>& multipath_aln_pair : multipath_aln_pairs_out) {
            topologically_order_subpaths(multipath_aln_pair.first);
            topologically_order_subpaths(multipath_aln_pair.second);
        }
        
        // put pairs in score sorted order and compute mapping quality of best pair using the score
        sort_and_compute_mapping_quality(multipath_aln_pairs_out, cluster_pairs, &duplicate_pairs_out);
        
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

    pair<bdsg::HashGraph*, bool> MultipathMapper::extract_maximal_graph(const Alignment& alignment, const memcluster_t& cluster) {
        
        // Figure out the aligner to use
        auto aligner = get_aligner(!alignment.quality().empty());
        
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
        
        auto cluster_graph = new bdsg::HashGraph();
        algorithms::extract_containing_graph(xindex, cluster_graph, positions, forward_max_dist, backward_max_dist,
                                             num_alt_alns > 1 ? reversing_walk_length : 0);
        
        return make_pair(cluster_graph, cluster.size() == 1);
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

    pair<bdsg::HashGraph*, bool> MultipathMapper::extract_restrained_graph(const Alignment& alignment, const memcluster_t& cluster) {
        
        // Figure out the aligner to use
        auto aligner = get_aligner(!alignment.quality().empty());
        
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
        
        bdsg::HashGraph* cluster_graph = nullptr;
        bool do_extract = true;
        bool connected = false;
        while (do_extract) {
            
            // get rid of the old graph (if there is one)
            delete cluster_graph;
            
            // extract according to the current search distances
            cluster_graph = new bdsg::HashGraph();
            algorithms::extract_containing_graph(xindex, cluster_graph, positions, forward_dist, backward_dist,
                                                 num_alt_alns > 1 ? reversing_walk_length : 0);
            
            // we can avoid a costly algorithm when the cluster was extracted from one position (and therefore
            // must be connected)
            if (cluster.size() == 1 || algorithms::is_weakly_connected(cluster_graph)) {
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
        return make_pair(cluster_graph, connected);
    }

    pair<bdsg::HashGraph*, bool> MultipathMapper::extract_cluster_graph(const Alignment& alignment, const memcluster_t& cluster) {
        if (restrained_graph_extraction) {
            return extract_restrained_graph(alignment, cluster);
        }
        else {
            return extract_maximal_graph(alignment, cluster);
        }
    }
    
    auto MultipathMapper::query_cluster_graphs(const Alignment& alignment,
                                               const vector<MaximalExactMatch>& mems,
                                               const vector<memcluster_t>& clusters) -> vector<clustergraph_t> {
        
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
        unordered_map<size_t, pair<bdsg::HashGraph*, bool>> cluster_graphs;
        
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
            auto cluster_graph = extract_cluster_graph(alignment, cluster);
            
            // check if this subgraph overlaps with any previous subgraph (indicates a probable clustering failure where
            // one cluster was split into multiple clusters)
            unordered_set<size_t> overlapping_graphs;
            
            if (!do_merge_suppression) {
                cluster_graph.first->for_each_handle([&](const handle_t& handle) {
                    id_t node_id = cluster_graph.first->get_id(handle);
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
                for (auto& mem_hit : cluster) {
                    hit_to_cluster[mem_hit] = i;
                }
            }
            
            if (overlapping_graphs.empty()) {
                // there is no overlap with any other graph, suggesting a new unique hit
                
#ifdef debug_multipath_mapper
                cerr << "cluster graph does not overlap with any other cluster graphs, adding as cluster " << i << endl;
#endif
                cluster_graphs[i] = cluster_graph;
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
                if (remaining_idx == i) {
                    // the new graph was chosen to remain, so add it to the record
                    cluster_graphs[i] = cluster_graph;
                    merging_graph = cluster_graph.first;
                    all_connected = cluster_graph.second;
                }
                else {
                    // the new graph will be merged into an existing graph
                    merging_graph = cluster_graphs[remaining_idx].first;
                    all_connected = cluster_graphs[remaining_idx].second && cluster_graph.second;
                    algorithms::extend(cluster_graph.first, merging_graph);
                    delete cluster_graph.first;
                }
                
                // merge any other chained graphs into the remaining graph
                for (size_t j : overlapping_graphs) {
                    if (j != remaining_idx) {
                        auto removing_graph = cluster_graphs[j];
                        algorithms::extend(removing_graph.first, merging_graph);
                        all_connected = all_connected && removing_graph.second;
                        delete removing_graph.first;
                        cluster_graphs.erase(j);
                    }
                }
                cluster_graphs[remaining_idx].second = all_connected;
                
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
            if (!cluster_graph.second.second) {
                vector<unordered_set<id_t>> connected_components = algorithms::weakly_connected_components(cluster_graph.second.first);
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
                cluster_graphs[max_graph_idx + i] = make_pair(new bdsg::HashGraph(), true);
            }
            
            // divvy up the nodes
            auto joined_graph = cluster_graphs[multicomponent_graph.first].first;
            joined_graph->for_each_handle([&](const handle_t& handle) {
                for (size_t j = 0; j < multicomponent_graph.second.size(); j++) {
                    if (multicomponent_graph.second[j].count(joined_graph->get_id(handle))) {
                        cluster_graphs[max_graph_idx + j].first->create_handle(joined_graph->get_sequence(handle),
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
                        auto comp_graph = cluster_graphs[max_graph_idx + j].first;
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
            delete cluster_graphs[multicomponent_graph.first].first;
            cluster_graphs.erase(multicomponent_graph.first);
            
            if (do_merge_suppression) {
                // we need to re-assign the hits to the new cluster graphs
                for (auto& mem_hit : clusters[multicomponent_graph.first]) {
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
        for (const auto& cluster_graph : cluster_graphs) {
#ifdef debug_multipath_mapper
            cerr << "adding cluster graph " << cluster_graph.first << " to return vector at index " << cluster_graphs_out.size() << endl;
#endif
            cluster_to_idx[cluster_graph.first] = cluster_graphs_out.size();
            cluster_graphs_out.emplace_back(cluster_graph.second.first, memcluster_t(), 0);
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
                        get<1>(cluster_graphs_out[cluster_idx]).push_back(make_pair(&mem, make_pos_t(hit)));
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
                auto cluster_graph = get<0>(cluster_graphs_out[i]);
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
                        get<1>(cluster_graphs_out[cluster_to_idx[iter->second]]).push_back(mem_hit);
#ifdef debug_multipath_mapper
                        cerr << "\tMEM " << mem.sequence() << " at " << mem_hit.second << " assigned as seed to cluster at index " << cluster_to_idx[iter->second] << endl;
#endif
                    }
                    else {
                        // also try to find other, unassigned MEMs in the clusters
                        auto id_iter = node_id_to_cluster_idxs.find(id(mem_hit.second));
                        if (id_iter != node_id_to_cluster_idxs.end()) {
                            for (size_t cluster_idx : id_iter->second) {
                                get<1>(cluster_graphs_out[cluster_idx]).push_back(mem_hit);
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
            auto& cluster_graph = cluster_graphs_out[i];
            get<2>(cluster_graph) = read_coverage(get<1>(cluster_graph));
#ifdef debug_multipath_mapper
            cerr << "compute read coverage of cluster at index " << i << " to be " << get<2>(cluster_graph) << endl;
#endif
            sort(get<1>(cluster_graph).begin(), get<1>(cluster_graph).end(),
                 [](const pair<const MaximalExactMatch*, pos_t>& hit_1,
                    const pair<const MaximalExactMatch*, pos_t>& hit_2) {
                return hit_1.first->length() > hit_2.first->length() ||
                       (hit_1.first->length() == hit_2.first->length() &&
                        (hit_1.first->begin < hit_2.first->begin ||
                         (hit_1.first->begin == hit_2.first->begin && hit_1.first->end < hit_2.first->end)));
            });
        }
            
        // find the node ID range for the cluster graphs to help set up a stable, system-independent ordering
        // note: technically this is not quite a total ordering, but it should be close to one
        unordered_map<bdsg::HashGraph*, pair<id_t, id_t>> node_range;
        node_range.reserve(cluster_graphs_out.size());
        for (const auto& cluster_graph : cluster_graphs_out) {
            node_range[get<0>(cluster_graph)] = make_pair(get<0>(cluster_graph)->min_node_id(),
                                                          get<0>(cluster_graph)->max_node_id());
        }
        
        // sort the cluster graphs descending by unique sequence coverage, breaking ties by scrambling according to a hash
        stable_sort(cluster_graphs_out.begin(), cluster_graphs_out.end(),
                    [&](const clustergraph_t& cluster_graph_1,
                        const clustergraph_t& cluster_graph_2) {
                        return (get<2>(cluster_graph_1) > get<2>(cluster_graph_2) ||
                                (get<2>(cluster_graph_1) == get<2>(cluster_graph_2) &&
                                 wang_hash<pair<id_t, id_t>>()(node_range[get<0>(cluster_graph_1)]) < wang_hash<pair<id_t, id_t>>()(node_range[get<0>(cluster_graph_2)])));
                    });
        
        return move(cluster_graphs_out);
        
        
    }
    
    void MultipathMapper::multipath_align(const Alignment& alignment, const bdsg::HashGraph* graph,
                                          memcluster_t& graph_mems,
                                          multipath_alignment_t& multipath_aln_out,
                                          const match_fanouts_t* fanouts) const {

#ifdef debug_multipath_mapper_alignment
        cerr << "constructing alignment graph" << endl;
#endif
        
        // the longest path we could possibly align to (full gap and a full sequence)
        auto aligner = get_aligner(!alignment.quality().empty());
        size_t target_length = alignment.sequence().size() + min(aligner->longest_detectable_gap(alignment), max_alignment_gap);
        
        // check if we can get away with using only one strand of the graph
        bool use_single_stranded = algorithms::is_single_stranded(graph);
        bool mem_strand = false;
        if (use_single_stranded) {
            mem_strand = is_rev(graph_mems[0].second);
            for (size_t i = 1; i < graph_mems.size(); i++) {
                if (is_rev(graph_mems[i].second) != mem_strand) {
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
        cerr << "initial alignment graph:" << endl;
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
        DagifiedGraph* dagified = nullptr;
        
        ExpandingOverlayGraph* align_dag = nullptr;
        if (algorithms::is_directed_acyclic(align_digraph)) {
            align_dag = &undagified;
        }
        else {
#ifdef debug_multipath_mapper_alignment
            cerr << "graph contains directed cycles, performing dagification" << endl;
#endif
            dagified = new DagifiedGraph(align_digraph, target_length);
            align_dag = dagified;
        }
        
#ifdef debug_multipath_mapper_alignment
        cerr << "final alignment graph:" << endl;
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
        
        function<pair<id_t, bool>(id_t)> translator = [&](const id_t& node_id) {
            handle_t original = align_digraph->get_underlying_handle(align_dag->get_underlying_handle(align_dag->get_handle(node_id)));
            return make_pair(graph->get_id(original), graph->get_is_reverse(original));
        };
        
#ifdef debug_multipath_mapper_alignment
        cerr << "making multipath alignment MEM graph" << endl;
#endif
        
        MultipathAlignmentGraph multi_aln_graph(*align_dag, graph_mems, translator, max_branch_trim_length, gcsa, fanouts);
        
        {
            // Compute a topological order over the graph
            vector<size_t> topological_order;
            multi_aln_graph.topological_sort(topological_order);
            
            // it's sometimes possible for transitive edges to survive the original construction algorithm, so remove them
            multi_aln_graph.remove_transitive_edges(topological_order);
            
            // prune this graph down the paths that have reasonably high likelihood
            multi_aln_graph.prune_to_high_scoring_paths(alignment, aligner,
                                                        max_suboptimal_path_score_ratio, topological_order);
        }
        
        if (snarl_manager) {
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
                              use_pessimistic_tail_alignment ? pessimistic_gap_multiplier : 0.0,
                              choose_band_padding, multipath_aln_out);
        
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
        
        // clean up (safe on nullptr if we didn't choose it)
        delete dagified;
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
                              use_pessimistic_tail_alignment ? pessimistic_gap_multiplier : 0.0,
                              choose_band_padding, multipath_aln_out);
        
        for (size_t j = 0; j < multipath_aln_out.subpath_size(); j++) {
            translate_oriented_node_ids(*multipath_aln_out.mutable_subpath(j)->mutable_path(), translator);
        }
        
        topologically_order_subpaths(multipath_aln_out);
        
#ifdef debug_multipath_mapper_alignment
        cerr << "completed multipath alignment: " << debug_string(multipath_aln_out) << endl;
#endif
    }
    
    int64_t MultipathMapper::read_coverage(const memcluster_t& mem_hits) {
        if (mem_hits.empty()) {
            return 0;
        }
        
        vector<pair<string::const_iterator, string::const_iterator>> mem_read_segments;
        mem_read_segments.reserve(mem_hits.size());
        for (auto& mem_hit : mem_hits) {
            mem_read_segments.emplace_back(mem_hit.first->begin, mem_hit.first->end);
        }
        std::sort(mem_read_segments.begin(), mem_read_segments.end());
        auto curr_begin = mem_read_segments[0].first;
        auto curr_end = mem_read_segments[0].second;
        
        int64_t total = 0;
        for (size_t i = 1; i < mem_read_segments.size(); i++) {
            if (mem_read_segments[i].first >= curr_end) {
                total += (curr_end - curr_begin);
                curr_begin = mem_read_segments[i].first;
                curr_end = mem_read_segments[i].second;
            }
            else if (mem_read_segments[i].second > curr_end) {
                curr_end = mem_read_segments[i].second;
            }
        }
        return total + (curr_end - curr_begin);
    }
    
    void MultipathMapper::strip_full_length_bonuses(multipath_alignment_t& multipath_aln) const {
        
        // TODO: this could technically be wrong if only one read in a pair has qualities
        int32_t full_length_bonus = get_aligner(!multipath_aln.quality().empty())->full_length_bonus;
        // strip bonus from source paths
        if (multipath_aln.start_size()) {
            // use the precomputed list of sources if we have it
            for (size_t i = 0; i < multipath_aln.start_size(); i++) {
                subpath_t* source_subpath = multipath_aln.mutable_subpath(multipath_aln.start(i));
                const edit_t& edit = source_subpath->path().mapping(0).edit(0);
                if (edit.to_length() > 0 && edit.from_length() == 0) {
                    source_subpath->set_score(source_subpath->score() - full_length_bonus);
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
                if (edit.to_length() > 0 && edit.from_length() == 0) {
                    source_subpath->set_score(source_subpath->score() - full_length_bonus);
                }
            }
        }
        // strip bonus from sink paths
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            subpath_t* subpath = multipath_aln.mutable_subpath(i);
            if (subpath->next_size() == 0) {
                const path_mapping_t& final_mapping = subpath->path().mapping(subpath->path().mapping_size() - 1);
                const edit_t& edit = final_mapping.edit(final_mapping.edit_size() - 1);
                if (edit.to_length() > 0 && edit.from_length() == 0) {
                    subpath->set_score(subpath->score() - full_length_bonus);
                }
            }
        }
    }
    
    int32_t MultipathMapper::compute_raw_mapping_quality_from_scores(const vector<double>& scores, MappingQualityMethod mapq_method,
                                                                     bool have_qualities, const vector<double>* multiplicities) const {
   
        // We should never actually compute a MAPQ with the None method. If we try, it means something has gone wrong.
        assert(mapq_method != None);
   
        auto aligner = get_aligner(have_qualities);
        int32_t raw_mapq;
        if (mapping_quality_method == Adaptive) {
            raw_mapq = aligner->compute_mapping_quality(scores, scores.size() < 2 ? true :
                                                        (scores[1] < scores[0] - get_aligner()->mapping_quality_score_diff(max_mapping_quality)),
                                                        multiplicities);
        }
        else {
            raw_mapq = aligner->compute_mapping_quality(scores, mapping_quality_method == Approx, multiplicities);
        }
        
        // arbitrary scaling, seems to help performance
        raw_mapq *= mapq_scaling_factor;
        
#ifdef debug_multipath_mapper
        cerr << "scores yield a raw MAPQ of " << raw_mapq << endl;
#endif

        return raw_mapq;

    }
    
    void MultipathMapper::sort_and_compute_mapping_quality(vector<multipath_alignment_t>& multipath_alns,
                                                           MappingQualityMethod mapq_method,
                                                           vector<size_t>* cluster_idxs) const {
        if (multipath_alns.empty()) {
            return;
        }
        
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
        vector<double> base_scores(multipath_alns.size(), 0.0);
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
            } else {
                // Just compute a single optimal alignment
                alignments = optimal_alignments(multipath_alns[i], 1);
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
            base_scores[i] = alignments.empty() ? 0 : alignments[0].score();
            
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
                double adjustment = best_linearization_total_score - base_scores[i];
                
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
        auto& scores = (include_population_component && all_multipaths_pop_consistent) ? pop_adjusted_scores : base_scores;
        
        // find the order of the scores
        vector<size_t> order(multipath_alns.size(), 0);
        for (size_t i = 1; i < multipath_alns.size(); i++) {
            order[i] = i;
        }
        // Sort, shuffling based on the aligned sequence to break ties.
        sort_shuffling_ties(order.begin(), order.end(),
            [&](const size_t i, const size_t j) { return scores[i] > scores[j]; },
            [&](const size_t seed_source) {return multipath_alns[seed_source].sequence(); });
        
        // translate the order to an index
        vector<size_t> index(multipath_alns.size());
        for (size_t i = 0; i < multipath_alns.size(); i++) {
            index[order[i]] = i;
        }
        
        // put the scores, clusters-of-origin, and alignments in order
        for (size_t i = 0; i < multipath_alns.size(); i++) {
            while (index[i] != i) {
                std::swap(scores[index[i]], scores[i]);
                std::swap(multipath_alns[index[i]], multipath_alns[i]);
                if (cluster_idxs) {
                    std::swap((*cluster_idxs)[index[i]], (*cluster_idxs)[i]);
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
            }
        }
        if (removed_so_far) {
            multipath_alns.resize(multipath_alns.size() - removed_so_far);
            scores.resize(scores.size() - removed_so_far);
            if (cluster_idxs) {
                cluster_idxs->resize(cluster_idxs->size() - removed_so_far);
            }
        }
        
#ifdef debug_multipath_mapper
        cerr << "scores obtained of multi-mappings:" << endl;
        for (size_t i = 0; i < scores.size(); i++) {
            Alignment aln;
            optimal_alignment(multipath_alns[i], aln);
            cerr << "\t" << scores[i] << " " << (aln.path().mapping_size() ? make_pos_t(aln.path().mapping(0).position()) : pos_t()) << endl;
        }
#endif
        
        if (mapq_method != None) {
            // Sometimes we are passed None, which means to not update the MAPQs at all. But otherwise, we do MAPQs.
            // Compute and set the mapping quality
            int32_t raw_mapq = compute_raw_mapping_quality_from_scores(scores, mapq_method, !multipath_alns.front().quality().empty());
            multipath_alns.front().set_mapping_quality(min<int32_t>(raw_mapq, max_mapping_quality));
        }
        
        if (report_group_mapq) {
            size_t num_reporting = min(multipath_alns.size(), max_alt_mappings);
            vector<size_t> reporting_idxs(num_reporting, 0);
            for (size_t i = 1; i < num_reporting; ++i) {
                reporting_idxs[i] = i;
            }
            double raw_mapq = get_aligner(!multipath_alns.front().quality().empty())->compute_group_mapping_quality(scores, reporting_idxs);
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
        
        if (multipath_aln_pairs.empty()) {
            return;
        }
        
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
        vector<double> base_scores(multipath_aln_pairs.size(), 0.0);
        
#ifdef debug_multipath_mapper
        // trackers that let us keep track of the alignment and population components separately
        vector<pair<double, double>> chosen_align_score(multipath_aln_pairs.size(),
                                                        make_pair(numeric_limits<double>::quiet_NaN(),
                                                                  numeric_limits<double>::quiet_NaN()));
        vector<pair<double, double>> chosen_population_score(multipath_aln_pairs.size(),
                                                             make_pair(numeric_limits<double>::quiet_NaN(),
                                                                       numeric_limits<double>::quiet_NaN()));
#endif
        
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
            vector<vector<Alignment>> alignments(2);
            
            if (query_population) {
                // We want to do population scoring
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
            } else {
                // Just compute a single optimal alignment
                alignments[0] = optimal_alignments(multipath_aln_pair.first, 1);
                alignments[1] = optimal_alignments(multipath_aln_pair.second, 1);
            }
            
#ifdef debug_multipath_mapper
            cerr << "Got " << alignments[0].size() << " and " << alignments[1].size() << " linearizations on each end" << endl;
#endif
            
            // We used to fail an assert if either list of optimal alignments
            // was empty, but now we handle it as if that side is an unmapped
            // read with score 0.
            
            // Compute the optimal alignment score ignoring population
            int32_t alignment_score = (alignments[0].empty() ? 0 : alignments[0][0].score()) +
                (alignments[1].empty() ? 0 : alignments[1][0].score());
            
            // This is the contribution to the alignment's score from the fragment length distribution
            double frag_score;
            if (alignments[0].empty() || alignments[1].empty()) {
                // Actually there should be no fragment score, because one or both ends are unmapped
                frag_score = 0;
            } else {
                // compute the fragment distribution's contribution to the score
                frag_score = fragment_length_log_likelihood(cluster_pairs[i].second) / log_base;
            }
            min_frag_score = min(frag_score, min_frag_score);
            
            // Record the base score, including fragment contribution
            base_scores[i] = alignment_score + frag_score;
            
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
#ifdef debug_multipath_mapper

                                if (end == 0) {
                                    chosen_align_score[i].first = alignments[end][j].score();
                                    chosen_population_score[i].first = pop_score.first / log_base;
                                }
                                else {
                                    chosen_align_score[i].second = alignments[end][j].score();
                                    chosen_population_score[i].second = pop_score.first / log_base;
                                }
#endif
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
        auto& scores = (include_population_component && all_multipaths_pop_consistent) ? pop_adjusted_scores : base_scores;
        
        for (auto& score : scores) {
            // Pull the min frag or extra score out of the score so it will be nonnegative
            score -= (include_population_component && all_multipaths_pop_consistent) ? min_extra_score : min_frag_score;
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
        
        // find the order of the scores
        vector<size_t> order(multipath_aln_pairs.size(), 0);
        for (size_t i = 1; i < multipath_aln_pairs.size(); i++) {
            order[i] = i;
        }
        sort_shuffling_ties(order.begin(), order.end(),
            [&](const size_t i, const size_t j) {
                return scores[i] > scores[j]; 
            },
            [&](const size_t seed_source) {
                return multipath_aln_pairs[seed_source].first.sequence() + multipath_aln_pairs[seed_source].second.sequence();
            });
        
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
#ifdef debug_multipath_mapper
                std::swap(chosen_align_score[index[i]], chosen_align_score[i]);
                std::swap(chosen_population_score[index[i]], chosen_population_score[i]);
#endif
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
            auto start1 = aln1.path().mapping(0).position().node_id();
            auto start2 = aln2.path().mapping(0).position().node_id();
        
            cerr << "\tpos:" << start1 << "(" << aln1.score() << ")-" << start2 << "(" << aln2.score() << ")"
                << " align:" << optimal_alignment_score(multipath_aln_pairs[i].first) + optimal_alignment_score(multipath_aln_pairs[i].second)
                << ", length: " << cluster_pairs[i].second;
            if (include_population_component && all_multipaths_pop_consistent) {
                cerr << ", pop adjusted aligns: " << chosen_align_score[i].first << " " << chosen_align_score[i].second << ", population: " << chosen_population_score[i].first << " " << chosen_population_score[i].second;
            }
            cerr << ", combined: " << scores[i];
            if (multiplicities) {
                cerr << ", multiplicity: " << multiplicities->at(i);
            }
            cerr << endl;
        }
#endif
        
        if (mapping_quality_method != None) {
            // Compute the raw mapping quality
            int32_t raw_mapq = compute_raw_mapping_quality_from_scores(scores, mapping_quality_method,
                                                                       !multipath_aln_pairs.front().first.quality().empty() &&
                                                                       !multipath_aln_pairs.front().second.quality().empty(),
                                                                       multiplicities);
            // Limit it to the max.
            int32_t mapq = min<int32_t>(raw_mapq, max_mapping_quality);
            multipath_aln_pairs.front().first.set_mapping_quality(mapq);
            multipath_aln_pairs.front().second.set_mapping_quality(mapq);
            
            if (multipath_aln_pairs.size() > 1) {
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
                    
                    // arbitrary scaling, seems to help performance
                    raw_mapq_1 *= mapq_scaling_factor;
                    raw_mapq_2 *= mapq_scaling_factor;
                    
#ifdef debug_multipath_mapper
                    cerr << "deduplicated raw MAPQs are " << raw_mapq_1 << " and " << raw_mapq_2 << endl;
#endif
                    
                    int32_t mapq_1 = min<int32_t>(raw_mapq_1, max_mapping_quality);
                    int32_t mapq_2 = min<int32_t>(raw_mapq_2, max_mapping_quality);
                    
                    multipath_aln_pairs.front().first.set_mapping_quality(mapq_1);
                    multipath_aln_pairs.front().second.set_mapping_quality(mapq_2);
                }
            }
        }
        
        if (report_group_mapq) {
            size_t num_reporting = min(multipath_aln_pairs.size(), max_alt_mappings);
            vector<size_t> reporting_idxs(num_reporting, 0);
            for (size_t i = 1; i < num_reporting; ++i) {
                reporting_idxs[i] = i;
            }
            auto aligner = get_aligner(!multipath_aln_pairs.front().first.quality().empty() &&
                                       !multipath_aln_pairs.front().second.quality().empty());
            double raw_mapq = aligner->compute_group_mapping_quality(scores, reporting_idxs);
            
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
    
    double MultipathMapper::read_coverage_z_score(int64_t coverage, const Alignment& alignment) const {
        /* algebraically equivalent to
         *
         *      Coverage - ReadLen / 4
         *  -------------------------------
         *  sqrt(ReadLen * 1/4 * (1 - 1/4))
         * 
         * from the Normal approximation to a Binomal(ReadLen, 1/4)
         */
        double root_len = sqrt(alignment.sequence().size());
        return 0.5773502691896258 * (4.0 * coverage / root_len - root_len);
    }
}




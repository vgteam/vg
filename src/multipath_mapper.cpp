//
//  multipath_mapper.cpp
//  
//
//

//#define debug_multipath_mapper_mapping
//#define debug_multipath_mapper_alignment
//#define debug_validate_multipath_alignments
//#define debug_report_startup_training

#include "multipath_mapper.hpp"

#include "algorithms/extract_containing_graph.hpp"
#include "algorithms/extract_connecting_graph.hpp"
#include "algorithms/extract_extending_graph.hpp"
#include "algorithms/topological_sort.hpp"
#include "algorithms/weakly_connected_components.hpp"

#include <structures/union_find.hpp>

namespace vg {
    
    //size_t MultipathMapper::PRUNE_COUNTER = 0;
    //size_t MultipathMapper::SUBGRAPH_TOTAL = 0;
    
    MultipathMapper::MultipathMapper(xg::XG* xg_index, gcsa::GCSA* gcsa_index, gcsa::LCPArray* lcp_array,
                                     haplo::ScoreProvider* haplo_score_provider, SnarlManager* snarl_manager) :
        BaseMapper(xg_index, gcsa_index, lcp_array, haplo_score_provider),
        snarl_manager(snarl_manager)
    {
        // nothing to do
    }
    
    MultipathMapper::~MultipathMapper() {
        
    }
    
    void MultipathMapper::multipath_map(const Alignment& alignment,
                                        vector<MultipathAlignment>& multipath_alns_out,
                                        size_t max_alt_mappings) {
        multipath_map_internal(alignment, mapping_quality_method, multipath_alns_out, max_alt_mappings);
    }
    
    void MultipathMapper::multipath_map_internal(const Alignment& alignment,
                                                 MappingQualityMethod mapq_method,
                                                 vector<MultipathAlignment>& multipath_alns_out,
                                                 size_t max_alt_mappings) {
        
#ifdef debug_multipath_mapper_mapping
        cerr << "multipath mapping read " << pb2json(alignment) << endl;
        cerr << "querying MEMs..." << endl;
#endif
    
        // query MEMs using GCSA2
        double dummy1; double dummy2;
        vector<MaximalExactMatch> mems = find_mems_deep(alignment.sequence().begin(), alignment.sequence().end(),
                                                        dummy1, dummy2, 0, min_mem_length, mem_reseed_length,
                                                        false, true, true, false);
        
#ifdef debug_multipath_mapper_mapping
        cerr << "obtained MEMs:" << endl;
        for (MaximalExactMatch mem : mems) {
            cerr << "\t" << mem << " (" << mem.nodes.size() << " hits)" << endl;
        }
        cerr << "clustering MEMs..." << endl;
#endif
        
        // TODO: use the automatic expected MEM length algorithm to restrict the MEMs used for clustering?
        
        // cluster the MEMs
        vector<memcluster_t> clusters;
        // memos for the results of expensive succinct operations that we may need to do multiple times
        OrientedDistanceClusterer::paths_of_node_memo_t paths_of_node_memo;
        OrientedDistanceClusterer::oriented_occurences_memo_t oriented_occurences_memo;
        OrientedDistanceClusterer::handle_memo_t handle_memo;
        // TODO: Making OrientedDistanceClusterers is the only place we actually
        // need to distinguish between regular_aligner and qual_adj_aligner
        if (adjust_alignments_for_base_quality) {
            OrientedDistanceClusterer clusterer(alignment, mems, *get_qual_adj_aligner(), xindex, max_expected_dist_approx_error,
                                                min_clustering_mem_length, unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
            clusters = clusterer.clusters(max_mapping_quality, log_likelihood_approx_factor);
        }
        else {
            OrientedDistanceClusterer clusterer(alignment, mems, *get_regular_aligner(), xindex, max_expected_dist_approx_error,
                                                min_clustering_mem_length, unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
            clusters = clusterer.clusters(max_mapping_quality, log_likelihood_approx_factor);
        }
        
#ifdef debug_multipath_mapper_mapping
        cerr << "obtained clusters:" << endl;
        for (int i = 0; i < clusters.size(); i++) {
            cerr << "\tcluster " << i << endl;
            for (pair<const MaximalExactMatch*, pos_t>  hit : clusters[i]) {
                cerr << "\t\t" << hit.second << " " <<  hit.first->sequence() << endl;
            }
        }
        cerr << "extracting subgraphs..." << endl;
#endif
        
        // extract graphs around the clusters
        auto cluster_graphs = query_cluster_graphs(alignment, mems, clusters);
        
        // actually perform the alignments and post-process to meeth MultipathAlignment invariants
        align_to_cluster_graphs(alignment, mapq_method, cluster_graphs, multipath_alns_out, max_alt_mappings);
        
        if (multipath_alns_out.empty()) {
            // add a null alignment so we know it wasn't mapped
            multipath_alns_out.emplace_back();
            to_multipath_alignment(alignment, multipath_alns_out.back());
            
            // in case we're realigning GAMs that have paths already
            multipath_alns_out.back().clear_subpath();
            multipath_alns_out.back().clear_start();
        }
        
        if (likely_mismapping(multipath_alns_out.front())) {
            // we can't distinguish this alignment from the longest MEM of a random sequence
#ifdef debug_multipath_mapper_mapping
            cerr << "mapping is not distinguishable from a random sequence, snapping MAPQ to 0" << endl;
#endif
            
            multipath_alns_out.front().set_mapping_quality(0);
        }
        
        // if we computed extra alignments to get a mapping quality, remove them
        if (multipath_alns_out.size() > max_alt_mappings) {
            multipath_alns_out.resize(max_alt_mappings);
        }
        
        if (strip_bonuses) {
            for (MultipathAlignment& multipath_aln : multipath_alns_out) {
                strip_full_length_bonuses(multipath_aln);
            }
        }
        
        // clean up the cluster graphs
        for (auto cluster_graph : cluster_graphs) {
            delete get<0>(cluster_graph);
        }
    }
    
    void MultipathMapper::align_to_cluster_graphs(const Alignment& alignment,
                                                  MappingQualityMethod mapq_method,
                                                  vector<clustergraph_t>& cluster_graphs,
                                                  vector<MultipathAlignment>& multipath_alns_out,
                                                  size_t max_alt_mappings) {
    
#ifdef debug_multipath_mapper_mapping
        cerr << "sorting subgraphs by read coverage..." << endl;
#endif
        
        
        
#ifdef debug_multipath_mapper_mapping
        cerr << "aligning to subgraphs..." << endl;
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
            if (get<2>(cluster_graph) < mem_coverage_min_ratio * get<2>(cluster_graphs[0])
                || num_mappings >= num_mappings_to_compute) {
#ifdef debug_multipath_mapper_mapping
                cerr << "halting further alignments, either because MEM coverage of " << get<2>(cluster_graph) << " is too far below optimum of " << get<2>(cluster_graphs[0]) << " or because already made " << num_mappings << " of " << num_mappings_to_compute << " mappings" << endl;
#endif
                
//#pragma omp atomic
//                PRUNE_COUNTER += cluster_graphs.size() - num_mappings;
                break;
            }
            
#ifdef debug_multipath_mapper_alignment
            cerr << "performing alignment to subgraph " << pb2json(get<0>(cluster_graph)->graph) << endl;
#endif
            
            multipath_alns_out.emplace_back();
            multipath_align(alignment, get<0>(cluster_graph), get<1>(cluster_graph), multipath_alns_out.back());
            
            num_mappings++;
        }
        
        // split up any alignments that ended up being disconnected
        split_multicomponent_alignments(multipath_alns_out);
        
#ifdef debug_multipath_mapper_mapping
        cerr << "topologically ordering " << multipath_alns_out.size() << " multipath alignments" << endl;
#endif
        for (MultipathAlignment& multipath_aln : multipath_alns_out) {
            topologically_order_subpaths(multipath_aln);
        }
        
#ifdef debug_multipath_mapper_mapping
        cerr << "computing mapping quality and sorting mappings" << endl;
#endif
        sort_and_compute_mapping_quality(multipath_alns_out, mapq_method);
        
        if (!multipath_alns_out.empty() ? likely_mismapping(multipath_alns_out.front()) : false) {
            multipath_alns_out.front().set_mapping_quality(0);
        }
        
        // for debugging: an expensive check for invariant validity that can be turned on
        // with a preprocessor flag
#ifdef debug_validate_multipath_alignments
        for (MultipathAlignment& multipath_aln : multipath_alns_out) {
#ifdef debug_multipath_mapper_mapping
            cerr << "validating multipath alignment:" << endl;
            cerr << pb2json(multipath_aln) << endl;
#endif
            if (!validate_multipath_alignment(multipath_aln)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln.name() << " failed to validate" << endl;
            }
        }
#endif
        
    }
    
    void MultipathMapper::attempt_unpaired_multipath_map_of_pair(const Alignment& alignment1, const Alignment& alignment2,
                                                                 vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                                 vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer) {
        
        // compute single ended mappings, and make sure we also compute mapping qualities to assess
        // mapping ambiguity
        vector<MultipathAlignment> multipath_alns_1, multipath_alns_2;
        multipath_map_internal(alignment1, mapping_quality_method == None ? Approx : mapping_quality_method,
                               multipath_alns_1, 1);
        multipath_map_internal(alignment2, mapping_quality_method == None ? Approx : mapping_quality_method,
                               multipath_alns_2, 1);
        
        bool is_ambiguous = true;
        
        if (!multipath_alns_1.empty() && !multipath_alns_2.empty()) {
            MultipathAlignment& multipath_aln_1 = multipath_alns_1.front();
            MultipathAlignment& multipath_aln_2 = multipath_alns_2.front();
            
            auto match_score = get_aligner()->match;
            auto full_length_bonus = get_aligner()->full_length_bonus;
            
            // score possible of a perfect match (at full base quality)
            int32_t max_score_1 = multipath_aln_1.sequence().size() * match_score + 2 * full_length_bonus * !strip_bonuses;
            int32_t max_score_2 = multipath_aln_2.sequence().size() * match_score + 2 * full_length_bonus * !strip_bonuses;
            
#ifdef debug_multipath_mapper_mapping
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
                
#ifdef debug_multipath_mapper_mapping
                cerr << "fragment length between mappings measured at " << fragment_length << endl;
#endif
                
                // can we obtain a distance between these positions?
                if (fragment_length != numeric_limits<int64_t>::max()) {
                    
                    // record the unambiguous mappings and the fragment length
                    
                    
#ifdef debug_multipath_mapper_mapping
                    cerr << "registering measurement, now have " << fragment_length_distr.curr_sample_size() << " of " << fragment_length_distr.max_sample_size() << endl;
#endif
                    
                    multipath_aln_pairs_out.emplace_back(move(multipath_aln_1), move(multipath_aln_2));
                    multipath_aln_pairs_out.front().first.set_paired_read_name(multipath_aln_pairs_out.front().second.name());
                    multipath_aln_pairs_out.front().second.set_paired_read_name(multipath_aln_pairs_out.front().first.name());
                    
                    fragment_length_distr.register_fragment_length(fragment_length);
                    
                    is_ambiguous = false;
                }
            }
        }
        
        if (is_ambiguous) {
            // we didn't find an unambiguous pairing in single-ended mode, buffer these for once
            // the paired mode is finalized
#ifdef debug_multipath_mapper_mapping
            cerr << "couldn't find unambiguous mapping, adding pair to ambiguous buffer" << endl;
#endif
            
            ambiguous_pair_buffer.emplace_back(alignment1, alignment2);
        }
        
        
        // for debugging:
        // we must have just finalized the distribution or else we wouldn't have entered this function
#ifdef debug_report_startup_training
        if (fragment_length_distr.is_finalized()) {
            cerr << "finalized read distribution with " << fragment_length_distr.max_sample_size() << " measurements on read pair " << alignment1.name() << ", " << alignment2.name() << endl;
            cerr << "mean: " << fragment_length_distr.mean() << endl;
            cerr << "std dev: " << fragment_length_distr.stdev() << endl;
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
    
    bool MultipathMapper::attempt_rescue(const MultipathAlignment& multipath_aln, const Alignment& other_aln,
                                         bool rescue_forward, MultipathAlignment& rescue_multipath_aln) {
        
#ifdef debug_multipath_mapper_mapping
        cerr << "attemping pair rescue in " << (rescue_forward ? "forward" : "backward") << " direction from " << pb2json(multipath_aln) << endl;
#endif
        
        // get the position to jump from and the distance to jump
        Alignment opt_anchoring_aln;
        optimal_alignment(multipath_aln, opt_anchoring_aln);
        pos_t pos_from = rescue_forward ? initial_position(opt_anchoring_aln.path()) : final_position(opt_anchoring_aln.path());
        int64_t jump_dist = rescue_forward ? fragment_length_distr.mean() - other_aln.sequence().size() : -fragment_length_distr.mean();
        
        // get the seed position(s) for the rescue by jumping along paths
        vector<pos_t> jump_positions = xindex->jump_along_closest_path(id(pos_from), is_rev(pos_from), offset(pos_from), jump_dist, 250);
        
#ifdef debug_multipath_mapper_mapping
        cerr << "found jump positions:" << endl;
        for (pos_t& pos : jump_positions) {
            cerr << "\t" << pos << endl;
        }
#endif
        if (jump_positions.empty()) {
            return false;
        }
        
        // pull out the graph around the position(s) we jumped to
        VG rescue_graph;
        vector<size_t> backward_dist(jump_positions.size(), 6 * fragment_length_distr.stdev());
        vector<size_t> forward_dist(jump_positions.size(), 6 * fragment_length_distr.stdev() + other_aln.sequence().size());
        algorithms::extract_containing_graph(xindex, rescue_graph.graph, jump_positions, backward_dist, forward_dist);
        rescue_graph.build_indexes();
        
#ifdef debug_multipath_mapper_mapping
        cerr << "got rescue graph " << pb2json(rescue_graph.graph) << endl;
#endif
        
        
        // TODO: repetitive code with multipath_align
        
        // the longest path we could possibly align to (full gap and a full sequence)
        size_t target_length = other_aln.sequence().size() + get_aligner()->longest_detectable_gap(other_aln);
        
        // convert from bidirected to directed
        unordered_map<id_t, pair<id_t, bool> > node_trans;
        VG align_graph = rescue_graph.split_strands(node_trans);
        // if necessary, convert from cyclic to acylic
        if (!rescue_graph.is_directed_acyclic()) {
            unordered_map<id_t, pair<id_t, bool> > dagify_trans;
            align_graph = align_graph.dagify(target_length, // high enough that num SCCs is never a limiting factor
                                             dagify_trans,
                                             target_length,
                                             0); // no maximum on size of component
            node_trans = align_graph.overlay_node_translations(dagify_trans, node_trans);
        }
        
        // put local alignment here
        Alignment aln = other_aln;
        // in case we're realigning a GAM, get rid of the path
        aln.clear_path();
        
        align_graph.lazy_sort();
        
        get_aligner()->align(aln, align_graph.graph, true, false);
        translate_oriented_node_ids(*aln.mutable_path(), node_trans);
        
        to_multipath_alignment(aln, rescue_multipath_aln);
        identify_start_subpaths(rescue_multipath_aln);
        
        vector<double> score(1, aln.score());
        int32_t raw_mapq = get_aligner()->compute_mapping_quality(score, mapping_quality_method == None || mapping_quality_method == Approx);
        int32_t adjusted_mapq = min(raw_mapq, min(max_mapping_quality, multipath_aln.mapping_quality()));
        rescue_multipath_aln.set_mapping_quality(adjusted_mapq);
        
#ifdef debug_multipath_mapper_mapping
        cerr << "rescued alignment is " << pb2json(rescue_multipath_aln) << endl;
        cerr << "rescued alignment has effective match length " << pseudo_length(rescue_multipath_aln) / 3 << ", which gives p-value " << random_match_p_value(pseudo_length(rescue_multipath_aln) / 3, rescue_multipath_aln.sequence().size()) << endl;
#endif
        return (raw_mapq >= min(25, max_mapping_quality)
                && random_match_p_value(
                    pseudo_length(rescue_multipath_aln) / 3,
                    rescue_multipath_aln.sequence().size()) < max_mapping_p_value * 0.1);
    }
    
    bool MultipathMapper::likely_mismapping(const MultipathAlignment& multipath_aln) {
#ifdef debug_multipath_mapper_mapping
        cerr << "effective match length of read " << multipath_aln.name() << " is " << pseudo_length(multipath_aln) / 3 << " in read length " << multipath_aln.sequence().size() << ", yielding p-value " << random_match_p_value(pseudo_length(multipath_aln) / 3, multipath_aln.sequence().size()) << endl;
#endif
        
        // empirically, we get better results by scaling the pseudo-length down, I have no good explanation for this probabilistically
        return random_match_p_value(pseudo_length(multipath_aln) / 3, multipath_aln.sequence().size()) > max_mapping_p_value;
    }
    
    size_t MultipathMapper::pseudo_length(const MultipathAlignment& multipath_aln) const {
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
            double p_value = 1.0 - pow(1.0 - exp(-(match_length * pseudo_length_multiplier)), xindex->seq_length * read_length);
            if (p_value_memo.size() < max_p_value_memo_size) {
                p_value_memo[make_pair(match_length, read_length)] = p_value;
            }
            return p_value;
        }
    }
    
    void MultipathMapper::calibrate_mismapping_detection(size_t num_simulations, size_t simulated_read_length) {
        // we don't want to do base quality adjusted alignments for this stage since we are just simulating random sequences
        // with no base qualities
        bool reset_quality_adjustments = adjust_alignments_for_base_quality;
        adjust_alignments_for_base_quality = false;
        
        // compute the pseudo length of a bunch of randomly generated sequences
        vector<double> lengths(num_simulations, 0.0);
        double length_sum = 0.0;
        double max_length = numeric_limits<double>::min();
#pragma omp parallel for
        for (size_t i = 0; i < num_simulations; i++) {
            
            Alignment alignment;
            alignment.set_sequence(random_sequence(simulated_read_length));
            vector<MultipathAlignment> multipath_alns;
            multipath_map(alignment, multipath_alns, 1);
            
            if (!multipath_alns.empty()) {
                lengths[i] = pseudo_length(multipath_alns.front());
#pragma omp critical
                {
                    length_sum += lengths[i];
                    max_length = max(lengths[i], max_length);
                }
            }
        }
        
        // reset the memo of p-values (which we are calibrating) for any updates using the default parameter during the null mappings
        p_value_memo.clear();
        
        // model the lengths as the maximum of genome_size * read_length exponential variables, which gives density function:
        //
        //   GLS exp(-Sx) (1 - exp(-Sx))^(GL - 1)
        //
        // where:
        //   G = genome size
        //   R = read length
        //   S = scale parameter (which we optimize below)
        
        
        // compute the log of the 1st and 2nd derivatives for the log likelihood (split up by positive and negative summands)
        // we have to do it this wonky way because the exponentiated numbers get very large and cause overflow otherwise
        
        double log_deriv_neg_part = log(length_sum);
        
        function<double(double)> log_deriv_pos_part = [&](double scale) {
            double accumulator = numeric_limits<double>::lowest();
            for (size_t i = 0; i < lengths.size(); i++) {
                double length = lengths[i];
                accumulator = add_log(accumulator, log(length) - scale * length - log(1.0 - exp(-scale * length)));
            }
            accumulator += log(xindex->seq_length * simulated_read_length - 1.0);
            return add_log(accumulator, log(num_simulations / scale));
        };
        
        function<double(double)> log_deriv2_neg_part = [&](double scale) {
            double accumulator = numeric_limits<double>::lowest();
            for (size_t i = 0; i < lengths.size(); i++) {
                double length = lengths[i];
                accumulator = add_log(accumulator, 2.0 * log(length) - scale * length - 2.0 * log(1.0 - exp(-scale * length)));
            }
            accumulator += log(xindex->seq_length * simulated_read_length - 1.0);
            return add_log(accumulator, log(num_simulations / (scale * scale)));
        };
        
        // use Newton's method to find the MLE
        double tolerance = 1e-10;
        double scale = 1.0 / max_length;
        double prev_scale = scale * (1.0 + 10.0 * tolerance);
        while (abs(prev_scale / scale - 1.0) > tolerance) {
            prev_scale = scale;
            double log_d2 = log_deriv2_neg_part(scale);
            double log_d_pos = log_deriv_pos_part(scale);
            double log_d_neg = log_deriv_neg_part;
            // determine if the value of the 1st deriv is positive or negative, and compute the
            // whole ratio to the 2nd deriv from the positive and negative parts accordingly
            if (log_d_pos > log_d_neg) {
                scale += exp(subtract_log(log_d_pos, log_d_neg) - log_d2);
            }
            else {
                scale -= exp(subtract_log(log_d_neg, log_d_pos) - log_d2);
            }
        }
        
#ifdef debug_report_startup_training
        cerr << "trained scale: " << scale << endl;
#endif
        
        // set the multipler to the maximimum likelihood
        pseudo_length_multiplier = scale;
        
        adjust_alignments_for_base_quality = reset_quality_adjustments;
    }
    
    int64_t MultipathMapper::distance_between(const MultipathAlignment& multipath_aln_1,
                                              const MultipathAlignment& multipath_aln_2,
                                              bool full_fragment, bool forward_strand) const {
        Alignment aln_1;
        optimal_alignment(multipath_aln_1, aln_1);
        pos_t pos_1 = initial_position(aln_1.path());
        
        Alignment aln_2;
        optimal_alignment(multipath_aln_2, aln_2);
        pos_t pos_2 = full_fragment ? final_position(aln_2.path()) : initial_position(aln_2.path());
#ifdef debug_multipath_mapper_mapping
        cerr << "measuring left-to-" << (full_fragment ? "right" : "left") << " end distance between " << pos_1 << " and " << pos_2 << endl;
#endif
        return xindex->closest_shared_path_oriented_distance(id(pos_1), offset(pos_1), is_rev(pos_1),
                                                             id(pos_2), offset(pos_2), is_rev(pos_2),
                                                             forward_strand);
    }
    
    bool MultipathMapper::is_consistent(int64_t distance) const {
        return (distance < fragment_length_distr.mean() + 10.0 * fragment_length_distr.stdev()
                && distance > fragment_length_distr.mean() - 10.0 * fragment_length_distr.stdev());
    }
    
    bool MultipathMapper::are_consistent(const MultipathAlignment& multipath_aln_1,
                                         const MultipathAlignment& multipath_aln_2) const {
        
        return is_consistent(distance_between(multipath_aln_1, multipath_aln_2, true));
    }
    
    bool MultipathMapper::share_start_position(const MultipathAlignment& multipath_aln_1,
                                               const MultipathAlignment& multipath_aln_2) const {
        
        unordered_set<pos_t> starts;
        for (size_t i = 0; i < multipath_aln_1.start_size(); i++) {
            starts.insert(make_pos_t(multipath_aln_1.subpath(multipath_aln_1.start(i)).path().mapping(0).position()));
        }
        
        for (size_t i = 0; i < multipath_aln_2.start_size(); i++) {
            if (starts.count(make_pos_t(multipath_aln_2.subpath(multipath_aln_2.start(i)).path().mapping(0).position()))) {
                return true;
            }
        }
        return false;
    }
    
    void MultipathMapper::establish_strand_consistency(vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs,
                                                       vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                       OrientedDistanceClusterer::paths_of_node_memo_t* paths_of_node_memo,
                                                       OrientedDistanceClusterer::oriented_occurences_memo_t* oriented_occurences_memo,
                                                       OrientedDistanceClusterer::handle_memo_t* handle_memo) {
        
#ifdef debug_multipath_mapper_mapping
        cerr << "establishing consistency between mapped pairs" << endl;
#endif
        
        int64_t search_dist = 0.5 * fragment_length_distr.mean() + 5.0 * fragment_length_distr.stdev();
        vector<pair<bool, bool>> strand_assignments;
        strand_assignments.reserve(multipath_aln_pairs.size());
        for (const pair<MultipathAlignment, MultipathAlignment>& multipath_aln_pair : multipath_aln_pairs) {
            Alignment optimal_aln_1, optimal_aln_2;
            optimal_alignment(multipath_aln_pair.first, optimal_aln_1);
            optimal_alignment(multipath_aln_pair.second, optimal_aln_2);
            pos_t pos_1 = initial_position(optimal_aln_1.path());
            pos_t pos_2 = initial_position(optimal_aln_2.path());
            
            strand_assignments.push_back(xindex->validate_strand_consistency(id(pos_1), offset(pos_1), is_rev(pos_1),
                                                                             id(pos_2), offset(pos_2), is_rev(pos_2),
                                                                             search_dist, paths_of_node_memo,
                                                                             oriented_occurences_memo, handle_memo));
            
#ifdef debug_multipath_mapper_mapping
            cerr << "pair has initial positions " << pos_1 << " and " << pos_2 << " on strands " << (strand_assignments.back().first ? "-" : "+") << " and " << (strand_assignments.back().second ? "-" : "+") << endl;
#endif
        }
        
        size_t end = multipath_aln_pairs.size();
        for (size_t i = 0; i < end; ) {
            // move strand inconsistent mappings to the end
            if (strand_assignments[i].first != strand_assignments[i].second) {
#ifdef debug_multipath_mapper_mapping
                cerr << "removing inconsistent strand at " << i << " and not advancing index" << endl;
#endif                
                std::swap(multipath_aln_pairs[i], multipath_aln_pairs[end - 1]);
                std::swap(strand_assignments[i], strand_assignments[end - 1]);
                std::swap(cluster_pairs[i], cluster_pairs[end - 1]);
                
                --end;
            }
            else {
#ifdef debug_multipath_mapper_mapping
                cerr << "identifying " << i << " as consistent" << endl;
#endif
                // reverse the distance if it's on the reverse strand
                if (strand_assignments[i].first) {
#ifdef debug_multipath_mapper_mapping
                    cerr << "\tinverting distance " << cluster_pairs[i].second << " because on negative strand" << endl;
#endif
                    cluster_pairs[i].second = -cluster_pairs[i].second;
                }
                ++i;
            }
        }
        
        // remove the inconsistent mappings
        if (end != multipath_aln_pairs.size()) {
#ifdef debug_multipath_mapper_mapping
            cerr << "found " << multipath_aln_pairs.size() - end << " strand inconsitent pairs, removing now" << endl;
#endif
            multipath_aln_pairs.resize(end);
            cluster_pairs.resize(end);
        }
    }
    
    bool MultipathMapper::align_to_cluster_graphs_with_rescue(const Alignment& alignment1, const Alignment& alignment2,
                                                              vector<clustergraph_t>& cluster_graphs1,
                                                              vector<clustergraph_t>& cluster_graphs2,
                                                              bool block_rescue_from_1, bool block_rescue_from_2,
                                                              vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                              vector<pair<pair<size_t, size_t>, int64_t>>& pair_distances,
                                                              size_t max_alt_mappings) {
        
        vector<MultipathAlignment> multipath_alns_1, multipath_alns_2;
        if (!block_rescue_from_1) {
            align_to_cluster_graphs(alignment1, mapping_quality_method == None ? Approx : mapping_quality_method,
                                    cluster_graphs1, multipath_alns_1, max_alt_mappings);
        }
        if (!block_rescue_from_2) {
            align_to_cluster_graphs(alignment2, mapping_quality_method == None ? Approx : mapping_quality_method,
                                    cluster_graphs2, multipath_alns_2, max_alt_mappings);
        }
        
        if (multipath_alns_1.empty() || multipath_alns_2.empty() ? false :
            ((multipath_alns_1.front().mapping_quality() >= min(60, max_mapping_quality)
              && multipath_alns_2.front().mapping_quality() >= min(60, max_mapping_quality)) ?
             are_consistent(multipath_alns_1.front(), multipath_alns_2.front()) : false)) {
            
            // we are able to obtain confident matches that satisfy the pairing constraints
#ifdef debug_multipath_mapper_mapping
            cerr << "found consistent, confident pair mapping from independent end mapping" << endl;
#endif
            multipath_aln_pairs_out.emplace_back(move(multipath_alns_1.front()), move(multipath_alns_2.front()));
            pair_distances.emplace_back(pair<size_t, size_t>(),
                                        distance_between(multipath_aln_pairs_out.back().first, multipath_aln_pairs_out.back().second, true));
            return true;
        }
        
        int32_t max_score_diff = get_aligner()->mapping_quality_score_diff(max_mapping_quality);
        
        int32_t top_score_1 = multipath_alns_1.empty() ? 0 : optimal_alignment_score(multipath_alns_1.front());
        int32_t top_score_2 = multipath_alns_2.empty() ? 0 : optimal_alignment_score(multipath_alns_2.front());
        
        size_t num_rescuable_alns_1 = block_rescue_from_1 ? 0 : min(multipath_alns_1.size(), max_rescue_attempts);
        size_t num_rescuable_alns_2 = block_rescue_from_2 ? 0 : min(multipath_alns_2.size(), max_rescue_attempts);
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
        
        vector<MultipathAlignment> rescue_multipath_alns_1(num_rescuable_alns_2), rescue_multipath_alns_2(num_rescuable_alns_1);
        
        unordered_set<size_t> rescued_from_1, rescued_from_2;
        
#ifdef debug_multipath_mapper_mapping
        cerr << "rescuing from " << num_rescuable_alns_1 << " read1's and " << num_rescuable_alns_2 << " read2's" << endl;
#endif
        
        for (size_t i = 0; i < num_rescuable_alns_1; i++) {
            MultipathAlignment rescue_multipath_aln;
            if (attempt_rescue(multipath_alns_1[i], alignment2, true, rescue_multipath_aln)) {
                rescued_from_1.insert(i);
                rescue_multipath_alns_2[i] = move(rescue_multipath_aln);
            }
        }
        
        for (size_t i = 0; i < num_rescuable_alns_2; i++) {
            MultipathAlignment rescue_multipath_aln;
            if (attempt_rescue(multipath_alns_2[i], alignment1, false, rescue_multipath_aln)) {
                rescued_from_2.insert(i);
                rescue_multipath_alns_1[i] = move(rescue_multipath_aln);
            }
        }
        
        bool found_consistent = false;
        if (!rescued_from_1.empty() && !rescued_from_2.empty()) {
#ifdef debug_multipath_mapper_mapping
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
                    
#ifdef debug_multipath_mapper_mapping
                    cerr << "checking duplication between mapped read1 " << i << " and rescued read1 " << j << endl;
#endif
                    if (abs(distance_between(multipath_alns_1[i], rescue_multipath_alns_1[j])) < 20) {
#ifdef debug_multipath_mapper_mapping
                        cerr << "found duplicate, now checking rescued read2 " << i << " and mapped read2 " << j << endl;
#endif
                        if (abs(distance_between(rescue_multipath_alns_2[i], multipath_alns_2[j])) < 20) {
#ifdef debug_multipath_mapper_mapping
                            cerr << "found duplicate, marking entire pair as duplicate" << endl;
#endif
                            // these two alignments found each other with their rescue, we don't want to add duplicate mappings
                            duplicate = true;
                            found_duplicate.insert(j);
                            
                            // move the original mappings
                            int64_t dist = distance_between(multipath_alns_1[i], multipath_alns_2[j], true);
                            if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
                                multipath_aln_pairs_out.emplace_back(move(multipath_alns_1[i]), move(multipath_alns_2[j]));
                                pair_distances.emplace_back(pair<size_t, size_t>(), dist);
                                found_consistent = true;
                            }
                            
                            break;
                        }
                    }
                }
                
                // if we haven't already moved the pair and marked it as a duplicate, move the rescued pair into the output vector
                if (!duplicate) {
#ifdef debug_multipath_mapper_mapping
                    cerr << "adding read1 and rescued read2 " << i << " to output vector" << endl;
#endif
                    int64_t dist = distance_between(multipath_alns_1[i], rescue_multipath_alns_2[i], true);
                    if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
                        multipath_aln_pairs_out.emplace_back(move(multipath_alns_1[i]), move(rescue_multipath_alns_2[i]));
                        pair_distances.emplace_back(pair<size_t, size_t>(), dist);
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
#ifdef debug_multipath_mapper_mapping
                cerr << "adding rescued read1 and read2 " << j << " to output vector" << endl;
#endif
                int64_t dist = distance_between(rescue_multipath_alns_1[j], multipath_alns_2[j], true);
                if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
                    multipath_aln_pairs_out.emplace_back(move(rescue_multipath_alns_1[j]), move(multipath_alns_2[j]));
                    pair_distances.emplace_back(pair<size_t, size_t>(), dist);
                    found_consistent = true;
                }
            }
        }
        else if (!rescued_from_1.empty()) {
#ifdef debug_multipath_mapper_mapping
            cerr << "successfully rescued from only read 1" << endl;
#endif
            for (size_t i: rescued_from_1) {
                int64_t dist = distance_between(multipath_alns_1[i], rescue_multipath_alns_2[i], true);
                if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
                    multipath_aln_pairs_out.emplace_back(move(multipath_alns_1[i]), move(rescue_multipath_alns_2[i]));
                    pair_distances.emplace_back(pair<size_t, size_t>(), dist);
                    found_consistent = true;
                }
            }
        }
        else if (!rescued_from_2.empty()) {
#ifdef debug_multipath_mapper_mapping
            cerr << "successfully rescued from only read 2" << endl;
#endif
            for (size_t i : rescued_from_2) {
                int64_t dist = distance_between(rescue_multipath_alns_1[i], multipath_alns_2[i], true);
                if (dist != numeric_limits<int64_t>::max() && dist >= 0) {
                    multipath_aln_pairs_out.emplace_back(move(rescue_multipath_alns_1[i]), move(multipath_alns_2[i]));
                    pair_distances.emplace_back(pair<size_t, size_t>(), dist);
                    found_consistent = true;
                }
            }
        }
        
        if (found_consistent) {
            // compute the paired mapping quality
            sort_and_compute_mapping_quality(multipath_aln_pairs_out, pair_distances);
        }
        else {
#ifdef debug_multipath_mapper_mapping
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
                    multipath_aln_pairs_out.emplace_back(move(multipath_alns_1[i]), MultipathAlignment());
                    to_multipath_alignment(alignment2, multipath_aln_pairs_out.back().second);
                    multipath_aln_pairs_out.back().second.clear_subpath();
                    multipath_aln_pairs_out.back().second.clear_start();
                }
                else {
                    multipath_aln_pairs_out.emplace_back(MultipathAlignment(), move(multipath_alns_2[i]));
                    to_multipath_alignment(alignment1, multipath_aln_pairs_out.back().first);
                    multipath_aln_pairs_out.back().first.clear_subpath();
                    multipath_aln_pairs_out.back().first.clear_start();
                }
            }
        }
        
#ifdef debug_validate_multipath_alignments
        for (pair<MultipathAlignment, MultipathAlignment>& multipath_aln_pair : multipath_aln_pairs_out) {
#ifdef debug_multipath_mapper_mapping
            cerr << "validating multipath alignments:" << endl;
            cerr << pb2json(multipath_aln_pair.first) << endl;
            cerr << pb2json(multipath_aln_pair.second) << endl;
#endif
            if (!validate_multipath_alignment(multipath_aln_pair.first)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.first.name() << " failed to validate" << endl;
            }
            if (!validate_multipath_alignment(multipath_aln_pair.second)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.second.name() << " failed to validate" << endl;
            }
        }
#endif
        
        if (multipath_aln_pairs_out.size() > max_alt_mappings) {
            multipath_aln_pairs_out.resize(max_alt_mappings);
        }
        
        if (mapping_quality_method == None) {
            for (pair<MultipathAlignment, MultipathAlignment>& multipath_aln_pair : multipath_aln_pairs_out) {
                multipath_aln_pair.first.clear_mapping_quality();
                multipath_aln_pair.second.clear_mapping_quality();
            }
        }
        
        return found_consistent;
    }
    
    void MultipathMapper::attempt_rescue_for_secondaries(const Alignment& alignment1, const Alignment& alignment2,
                                                         vector<clustergraph_t>& cluster_graphs1,
                                                         vector<clustergraph_t>& cluster_graphs2,
                                                         vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                         vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs) {
        
#ifdef debug_multipath_mapper_mapping
        cerr << "using rescue to find secondary mappings" << endl;
#endif

        
        unordered_set<size_t> paired_clusters_1, paired_clusters_2;
        
        for (size_t i = 0; i < multipath_aln_pairs_out.size(); i++) {
            // keep track of which clusters already have consistent pairs
            paired_clusters_1.insert(cluster_pairs[i].first.first);
            paired_clusters_2.insert(cluster_pairs[i].first.second);
        }
        
        int32_t cluster_score_1 = get_aligner()->match * get<2>(cluster_graphs1[cluster_pairs.front().first.first]);
        int32_t cluster_score_2 = get_aligner()->match * get<2>(cluster_graphs2[cluster_pairs.front().first.second]);
        int32_t max_score_diff = secondary_rescue_score_diff * get_aligner()->mapping_quality_score_diff(max_mapping_quality);
        
        vector<pair<MultipathAlignment, MultipathAlignment>> rescued_secondaries;
        vector<pair<pair<size_t, size_t>, int64_t>> rescued_distances;
        
        auto align_and_rescue = [&](const Alignment& anchor_aln, const Alignment& rescue_aln,
                                    vector<clustergraph_t>& cluster_graphs, unordered_set<size_t>& paired_clusters,
                                    int32_t max_score, bool anchor_is_read_1) {
            
#ifdef debug_multipath_mapper_mapping
            cerr << "checking for rescues from read " << (anchor_is_read_1 ? 1 : 2) << endl;
#endif

            
            size_t num_rescues = 0;
            for (size_t i = 0; i < cluster_graphs.size() && num_rescues < secondary_rescue_attempts; i++) {
                if (paired_clusters.count(i)) {
                    // we already have a consistent pair from this cluster
#ifdef debug_multipath_mapper_mapping
                    cerr << "cluster " << i << " is already in a pair" << endl;
#endif

                    continue;
                }
                
#ifdef debug_multipath_mapper_mapping
                cerr << "cluster " << i << "'s approximate score is " << get<2>(cluster_graphs[i]) * get_aligner()->match << ", looking for " << max_score - max_score_diff << endl;
#endif
                
                if (get<2>(cluster_graphs[i]) * get_aligner()->match < max_score - max_score_diff) {
                    // the approximate score of the remaining is too low to consider
                    break;
                }
                
                // TODO: repetitive with align_to_cluster_graphs
                
                // make the alignment
                vector<MultipathAlignment> cluster_multipath_alns;
                cluster_multipath_alns.emplace_back();
                multipath_align(anchor_aln, get<0>(cluster_graphs[i]), get<1>(cluster_graphs[i]), cluster_multipath_alns.back());
                
                // split it up if it turns out to be multiple components
                split_multicomponent_alignments(cluster_multipath_alns);
                
                // order the subpaths
                for (MultipathAlignment& multipath_aln : cluster_multipath_alns) {
                    topologically_order_subpaths(multipath_aln);
                }
                
                // if we split it up, move the best one to the front
                if (cluster_multipath_alns.size() > 1) {
                    sort_and_compute_mapping_quality(cluster_multipath_alns, None);
                }
                
                // rescue from the alignment
                MultipathAlignment rescue_multipath_aln;
                if (!likely_mismapping(cluster_multipath_alns.front())) {
                    bool rescued = attempt_rescue(cluster_multipath_alns.front(), rescue_aln, anchor_is_read_1, rescue_multipath_aln);
#ifdef debug_multipath_mapper_mapping
                    cerr << "rescued alignment is " << pb2json(rescue_multipath_aln) << endl;
#endif
                    if (rescued) {
#ifdef debug_multipath_mapper_mapping
                        cerr << "rescue succeeded, adding to rescue pair vector" << endl;
#endif
                        if (anchor_is_read_1) {
                            int64_t dist = distance_between(cluster_multipath_alns.front(), rescue_multipath_aln, true);
                            if (dist >= 0 && dist != numeric_limits<int64_t>::max()) {
                                rescued_secondaries.emplace_back(move(cluster_multipath_alns.front()), move(rescue_multipath_aln));
                                rescued_distances.emplace_back(pair<size_t, size_t>(), dist);
                                
                            }
                        }
                        else {
                            int64_t dist = distance_between(rescue_multipath_aln, cluster_multipath_alns.front(), true);
                            if (dist >= 0 && dist != numeric_limits<int64_t>::max()) {
                                rescued_secondaries.emplace_back(move(rescue_multipath_aln), move(cluster_multipath_alns.front()));
                                rescued_distances.emplace_back(pair<size_t, size_t>(), dist);
                                
                            }
                        }
                    }
                }
                
                num_rescues++;
            }
        };
        
        // perform routine for both read ends
        align_and_rescue(alignment1, alignment2, cluster_graphs1, paired_clusters_1, cluster_score_1, true);
        align_and_rescue(alignment2, alignment1, cluster_graphs2, paired_clusters_2, cluster_score_2, false);
        
        if (!rescued_secondaries.empty()) {
            // we found mappings that could be rescues of each other
            
#ifdef debug_multipath_mapper_mapping
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
            }
            
            // merge the rescued secondaries into the return vector
            merge_rescued_mappings(multipath_aln_pairs_out, cluster_pairs, rescued_secondaries, rescued_distances);
        }
    }
    
    void MultipathMapper::multipath_map_paired(const Alignment& alignment1, const Alignment& alignment2,
                                               vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                               vector<pair<Alignment, Alignment>>& ambiguous_pair_buffer,
                                               size_t max_alt_mappings) {
        
#ifdef debug_multipath_mapper_mapping
        cerr << "multipath mapping paired reads " << pb2json(alignment1) << " and " << pb2json(alignment2) << endl;
#endif
        
        // empty the output vector (just for safety)
        multipath_aln_pairs_out.clear();
        
        if (!fragment_length_distr.is_finalized()) {
            // we have not estimated a fragment length distribution yet, so we revert to single ended mode and look
            // for unambiguous pairings
            
#ifdef debug_multipath_mapper_mapping
            cerr << "no fragment length distribution yet, looking for unambiguous single ended pairs" << endl;
#endif
            
            attempt_unpaired_multipath_map_of_pair(alignment1, alignment2, multipath_aln_pairs_out, ambiguous_pair_buffer);
            
            return;
        }
        
        // the fragment length distribution has been estimated, so we can do full-fledged paired mode
    
        // query MEMs using GCSA2
        double dummy1, dummy2;
        vector<MaximalExactMatch> mems1 = find_mems_deep(alignment1.sequence().begin(), alignment1.sequence().end(), dummy1, dummy2,
                                                         0, min_mem_length, mem_reseed_length, false, true, true, false);
        vector<MaximalExactMatch> mems2 = find_mems_deep(alignment2.sequence().begin(), alignment2.sequence().end(), dummy1, dummy2,
                                                         0, min_mem_length, mem_reseed_length, false, true, true, false);
        
#ifdef debug_multipath_mapper_mapping
        cerr << "obtained read1 MEMs:" << endl;
        for (MaximalExactMatch mem : mems1) {
            cerr << "\t" << mem << " (" << mem.nodes.size() << " hits filled out of " << mem.match_count << ")" << endl;
        }
        cerr << "obtained read2 MEMs:" << endl;
        for (MaximalExactMatch mem : mems2) {
            cerr << "\t" << mem << " (" << mem.nodes.size() << " hits filled out of " << mem.match_count << ")" << endl;
        }
#endif
        
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
        
        // intialize memos for the results of expensive succinct operations that we may need to do multiple times
        OrientedDistanceClusterer::paths_of_node_memo_t paths_of_node_memo;
        OrientedDistanceClusterer::oriented_occurences_memo_t oriented_occurences_memo;
        OrientedDistanceClusterer::handle_memo_t handle_memo;
        
        // do we want to try to only cluster one read end and rescue the other?
        bool do_repeat_rescue_from_1 = min_match_count_2 > rescue_only_min && min_match_count_1 <= rescue_only_anchor_max;
        bool do_repeat_rescue_from_2 = min_match_count_1 > rescue_only_min && min_match_count_2 <= rescue_only_anchor_max;
        
        bool rescued_order_length_runs_1 = false, rescued_order_length_runs_2 = false;
        
#ifdef debug_multipath_mapper_mapping
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
                                                     cluster_pairs, max_alt_mappings, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
            
            if (multipath_aln_pairs_out.empty() && do_repeat_rescue_from_1 && !do_repeat_rescue_from_2) {
                // we've clustered and extracted read 1, but rescue failed, so do the same for read 2 to prepare for the
                // normal pair clustering routine
                
#ifdef debug_multipath_mapper_mapping
                cerr << "repeat rescue failed from read 1, extracting clusters for read 2 and transitioning to standard clustering approach" << endl;
#endif
                
                // rescue high count runs of order-length MEMs now that we're going to cluster here
                if (order_length_repeat_hit_max && !rescued_order_length_runs_2) {
                    rescue_high_count_order_length_mems(mems2, order_length_repeat_hit_max);
                    rescued_order_length_runs_2 = true;
                }
                
                // do the clustering
                if (adjust_alignments_for_base_quality) {
                    OrientedDistanceClusterer clusterer2(alignment2, mems2, *get_qual_adj_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                         unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
                    clusters2 = clusterer2.clusters(max_mapping_quality, log_likelihood_approx_factor);
                }
                else {
                    OrientedDistanceClusterer clusterer2(alignment2, mems2, *get_regular_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                         unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
                    clusters2 = clusterer2.clusters(max_mapping_quality, log_likelihood_approx_factor);
                }
                
                cluster_graphs2 = query_cluster_graphs(alignment2, mems2, clusters2);
            }
            
            if (multipath_aln_pairs_out.empty() && do_repeat_rescue_from_2 && !do_repeat_rescue_from_1) {
                // we've clustered and extracted read 2, but rescue failed, so do the same for read 1 to prepare for the
                // normal pair clustering routine
                
#ifdef debug_multipath_mapper_mapping
                cerr << "repeat rescue failed from read 2, extracting clusters for read 1 and transitioning to standard clustering approach" << endl;
#endif
                
                // rescue high count runs of order-length MEMs now that we're going to cluster here
                if (order_length_repeat_hit_max && !rescued_order_length_runs_1) {
                    rescue_high_count_order_length_mems(mems1, order_length_repeat_hit_max);
                    rescued_order_length_runs_1 = true;
                }
                
                // do the clustering
                if (adjust_alignments_for_base_quality) {
                    OrientedDistanceClusterer clusterer1(alignment1, mems1, *get_qual_adj_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                         unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
                    clusters1 = clusterer1.clusters(max_mapping_quality, log_likelihood_approx_factor);
                }
                else {
                    OrientedDistanceClusterer clusterer1(alignment1, mems1, *get_regular_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                         unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
                    clusters1 = clusterer1.clusters(max_mapping_quality, log_likelihood_approx_factor);
                }
                
                cluster_graphs1 = query_cluster_graphs(alignment1, mems1, clusters1);
            }
        }
        else {
            // we have reasonably unique hits on both reads, so cluster them both and extract the graphs (i.e. don't try
            // to rely on rescue for either end yet)
            
#ifdef debug_multipath_mapper_mapping
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
            if (adjust_alignments_for_base_quality) {
                OrientedDistanceClusterer clusterer1(alignment1, mems1, *get_qual_adj_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                     unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
                clusters1 = clusterer1.clusters(max_mapping_quality, log_likelihood_approx_factor);
                OrientedDistanceClusterer clusterer2(alignment2, mems2, *get_qual_adj_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                     unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
                clusters2 = clusterer2.clusters(max_mapping_quality, log_likelihood_approx_factor);
            }
            else {
                OrientedDistanceClusterer clusterer1(alignment1, mems1, *get_regular_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                     unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
                clusters1 = clusterer1.clusters(max_mapping_quality, log_likelihood_approx_factor);
                OrientedDistanceClusterer clusterer2(alignment2, mems2, *get_regular_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                     unstranded_clustering, &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
                clusters2 = clusterer2.clusters(max_mapping_quality, log_likelihood_approx_factor);
            }
            
            // extract graphs around the clusters and get the assignments of MEMs to these graphs
            cluster_graphs1 = query_cluster_graphs(alignment1, mems1, clusters1);
            cluster_graphs2 = query_cluster_graphs(alignment2, mems2, clusters2);
        }
        
#ifdef debug_multipath_mapper_mapping
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
            
            // make vectors of cluster pointers to shim into the cluster pairing function
            vector<memcluster_t*> cluster_mems_1(cluster_graphs1.size()), cluster_mems_2(cluster_graphs2.size());
            for (size_t i = 0; i < cluster_mems_1.size(); i++) {
                cluster_mems_1[i] = &(get<1>(cluster_graphs1[i]));
            }
            for (size_t i = 0; i < cluster_mems_2.size(); i++) {
                cluster_mems_2[i] = &(get<1>(cluster_graphs2[i]));
            }
            
            // Chebyshev bound for 99% of all fragments regardless of distribution
            // TODO: I don't love having this internal aspect of the stranded/unstranded clustering outside the clusterer...
            int64_t max_separation, min_separation;
            if (unstranded_clustering) {
                max_separation = (int64_t) ceil(abs(fragment_length_distr.mean()) + 10.0 * fragment_length_distr.stdev());
                min_separation = -max_separation;
            }
            else {
                max_separation = (int64_t) ceil(fragment_length_distr.mean() + 10.0 * fragment_length_distr.stdev());
                min_separation = (int64_t) fragment_length_distr.mean() - 10.0 * fragment_length_distr.stdev();
            }
            
            // Compute the pairs of cluster graphs and their approximate distances from each other
            cluster_pairs = OrientedDistanceClusterer::pair_clusters(alignment1, alignment2, cluster_mems_1, cluster_mems_2,
                                                                     xindex, min_separation, max_separation, unstranded_clustering,
                                                                     &paths_of_node_memo, &oriented_occurences_memo, &handle_memo);
#ifdef debug_multipath_mapper_mapping
            cerr << "obtained cluster pairs:" << endl;
            for (int i = 0; i < cluster_pairs.size(); i++) {
                pos_t pos_1 = get<1>(cluster_graphs1[cluster_pairs[i].first.first]).front().second;
                pos_t pos_2 = get<1>(cluster_graphs2[cluster_pairs[i].first.second]).back().second;
                int64_t dist = xindex->closest_shared_path_oriented_distance(id(pos_1), offset(pos_1), is_rev(pos_1),
                                                                             id(pos_2), offset(pos_2), is_rev(pos_2));
                cerr << "\tpair "  << i << " at distance " << dist << endl;
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
                // only perform the mappings that satisfy the expectations on distance
                
                align_to_cluster_graph_pairs(alignment1, alignment2, cluster_graphs1, cluster_graphs2, cluster_pairs,
                                             multipath_aln_pairs_out, max_alt_mappings, &paths_of_node_memo,
                                             &oriented_occurences_memo, &handle_memo);
                
                // do we produce at least one good looking pair alignments from the clustered clusters?
                if (multipath_aln_pairs_out.empty() ? true : (likely_mismapping(multipath_aln_pairs_out.front().first) ||
                                                              likely_mismapping(multipath_aln_pairs_out.front().second))) {
                    
#ifdef debug_multipath_mapper_mapping
                    cerr << "one end of the pair may be mismapped, attempting individual end mappings" << endl;
#endif
                    // we're not happy with the pairs we got, try to get a good pair by rescuing from single ended alignments
                    // but block rescue from any sides that we already tried rescue from in the repeat rescue routine
                    
                    vector<pair<MultipathAlignment, MultipathAlignment>> rescue_aln_pairs;
                    vector<pair<pair<size_t, size_t>, int64_t>> rescue_distances;
                    bool rescued = align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                                       do_repeat_rescue_from_1, do_repeat_rescue_from_2,
                                                                       rescue_aln_pairs, rescue_distances, max_alt_mappings);
                    
                    // if we find consistent pairs by rescue, merge the two lists
                    if (rescued) {
#ifdef debug_multipath_mapper_mapping
                        cerr << "found some rescue pairs, merging into current list of consistent mappings" << endl;
#endif
                        
                        merge_rescued_mappings(multipath_aln_pairs_out, cluster_pairs, rescue_aln_pairs, rescue_distances);
                        
                        // if we still haven't found mappings that are distinguishable from matches to random sequences,
                        // don't let them have any mapping quality
                        if (likely_mismapping(multipath_aln_pairs_out.front().first) ||
                            likely_mismapping(multipath_aln_pairs_out.front().second)) {
                            multipath_aln_pairs_out.front().first.set_mapping_quality(0);
                            multipath_aln_pairs_out.front().second.set_mapping_quality(0);
                        }
                    }
                    else {
                        // rescue didn't find any consistent mappings, revert to the single ended mappings
                        std::swap(multipath_aln_pairs_out, rescue_aln_pairs);
                    }
                }
                else if (multipath_aln_pairs_out.front().first.mapping_quality() == max_mapping_quality &&
                         multipath_aln_pairs_out.front().second.mapping_quality() == max_mapping_quality) {
                    
                    // we're very confident about this pair, but it might be because we over-pruned at the clustering stage
                    // so we use this routine to use rescue on other very good looking independent end clusters
                    attempt_rescue_for_secondaries(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                   multipath_aln_pairs_out, cluster_pairs);
                }
            }
            else {
                // revert to independent single ended mappings, but skip any rescues that we already tried
                
#ifdef debug_multipath_mapper_mapping
                cerr << "could not find a consistent pair, reverting to single ended mapping" << endl;
#endif
                align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2, do_repeat_rescue_from_1,
                                                    do_repeat_rescue_from_2, multipath_aln_pairs_out, cluster_pairs, max_alt_mappings);
                
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
        
        // remove the full length bonus if we don't want it in the final score
        if (strip_bonuses) {
            for (pair<MultipathAlignment, MultipathAlignment>& multipath_aln_pair : multipath_aln_pairs_out) {
                strip_full_length_bonuses(multipath_aln_pair.first);
                strip_full_length_bonuses(multipath_aln_pair.second);
            }
        }
        
        // add pair names
        for (pair<MultipathAlignment, MultipathAlignment>& multipath_aln_pair : multipath_aln_pairs_out) {
            multipath_aln_pair.first.set_paired_read_name(multipath_aln_pair.second.name());
            multipath_aln_pair.second.set_paired_read_name(multipath_aln_pair.first.name());
        }
        
        // clean up the VG objects on the heap
        for (auto cluster_graph : cluster_graphs1) {
            delete get<0>(cluster_graph);
        }
        for (auto cluster_graph : cluster_graphs2) {
            delete get<0>(cluster_graph);
        }
    }
    
    void MultipathMapper::split_multicomponent_alignments(vector<MultipathAlignment>& multipath_alns_out) const {
        
        for (size_t i = 0; i < multipath_alns_out.size(); i++) {
            
            vector<vector<int64_t>> comps = connected_components(multipath_alns_out[i]);
            
            if (comps.size() > 1) {
#ifdef debug_multipath_mapper_mapping
                cerr << "splitting multicomponent alignment " << pb2json(multipath_alns_out[i]) << endl;
#endif
                // split this multipath alignment into its connected components
                for (size_t j = 1; j < comps.size(); j++) {
                    multipath_alns_out.emplace_back();
                    extract_sub_multipath_alignment(multipath_alns_out[i], comps[j],
                                                    multipath_alns_out.back());;
                }
                // put the first component into the original location
                MultipathAlignment last_component;
                extract_sub_multipath_alignment(multipath_alns_out[i], comps[0], last_component);
                multipath_alns_out[i] = last_component;
            }
        }
    }
    
    void MultipathMapper::attempt_rescue_of_repeat_from_non_repeat(const Alignment& alignment1, const Alignment& alignment2,
                                                                   const vector<MaximalExactMatch>& mems1, const vector<MaximalExactMatch>& mems2,
                                                                   bool do_repeat_rescue_from_1, bool do_repeat_rescue_from_2,
                                                                   vector<memcluster_t>& clusters1, vector<memcluster_t>& clusters2,
                                                                   vector<clustergraph_t>& cluster_graphs1, vector<clustergraph_t>& cluster_graphs2,
                                                                   vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                                   vector<pair<pair<size_t, size_t>, int64_t>>& pair_distances, size_t max_alt_mappings,
                                                                   OrientedDistanceClusterer::paths_of_node_memo_t* paths_of_node_memo,
                                                                   OrientedDistanceClusterer::oriented_occurences_memo_t* oriented_occurences_memo,
                                                                   OrientedDistanceClusterer::handle_memo_t* handle_memo) {
        
        bool rescue_succeeded_from_1 = false, rescue_succeeded_from_2 = false;
        
        if (do_repeat_rescue_from_1) {
            
#ifdef debug_multipath_mapper_mapping
            cerr << "attempting repeat rescue from read 1" << endl;
#endif
            
            // get the clusters for the non repeat
            if (adjust_alignments_for_base_quality) {
                OrientedDistanceClusterer clusterer1(alignment1, mems1, *get_qual_adj_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                     unstranded_clustering, paths_of_node_memo, oriented_occurences_memo, handle_memo);
                clusters1 = clusterer1.clusters(max_mapping_quality, log_likelihood_approx_factor);
            }
            else {
                OrientedDistanceClusterer clusterer1(alignment1, mems1, *get_regular_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                     unstranded_clustering, paths_of_node_memo, oriented_occurences_memo, handle_memo);
                clusters1 = clusterer1.clusters(max_mapping_quality, log_likelihood_approx_factor);
            }
            
            // extract the graphs around the clusters
            cluster_graphs1 = query_cluster_graphs(alignment1, mems1, clusters1);
            
            // attempt rescue from these graphs
            vector<pair<MultipathAlignment, MultipathAlignment>> rescued_pairs;
            vector<pair<pair<size_t, size_t>, int64_t>> rescued_distances;
            rescue_succeeded_from_1 = align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                                          false, true, rescued_pairs, pair_distances, max_alt_mappings);
            
            // move the rescued pairs to the output vectors
            if (rescue_succeeded_from_1) {
                
#ifdef debug_multipath_mapper_mapping
                cerr << "repeat rescue succeeded" << endl;
#endif
                
                for (auto& multipath_aln_pair : rescued_pairs) {
                    multipath_aln_pairs_out.emplace_back(move(multipath_aln_pair));
                }
                for (auto& pair_distance : rescued_distances) {
                    pair_distances.emplace_back(pair_distance);
                }
            }
        }
        
        if (do_repeat_rescue_from_2) {
            // TODO: duplicative code
            
#ifdef debug_multipath_mapper_mapping
            cerr << "attempting repeat rescue from read 2" << endl;
#endif
            
            // get the clusters for the non repeat
            if (adjust_alignments_for_base_quality) {
                OrientedDistanceClusterer clusterer2(alignment2, mems2, *get_qual_adj_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                     unstranded_clustering, paths_of_node_memo, oriented_occurences_memo, handle_memo);
                clusters2 = clusterer2.clusters(max_mapping_quality, log_likelihood_approx_factor);
            }
            else {
                OrientedDistanceClusterer clusterer2(alignment2, mems2, *get_regular_aligner(), xindex, max_expected_dist_approx_error, min_clustering_mem_length,
                                                     unstranded_clustering, paths_of_node_memo, oriented_occurences_memo, handle_memo);
                clusters2 = clusterer2.clusters(max_mapping_quality, log_likelihood_approx_factor);
            }
            
            // extract the graphs around the clusters
            cluster_graphs2 = query_cluster_graphs(alignment2, mems2, clusters2);
            
            // attempt rescue from these graphs
            vector<pair<MultipathAlignment, MultipathAlignment>> rescued_pairs;
            vector<pair<pair<size_t, size_t>, int64_t>> rescued_distances;
            rescue_succeeded_from_2 = align_to_cluster_graphs_with_rescue(alignment1, alignment2, cluster_graphs1, cluster_graphs2,
                                                                          true, false, rescued_pairs, pair_distances, max_alt_mappings);
            
            // move the rescued pairs to the output vectors
            if (rescue_succeeded_from_2) {
                
#ifdef debug_multipath_mapper_mapping
                cerr << "repeat rescue succeeded" << endl;
#endif
                for (auto& multipath_aln_pair : rescued_pairs) {
                    multipath_aln_pairs_out.emplace_back(move(multipath_aln_pair));
                }
                for (auto& pair_distance : rescued_distances) {
                    pair_distances.emplace_back(pair_distance);
                }
            }
        }
        
        // re-sort the rescued alignments if we actually did it from both sides
        if (rescue_succeeded_from_1 && rescue_succeeded_from_2) {
            sort_and_compute_mapping_quality(multipath_aln_pairs_out, pair_distances);
        }
    }
    
    void MultipathMapper::merge_rescued_mappings(vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                 vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                 vector<pair<MultipathAlignment, MultipathAlignment>>& rescued_multipath_aln_pairs,
                                                 vector<pair<pair<size_t, size_t>, int64_t>>& rescued_cluster_pairs) const {
                
        size_t num_unrescued_pairs = multipath_aln_pairs_out.size();
        for (size_t j = 0; j < rescued_multipath_aln_pairs.size(); j++) {
            
            // make sure this pair isn't a duplicate with any of the original pairs
            bool duplicate = false;
            for (size_t i = 0; i < num_unrescued_pairs; i++) {
#ifdef debug_multipath_mapper_mapping
                cerr << "checking if rescue pair " << j << " is duplicate of original pair " << i << endl;
#endif
                if (abs(distance_between(multipath_aln_pairs_out[i].first, rescued_multipath_aln_pairs[j].first)) < 20) {
                    if (abs(distance_between(multipath_aln_pairs_out[i].second, rescued_multipath_aln_pairs[j].second)) < 20) {
#ifdef debug_multipath_mapper_mapping
                        cerr << "found a duplicate" << endl;
#endif
                        duplicate = true;
                        break;
                    }
                }
            }
            
            if (!duplicate) {
#ifdef debug_multipath_mapper_mapping
                cerr << "no duplicate, adding to return vector if distance is finite and positive" << endl;
#endif
                // add a dummy pair to hold the distance
                cluster_pairs.emplace_back(rescued_cluster_pairs[j]);
                multipath_aln_pairs_out.emplace_back(move(rescued_multipath_aln_pairs[j]));
            }
        }
        
        sort_and_compute_mapping_quality(multipath_aln_pairs_out, cluster_pairs);
    }
    
    void MultipathMapper::split_multicomponent_alignments(vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                          vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs) const {
        
        for (size_t i = 0; i < multipath_aln_pairs_out.size(); i++) {
            vector<vector<int64_t>> connected_components_1 = connected_components(multipath_aln_pairs_out[i].first);
            vector<vector<int64_t>> connected_components_2 = connected_components(multipath_aln_pairs_out[i].second);
            
#ifdef debug_multipath_mapper_mapping
            cerr << "finding connected components for mapping:" << endl;
            cerr  << pb2json(multipath_aln_pairs_out[i].first) << endl;
            cerr  << pb2json(multipath_aln_pairs_out[i].second) << endl;
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
            
            
            if (connected_components_1.size() > 1 && connected_components_2.size() > 1) {
#ifdef debug_multipath_mapper_mapping
                cerr << "splitting both multicomponent alignments" << endl;
#endif
                // need to split both ends
                for (size_t j = 0; j < connected_components_1.size(); j++) {
                    for (size_t k = 0; k < connected_components_2.size(); k++) {
                        if (j == 0 && k == 0) {
                            // save the first component to go back into the original index in the return vector
                            continue;
                        }
                        multipath_aln_pairs_out.emplace_back();
                        extract_sub_multipath_alignment(multipath_aln_pairs_out[i].first, connected_components_1[j],
                                                        multipath_aln_pairs_out.back().first);
                        extract_sub_multipath_alignment(multipath_aln_pairs_out[i].second, connected_components_2[k],
                                                        multipath_aln_pairs_out.back().second);
                        // add a distance for scoring purposes (possibly unstranded since this comes before enforcing strand consistency)
                        cluster_pairs.emplace_back(cluster_pairs[i].first,
                                                   distance_between(multipath_aln_pairs_out.back().first, multipath_aln_pairs_out.back().second,
                                                                    true, unstranded_clustering));
#ifdef debug_multipath_mapper_mapping
                        cerr << "added component pair at distance " << cluster_pairs.back().second << ":" << endl;
                        cerr  << pb2json(multipath_aln_pairs_out.back().first) << endl;
                        cerr  << pb2json(multipath_aln_pairs_out.back().second) << endl;
#endif                    
                    }
                }
                // put the first component pair into the original location
                MultipathAlignment last_component;
                extract_sub_multipath_alignment(multipath_aln_pairs_out[i].first, connected_components_1[0], last_component);
                multipath_aln_pairs_out[i].first = last_component;
                last_component.Clear();
                extract_sub_multipath_alignment(multipath_aln_pairs_out[i].second, connected_components_2[0], last_component);
                multipath_aln_pairs_out[i].second = last_component;
                cluster_pairs[i].second = distance_between(multipath_aln_pairs_out[i].first, multipath_aln_pairs_out[i].second,
                                                           true, unstranded_clustering);
                
#ifdef debug_multipath_mapper_mapping
                cerr << "added component pair at distance " << cluster_pairs[i].second << ":" << endl;
                cerr  << pb2json(multipath_aln_pairs_out[i].first) << endl;
                cerr  << pb2json(multipath_aln_pairs_out[i].second) << endl;
#endif
                
            }
            else if (connected_components_1.size() > 1) {
#ifdef debug_multipath_mapper_mapping
                cerr << "splitting read 1 multicomponent alignments" << endl;
#endif
                // only need to split first end
                for (size_t j = 1; j < connected_components_1.size(); j++) {
                    multipath_aln_pairs_out.emplace_back(MultipathAlignment(), multipath_aln_pairs_out[i].second);
                    extract_sub_multipath_alignment(multipath_aln_pairs_out[i].first, connected_components_1[j],
                                                    multipath_aln_pairs_out.back().first);
                    // add a distance for scoring purposes (possibly unstranded since this comes before enforcing strand consistency)
                    cluster_pairs.emplace_back(cluster_pairs[i].first,
                                               distance_between(multipath_aln_pairs_out.back().first, multipath_aln_pairs_out.back().second,
                                                                true, unstranded_clustering));
#ifdef debug_multipath_mapper_mapping
                    cerr << "added component pair at distance " << cluster_pairs.back().second << ":" << endl;
                    cerr  << pb2json(multipath_aln_pairs_out.back().first) << endl;
                    cerr  << pb2json(multipath_aln_pairs_out.back().second) << endl;
#endif
                }
                // put the first component into the original location
                MultipathAlignment last_component;
                extract_sub_multipath_alignment(multipath_aln_pairs_out[i].first, connected_components_1[0], last_component);
                multipath_aln_pairs_out[i].first = last_component;
                cluster_pairs[i].second = distance_between(multipath_aln_pairs_out[i].first, multipath_aln_pairs_out[i].second,
                                                           true, unstranded_clustering);
                
#ifdef debug_multipath_mapper_mapping
                cerr << "added component pair at distance " << cluster_pairs[i].second << ":" << endl;
                cerr  << pb2json(multipath_aln_pairs_out[i].first) << endl;
                cerr  << pb2json(multipath_aln_pairs_out[i].second) << endl;
#endif
            }
            else if (connected_components_2.size() > 1) {
#ifdef debug_multipath_mapper_mapping
                cerr << "splitting read 2 multicomponent alignments" << endl;
#endif
                // only need to split second end
                for (size_t j = 1; j < connected_components_2.size(); j++) {
                    multipath_aln_pairs_out.emplace_back(multipath_aln_pairs_out[i].first, MultipathAlignment());
                    extract_sub_multipath_alignment(multipath_aln_pairs_out[i].second, connected_components_2[j],
                                                    multipath_aln_pairs_out.back().second);
                    // add a distance for scoring purposes (possibly unstranded since this comes before enforcing strand consistency)
                    cluster_pairs.emplace_back(cluster_pairs[i].first,
                                               distance_between(multipath_aln_pairs_out.back().first, multipath_aln_pairs_out.back().second,
                                                                true, unstranded_clustering));
#ifdef debug_multipath_mapper_mapping
                    cerr << "added component pair at distance " << cluster_pairs.back().second << ":" << endl;
                    cerr  << pb2json(multipath_aln_pairs_out.back().first) << endl;
                    cerr  << pb2json(multipath_aln_pairs_out.back().second) << endl;
#endif
                }
                // put the first component into the original location
                MultipathAlignment last_component;
                extract_sub_multipath_alignment(multipath_aln_pairs_out[i].second, connected_components_2[0], last_component);
                multipath_aln_pairs_out[i].second = last_component;
                cluster_pairs[i].second = distance_between(multipath_aln_pairs_out[i].first, multipath_aln_pairs_out[i].second,
                                                           true, unstranded_clustering);
                
#ifdef debug_multipath_mapper_mapping
                cerr << "added component pair at distance " << cluster_pairs[i].second << ":" << endl;
                cerr  << pb2json(multipath_aln_pairs_out[i].first) << endl;
                cerr  << pb2json(multipath_aln_pairs_out[i].second) << endl;
#endif
            }
        }
    }
    
    void MultipathMapper::align_to_cluster_graph_pairs(const Alignment& alignment1, const Alignment& alignment2,
                                                       vector<clustergraph_t>& cluster_graphs1,
                                                       vector<clustergraph_t>& cluster_graphs2,
                                                       vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs,
                                                       vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs_out,
                                                       size_t max_alt_mappings,
                                                       OrientedDistanceClusterer::paths_of_node_memo_t* paths_of_node_memo,
                                                       OrientedDistanceClusterer::oriented_occurences_memo_t* oriented_occurences_memo,
                                                       OrientedDistanceClusterer::handle_memo_t* handle_memo) {
        
        assert(multipath_aln_pairs_out.empty());
        
        auto get_pair_coverage = [&](const pair<size_t, size_t>& cluster_pair) {
            return get<2>(cluster_graphs1[cluster_pair.first]) + get<2>(cluster_graphs2[cluster_pair.second]);
        };
        
        // sort the pairs descending by total unique sequence coverage
        std::sort(cluster_pairs.begin(), cluster_pairs.end(),
                  [&](const pair<pair<size_t, size_t>, int64_t>& a, const pair<pair<size_t, size_t>, int64_t>& b) {
                      // We need to be able to look up the coverage for the graph an input cluster went into.
                      // Compute total coverage following all the redirects and see if
                      // it's in the right order.
                      return get_pair_coverage(a.first) > get_pair_coverage(b.first);
                  });
        
#ifdef debug_multipath_mapper_mapping
        cerr << "aligning to cluster pairs..." << endl;
#endif
        
        // we may need to compute an extra mapping above the one we'll report if we're computing mapping quality
        size_t num_mappings_to_compute = mapping_quality_method != None ? max(num_mapping_attempts, (size_t) 2) : num_mapping_attempts;
        
        // TODO: some cluster pairs will produce redundant subgraph pairs.
        // We'll end up with redundant pairs being output.
        
        // align to each cluster pair
        multipath_aln_pairs_out.reserve(min(num_mappings_to_compute, cluster_pairs.size()));
        size_t num_mappings = 0;
        for (size_t i = 0; i < cluster_pairs.size(); ++i) {
            // For each cluster pair
            const pair<pair<size_t, size_t>, int64_t>& cluster_pair = cluster_pairs[i];
            
            // if we have a cluster graph pair with small enough MEM coverage
            // compared to the best one or we've made the maximum number of
            // alignments we stop producing alternate alignments
            if (get_pair_coverage(cluster_pair.first) < mem_coverage_min_ratio * get_pair_coverage(cluster_pairs[0].first)
                || num_mappings >= num_mappings_to_compute) {
                
                // remove the rest of the cluster pairs so we maintain the invariant that there are the
                // same number of cluster pairs as alternate mappings
                cluster_pairs.resize(i);
                
                break;
            }
            
            VG* vg1 = get<0>(cluster_graphs1[cluster_pair.first.first]);
            VG* vg2 = get<0>(cluster_graphs2[cluster_pair.first.second]);
            
            memcluster_t& graph_mems1 = get<1>(cluster_graphs1[cluster_pair.first.first]);
            memcluster_t& graph_mems2 = get<1>(cluster_graphs2[cluster_pair.first.second]);
            
#ifdef debug_multipath_mapper_mapping
            cerr << "doing pair " << cluster_pair.first.first << " " << cluster_pair.first.second << endl;
            cerr << "performing alignments to subgraphs " << pb2json(vg1->graph) << " and " << pb2json(vg2->graph) << endl;
#endif
            
            // Do the two alignments
            multipath_aln_pairs_out.emplace_back();
            multipath_align(alignment1, vg1, graph_mems1, multipath_aln_pairs_out.back().first);
            multipath_align(alignment2, vg2, graph_mems2, multipath_aln_pairs_out.back().second);
            
            num_mappings++;
        }
        
        // split up any multi-component multipath alignments
        split_multicomponent_alignments(multipath_aln_pairs_out, cluster_pairs);
        
        // downstream algorithms assume multipath alignments are topologically sorted (including the scoring
        // algorithm in the next step)
        for (pair<MultipathAlignment, MultipathAlignment>& multipath_aln_pair : multipath_aln_pairs_out) {
            topologically_order_subpaths(multipath_aln_pair.first);
            topologically_order_subpaths(multipath_aln_pair.second);
        }
        
        // if we haven't been checking strand consistency, enforce it now at the end
        if (unstranded_clustering) {
            establish_strand_consistency(multipath_aln_pairs_out, cluster_pairs, paths_of_node_memo, oriented_occurences_memo, handle_memo);
        }
        
        // put pairs in score sorted order and compute mapping quality of best pair using the score
        sort_and_compute_mapping_quality(multipath_aln_pairs_out, cluster_pairs);
        
#ifdef debug_validate_multipath_alignments
        for (pair<MultipathAlignment, MultipathAlignment>& multipath_aln_pair : multipath_aln_pairs_out) {
#ifdef debug_multipath_mapper_mapping
            cerr << "validating multipath alignments:" << endl;
            cerr << pb2json(multipath_aln_pair.first) << endl;
            cerr << pb2json(multipath_aln_pair.second) << endl;
#endif
            if (!validate_multipath_alignment(multipath_aln_pair.first)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.first.name() << " failed to validate" << endl;
            }
            if (!validate_multipath_alignment(multipath_aln_pair.second)) {
                cerr << "### WARNING ###" << endl;
                cerr << "multipath alignment of read " << multipath_aln_pair.second.name() << " failed to validate" << endl;
            }
        }
#endif
        
    }
    
    auto MultipathMapper::query_cluster_graphs(const Alignment& alignment,
                                               const vector<MaximalExactMatch>& mems,
                                               const vector<memcluster_t>& clusters) -> vector<clustergraph_t> {
        
        // Figure out the aligner to use
        BaseAligner* aligner = get_aligner();
        
        // We populate this with all the cluster graphs.
        vector<clustergraph_t> cluster_graphs_out;
        
        // we will ensure that nodes are in only one cluster, use this to record which one
        unordered_map<id_t, size_t> node_id_to_cluster;
        
        // to hold the clusters as they are (possibly) merged
        unordered_map<size_t, VG*> cluster_graphs;
        
        // to keep track of which clusters have been merged
        UnionFind union_find(clusters.size());
        
        for (size_t i = 0; i < clusters.size(); i++) {
            
#ifdef debug_multipath_mapper_mapping
            cerr << "extracting subgraph for cluster " << i << endl;
#endif
            
            // gather the parameters for subgraph extraction from the MEM hits
            
            const memcluster_t& cluster = clusters[i];
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
                forward_max_dist.push_back(aligner->longest_detectable_gap(alignment, mem_hit.first->end)
                                           + (alignment.sequence().end() - mem_hit.first->begin));
                backward_max_dist.push_back(aligner->longest_detectable_gap(alignment, mem_hit.first->begin)
                                            + (mem_hit.first->begin - alignment.sequence().begin()));
            }
            
            
            // TODO: a progressive expansion of the subgraph if the MEM hit is already contained in
            // a cluster graph somewhere?
            
            // extract the subgraph within the search distance
            
            VG* cluster_graph = new VG();
            Graph& graph = cluster_graph->graph;
            
            // extract the protobuf Graph in place in the VG
            algorithms::extract_containing_graph(xindex, graph, positions, forward_max_dist,
                                                 backward_max_dist);
            
            // check if this subgraph overlaps with any previous subgraph (indicates a probable clustering failure where
            // one cluster was split into multiple clusters)
            unordered_set<size_t> overlapping_graphs;
            for (size_t j = 0; j < graph.node_size(); j++) {
                id_t node_id = graph.node(j).id();
                if (node_id_to_cluster.count(node_id)) {
                    overlapping_graphs.insert(node_id_to_cluster[node_id]);
                }
                else {
                    node_id_to_cluster[node_id] = i;
                }
            }
            
            if (overlapping_graphs.empty()) {
                // there is no overlap with any other graph, suggesting a new unique hit
                
#ifdef debug_multipath_mapper_mapping
                cerr << "cluster graph does not overlap with any other cluster graphs, adding as cluster " << i << endl;
#endif
                cluster_graphs[i] = cluster_graph;
                
                // now that we know we're going to save the graph, manually trigger the index building since
                // we circumvented the constructors
                cluster_graph->build_indexes();
            }
            else {
                // this graph overlaps at least one other graph, so we merge them into one
                
#ifdef debug_multipath_mapper_mapping
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
                
#ifdef debug_multipath_mapper_mapping
                cerr << "merging as cluster " << remaining_idx << endl;
#endif
                
                VG* merging_graph;
                if (remaining_idx == i) {
                    // the new graph was chosen to remain, so add it to the record
                    cluster_graphs[i] = cluster_graph;
                    merging_graph = cluster_graph;
                }
                else {
                    // the new graph will be merged into an existing graph
                    merging_graph = cluster_graphs[remaining_idx];
                    merging_graph->extend(graph);
                    delete cluster_graph;
                }
                
                // merge any other chained graphs into the remaining graph
                for (size_t j : overlapping_graphs) {
                    if (j != remaining_idx) {
                        VG* removing_graph = cluster_graphs[j];
                        merging_graph->extend(removing_graph->graph);
                        delete removing_graph;
                        cluster_graphs.erase(j);
                    }
                }
                
                Graph& merged_graph = merging_graph->graph;
                for (size_t j = 0; j < merged_graph.node_size(); j++) {
                    node_id_to_cluster[merged_graph.node(j).id()] = remaining_idx;
                }
            }
        }
        
        // check if any cluster graphs pulled disconnected components (as a result of a clustering failure)
        // and if so split them up
        
        // keeps track of the connected components of any multicomponent graph
        vector<pair<size_t, vector<unordered_set<id_t>>>> multicomponent_graphs;
        unordered_map<size_t, vector<size_t>> multicomponent_splits;
        
        size_t max_graph_idx = 0;
        for (const pair<size_t, VG*> cluster_graph : cluster_graphs) {
            vector<unordered_set<id_t>> connected_components = algorithms::weakly_connected_components(cluster_graph.second);
            if (connected_components.size() > 1) {
                multicomponent_graphs.emplace_back(cluster_graph.first, std::move(connected_components));
            }
            max_graph_idx = max(cluster_graph.first, max_graph_idx);
        }
            
        // did we find any graphs that consist of disconnected components?
        for (pair<size_t, vector<unordered_set<id_t>>>& multicomponent_graph : multicomponent_graphs) {
            max_graph_idx++;
            // make a new graph for each of the components
#ifdef debug_multipath_mapper_mapping
            cerr << "cluster graph " << multicomponent_graph.first << " has multiple connected components, splitting now" << endl;
            for (size_t i = 0; i < multicomponent_graph.second.size(); i++) {
                cerr << "component " << max_graph_idx + i << ":" << endl;
                
                for (auto node_id : multicomponent_graph.second[i]) {
                    cerr << "\t" << node_id << endl;
                }
            }
#endif
            
            for (size_t i = 0; i < multicomponent_graph.second.size(); i++) {
                cluster_graphs[max_graph_idx + i] = new VG();
            }
            
            Graph& joined_graph = cluster_graphs[multicomponent_graph.first]->graph;
            
            // divvy up the nodes
            for (size_t i = 0; i < joined_graph.node_size(); i++) {
                const Node& node = joined_graph.node(i);
                for (size_t j = 0; j < multicomponent_graph.second.size(); j++) {
                    if (multicomponent_graph.second[j].count(node.id())) {
                        cluster_graphs[max_graph_idx + j]->add_node(node);
                        node_id_to_cluster[node.id()] = max_graph_idx + j;
                        break;
                        }
                    }
                }
            
            // divvy up the edges
            for (size_t i = 0; i < joined_graph.edge_size(); i++) {
                const Edge& edge = joined_graph.edge(i);
                for (size_t j = 0; j < multicomponent_graph.second.size(); j++) {
                    if (multicomponent_graph.second[j].count(edge.from())) {
                        cluster_graphs[max_graph_idx + j]->add_edge(edge);
                        break;
                    }
                }
            }
            
#ifdef debug_multipath_mapper_mapping
            cerr << "split graphs:" << endl;
            for (size_t i = 0; i < multicomponent_graph.second.size(); i++) {
                cerr << "component " << max_graph_idx + i << ":" << endl;
                cerr << pb2json(cluster_graphs[max_graph_idx + i]->graph) << endl;
            }
#endif
            
            // remove the old graph
            delete cluster_graphs[multicomponent_graph.first];
            cluster_graphs.erase(multicomponent_graph.first);
            
            max_graph_idx += multicomponent_graph.second.size();
        }
        
        // move the pointers to the return vector and figure out which graph in the return
        // vector each MEM cluster ended up in
        cluster_graphs_out.reserve(cluster_graphs.size());
        unordered_map<size_t, size_t> cluster_to_idx;
        for (const auto& cluster_graph : cluster_graphs) {
#ifdef debug_multipath_mapper_mapping
            cerr << "adding cluster graph " << cluster_graph.first << " to return vector at index " << cluster_graphs_out.size() << endl;
#endif
            cluster_to_idx[cluster_graph.first] = cluster_graphs_out.size();
            cluster_graphs_out.emplace_back(cluster_graph.second, memcluster_t(), 0);
        }

        
#ifdef debug_multipath_mapper_mapping
        cerr << "computing MEM assignments to cluster graphs" << endl;
#endif
        // which MEMs are in play for which cluster?
        for (const MaximalExactMatch& mem : mems) {
            for (gcsa::node_type hit : mem.nodes) {
                id_t node_id = gcsa::Node::id(hit);
                if (node_id_to_cluster.count(node_id)) {
                    size_t cluster_idx = cluster_to_idx[node_id_to_cluster[node_id]];
                    get<1>(cluster_graphs_out[cluster_idx]).push_back(make_pair(&mem, make_pos_t(hit)));
#ifdef debug_multipath_mapper_mapping
                    cerr << "\tMEM " << mem.sequence() << " at " << make_pos_t(hit) << " found in cluster " << node_id_to_cluster[node_id] << " at index " << cluster_idx << endl;
#endif
                }
            }
        }
        
        // compute the read coverage of each cluster graph and sort the assigned MEMs by length
        // and then lexicographically by read index
        for (size_t i = 0; i < cluster_graphs_out.size(); i++) {
            auto& cluster_graph = cluster_graphs_out[i];
            get<2>(cluster_graph) = read_coverage(get<1>(cluster_graph));
#ifdef debug_multipath_mapper_mapping
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
        
            
        // sort the cluster graphs descending by unique sequence coverage
        std::sort(cluster_graphs_out.begin(), cluster_graphs_out.end(),
                  [](const clustergraph_t& cluster_graph_1,
                     const clustergraph_t& cluster_graph_2) {
                      return get<2>(cluster_graph_1) > get<2>(cluster_graph_2);
                  });
        
        return move(cluster_graphs_out);
        
        
    }
    
    void MultipathMapper::multipath_align(const Alignment& alignment, VG* vg,
                                          memcluster_t& graph_mems,
                                          MultipathAlignment& multipath_aln_out) const {

#ifdef debug_multipath_mapper_alignment
        cerr << "constructing alignment graph" << endl;
#endif
        
        // the longest path we could possibly align to (full gap and a full sequence)
        size_t target_length = alignment.sequence().size() + get_aligner()->longest_detectable_gap(alignment);
        
        // convert from bidirected to directed
        unordered_map<id_t, pair<id_t, bool> > node_trans;
        
        // check if we can get away with using only one strand of the graph
        bool use_single_stranded = vg->is_single_stranded();
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
        // TODO: I do this without the copy constructor for the forward strand?
        VG align_graph = use_single_stranded ? (mem_strand ? vg->reverse_complement_graph(node_trans) : *vg) : vg->split_strands(node_trans);
        
        // if we are using only the forward strand of the current graph, a make trivial node translation so
        // the later code's expectations are met
        if (use_single_stranded && !mem_strand) {
            vg->identity_translation(node_trans);
        }
        
        // if necessary, convert from cyclic to acylic
        if (!vg->is_directed_acyclic()) {
            unordered_map<id_t, pair<id_t, bool> > dagify_trans;
            align_graph = align_graph.dagify(target_length, // high enough that num SCCs is never a limiting factor
                                             dagify_trans,
                                             target_length,
                                             0); // no maximum on size of component
            node_trans = align_graph.overlay_node_translations(dagify_trans, node_trans);
        }
        
        // put the internal graph in topological order for the MultipathAlignmentGraph algorithm
        align_graph.lazy_sort();
        
#ifdef debug_multipath_mapper_alignment
        cerr << "making multipath alignment MEM graph" << endl;
#endif
        
        // construct a graph that summarizes reachability between MEMs
        MultipathAlignmentGraph multi_aln_graph(align_graph, graph_mems, node_trans, gcsa, snarl_manager, max_snarl_cut_size);
        
        vector<size_t> topological_order;
        multi_aln_graph.topological_sort(topological_order);
        
#ifdef debug_multipath_mapper_alignment
        cerr << "computed topological sort, topology is: " << endl;
        for (size_t i = 0; i < multi_aln_graph.match_nodes.size(); i++) {
            cerr << "node " << i << ", " << pb2json(multi_aln_graph.match_nodes[i].path.mapping(0).position()) << " ";
            for (auto iter = multi_aln_graph.match_nodes[i].begin; iter != multi_aln_graph.match_nodes[i].end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            for (pair<size_t, size_t> edge : multi_aln_graph.match_nodes[i].edges) {
                cerr << "\tto " << edge.first << ", dist " << edge.second << endl;
            }
        }
#endif
        
        // it's sometimes possible for transitive edges to survive the original construction algorithm, so remove them
        multi_aln_graph.remove_transitive_edges(topological_order);
        
#ifdef debug_multipath_mapper_alignment
        cerr << "removed transitive edges, topology is:" << endl;
        for (size_t i = 0; i < multi_aln_graph.match_nodes.size(); i++) {
            cerr << "node " << i << ", " << pb2json(multi_aln_graph.match_nodes[i].path.mapping(0).position()) << " ";
            for (auto iter = multi_aln_graph.match_nodes[i].begin; iter != multi_aln_graph.match_nodes[i].end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            for (pair<size_t, size_t> edge : multi_aln_graph.match_nodes[i].edges) {
                cerr << "\tto " << edge.first << ", dist " << edge.second << endl;
            }
        }
#endif
        
        // prune this graph down the paths that have reasonably high likelihood
        multi_aln_graph.prune_to_high_scoring_paths(alignment, get_aligner(),
                                                    max_suboptimal_path_score_ratio, topological_order);
        
#ifdef debug_multipath_mapper_alignment
        cerr << "pruned to high scoring paths, topology is:" << endl;
        for (size_t i = 0; i < multi_aln_graph.match_nodes.size(); i++) {
            cerr << "node " << i << ", " << pb2json(multi_aln_graph.match_nodes[i].path.mapping(0).position()) << " ";
            for (auto iter = multi_aln_graph.match_nodes[i].begin; iter != multi_aln_graph.match_nodes[i].end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            for (pair<size_t, size_t> edge : multi_aln_graph.match_nodes[i].edges) {
                cerr << "\tto " << edge.first << ", dist " << edge.second << endl;
            }
        }
#endif
        
        // create a new multipath alignment object and transfer over data from alignment
        multipath_aln_out.set_sequence(alignment.sequence());
        multipath_aln_out.set_quality(alignment.quality());
        multipath_aln_out.set_name(alignment.name());
        multipath_aln_out.set_sample_name(alignment.sample_name());
        multipath_aln_out.set_read_group(alignment.read_group());
        
#ifdef debug_multipath_mapper_alignment
        cerr << "transferred over read information" << endl;
#endif
        
        // add a subpath for each of the exact match nodes
        for (int64_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
            ExactMatchNode& match_node = multi_aln_graph.match_nodes[j];
            Subpath* subpath = multipath_aln_out.add_subpath();
            *subpath->mutable_path() = match_node.path;
            int32_t match_score = get_aligner()->score_exact_match(match_node.begin, match_node.end,
                alignment.quality().begin() + (match_node.begin - alignment.sequence().begin()));
                
            subpath->set_score(match_score + get_aligner()->full_length_bonus *
                ((match_node.begin == alignment.sequence().begin()) +
                (match_node.end == alignment.sequence().end())));
        }
        
#ifdef debug_multipath_mapper_alignment
        cerr << "doing DP between MEMs" << endl;
#endif
        
        // perform alignment in the intervening sections
        for (int64_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
#ifdef debug_multipath_mapper_alignment
            cerr << "checking for intervening alignments from match node " << j << endl;
#endif
            
            ExactMatchNode& src_match_node = multi_aln_graph.match_nodes[j];
            Subpath* src_subpath = multipath_aln_out.mutable_subpath(j);
            
            // check whether this is a simple split MEM, in which case we can make the connection directly
            // and avoid making an intervening alignment
            if (src_match_node.edges.size() == 1) {
                pair<size_t, size_t>& edge = src_match_node.edges[0];
                // check whether the positions are right next to each other in both the read and the graph
                //
                // strictly speaking, there could be other paths between the match nodes if the ends of the paths
                // are at node boundaries, but even then the best alignment would be from the empty sequence to
                // the path of length 0
                if (src_match_node.end == multi_aln_graph.match_nodes[edge.first].begin && edge.second == 0) {
#ifdef debug_multipath_mapper_alignment
                    cerr << "match node has a simple split connection onto node " << edge.first << endl;
#endif
                    // connect the subpath directly without an intervening alignment
                    src_subpath->add_next(edge.first);
                    continue;
                }
            }
            
            const Path& path = src_subpath->path();
            const Mapping& final_mapping = path.mapping(path.mapping_size() - 1);
            const Position& final_mapping_position = final_mapping.position();
            // make a pos_t that points to the final base in the match
            pos_t src_pos = make_pos_t(final_mapping_position.node_id(),
                                       final_mapping_position.is_reverse(),
                                       final_mapping_position.offset() + mapping_from_length(final_mapping) - 1);
            
            // the longest gap that could be detected at this position in the read
            size_t src_max_gap = get_aligner()->longest_detectable_gap(alignment, src_match_node.end);
            
            unordered_set<pair<size_t, size_t>> edges_for_removal;
            
            for (const pair<size_t, size_t>& edge : src_match_node.edges) {
                ExactMatchNode& dest_match_node = multi_aln_graph.match_nodes[edge.first];
                pos_t dest_pos = make_pos_t(multipath_aln_out.subpath(edge.first).path().mapping(0).position());
                
#ifdef debug_multipath_mapper_alignment
                cerr << "forming intervening alignment for edge to node " << edge.first << endl;
#endif
                
                size_t intervening_length = dest_match_node.begin - src_match_node.end;
                size_t max_dist = intervening_length +
                    std::min(src_max_gap, get_aligner()->longest_detectable_gap(alignment, dest_match_node.begin)) + 1;
                
#ifdef debug_multipath_mapper_alignment
                cerr << "read dist: " << intervening_length << ", source max gap: "
                     << src_max_gap << ", dest max gap "
                     << get_aligner()->longest_detectable_gap(alignment, dest_match_node.begin) << endl;
#endif
                
                // extract the graph between the matches
                Graph connecting_graph;
                unordered_map<id_t, id_t> connect_trans = algorithms::extract_connecting_graph(&align_graph,     // DAG with split strands
                                                                                               connecting_graph, // graph to extract into
                                                                                               max_dist,         // longest distance necessary
                                                                                               src_pos,          // end of earlier match
                                                                                               dest_pos,         // beginning of later match
                                                                                               false,            // do not extract the end positions in the matches
                                                                                               false,            // do not bother finding all cycles (it's a DAG)
                                                                                               true,             // remove tips
                                                                                               true,             // only include nodes on connecting paths
                                                                                               true);            // enforce max distance strictly
                
                
                if (connecting_graph.node_size() == 0) {
                    // the MEMs weren't connectable with a positive score after all, mark the edge for removal
                    edges_for_removal.insert(edge);
                    continue;
                }
                
                // transfer the substring between the matches to a new alignment
                Alignment intervening_sequence;
                intervening_sequence.set_sequence(alignment.sequence().substr(src_match_node.end - alignment.sequence().begin(),
                                                                              dest_match_node.begin - src_match_node.end));
                if (!alignment.quality().empty()) {
                    intervening_sequence.set_quality(alignment.quality().substr(src_match_node.end - alignment.sequence().begin(),
                                                                                dest_match_node.begin - src_match_node.end));
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << "aligning sequence " << intervening_sequence.sequence() << " to connecting graph: " << pb2json(connecting_graph) << endl;
#endif
                
                bool added_direct_connection = false;
                // TODO a better way of choosing the number of alternate alignments
                // TODO alternate alignments restricted only to distinct node paths?
                vector<Alignment> alt_alignments;
                get_aligner()->align_global_banded_multi(intervening_sequence, alt_alignments,
                                                         connecting_graph, num_alt_alns, band_padding, true);
                
                for (Alignment& connecting_alignment : alt_alignments) {
#ifdef debug_multipath_mapper_alignment
                    cerr << "translating connecting alignment: " << pb2json(connecting_alignment) << endl;
#endif
                    
                    const Path& aligned_path = connecting_alignment.path();
                    const Mapping& first_mapping = aligned_path.mapping(0);
                    const Mapping& last_mapping = aligned_path.mapping(aligned_path.mapping_size() - 1);
                    
                    bool add_first_mapping = mapping_from_length(first_mapping) != 0 || mapping_to_length(first_mapping) != 0;
                    bool add_last_mapping = ((mapping_from_length(last_mapping) != 0 || mapping_to_length(last_mapping) != 0)
                                             && aligned_path.mapping_size() > 1);
                    
                    if (!(add_first_mapping || add_last_mapping) && aligned_path.mapping_size() <= 2) {
                        if (!added_direct_connection) {
                            // edge case where there is a simple split but other non-simple edges intersect the target
                            // at the same place (so it passes the previous filter)
                            // it actually doesn't need an alignment, just a connecting edge
                            src_subpath->add_next(edge.first);
                            added_direct_connection = true;
                        }
                        continue;
                    }
                    
                    // create a subpath between the matches for this alignment
                    Subpath* connecting_subpath = multipath_aln_out.add_subpath();
                    connecting_subpath->set_score(connecting_alignment.score());
                    Path* subpath_path = connecting_subpath->mutable_path();
                    
                    int32_t rank = 1;
                    
                    // check to make sure the first is not an empty anchoring mapping
                    if (add_first_mapping) {
                        Mapping* mapping = subpath_path->add_mapping();
                        *mapping = first_mapping;
                        mapping->set_rank(rank);
#ifdef debug_multipath_mapper_alignment
                        cerr << "first mapping is not empty, formed mapping: " << pb2json(*mapping) << endl;
#endif
                        rank++;
                    }
                    // add all mapping in between the ends
                    for (size_t j = 1; j < aligned_path.mapping_size() - 1; j++) {
                        Mapping* mapping = subpath_path->add_mapping();
                        *mapping = aligned_path.mapping(j);
                        mapping->set_rank(rank);
#ifdef debug_multipath_mapper_alignment
                        cerr << "added middle mapping: " << pb2json(*mapping) << endl;
#endif
                        rank++;
                    }
                    // check to make sure the last is not an empty anchoring mapping or the same as the first
                    if (add_last_mapping) {
                        Mapping* mapping = subpath_path->add_mapping();
                        *mapping = last_mapping;
                        mapping->set_rank(rank);
#ifdef debug_multipath_mapper_alignment
                        cerr << "final mapping is not empty, formed mapping: " << pb2json(*mapping) << endl;
#endif
                    }
                
                    // add the appropriate connections
                    src_subpath->add_next(multipath_aln_out.subpath_size() - 1);
                    connecting_subpath->add_next(edge.first);
                    
                    // translate the path into the space of the main graph unless the path is null
                    if (connecting_subpath->path().mapping_size() != 0) {
                        translate_node_ids(*connecting_subpath->mutable_path(), connect_trans);
                        Mapping* first_subpath_mapping = connecting_subpath->mutable_path()->mutable_mapping(0);
                        if (first_subpath_mapping->position().node_id() == final_mapping.position().node_id()) {
                            first_subpath_mapping->mutable_position()->set_offset(offset(src_pos) + 1);
                        }
                    }
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "subpath from " << j << " to " << edge.first << ":" << endl;
                    cerr << pb2json(*connecting_subpath) << endl;
#endif
                }
            }
            
            if (!edges_for_removal.empty()) {
                auto new_end = std::remove_if(src_match_node.edges.begin(), src_match_node.edges.end(),
                                              [&](const pair<size_t, size_t>& edge) {
                                                  return edges_for_removal.count(edge);
                                              });
                src_match_node.edges.resize(new_end - src_match_node.edges.begin());
            }
        }
        
        // function to reorder the nodes of a Protobuf graph in topological order, flip doubly reversing edges,
        // and remove empty sequence nodes (invariants required for gssw alignment)
        // TODO: this is duplicative with VG::sort, but I don't want to construct a VG here
        auto groom_graph_for_gssw = [](Graph& graph) {
            // remove empty nodes
            size_t end = graph.node_size();
            size_t idx = 0;
            unordered_set<id_t> removed_nodes;
            while (idx < end) {
                if (graph.node(idx).sequence().empty()) {
                    end--;
                    removed_nodes.insert(graph.node(idx).id());
                    std::swap(*graph.mutable_node(idx), *graph.mutable_node(end));
                }
                else {
                    idx++;
                }
            }
            if (end != graph.node_size()) {
                graph.mutable_node()->DeleteSubrange(end, graph.node_size() - end);
                
                // look for any edges connecting them and remove these too
                end = graph.edge_size();
                idx = 0;
                while (idx < end) {
                    Edge* edge = graph.mutable_edge(idx);
                    if (removed_nodes.count(edge->from()) || removed_nodes.count(edge->to())) {
                        end--;
                        std::swap(*edge, *graph.mutable_edge(end));
                    }
                    else {
                        idx++;
                    }
                }
                graph.mutable_edge()->DeleteSubrange(end, graph.edge_size() - end);
            }
            
            // flip doubly reversing edges
            for (size_t i = 0; i < graph.edge_size(); i++) {
                Edge* edge = graph.mutable_edge(i);
                if (edge->from_start() && edge->to_end()) {
                    id_t tmp = edge->from();
                    edge->set_from(edge->to());
                    edge->set_to(tmp);
                    edge->set_from_start(false);
                    edge->set_to_end(false);
                }
            }
            
            // associate node ids with their index
            unordered_map<id_t, size_t> node_idx;
            for (size_t i = 0; i < graph.node_size(); i++) {
                node_idx[graph.node(i).id()] = i;
            }
            
            // construct adjacency list and compute in degrees
            vector<size_t> in_degree(graph.node_size(), 0);
            vector<vector<size_t>> adj_list(graph.node_size());
            for (size_t i = 0; i < graph.edge_size(); i++) {
                const Edge& edge = graph.edge(i);
                size_t to_idx = node_idx[edge.to()];
                adj_list[node_idx[edge.from()]].push_back(to_idx);
                in_degree[to_idx]++;
            }
            
            // get the topological ordering of the graph (Kahn's algorithm)
            vector<size_t> source_stack;
            for (size_t i = 0; i < graph.node_size(); i++) {
                if (in_degree[i] == 0) {
                    source_stack.push_back(i);
                }
            }
            
            vector<id_t> order(graph.node_size());
            size_t next = 0;
            while (!source_stack.empty()) {
                size_t src = source_stack.back();
                source_stack.pop_back();
                
                for (size_t dest : adj_list[src]) {
                    in_degree[dest]--;
                    if (in_degree[dest] == 0) {
                        source_stack.push_back(dest);
                    }
                }
                
                order[next] = src;
                next++;
            }
            
            // identify the index that we want each node to end up at
            vector<size_t> index(order.size());
            for (size_t i = 0; i < order.size(); i++) {
                index[order[i]] = i;
            }

            // in place permutation according to the topological order
            for (size_t i = 0; i < graph.node_size(); i++) {
                while (index[i] != i) {
                    std::swap(*graph.mutable_node(i), *graph.mutable_node(index[i]));
                    std::swap(index[i], index[index[i]]);
                }
            }
        };
        
        vector<bool> is_source_node(multi_aln_graph.match_nodes.size(), true);
        for (size_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
            ExactMatchNode& match_node = multi_aln_graph.match_nodes[j];
            if (match_node.edges.empty()) {
                if (match_node.end != alignment.sequence().end()) {
                    
                    Subpath* sink_subpath = multipath_aln_out.mutable_subpath(j);
                    
                    int64_t target_length = (alignment.sequence().end() - match_node.end) +
                                            get_aligner()->longest_detectable_gap(alignment, match_node.end);
                    pos_t end_pos = final_position(match_node.path);
                    // want past-the-last instead of last index here
                    get_offset(end_pos)++;
                    
                    Graph tail_graph;
                    unordered_map<id_t, id_t> tail_trans = algorithms::extract_extending_graph(&align_graph,
                                                                                               tail_graph,
                                                                                               target_length,
                                                                                               end_pos,
                                                                                               false,         // search forward
                                                                                               false);        // no need to preserve cycles (in a DAG)
                    
                    // ensure invariants that gssw-based alignment expects
                    groom_graph_for_gssw(tail_graph);
                    
                    // get the sequence remaining in the right tail
                    Alignment right_tail_sequence;
                    right_tail_sequence.set_sequence(alignment.sequence().substr(match_node.end - alignment.sequence().begin(),
                                                                                 alignment.sequence().end() - match_node.end));
                    if (!alignment.quality().empty()) {
                        right_tail_sequence.set_quality(alignment.quality().substr(match_node.end - alignment.sequence().begin(),
                                                                                   alignment.sequence().end() - match_node.end));
                    }
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "aligning sequence: " << right_tail_sequence.sequence() << endl << "to right tail graph: " << pb2json(tail_graph) << endl;
#endif
                    
                    vector<Alignment> alt_alignments;
                    if (tail_graph.node_size() == 0) {
                        // edge case for when a read keeps going past the end of a graph
                        alt_alignments.emplace_back();
                        Alignment& tail_alignment = alt_alignments.back();
                        tail_alignment.set_score(get_aligner()->score_gap(right_tail_sequence.sequence().size()));
                        Mapping* insert_mapping = tail_alignment.mutable_path()->add_mapping();
                        
                        // add a soft clip
                        Edit* edit = insert_mapping->add_edit();
                        edit->set_to_length(right_tail_sequence.sequence().size());
                        edit->set_sequence(right_tail_sequence.sequence());
                        
                        // make it at the correct position
                        const Path& anchoring_path = multi_aln_graph.match_nodes[j].path;
                        const Mapping& anchoring_mapping = anchoring_path.mapping(anchoring_path.mapping_size() - 1);
                        Position* anchoring_position = insert_mapping->mutable_position();
                        anchoring_position->set_node_id(anchoring_mapping.position().node_id());
                        anchoring_position->set_is_reverse(anchoring_mapping.position().is_reverse());
                        anchoring_position->set_offset(anchoring_mapping.position().offset() + mapping_from_length(anchoring_mapping));
#ifdef debug_multipath_mapper_alignment
                        cerr << "read overhangs end of graph, manually added softclip: " << pb2json(tail_alignment) << endl;
#endif
                        // the ID translator is empty, so add this ID here so it doesn't give an out of index error
                        id_t node_id = insert_mapping->position().node_id();
                        tail_trans[node_id] = node_id;
                    }
                    else {
                        // align against the graph
                        get_aligner()->align_pinned_multi(right_tail_sequence, alt_alignments, tail_graph, true, num_alt_alns);
                    }
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "made " << alt_alignments.size() << " tail alignments" << endl;
#endif
                    
                    const Mapping& final_mapping = match_node.path.mapping(match_node.path.mapping_size() - 1);
                    
                    for (Alignment& tail_alignment : alt_alignments) {
                        
                        sink_subpath->add_next(multipath_aln_out.subpath_size());
                        
                        Subpath* tail_subpath = multipath_aln_out.add_subpath();
                        *tail_subpath->mutable_path() = tail_alignment.path();
                        tail_subpath->set_score(tail_alignment.score());
                        
                        translate_node_ids(*tail_subpath->mutable_path(), tail_trans);
                        Mapping* first_mapping = tail_subpath->mutable_path()->mutable_mapping(0);
                        if (first_mapping->position().node_id() == final_mapping.position().node_id()) {
                            first_mapping->mutable_position()->set_offset(offset(end_pos));
                        }
                        else if (tail_alignment.path().mapping_size() == 1 && first_mapping->edit_size() == 1
                                 && first_mapping->edit(0).from_length() == 0 && first_mapping->edit(0).to_length() > 0
                                 && first_mapping->position().node_id() != final_mapping.position().node_id()) {
                            // this is a pure soft-clip on the beginning of the next node, we'll move it to the end
                            // of the match node to match invariants expected by other parts of the code base
                            Position* pos = first_mapping->mutable_position();
                            pos->set_node_id(final_mapping.position().node_id());
                            pos->set_is_reverse(final_mapping.position().is_reverse());
                            pos->set_offset(final_mapping.position().offset() + mapping_from_length(final_mapping));
                        }
#ifdef debug_multipath_mapper_alignment
                        cerr << "subpath from " << j << " to right tail:" << endl;
                        cerr << pb2json(*tail_subpath) << endl;
#endif
                    }
                }
            }
            else {
                for (const pair<size_t, size_t>& edge : match_node.edges) {
                    is_source_node[edge.first] = false;
                }
            }
        }
        
        for (size_t j = 0; j < multi_aln_graph.match_nodes.size(); j++) {
            if (is_source_node[j]) {
                ExactMatchNode& match_node = multi_aln_graph.match_nodes[j];
                if (match_node.begin != alignment.sequence().begin()) {
                    
                    int64_t target_length = (match_node.begin - alignment.sequence().begin()) +
                                            get_aligner()->longest_detectable_gap(alignment, match_node.begin);
                    pos_t begin_pos = initial_position(match_node.path);

                    
                    Graph tail_graph;
                    unordered_map<id_t, id_t> tail_trans = algorithms::extract_extending_graph(&align_graph,
                                                                                               tail_graph,
                                                                                               target_length,
                                                                                               begin_pos,
                                                                                               true,          // search backward
                                                                                               false);        // no need to preserve cycles (in a DAG)
                    
                    // ensure invariants that gssw-based alignment expects
                    groom_graph_for_gssw(tail_graph);
                    
                    Alignment left_tail_sequence;
                    left_tail_sequence.set_sequence(alignment.sequence().substr(0, match_node.begin - alignment.sequence().begin()));
                    if (!alignment.quality().empty()) {
                        left_tail_sequence.set_quality(alignment.quality().substr(0, match_node.begin - alignment.sequence().begin()));
                    }
                    
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "aligning sequence: " << left_tail_sequence.sequence() << endl << "to left tail graph: " << pb2json(tail_graph) << endl;
#endif
                    vector<Alignment> alt_alignments;
                    if (tail_graph.node_size() == 0) {
                        // edge case for when a read keeps going past the end of a graph
                        alt_alignments.emplace_back();
                        Alignment& tail_alignment = alt_alignments.back();
                        tail_alignment.set_score(get_aligner()->score_gap(left_tail_sequence.sequence().size()));
                        Mapping* insert_mapping = tail_alignment.mutable_path()->add_mapping();
                        
                        // add a soft clip
                        Edit* edit = insert_mapping->add_edit();
                        edit->set_to_length(left_tail_sequence.sequence().size());
                        edit->set_sequence(left_tail_sequence.sequence());
                        
                        // make it at the correct position
                        *insert_mapping->mutable_position() = multi_aln_graph.match_nodes[j].path.mapping(0).position();
#ifdef debug_multipath_mapper_alignment
                        cerr << "read overhangs end of graph, manually added softclip: " << pb2json(tail_alignment) << endl;
#endif
                        // the ID translator is empty, so add this ID here so it doesn't give an out of index error
                        id_t node_id = insert_mapping->position().node_id();
                        tail_trans[node_id] = node_id;
                    }
                    else {
                        get_aligner()->align_pinned_multi(left_tail_sequence, alt_alignments, tail_graph, false, num_alt_alns);
                    }
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "made " << alt_alignments.size() << " tail alignments" << endl;
#endif
                    const Mapping& first_mapping = match_node.path.mapping(0);
                    for (Alignment& tail_alignment : alt_alignments) {
                        Subpath* tail_subpath = multipath_aln_out.add_subpath();
                        *tail_subpath->mutable_path() = tail_alignment.path();
                        tail_subpath->set_score(tail_alignment.score());
                        
                        tail_subpath->add_next(j);
                        multipath_aln_out.add_start(multipath_aln_out.subpath_size() - 1);
                        
                        translate_node_ids(*tail_subpath->mutable_path(), tail_trans);
#ifdef debug_multipath_mapper_alignment
                        cerr << "subpath from " << j << " to left tail:" << endl;
                        cerr << pb2json(*tail_subpath) << endl;
#endif
                        Mapping* final_mapping = tail_subpath->mutable_path()->mutable_mapping(tail_subpath->path().mapping_size() - 1);
                        if (tail_subpath->path().mapping_size() == 1 && final_mapping->edit_size() == 1
                            && final_mapping->edit(0).from_length() == 0 && final_mapping->edit(0).to_length() > 0
                            && final_mapping->position().node_id() != first_mapping.position().node_id()) {
                            // this is a pure soft clip on the end of the previous node, so we move it over to the
                            // beginning of the match node to match invariants in rest of code base
                            *final_mapping->mutable_position() = first_mapping.position();
                        }
                    }
                }
                else {
                    multipath_aln_out.add_start(j);
                }
            }
        }
        
#ifdef debug_multipath_mapper_alignment
        cerr << "multipath alignment before translation: " << pb2json(multipath_aln_out) << endl;
#endif
        for (size_t j = 0; j < multipath_aln_out.subpath_size(); j++) {
            translate_oriented_node_ids(*multipath_aln_out.mutable_subpath(j)->mutable_path(), node_trans);
        }
#ifdef debug_multipath_mapper_alignment
        cerr << "completed multipath alignment: " << pb2json(multipath_aln_out) << endl;
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
    
    void MultipathMapper::strip_full_length_bonuses(MultipathAlignment& mulipath_aln) const {
        
        int32_t full_length_bonus = get_aligner()->full_length_bonus;
        // strip bonus from source paths
        if (mulipath_aln.start_size()) {
            // use the precomputed list of sources if we have it
            for (size_t i = 0; i < mulipath_aln.start_size(); i++) {
                Subpath* source_subpath = mulipath_aln.mutable_subpath(mulipath_aln.start(i));
                if (edit_is_insertion(source_subpath->path().mapping(0).edit(0))) {
                    source_subpath->set_score(source_subpath->score() - full_length_bonus);
                }
            }
        }
        else {
            // find sources
            vector<bool> is_source(mulipath_aln.subpath_size(), true);
            for (size_t i = 0; i < mulipath_aln.subpath_size(); i++) {
                const Subpath& subpath = mulipath_aln.subpath(i);
                for (size_t j = 0; j < subpath.next_size(); j++) {
                    is_source[subpath.next(j)] = false;
                }
            }
            // strip the bonus from the sources
            for (size_t i = 0; i < mulipath_aln.subpath_size(); i++) {
                if (!is_source[i]) {
                    continue;
                }
                Subpath* source_subpath = mulipath_aln.mutable_subpath(i);
                if (edit_is_insertion(source_subpath->path().mapping(0).edit(0))) {
                    source_subpath->set_score(source_subpath->score() - full_length_bonus);
                }
            }
        }
        // strip bonus from sink paths
        for (size_t i = 0; i < mulipath_aln.subpath_size(); i++) {
            Subpath* subpath = mulipath_aln.mutable_subpath(i);
            if (subpath->next_size() == 0) {
                const Mapping& final_mapping = subpath->path().mapping(subpath->path().mapping_size() - 1);
                if (edit_is_insertion(final_mapping.edit(final_mapping.edit_size() - 1))) {
                    subpath->set_score(subpath->score() - full_length_bonus);
                }
            }
        }
    }
    
    void MultipathMapper::sort_and_compute_mapping_quality(vector<MultipathAlignment>& multipath_alns,
                                                           MappingQualityMethod mapq_method) const {
        if (multipath_alns.empty()) {
            return;
        }
        
        // only do the population MAPQ if it might disambiguate two paths (since it's not
        // as cheap as just using the score)
        bool include_population_component = (use_population_mapqs && multipath_alns.size() > 1);
        // records whether all of the pathsenumerated across all multipath alignments followed the edges in the index
        bool all_paths_pop_consistent = true;
        
        double log_base = get_aligner()->log_base;
        
        // The score of the optimal Alignment for each MultipathAlignment, not adjusted for population
        vector<double> base_scores(multipath_alns.size(), 0.0);
        // The scores of the best Alignment for each MultipathAlignment, adjusted for population.
        // These can be negativem but will be bumped up to all be positive later.
        vector<double> pop_adjusted_scores;
        if (include_population_component) {
            pop_adjusted_scores.resize(multipath_alns.size());
        }
        
        // We need to track the min population score by itself (not converted
        // to the alignment score log base) so we can adjust the
        // pop_adjusted_scores up, turning the largest penalty into a 0 bonus.
        double min_pop_score = numeric_limits<double>::max();

        
        for (size_t i = 0; i < multipath_alns.size(); i++) {
            // Score all the multipath alignment candidates, optionally using
            // population adjustment
            
            // We will query the population database for this alignment if it
            // is turned on and it succeeded for the others.
            bool query_population = include_population_component && all_paths_pop_consistent;
            
            // Generate the top alignment, or the top population_max_paths
            // alignments if we are doing multiple alignments for population
            // scoring.
            auto wanted_alignments = query_population ? population_max_paths : 1;
            auto alignments = optimal_alignments(multipath_alns[i], wanted_alignments);
            assert(!alignments.empty());
            
#ifdef debug_multipath_mapper_mapping
            cerr << "Got " << alignments.size() << " / " << wanted_alignments << " tracebacks for multipath " << i << endl;
#endif
#ifdef debug_multipath_mapper_alignment
            cerr << pb2json(multipath_alns[i]) << endl;
#endif
           
            // Collect the score of the optimal alignment, to use if population
            // scoring fails for a multipath alignment. Put it in the optimal
            // base score.
            base_scores[i] = alignments[0].score();
            
            if (query_population) {
                
                // Make sure to grab the memo
                auto& memo = get_rr_memo(recombination_penalty, xindex->get_haplotype_count());
                
                // Now compute population scores for all the top paths
                vector<double> alignment_pop_scores(alignments.size(), 0.0);
                for (size_t j = 0; j < alignments.size(); j++) {
                    // Score each alignment if possible
                    auto pop_score = haplo_score_provider->score(alignments[j].path(), memo);
                    
#ifdef debug_multipath_mapper_mapping
                    cerr << "Got pop score " << pop_score.first << ", " << pop_score.second << " for alignment " << j
                        << " score " << alignments[j].score() << " of multipath " << i << endl;
#endif
#ifdef debug_multipath_mapper_alignment
                    cerr << pb2json(alignments[j]) << endl;
#endif
                   
                    alignment_pop_scores[j] = pop_score.first;
                    
                    all_paths_pop_consistent &= pop_score.second;
                }
                
                if (!all_paths_pop_consistent) {
                    // If we failed, bail out on population score correction for the whole MultipathAlignment.
                    
                    // Go and do the next MultipathAlignment since we have the base score for this one
                    continue;
                }
                
                // Otherwise, pick the best adjusted score and its pop score
                pop_adjusted_scores[i] = numeric_limits<double>::min();
                double best_alignment_pop_score;
                for (size_t j = 0; j < alignments.size(); j++) {
                    // Compute the adjusted score for each alignment
                    auto adjusted_score = alignments[j].score() + alignment_pop_scores[j] / log_base;
                    
                    if (adjusted_score > pop_adjusted_scores[i]) {
                        // It is the best, so use it.
                        // TODO: somehow know we want this Alignment when collapsing the MultipathAlignment later.
                        pop_adjusted_scores[i] = adjusted_score;
                        best_alignment_pop_score = alignment_pop_scores[j];
                    }
                }
                
                // Save the actual population score for it in the min pop score, if it is lower.
                min_pop_score = min(min_pop_score, best_alignment_pop_score);
            }
        }
        
        if (include_population_component && all_paths_pop_consistent) {
            for (auto& score : pop_adjusted_scores) {
                // Adjust the adjusted scores up/down by the contributoon of the min pop score to ensure no scores are negative
                score -= min_pop_score / log_base;
            }
        }
        
        // Select whether to use base or adjusted scores depending on whether
        // we did population-aware alignment and succeeded for all the
        // multipath alignments.
        auto& scores = (include_population_component && all_paths_pop_consistent) ? pop_adjusted_scores : base_scores;
        
        // find the order of the scores
        vector<size_t> order(multipath_alns.size(), 0);
        for (size_t i = 1; i < multipath_alns.size(); i++) {
            order[i] = i;
        }
        sort(order.begin(), order.end(), [&](const size_t i, const size_t j) { return scores[i] > scores[j]; });
        
        // translate the order to an index
        vector<size_t> index(multipath_alns.size());
        for (size_t i = 0; i < multipath_alns.size(); i++) {
            index[order[i]] = i;
        }
        
        // put the scores and alignments in order
        for (size_t i = 0; i < multipath_alns.size(); i++) {
            while (index[i] != i) {
                std::swap(scores[index[i]], scores[i]);
                std::swap(multipath_alns[index[i]], multipath_alns[i]);
                std::swap(index[index[i]], index[i]);
                
            }
        }
        
#ifdef debug_multipath_mapper_mapping
        cerr << "scores obtained of multi-mappings:" << endl;
        for (auto score : scores) {
            cerr << "\t" << score << endl;
        }
#endif
        
        if (mapq_method != None) {
            int32_t raw_mapq;
            if (mapq_method == Adaptive) {
                raw_mapq = get_aligner()->compute_mapping_quality(scores, scores.size() < 2 ? true :
                                                                  (scores[1] < scores[0] - get_aligner()->mapping_quality_score_diff(max_mapping_quality)));
            }
            else {
                raw_mapq = get_aligner()->compute_mapping_quality(scores, mapq_method == Approx);
            }
            
            // arbitrary scaling, seems to help performance
            raw_mapq *= mapq_scaling_factor;
            
#ifdef debug_multipath_mapper_mapping
            cerr << "scores yield a raw MAPQ of " << raw_mapq << endl;
#endif
            multipath_alns.front().set_mapping_quality(min<int32_t>(raw_mapq, max_mapping_quality));
        }
    }
    
    // TODO: pretty duplicative with the unpaired version
    void MultipathMapper::sort_and_compute_mapping_quality(vector<pair<MultipathAlignment, MultipathAlignment>>& multipath_aln_pairs,
                                                           vector<pair<pair<size_t, size_t>, int64_t>>& cluster_pairs) const {
        assert(multipath_aln_pairs.size() == cluster_pairs.size());
        
        if (multipath_aln_pairs.empty()) {
            return;
        }
        
        // only do the population MAPQ if it might disambiguate two paths (since it's not
        // as cheap as just using the score)
        bool include_population_component = (use_population_mapqs && multipath_aln_pairs.size() > 1);
        // records whether of the paths followed the edges in the index
        bool all_paths_pop_consistent = true;
        
        double log_base = get_aligner()->log_base;
        
        // the scores of the optimal alignments and fragments, ignoring population
        vector<double> base_scores(multipath_aln_pairs.size(), 0.0);
        
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
            pair<MultipathAlignment, MultipathAlignment>& multipath_aln_pair = multipath_aln_pairs[i];
            
            // We will query the population database for this alignment if it
            // is turned on and it succeeded for the others.
            bool query_population = include_population_component && all_paths_pop_consistent;
            
            // Generate the top alignments on each side, or the top
            // population_max_paths alignments if we are doing multiple
            // alignments for population scoring.
            auto alignments1 = optimal_alignments(multipath_aln_pair.first, query_population ? population_max_paths : 1);
            auto alignments2 = optimal_alignments(multipath_aln_pair.second, query_population ? population_max_paths : 1);
            assert(!alignments1.empty());
            assert(!alignments2.empty());
            
            // Compute the optimal alignment score ignoring population
            int32_t alignment_score = alignments1[0].score() + alignments2[0].score();
            
            // compute the fragment distribution's contribution to the score
            double frag_score = fragment_length_log_likelihood(cluster_pairs[i].second) / log_base;
            min_frag_score = min(frag_score, min_frag_score);
            
            // Record the base score, including fragment contribution
            base_scores[i] = alignment_score + frag_score;
            
            if (query_population) {
                // We also want to select the optimal population-scored alignment on each side and compute a pop-adjusted score.
                
                // Make sure to grab the memo
                auto& memo = get_rr_memo(recombination_penalty, xindex->get_haplotype_count());
                
                // What's the population score for each alignment?
                vector<double> pop_scores1(alignments1.size());
                vector<double> pop_scores2(alignments2.size());
                
                for (size_t j = 0; j < alignments1.size(); j++) {
                    // Pop score the first alignments
                    auto pop_score = haplo_score_provider->score(alignments1[j].path(), memo);
                    pop_scores1[j] = pop_score.first;
                    all_paths_pop_consistent &= pop_score.second;
                }
                
                for (size_t j = 0; j < alignments2.size(); j++) {
                    // Pop score the second alignments
                    auto pop_score = haplo_score_provider->score(alignments2[j].path(), memo);
                    pop_scores2[j] = pop_score.first;
                    all_paths_pop_consistent &= pop_score.second;
                }
                
                if (!all_paths_pop_consistent) {
                    // If we couldn't score everything, bail
                    continue;
                }
                
                // Pick the best alignment on each side
                auto best_index1 = max_element(pop_scores1.begin(), pop_scores1.end()) - pop_scores1.begin();
                auto best_index2 = max_element(pop_scores2.begin(), pop_scores2.end()) - pop_scores2.begin();
                
                // And compute base score for those alignments and pop adjustment.
                // TODO: Synchronize with single-end case about whether we store things base-adjusted or not!
                alignment_score = alignments1[0].score() + alignments2[0].score();
                double pop_score = (pop_scores1[best_index1] + pop_scores2[best_index2]) / log_base;
                
                // Compute the total pop adjusted score for this MultipathAlignment
                pop_adjusted_scores[i] = alignment_score + frag_score + pop_score;
                
                // Record our extra score if it was a new minimum
                min_extra_score = min(frag_score + pop_score, min_extra_score);
            }
        }
        
        // Decide which scores to use depending on whether we have pop adjusted scores we want to use
        auto& scores = (include_population_component && all_paths_pop_consistent) ? pop_adjusted_scores : base_scores;
        
        for (auto& score : scores) {
            // Pull the min frag or extra score out of the score so it will be nonnegative
            score -= (include_population_component && all_paths_pop_consistent) ? min_extra_score : min_frag_score;
        }
        
        // find the order of the scores
        vector<size_t> order(multipath_aln_pairs.size(), 0);
        for (size_t i = 1; i < multipath_aln_pairs.size(); i++) {
            order[i] = i;
        }
        sort(order.begin(), order.end(), [&](const size_t i, const size_t j) { return scores[i] > scores[j]; });
        
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
                std::swap(index[index[i]], index[i]);
                
            }
        }
        
#ifdef debug_multipath_mapper_mapping
        cerr << "scores and distances obtained of multi-mappings:" << endl;
        for (int i = 0; i < multipath_aln_pairs.size(); i++) {
            cerr << "\talign:" << optimal_alignment_score(multipath_aln_pairs[i].first) + optimal_alignment_score(multipath_aln_pairs[i].second) << ", length: " << cluster_pairs[i].second << ", combined: " << scores[i] << endl;
        }
#endif
        
        if (mapping_quality_method != None) {
            int32_t raw_mapq;
            if (mapping_quality_method == Adaptive) {
                raw_mapq = get_aligner()->compute_mapping_quality(scores, scores.size() < 2 ? true :
                                                                  (scores[1] < scores[0] - get_aligner()->mapping_quality_score_diff(max_mapping_quality)));
            }
            else {
                raw_mapq = get_aligner()->compute_mapping_quality(scores, mapping_quality_method == Approx);
            }
            
            // arbitrary scaling, seems to help performance
            raw_mapq *= mapq_scaling_factor;
            
#ifdef debug_multipath_mapper_mapping
            cerr << "scores yield a raw MAPQ of " << raw_mapq << endl;
#endif
            int32_t mapq = min<int32_t>(raw_mapq, max_mapping_quality);
            multipath_aln_pairs.front().first.set_mapping_quality(mapq);
            multipath_aln_pairs.front().second.set_mapping_quality(mapq);
            
            if (multipath_aln_pairs.size() > 1) {
                // find the duplicates of the optimal pair (initially mark with only the pair itself)
                vector<size_t> duplicates_1(1, 0);
                vector<size_t> duplicates_2(1, 0);
                vector<size_t> to_remove;
                for (size_t i = 1; i < multipath_aln_pairs.size(); i++) {
                    bool duplicate_1 = share_start_position(multipath_aln_pairs[0].first, multipath_aln_pairs[i].first);
                    bool duplicate_2 = share_start_position(multipath_aln_pairs[0].second, multipath_aln_pairs[i].second);
                    if (duplicate_1 && duplicate_2) {
#ifdef debug_multipath_mapper_mapping
                        cerr << "found double end duplication at index " << i << endl;
#endif
                        
                        // this pair is a complete duplication (not just one end) we want it gone
                        to_remove.push_back(i);
                    }
                    else if (duplicate_1) {
#ifdef debug_multipath_mapper_mapping
                        cerr << "found left end duplication at index " << i << endl;
#endif
                        duplicates_1.push_back(i);
                    }
                    else if (duplicate_2) {
#ifdef debug_multipath_mapper_mapping
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
                            scores[i - removed_so_far] = move(scores[i]);
                            cluster_pairs[i - removed_so_far] = move(cluster_pairs[i]);
                        }
                    }
                    
                    // remove the end positions that are now empty
                    multipath_aln_pairs.resize(multipath_aln_pairs.size() - to_remove.size());
                    scores.resize(scores.size() - to_remove.size());
                    cluster_pairs.resize(cluster_pairs.size() - to_remove.size());
                    
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
                    int32_t raw_mapq_1 = get_aligner()->compute_group_mapping_quality(scores, duplicates_1);
                    int32_t raw_mapq_2 = get_aligner()->compute_group_mapping_quality(scores, duplicates_2);
                    
                    // arbitrary scaling, seems to help performance
                    raw_mapq_1 *= mapq_scaling_factor;
                    raw_mapq_2 *= mapq_scaling_factor;
                    
#ifdef debug_multipath_mapper_mapping
                    cerr << "deduplicated raw MAPQs are " << raw_mapq_1 << " and " << raw_mapq_2 << endl;
#endif
                    
                    int32_t mapq_1 = min<int32_t>(raw_mapq_1, max_mapping_quality);
                    int32_t mapq_2 = min<int32_t>(raw_mapq_2, max_mapping_quality);
                    
                    multipath_aln_pairs.front().first.set_mapping_quality(mapq_1);
                    multipath_aln_pairs.front().second.set_mapping_quality(mapq_2);
                }
            }
        }
    }
            
    double MultipathMapper::fragment_length_log_likelihood(int64_t length) const {
        double dev = length - fragment_length_distr.mean();
        return -dev * dev / (2.0 * fragment_length_distr.stdev() * fragment_length_distr.stdev());
    }
    
    void MultipathMapper::set_automatic_min_clustering_length(double random_mem_probability) {
        min_clustering_mem_length = max<int>(log(1.0 - pow(random_mem_probability, 1.0 / xindex->seq_length)) / log(0.25), 1);
    }
    
    bool MultipathMapper::validate_multipath_alignment(const MultipathAlignment& multipath_aln) const {
        
        // are the subpaths in topological order?
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            for (size_t j = 0; j < subpath.next_size(); j++) {
                if (subpath.next(j) <= i) {
#ifdef debug_multipath_mapper_mapping
                    cerr << "validation failure on topological order" << endl;
#endif
                    return false;
                }
            }
        }
        
        // are the start subpaths properly labeled (if they are included)?
        
        if (multipath_aln.start_size()) {
            vector<bool> is_source(multipath_aln.subpath_size(), true);
            for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
                const Subpath& subpath = multipath_aln.subpath(i);
                for (size_t j = 0; j < subpath.next_size(); j++) {
                    is_source[subpath.next(j)] = false;
                }
            }
            
            size_t num_starts = 0;
            for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
                num_starts += is_source[i];
            }
            
            if (num_starts != multipath_aln.start_size()) {
#ifdef debug_multipath_mapper_mapping
                cerr << "validation failure on correct number of starts" << endl;
                for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
                    if (is_source[i]) {
                        cerr << i << " ";
                    }
                }
                cerr << endl;
#endif
                return false;
            }
            
            for (size_t i = 0; i < multipath_aln.start_size(); i++) {
                if (!is_source[multipath_aln.start(i)]) {
#ifdef debug_multipath_mapper_mapping
                    cerr << "validation failure on correctly identified starts" << endl;
                    for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
                        if (is_source[i]) {
                            cerr << i << " ";
                        }
                        cerr << endl;
                    }
#endif
                    return false;
                }
            }
        }
        
        // are the subpaths contiguous along the read?
        
        vector<pair<int64_t, int64_t>> subpath_read_interval(multipath_aln.subpath_size(), make_pair<int64_t, int64_t>(-1, -1));
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            
            if (subpath_read_interval[i].first < 0) {
                subpath_read_interval[i].first = 0;
            }
            
            const Subpath& subpath = multipath_aln.subpath(i);
            int64_t subsequence_length = path_to_length(subpath.path());
            subpath_read_interval[i].second = subpath_read_interval[i].first + subsequence_length;
            
            if (!subpath.next_size()) {
                if (subpath_read_interval[i].second != multipath_aln.sequence().size()) {
#ifdef debug_multipath_mapper_mapping
                    cerr << "validation failure on using complete read" << endl;
                    cerr << "subpath " <<  i << " ends on sequence index " << subpath_read_interval[i].second << " of " << multipath_aln.sequence().size() << endl;
                    cerr << pb2json(subpath) << endl;
                    for (size_t j = 0; j < multipath_aln.subpath_size(); j++) {
                        cerr << j << " (" << subpath_read_interval[j].first << ", " << subpath_read_interval[j].second << "): ";
                        for (size_t k = 0; k < multipath_aln.subpath(j).next_size(); k++) {
                            cerr << multipath_aln.subpath(j).next(k) << " ";
                        }
                        cerr << endl;
                    }
#endif
                    return false;
                }
            }
            else {
                for (size_t j = 0; j < subpath.next_size(); j++) {
                    if (subpath_read_interval[subpath.next(j)].first >= 0) {
                        if (subpath_read_interval[subpath.next(j)].first != subpath_read_interval[i].second) {
#ifdef debug_multipath_mapper_mapping
                            cerr << "validation failure on read contiguity" << endl;
#endif
                            return false;
                        }
                    }
                    else {
                        subpath_read_interval[subpath.next(j)].first = subpath_read_interval[i].second;
                    }
                }
            }
        }
        
        // are all of the subpaths nonempty?
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            if (multipath_aln.subpath(i).path().mapping_size() == 0) {
#ifdef debug_multipath_mapper_mapping
                cerr << "validation failure on containing only nonempty paths" << endl;
                cerr << "subpath " << i << ": " << pb2json(multipath_aln.subpath(i)) << endl;
#endif
                return false;
            }
            for (size_t j = 0; j < multipath_aln.subpath(i).path().mapping_size(); j++) {
                if (multipath_aln.subpath(i).path().mapping(j).edit_size() == 0) {
#ifdef debug_multipath_mapper_mapping
                    cerr << "validation failure on containing only nonempty mappings" << endl;
                    cerr << "subpath " << i << ": " << pb2json(multipath_aln.subpath(i)) << endl;
#endif
                    return false;
                }
            }
        }
        
        
        // are the subpaths contiguous within the graph?
        
        auto validate_adjacent_mappings = [&](const Mapping& mapping_from, const Mapping& mapping_to) {
            size_t mapping_from_end_offset = mapping_from.position().offset() + mapping_from_length(mapping_from);
            if (mapping_from.position().node_id() == mapping_to.position().node_id() &&
                mapping_from.position().is_reverse() == mapping_to.position().is_reverse()) {
                if (mapping_to.position().offset() != mapping_from_end_offset) {
#ifdef debug_multipath_mapper_mapping
                    cerr << "validation failure on within-node adjacency" << endl;
                    cerr << pb2json(mapping_from) << "->" << pb2json(mapping_to) << endl;
#endif
                    return false;
                }
            }
            else {
                if (mapping_from_end_offset != xindex->node_length(mapping_from.position().node_id())) {
#ifdef debug_multipath_mapper_mapping
                    cerr << "validation failure on using edge at middle of node" << endl;
                    cerr << pb2json(mapping_from) << "->" << pb2json(mapping_to) << endl;
#endif
                    return false;
                }
                
                vector<Edge> edges = xindex->edges_of(mapping_from.position().node_id());
                bool found_edge = false;
                for (Edge& edge : edges) {
                    if (edge.from() == mapping_from.position().node_id() &&
                        edge.from_start() == mapping_from.position().is_reverse()) {
                        if (edge.to() == mapping_to.position().node_id() &&
                            edge.to_end() == mapping_to.position().is_reverse()) {
                            found_edge = true;
                            break;
                        }
                    }
                    if (edge.to() == mapping_from.position().node_id() &&
                        edge.to_end() != mapping_from.position().is_reverse()) {
                        if (edge.from() == mapping_to.position().node_id() &&
                            edge.from_start() != mapping_to.position().is_reverse()) {
                            found_edge = true;
                            break;
                        }
                    }
                }
                
                if (!found_edge) {
#ifdef debug_multipath_mapper_mapping
                    cerr << "validation failure on nodes not connected by an edge" << endl;
                    cerr << pb2json(mapping_from) << "->" << pb2json(mapping_to) << endl;
#endif
                    return false;
                }
            }
            return true;
        };
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            const Path& path = subpath.path();
            for (size_t j = 1; j < path.mapping_size(); j++) {
                if (!validate_adjacent_mappings(path.mapping(j - 1), path.mapping(j))) {
                    return false;
                }
            }
            const Mapping& final_mapping = path.mapping(path.mapping_size() - 1);
            for (size_t j = 0; j < subpath.next_size(); j++) {
                if (!validate_adjacent_mappings(final_mapping, multipath_aln.subpath(subpath.next(j)).path().mapping(0))) {
                    return false;
                }
            }
        }
        
        
        // do the paths represent valid alignments of the associated read string and graph path?
        
        auto validate_mapping_edits = [&](const Mapping& mapping, const string& subseq) {
            string node_seq = xindex->node_sequence(mapping.position().node_id());
            string rev_node_seq = reverse_complement(node_seq);
            size_t node_idx = mapping.position().offset();
            size_t seq_idx = 0;
            for (size_t i = 0; i < mapping.edit_size(); i++) {
                const Edit& edit = mapping.edit(i);
                if (edit_is_match(edit)) {
                    for (size_t j = 0; j < edit.from_length(); j++, node_idx++, seq_idx++) {
                        if ((mapping.position().is_reverse() ? rev_node_seq[node_idx] : node_seq[node_idx]) != subseq[seq_idx]) {
#ifdef debug_multipath_mapper_mapping
                            cerr << "validation failure on match that does not match" << endl;
                            cerr << pb2json(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                    }
                }
                else if (edit_is_sub(edit)) {
                    for (size_t j = 0; j < edit.from_length(); j++, node_idx++, seq_idx++) {
                        if ((mapping.position().is_reverse() ? rev_node_seq[node_idx] : node_seq[node_idx]) == subseq[seq_idx]) {
#ifdef debug_multipath_mapper_mapping
                            cerr << "validation failure on mismatch that matches" << endl;
                            cerr << pb2json(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                        if (edit.sequence()[j] != subseq[seq_idx]) {
#ifdef debug_multipath_mapper_mapping
                            cerr << "validation failure on substitution sequence that does not match read" << endl;
                            cerr << pb2json(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                    }
                }
                else if (edit_is_insertion(edit)) {
                    for (size_t j = 0; j < edit.to_length(); j++, seq_idx++) {
                        if (edit.sequence()[j] != subseq[seq_idx]) {
#ifdef debug_multipath_mapper_mapping
                            cerr << "validation failure on insertion sequence that does not match read" << endl;
                            cerr << pb2json(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                    }
                }
                else if (edit_is_deletion(edit)) {
                    node_idx += edit.from_length();
                }
            }
            return true;
        };
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            const Path& path = subpath.path();
            size_t read_start = subpath_read_interval[i].first;
            for (size_t j = 0; j < path.mapping_size(); j++) {
                size_t read_mapping_len = mapping_to_length(path.mapping(j));
                if (!validate_mapping_edits(path.mapping(j), multipath_aln.sequence().substr(read_start, read_mapping_len))) {
                    return false;
                }
                read_start += read_mapping_len;
            }
        }
        
        // do the scores match the alignments?
        
        // TODO: this really deserves a test, but there's a factoring problem because the qual adj aligner needs to know
        // the node sequence to score mismatches but the node sequence is not stored in the Alignment object
        
//        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
//            const Subpath& subpath = multipath_aln.subpath(i);
//            Alignment& alignment;
//            *alignment.mutable_sequence() = multipath_aln.sequence().substr(subpath_read_interval[i].first,
//                                                                            subpath_read_interval[i].second - subpath_read_interval[i].first);
//            *alignment.mutable_quality() = multipath_aln.quality().substr(subpath_read_interval[i].first,
//                                                                          subpath_read_interval[i].second - subpath_read_interval[i].first);
//            *alignment.mutable_path() = subpath.path();
//        }
        
        
        return true;
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
    
    MultipathAlignmentGraph::MultipathAlignmentGraph(VG& vg, const MultipathMapper::memcluster_t& hits,
                                                     const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                                     gcsa::GCSA* gcsa, SnarlManager* cutting_snarls, int64_t max_snarl_cut_size) {

        // create the injection translator, which maps a node in the original graph to every one of its occurrences
        // in the dagified graph
        unordered_multimap<id_t, pair<id_t, bool> > injection_trans;
        for (const auto& trans_record : projection_trans) {
#ifdef debug_multipath_mapper_alignment
            cerr << trans_record.second.first << "->" << trans_record.first << (trans_record.second.second ? "-" : "+") << endl;
#endif
            injection_trans.emplace(trans_record.second.first, make_pair(trans_record.first, trans_record.second.second));
        }
        
        // initialize the match nodes
        create_match_nodes(vg, hits, projection_trans, injection_trans);
        
        if (gcsa) {
            // we indicated that these MEMs came from a GCSA, so there might be order-length MEMs that we can combine
            collapse_order_length_runs(vg, gcsa);
        }
        
#ifdef debug_multipath_mapper_alignment
        cerr << "nodes after adding and collapsing:" << endl;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            ExactMatchNode& match_node = match_nodes[i];
            cerr << i << " " << pb2json(match_node.path) << " ";
            for (auto iter = match_node.begin; iter != match_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
        }
#endif
        
        if (cutting_snarls) {
            // we indicated a snarl manager that owns the snarls we want to cut out of exact matches
            resect_snarls_from_paths(cutting_snarls, projection_trans, max_snarl_cut_size);
        }
        
        // compute reachability and add edges
        add_reachability_edges(vg, hits, projection_trans, injection_trans);
        
#ifdef debug_multipath_mapper_alignment
        cerr << "final mem graph:" << endl;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            ExactMatchNode& match_node = match_nodes[i];
            cerr << i << " " << pb2json(match_node.path) << " ";
            for (auto iter = match_node.begin; iter != match_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            cerr << "\t";
            for (auto edge : match_node.edges) {
                cerr << "(to:" << edge.first << ", graph dist:" << edge.second << ", read dist: " << (match_nodes[edge.first].begin - match_node.end) << ") ";
            }
            cerr << endl;
        }
#endif
    }
            
    void MultipathAlignmentGraph::create_match_nodes(VG& vg, const MultipathMapper::memcluster_t& hits,
                                                     const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                                     const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans) {
        
        
#ifdef debug_multipath_mapper_alignment
        cerr << "walking out MEMs in graph" << endl;
#endif
        
        // map of node ids in the dagified graph to the indices in the matches that contain them
        unordered_map<int64_t, vector<int64_t>> node_matches;
        
        // walk the matches and filter out redundant sub-MEMs
        for (int64_t i = 0; i < hits.size(); i++) {
            
            const pair<const MaximalExactMatch*, pos_t>& hit = hits[i];
            
            // the part of the read we're going to match
            string::const_iterator begin = hit.first->begin;
            string::const_iterator end = hit.first->end;
            int64_t mem_length = end - begin;
            // the start of the hit in the original graph
            const pos_t& hit_pos = hit.second;
            
#ifdef debug_multipath_mapper_alignment
            cerr << "walking MEM hit " << hit_pos << " " << hit.first->sequence() << endl;
#endif
            
            auto hit_range = injection_trans.equal_range(id(hit_pos));
            for (auto iter = hit_range.first; iter != hit_range.second; iter++) {
                // this graph is unrolled/dagified, so all orientations should match
                if ((*iter).second.second != is_rev(hit_pos)) {
                    continue;
                }
                
                // an id that corresponds to the original node
                id_t injected_id = (*iter).second.first;
                
#ifdef debug_multipath_mapper_alignment
                cerr << "hit node exists in graph as " << injected_id << endl;
#endif
                
                // check all MEMs that traversed this node to see if this is a redundant sub-MEM
                bool is_partial_mem = false;
                if (node_matches.count(injected_id)) {
#ifdef debug_multipath_mapper_alignment
                    cerr << "we need to check if this is a redundant sub MEM, there are previous that visited this hit" << endl;
#endif
                    
                    for (int64_t j : node_matches[injected_id]) {
                        ExactMatchNode& match_node = match_nodes[j];
                        
                        if (begin < match_node.begin || end > match_node.end) {
#ifdef debug_multipath_mapper_alignment
                            if (begin < match_node.begin) {
                                cerr << "this MEM is earlier in the read than the other, so this is not redundant" << endl;
                            }
                            else if (end > match_node.end) {
                                cerr << "this MEM is later in the read than the other, so this is not redundant" << endl;
                            }
#endif
                            // the hit does not fall on the same section of the read as the other match, so
                            // it cannot be contained in it
                            continue;
                        }
                        
                        int64_t relative_offset = begin - match_node.begin;
#ifdef debug_multipath_mapper_alignment
                        cerr << "the match on node " << j << " has an relative offset of " << relative_offset << " to the this MEM in the read" << endl;
#endif
                        
                        Path& path = match_node.path;
                        
                        // if this is a partial MEM, we should be able to predict its hit location by traversing the path
                        // of the parent MEM by a distance equal to the relative offset
                        
#ifdef debug_multipath_mapper_alignment
                        cerr << "traversing putative parent MEM with path " << pb2json(path) << endl;
#endif
                        
                        int64_t prefix_length = 0;
                        for (size_t k = 0; k < path.mapping_size(); k++) {
                            if (prefix_length > relative_offset) {
#ifdef debug_multipath_mapper_alignment
                                cerr << "we have passed where the location would be, breaking out of loop" << endl;
#endif
                                break;
                            }
                            const Mapping& mapping = path.mapping(k);
                            // the length through this mapping
                            int64_t prefix_through_length = prefix_length + mapping_from_length(mapping);
#ifdef debug_multipath_mapper_alignment
                            cerr << "after traversing the " << k << "-th step, we have covered a distance of " << prefix_through_length << endl;
#endif
                            if (prefix_through_length > relative_offset) {
                                // we cross the relative offset on this node, so check if the path is in the predicted
                                // position for a redundant sub-MEM
                                id_t node_id_here = mapping.position().node_id();
                                is_partial_mem = is_partial_mem || (injected_id == node_id_here
                                                                    && offset(hit_pos) == mapping.position().offset() + relative_offset - prefix_length
                                                                    && projection_trans.at(node_id_here).second == is_rev(hit_pos));
#ifdef debug_multipath_mapper_alignment
                                cerr << "this mapping crosses where we would expect a child to be: " << node_id_here << (projection_trans.at(node_id_here).second ? "-" : "+") << ":" << mapping.position().offset() + relative_offset - prefix_length << endl;
                                cerr << "this MEM is actually at: " << injected_id << (is_rev(hit_pos) ? "-" : "+") << ":" << offset(hit_pos) << endl;
#endif
                                
                            }
                            prefix_length = prefix_through_length;
                        }
                        if (is_partial_mem) {
                            break;
                        }
                    }
                }
                
                // don't walk the match of false partial hits
                if (is_partial_mem) {
#ifdef debug_multipath_mapper_alignment
                    cerr << "this MEM is identified as a redundant sub-MEM, so we skip it" << endl;
#endif
                    continue;
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << "performing DFS to walk out match" << endl;
#endif
                
                // stack for DFS, each record contains tuples of (read begin, node offset, next node index, next node ids)
                vector<tuple<string::const_iterator, size_t, size_t, vector<NodeTraversal>>> stack;
                stack.emplace_back(begin, offset(hit_pos), 0,
                                   vector<NodeTraversal>{NodeTraversal(vg.get_node(injected_id))});
                
                while (!stack.empty()) {
                    auto& back = stack.back();
                    if (get<2>(back) == get<3>(back).size()) {
#ifdef debug_multipath_mapper_alignment
                        cerr << "traversed all edges out of current traversal" << endl;
#endif
                        stack.pop_back();
                        continue;
                    }
                    NodeTraversal trav = get<3>(back)[get<2>(back)];
                    get<2>(back)++;
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "checking node " << trav.node->id() << endl;
#endif
                    
                    const string& node_seq = trav.node->sequence();
                    size_t node_idx = get<1>(back);
                    string::const_iterator read_iter = get<0>(back);
                    
                    // look for a match along the entire node sequence
                    for (; node_idx < node_seq.size() && read_iter != end; node_idx++, read_iter++) {
                        if (node_seq[node_idx] != *read_iter) {
#ifdef debug_multipath_mapper_alignment
                            cerr << "node sequence does not match read" << endl;
#endif
                            break;
                        }
                    }
                    
                    if (read_iter == end) {
                        // finished walking match
#ifdef debug_multipath_mapper_alignment
                        cerr << "reached end of read sequence, finished walking match" << endl;
#endif
                        break;
                    }
                    else if (node_idx == node_seq.size()) {
                        // matched entire node
                        stack.emplace_back(read_iter, 0, 0, vector<NodeTraversal>());
                        vg.nodes_next(trav, get<3>(stack.back()));
                    }
                }
                
                // if we left a trace in the stack we found a complete match, but sometimes MEMs that overhang
                // the edge of the subraph find their way in (only when they are not part of the alignment
                // represented by the cluster used to query the subgraph) in which case we just skip this MEM
                if (stack.empty()) {
#ifdef debug_multipath_mapper_alignment
                    cerr << "this MEM overhangs the end of the graph" << endl;
#endif
                    continue;
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << "converting into a Path at idx " << match_nodes.size() << endl;
#endif
                
                int64_t match_node_idx = match_nodes.size();
                match_nodes.emplace_back();
                ExactMatchNode& match_node = match_nodes.back();
                Path& path = match_node.path;
                match_node.begin = begin;
                match_node.end = end;
                int64_t length_remaining = end - begin;
                
                // walk out the match
                int32_t rank = 1;
                for (auto search_record : stack) {
                    int64_t offset = get<1>(search_record);
                    Node* node = get<3>(search_record)[get<2>(search_record) - 1].node;
                    int64_t length = std::min((int64_t) node->sequence().size() - offset, length_remaining);
                    
                    Mapping* mapping = path.add_mapping();
                    mapping->set_rank(rank);
                    
                    Edit* edit = mapping->add_edit();
                    edit->set_from_length(length);
                    edit->set_to_length(length);
                    
                    // note: the graph is dagified and unrolled, so all hits should be on the forward strand
                    Position* position = mapping->mutable_position();
                    position->set_node_id(node->id());
                    position->set_offset(offset);
                    
                    // record that each node occurs in this match so we can filter out sub-MEMs
                    node_matches[node->id()].push_back(match_node_idx);
#ifdef debug_multipath_mapper_alignment
                    cerr << "associating node " << node->id() << " with a match at idx " << match_node_idx << endl;
#endif
                    
                    rank++;
                    length_remaining -= length;
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << pb2json(path) << endl;
#endif
            }
        }
    }
            
    void MultipathAlignmentGraph::collapse_order_length_runs(VG& vg, gcsa::GCSA* gcsa) {
        
#ifdef debug_multipath_mapper_alignment
        cerr << "looking for runs of order length MEMs to collapse" << endl;
#endif
        
        vector<vector<size_t>> merge_groups;
        
        size_t num_order_length_mems = 0;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            
            ExactMatchNode& match_node = match_nodes[i];
            
            if (match_node.end - match_node.begin < gcsa->order()) {
                // we have passed all of the order length MEMs, bail out of loop
#ifdef debug_multipath_mapper_alignment
                cerr << "found " << i << " order length MEMs" << endl;
#endif
                
                num_order_length_mems = i;
                break;
            }
        }
        
        // find the order of the order-length MEMs lexicographically along the read
        vector<size_t> order(num_order_length_mems, 0);
        for (size_t i = 1; i < order.size(); i++) {
            order[i] = i;
        }
        sort(order.begin(), order.end(), [&](size_t i, size_t j) {
            return match_nodes[i].begin < match_nodes[j].begin || (match_nodes[i].begin == match_nodes[j].begin &&
                                                                   match_nodes[i].end < match_nodes[j].end);
        });
        
        for (size_t i : order) {
            
            ExactMatchNode& match_node = match_nodes[i];
            
#ifdef debug_multipath_mapper_alignment
            cerr << "## checking if MEM " << i << " can be an extension: " << endl;
            cerr << "\t";
            for (auto iter = match_node.begin; iter != match_node.end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
            cerr << "\t" << pb2json(match_node.path) << endl;
#endif
            
            // try to find any run of MEMs that could be merged with this MEM
            bool found_merge_group = false;
            for (size_t j = 0; j < merge_groups.size(); j++) {
                
                vector<size_t>& merge_group = merge_groups[j];
                
                // because of the sort order, the last node in this run should overlap the current MEM
                // if any of them can
                ExactMatchNode& last_run_node = match_nodes[merge_group[merge_group.size() - 1]];
                
#ifdef debug_multipath_mapper_alignment
                cerr << "checking against extending MEM " << merge_group[merge_group.size() - 1] << " in merge group " << j << ":" << endl;
                cerr << "\t";
                for (auto iter = last_run_node.begin; iter != last_run_node.end; iter++) {
                    cerr << *iter;
                }
                cerr << endl;
                cerr << "\t" << pb2json(last_run_node.path) << endl;
#endif
                
                // do they overhang an amount on the read that indicates they overlap and could be merged?
                int64_t overhang = last_run_node.end - match_node.begin;
                if (last_run_node.begin < match_node.begin && match_node.end > last_run_node.end && overhang >= 0) {
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "MEMs overlap on the read, checking for consistency with overhang on the path" << endl;
#endif
                    
                    // get the initial position of the node further to the right
                    pos_t match_node_initial_pos = make_pos_t(match_node.path.mapping(0).position());
                    
                    // get the position at the overhang back from the end of the node further to the left
                    int64_t remaining = last_run_node.end - last_run_node.begin;
                    pos_t last_run_node_internal_pos;
                    for (size_t k = 0; k < last_run_node.path.mapping_size(); k++) {
                        
                        int64_t mapping_length = mapping_from_length(last_run_node.path.mapping(k));
                        
                        if (remaining - mapping_length < overhang) {
                            // we will cross the position that should line up with the initial position on this mapping
                            
                            const Position& overhang_position = last_run_node.path.mapping(k).position();
                            
                            get_id(last_run_node_internal_pos) = overhang_position.node_id();
                            get_is_rev(last_run_node_internal_pos) = overhang_position.is_reverse();
                            get_offset(last_run_node_internal_pos) = overhang_position.offset() + remaining - overhang;
                            
                            break;
                        }
                        
                        remaining -= mapping_length;
                    }
                    
                    // get the final position of the node further to the left
                    const Mapping& final_mapping = last_run_node.path.mapping(last_run_node.path.mapping_size() - 1);
                    const Position& final_mapping_position = final_mapping.position();
                    pos_t last_run_node_final_pos = make_pos_t(final_mapping_position.node_id(),
                                                               final_mapping_position.is_reverse(),
                                                               final_mapping_position.offset() + mapping_from_length(final_mapping));
                    
                    // get the position at the overhang into the node further to the right
                    remaining = match_node.end - match_node.begin;
                    pos_t match_node_internal_pos;
                    for (int64_t k = match_node.path.mapping_size() - 1; k >= 0; k--) {
                        
                        int64_t mapping_length = mapping_from_length(match_node.path.mapping(k));
                        remaining -= mapping_length;
                        
                        if (remaining < overhang) {
                            // we will cross the position that should line up with the initial position on this mapping
                            
                            const Position& overhang_position = match_node.path.mapping(k).position();
                            
                            get_id(match_node_internal_pos) = overhang_position.node_id();
                            get_is_rev(match_node_internal_pos) = overhang_position.is_reverse();
                            get_offset(match_node_internal_pos) = overhang_position.offset() + overhang - remaining;
                            
                            break;
                        }
                    }
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "overhang segments on extension node " << match_node_initial_pos << ", " << match_node_internal_pos << " and extending node " << last_run_node_internal_pos << ", " << last_run_node_final_pos << endl;
#endif
                    
                    // do the positions match up as we would expect if these are actually part of the same match?
                    if (match_node_initial_pos == last_run_node_internal_pos && last_run_node_final_pos == match_node_internal_pos) {
                        
#ifdef debug_multipath_mapper_alignment
                        cerr << "found matching overhang" << endl;
#endif
                        
                        // add the match node to this merge group
                        merge_group.push_back(i);
                        found_merge_group = true;
                        
                        break;
                    }
                    else if (overhang == 0 && offset(match_node_initial_pos) == 0) {
                        // it could still be that these are two end-to-end matches that got assigned to the beginning
                        // and end of two nodes connected by an edge
                        
                        if (offset(last_run_node_final_pos) == vg.get_node(final_mapping_position.node_id())->sequence().size()
                            && vg.has_edge(NodeSide(id(last_run_node_final_pos), !is_rev(last_run_node_final_pos)),
                                           NodeSide(id(match_node_initial_pos), is_rev(match_node_initial_pos)))) {
                            
#ifdef debug_multipath_mapper_alignment
                            cerr << "found end to end connection over an edge" << endl;
#endif
                            
                            // add the match node to this merge group
                            merge_group.push_back(i);
                            found_merge_group = true;
                            
                            break;
                        }
                    }
                }
            }
            
            if (!found_merge_group) {
                // make a new merge group consisting of only this MEM
                merge_groups.emplace_back(1, i);
            }
        }
        
#ifdef debug_multipath_mapper_alignment
        cerr << "merge groups among order length MEMs:" << endl;
        for (auto& group : merge_groups) {
            cerr << "\t";
            for (auto i : group) {
                cerr << i << " ";
            }
            cerr << endl;
        }
#endif
        
        if (merge_groups.size() != num_order_length_mems) {
            // we found at least one merge to do, now we need to actually do the merges
            
            unordered_set<size_t> to_remove;
            
            for (const vector<size_t>& merge_group : merge_groups) {
                
                // merge the paths into the first node in the group (arbitrarily)
                ExactMatchNode& merge_into_node = match_nodes[merge_group[0]];
                
                for (size_t i = 1; i < merge_group.size(); i++) {
                    
                    // mark the node we're merging from for removal
                    to_remove.insert(merge_group[i]);
                    
                    ExactMatchNode& merge_from_node = match_nodes[merge_group[i]];
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "merging into node " << merge_group[0] << " path " << pb2json(merge_into_node.path) << endl;
                    cerr << "from node " << merge_group[i] << " path " << pb2json(merge_from_node.path) << endl;
#endif
                    
                    // walk backwards until we find the first mapping to add
                    int64_t to_add_length = merge_from_node.end - merge_into_node.end;
                    int64_t remaining = to_add_length;
                    int64_t first_mapping_to_add_idx = 0;
                    for (int64_t j = merge_from_node.path.mapping_size() - 1; j >= 0; j--) {
                        remaining -= mapping_from_length(merge_from_node.path.mapping(j));
                        if (remaining <= 0) {
                            first_mapping_to_add_idx = j;
                            break;
                        }
                    }
                    
                    //cerr << "first mapping to add is idx " << first_mapping_to_add_idx << ": " << pb2json(merge_from_node.path.mapping(first_mapping_to_add_idx)) << endl;
                    //cerr << "adding total length " << to_add_length << " remaining length is " << remaining << endl;
                    
                    // handle the first mapping we add as a special case
                    const Mapping& first_mapping_to_add = merge_from_node.path.mapping(first_mapping_to_add_idx);
                    Mapping* final_merging_mapping = merge_into_node.path.mutable_mapping(merge_into_node.path.mapping_size() - 1);
                    //cerr << "merging it into mapping " << pb2json(*final_merging_mapping) << endl;
                    if (final_merging_mapping->position().node_id() == first_mapping_to_add.position().node_id() &&
                        final_merging_mapping->position().is_reverse() == first_mapping_to_add.position().is_reverse() &&
                        final_merging_mapping->position().offset() == first_mapping_to_add.position().offset() &&
                        mapping_from_length(*final_merging_mapping) == -remaining) {
                        
                        // the mappings are on the same node, so they can be combined
                        int64_t mapping_to_add_length = mapping_from_length(first_mapping_to_add) + remaining;
                        Edit* final_edit = final_merging_mapping->mutable_edit(final_merging_mapping->edit_size() - 1);
                        final_edit->set_from_length(final_edit->from_length() + mapping_to_add_length);
                        final_edit->set_to_length(final_edit->to_length() + mapping_to_add_length);
                        
                        //cerr << "blended mapping is " << pb2json(*final_merging_mapping) << endl;
                    }
                    else {
                        // we need to add this as a new mapping
                        Mapping* new_mapping = merge_into_node.path.add_mapping();
                        *new_mapping = first_mapping_to_add;
                        new_mapping->set_rank(final_merging_mapping->rank() + 1);
                        //cerr << "new separate mapping is " << pb2json(*new_mapping) << endl;
                        
                    }
                    
                    // add the remaining mappings as new mappings
                    for (size_t j = first_mapping_to_add_idx + 1; j < merge_from_node.path.mapping_size(); j++) {
                        Mapping* new_mapping = merge_into_node.path.add_mapping();
                        *new_mapping = merge_from_node.path.mapping(j);
                        new_mapping->set_rank(merge_into_node.path.mapping(merge_into_node.path.mapping_size() - 2).rank() + 1);
                        //cerr << "added new separate mapping " << pb2json(*new_mapping) << endl;
                    }
                    
                    // merge the substrings on the node
                    merge_into_node.end = merge_from_node.end;
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "merged path is " << pb2json(merge_from_node.path) << endl;
#endif
                }
            }
            
            // remove all of the nodes we merged into other nodes
            size_t removed_so_far = 0;
            for (size_t i = 0; i < match_nodes.size(); i++) {
                if (to_remove.count(i)) {
                    removed_so_far++;
                }
                else if (removed_so_far > 0) {
                    match_nodes[i - removed_so_far] = move(match_nodes[i]);
                }
            }
            
            match_nodes.resize(match_nodes.size() - to_remove.size());
        }
    }
            
    void MultipathAlignmentGraph::resect_snarls_from_paths(SnarlManager* cutting_snarls,
                                                           const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                                           int64_t max_snarl_cut_size) {
#ifdef debug_multipath_mapper_alignment
        cerr << "cutting with snarls" << endl;
#endif
        
        size_t num_original_match_nodes = match_nodes.size();
        for (size_t i = 0; i < num_original_match_nodes; i++) {
            
            // first compute the segments we want to cut out
            
            ExactMatchNode* match_node = &match_nodes[i];
            Path* path = &match_node->path;
            
#ifdef debug_multipath_mapper_alignment
            cerr << "cutting node at index " << i << " with path " << pb2json(*path) << endl;
#endif
            
            // this list holds the beginning of the current segment at each depth in the snarl hierarchy
            // as we traverse the exact match, the beginning is recorded in both sequence distance and node index
            list<pair<size_t, size_t>> level_segment_begin;
            level_segment_begin.emplace_back(0, 0);
            
            // we record which segments we are going to cut out of the match here
            vector<pair<size_t, size_t>> cut_segments;
            
            auto curr_level = level_segment_begin.begin();
            size_t prefix_length = 0;
            for (size_t j = 0, last = path->mapping_size() - 1; j <= last; j++) {
                const Position& position = path->mapping(j).position();
                const auto& projection = projection_trans.at(position.node_id());
                id_t projected_id = projection.first;
                bool projected_rev = (projection.second != position.is_reverse());
                
                if (j > 0) {
                    // we have entered this node on this iteration
                    if (cutting_snarls->into_which_snarl(projected_id, !projected_rev)) {
                        // as we enter this node, we are leaving the snarl we were in
                        
                        // since we're going up a level, we need to check whether we need to cut out the segment we've traversed
                        if (prefix_length - curr_level->first <= max_snarl_cut_size || !max_snarl_cut_size) {
                            cut_segments.emplace_back(curr_level->second, j);
                        }
                        
                        curr_level++;
                        if (curr_level == level_segment_begin.end()) {
                            // we were already at the highest level seen so far, so we need to add a new one
                            // the entire previous part of the match is contained in this level, so we start
                            // the segment from 0
                            curr_level = level_segment_begin.insert(level_segment_begin.end(), make_pair(0, 0));
                        }
                    }
                }
                
                // cross to the other side of the node
                prefix_length += mapping_from_length(path->mapping(j));
                
                if (j < last) {
                    // we are going to leave this node next iteration
                    if (cutting_snarls->into_which_snarl(projected_id, projected_rev)) {
                        // as we leave this node, we are entering a new deeper snarl
                        
                        // the segment in the new level will begin at the end of the current node
                        if (curr_level == level_segment_begin.begin()) {
                            // we are already at the lowest level seen so far, so we need to add a new one
                            level_segment_begin.emplace_front(prefix_length, j + 1);
                            curr_level--;
                        }
                        else {
                            // the lower level is in the record already, so we update its segment start
                            curr_level--;
                            *curr_level = make_pair(prefix_length, j + 1);
                        }
                    }
                }
            }
            
            // check the final segment for a cut unless we're at the highest level in the match
            auto last = level_segment_begin.end();
            last--;
            if ((prefix_length - curr_level->first <= max_snarl_cut_size || !max_snarl_cut_size) && curr_level != last) {
                cut_segments.emplace_back(curr_level->second, path->mapping_size());
            }
            
#ifdef debug_multipath_mapper_alignment
                cerr << "found " << cut_segments.size() << " cut segments:" << endl;
                for (auto seg : cut_segments) {
                    cerr << "\t" << seg.first << ":" << seg.second << endl;
                }
#endif
            
            // did we cut out any segments?
            if (!cut_segments.empty()) {
                
                // we may have decided to cut the segments of both a parent and child snarl, so now we
                // collapse the list of intervals, which is sorted on the end index by construction
                //
                // snarl nesting properties guarantee that there will be at least one node between any
                // cut segments that are not nested, so we don't need to deal with the case where the
                // segments are partially overlapping (i.e. it's a bit easier than the general interval
                // intersection problem)
                vector<pair<size_t, size_t>> keep_segments;
                size_t curr_keep_seg_end = path->mapping_size();
                auto riter = cut_segments.rbegin();
                if (riter->second == curr_keep_seg_end) {
                    // don't add an empty keep segment in the first position
                    curr_keep_seg_end = riter->first;
                    riter++;
                }
                for (; riter != cut_segments.rend(); riter++) {
                    if (riter->second < curr_keep_seg_end) {
                        // this is a new interval
                        keep_segments.emplace_back(riter->second, curr_keep_seg_end);
                        curr_keep_seg_end = riter->first;
                    }
                }
                if (curr_keep_seg_end > 0) {
                    // we are not cutting off the left tail, so add a keep segment for it
                    keep_segments.emplace_back(0, curr_keep_seg_end);
                }
                
                // make a new node for all but one of the keep segments
                size_t prefix_length = 0;
                size_t prefix_idx = 0;
                for (auto iter = keep_segments.end() - 1; iter != keep_segments.begin(); iter--) {
#ifdef debug_multipath_mapper_alignment
                    cerr << "making path for keep segment " << iter->first << ":" << iter->second << " at idx " << match_nodes.size() << endl;
#endif
                    match_nodes.emplace_back();
                    
                    // update pointers in case the vector reallocates
                    match_node = &match_nodes[i];
                    path = &match_node->path;
                    
                    // measure the length of path between this keep segment and the last one
                    while (prefix_idx < iter->first) {
                        prefix_length += mapping_from_length(path->mapping(prefix_idx));
                        prefix_idx++;
                    }
                    
                    size_t keep_segment_length = 0;
                    ExactMatchNode& cut_node = match_nodes.back();
                    Path& cut_path = cut_node.path;
                    // transfer over the keep segment from the main path and measure the length
                    int32_t rank = 1;
                    for (size_t j = iter->first; j < iter->second; j++, rank++) {
                        Mapping* mapping = cut_path.add_mapping();
                        *mapping = path->mapping(j);
                        mapping->set_rank(rank);
                        keep_segment_length += mapping_from_length(*mapping);
                    }
                    
                    // identify the substring of the MEM that stays on this node
                    cut_node.begin = match_node->begin + prefix_length;
                    cut_node.end = cut_node.begin + keep_segment_length;
                    
                    prefix_length += keep_segment_length;
                    prefix_idx = iter->second;
#ifdef debug_multipath_mapper_alignment
                    cerr << "new cut path: " << pb2json(cut_path) << endl;
#endif
                    
                }
                
                while (prefix_idx < keep_segments.front().first) {
                    prefix_length += mapping_from_length(path->mapping(prefix_idx));
                    prefix_idx++;
                }
                
                // replace the path of the original node with the final keep segment
                size_t keep_segment_length = 0;
                Path new_path;
                int32_t rank = 1;
                for (size_t j = keep_segments.front().first; j < keep_segments.front().second; j++, rank++) {
                    Mapping* mapping = new_path.add_mapping();
                    *mapping = path->mapping(j);
                    mapping->set_rank(rank);
                    keep_segment_length += mapping_from_length(*mapping);
                }
                *path = new_path;
                
                // update the substring of MEM to match the path
                match_node->begin += prefix_length;
                match_node->end = match_node->begin + keep_segment_length;
#ifdef debug_multipath_mapper_alignment
                cerr << "new in place cut path for keep segment " << keep_segments.front().first << ":" << keep_segments.front().second << " " << pb2json(new_path) << endl;
#endif
            }
        }
    }
            
    void MultipathAlignmentGraph::add_reachability_edges(VG& vg, const MultipathMapper::memcluster_t& hits,
                                                         const unordered_map<id_t, pair<id_t, bool>>& projection_trans,
                                                         const unordered_multimap<id_t, pair<id_t, bool>>& injection_trans) {
        
        
#ifdef debug_multipath_mapper_alignment
        cerr << "computing reachability" << endl;
#endif
        
        // now we calculate reachability between the walked matches so we know which ones
        // to connect with intervening alignments
        
        auto start_offset = [&](size_t idx) {
            return match_nodes[idx].path.mapping(0).position().offset();
        };
        
        auto end_offset = [&](size_t idx) {
            Path& path = match_nodes[idx].path;
            const Mapping& mapping = path.mapping(path.mapping_size() - 1);
            return mapping.position().offset() + mapping_from_length(mapping);
        };
        
        auto start_node_id = [&](size_t idx) {
            return match_nodes[idx].path.mapping(0).position().node_id();
        };
        
        auto end_node_id = [&](size_t idx) {
            Path& path = match_nodes[idx].path;
            return path.mapping(path.mapping_size() - 1).position().node_id();
        };
        
        auto endpoint_offset = [&](size_t idx, bool end) {
            return end ? end_offset(idx) : start_offset(idx);
        };
        
        auto endpoint_node_id = [&](size_t idx, bool end) {
            return end ? end_node_id(idx) : start_node_id(idx);
        };
        
        // record the start and end node ids of every exact match
        unordered_map<id_t, vector<size_t>> exact_match_starts;
        unordered_map<id_t, vector<size_t>> exact_match_ends;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            Path& path = match_nodes[i].path;
            exact_match_starts[path.mapping(0).position().node_id()].push_back(i);
            exact_match_ends[path.mapping(path.mapping_size() - 1).position().node_id()].push_back(i);
        }
        
#ifdef debug_multipath_mapper_alignment
        cerr << "recorded starts: " << endl;
        for (const auto& rec : exact_match_starts) {
            cerr << "\t" << rec.first << ": ";
            for (auto l : rec.second) {
                cerr << l << " ";
            }
            cerr << endl;
        }
        
        cerr << "recorded ends: " << endl;
        for (const auto& rec : exact_match_ends) {
            cerr << "\t" << rec.first << ": ";
            for (auto l : rec.second) {
                cerr << l << " ";
            }
            cerr << endl;
        }
#endif
        
        
        // sort the MEMs starting and ending on each node in node sequence order
        for (pair<const id_t, vector<size_t>>& node_match_starts : exact_match_starts) {
            std::sort(node_match_starts.second.begin(), node_match_starts.second.end(),
                      [&](const size_t idx_1, const size_t idx_2) {
                          return start_offset(idx_1) < start_offset(idx_2);
                      });
        }
        for (pair<const id_t, vector<size_t>>& node_match_ends : exact_match_ends) {
            std::sort(node_match_ends.second.begin(), node_match_ends.second.end(),
                      [&](const size_t idx_1, const size_t idx_2) {
                          return end_offset(idx_1) < end_offset(idx_2);
                      });
        }
        
        // some structures we will fill out with DP:
        
        // for each node, the starts and ends of MEMs that can reach this node without passing any other
        // MEM starts or ends
        unordered_map<id_t, unordered_map<size_t, size_t>> reachable_ends;
        unordered_map<id_t, unordered_map<size_t, size_t>> reachable_starts;
        
        // for the start of each MEM, the starts and ends of other MEMs that can reach it without passing any
        // other start or end
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_starts_from_start;
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_ends_from_start;
        
        // for the end of each MEM, the ends of other MEMs that can reach it without passing any
        // other start or end
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_ends_from_end;
        unordered_map<size_t, vector<pair<size_t, size_t>>> reachable_starts_from_end;
        
        // note: graph has been sorted into topological order
        
        Graph& graph = vg.graph;
        for (int64_t i = 0; i < graph.node_size(); i++) {
            Node* node = graph.mutable_node(i);
            id_t node_id = node->id();
            
#ifdef debug_multipath_mapper_alignment
            cerr << "DP step for node " << node_id << endl;
#endif
            
            size_t node_length = node->sequence().size();
            
            // do any MEMs start or end on this node?
            bool contains_starts = exact_match_starts.count(node_id);
            bool contains_ends = exact_match_ends.count(node_id);
            
            // we will use DP to carry reachability information forward onto the next nodes
            vector<NodeTraversal> nexts;
            vg.nodes_next(NodeTraversal(node), nexts);
            
            if (contains_starts && contains_ends) {
                // since there are both starts and ends on this node, we have to traverse both lists simultaneously
                // to assess reachability within the same node
                
#ifdef debug_multipath_mapper_alignment
                cerr << "node " << node_id << " contains both starts and ends of MEMs" << endl;
#endif
                
                vector<size_t>& ends = exact_match_ends[node_id];
                vector<size_t>& starts = exact_match_starts[node_id];
                
                
                // find the range of starts and ends in the list with the same offset
                
                size_t start_range_begin = 0;
                size_t start_range_end = 1;
                size_t end_range_begin = 0;
                size_t end_range_end = 1;
                
                size_t curr_start_offset = start_offset(starts[start_range_begin]);
                size_t curr_end_offset = end_offset(ends[end_range_begin]);
                size_t prev_offset = 0;
                
                while (end_range_end == ends.size() ? false : end_offset(ends[end_range_end]) == curr_end_offset) {
                    end_range_end++;
                }
                while (start_range_end == starts.size() ? false : start_offset(starts[start_range_end]) == curr_start_offset) {
                    start_range_end++;
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << "first starts in range " << start_range_begin << ":" << start_range_end << " are located at initial offset " << curr_start_offset << endl;
                cerr << "first ends in range " << end_range_begin << ":" << end_range_end << " are located at initial offset " << curr_end_offset << endl;
#endif
                
                // connect the first range of starts or ends to the incoming starts and ends
                
                size_t prev_end_range_begin = end_range_begin;
                size_t prev_start_range_begin = start_range_begin;
                
                bool at_end = (curr_end_offset <= curr_start_offset);
                unordered_map<id_t, unordered_map<size_t, size_t>>* reachable_endpoints;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_ends_from_endpoint;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_starts_from_endpoint;
                vector<size_t>* endpoints;
                size_t* range_begin;
                size_t* range_end;
                size_t* prev_range_begin;
                size_t* curr_offset;
                if (at_end) {
                    reachable_endpoints = &reachable_ends;
                    reachable_starts_from_endpoint = &reachable_starts_from_end;
                    reachable_ends_from_endpoint = &reachable_ends_from_end;
                    endpoints = &ends;
                    range_begin = &end_range_begin;
                    range_end = &end_range_end;
                    prev_range_begin = &prev_end_range_begin;
                    curr_offset = &curr_end_offset;
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "first endpoint is an end" << endl;
#endif
                }
                else {
                    reachable_endpoints = &reachable_starts;
                    reachable_starts_from_endpoint = &reachable_starts_from_start;
                    reachable_ends_from_endpoint = &reachable_ends_from_start;
                    endpoints = &starts;
                    range_begin = &start_range_begin;
                    range_end = &start_range_end;
                    prev_range_begin = &prev_start_range_begin;
                    curr_offset = &curr_start_offset;
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "first endpoint is a start" << endl;
#endif
                }
                
                for (size_t j = *range_begin; j < *range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
#ifdef debug_multipath_mapper_alignment
                        cerr << "identifying end " << incoming_end.first << " as reachable from endpoint " << endpoints->at(j) << endl;
#endif
                        (*reachable_ends_from_endpoint)[endpoints->at(j)].emplace_back(incoming_end.first, incoming_end.second + *curr_offset);
                    }
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
#ifdef debug_multipath_mapper_alignment
                        cerr << "identifying start " << incoming_start.first << " as reachable from endpoint " << endpoints->at(j) << endl;
#endif
                        (*reachable_starts_from_endpoint)[endpoints->at(j)].emplace_back(incoming_start.first, incoming_start.second + *curr_offset);
                    }
                }
                
                
                bool prev_is_end = at_end;
                *range_begin = *range_end;
                prev_offset = *curr_offset;
                if (*range_begin != endpoints->size()) {
                    *curr_offset = endpoint_offset(endpoints->at(*range_begin), at_end);
                    while (*range_end == endpoints->size() ? false : endpoint_offset(endpoints->at(*range_end), at_end) == *curr_offset) {
                        (*range_end)++;
                    }
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << "next range is " << *range_begin << ":" << *range_end << " at offset " << *curr_offset << endl;
#endif
                
                // iterate along ranges of starts or ends in order of their offsets
                
                while (start_range_begin < starts.size() && end_range_begin < ends.size()) {
                    at_end = (curr_end_offset <= curr_start_offset);
                    if (at_end) {
                        reachable_endpoints = &reachable_ends;
                        reachable_starts_from_endpoint = &reachable_starts_from_end;
                        reachable_ends_from_endpoint = &reachable_ends_from_end;
                        endpoints = &ends;
                        range_begin = &end_range_begin;
                        range_end = &end_range_end;
                        prev_range_begin = &prev_end_range_begin;
                        curr_offset = &curr_end_offset;
                    }
                    else {
                        reachable_endpoints = &reachable_starts;
                        reachable_starts_from_endpoint = &reachable_starts_from_start;
                        reachable_ends_from_endpoint = &reachable_ends_from_start;
                        endpoints = &starts;
                        range_begin = &start_range_begin;
                        range_end = &start_range_end;
                        prev_range_begin = &prev_start_range_begin;
                        curr_offset = &curr_start_offset;
                    }
#ifdef debug_multipath_mapper_alignment
                    cerr << "at " << (at_end ? "end" : "start") << "s in range " << *range_begin << ":" << *range_end << endl;
#endif
                    
                    size_t dist_between = *curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    if (prev_is_end) {
#ifdef debug_multipath_mapper_alignment
                        cerr << "looking backwards to ends in range " << prev_end_range_begin << ":" << end_range_begin << endl;
#endif
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_mapper_alignment
                                cerr << "identifying end " << ends[j] << " as reachable from " << (at_end ? "end" : "start") << " " << endpoints->at(k) << endl;
#endif
                                (*reachable_ends_from_endpoint)[endpoints->at(k)].emplace_back(ends[j], dist_between);
                            }
                        }
                    }
                    else {
#ifdef debug_multipath_mapper_alignment
                        cerr << "looking backwards to starts in range " << prev_start_range_begin << ":" << start_range_begin << endl;
#endif
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_mapper_alignment
                                cerr << "identifying start " << starts[j] << " as reachable from " << (at_end ? "end" : "start") << " " << endpoints->at(k) << endl;
#endif
                                (*reachable_starts_from_endpoint)[endpoints->at(k)].emplace_back(starts[j], dist_between);
                            }
                        }
                    }
                    
                    // record the properties of this range
                    *prev_range_begin = *range_begin;
                    prev_is_end = at_end;
                    
                    // advance to the next range
                    *range_begin = *range_end;
                    prev_offset = *curr_offset;
                    if (*range_begin != endpoints->size()) {
                        *curr_offset = endpoint_offset(endpoints->at(*range_begin), at_end);
                        while (*range_end == endpoints->size() ? false : endpoint_offset(endpoints->at(*range_end), at_end) == *curr_offset) {
                            (*range_end)++;
                        }
                    }
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "next range is " << *range_begin << ":" << *range_end << " at offset " << *curr_offset << endl;
#endif
                }
                
                // finish off the list of starts or ends on this node
                
                at_end = (end_range_begin < ends.size());
                if (at_end) {
                    reachable_endpoints = &reachable_ends;
                    reachable_starts_from_endpoint = &reachable_starts_from_end;
                    reachable_ends_from_endpoint = &reachable_ends_from_end;
                    endpoints = &ends;
                    range_begin = &end_range_begin;
                    range_end = &end_range_end;
                    prev_range_begin = &prev_end_range_begin;
                    curr_offset = &curr_end_offset;
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "final endpoint(s) are end(s)" << endl;
#endif
                }
                else {
                    reachable_endpoints = &reachable_starts;
                    reachable_starts_from_endpoint = &reachable_starts_from_start;
                    reachable_ends_from_endpoint = &reachable_ends_from_start;
                    endpoints = &starts;
                    range_begin = &start_range_begin;
                    range_end = &start_range_end;
                    prev_range_begin = &prev_start_range_begin;
                    curr_offset = &curr_start_offset;
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "final endpoint(s) are start(s)" << endl;
#endif
                }
                
                while (*range_begin < endpoints->size()) {
                    
                    size_t dist_between = *curr_offset - prev_offset;
                    
                    if (prev_is_end) {
#ifdef debug_multipath_mapper_alignment
                        cerr << "looking backwards to ends in range " << prev_end_range_begin << ":" << end_range_begin << endl;
#endif
                        for (size_t j = prev_end_range_begin; j < end_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_mapper_alignment
                                cerr << "identifying end " << ends[j] << " as reachable from endpoint " << endpoints->at(k) << endl;
#endif
                                (*reachable_ends_from_endpoint)[endpoints->at(k)].push_back(make_pair(ends[j], dist_between));
                            }
                        }
                    }
                    else {
#ifdef debug_multipath_mapper_alignment
                        cerr << "looking backwards to starts in range " << prev_start_range_begin << ":" << start_range_begin << endl;
#endif
                        for (size_t j = prev_start_range_begin; j < start_range_begin; j++) {
                            for (size_t k = *range_begin; k < *range_end; k++) {
#ifdef debug_multipath_mapper_alignment
                                cerr << "identifying start " << starts[j] << " as reachable from endpoint " << endpoints->at(k) << endl;
#endif
                                (*reachable_starts_from_endpoint)[endpoints->at(k)].push_back(make_pair(starts[j], dist_between));
                            }
                        }
                    }
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "moving to next endpoint range" << endl;
#endif
                    
                    *prev_range_begin = *range_begin;
                    *range_begin = *range_end;
                    prev_offset = *curr_offset;
                    prev_is_end = at_end;
                    
                    if (*range_begin != endpoints->size()) {
                        *curr_offset = endpoint_offset(endpoints->at(*range_begin), at_end);
                        while (*range_end == endpoints->size() ? false : endpoint_offset(endpoints->at(*range_end), at_end) == *curr_offset) {
                            (*range_end)++;
                        }
                    }
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "next range is " << *range_begin << ":" << *range_end << " at offset " << *curr_offset << endl;
#endif
                }
                
                // carry forward the reachability of the last range onto the next nodes
                size_t dist_thru = node_length - *curr_offset;
                
#ifdef debug_multipath_mapper_alignment
                cerr << "carrying forward reachability onto next nodes at distance " << dist_thru << endl;
#endif
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t>& reachable_endpoints_next = (*reachable_endpoints)[next.node->id()];
                    for (size_t j = *prev_range_begin; j < endpoints->size(); j++) {
                        if (reachable_endpoints_next.count(endpoints->at(j))) {
                            reachable_endpoints_next[endpoints->at(j)] = std::min(reachable_endpoints_next[endpoints->at(j)], dist_thru);
                        }
                        else {
                            reachable_endpoints_next[endpoints->at(j)] = dist_thru;
                        }
                    }
                }
            }
            else if (contains_starts || contains_ends) {
                // this nodes contains at least one start or end, but either all starts or all ends
                // record the incoming starts/ends for the starts/ends on this node
                
                unordered_map<id_t, unordered_map<size_t, size_t>>* reachable_endpoints;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_ends_from_endpoint;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_starts_from_endpoint;
                unordered_map<size_t, vector<pair<size_t, size_t>>>* reachable_endpoints_from_endpoint;
                vector<size_t>* endpoints;
                if (contains_ends) {
#ifdef debug_multipath_mapper_alignment
                    cerr << "node " << node_id << " contains only ends of MEMs" << endl;
#endif
                    reachable_endpoints = &reachable_ends;
                    reachable_starts_from_endpoint = &reachable_starts_from_end;
                    reachable_ends_from_endpoint = &reachable_ends_from_end;
                    reachable_endpoints_from_endpoint = &reachable_ends_from_end;
                    endpoints = &exact_match_ends[node_id];
                }
                else {
#ifdef debug_multipath_mapper_alignment
                    cerr << "node " << node_id << " contains only starts of MEMs" << endl;
#endif
                    reachable_endpoints = &reachable_starts;
                    reachable_starts_from_endpoint = &reachable_starts_from_start;
                    reachable_ends_from_endpoint = &reachable_ends_from_start;
                    reachable_endpoints_from_endpoint = &reachable_starts_from_start;
                    endpoints = &exact_match_starts[node_id];
                }
                
                // the starts/ends coming into this node from outside
                size_t range_begin = 0;
                size_t range_end = 1;
                size_t curr_offset = endpoint_offset(endpoints->at(range_begin), contains_ends);
                size_t prev_offset = curr_offset;
                // find the range of endpoints that are at the first offset
                while (range_end == endpoints->size() ? false : endpoint_offset(endpoints->at(range_end), contains_ends) == curr_offset) {
                    range_end++;
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << "initial range " << range_begin << ":" << range_end << " is at offset " << curr_offset << endl;
#endif
                
                // connect the range to the incoming starts/ends
                for (size_t j = range_begin; j < range_end; j++) {
                    for (const pair<size_t, size_t>& incoming_start : reachable_starts[node_id]) {
#ifdef debug_multipath_mapper_alignment
                        cerr << "identifying start " << incoming_start.first << " as reachable from " << (contains_ends ? "end" : "start") << " " << endpoints->at(j) << endl;
#endif
                        (*reachable_starts_from_endpoint)[endpoints->at(j)].emplace_back(incoming_start.first, incoming_start.second + curr_offset);
                    }
                    for (const pair<size_t, size_t>& incoming_end : reachable_ends[node_id]) {
#ifdef debug_multipath_mapper_alignment
                        cerr << "identifying end " << incoming_end.first << " as reachable from " << (contains_ends ? "end" : "start") << " " << endpoints->at(j) << endl;
#endif
                        (*reachable_ends_from_endpoint)[endpoints->at(j)].emplace_back(incoming_end.first, incoming_end.second + curr_offset);
                    }
                }
                
                // the reachable endpoints internal to this node
                size_t prev_range_begin = range_begin;
                range_begin = range_end;
                while (range_begin < endpoints->size()) {
                    // find the range of endpoints at this offset
                    prev_offset = curr_offset;
                    curr_offset = endpoint_offset(endpoints->at(range_begin), contains_ends);
                    while (range_end == endpoints->size() ? false : endpoint_offset(endpoints->at(range_end), contains_ends) == curr_offset) {
                        range_end++;
                    }
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "next range " << range_begin << ":" << range_end << " is at offset " << curr_offset << endl;
#endif
                    
                    size_t dist_between = curr_offset - prev_offset;
                    
                    // connect this range to the previous range
                    for (size_t j = range_begin; j < range_end; j++) {
                        for (size_t k = prev_range_begin; k < range_begin; k++) {
#ifdef debug_multipath_mapper_alignment
                            cerr << "identifying " << (contains_ends ? "end" : "start") << " " << endpoints->at(k) << " as reachable from " << (contains_ends ? "end" : "start") << " " << endpoints->at(j) << endl;
#endif
                            (*reachable_endpoints_from_endpoint)[endpoints->at(j)].push_back(make_pair(endpoints->at(k), dist_between));
                        }
                    }
                    prev_range_begin = range_begin;
                    range_begin = range_end;
                }
                
                // this node contains at least one endpoint of a MEM, so carry forward the reachability of all
                // endpoints at the final offset onto the next nodes
                
                size_t dist_thru = node_length - curr_offset;
                
#ifdef debug_multipath_mapper_alignment
                cerr << "carrying forward reachability onto next nodes at distance " << dist_thru << endl;
#endif
                
                for (NodeTraversal next : nexts) {
                    
                    unordered_map<size_t, size_t>& reachable_endpoints_next = (*reachable_endpoints)[next.node->id()];
                    for (size_t j = prev_range_begin; j < endpoints->size(); j++) {
                        if (reachable_endpoints_next.count(endpoints->at(j))) {
                            reachable_endpoints_next[endpoints->at(j)] = std::min(reachable_endpoints_next[endpoints->at(j)], dist_thru);
                        }
                        else {
                            reachable_endpoints_next[endpoints->at(j)] = dist_thru;
                        }
                    }
                }
            }
            else {
                // this node doesn't contain the start or end of any MEM, so we carry forward the reachability
                // into this node onto the next nodes
                
#ifdef debug_multipath_mapper_alignment
                cerr << "node " << node_id << " does not contain starts or ends of MEMs, carrying forward reachability" << endl;
#endif
                
                for (NodeTraversal next : nexts) {
                    unordered_map<size_t, size_t>& reachable_ends_next = reachable_ends[next.node->id()];
                    for (const pair<size_t, size_t>& reachable_end : reachable_ends[node_id]) {
                        size_t dist_thru = reachable_end.second + node_length;
#ifdef debug_multipath_mapper_alignment
                        cerr << "\tend " << reachable_end.first << " at dist " << dist_thru << " to node " << next.node->id() << endl;
#endif
                        if (reachable_ends_next.count(reachable_end.first)) {
                            reachable_ends_next[reachable_end.first] = std::min(reachable_ends_next[reachable_end.first],
                                                                                dist_thru);
                        }
                        else {
                            reachable_ends_next[reachable_end.first] = dist_thru;
                        }
                    }
                    
                    unordered_map<size_t, size_t>& reachable_starts_next = reachable_starts[next.node->id()];
                    for (const pair<size_t, size_t>& reachable_start : reachable_starts[node_id]) {
                        size_t dist_thru = reachable_start.second + node_length;
#ifdef debug_multipath_mapper_alignment
                        cerr << "\tstart " << reachable_start.first << " at dist " << dist_thru << " to node " << next.node->id() << endl;
#endif
                        if (reachable_starts_next.count(reachable_start.first)) {
                            reachable_starts_next[reachable_start.first] = std::min(reachable_starts_next[reachable_start.first],
                                                                                    dist_thru);
                        }
                        else {
                            reachable_starts_next[reachable_start.first] = dist_thru;
                        }
                    }
                }
            }
        }
        
#ifdef debug_multipath_mapper_alignment
        cerr << "final reachability:" << endl;
        cerr << "\tstarts from starts" << endl;
        for (const auto& record : reachable_starts_from_start) {
            cerr << "\t\t" << record.first << ":";
            for (const auto& endpoint : record.second) {
                cerr << " " << endpoint.first << "(" << endpoint.second << ")";
            }
            cerr << endl;
        }
        cerr << "\tends from starts" << endl;
        for (const auto& record : reachable_ends_from_start) {
            cerr << "\t\t" << record.first << ":";
            for (const auto& endpoint : record.second) {
                cerr << " " << endpoint.first << "(" << endpoint.second << ")";
            }
            cerr << endl;
        }
        cerr << "\tstarts from ends" << endl;
        for (const auto& record : reachable_starts_from_end) {
            cerr << "\t\t" << record.first << ":";
            for (const auto& endpoint : record.second) {
                cerr << " " << endpoint.first << "(" << endpoint.second << ")";
            }
            cerr << endl;
        }
        cerr << "\tends from ends" << endl;
        for (const auto& record : reachable_ends_from_end) {
            cerr << "\t\t" << record.first << ":";
            for (const auto& endpoint : record.second) {
                cerr << " " << endpoint.first << "(" << endpoint.second << ")";
            }
            cerr << endl;
        }
        cerr << "setting up structure of MEM graph" << endl;
#endif
        
        // now we have the reachability information for the start and end of every MEM in the graph. we
        // will use this to navigate between the MEMs in a way that respects graph reachability so that this
        // phase of the algorithm only needs to pay attention to read colinearity and transitive reducibility
        
        vector<unordered_map<size_t, size_t>> noncolinear_shells(match_nodes.size());
        
        // tuples of (overlap size, index onto, index from, dist)
        vector<tuple<size_t, size_t, size_t, size_t>> confirmed_overlaps;
        
        for (size_t i = 0; i < graph.node_size(); i++) {
            id_t node_id = graph.node(i).id();
            
#ifdef debug_multipath_mapper_alignment
            cerr << "looking for edges for starts on node " << node_id << endl;
#endif
            
            if (!exact_match_starts.count(node_id)) {
#ifdef debug_multipath_mapper_alignment
                cerr << "there are no starts on this node" << endl;
#endif
                continue;
            }
            
            // keep track of the starts that are located at the same offset at earlier positions in the vector of
            // of starts (used later in the overlap finding step)
            vector<size_t> colocated_starts;
            
            vector<size_t>& starts = exact_match_starts[node_id];
            vector<size_t>& ends = exact_match_ends[node_id];
            // index of the next end that is past the start we are on
            size_t next_end_idx = 0;
            // sentinel that will never be equal to the first offset
            size_t curr_start_offset = numeric_limits<size_t>::max();
            
            for (size_t start_idx = 0; start_idx < starts.size(); start_idx++) {
                // traverse all of the reachable starts to find the adjacent ends that might be colinear
                
                size_t start = starts[start_idx];
#ifdef debug_multipath_mapper_alignment
                cerr << "searching backward from start " << start << endl;
#endif
                
                ExactMatchNode& start_node = match_nodes[start];
                unordered_map<size_t, size_t>& noncolinear_shell = noncolinear_shells[start];
                
                // pairs of (dist, index)
                priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> start_queue;
                priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> end_queue;
                start_queue.emplace(0, start);
                
                unordered_set<size_t> traversed_start;
                
                while (!start_queue.empty()) {
                    pair<size_t, size_t> start_here = start_queue.top();
                    start_queue.pop();
                    if (traversed_start.count(start_here.second)) {
                        continue;
                    }
                    traversed_start.insert(start_here.second);
#ifdef debug_multipath_mapper_alignment
                    cerr << "traversing initial start " << start_here.second << " at distance " << start_here.first << endl;
#endif
                    
                    // the minimum distance to each of the starts or ends this can reach is the sum of the min distance
                    // between them and the distance already traversed
                    for (const pair<size_t, size_t>& end : reachable_ends_from_start[start_here.second]) {
                        end_queue.emplace(start_here.first + end.second, end.first);
#ifdef debug_multipath_mapper_alignment
                        cerr << "found reachable end " << end.first << " at distance " << start_here.first + end.second << endl;
#endif
                    }
                    
                    for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here.second]) {
                        start_queue.emplace(start_here.first + start_next.second, start_next.first);
                    }
                }
                
                // now we've traversed all of the starts, we have the set of ends that can be reached
                // without passing another end
                
                unordered_set<size_t> traversed_end;
                
                while (!end_queue.empty()) {
                    size_t candidate_end = end_queue.top().second;
                    size_t candidate_dist = end_queue.top().first;
                    end_queue.pop();
                    
                    if (traversed_end.count(candidate_end)) {
                        continue;
                    }
#ifdef debug_multipath_mapper_alignment
                    cerr << "considering end " << candidate_end << " as candidate for edge of dist " << candidate_dist << endl;
#endif
                    traversed_end.insert(candidate_end);
                    
                    ExactMatchNode& candidate_end_node = match_nodes[candidate_end];
                    
                    if (candidate_end_node.end <= start_node.begin) {
                        // these MEMs are read colinear and graph reachable, so connect them
                        candidate_end_node.edges.push_back(make_pair(start, candidate_dist));
                        
#ifdef debug_multipath_mapper_alignment
                        cerr << "connection is read colinear, adding edge" << endl;
#endif
                        
                        // skip to the predecessor's noncolinear shell, whose connections might not be blocked by
                        // this connection
                        for (const pair<size_t, size_t>& shell_pred : noncolinear_shells[candidate_end]) {
#ifdef debug_multipath_mapper_alignment
                            cerr << "enqueueing " << shell_pred.first << " at dist " << shell_pred.second + candidate_dist << " from noncolinear shell" << endl;
#endif
                            end_queue.emplace(candidate_dist + shell_pred.second, shell_pred.first);
                        }
                    }
                    else if (start_node.end > candidate_end_node.end && start_node.begin > candidate_end_node.begin) {
                        // the MEM can be made colinear by removing an overlap, which will not threaten reachability
                        size_t overlap = candidate_end_node.end - start_node.begin;
                        confirmed_overlaps.emplace_back(overlap, start, candidate_end, candidate_dist + overlap);
                        
#ifdef debug_multipath_mapper_alignment
                        cerr << "connection is overlap colinear, recording to add edge later" << endl;
#endif
                        
                        // the end of this node might not actually block connections since it's going to intersect the middle of the node
                        // so we need to find predecessors to this end too
                        
                        // add any ends directly reachable from the end
                        for (const pair<size_t, size_t>& exposed_end : reachable_ends_from_end[candidate_end]) {
                            end_queue.emplace(candidate_dist + exposed_end.second, exposed_end.first);
#ifdef debug_multipath_mapper_alignment
                            cerr << "found reachable exposed end " << exposed_end.first << " at distance " << candidate_dist + exposed_end.second << endl;
#endif
                        }
                        
                        // traverse through any exposes starts to see if we can find other exposed ends
                        priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> exposed_start_queue;
                        unordered_set<size_t> traversed_exposed_start;
                        
                        // inialize the queue with the directly reachable exposed starts
                        for (const pair<size_t, size_t>& exposed_start : reachable_starts_from_end[candidate_end]) {
#ifdef debug_multipath_mapper_alignment
                            cerr << "initializing exposed start traversal with " << exposed_start.first << " at distance " << candidate_dist + exposed_start.second << endl;
#endif
                            exposed_start_queue.emplace(candidate_dist + exposed_start.second, exposed_start.first);
                        }
                        
                        while (!exposed_start_queue.empty()) {
                            pair<size_t, size_t> start_here = exposed_start_queue.top();
                            exposed_start_queue.pop();
                            if (traversed_exposed_start.count(start_here.second)) {
                                continue;
                            }
                            traversed_exposed_start.insert(start_here.second);
#ifdef debug_multipath_mapper_alignment
                            cerr << "traversing exposed start " << start_here.second << " at distance " << start_here.first << endl;
#endif
                            
                            // the minimum distance to each of the starts or ends this can reach is the sum of the min distance
                            // between them and the distance already traversed
                            for (const pair<size_t, size_t>& end : reachable_ends_from_start[start_here.second]) {
                                end_queue.emplace(start_here.first + end.second, end.first);
#ifdef debug_multipath_mapper_alignment
                                cerr << "found reachable exposed end " << end.first << " at distance " << start_here.first + end.second << endl;
#endif
                            }
                            
                            for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here.second]) {
                                exposed_start_queue.emplace(start_here.first + start_next.second, start_next.first);
                            }
                        }
                        
                        // also skip to the predecessor's noncolinear shell, whose connections might not be blocked by
                        // this connection
                        for (const pair<size_t, size_t>& shell_pred : noncolinear_shells[candidate_end]) {
#ifdef debug_multipath_mapper_alignment
                            cerr << "enqueueing " << shell_pred.first << " at dist " << candidate_dist + shell_pred.second << " from noncolinear shell" << endl;
#endif
                            end_queue.emplace(candidate_dist + shell_pred.second, shell_pred.first);
                        }
                    }
                    else {
                        // these MEMs are noncolinear, so add this predecessor to the noncolinear shell
                        if (noncolinear_shell.count(candidate_end)) {
                            noncolinear_shell[candidate_end] = std::min(candidate_dist + (start_node.end - start_node.begin),
                                                                        noncolinear_shell[candidate_end]);
                        }
                        else {
                            noncolinear_shell[candidate_end] = candidate_dist + (start_node.end - start_node.begin);
                        }
                        
#ifdef debug_multipath_mapper_alignment
                        cerr << "connection is noncolinear, add to shell at dist " << candidate_dist + (start_node.end - start_node.begin) << " and continue to search backwards" << endl;
#endif
                        
                        // there is no connection to block further connections back, so any of this MEMs
                        // predecessors could still be colinear
                        
                        // find the ends that can reach it directly
                        for (const pair<size_t, size_t>& pred_end : reachable_ends_from_end[candidate_end]) {
                            end_queue.emplace(candidate_dist + pred_end.second, pred_end.first);
#ifdef debug_multipath_mapper_alignment
                            cerr << "found reachable end " << pred_end.first << " at distance " << candidate_dist + pred_end.second << endl;
#endif
                        }
                        
                        // traverse backward through any starts to find more ends that can reach this MEM
                        priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>, std::greater<pair<size_t, size_t>>> pred_start_queue;
                        
                        // initialize the queue with the immediate start neighbors
                        for (const pair<size_t, size_t>& pred_start : reachable_starts_from_end[candidate_end]) {
                            pred_start_queue.emplace(candidate_dist + pred_start.second, pred_start.first);
                        }
                        
                        unordered_set<size_t> pred_traversed;
                        
                        // traverse backwards through starts, stopping at any ends
                        while (!pred_start_queue.empty()) {
                            size_t start_here = pred_start_queue.top().second;
                            size_t start_dist = pred_start_queue.top().first;
                            pred_start_queue.pop();
                            if (pred_traversed.count(start_here)) {
                                continue;
                            }
                            pred_traversed.insert(start_here);
                            
#ifdef debug_multipath_mapper_alignment
                            cerr << "traversing predecessor start " << start_here << " at distance " << start_dist << endl;
#endif
                            
                            for (const pair<size_t, size_t>& pred_end : reachable_ends_from_start[start_here]) {
                                end_queue.emplace(start_dist + pred_end.second, pred_end.first);
#ifdef debug_multipath_mapper_alignment
                                cerr << "found reachable end " << pred_end.first << " at distance " << candidate_dist + pred_end.second << endl;
#endif
                            }
                            for (const pair<size_t, size_t>& start_next : reachable_starts_from_start[start_here]) {
                                pred_start_queue.emplace(start_dist + start_next.second, start_next.first);
#ifdef debug_multipath_mapper_alignment
                                cerr << "found intermediate start " << start_next.first << " at distance " << start_dist + start_next.second << endl;
#endif
                            }
                        }
                    }
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << "walking path to look for overlaps" << endl;
#endif
                
                size_t prev_start_offset = curr_start_offset;
                curr_start_offset = start_offset(start);
                // update the list of starts at this offset earlier in the starts vector
                if (curr_start_offset != prev_start_offset) {
                    colocated_starts.clear();
                }
                colocated_starts.push_back(start);
                
                // move the next end pointer to the one immediately following this start on the node
                while (next_end_idx >= ends.size() ? false : end_offset(ends[next_end_idx]) <= curr_start_offset) {
                    next_end_idx++;
                }
                
                Path& match_path = match_nodes[start].path;
                // the starts that are on this path
                unordered_set<size_t> path_starts(colocated_starts.begin(), colocated_starts.end());
                // records of (node_idx, overlap length)
                vector<pair<size_t, size_t>> overlap_candidates;
                
                if (match_path.mapping_size() == 1) {
                    // TODO: this edge case is a little duplicative, probably could merge
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "path is one mapping long" << endl;
#endif
                    
                    size_t final_offset = end_offset(start);
                    // record which starts are on the path on this node
                    for (size_t path_start_idx = start_idx + 1;
                         path_start_idx >= starts.size() ? false : start_offset(starts[path_start_idx]) < final_offset;
                         path_start_idx++) {
                        
                        path_starts.insert(starts[path_start_idx]);
                        
                    }
                    // record which ends are on the path on this node
                    for (size_t path_end_idx = next_end_idx; path_end_idx < ends.size(); path_end_idx++) {
                        size_t end_offset_here = end_offset(ends[path_end_idx]);
                        if (end_offset_here < final_offset) {
                            overlap_candidates.emplace_back(ends[path_end_idx], end_offset_here - curr_start_offset);
                        }
                        else {
                            break;
                        }
                    }
                }
                else {
#ifdef debug_multipath_mapper_alignment
                    cerr << "path is multiple mappings long" << endl;
#endif
                    
                    // record which starts are on the path on the first node
                    for (size_t path_start_idx = start_idx + 1; path_start_idx < starts.size(); path_start_idx++) {
                        path_starts.insert(starts[path_start_idx]);
                    }
                    // record which ends are on the path on the first node
                    for (size_t path_end_idx = next_end_idx; path_end_idx < ends.size(); path_end_idx++) {
                        overlap_candidates.emplace_back(ends[path_end_idx], end_offset(ends[path_end_idx]) - curr_start_offset);
                    }
                    size_t traversed_length = mapping_from_length(match_path.mapping(0));
                    
                    for (size_t j = 1; j + 1 < match_path.mapping_size(); j++) {
                        id_t path_node_id = match_path.mapping(j).position().node_id();
                        // record which starts are on the path on this node
                        for (size_t path_start : exact_match_starts[path_node_id]) {
                            path_starts.insert(path_start);
                        }
                        // record which ends are on the path on this node
                        for (size_t path_end : exact_match_ends[path_node_id]) {
                            overlap_candidates.emplace_back(path_end, end_offset(path_end) + traversed_length);
                        }
                        
                        traversed_length += mapping_from_length(match_path.mapping(j));
                    }
                    
                    id_t final_node_id = match_path.mapping(match_path.mapping_size() - 1).position().node_id();
                    vector<size_t>& final_starts = exact_match_starts[final_node_id];
                    vector<size_t>& final_ends = exact_match_ends[final_node_id];
                    
                    size_t final_offset = end_offset(start);
                    // record which starts are on the path on the last node
                    for (size_t path_start_idx = 0;
                         path_start_idx >= final_starts.size() ? false : start_offset(final_starts[path_start_idx]) < final_offset;
                         path_start_idx++) {
                        
                        path_starts.insert(final_starts[path_start_idx]);
                        
                    }
                    // record which ends are on the path on the last node
                    for (size_t path_end_idx = 0; path_end_idx < final_ends.size(); path_end_idx++) {
                        size_t end_offset_here = end_offset(final_ends[path_end_idx]);
                        if (end_offset_here < final_offset) {
                            overlap_candidates.emplace_back(final_ends[path_end_idx], end_offset_here + traversed_length);
                        }
                        else {
                            break;
                        }
                    }
                }
                
                for (const pair<size_t, size_t>& overlap_candidate : overlap_candidates) {
#ifdef debug_multipath_mapper_alignment
                    cerr << "considering candidate overlap from " << overlap_candidate.first << " at dist " << overlap_candidate.second << endl;
#endif
                    
                    if (path_starts.count(overlap_candidate.first)) {
                        // the start of this MEM is also on the path, so this can't be an overhanging overlap
                        continue;
                    }
                    
                    ExactMatchNode& overlap_node = match_nodes[overlap_candidate.first];
                    
                    // how much do the paths overlap?
                    size_t overlap = overlap_candidate.second;
                    
                    // are the matches read colinear after removing the overlap?
                    if (start_node.begin + overlap >= overlap_node.end) {
#ifdef debug_multipath_mapper_alignment
                        cerr << "confirmed overlap colinear with overlap of " << overlap << endl;
#endif
                        confirmed_overlaps.emplace_back(overlap, start, overlap_candidate.first, 0);
                    }
                    else if (overlap_node.begin < start_node.begin && overlap_node.end < start_node.end) {
#ifdef debug_multipath_mapper_alignment
                        cerr << "confirmed overlap colinear with longer read overlap of " << overlap_node.end - start_node.begin << endl;
#endif
                        // there is still an even longer read overlap we need to remove
                        size_t read_overlap = overlap_node.end - start_node.begin;
                        confirmed_overlaps.emplace_back(read_overlap, start, overlap_candidate.first, read_overlap - overlap);
                    }
                    else {
#ifdef debug_multipath_mapper_alignment
                        cerr << "not colinear even with overlap, adding to non-colinear shell at distance " << overlap_candidate.second << endl;
#endif
                        // the overlapping node is still not reachable so it is in the noncolinear shell of this node
                        noncolinear_shell[overlap_candidate.first] = (start_node.end - start_node.begin) - overlap_candidate.second;
                    }
                }
            }
        }
        
#ifdef debug_multipath_mapper_alignment
        cerr << "breaking nodes at overlap edges (" << confirmed_overlaps.size() << " times)" << endl;
#endif
        
        // now we've found all overlap edges, so we can add them into the graph in an order such that they don't
        // conflict (note that all overlap are from an earlier node onto a later one, so we don't need to worry
        // about overlaps coming in from both directions)
        
        // sort in descending order of overlap length and group by the node that is being cut among overlaps of same length
        std::sort(confirmed_overlaps.begin(), confirmed_overlaps.end(),
                  std::greater<tuple<size_t, size_t, size_t, size_t>>());
        
        // keep track of whether another node is holding the suffix of one of the original nodes because of a split
        unordered_map<size_t, size_t> node_with_suffix;
        
        // split up each node with an overlap edge onto it
        auto iter = confirmed_overlaps.begin();
        while (iter != confirmed_overlaps.end()) {
            // find the range of overlaps that want to cut this node at the same place
            auto iter_range_end = iter;
            while (get<0>(*iter_range_end) == get<0>(*iter) && get<1>(*iter_range_end) == get<1>(*iter)) {
                iter_range_end++;
                if (iter_range_end == confirmed_overlaps.end()) {
                    break;
                }
            }
            
#ifdef debug_multipath_mapper_alignment
            cerr << "performing an overlap split onto " << get<1>(*iter) << " of length " << get<0>(*iter) << endl;
#endif
            
            
            ExactMatchNode* onto_node = &match_nodes[get<1>(*iter)];
            
#ifdef debug_multipath_mapper_alignment
            cerr << "before splitting:" << endl;
            cerr << "onto node:" << endl << "\t";
            for (auto node_iter = onto_node->begin; node_iter != onto_node->end; node_iter++) {
                cerr << *node_iter;
            }
            cerr << endl << "\t" << pb2json(onto_node->path) << endl;
#endif
            
            // store the full path and remove it from the node
            Path full_path = std::move(onto_node->path);
            onto_node->path.Clear();
            
            // add mappings from the path until reaching the overlap point
            int64_t remaining = get<0>(*iter);
            int64_t mapping_idx = 0;
            int64_t mapping_len = mapping_from_length(full_path.mapping(mapping_idx));
            while (remaining >= mapping_len) {
                *onto_node->path.add_mapping() = full_path.mapping(mapping_idx);
                remaining -= mapping_len;
                mapping_idx++;
                if (mapping_idx == full_path.mapping_size()) {
                    break;
                }
                mapping_len = mapping_from_length(full_path.mapping(mapping_idx));
            }
            
            if (mapping_idx == full_path.mapping_size() && !remaining) {
                // TODO: isn't this case covered by taking the entire range of splits at the same place?
                
                // the overlap covered the match, so connect it to the onto node's successors
                // rather than splitting it into two nodes
                
                while (iter != iter_range_end) {
                    for (const pair<size_t, size_t> edge : onto_node->edges) {
                        match_nodes[get<2>(*iter)].edges.emplace_back(edge.first, edge.second + get<3>(*iter));
                    }
                    iter++;
                }
            }
            else {
                // the overlap was in the middle of the match, so split the onto node into a
                // prefix and suffix
                
                // make a new node to hold the suffix of the match
                size_t suffix_idx = match_nodes.size();
                match_nodes.emplace_back();
                ExactMatchNode& suffix_node = match_nodes.back();
                
                // get the pointer from the onto node back in case the vector reallocated
                onto_node = &match_nodes[get<1>(*iter)];
                
                // divide up the match on the read
                suffix_node.end = onto_node->end;
                suffix_node.begin = onto_node->begin + get<0>(*iter);
                onto_node->end = suffix_node.begin;
                
                // transfer the outgoing edges onto the new node
                suffix_node.edges = std::move(onto_node->edges);
                
                // clear the old edges and add a single edge to the suffix
                onto_node->edges.clear();
                onto_node->edges.emplace_back(suffix_idx, 0);
                
                // keep track of the relationship of suffix nodes to original nodes
                if (!node_with_suffix.count(get<1>(*iter))) {
                    // since we take longest overlaps first, only the first split onto a node will change
                    // which node contains the suffix of the original node
                    node_with_suffix[get<1>(*iter)] = suffix_idx;
                }
                
                if (remaining) {
                    // the overlap point is in the middle of a node, need to split a mapping
                    
                    const Mapping& split_mapping = full_path.mapping(mapping_idx);
                    
                    // add the prefix of the mapping to the original node
                    Mapping* prefix_split = onto_node->path.add_mapping();
                    prefix_split->mutable_position()->set_node_id(split_mapping.position().node_id());
                    prefix_split->mutable_position()->set_offset(split_mapping.position().offset());
                    Edit* prefix_edit = prefix_split->add_edit();
                    prefix_edit->set_from_length(remaining);
                    prefix_edit->set_to_length(remaining);
                    
                    // add the suffix of the mapping to the new node
                    Mapping* suffix_split = suffix_node.path.add_mapping();
                    suffix_split->mutable_position()->set_node_id(split_mapping.position().node_id());
                    suffix_split->mutable_position()->set_offset(split_mapping.position().offset() + remaining);
                    Edit* suffix_edit = suffix_split->add_edit();
                    suffix_edit->set_from_length(mapping_len - remaining);
                    suffix_edit->set_to_length(mapping_len - remaining);
                    
                    mapping_idx++;
                }
                
                // add the remaining mappings to the suffix node
                for (; mapping_idx < full_path.mapping_size(); mapping_idx++) {
                    *suffix_node.path.add_mapping() = full_path.mapping(mapping_idx);
                }
                
#ifdef debug_multipath_mapper_alignment
                cerr << "after splitting:" << endl;
                cerr << "onto node:" << endl << "\t";
                for (auto node_iter = onto_node->begin; node_iter != onto_node->end; node_iter++) {
                    cerr << *node_iter;
                }
                cerr << endl << "\t" << pb2json(onto_node->path) << endl;
                cerr << "suffix node:" << endl << "\t";
                for (auto node_iter = suffix_node.begin; node_iter != suffix_node.end; node_iter++) {
                    cerr << *node_iter;
                }
                cerr << endl << "\t" << pb2json(suffix_node.path) << endl;
#endif
                
                while (iter != iter_range_end) {
                    // index of the node that contains the end of the original node we recorded the overlap from
                    size_t splitting_idx = node_with_suffix.count(get<2>(*iter)) ? node_with_suffix[get<2>(*iter)] : get<2>(*iter);
                    
#ifdef debug_multipath_mapper_alignment
                    cerr << "adding an overlap edge from node " << splitting_idx << " at distance " << get<3>(*iter) << endl;
                    cerr << "\t";
                    for (auto node_iter = match_nodes[splitting_idx].begin; node_iter != match_nodes[splitting_idx].end; node_iter++) {
                        cerr << *node_iter;
                    }
                    cerr << endl;
#endif
                    
                    // get the next node that overlaps onto the other node at this index and add the overlap edge
                    match_nodes[splitting_idx].edges.emplace_back(suffix_idx, get<3>(*iter));
                    
                    iter++;
                }
            }
        }
    }
    
    // Kahn's algorithm
    void MultipathAlignmentGraph::topological_sort(vector<size_t>& order_out) {
        order_out.resize(match_nodes.size());
        
        vector<size_t> in_degree(match_nodes.size());
        
        for (const ExactMatchNode& match_node : match_nodes) {
            for (const pair<size_t, size_t>& edge : match_node.edges) {
                in_degree[edge.first]++;
            }
        }
        
        list<size_t> source_queue;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            if (in_degree[i] == 0) {
                source_queue.push_back(i);
            }
        }
        
        size_t next = 0;
        while (!source_queue.empty()) {
            size_t src = source_queue.front();
            source_queue.pop_front();
            
            for (const pair<size_t, size_t>& edge : match_nodes[src].edges) {
                in_degree[edge.first]--;
                if (in_degree[edge.first] == 0) {
                    source_queue.push_back(edge.first);
                }
            }
            
            order_out[next] = src;
            next++;
        }
    }
    
    void MultipathAlignmentGraph::reorder_adjacency_lists(const vector<size_t>& order) {
        vector<vector<pair<size_t, size_t>>> reverse_graph(match_nodes.size());
        for (size_t i = 0; i < match_nodes.size(); i++) {
            for (const pair<size_t, size_t>& edge : match_nodes[i].edges) {
                reverse_graph[edge.first].emplace_back(i, edge.second);
            }
        }
        for (ExactMatchNode& match_node : match_nodes) {
            size_t out_degree = match_node.edges.size();
            match_node.edges.clear();
            match_node.edges.reserve(out_degree);
        }
        for (size_t i : order) {
            for (const pair<size_t, size_t>& edge : reverse_graph[i]) {
                match_nodes[edge.first].edges.emplace_back(i, edge.second);
            }
        }
    }
    
    void MultipathAlignmentGraph::remove_transitive_edges(const vector<size_t>& topological_order) {
        // algorithm assumes edges are also sorted in topological order, which guarantees that we will
        // traverse a path that reveals an edge as transitive before actually traversing the transitive edge
        reorder_adjacency_lists(topological_order);
        
        for (size_t i : topological_order) {
            vector<pair<size_t, size_t>>& edges = match_nodes[i].edges;
            
            // if there is only one edge out of a node, that edge can never be transitive
            // (this optimization covers most cases)
            if (edges.size() <= 1) {
                continue;
            }
            
            vector<bool> keep(edges.size(), true);
            unordered_set<size_t> traversed;
            
            for (size_t j = 0; j < edges.size(); j++) {
                const pair<size_t, size_t>& edge = edges[j];
                if (traversed.count(edge.first)) {
                    // we can reach the target of this edge by another path, so it is transitive
                    keep[j] = false;
                    continue;
                }
                
                // DFS to mark all reachable nodes from this edge
                vector<size_t> stack{edge.first};
                traversed.insert(edge.first);
                while (!stack.empty()) {
                    size_t idx = stack.back();
                    stack.pop_back();
                    for (const pair<size_t, size_t>& edge_from : match_nodes[idx].edges) {
                        if (!traversed.count(edge_from.first)) {
                            stack.push_back(edge_from.first);
                            traversed.insert(edge_from.first);
                        }
                    }
                }
            }
            
            // remove the transitive edges we found
            size_t next_idx = 0;
            for (size_t j = 0; j < edges.size(); j++) {
                if (keep[j] && j != next_idx) {
                    edges[next_idx] = edges[j];
                    next_idx++;
                }
                else if (keep[j]) {
                    next_idx++;
                }
            }
            edges.resize(next_idx);
        }
    }
    
    MultipathAlignmentGraph::~MultipathAlignmentGraph() {
        
    }
    
    void MultipathAlignmentGraph::prune_to_high_scoring_paths(const Alignment& alignment, const BaseAligner* aligner,
                                                              double max_suboptimal_score_ratio, const vector<size_t>& topological_order) {
        
        if (match_nodes.empty()) {
            return;
        }
        
        unordered_map<pair<size_t, size_t>, int32_t> edge_weights;

        vector<int32_t> node_weights(match_nodes.size());
    
        // TODO: is the lower bound too strict?
        
        // compute the weight of edges and node matches
        for (size_t i = 0; i < match_nodes.size(); i++) {
            ExactMatchNode& from_node = match_nodes[i];
            node_weights[i] = aligner->match * (from_node.end - from_node.begin)
                              + aligner->full_length_bonus * ((from_node.begin == alignment.sequence().begin())
                                                             + (from_node.end == alignment.sequence().end()));
            
            for (const pair<size_t, size_t>& edge : from_node.edges) {
                ExactMatchNode& to_node = match_nodes[edge.first];
                
                int64_t graph_dist = edge.second;
                int64_t read_dist = to_node.begin - from_node.end;
                
                if (read_dist > graph_dist) {
                    // the read length in between the MEMs is longer than the distance, suggesting a read insert
                    // and potentially another mismatch on the other end
                    int64_t gap_length = read_dist - graph_dist;
                    edge_weights[make_pair(i, edge.first)] = -(gap_length - 1) * aligner->gap_extension - aligner->gap_open
                                                             - (graph_dist > 0) * aligner->mismatch;
                }
                else if (read_dist < graph_dist) {
                    // the read length in between the MEMs is shorter than the distance, suggesting a read deletion
                    // and potentially another mismatch on the other end
                    int64_t gap_length = graph_dist - read_dist;
                    edge_weights[make_pair(i, edge.first)] = -(gap_length - 1) * aligner->gap_extension - aligner->gap_open
                                                             - (read_dist > 0) * aligner->mismatch;
                }
                else {
                    // the read length in between the MEMs is the same as the distance, suggesting a pure mismatch
                    edge_weights[make_pair(i, edge.first)] = -((graph_dist > 0) + (graph_dist > 1)) * aligner->mismatch;
                }
            }
        }
        
        vector<int32_t> forward_scores = node_weights;
        vector<int32_t> backward_scores = node_weights;
        
        // forward DP
        for (int64_t i = 0; i < topological_order.size(); i++) {
            size_t idx = topological_order[i];
            int32_t from_score = forward_scores[idx];
            for (const pair<size_t, size_t>& edge : match_nodes[idx].edges) {
                forward_scores[edge.first] = std::max(forward_scores[edge.first],
                                                      node_weights[edge.first] + from_score + edge_weights[make_pair(idx, edge.first)]);
            }
        }
        
        // backward DP
        for (int64_t i = topological_order.size() - 1; i >= 0; i--) {
            size_t idx = topological_order[i];
            int32_t score_here = node_weights[idx];
            for (const pair<size_t, size_t>& edge : match_nodes[idx].edges) {
                backward_scores[idx] = std::max(backward_scores[idx],
                                                score_here + backward_scores[edge.first] + edge_weights[make_pair(idx, edge.first)]);
            }
        }
        
        // compute the minimum score we will require of a node or edge
        int32_t min_path_score = *std::max_element(forward_scores.begin(), forward_scores.end()) / max_suboptimal_score_ratio;
        
        // use forward-backward to find nodes/edges on some path with a score above the minimum
        unordered_set<size_t> keep_nodes;
        unordered_set<pair<size_t, size_t>> keep_edges;
        vector<size_t> removed_in_prefix(match_nodes.size() + 1, 0);
        for (size_t i = 0; i < match_nodes.size(); i++) {
            if (forward_scores[i] + backward_scores[i] - node_weights[i] >= min_path_score) {
                keep_nodes.insert(i);
                for (const pair<size_t, size_t>& edge : match_nodes[i].edges) {
                    if (forward_scores[i] + backward_scores[edge.first] + edge_weights[make_pair(i, edge.first)] >= min_path_score) {
                        keep_edges.emplace(i, edge.first);
                    }
                }
                removed_in_prefix[i + 1] = removed_in_prefix[i];
            }
            else {
                removed_in_prefix[i + 1] = removed_in_prefix[i] + 1;
            }
        }
        
        // prune down to these nodes and edges
        size_t next = 0;
        for (size_t i = 0; i < match_nodes.size(); i++) {
            if (keep_nodes.count(i)) {
                if (i != next) {
                    match_nodes[next] = std::move(match_nodes[i]);
                }
                vector<pair<size_t, size_t>>& edges = match_nodes[next].edges;
                
                size_t new_end = edges.size();
                for (size_t j = 0; j < new_end;) {
                    pair<size_t, size_t>& edge = edges[j];
                    if (!keep_edges.count(make_pair(i, edge.first))) {
                        new_end--;
                        edge = edges[new_end];
                    }
                    else {
                        edge.first -= removed_in_prefix[edge.first];
                        j++;
                    }
                }
                edges.resize(new_end);
                
                next++;
            }
        }
        match_nodes.resize(next);
    }
}




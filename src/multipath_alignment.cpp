//
//  multipath_alignment.cpp
//

#include "multipath_alignment.hpp"
#include <structures/immutable_list.hpp>
#include <structures/min_max_heap.hpp>

#include <type_traits>

//#define debug_multiple_tracebacks
//#define debug_verbose_validation

using namespace std;
using namespace structures;

namespace vg {
    
    void topologically_order_subpaths(MultipathAlignment& multipath_aln) {
        // Kahn's algorithm
        
        vector<size_t> index(multipath_aln.subpath_size(), 0);
        size_t order_idx = 0;
        
        vector<size_t> stack;
        vector<size_t> in_degree(multipath_aln.subpath_size(), 0);
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            for (size_t j = 0; j < subpath.next_size(); j++) {
                in_degree[subpath.next(j)]++;
            }
        }
        
        // identify the source nodes and add them to the stack
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            if (!in_degree[i]) {
                stack.push_back(i);
            }
        }
        
        while (!stack.empty()) {
            // pop a source node and add it to the topological order
            size_t here = stack.back();
            stack.pop_back();
            
            index[here] = order_idx;
            order_idx++;
            
            // remove the node's edges
            const Subpath& subpath = multipath_aln.subpath(here);
            for (size_t i = 0; i < subpath.next_size(); i++) {
                size_t next = subpath.next(i);
                in_degree[next]--;
                // if a node is now a source, stack it up
                if (!in_degree[next]) {
                    stack.push_back(next);
                }
            }
        }
        
        // translate the edges to the new indices
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            Subpath* subpath = multipath_aln.mutable_subpath(i);
            for (size_t j = 0; j < subpath->next_size(); j++) {
                subpath->set_next(j, index[subpath->next(j)]);
            }
        }
        
        // translate the start nodes
        for (size_t i = 0; i < multipath_aln.start_size(); i++) {
            multipath_aln.set_start(i, index[multipath_aln.start(i)]);
        }
        
        // in place permutation according to the topological order
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            while (index[i] != i) {
                std::swap(*multipath_aln.mutable_subpath(i), *multipath_aln.mutable_subpath(index[i]));
                std::swap(index[i], index[index[i]]);
            }
        }
    }
    
    void identify_start_subpaths(MultipathAlignment& multipath_aln) {
        
        // remove start nodes if there are any (avoids doubling them if this function is used liberally)
        multipath_aln.clear_start();
        
        // label nodes with incoming edges
        vector<bool> has_incoming_edge(multipath_aln.subpath_size(), false);
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            for (size_t j = 0; j < subpath.next_size(); j++) {
                has_incoming_edge[subpath.next(j)] = true;
            }
        }
        
        // construct list of starts
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            if (!has_incoming_edge[i]) {
                multipath_aln.add_start(i);
            }
        }
    }
    
    /// We define this struct for holding the dynamic programming problem for a
    /// multipath alignment, which we use for finding the optimal alignment,
    /// scoring the optimal alignment, and enumerating the top alignments.
    struct MultipathProblem {
        // Score of the optimal alignment ending immediately before this
        // subpath. To get the score of the optimal alignment ending with the
        // subpath, add the subpath's score.
        vector<int32_t> prefix_score;
        // previous subpath for traceback (we refer to subpaths by their index)
        vector<int64_t> prev_subpath;
        // the length of read sequence preceding this subpath
        vector<int64_t> prefix_length;
        
        /// Make a new MultipathProblem over the given number of subpaths with scores
        /// initialized according to whether we're doing a local or global traceback
        MultipathProblem(const MultipathAlignment& multipath_aln, bool subpath_global)
            : prefix_score(multipath_aln.subpath_size(), subpath_global ? numeric_limits<int32_t>::min() / 2 : 0),
              prev_subpath(multipath_aln.subpath_size(), -1), prefix_length(multipath_aln.subpath_size(), 0) {
            
            if (subpath_global) {
                // set the starting score at sources to 0 so that alignments can only start there
                for (const auto& i : multipath_aln.start()) {
                    prefix_score[i] = 0;
                }
            }
        }
    };
    
    /// Internal helper function for running the dynamic programming problem
    /// represented by a multipath alignment. Returns the filled DP problem,
    /// the optimal ending subpath, or -1 if no subpath is optimal, and the
    /// optimal score, or 0 if no score is optimal. An option toggles whether
    /// the traceback should be global (a source to a sink in the multipath DAG)
    /// or local (starting and ending at any subpath)
    tuple<MultipathProblem, int64_t, int32_t> run_multipath_dp(const MultipathAlignment& multipath_aln,
                                                               bool subpath_global = false) {
        
        // Create and unpack the return value (including setting up the DP table)
        tuple<MultipathProblem, int64_t, int32_t> to_return(MultipathProblem(multipath_aln, subpath_global),
                                                            -1, 0);
        auto& problem = get<0>(to_return);
        auto& opt_subpath = get<1>(to_return);
        auto& opt_score = get<2>(to_return);
    
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            int32_t extended_score = problem.prefix_score[i] + subpath.score();
            // carry DP forward
            if (subpath.next_size() > 0) {
                int64_t thru_length = path_to_length(subpath.path()) + problem.prefix_length[i];
                for (size_t j = 0; j < subpath.next_size(); j++) {
                    int64_t next = subpath.next(j);
                    problem.prefix_length[next] = thru_length;
                    
                    // can we improve prefix score on following subpath through this one?
                    if (extended_score >= problem.prefix_score[next]) {
                        problem.prev_subpath[next] = i;
                        problem.prefix_score[next] = extended_score;
                    }
                }
            }
            // check if an alignment is allowed to end here according to global/local rules and
            // if so whether it's optimal
            if (extended_score >= opt_score && (!subpath_global || subpath.next_size() == 0)) {
                // We have a better optimal subpath
                opt_score = extended_score;
                opt_subpath = i;
            }
        }
                
        return to_return;
    
    }
    
    /// We define this helper to turn tracebacks through a DP problem into
    /// Paths that we can put in an Alignment. We use iterators to the start
    /// and past-the-end of the traceback (in some kind of list of int64_t
    /// subpath indexes) to define it.
    template<typename TracebackIterator>
    void populate_path_from_traceback(const MultipathAlignment& multipath_aln, const MultipathProblem& problem,
        TracebackIterator traceback_start, TracebackIterator traceback_end, Path* output) {
        
        static_assert(is_convertible<decltype(*traceback_start), int64_t>::value, "traceback must contain int64_t items");
        
        if (traceback_start == traceback_end) {
            // We have been given an empty range. Do nothing.
            return;
        }
        
        // We need to maintain a persistent iterator so we can easily get whatever's before traceback_end.
        auto current_subpath = traceback_start;
    
        // check for a softclip of entire subpaths on the beginning
        if (problem.prefix_length[*current_subpath]) {
            Mapping* soft_clip_mapping = output->add_mapping();
            
            soft_clip_mapping->set_rank(1);
            
            Edit* edit = soft_clip_mapping->add_edit();
            edit->set_to_length(problem.prefix_length[*current_subpath]);
            edit->set_sequence(multipath_aln.sequence().substr(0, problem.prefix_length[*current_subpath]));
            
            *soft_clip_mapping->mutable_position() = multipath_aln.subpath(*current_subpath).path().mapping(0).position();
        }
        
        // merge the subpaths into one optimal path in the Alignment object
        for (auto next_subpath = current_subpath; next_subpath != traceback_end; ++next_subpath) {
            // For all subpaths in the traceback
            
            // Advance only if we don't hit the end
            current_subpath = next_subpath;
            
            if (output->mapping_size() == 0) {
                // There's nothing in the output yet, so just copy all the mappings from this subpath
                *output = multipath_aln.subpath(*current_subpath).path();
                
                for(size_t i = 0; i < output->mapping_size(); i++) {
                    // Set all the ranks
                    output->mutable_mapping(i)->set_rank(i + 1);
                }
            } else {
                // There's already content in the output so we have to merge stuff
                Mapping* curr_end_mapping = output->mutable_mapping(output->mapping_size() - 1);
            
                // get the first mapping of the next path
                const Path& next_path = multipath_aln.subpath(*current_subpath).path();
                const Mapping& next_start_mapping = next_path.mapping(0);
                
                size_t mapping_start_idx = 0;
                // merge mappings if they occur on the same node and same strand
                if (curr_end_mapping->position().node_id() == next_start_mapping.position().node_id()
                    && curr_end_mapping->position().is_reverse() == next_start_mapping.position().is_reverse()) {
                    
                    Edit* last_edit = curr_end_mapping->mutable_edit(curr_end_mapping->edit_size() - 1);
                    const Edit& first_edit = next_start_mapping.edit(0);
                    
                    // merge the first edit if it is the same type
                    size_t edit_start_idx = 0;
                    if ((last_edit->from_length() > 0) == (first_edit.from_length() > 0)
                            && (last_edit->to_length() > 0) == (first_edit.to_length() > 0)
                            && (last_edit->sequence().empty()) == (first_edit.sequence().empty())) {
                            
                            last_edit->set_from_length(last_edit->from_length() + first_edit.from_length());
                            last_edit->set_to_length(last_edit->to_length() + first_edit.to_length());
                            last_edit->set_sequence(last_edit->sequence() + first_edit.sequence());
                        
                        edit_start_idx++;
                    }
                        
                    // append the rest of the edits
                    for (size_t j = edit_start_idx; j < next_start_mapping.edit_size(); j++) {
                        *curr_end_mapping->add_edit() = next_start_mapping.edit(j);
                    }
                    
                    mapping_start_idx++;
                }
                
                // append the rest of the mappings
                for (size_t j = mapping_start_idx; j < next_path.mapping_size(); j++) {
                    Mapping* next_mapping = output->add_mapping();
                    *next_mapping = next_path.mapping(j);
                    next_mapping->set_rank(output->mapping_size());
                }
            }
        }
        
        // Now current_subpath is right before traceback_end
        
        // check for a softclip of entire subpaths on the end
        int64_t seq_thru_length = problem.prefix_length[*current_subpath] + path_to_length(multipath_aln.subpath(*current_subpath).path());
        if (seq_thru_length < multipath_aln.sequence().size()) {
            
            if (output->mapping_size() == 0) {
                output->add_mapping();
            }
            
            Mapping* final_mapping = output->mutable_mapping(output->mapping_size() - 1);
            
            Edit* edit = final_mapping->add_edit();
            edit->set_to_length(multipath_aln.sequence().size() - seq_thru_length);
            edit->set_sequence(multipath_aln.sequence().substr(seq_thru_length,
                                                               multipath_aln.sequence().size() - seq_thru_length));
        }
    }
    
    int32_t optimal_alignment_internal(const MultipathAlignment& multipath_aln, Alignment* aln_out,
                                       bool subpath_global) {
        
        // Run the dynamic programming
        auto dp_result = run_multipath_dp(multipath_aln, subpath_global);
        
        // C++17 finally gets http://en.cppreference.com/w/cpp/language/structured_binding
        // Until then we have to unpack tuples like this.
        
        // Get the filled DP problem
        MultipathProblem& problem = get<0>(dp_result);
        // And the optimal final subpath
        int64_t& opt_subpath = get<1>(dp_result);
        // And the optimal score
        int32_t& opt_score = get<2>(dp_result);
        
        // are we constructing the alignment, or just getting the score?
        if (aln_out && opt_subpath >= 0) {
            
            // traceback the optimal subpaths until hitting sentinel (-1)
            list<int64_t> opt_traceback;
            int64_t curr = opt_subpath;
            while (curr >= 0) {
                opt_traceback.push_front(curr);
                curr = problem.prev_subpath[curr];
            }
            
            Path* opt_path = aln_out->mutable_path();
            
            // Fill in the path in the alignment with the alignment represented
            // by this traceback in this DP problem for this multipath
            // alignment.
            populate_path_from_traceback(multipath_aln, problem, opt_traceback.begin(), opt_traceback.end(), opt_path);
            
            
        }
        
        // Return the optimal score, or 0 if unaligned.
        return opt_score;
    }
    
    void optimal_alignment(const MultipathAlignment& multipath_aln, Alignment& aln_out, bool subpath_global) {

        // transfer read information over to alignment
        transfer_read_metadata(multipath_aln, aln_out);
        aln_out.set_mapping_quality(multipath_aln.mapping_quality());
        
        // do dynamic programming and traceback the optimal alignment
        int32_t score = optimal_alignment_internal(multipath_aln, &aln_out, subpath_global);
        
        aln_out.set_score(score);
    }
    
    int32_t optimal_alignment_score(const MultipathAlignment& multipath_aln, bool subpath_global){
        // do dynamic programming without traceback
        return optimal_alignment_internal(multipath_aln, nullptr, subpath_global);
    }
    
    vector<Alignment> optimal_alignments(const MultipathAlignment& multipath_aln, size_t count) {
        
#ifdef debug_multiple_tracebacks
        cerr << "Computing top " << count << " alignments" << endl;
#endif
        
        // Keep a list of what we're going to emit.
        vector<Alignment> to_return;
        
        // Fill out the dynamic programming problem
        auto dp_result = run_multipath_dp(multipath_aln);
        // Get the filled DP problem
        MultipathProblem& problem = get<0>(dp_result);
        // And the optimal final subpath
        int64_t& opt_subpath = get<1>(dp_result);
        // And the optimal score
        int32_t& opt_score = get<2>(dp_result);
        
        // Keep lists of DP steps, which are subpath numbers to visit.
        // Even going to the end subpath (where prefix length + subpath length = read length) is a DP step
        // We never deal with empty lists; we always seed with the traceback start node.
        using step_list_t = ImmutableList<int64_t>;
        
        // Put them in a size-limited priority queue by score difference (positive) from optimal
        MinMaxHeap<pair<int32_t, step_list_t>> queue;
        
        // We define a function to put stuff in the queue and limit its size to
        // the (count - to_return.size()) items with lowest penalty.
        auto try_enqueue = [&](const pair<int32_t, step_list_t>& item) {
            auto max_size = count - to_return.size();
            if (queue.size() < max_size || item < queue.max()) {
                // The item belongs in the queue because it fits or it beats
                // the current worst thing.
                queue.push(item);
            } else {
#ifdef debug_multiple_tracebacks
                cerr << "Rejected! Queue is full!" << endl;
#endif
            }
            
            while(queue.size() > max_size) {
                // We have more possibilities than we need to consider to emit
                // the top count alignments. Get rid of the worst one.
                queue.pop_max();
#ifdef debug_multiple_tracebacks
                cerr << "Existing item displaced!" << endl;
#endif
            }
        };
        
        // Also, subpaths only keep track of their nexts, so we need to invert
        // that so we can get all valid prev subpaths.
        vector<vector<int64_t>> prev_subpaths;
        
        // We want to be able to start the traceback only from places where we
        // won't get shorter versions of same- or higher-scoring alignments.
        // This means that we want exactly the subpaths that have no successors
        // with nonnegative subpath score.
        
        // So we go through all the subpaths, check all their successors, and
        // add starting points for all satisfactory subpaths to the queue.
        // Sinks in the graphs will always be in here, as will the starting
        // point for tracing back the optimal alignment.
        
        // We know what the penalty from optimal is for each, because we know
        // the optimal score overall and the score we would get for the optimal
        // alignment ending at each.
        
        prev_subpaths.resize(multipath_aln.subpath_size());
        for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
            // For each subpath
            
            // If it has no successors, we can start a traceback here
            bool valid_traceback_start = true;
            
            for (auto& next_subpath : multipath_aln.subpath(i).next()) {
                // For each next subpath it lists
                
                // Register this subpath as a predecessor of the next
                prev_subpaths[next_subpath].push_back(i);
                
                if (multipath_aln.subpath(next_subpath).score() >= 0) {
                    // This successor has a nonnegative score, so taking it
                    // after us would generate a longer, same- or
                    // higher-scoring alignment. So we shouldn't start a
                    // traceback from subpath i.
                    valid_traceback_start = false;
                }
            }
            
            if (valid_traceback_start) {
                // We can start a traceback here.
                
                // The score penalty for starting here is the optimal score minus the optimal score starting here
                auto penalty = opt_score - (problem.prefix_score[i] + multipath_aln.subpath(i).score());
                
                // The path is just to be here
                step_list_t starting_path{i};
                
#ifdef debug_multiple_tracebacks
                cerr << "Could end at subpath " << i << " with penalty " << penalty << endl;
#endif
                
                try_enqueue(make_pair(penalty, starting_path));
            }
        }
        
        while (!queue.empty() && to_return.size() < count) {
            // Each iteration
            
            // Grab the best list as our basis
            int32_t basis_score_difference;
            step_list_t basis;
            tie(basis_score_difference, basis) = queue.min();
            queue.pop_min();
            
            assert(!basis.empty());
            
#ifdef debug_multiple_tracebacks
            size_t basis_size = 0;
            for (auto& i : basis) {
                basis_size++;
            }
            cerr << "Consider " << basis_size << " element traceback to " << basis.front() << " with penalty "
                 << basis_score_difference << endl;
            cerr << "\t" << pb2json(multipath_aln.subpath(basis.front()).path()) << endl;
#endif
            
            if (problem.prev_subpath[basis.front()] == -1) {
                // If it leads all the way to a subpath that is optimal as a start
                
                // Make an Alignment to emit it in
                to_return.emplace_back();
                Alignment& aln_out = to_return.back();
                
                // Set up read info and MAPQ
                // TODO: MAPQ on secondaries?
                transfer_read_metadata(multipath_aln, aln_out);
                aln_out.set_mapping_quality(multipath_aln.mapping_quality());
                
                // Populate path
                populate_path_from_traceback(multipath_aln, problem, basis.begin(), basis.end(), aln_out.mutable_path());
                
                // Set score
                aln_out.set_score(opt_score - basis_score_difference);
                
#ifdef debug_multiple_tracebacks
                cerr << "Traceback reaches start; emit with score " << aln_out.score() << endl;
#endif
                
            } else {
                // The path does not lead all the way to a source
                
                // Find out all the places we can come from, and the score
                // penalties, relative to the optimal score for an alignment
                // visiting the basis's lead subpath, that we would take if we
                // came from each.
                // Note that we will only do this once per subpath, when we are
                // working on the optimal alignment going through that subpath.
                list<pair<int64_t, int32_t>> destinations;
                
                // The destinations will be all places we could have arrived here from
                auto& here = basis.front();
                
                // To compute the additional score difference, we need to know what our optimal prefix score was.
                auto& best_prefix_score = problem.prefix_score[here];
                
                for (auto& prev : prev_subpaths[here]) {
                    // For each, compute the score of the optimal alignment ending at that predecessor
                    auto prev_opt_score = problem.prefix_score[prev] + multipath_aln.subpath(prev).score();
                    
                    // What's the difference we would take if we went with this predecessor?
                    auto additional_penalty = best_prefix_score - prev_opt_score;
                    
                    destinations.emplace_back(prev, additional_penalty);
                }
                
                // TODO: unify loops!
                
                for (auto& destination : destinations) {
                    // Prepend each of the things that can be prepended
                    
                    // Unpack
                    auto& prev = destination.first;
                    auto& additional_penalty = destination.second;
                    
                    // Make an extended path
                    auto extended_path = basis.push_front(prev);
                    
                    // Calculate the score differences from optimal
                    auto total_penalty = basis_score_difference + additional_penalty;
                    
#ifdef debug_multiple_tracebacks
                    cerr << "\tAugment with " << prev << " to penalty " << total_penalty << endl;
#endif
                    
                    // Put them in the priority queue
                    try_enqueue(make_pair(total_penalty, extended_path));
                }
                
            }
        }
        
        return to_return;
        
    }
    
    vector<Alignment> optimal_alignments_with_disjoint_subpaths(const MultipathAlignment& multipath_aln, size_t count) {
        
#ifdef debug_multiple_tracebacks
        cerr << "Computing top " << count << " alignments with disjoint subpaths" << endl;
#endif
        
        // Keep a list of what we're going to emit.
        vector<Alignment> to_return;
        
        // Fill out the dynamic programming problem
        auto dp_result = run_multipath_dp(multipath_aln);
        // Get the filled DP problem
        MultipathProblem& problem = get<0>(dp_result);
        // And the optimal final subpath
        int64_t& opt_subpath = get<1>(dp_result);
        // And the optimal score
        int32_t& opt_score = get<2>(dp_result);
        
        // Keep lists of DP steps
        using step_list_t = ImmutableList<int64_t>;
        
        // Have a queue just for end positions
        MinMaxHeap<pair<int32_t, step_list_t>> end_queue;
        
        // Also, subpaths only keep track of their nexts, so we need to invert
        // that so we can get all valid prev subpaths.
        vector<vector<int64_t>> prev_subpaths;
        
        prev_subpaths.resize(multipath_aln.subpath_size());
        for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
            // For each subpath
            
            // If it has no successors, we can start a traceback here
            bool valid_traceback_start = true;
            
            for (auto& next_subpath : multipath_aln.subpath(i).next()) {
                // For each next subpath it lists
                
                // Register this subpath as a predecessor of the next
                prev_subpaths[next_subpath].push_back(i);
                
                if (multipath_aln.subpath(next_subpath).score() >= 0) {
                    // This successor has a nonnegative score, so taking it
                    // after us would generate a longer, same- or
                    // higher-scoring alignment. So we shouldn't start a
                    // traceback from subpath i.
                    valid_traceback_start = false;
                }
            }
            
            if (valid_traceback_start) {
                // We can start a traceback here.
                
                // The score penalty for starting here is the optimal score minus the optimal score starting here
                auto penalty = opt_score - (problem.prefix_score[i] + multipath_aln.subpath(i).score());
                
                // The path is just to be here
                step_list_t starting_path{i};
                
#ifdef debug_multiple_tracebacks
                cerr << "Could end at subpath " << i << " with penalty " << penalty << endl;
#endif
                
                end_queue.push(make_pair(penalty, starting_path));
            }
        }
        
        // Keep a bit vector of the subpaths that have been used, so we can reject
        // them. TODO: We get the optimal alignment for each end, subject to
        // the constraint, but any other subpath may be used in a suboptimal
        // alignment for that subpath, and we may never see its optimal
        // alignment.
        vector<bool> subpath_is_used(multipath_aln.subpath_size(), false);
        
        while (!end_queue.empty() && to_return.size() < count) {
            // For each distinct ending subpath in the multipath
            
#ifdef debug_multiple_tracebacks
            cerr << "Look for alignment " << to_return.size() << " ending with " << end_queue.min().second.front() << endl;
#endif
            
            // Make a real queue for starting from it
            MinMaxHeap<pair<int32_t, step_list_t>> queue;
            queue.push(end_queue.min());
            end_queue.pop_min();
            
            if (subpath_is_used[queue.min().second.front()]) {
                // We shouldn't ever have the place we want to trace back from already used, but if it is already used we don't want to use it.
                
#ifdef debug_multiple_tracebacks
                cerr << "Skip " << queue.min().second.front() << " because it was already emitted in a previous traceback" << endl;
#endif
                
                continue;
            }
            
            // We also want to remember the lowest penalty with which each
            // subpath has been queued, so we can do a real Dijkstra traversal
            // and not waste all our time on combinatorial paths to get places
            // with the same or higher penalty.
            vector<size_t> min_penalty_for_subpath(multipath_aln.subpath_size(), numeric_limits<size_t>::max());
            // Seed with the end we are starting with.
            min_penalty_for_subpath[queue.min().second.front()] = queue.min().first;
            
            // We also track visited-ness, so we don;t query edges for the same thing twice.
            // TODO: This is the world's most hacky Dijkstra and needs to be rewritten from the top with an understanding of what it is supposed to be doing.
            vector<bool> subpath_is_visited(multipath_aln.subpath_size(), false);
        
            while (!queue.empty() && to_return.size() < count) {
                // Each iteration
                
                // Grab the best list as our basis
                int32_t basis_score_difference;
                step_list_t basis;
                tie(basis_score_difference, basis) = queue.min();
                queue.pop_min();
                
                assert(!basis.empty());
                
#ifdef debug_multiple_tracebacks
                size_t basis_size = 0;
                for (auto& i : basis) {
                    basis_size++;
                }
                cerr << "Consider " << basis_size << " element traceback to " << basis.front() << " with penalty "
                     << basis_score_difference << endl;
                cerr << "\t" << pb2json(multipath_aln.subpath(basis.front()).path()) << endl;
#endif

                if (subpath_is_used[basis.front()]) {
                    // We already used this start and can't use it again. Try something else.
                    // TODO: This shouldn't happen; we are also catching this case on enqueue.
#ifdef debug_multiple_tracebacks
                    cerr << "Traceback reaches already used subpath; skip" << endl;
#endif
                    
                    continue;
                }
                
                if (subpath_is_visited[basis.front()]) {
                    // We already processed this; this must be a higher cost version of the same thing.
                    assert(basis_score_difference >= min_penalty_for_subpath[basis.front()]);
                    
#ifdef debug_multiple_tracebacks
                    cerr << "Found more expensive version of subpath that has already been processed; skip" << endl;
#endif
                    
                    continue;
                }
                subpath_is_visited[basis.front()] = true;
                
                if (problem.prev_subpath[basis.front()] == -1) {
                    // If it leads all the way to a subpath that is optimal as a start
                    
                    // Make an Alignment to emit it in
                    to_return.emplace_back();
                    Alignment& aln_out = to_return.back();
                    
                    // Set up read info and MAPQ
                    // TODO: MAPQ on secondaries?
                    transfer_read_metadata(multipath_aln, aln_out);
                    aln_out.set_mapping_quality(multipath_aln.mapping_quality());
                    
                    // Populate path
                    populate_path_from_traceback(multipath_aln, problem, basis.begin(), basis.end(), aln_out.mutable_path());
                    
                    // Set score
                    aln_out.set_score(opt_score - basis_score_difference);
                    
#ifdef debug_multiple_tracebacks
                    cerr << "Traceback reaches start; emit with score " << aln_out.score() << endl;
#endif

                    for (auto& subpath : basis) {
                        // Record the used-ness of all the subpaths
                        subpath_is_used[subpath] = true;
                    }

                    // Break out of this loop and try a different starting position
                    break;
                    
                } else {
                    // The path does not lead all the way to a source
                    
                    // The destinations will be all places we could have arrived here from
                    auto& here = basis.front();
                    
                    // To compute the additional score difference, we need to know what our optimal prefix score was.
                    auto& best_prefix_score = problem.prefix_score[here];
                    
                    for (auto& prev : prev_subpaths[here]) {
                        // For each candidate previous subpath
                        
                        if (subpath_is_used[prev]) {
                            // This subpath has already been used in an emitted alignment, so we can't use it.
                            
#ifdef debug_multiple_tracebacks
                            cerr << "\tSkip " << prev << " which is already used" << endl;
#endif
                            
                            continue;
                        }
                        
                        // For each, compute the score of the optimal alignment ending at that predecessor
                        auto prev_opt_score = problem.prefix_score[prev] + multipath_aln.subpath(prev).score();
                        
                        // What's the difference we would take if we went with this predecessor?
                        auto additional_penalty = best_prefix_score - prev_opt_score;
                        
                        // Calculate the score differences from optimal
                        auto total_penalty = basis_score_difference + additional_penalty;
                        
                        if (total_penalty >= min_penalty_for_subpath[prev]) {
                            // This previous subpath is already reachable with a penalty as good or better.
                            // Don't bother with it again
                            
#ifdef debug_multiple_tracebacks
                            cerr << "\tSkip " << prev << " with penalty " << total_penalty << " >= " << min_penalty_for_subpath[prev] << endl;
#endif
                            
                            continue;
                        }
                        
                        // Record that this is the cheapest we managed to get here
                        min_penalty_for_subpath[prev] = total_penalty;
                        
                        // Make an extended path
                        auto extended_path = basis.push_front(prev);
                        
#ifdef debug_multiple_tracebacks
                        cerr << "\tAugment with " << prev << " to penalty " << total_penalty << endl;
#endif
                        
                        // Put them in the priority queue
                        queue.push(make_pair(total_penalty, extended_path));
                    }
                    
                    
                }
            }
        }
        
        return to_return;
        
    }
    
    /// Stores the reverse complement of a Subpath in another Subpath
    ///
    /// note: this is not included in the header because reversing a subpath without going through
    /// the multipath alignment can break invariants related to the edge lists
    ///
    ///  Args:
    ///    subpath           subpath to reverse complement
    ///    node_length       a function that returns the length of a node sequence from its node ID
    ///    rev_comp_out      empty subpath to store reverse complement in (data will be overwritten
    ///                      if not empty)
    ///
    inline void rev_comp_subpath(const Subpath& subpath, const function<int64_t(int64_t)>& node_length,
                          Subpath& rev_comp_out) {
        
        *(rev_comp_out.mutable_path()) = reverse_complement_path(subpath.path(), node_length);
        rev_comp_out.set_score(subpath.score());
        // leave reversing the edges to the multipath alignment
    }
    
    void rev_comp_multipath_alignment(const MultipathAlignment& multipath_aln, const function<int64_t(int64_t)>& node_length,
                                      MultipathAlignment& rev_comp_out) {
        
        // reverse complement sequence
        rev_comp_out.set_sequence(reverse_complement(multipath_aln.sequence()));
        // reverse base qualities
        rev_comp_out.set_quality(string(multipath_aln.quality().rbegin(), multipath_aln.quality().rend()));
        
        // transfer the rest of the metadata directly
        rev_comp_out.set_read_group(multipath_aln.read_group());
        rev_comp_out.set_name(multipath_aln.name());
        rev_comp_out.set_sample_name(multipath_aln.sample_name());
        rev_comp_out.set_paired_read_name(multipath_aln.paired_read_name());
        
        vector< vector<size_t> > reverse_edge_lists(multipath_aln.subpath_size());
        vector<size_t> reverse_starts;
        
        // remove subpaths to avoid duplicating
        rev_comp_out.clear_subpath();
        
        // add subpaths in reverse order to maintain topological ordering
        for (int64_t i = multipath_aln.subpath_size() - 1; i >= 0; i--) {
            const Subpath& subpath = multipath_aln.subpath(i);
            Subpath* rc_subpath = rev_comp_out.add_subpath();
            rev_comp_subpath(subpath, node_length, *rc_subpath);
            
            if (subpath.next_size() > 0) {
                // collect edges by their target (for reversing)
                for (size_t j = 0; j < subpath.next_size(); j++) {
                    reverse_edge_lists[subpath.next(j)].push_back(i);
                }
            }
            else {
                // sink subpaths become sources in reverse
                reverse_starts.push_back(i);
            }
        }
        
        // add reversed edges
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            Subpath* rc_subpath = rev_comp_out.mutable_subpath(i);
            vector<size_t>& reverse_edge_list = reverse_edge_lists[multipath_aln.subpath_size() - i - 1];
            for (size_t j = 0; j < reverse_edge_list.size(); j++) {
                rc_subpath->add_next(multipath_aln.subpath_size() - reverse_edge_list[j] - 1);
            }
        }
        
        // remove start nodes that are invalid in reverse
        rev_comp_out.clear_start();
        
        // assume that if the original multipath alignment had its starts labeled they want them
        // labeled in the reverse complement too
        if (multipath_aln.start_size() > 0) {
            for (size_t i = 0; i < reverse_starts.size(); i++) {
                rev_comp_out.add_start(multipath_aln.subpath_size() - reverse_starts[i] - 1);
            }
        }
    }
    
    void rev_comp_multipath_alignment_in_place(MultipathAlignment* multipath_aln,
                                               const function<int64_t(int64_t)>& node_length) {
        
        // reverse complement sequence
        reverse_complement_in_place(*multipath_aln->mutable_sequence());
        // reverse base qualities
        string* quality = multipath_aln->mutable_quality();
        std::reverse(quality->begin(), quality->end());
        
        // current reverse edges
        vector< vector<int64_t> > reverse_edge_lists(multipath_aln->subpath_size());
        // current sink nodes (will be starts)
        vector<int64_t> reverse_starts;
        
        int64_t subpath_swap_size = multipath_aln->subpath_size() / 2;
        int64_t last = multipath_aln->subpath_size() - 1;
        for (int64_t i = 0, j = last; i < subpath_swap_size; i++, j--) {
            Subpath* subpath_1 = multipath_aln->mutable_subpath(i);
            Subpath* subpath_2 = multipath_aln->mutable_subpath(j);
            
            // add reverse edges for first subpath
            if (subpath_1->next_size() > 0) {
                for (int64_t k = 0; k < subpath_1->next_size(); k++) {
                    reverse_edge_lists[subpath_1->next(k)].push_back(i);
                }
            }
            else {
                reverse_starts.push_back(i);
            }
            
            // add reverse edges for second subpath
            if (subpath_2->next_size() > 0) {
                for (int64_t k = 0; k < subpath_2->next_size(); k++) {
                    reverse_edge_lists[subpath_2->next(k)].push_back(j);
                }
            }
            else {
                reverse_starts.push_back(j);
            }
            
            // clear current edges
            subpath_1->clear_next();
            subpath_2->clear_next();
            
            // reverse complement the paths
            reverse_complement_path_in_place(subpath_1->mutable_path(), node_length);
            reverse_complement_path_in_place(subpath_2->mutable_path(), node_length);
            
            // swap their positions (to maintain topological ordering)
            std::swap(*subpath_1, *subpath_2);
        }
        
        // repeat process for the middle subpath if there is an odd number
        if (multipath_aln->subpath_size() % 2) {
            Subpath* subpath = multipath_aln->mutable_subpath(subpath_swap_size);
            if (subpath->next_size() > 0) {
                for (int64_t k = 0; k < subpath->next_size(); k++) {
                    reverse_edge_lists[subpath->next(k)].push_back(subpath_swap_size);
                }
            }
            else {
                reverse_starts.push_back(subpath_swap_size);
            }
            
            subpath->clear_next();
            reverse_complement_path_in_place(subpath->mutable_path(), node_length);
        }
        
        // add reversed edges
        for (int64_t i = 0, j = last; i < multipath_aln->subpath_size(); i++, j--) {
            vector<int64_t> edges = reverse_edge_lists[j];
            Subpath* subpath = multipath_aln->mutable_subpath(i);
            for (int64_t k : edges) {
                subpath->add_next(last - k);
            }
        }
        
        // if we had starts labeled before, label them again
        if (multipath_aln->start_size() > 0) {
            multipath_aln->clear_start();
            for (int64_t i : reverse_starts) {
                multipath_aln->add_start(last - i);
            }
        }
    }
    
    void to_multipath_alignment(const Alignment& aln, MultipathAlignment& multipath_aln_out) {
        
        // clear repeated fields
        multipath_aln_out.clear_subpath();
        multipath_aln_out.clear_start();
        
        // transfer read and alignment metadata
        transfer_read_metadata(aln, multipath_aln_out);
        multipath_aln_out.set_mapping_quality(aln.mapping_quality());
        
        // transfer alignment and score
        if (aln.has_path() || aln.score()) {
            Subpath* subpath = multipath_aln_out.add_subpath();
            subpath->set_score(aln.score());
            *(subpath->mutable_path()) = aln.path();
        }
        
    }
    
    void transfer_read_metadata(const MultipathAlignment& from, MultipathAlignment& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        to.set_read_group(from.read_group());
        to.set_name(from.name());
        to.set_sample_name(from.sample_name());
        to.set_paired_read_name(from.paired_read_name());
    }
    
    void transfer_read_metadata(const Alignment& from, MultipathAlignment& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        to.set_read_group(from.read_group());
        to.set_name(from.name());
        to.set_sample_name(from.sample_name());
        
        // no difference in these fields for MultipathAlignments
        if (from.has_fragment_prev()) {
            to.set_paired_read_name(from.fragment_prev().name());
        }
        else if (from.has_fragment_next()) {
            to.set_paired_read_name(from.fragment_next().name());
        }
    }
    
    void transfer_read_metadata(const MultipathAlignment& from, Alignment& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        to.set_read_group(from.read_group());
        to.set_name(from.name());
        to.set_sample_name(from.sample_name());
        
        // note: not transferring paired_read_name because it is unclear whether
        // it should go into fragment_prev or fragment_next
    }
    
    void merge_non_branching_subpaths(MultipathAlignment& multipath_aln) {
        
        vector<size_t> in_degree(multipath_aln.subpath_size(), 0);
        for (const Subpath& subpath : multipath_aln.subpath()) {
            for (int64_t next : subpath.next()) {
                in_degree[next]++;
            }
        }
        
        vector<bool> removed(multipath_aln.subpath_size(), false);
        vector<size_t> removed_so_far(multipath_aln.subpath_size(), 0);
        
        auto get_mergeable_next = [&](const Subpath& subpath) {
            if (subpath.next_size() == 1) {
                if (in_degree[subpath.next(0)] == 1) {
                    return int64_t(subpath.next(0));
                }
            }
            return int64_t(-1);
        };
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            
            if (i > 0) {
                removed_so_far[i] = removed_so_far[i - 1];
            }
            
            // this one has been marked for removal,
            if (removed[i]) {
                removed_so_far[i]++;
                continue;
            }
            
            // the subpath we might merge into
            Subpath* subpath = multipath_aln.mutable_subpath(i);
            
            // move it up in the vector if we've removed earlier subpaths
            if (removed_so_far[i] > 0) {
                *multipath_aln.mutable_subpath(i - removed_so_far[i]) = move(*subpath);
                subpath = multipath_aln.mutable_subpath(i - removed_so_far[i]);
            }
            
            int64_t last = -1;
            // iterate through non-branching subpaths
            for (int64_t j = get_mergeable_next(*subpath); j >= 0; j = get_mergeable_next(multipath_aln.subpath(j))) {
                
                // mark the next one for removal
                removed[j] = true;
                
                const Subpath& merge_subpath = multipath_aln.subpath(j);
                
                subpath->set_score(subpath->score() + merge_subpath.score());
                
                const Path& merge_path = merge_subpath.path();
                if (merge_path.mapping_size() == 0) {
                    continue;
                }
                
                Path* path = subpath->mutable_path();
                Mapping* final_mapping = path->mutable_mapping(path->mapping_size() - 1);
                const Position& final_position = final_mapping->position();
                
                const Mapping& first_mapping = merge_path.mapping(0);
                const Position& first_position = first_mapping.position();
                
                int64_t mapping_idx = 0;
                
                // do we need to merge the abutting mappings?
                if (first_position.node_id() == final_position.node_id() &&
                    first_position.is_reverse() == final_position.is_reverse() &&
                    first_position.offset() == final_position.offset() + mapping_from_length(*final_mapping)) {
                    
                    // do we need to merge the abutting edits?
                    int64_t edit_idx = 0;
                    if (final_mapping->edit_size() && first_mapping.edit_size()) {
                        Edit* final_edit = final_mapping->mutable_edit(0);
                        const Edit& first_edit = first_mapping.edit(0);
                        if ((first_edit.from_length() > 0) == (final_edit->from_length() > 0) &&
                            (first_edit.to_length() > 0) == (final_edit->to_length() > 0) &&
                            first_edit.sequence().empty() == final_edit->sequence().empty()) {
                            
                            final_edit->set_from_length(final_edit->from_length() + first_edit.from_length());
                            final_edit->set_to_length(final_edit->to_length() + first_edit.to_length());
                            final_edit->set_sequence(final_edit->sequence() + first_edit.sequence());
                            
                            edit_idx++;
                        }
                    }
                    
                    // append rest of the edits
                    for (; edit_idx < first_mapping.edit_size(); edit_idx++) {
                        *final_mapping->add_edit() = first_mapping.edit(edit_idx);
                    }
                    
                    mapping_idx++;
                }
                
                // append rest of the mappings
                for (; mapping_idx < merge_path.mapping_size(); mapping_idx++) {
                    *path->add_mapping() = merge_path.mapping(mapping_idx);
                }
                
                last = j;
            }
            
            // move the adjacencies over from the last one we merged in
            if (last >= 0) {
                subpath->clear_next();
                for (int64_t next : multipath_aln.subpath(last).next()) {
                    subpath->add_next(next);
                }
            }
        }
        
        // did we merge and remove any subpaths?
        if (removed_so_far.back()) {
            // trim the vector of subpaths
            multipath_aln.mutable_subpath()->DeleteSubrange(multipath_aln.subpath_size() - removed_so_far.back(),
                                                            removed_so_far.back());
            
            // update the indexes of the adjacencies
            for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
                Subpath* subpath = multipath_aln.mutable_subpath(i);
                for (size_t j = 0; j < subpath->next_size(); j++) {
                    subpath->set_next(j, subpath->next(j) - removed_so_far[subpath->next(j)]);
                }
            }
        }
    }
    
    vector<vector<int64_t>> connected_components(const MultipathAlignment& multipath_aln) {
        
        int64_t comps = 0;
        
        vector<vector<int64_t>> reverse_edge_lists(multipath_aln.subpath_size());
        
        for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            // collect edges by their target
            for (size_t j = 0; j < subpath.next_size(); j++) {
                reverse_edge_lists[subpath.next(j)].push_back(i);
            }
        }
        
        vector<bool> collected(multipath_aln.subpath_size(), false);
        
        vector<vector<int64_t>> components;
        
        for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
            if (collected[i]) {
                continue;
            }
            
            components.emplace_back();
            
            vector<int64_t> stack{i};
            collected[i] = true;
            while (!stack.empty()) {
                int64_t at = stack.back();
                stack.pop_back();
                
                components.back().push_back(at);
                
                const Subpath& subpath = multipath_aln.subpath(at);
                for (int64_t j = 0; j < subpath.next_size(); j++) {
                    int64_t idx = subpath.next(j);
                    if (!collected[idx]) {
                        collected[idx] = true;
                        stack.push_back(idx);
                    }
                }
                for (int64_t idx : reverse_edge_lists[at]) {
                    if (!collected[idx]) {
                        collected[idx] = true;
                        stack.push_back(idx);
                    }
                }
            }
        }
        
        return std::move(components);
    }
    
    void extract_sub_multipath_alignment(const MultipathAlignment& multipath_aln,
                                         const vector<int64_t>& subpath_indexes,
                                         MultipathAlignment& sub_multipath_aln) {
        sub_multipath_aln.Clear();
        transfer_read_metadata(multipath_aln, sub_multipath_aln);
        
        // create subpaths for each of the ones we're retaining and record the translation
        unordered_map<int64_t, int64_t> new_index;
        for (int64_t i = 0; i < subpath_indexes.size(); i++) {
            int64_t old_idx = subpath_indexes[i];
            const Subpath& old_subpath = multipath_aln.subpath(old_idx);
            
            Subpath* subpath = sub_multipath_aln.add_subpath();
            *subpath->mutable_path() = old_subpath.path();
            subpath->set_score(old_subpath.score());
            
            new_index[old_idx] = i;
        }
        
        // add edges according to the translation
        for (int64_t i = 0; i < subpath_indexes.size(); i++) {
            const Subpath& old_subpath = multipath_aln.subpath(subpath_indexes[i]);
            Subpath* new_subpath = sub_multipath_aln.mutable_subpath(i);
            for (int64_t j = 0; j < old_subpath.next_size(); j++) {
                if (new_index.count(old_subpath.next(j))) {
                    new_subpath->add_next(new_index[old_subpath.next(j)]);
                }
            }
        }
        
        // assume that if we had starts labeled before, we want them again
        if (multipath_aln.start_size() > 0) {
            identify_start_subpaths(sub_multipath_aln);
        }
    }
    
    
    bool validate_multipath_alignment(const MultipathAlignment& multipath_aln, const HandleGraph& handle_graph) {
        
        // are the subpaths in topological order?
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            for (size_t j = 0; j < subpath.next_size(); j++) {
                if (subpath.next(j) <= i) {
#ifdef debug_verbose_validation
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
#ifdef debug_verbose_validation
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
#ifdef debug_verbose_validation
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
#ifdef debug_verbose_validation
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
#ifdef debug_verbose_validation
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
#ifdef debug_verbose_validation
                cerr << "validation failure on containing only nonempty paths" << endl;
                cerr << "subpath " << i << ": " << pb2json(multipath_aln.subpath(i)) << endl;
#endif
                return false;
            }
            for (size_t j = 0; j < multipath_aln.subpath(i).path().mapping_size(); j++) {
                if (multipath_aln.subpath(i).path().mapping(j).edit_size() == 0) {
#ifdef debug_verbose_validation
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
#ifdef debug_verbose_validation
                    cerr << "validation failure on within-node adjacency" << endl;
                    cerr << pb2json(mapping_from) << "->" << pb2json(mapping_to) << endl;
#endif
                    return false;
                }
            }
            else {
                if (mapping_from_end_offset != handle_graph.get_length(handle_graph.get_handle(mapping_from.position().node_id()))) {
#ifdef debug_verbose_validation
                    cerr << "validation failure on using edge at middle of node" << endl;
                    cerr << pb2json(mapping_from) << "->" << pb2json(mapping_to) << endl;
#endif
                    return false;
                }
                
                handle_t handle_from = handle_graph.get_handle(mapping_from.position().node_id(), mapping_from.position().is_reverse());
                handle_t handle_to = handle_graph.get_handle(mapping_to.position().node_id(), mapping_to.position().is_reverse());
                
                bool found_edge = false;
                function<bool(const handle_t&)> check_for_edge = [&](const handle_t& next_handle) {
                    found_edge = (next_handle == handle_to);
                    return !found_edge;
                };
                handle_graph.follow_edges(handle_from, false, check_for_edge);
                
                if (!found_edge) {
#ifdef debug_verbose_validation
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
            string node_seq = handle_graph.get_sequence(handle_graph.get_handle(mapping.position().node_id()));
            string rev_node_seq = reverse_complement(node_seq);
            size_t node_idx = mapping.position().offset();
            size_t seq_idx = 0;
            for (size_t i = 0; i < mapping.edit_size(); i++) {
                const Edit& edit = mapping.edit(i);
                if (edit_is_match(edit)) {
                    for (size_t j = 0; j < edit.from_length(); j++, node_idx++, seq_idx++) {
                        if ((mapping.position().is_reverse() ? rev_node_seq[node_idx] : node_seq[node_idx]) != subseq[seq_idx]) {
#ifdef debug_verbose_validation
                            cerr << "validation failure on match that does not match for read " << multipath_aln.name() << endl;
                            cerr << pb2json(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                    }
                }
                else if (edit_is_sub(edit)) {
                    for (size_t j = 0; j < edit.from_length(); j++, node_idx++, seq_idx++) {
                        if ((mapping.position().is_reverse() ? rev_node_seq[node_idx] : node_seq[node_idx]) == subseq[seq_idx]) {
#ifdef debug_verbose_validation
                            cerr << "validation failure on mismatch that matches" << endl;
                            cerr << pb2json(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                        if (edit.sequence()[j] != subseq[seq_idx]) {
#ifdef debug_verbose_validation
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
#ifdef debug_verbose_validation
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
    
    void view_multipath_alignment(ostream& out, const MultipathAlignment& multipath_aln, const HandleGraph& handle_graph) {
        
        size_t max_line_length = 128;
        
        vector<pair<int64_t, int64_t>> subpath_read_interval(multipath_aln.subpath_size(), pair<int64_t, int64_t>(0, 0));
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            
            const Subpath& subpath = multipath_aln.subpath(i);
            subpath_read_interval[i].second = subpath_read_interval[i].first + path_to_length(subpath.path());
            
            for (int64_t j : subpath.next()) {
                subpath_read_interval[j].first = subpath_read_interval[i].second;
            }
        }
        
        auto format_position = [](const Position& pos) {
            stringstream strm;
            strm << pos.node_id() << (pos.is_reverse() ? "-" : "+") << (pos.offset() ? (":" + to_string(pos.offset())) : "");
            return strm.str();
        };
        
        for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const Subpath& subpath = multipath_aln.subpath(i);
            
            stringstream read_strm;
            stringstream node_strm;
            
            out << "subpath " << i << " (score " << subpath.score() << ")" << endl;
            
            string read_header = "read[" + to_string(subpath_read_interval[i].first) + ":" + to_string(subpath_read_interval[i].second) + "]";
            string node_header = "graph" + string(max(int(read_header.size()) - 5, 0), ' ');
            
            read_strm << "\t" << read_header;
            node_strm << "\t" << node_header;
            
            size_t line_length = 8 + read_header.size();
            
            int64_t read_at = subpath_read_interval[i].first;
            
            bool first_mapping = true;
            for (size_t j = 0; j < subpath.path().mapping_size(); j++) {
                const Mapping& mapping = subpath.path().mapping(j);
                string pos_string = format_position(mapping.position());
                
                
                stringstream mapping_read_strm;
                stringstream mapping_node_strm;
                mapping_read_strm << pos_string << " ";
                mapping_node_strm << string(pos_string.size() + 1, ' ');
                
                string node_seq = handle_graph.get_sequence(handle_graph.get_handle(mapping.position().node_id(),
                                                                                    mapping.position().is_reverse()));
                
                int64_t node_at = mapping.position().offset();
                for (const Edit& edit : mapping.edit()) {
                    if (edit.from_length() > 0 && edit.to_length() > 0) {
                        mapping_read_strm << multipath_aln.sequence().substr(read_at, edit.to_length());
                        mapping_node_strm << node_seq.substr(node_at, edit.from_length());
                    }
                    else if (edit.from_length() > 0) {
                        mapping_read_strm << string(edit.from_length(), '-');
                        mapping_node_strm << node_seq.substr(node_at, edit.from_length());
                    }
                    else if (edit.to_length() > 0) {
                        mapping_read_strm << multipath_aln.sequence().substr(read_at, edit.to_length());
                        mapping_node_strm << string(edit.to_length(), '-');
                    }
                    
                    read_at += edit.to_length();
                    node_at += edit.from_length();
                }
                
                string mapping_read_str = mapping_read_strm.str();
                string mapping_node_str = mapping_node_strm.str();
                
                if (line_length + mapping_node_str.size() + 1 <= max_line_length || first_mapping) {
                    read_strm << " " << mapping_read_str;
                    node_strm << " " << mapping_node_str;
                    line_length += mapping_node_str.size() + 1;
                    
                    first_mapping = false;
                }
                else {
                    out << read_strm.str() << endl;
                    out << node_strm.str() << endl << endl;
                    
                    read_strm.str("");
                    node_strm.str("");
                    read_strm << "\tread" << string(read_header.size() - 4, ' ') << " " << mapping_read_str;
                    node_strm << "\tgraph" << string(read_header.size() - 5, ' ') << " " << mapping_node_str;
                    
                    line_length = 9 + read_header.size() + mapping_node_str.size();
                }
                
            }
            
            out << read_strm.str() << endl;
            out << node_strm.str() << endl;
            
            for (size_t j = 0; j < subpath.next_size(); j++) {
                out << "\t-> " << subpath.next(j) << endl;
            }
        }
    }
}







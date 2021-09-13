//
//  multipath_alignment.cpp
//

#include "multipath_alignment.hpp"
#include "haplotypes.hpp"
#include <structures/immutable_list.hpp>
#include <structures/min_max_heap.hpp>

#include <type_traits>

//#define debug_multiple_tracebacks
//#define debug_search
//#define debug_trace
//#define debug_verbose_validation
//#define debug_cigar
//#define debug_find_match
//#define debug_remove_empty

using namespace std;
using namespace structures;
using namespace vg::io;

namespace vg {

    multipath_alignment_t::multipath_alignment_t() : _mapping_quality(0) {
        // i don't understand why this is necessary, but somehow mapq is getting defaulted to 160?
    }

    multipath_alignment_t::multipath_alignment_t(const multipath_alignment_t& other) {
        *this = other;
    }

    multipath_alignment_t::multipath_alignment_t(multipath_alignment_t&& other) {
        *this = move(other);
    }

    multipath_alignment_t& multipath_alignment_t::operator=(const multipath_alignment_t& other) {
        _sequence = other._sequence;
        _quality = other._quality;
        _subpath = other._subpath;
        _mapping_quality = other._mapping_quality;
        _start = other._start;
        other.for_each_annotation([&](const string& name, anno_type_t type, const void* annotation) {
            switch (type) {
                case Null:
                    set_annotation(name);
                    break;
                case Double:
                    set_annotation(name, *((const double*) annotation));
                    break;
                case Bool:
                    set_annotation(name, *((const bool*) annotation));
                    break;
                case String:
                    set_annotation(name, *((const string*) annotation));
                    break;
                default:
                    cerr << "error: unrecognized annotation type" << endl;
                    exit(1);
                    break;
            }
        });
        return *this;
    }

    multipath_alignment_t& multipath_alignment_t::operator=(multipath_alignment_t&& other) {
        if (this != &other) {
            _sequence = move(other._sequence);
            _quality = move(other._quality);
            _subpath = move(other._subpath);
            _mapping_quality = move(other._mapping_quality);
            _start = move(other._start);
            _annotation = move(other._annotation);
            other._annotation.clear();
        }
        return *this;
    }

    multipath_alignment_t::~multipath_alignment_t() {
        while (!_annotation.empty()) {
            clear_annotation(_annotation.begin()->first);
        }
    }

    void multipath_alignment_t::clear_annotation(const string& annotation_name) {
        auto iter = _annotation.find(annotation_name);
        if (iter != _annotation.end()) {
            switch (iter->second.first) {
                case Null:
                    break;
                case Double:
                    free((double*) iter->second.second);
                    break;
                case Bool:
                    free((bool*) iter->second.second);
                    break;
                case String:
                    delete ((string*) iter->second.second);
                    break;
                default:
                    cerr << "error: unrecognized annotation type" << endl;
                    exit(1);
                    break;
            }
            _annotation.erase(iter);
        }
    }

    void multipath_alignment_t::set_annotation(const string& annotation_name) {
        clear_annotation(annotation_name);
        _annotation[annotation_name] = make_pair(Null, (void*) nullptr);
    }

    void multipath_alignment_t::set_annotation(const string& annotation_name, double value) {
        clear_annotation(annotation_name);
        auto ptr = (double*) malloc(sizeof(double));
        *ptr = value;
        _annotation[annotation_name] = make_pair(Double, (void*) ptr);
    }

    void multipath_alignment_t::set_annotation(const string& annotation_name, bool value) {
        clear_annotation(annotation_name);
        auto ptr = (bool*) malloc(sizeof(bool));
        *ptr = value;
        _annotation[annotation_name] = make_pair(Bool, (void*) ptr);
    }

    void multipath_alignment_t::set_annotation(const string& annotation_name, const string& value) {
        clear_annotation(annotation_name);
        auto ptr = new string();
        *ptr = value;
        _annotation[annotation_name] = make_pair(String, (void*) ptr);
    }

    pair<multipath_alignment_t::anno_type_t, const void*>
    multipath_alignment_t::get_annotation(const string& annotation_name) const {
        auto iter = _annotation.find(annotation_name);
        if (iter != _annotation.end()) {
            return iter->second;
        }
        else {
            return pair<anno_type_t, void*>(Null, nullptr);
        }
    }

    void multipath_alignment_t::for_each_annotation(function<void(const string&, anno_type_t, const void*)> lambda) const {
        for (const auto& annotation : _annotation) {
            lambda(annotation.first, annotation.second.first, annotation.second.second);
        }
    }
    
    /// Return either the vector of topological order by index or the vector of indexes within the topological order
    vector<size_t> subpath_topological_order(const multipath_alignment_t& multipath_aln,
                                             bool do_index) {
        // Kahn's algorithm
        
        vector<size_t> return_val(multipath_aln.subpath_size(), 0);
        size_t idx = 0;
        
        vector<size_t> stack;
        vector<size_t> in_degree(multipath_aln.subpath_size(), 0);
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
            const subpath_t& subpath = multipath_aln.subpath(i);
            for (size_t j = 0; j < subpath.next_size(); ++j) {
                ++in_degree[subpath.next(j)];
            }
            for (const auto& connection : subpath.connection()) {
                ++in_degree[connection.next()];
            }
        }
        
        // identify the source nodes and add them to the stack
        for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
            if (!in_degree[i]) {
                stack.push_back(i);
            }
        }
        
        while (!stack.empty()) {
            // pop a source node and add it to the topological order
            size_t here = stack.back();
            stack.pop_back();
            
            if (do_index) {
                return_val[here] = idx;
            }
            else {
                return_val[idx] = here;
            }
            ++idx;
            
            // remove the node's edges
            const subpath_t& subpath = multipath_aln.subpath(here);
            for (size_t i = 0; i < subpath.next_size(); ++i) {
                size_t next = subpath.next(i);
                --in_degree[next];
                // if a node is now a source, stack it up
                if (!in_degree[next]) {
                    stack.push_back(next);
                }
            }
            // repeat for connections
            for (const auto& connection : subpath.connection()) {
                --in_degree[connection.next()];
                if (!in_degree[connection.next()]) {
                    stack.push_back(connection.next());
                }
            }
        }
        return return_val;
    }
    
    void topologically_order_subpaths(multipath_alignment_t& multipath_aln) {
        
        vector<size_t> index = subpath_topological_order(multipath_aln, true);
        
        // translate the edges to the new indices
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            subpath_t* subpath = multipath_aln.mutable_subpath(i);
            for (size_t j = 0; j < subpath->next_size(); j++) {
                subpath->set_next(j, index[subpath->next(j)]);
            }
            for (size_t j = 0; j < subpath->connection_size(); ++j) {
                connection_t* connection = subpath->mutable_connection(j);
                connection->set_next(index[connection->next()]);
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

    void remove_empty_alignment_sections(multipath_alignment_t& multipath_aln) {
        
        vector<bool> is_empty(multipath_aln.subpath_size(), false);
        bool any_empty_subpaths = false;
        
        // find subpaths that don't have any aligned bases
        for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
            auto path = multipath_aln.mutable_subpath(i)->mutable_path();
            size_t mappings_removed = 0;
            for (size_t j = 0; j < path->mapping_size(); ++j) {
                auto mapping = path->mutable_mapping(j);
                size_t edits_removed = 0;
                for (size_t k = 0; k < mapping->edit_size(); ++k) {
                    auto edit = mapping->mutable_edit(k);
                    if (edit->from_length() == 0 && edit->to_length() == 0) {
                        ++edits_removed;
                    }
                    else if (edits_removed != 0) {
                        *mapping->mutable_edit(k - edits_removed) = move(*edit);
                    }
                }
                mapping->mutable_edit()->resize(mapping->edit_size() - edits_removed);
                if (mapping->edit().empty()) {
                    ++mappings_removed;
                }
                else if (mappings_removed != 0) {
                    *path->mutable_mapping(j - mappings_removed) = move(*mapping);
                }
            }
            path->mutable_mapping()->resize(path->mapping_size() - mappings_removed);
            is_empty[i] = path->mapping().empty();
            any_empty_subpaths = any_empty_subpaths || path->mapping().empty();
        }
        
#ifdef debug_remove_empty
        cerr << "subpaths empty:" << endl;
        for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
            cerr << "\t" << i << " " << is_empty[i] << endl;
        }
#endif
        
        if (any_empty_subpaths) {
            // there's at least one empty subpath
            
            // compute the transitive edges through empty subpaths
            for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
                if (is_empty[i]) {
                    continue;
                }
                // check if any of this subpaths nexts are empty
                subpath_t& subpath = *multipath_aln.mutable_subpath(i);
                
                set<size_t> transitive_nexts;
                map<size_t, int32_t> transitive_connections;
                
                unordered_set<size_t> dequeued;
                // records of (index, no connections, score)
                priority_queue<tuple<int32_t, bool, size_t>> queue;
                for (auto n : subpath.next()) {
                    if (is_empty[n]) {
                        queue.emplace(0, true, n);
                    }
                }
                for (const auto& connection : subpath.connection()) {
                    if (is_empty[connection.next()]) {
                        queue.emplace(connection.score(), false, connection.next());
                    }
                }
                // dijkstra works to compute max here because it's a DAG
                while (!queue.empty()) {
                    
                    int32_t score;
                    bool not_connection;
                    size_t idx;
                    tie(score, not_connection, idx) = queue.top();
                    queue.pop();
                    
                    if (dequeued.count(idx)) {
                        continue;
                    }
                    dequeued.insert(idx);
                    
                    const auto& subpath_here = multipath_aln.subpath(idx);
                    for (auto next : subpath_here.next()) {
                        if (is_empty[next]) {
                            queue.emplace(score, not_connection, next);
                        }
                        else if (not_connection) {
                            transitive_nexts.insert(next);
                        }
                        else {
                            auto it = transitive_connections.find(next);
                            if (it == transitive_connections.end()) {
                                transitive_connections[next] = score;
                            }
                            else {
                                it->second = max(score, it->second);
                            }
                        }
                    }
                    for (const auto& connection : subpath_here.connection()) {
                        int32_t score_thru = connection.score() + score;
                        if (is_empty[connection.next()]) {
                            queue.emplace(score_thru, false, connection.next());
                        }
                        else {
                            auto it = transitive_connections.find(connection.next());
                            if (it == transitive_connections.end()) {
                                transitive_connections[connection.next()] = score_thru;
                            }
                            else {
                                it->second = max(score_thru, it->second);
                            }
                        }
                    }
                }
#ifdef debug_remove_empty
                cerr << "transitive nexts for subpath " << i << endl;
                for (size_t j : transitive_nexts) {
                    cerr << "\t" << j << endl;
                }
                cerr << "transitive connections for subpath " << i << endl;
                for (auto c : transitive_connections) {
                    cerr << "\t" << c.first << " " << c.second << endl;
                }
#endif
                
                // add edges for the nexts reachable through empty subpaths
                for (size_t j : transitive_nexts) {
                    bool found = false;
                    for (size_t k = 0; k < subpath.next_size() && !found; ++k) {
                        found = (subpath.next(k) == j);
                    }
                    if (!found) {
                        subpath.add_next(j);
                    }
                }
                
                for (const pair<size_t, int32_t>& cnxn : transitive_connections) {
                    bool found = false;
                    size_t k = 0;
                    for (; k < subpath.connection_size() && !found; ++k) {
                        found = (subpath.connection(k).next() == cnxn.first);
                    }
                    if (!found) {
                        auto connection = subpath.add_connection();
                        connection->set_next(cnxn.first);
                        connection->set_score(cnxn.second);
                    }
                    else {
                        subpath.mutable_connection(k - 1)->set_score(max(subpath.connection(k - 1).score(), cnxn.second));
                    }
                }
            }
            
#ifdef debug_remove_empty
            cerr << "before removing subpaths" << endl;
            cerr << debug_string(multipath_aln) << endl;
#endif
            
            // relocate the subpaths we're keeping at the front of the vector
            vector<size_t> removed_so_far(multipath_aln.subpath_size() + 1, 0);
            for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
                if (is_empty[i]) {
                    removed_so_far[i + 1] = removed_so_far[i] + 1;
                    continue;
                }
                else {
                    removed_so_far[i + 1] = removed_so_far[i];
                }
                
                if (removed_so_far[i] > 0) {
                    *multipath_aln.mutable_subpath(i - removed_so_far[i]) = move(*multipath_aln.mutable_subpath(i));
                }
            }
            
            // delete the end of the subpaths
            multipath_aln.mutable_subpath()->resize(multipath_aln.subpath_size() - removed_so_far.back());
            
#ifdef debug_remove_empty
            cerr << "before updating edges" << endl;
            cerr << debug_string(multipath_aln) << endl;
#endif
            
            // reassign the next and connection indexes
            for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
                
                subpath_t& subpath = *multipath_aln.mutable_subpath(i);
                
                size_t nexts_removed_so_far = 0;
                unordered_set<int64_t> nexts_seen;
                for (size_t j = 0; j < subpath.next_size(); ++j) {
                    int64_t updated_next = subpath.next(j) - removed_so_far[subpath.next(j)];
                    if (is_empty[subpath.next(j)] || nexts_seen.count(updated_next)) {
                        ++nexts_removed_so_far;
                    }
                    else {
                        subpath.set_next(j - nexts_removed_so_far, updated_next);
                        nexts_seen.insert(updated_next);
                    }
                }
                if (nexts_removed_so_far) {
                    subpath.mutable_next()->resize(subpath.next_size() - nexts_removed_so_far);
                }
                
                size_t connections_removed_so_far = 0;
                unordered_set<pair<size_t, int64_t>> connections_seen;
                for (size_t j = 0; j < subpath.connection_size(); ++j) {
                    auto connection = subpath.mutable_connection(j);
                    auto updated_connection = pair<size_t, int64_t>(connection->next() - removed_so_far[connection->next()],
                                                                    connection->score());

                    if (is_empty[connection->next()] || connections_seen.count(updated_connection)) {
                        ++connections_removed_so_far;
                    }
                    else {
                        connection->set_next(updated_connection.first);
                        if (connections_removed_so_far) {
                            *subpath.mutable_connection(j - connections_removed_so_far) = *connection;
                        }
                        connections_seen.insert(updated_connection);
                    }
                }
                if (connections_removed_so_far) {
                    subpath.mutable_connection()->resize(subpath.connection_size() - connections_removed_so_far);
                }
            }
            
#ifdef debug_remove_empty
            cerr << "before updating starts" << endl;
            cerr << debug_string(multipath_aln) << endl;
#endif
            
            // update the starts
            bool found_deleted_start = false;
            for (size_t i = 0; !found_deleted_start && i < multipath_aln.start_size(); ++i) {
                found_deleted_start = is_empty[multipath_aln.start(i)];
                multipath_aln.set_start(i, multipath_aln.start(i) - removed_so_far[multipath_aln.start(i)]);
            }
            
            if (found_deleted_start) {
                // recompute the edges from scratch (could be done faster, but
                // this is easy);
                identify_start_subpaths(multipath_aln);
            }
        }
    }
    
    void identify_start_subpaths(multipath_alignment_t& multipath_aln) {
        
        // remove start nodes if there are any (avoids doubling them if this function is used liberally)
        multipath_aln.clear_start();
        
        // label nodes with incoming edges
        vector<bool> has_incoming_edge(multipath_aln.subpath_size(), false);
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const subpath_t& subpath = multipath_aln.subpath(i);
            for (size_t j = 0; j < subpath.next_size(); j++) {
                has_incoming_edge[subpath.next(j)] = true;
            }
            for (const connection_t& connection : subpath.connection()) {
                has_incoming_edge[connection.next()] = true;
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
        MultipathProblem(const multipath_alignment_t& multipath_aln, bool subpath_global)
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
    tuple<MultipathProblem, int64_t, int32_t> run_multipath_dp(const multipath_alignment_t& multipath_aln,
                                                               bool subpath_global = false) {
        
        // Create and unpack the return value (including setting up the DP table). Initialise score according
        // to whether the alignment is local or global 
        tuple<MultipathProblem, int64_t, int32_t> to_return(MultipathProblem(multipath_aln, subpath_global),
                                                            -1, subpath_global ? numeric_limits<int32_t>::min() : 0);
        auto& problem = get<0>(to_return);
        auto& opt_subpath = get<1>(to_return);
        auto& opt_score = get<2>(to_return);
    
        for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
            const subpath_t& subpath = multipath_aln.subpath(i);
            int32_t extended_score = problem.prefix_score[i] + subpath.score();
            // carry DP forward
            if (subpath.next_size() != 0 || subpath.connection_size() != 0) {
                int64_t thru_length = path_to_length(subpath.path()) + problem.prefix_length[i];
                for (size_t j = 0; j < subpath.next_size(); ++j) {
                    int64_t next = subpath.next(j);
                    problem.prefix_length[next] = thru_length;
                    
                    // can we improve prefix score on following subpath through this one?
                    if (extended_score >= problem.prefix_score[next]) {
                        problem.prev_subpath[next] = i;
                        problem.prefix_score[next] = extended_score;
                    }
                }
                // repeat DP across scored connections
                for (size_t j = 0; j < subpath.connection_size(); ++j) {
                    const connection_t& connection = subpath.connection(j);
                    
                    problem.prefix_length[connection.next()] = thru_length;
                    if (extended_score + connection.score() >= problem.prefix_score[connection.next()]) {
                        problem.prev_subpath[connection.next()] = i;
                        problem.prefix_score[connection.next()] = extended_score + connection.score();
                    }
                }
            }
            // check if an alignment is allowed to end here according to global/local rules and
            // if so whether it's optimal
            if (extended_score >= opt_score && (!subpath_global || (subpath.next_size() == 0 &&
                                                                    subpath.connection_size() == 0))) {
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
    void populate_path_from_traceback(const multipath_alignment_t& multipath_aln, const MultipathProblem& problem,
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
            
            Position* position = soft_clip_mapping->mutable_position();
            const position_t& pos_from = multipath_aln.subpath(*current_subpath).path().mapping(0).position();
            position->set_node_id(pos_from.node_id());
            position->set_is_reverse(pos_from.is_reverse());
            position->set_offset(pos_from.offset());
        }
        
        // merge the subpaths into one optimal path in the Alignment object
        for (auto next_subpath = current_subpath; next_subpath != traceback_end; ++next_subpath) {
            // For all subpaths in the traceback
            
            // Advance only if we don't hit the end
            current_subpath = next_subpath;
            
            if (output->mapping_size() == 0) {
                // There's nothing in the output yet, so just copy all the mappings from this subpath
                to_proto_path(multipath_aln.subpath(*current_subpath).path(), *output);
                
                for(size_t i = 0; i < output->mapping_size(); i++) {
                    // Set all the ranks
                    output->mutable_mapping(i)->set_rank(i + 1);
                }
            } else {
                // There's already content in the output so we have to merge stuff
                Mapping* curr_end_mapping = output->mutable_mapping(output->mapping_size() - 1);
            
                // get the first mapping of the next path
                const path_t& next_path = multipath_aln.subpath(*current_subpath).path();
                const path_mapping_t& next_start_mapping = next_path.mapping(0);
                
                size_t mapping_start_idx = 0;
                // merge mappings if they occur on the same node and same strand
                if (curr_end_mapping->position().node_id() == next_start_mapping.position().node_id()
                    && curr_end_mapping->position().is_reverse() == next_start_mapping.position().is_reverse()) {
                    
                    Edit* last_edit = curr_end_mapping->mutable_edit(curr_end_mapping->edit_size() - 1);
                    const edit_t& first_edit = next_start_mapping.edit(0);
                    
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
                        to_proto_edit(next_start_mapping.edit(j), *curr_end_mapping->add_edit());
                    }
                    
                    mapping_start_idx++;
                }
                
                // append the rest of the mappings
                for (size_t j = mapping_start_idx; j < next_path.mapping_size(); j++) {
                    Mapping* next_mapping = output->add_mapping();
                    to_proto_mapping(next_path.mapping(j), *next_mapping);
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
    
    int32_t optimal_alignment_internal(const multipath_alignment_t& multipath_aln, Alignment* aln_out,
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
            deque<int64_t> opt_traceback;
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
    
    void optimal_alignment(const multipath_alignment_t& multipath_aln, Alignment& aln_out, bool subpath_global) {

        // transfer read information over to alignment
        transfer_read_metadata(multipath_aln, aln_out);
        aln_out.set_mapping_quality(multipath_aln.mapping_quality());
        
        // do dynamic programming and traceback the optimal alignment
        int32_t score = optimal_alignment_internal(multipath_aln, &aln_out, subpath_global);
        
        aln_out.set_score(score);
    }
    
    int32_t optimal_alignment_score(const multipath_alignment_t& multipath_aln, bool subpath_global){
        // do dynamic programming without traceback
        return optimal_alignment_internal(multipath_aln, nullptr, subpath_global);
    }

    int32_t worst_alignment_score(const multipath_alignment_t& multipath_aln) {
        
        if (multipath_aln.subpath().empty()) {
            return 0;
        }
        
        // initialize a DP table
        vector<int32_t> dp(multipath_aln.subpath_size(), numeric_limits<int32_t>::max());
        
        // initial conditions, allow alignments to begin at the starts
        for (auto i : multipath_aln.start()) {
            dp[i] = 0;
        }
        
        int32_t opt = numeric_limits<int32_t>::max();
        for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
            const auto& subpath = multipath_aln.subpath(i);
            int32_t score_thru = dp[i] + subpath.score();
            if (subpath.next().empty()) {
                // this is a sink, check for optimality
                opt = min(opt, score_thru);
            }
            else {
                // carry the DP through to the next subpaths
                for (auto j : subpath.next()) {
                    dp[j] = min(dp[j], score_thru);
                }
            }
        }
        return max<int32_t>(opt, 0);
    }
    
    vector<Alignment> optimal_alignments(const multipath_alignment_t& multipath_aln, size_t count) {
        
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
        vector<vector<pair<int64_t, int32_t>>> prev_subpaths(multipath_aln.subpath_size());

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
        
        for (int64_t i = 0; i < multipath_aln.subpath_size(); ++i) {
            // For each subpath
            
            // If it has no successors, we can start a traceback here
            bool valid_traceback_start = true;
            
            for (auto& next_subpath : multipath_aln.subpath(i).next()) {
                // For each next subpath it lists
                
                // Register this subpath as a predecessor of the next
                prev_subpaths[next_subpath].emplace_back(i, 0);
                
                if (multipath_aln.subpath(next_subpath).score() >= 0) {
                    // This successor has a nonnegative score, so taking it
                    // after us would generate a longer, same- or
                    // higher-scoring alignment. So we shouldn't start a
                    // traceback from subpath i.
                    valid_traceback_start = false;
                }
            }
            
            for (const auto& connection : multipath_aln.subpath(i).connection()) {
                
                // register the connection
                prev_subpaths[connection.next()].emplace_back(i, connection.score());
                
                if (multipath_aln.subpath(connection.next()).score() + connection.score() >= 0) {
                    // Taking the connection would lead to a longer or better alignment
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
                
            }
            else {
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
                    auto prev_opt_score = problem.prefix_score[prev.first] + multipath_aln.subpath(prev.first).score() + prev.second;
                    
                    // What's the difference we would take if we went with this predecessor?
                    auto additional_penalty = best_prefix_score - prev_opt_score;
                    
                    destinations.emplace_back(prev.first, additional_penalty);
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
    
    vector<Alignment> optimal_alignments_with_disjoint_subpaths(const multipath_alignment_t& multipath_aln, size_t count) {
        
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
        vector<vector<pair<int64_t, int32_t>>> prev_subpaths(multipath_aln.subpath_size());
        
        for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
            // For each subpath
            
            // If it has no successors, we can start a traceback here
            bool valid_traceback_start = true;
            
            for (auto& next_subpath : multipath_aln.subpath(i).next()) {
                // For each next subpath it lists
                
                // Register this subpath as a predecessor of the next
                prev_subpaths[next_subpath].emplace_back(i, 0);
                
                if (multipath_aln.subpath(next_subpath).score() >= 0) {
                    // This successor has a nonnegative score, so taking it
                    // after us would generate a longer, same- or
                    // higher-scoring alignment. So we shouldn't start a
                    // traceback from subpath i.
                    valid_traceback_start = false;
                }
            }
            
            for (const auto& connection : multipath_aln.subpath(i).connection()) {
                
                // register the connection
                prev_subpaths[connection.next()].emplace_back(i, connection.score());
                
                if (multipath_aln.subpath(connection.next()).score() + connection.score() >= 0) {
                    // Taking the connection would lead to a longer or better alignment
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
                        
                        if (subpath_is_used[prev.first]) {
                            // This subpath has already been used in an emitted alignment, so we can't use it.
                            
#ifdef debug_multiple_tracebacks
                            cerr << "\tSkip " << prev << " which is already used" << endl;
#endif
                            
                            continue;
                        }
                        
                        // For each, compute the score of the optimal alignment ending at that predecessor
                        auto prev_opt_score = problem.prefix_score[prev.first] + multipath_aln.subpath(prev.first).score() + prev.second;
                        
                        // What's the difference we would take if we went with this predecessor?
                        auto additional_penalty = best_prefix_score - prev_opt_score;
                        
                        // Calculate the score differences from optimal
                        auto total_penalty = basis_score_difference + additional_penalty;
                        
                        if (total_penalty >= min_penalty_for_subpath[prev.first]) {
                            // This previous subpath is already reachable with a penalty as good or better.
                            // Don't bother with it again
                            
#ifdef debug_multiple_tracebacks
                            cerr << "\tSkip " << prev.first << " with penalty " << total_penalty << " >= " << min_penalty_for_subpath[prev.first] << endl;
#endif
                            
                            continue;
                        }
                        
                        // Record that this is the cheapest we managed to get here
                        min_penalty_for_subpath[prev.first] = total_penalty;
                        
                        // Make an extended path
                        auto extended_path = basis.push_front(prev.first);
                        
#ifdef debug_multiple_tracebacks
                        cerr << "\tAugment with " << prev.first << " to penalty " << total_penalty << endl;
#endif
                        
                        // Put them in the priority queue
                        queue.push(make_pair(total_penalty, extended_path));
                    }
                    
                    
                }
            }
        }
        
        return to_return;
        
    }
   
    vector<Alignment> haplotype_consistent_alignments(const multipath_alignment_t& multipath_aln, const haplo::ScoreProvider& score_provider,
        size_t soft_count, size_t hard_count, bool optimal_first) {
        
#ifdef debug_multiple_tracebacks
        cerr << "Computing haplotype consistent alignments" << endl;
#endif

        // We can only work with a score provider that supports incremental search.
        assert(score_provider.has_incremental_search());
        
        // Keep a list of what we're going to emit.
        vector<Alignment> to_return;
        
        // Fill out the dynamic programming problem
        // TODO: are we duplicating work if we also get the top alignment?
        auto dp_result = run_multipath_dp(multipath_aln);
        // Get the filled DP problem
        MultipathProblem& problem = get<0>(dp_result);
        // And the optimal final subpath
        int64_t& opt_subpath = get<1>(dp_result);
        // And the optimal score
        int32_t& opt_score = get<2>(dp_result);
        
        if (optimal_first) {
            // Compute the optimal alignment and put it first.
            // TODO: It will also appear later if it is haplotype-consistent.
            // But we are allowed to produce duplicates so that's OK.
            to_return.emplace_back();
            Alignment& opt_aln = to_return.back();
            
            opt_aln.set_score(opt_score);
            
            // traceback the optimal subpaths until hitting sentinel (-1)
            list<int64_t> opt_traceback;
            int64_t curr = opt_subpath;
            while (curr >= 0) {
                opt_traceback.push_front(curr);
                curr = problem.prev_subpath[curr];
            }
            
            Path* opt_path = opt_aln.mutable_path();
            
            // Fill in the path in the alignment with the alignment represented
            // by this traceback in this DP problem for this multipath
            // alignment.
            populate_path_from_traceback(multipath_aln, problem, opt_traceback.begin(), opt_traceback.end(), opt_path);
            
#ifdef debug_multiple_tracebacks
            cerr << "Produced optimal alignment with score " << opt_aln.score() << endl;
#endif
        }
        
        // Keep lists of traceback steps as multipath subpath numbers
        using step_list_t = ImmutableList<int64_t>;
        
        // We define our own search state, which includes a haplotype search
        // state and a flag for whether all the edges crossed so far have been
        // present in the GBWT. This flag can remain true, and we can keep
        // searching, after hour haplotype search state becomes empty.
        struct SearchState {
            haplo::IncrementalSearchState haplo_state;
            bool all_edges_exist = true;
            bool started = false;
            
            inline bool empty() const {
                return haplo_state.empty();
            }
            
            inline size_t size() const {
                return haplo_state.size();
            }
            
            // More consistent things should be smaller, for good queueing
            inline bool operator<(const SearchState& other) const {
                return haplo_state.size() > other.haplo_state.size() ||
                    (haplo_state.size() == other.haplo_state.size() && all_edges_exist && !other.all_edges_exist);
            }
        };
        
        // This function advances an incremental haplotype search state with the edges along a subpath, if any.
        // A non-started input SearchState means to start the search.
        // Note that we interpret the path IN REVERSE, because we're doing a traceback.
        auto extend_with_subpath = [&](const SearchState& initial, int64_t subpath) {
            // Get the Path from the subpath.
            const path_t& path = multipath_aln.subpath(subpath).path();
            
            // No empty paths are allowed.
            assert(path.mapping_size() > 0);
            
            // Set up a search state scratch
            SearchState state = initial;
            if (!state.started) {
                // If our input state says we need to start a new search, start with the node from the last mapping.
                auto& pos = path.mapping(path.mapping_size() - 1).position();
                // We require everything to have mappings to actual places, even pure inserts.
                assert(pos.node_id() != 0);
                // Make sure to search in the orientation we are actually going
                state.haplo_state = score_provider.incremental_find(make_position(pos.node_id(), !pos.is_reverse(), 0));
                state.started = true;
                
#ifdef debug_multiple_tracebacks
                cerr << "New haplotype search for " << pb2json(pos) << " finds " << state.size() << " matching haplotypes" << endl;
#endif
            }
            
            // Otherwise we have already selected the last Mapping when we crossed the edge back into here.
            for (size_t i = path.mapping_size() - 1; i != 0; i--) {
                // For each transition between Mappings, we assume we are going between distinct node visits because of our precondition.
                // So find the next position looking left
                auto& pos = path.mapping(i - 1).position();
                
                if (!state.empty()) {
                    // Search the transition to it in the reverse orientation.
                    state.haplo_state = score_provider.incremental_extend(state.haplo_state,
                        make_position(pos.node_id(), !pos.is_reverse(), 0));
#ifdef debug_multiple_tracebacks
                    cerr << "Extend within subpath " << subpath
                        << " by going to node " << pos.node_id()
                        << " matches " << state.size() << " haplotypes" << endl;
#endif
                } else {
#ifdef debug_multiple_tracebacks
                    cerr << "Extend within subpath " << subpath
                        << " by going to node " << pos.node_id()
                        << " but haplotype match set is already empty" << endl;
#endif
                }
                
                if (state.empty() && state.all_edges_exist) {
                    // We have run out of haplotypes. It may be because we have traversed an edge not in the haplotype index.
                    // Check for that.
                    
                    // Look up where we came from
                    auto& prev_pos = path.mapping(i).position();
                    auto scratch = score_provider.incremental_find(make_position(prev_pos.node_id(), !prev_pos.is_reverse(), 0));
                    
                    // Search where we go to
                    scratch = score_provider.incremental_extend(scratch, make_position(pos.node_id(), !pos.is_reverse(), 0));
                    
                    if (scratch.empty()) {
                        // If no haplotypes go there, mark the state as having crossed an unused edge
                        state.all_edges_exist = false;
                        
#ifdef debug_multiple_tracebacks
                        cerr << "Extend within subpath " << subpath
                            << " by going node " << prev_pos.node_id() << " -> " << pos.node_id()
                            << " has no haplotypes crossing that edge" << endl;
#endif
                    }
                }
            
                if (state.empty() && !state.all_edges_exist) {
                    // If we enter this state we can never imporve, so return
                    return state;
                }
                
                // Otherwise loop until we run out of transitions
            }
            
            return state;
        };
        
        // This function advances an incremental haplotype search state with the edge between two subpaths, if any.
        // An empty input subpath means to start the search.
        // Note that we interpret the path IN REVERSE, because we're doing a traceback.
        // Even though old_subpath comes before new_subpath, since we're going backward, new_subpath comes in on the left.
        auto extend_between_subpaths = [&](const SearchState& initial, int64_t old_subpath, int64_t new_subpath) {
            // We can't have an un-started input state
            assert(initial.started);
            
            // See if the transition from the previous subpath to the next subpath is just two mappings abutting on the same node
            const path_t& old_path = multipath_aln.subpath(old_subpath).path();
            const path_t& new_path = multipath_aln.subpath(new_subpath).path();
            assert(old_path.mapping_size() > 0);
            assert(new_path.mapping_size() > 0);
            // We're going left from the first mapping on the old path
            const path_mapping_t& old_mapping = old_path.mapping(0);
            // And into the last mapping on the new path.
            const path_mapping_t& new_mapping = new_path.mapping(new_path.mapping_size() - 1);
            
            // Look up the positions
            auto& new_pos = new_mapping.position();
            auto& old_pos = old_mapping.position();
            
            if (new_pos.node_id() == old_pos.node_id() &&
                new_pos.is_reverse() == old_pos.is_reverse() &&
                new_pos.offset() + mapping_from_length(new_mapping) == old_pos.offset()) {
                // We actually are transitioning just within a node. No more state updates to do.
                
#ifdef debug_multiple_tracebacks
                cerr << "Extend between subpaths " << old_subpath << " and " << new_subpath
                    << " is just abutment on node " << new_pos.node_id() << endl;
#endif
                return initial;
            }
            
            // Otherwise there's an edge we have to look for.
            SearchState result = initial;
            
            if (!result.empty()) {
                // Try searching on the contained interval
                // Make sure to flip the orientation because we're searching left.
                result.haplo_state = score_provider.incremental_extend(initial.haplo_state,
                    make_position(new_mapping.position().node_id(), !new_mapping.position().is_reverse(), 0));
                    
#ifdef debug_multiple_tracebacks
                cerr << "Extend between subpaths " << old_subpath << " and " << new_subpath
                    << " by going node " << old_pos.node_id() << " (" << initial.size() << ") -> "
                    << new_pos.node_id() << " (" << result.size() << ")" << endl;
#endif
            } else {
#ifdef debug_multiple_tracebacks
                cerr << "Extend between subpaths " << old_subpath << " and " << new_subpath
                    << " by going node " << old_pos.node_id() << " -> " << new_pos.node_id()
                    << " starts from empty search" << endl;
#endif
            }
            
            if (result.empty() && result.all_edges_exist) {
                // We need to see if we not only ran out of agreeing haplotypes but also took an unused edge
                // Look up where we came from
                auto scratch = score_provider.incremental_find(make_position(old_pos.node_id(), !old_pos.is_reverse(), 0));
                
                // Search where we go to
                scratch = score_provider.incremental_extend(scratch, make_position(new_pos.node_id(), !new_pos.is_reverse(), 0));
                
                if (scratch.empty()) {
                    // If no haplotypes go there, mark the state as having crossed an unused edge
                    result.all_edges_exist = false;
                    
#ifdef debug_multiple_tracebacks
                    cerr << "Extend between subpaths " << old_subpath << " and " << new_subpath
                        << " by going node " << old_pos.node_id() << " -> " << new_pos.node_id()
                        << " has no haplotypes crossing that edge" << endl;
#endif
                }
                
            }
            
            return result;
        };
        
        // Our search is kind of complicated, because we want to enumerate all
        // haplotype-consistent linearizations, but pad out to n merely
        // scorable linearizations, in alignment score order.
        
        // So we keep our search queue in a min-max heap, where lower is a
        // better thing to extend.
        //
        // We sort by search state, which we define an order on where more
        // consistent haplotypes come before fewer, and then unscorable things
        // come last.
        // 
        // And then after that we search by score penalty from optimal.
        
        // Put them in a size-limited priority queue by search state, and then
        // score difference (positive) from optimal
        MinMaxHeap<tuple<SearchState, int32_t, step_list_t>> queue;
        
        // We define a function to put stuff in the queue and limit its size to
        // the (count - to_return.size()) items with lowest penalty.
        auto try_enqueue = [&](const tuple<SearchState, int32_t, step_list_t>& item) {
            // Work out how many things can be in the queue to compete to pad out the remaining return slots.
            // Make sure it doesn't try to go negative.
            size_t soft_max_size = soft_count - std::min(soft_count, to_return.size());
            size_t hard_max_size = hard_count ? hard_count - to_return.size() : numeric_limits<size_t>::max();

#ifdef debug_multiple_tracebacks
            if (queue.size() >= hard_max_size) {
                cerr << "We've reached the hard cap on queue size -- even haplotype consistent are rejected" << endl;
            } else if (!get<0>(item).empty()) {
                cerr << "Item is haplotype-consistent and must be queued" << endl;
            } else if (queue.size() < soft_max_size) {
                cerr << "Item fits in queue" << endl;
            } else if (!queue.empty() && item < queue.max()) {
                cerr << "Item beats worst thing in queue" << endl;
            }
#endif
            
            
            if ((!get<0>(item).empty() && queue.size() < hard_max_size)
                || (queue.size() < soft_max_size && queue.size() < hard_max_size)
                || (!queue.empty() && item < queue.max())) {
                // The item belongs in the queue because it fits or it beats
                // the current worst thing if present, or it's not eligible for removal.
                queue.push(item);
#ifdef debug_multiple_tracebacks
                cerr << "Allow item into queue (" << queue.size() << "/" << soft_max_size << "," << hard_max_size << ")" << endl;
#endif
                
            }
            
            while (!queue.empty()
                   && (queue.size() > hard_max_size || queue.size() > soft_max_size)
                   && (get<0>(queue.max()).empty() || queue.size() > hard_max_size)) {
                // We have more possibilities than we need to consider, and
                // some are eligible for removal. Get rid of the worst one.
                queue.pop_max();
#ifdef debug_multiple_tracebacks
                cerr << "Remove worst from queue (" << queue.size() << "/" << soft_max_size << "," << hard_max_size << ")" << endl;
#endif
            }
        };
        
        // Also, subpaths only keep track of their nexts, so we need to invert
        // that so we can get all valid prev subpaths.
        // TODO: This code is also duplicated
        vector<vector<pair<int64_t, int32_t>>> prev_subpaths;
        
        prev_subpaths.resize(multipath_aln.subpath_size());
        for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
            // For each subpath
            
            // If it has no successors, we can start a traceback here
            bool valid_traceback_start = true;
            
            for (auto& next_subpath : multipath_aln.subpath(i).next()) {
                // For each next subpath it lists
                
                // Register this subpath as a predecessor of the next
                prev_subpaths[next_subpath].emplace_back(i, 0);
                
                if (multipath_aln.subpath(next_subpath).score() >= 0) {
                    // This successor has a nonnegative score, so taking it
                    // after us would generate a longer, same- or
                    // higher-scoring alignment. So we shouldn't start a
                    // traceback from subpath i.
                    valid_traceback_start = false;
                }
            }
            
            for (const auto& connection : multipath_aln.subpath(i).connection()) {
                
                // register the connection
                prev_subpaths[connection.next()].emplace_back(i, connection.score());
                
                if (multipath_aln.subpath(connection.next()).score() + connection.score() >= 0) {
                    // Taking the connection would lead to a longer or better alignment
                    valid_traceback_start = false;
                }
            }
            
            if (valid_traceback_start) {
                // We can start a traceback here.
                
                // The path is just to be here
                step_list_t starting_path{i};
                
                // The search state is what we get just starting with this subpath
                SearchState state = extend_with_subpath(SearchState(), i);
                
                // The score penalty for starting here is the optimal score minus the optimal score starting here
                auto penalty = opt_score - (problem.prefix_score[i] + multipath_aln.subpath(i).score());
                
                
#ifdef debug_multiple_tracebacks
                cerr << "Could end at subpath " << i << " with " << state.size()
                    << " matches and scorability " << state.all_edges_exist << endl;
#endif
                
                try_enqueue(make_tuple(state, penalty, starting_path)); 
                
            }
        }
        
        while(!queue.empty()) {
            // Grab a traceback to try extending
            auto frame = queue.min();
            queue.pop_min();
            
            // Unpack the stack frame
            auto& state = get<0>(frame);
            auto& base_penalty = get<1>(frame);
            auto& basis = get<2>(frame);
            
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
                
                // Compute the score from the penalty
                aln_out.set_score(opt_score - base_penalty);
                
#ifdef debug_multiple_tracebacks
                cerr << "Traceback reaches start at " << basis.front() << " with " << state.size()
                    << " consistent haplotypes and scorability " << state.all_edges_exist << "; emit linearization "
                    << (to_return.size() - 1) << " with score " << aln_out.score() << endl;
                    
                for (auto& m : aln_out.path().mapping()) {
                    cerr << m.position().node_id() << " ";
                }
                cerr << endl;
#endif
            } else {
                // We can't optimally stop the traceback here. We have to come from somewhere.
                
                auto& here = basis.front();
                
                // To compute the additional score difference, we need to know what our optimal prefix score was.
                auto& best_prefix_score = problem.prefix_score[here];
                
                for (auto& prev : prev_subpaths[here]) {
                    // For each possible previous location
                    
                    // Compute the score of the optimal alignment ending at that predecessor
                    auto prev_opt_score = problem.prefix_score[prev.first] + multipath_aln.subpath(prev.first).score() + prev.second;
                    
                    // What's the difference we would take if we went with this predecessor?
                    auto additional_penalty = best_prefix_score - prev_opt_score;
                    
                    // Try extending into the subpath
                    auto extended_state = extend_between_subpaths(state, here, prev.first);
                    // And then with all of it
                    extended_state = extend_with_subpath(extended_state, prev.first);
                    
#ifdef debug_multiple_tracebacks
                    cerr << "Extending traceback from subpath " << here << " to and through " << prev.first << " keeps "
                        << extended_state.size() << " / " << state.size() << " haplotype matches, with scorability "
                        << extended_state.all_edges_exist << " and penalty " << base_penalty + additional_penalty << endl;
#endif
                    
                    // Save the result
                    try_enqueue(make_tuple(extended_state, base_penalty + additional_penalty, basis.push_front(prev.first)));
                }
            }
            
            
            
        }
            
        return to_return;
    }


    pair<int64_t, int64_t> aligned_interval(const multipath_alignment_t& multipath_aln) {
        
        if (multipath_aln.subpath().empty()) {
            return pair<int64_t, int64_t>(0, 0);
        }
        
        int64_t min_softclip_left = numeric_limits<int64_t>::max();
        int64_t min_softclip_right = numeric_limits<int64_t>::max();
        
        for (auto i : multipath_aln.start()) {
            const auto& edit = multipath_aln.subpath(i).path().mapping(0).edit(0);
            if (edit.from_length() == 0 && edit.to_length() != 0) {
                min_softclip_left = min<int64_t>(min_softclip_left, edit.to_length());
            }
            else {
                min_softclip_left = 0;
            }
        }
        
        vector<bool> is_sink(multipath_aln.subpath_size(), false);
        for (const auto& subpath : multipath_aln.subpath()) {
            if (subpath.next_size() == 0 && subpath.connection_size() == 0) {
                
                const auto& path = subpath.path();
                const auto& mapping = path.mapping(path.mapping_size() - 1);
                const auto& edit = mapping.edit(mapping.edit_size() - 1);
                if (edit.from_length() == 0 && edit.to_length() != 0) {
                    min_softclip_right = min<int64_t>(min_softclip_right, edit.to_length());
                }
                else {
                    min_softclip_right = 0;
                }
            }
        }
        
        if (min_softclip_left == numeric_limits<int64_t>::max()) {
            min_softclip_left = 0;
        }
        if (min_softclip_right == numeric_limits<int64_t>::max()) {
            min_softclip_right = 0;
        }
        return pair<int64_t, int64_t>(min_softclip_left,
                                      multipath_aln.sequence().size() - min_softclip_right);
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
    inline void rev_comp_subpath(const subpath_t& subpath, const function<int64_t(int64_t)>& node_length,
                                 subpath_t& rev_comp_out) {
        
        *(rev_comp_out.mutable_path()) = reverse_complement_path(subpath.path(), node_length);
        rev_comp_out.set_score(subpath.score());
        // leave reversing the edges to the multipath alignment
    }
    
    void rev_comp_multipath_alignment(const multipath_alignment_t& multipath_aln, const function<int64_t(int64_t)>& node_length,
                                      multipath_alignment_t& rev_comp_out) {
        
        // reverse complement sequence
        rev_comp_out.set_sequence(reverse_complement(multipath_aln.sequence()));
        // reverse base qualities
        rev_comp_out.set_quality(string(multipath_aln.quality().rbegin(), multipath_aln.quality().rend()));
        
        // transfer the rest of the metadata directly
        rev_comp_out.set_mapping_quality(multipath_aln.mapping_quality());
        
        vector<vector<size_t>> reverse_edge_lists(multipath_aln.subpath_size());
        vector<vector<pair<size_t, int32_t>>> reverse_connection_lists(multipath_aln.subpath_size());
        vector<size_t> reverse_starts;
        
        // remove subpaths to avoid duplicating
        rev_comp_out.clear_subpath();
        
        // add subpaths in reverse order to maintain topological ordering
        for (int64_t i = multipath_aln.subpath_size() - 1; i >= 0; i--) {
            const subpath_t& subpath = multipath_aln.subpath(i);
            subpath_t* rc_subpath = rev_comp_out.add_subpath();
            rev_comp_subpath(subpath, node_length, *rc_subpath);
            
            if (subpath.next_size() > 0 || subpath.connection_size() > 0) {
                // collect edges by their target (for reversing)
                for (size_t j = 0; j < subpath.next_size(); j++) {
                    reverse_edge_lists[subpath.next(j)].push_back(i);
                }
                for (const connection_t& connection : subpath.connection()) {
                    reverse_connection_lists[connection.next()].emplace_back(i, connection.score());
                }
            }
            else {
                // sink subpaths become sources in reverse
                reverse_starts.push_back(i);
            }
        }
        
        // add reversed edges
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            subpath_t* rc_subpath = rev_comp_out.mutable_subpath(i);
            vector<size_t>& reverse_edge_list = reverse_edge_lists[multipath_aln.subpath_size() - i - 1];
            for (size_t j = 0; j < reverse_edge_list.size(); j++) {
                rc_subpath->add_next(multipath_aln.subpath_size() - reverse_edge_list[j] - 1);
            }
            vector<pair<size_t, int32_t>>& reverse_connection_list = reverse_connection_lists[multipath_aln.subpath_size() - i - 1];
            for (size_t j = 0; j < reverse_connection_list.size(); ++j) {
                connection_t* connection = rc_subpath->add_connection();
                connection->set_next(multipath_aln.subpath_size() - reverse_connection_list[j].first - 1);
                connection->set_score(reverse_connection_list[j].second);
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
    
    void rev_comp_multipath_alignment_in_place(multipath_alignment_t* multipath_aln,
                                               const function<int64_t(int64_t)>& node_length) {
        
        // reverse complement sequence
        reverse_complement_in_place(*multipath_aln->mutable_sequence());
        // reverse base qualities
        string* quality = multipath_aln->mutable_quality();
        std::reverse(quality->begin(), quality->end());
        
        // current reverse edges
        vector<vector<uint32_t>> reverse_edge_lists(multipath_aln->subpath_size());
        vector<vector<connection_t>> reverse_connection_lists(multipath_aln->subpath_size());
        // current sink nodes (will be starts)
        vector<uint32_t> reverse_starts;
        
        uint32_t subpath_swap_size = multipath_aln->subpath_size() / 2;
        uint32_t last = multipath_aln->subpath_size() - 1;
        for (uint32_t i = 0, j = last; i < subpath_swap_size; i++, j--) {
            subpath_t* subpath_1 = multipath_aln->mutable_subpath(i);
            subpath_t* subpath_2 = multipath_aln->mutable_subpath(j);
            
            // add reverse edges for first subpath
            if (subpath_1->next_size() > 0 || subpath_1->connection_size() > 0) {
                for (uint32_t k = 0; k < subpath_1->next_size(); k++) {
                    reverse_edge_lists[subpath_1->next(k)].push_back(last - i);
                }
                for (uint32_t k = 0; k < subpath_1->connection_size(); k++) {
                    const connection_t& connection = subpath_1->connection(k);
                    reverse_connection_lists[connection.next()].emplace_back();
                    connection_t& new_connection = reverse_connection_lists[connection.next()].back();
                    new_connection.set_next(last - i);
                    new_connection.set_score(connection.score());
                }
            }
            else {
                reverse_starts.push_back(last - i);
            }
            
            // add reverse edges for second subpath
            if (subpath_2->next_size() > 0 || subpath_2->connection_size() > 0) {
                for (uint32_t k = 0; k < subpath_2->next_size(); k++) {
                    reverse_edge_lists[subpath_2->next(k)].push_back(last - j);
                }
                for (uint32_t k = 0; k < subpath_2->connection_size(); k++) {
                    const connection_t& connection = subpath_2->connection(k);
                    reverse_connection_lists[connection.next()].emplace_back();
                    connection_t& new_connection = reverse_connection_lists[connection.next()].back();
                    new_connection.set_next(last - j);
                    new_connection.set_score(connection.score());
                }
            }
            else {
                reverse_starts.push_back(last - j);
            }
            
            // clear current edges
            subpath_1->clear_next();
            subpath_2->clear_next();
            subpath_1->clear_connection();
            subpath_2->clear_connection();
            
            // reverse complement the paths
            reverse_complement_path_in_place(subpath_1->mutable_path(), node_length);
            reverse_complement_path_in_place(subpath_2->mutable_path(), node_length);
            
            // swap their positions (to maintain topological ordering)
            std::swap(*subpath_1, *subpath_2);
        }
        
        // repeat process for the middle subpath if there is an odd number
        if (multipath_aln->subpath_size() % 2) {
            subpath_t* subpath = multipath_aln->mutable_subpath(subpath_swap_size);
            if (subpath->next_size() > 0 || subpath->connection_size() > 0) {
                for (uint32_t k = 0; k < subpath->next_size(); k++) {
                    reverse_edge_lists[subpath->next(k)].push_back(subpath_swap_size);
                }
                for (uint32_t k = 0; k < subpath->connection_size(); k++) {
                    const connection_t& connection = subpath->connection(k);
                    reverse_connection_lists[connection.next()].emplace_back();
                    connection_t& new_connection = reverse_connection_lists[connection.next()].back();
                    new_connection.set_next(subpath_swap_size);
                    new_connection.set_score(connection.score());
                }
            }
            else {
                reverse_starts.push_back(subpath_swap_size);
            }
            
            subpath->clear_next();
            subpath->clear_connection();
            
            reverse_complement_path_in_place(subpath->mutable_path(), node_length);
        }
        
        // add reversed edges
        for (uint32_t i = 0, j = last; i < multipath_aln->subpath_size(); i++, j--) {
            subpath_t* subpath = multipath_aln->mutable_subpath(i);
            *subpath->mutable_next() = move(reverse_edge_lists[j]);
            *subpath->mutable_connection() = move(reverse_connection_lists[j]);
        }
        
        // if we had starts labeled before, label them again
        if (multipath_aln->start_size() > 0) {
            *multipath_aln->mutable_start() = move(reverse_starts);
        }
    }

    void convert_multipath_alignment_char(multipath_alignment_t& multipath_aln, char from, char to) {
        auto& seq = *multipath_aln.mutable_sequence();
        for (size_t i = 0; i < seq.size(); ++i) {
            if (seq[i] == from) {
                seq[i] = to;
            }
        }
        for (subpath_t& subpath : *multipath_aln.mutable_subpath()) {
            for (path_mapping_t& mapping : *subpath.mutable_path()->mutable_mapping()) {
                for (edit_t& edit : *mapping.mutable_edit()) {
                    if (!edit.sequence().empty()) {
                        auto& eseq = *edit.mutable_sequence();
                        for (size_t i = 0; i < eseq.size(); ++i) {
                            if (eseq[i] == from) {
                                eseq[i] = to;
                            }
                        }
                    }
                }
            }
        }
    }

    void convert_Us_to_Ts(multipath_alignment_t& multipath_aln) {
        convert_multipath_alignment_char(multipath_aln, 'U', 'T');
    }

    void convert_Ts_to_Us(multipath_alignment_t& multipath_aln) {
        convert_multipath_alignment_char(multipath_aln, 'T', 'U');
    }

    void to_proto_multipath_alignment(const multipath_alignment_t& multipath_aln,
                                      MultipathAlignment& proto_multipath_aln_out) {
        proto_multipath_aln_out.clear_subpath();
        proto_multipath_aln_out.clear_start();
        transfer_read_metadata(multipath_aln, proto_multipath_aln_out);
        proto_multipath_aln_out.set_mapping_quality(multipath_aln.mapping_quality());
        for (const auto& subpath : multipath_aln.subpath()) {
            auto subpath_copy = proto_multipath_aln_out.add_subpath();
            subpath_copy->set_score(subpath.score());
            for (auto next : subpath.next()) {
                subpath_copy->add_next(next);
            }
            for (const auto& connection : subpath.connection()) {
                auto connection_copy = subpath_copy->add_connection();
                connection_copy->set_next(connection.next());
                connection_copy->set_score(connection.score());
            }
            if (subpath.has_path()) {
                const auto& path = subpath.path();
                auto path_copy = subpath_copy->mutable_path();
                to_proto_path(path, *path_copy);
            }
        }
        for (auto start : multipath_aln.start()) {
            proto_multipath_aln_out.add_start(start);
        }
    }

    void from_proto_multipath_alignment(const MultipathAlignment& proto_multipath_aln,
                                        multipath_alignment_t& multipath_aln_out) {
        multipath_aln_out.clear_subpath();
        multipath_aln_out.clear_start();
        transfer_read_metadata(proto_multipath_aln, multipath_aln_out);
        multipath_aln_out.set_mapping_quality(proto_multipath_aln.mapping_quality());
        for (auto subpath : proto_multipath_aln.subpath()) {
            auto subpath_copy = multipath_aln_out.add_subpath();
            subpath_copy->set_score(subpath.score());
            for (auto next : subpath.next()) {
                subpath_copy->add_next(next);
            }
            for (const auto& connection : subpath.connection()) {
                auto connection_copy = subpath_copy->add_connection();
                connection_copy->set_next(connection.next());
                connection_copy->set_score(connection.score());
            }
            if (subpath.has_path()) {
                auto path = subpath.path();
                auto path_copy = subpath_copy->mutable_path();
                from_proto_path(path, *path_copy);
            }
        }
        
        for (auto start : proto_multipath_aln.start()) {
            multipath_aln_out.add_start(start);
        }
    }
    
    void to_multipath_alignment(const Alignment& aln, multipath_alignment_t& multipath_aln_out) {
        
        // clear repeated fields
        multipath_aln_out.clear_subpath();
        multipath_aln_out.clear_start();
        
        // transfer read and alignment metadata
        transfer_read_metadata(aln, multipath_aln_out);
        multipath_aln_out.set_mapping_quality(aln.mapping_quality());
        
        // transfer alignment and score
        if (aln.has_path() || aln.score()) {
            subpath_t* subpath = multipath_aln_out.add_subpath();
            subpath->set_score(aln.score());
            from_proto_path(aln.path(), *subpath->mutable_path());
        }
        identify_start_subpaths(multipath_aln_out);
    }

    template<class ProtoAlignment>
    void transfer_from_proto_annotation(const ProtoAlignment& from, multipath_alignment_t& to) {
        for_each_basic_annotation(from,
                                  [&](const string& anno_name) { to.set_annotation(anno_name); },
                                  [&](const string& anno_name, double value) { to.set_annotation(anno_name, value); },
                                  [&](const string& anno_name, bool value) { to.set_annotation(anno_name, value); },
                                  [&](const string& anno_name, const string& value) { to.set_annotation(anno_name, value); });
    }

    template<class ProtoAlignment>
    void transfer_to_proto_annotation(const multipath_alignment_t& from, ProtoAlignment& to) {
        from.for_each_annotation([&](const string& anno_name, multipath_alignment_t::anno_type_t type, const void* value) {
            switch (type) {
                case multipath_alignment_t::Null:
                    break;
                case multipath_alignment_t::Double:
                    set_annotation(to, anno_name, *((const double*) value));
                    break;
                case multipath_alignment_t::Bool:
                    set_annotation(to, anno_name, *((const bool*) value));
                    break;
                case multipath_alignment_t::String:
                    set_annotation(to, anno_name, *((const string*) value));
                    break;
                default:
                    cerr << "error: unrecognized annotation type" << endl;
                    exit(1);
                    break;
            }
        });
    }

    void transfer_read_metadata(const MultipathAlignment& from, multipath_alignment_t& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        transfer_from_proto_annotation(from, to);
    }

    void transfer_read_metadata(const multipath_alignment_t& from, MultipathAlignment& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        transfer_to_proto_annotation(from, to);
    }
    
    void transfer_read_metadata(const multipath_alignment_t& from, multipath_alignment_t& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        from.for_each_annotation([&](const string& anno_name, multipath_alignment_t::anno_type_t type, const void* value) {
            switch (type) {
                case multipath_alignment_t::Null:
                    break;
                case multipath_alignment_t::Double:
                    to.set_annotation(anno_name, *((const double*) value));
                    break;
                case multipath_alignment_t::Bool:
                    to.set_annotation(anno_name, *((const bool*) value));
                    break;
                case multipath_alignment_t::String:
                    to.set_annotation(anno_name, *((const string*) value));
                    break;
                default:
                    cerr << "error: unrecognized annotation type" << endl;
                    exit(1);
                    break;
            }
        });
    }
    
    void transfer_read_metadata(const Alignment& from, multipath_alignment_t& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        transfer_from_proto_annotation(from, to);
    }
    
    void transfer_read_metadata(const multipath_alignment_t& from, Alignment& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        transfer_to_proto_annotation(from, to);
    }

    void transfer_read_metadata(const Alignment& from, Alignment& to) {
        to.set_sequence(from.sequence());
        to.set_quality(from.quality());
        // TODO: do I still care about these fields now that they're taken out
        // of multipath_alignment_t?
        to.set_read_group(from.read_group());
        to.set_name(from.name());
        to.set_sample_name(from.sample_name());
        if (from.has_fragment_prev()) {
            *to.mutable_fragment_prev() = from.fragment_prev();
        }
        if (from.has_fragment_next()) {
            *to.mutable_fragment_next() = from.fragment_next();
        }
        if (from.has_annotation()) {
            *to.mutable_annotation() = from.annotation();
        }
    }

    void transfer_proto_metadata(const Alignment& from, MultipathAlignment& to) {
        // transfer over the fields that are included only in the protobuf object
        to.set_name(from.name());
        to.set_read_group(from.read_group());
        to.set_sample_name(from.sample_name());
        if (from.has_fragment_prev()) {
            to.set_paired_read_name(from.fragment_prev().name());
        }
        else if (from.has_fragment_next()) {
            to.set_paired_read_name(from.fragment_next().name());
        }
    }

    void transfer_proto_metadata(const MultipathAlignment& from, Alignment& to) {
        // transfer over the fields that are included only in the protobuf object
        to.set_name(from.name());
        to.set_read_group(from.read_group());
        to.set_sample_name(from.sample_name());
        
        // not doing paired name because need extra logic to decide if it's prev or next
    }
    
    void merge_non_branching_subpaths(multipath_alignment_t& multipath_aln,
                                      const unordered_set<size_t>* prohibited_merges) {
        
        vector<size_t> in_degree(multipath_aln.subpath_size(), 0);
        vector<bool> has_inward_connection(multipath_aln.subpath_size());
        for (const subpath_t& subpath : multipath_aln.subpath()) {
            for (auto next : subpath.next()) {
                in_degree[next]++;
            }
            for (const auto& connection : subpath.connection()) {
                has_inward_connection[connection.next()] = true;
            }
        }
        
        auto get_mergeable_next = [&](size_t idx) -> int64_t {
            const subpath_t& subpath = multipath_aln.subpath(idx);
            bool prohibited = false;
            if (prohibited_merges) {
                prohibited = prohibited_merges->count(idx);
            }
            if (!prohibited && subpath.next_size() == 1 && subpath.connection_size() == 0
                && in_degree[subpath.next(0)] == 1 && !has_inward_connection[subpath.next(0)]) {
                return subpath.next(0);
            }
            return -1;
        };
        
        vector<bool> removed(multipath_aln.subpath_size(), false);
        
        for (auto i : subpath_topological_order(multipath_aln, false)) {
            
            // this one has been marked for removal,
            if (removed[i]) {
                continue;
            }
            
            // the subpath we might merge into
            subpath_t* subpath = multipath_aln.mutable_subpath(i);
            
            int64_t last = -1;
            // iterate through non-branching subpaths
            for (int64_t j = get_mergeable_next(i); j >= 0; j = get_mergeable_next(j)) {
                
                // mark the next one for removal
                removed[j] = true;
                                
                subpath_t* merge_subpath = multipath_aln.mutable_subpath(j);
                
                subpath->set_score(subpath->score() + merge_subpath->score());
                
                path_t* merge_path = merge_subpath->mutable_path();
                if (merge_path->mapping_size() == 0) {
                    continue;
                }
                
                path_t* path = subpath->mutable_path();
                path_mapping_t* final_mapping = path->mutable_mapping(path->mapping_size() - 1);
                const position_t& final_position = final_mapping->position();
                
                path_mapping_t* first_mapping = merge_path->mutable_mapping(0);
                const position_t& first_position = first_mapping->position();
                
                int64_t mapping_idx = 0;
                
                // do we need to merge the abutting mappings?
                if (first_position.node_id() == final_position.node_id() &&
                    first_position.is_reverse() == final_position.is_reverse() &&
                    first_position.offset() == final_position.offset() + mapping_from_length(*final_mapping)) {
                    // do we need to merge the abutting edits?
                    int64_t edit_idx = 0;
                    if (final_mapping->edit_size() && first_mapping->edit_size()) {
                        edit_t* final_edit = final_mapping->mutable_edit(final_mapping->edit_size() - 1);
                        const edit_t& first_edit = first_mapping->edit(0);
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
                    for (; edit_idx < first_mapping->edit_size(); edit_idx++) {
                        *final_mapping->add_edit() = move(*first_mapping->mutable_edit(edit_idx));
                    }
                    
                    mapping_idx++;
                }
                
                for (; mapping_idx < merge_path->mapping_size(); mapping_idx++) {
                    *path->add_mapping() = move(*merge_path->mutable_mapping(mapping_idx));
                }
                
                last = j;
            }
            
            // move the adjacencies over from the last one we merged in
            if (last >= 0) {
                subpath->clear_next();
                subpath->clear_connection();
                for (int64_t next : multipath_aln.subpath(last).next()) {
                    subpath->add_next(next);
                }
                for (const auto& connection : multipath_aln.subpath(last).connection()) {
                    *subpath->add_connection() = connection;
                }
            }
        }
        
        // go back and do the removals
        vector<size_t> removed_so_far(multipath_aln.subpath_size(), 0);
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            if (i > 0) {
                removed_so_far[i] = removed_so_far[i - 1];
            }
            
            if (removed[i]) {
                // this one has been marked for removal
                removed_so_far[i]++;
                continue;
            }
            
            if (removed_so_far[i]) {
                // move it up in the vector past the removed subpaths
                *multipath_aln.mutable_subpath(i - removed_so_far[i]) = move(*multipath_aln.mutable_subpath(i));
            }
        }
        
        // did we merge and remove any subpaths?
        if (!removed_so_far.empty() && removed_so_far.back()) {
            // trim the vector of subpaths
            multipath_aln.mutable_subpath()->resize(multipath_aln.subpath_size() - removed_so_far.back());
            
            // update the indexes of the adjacencies
            for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
                subpath_t* subpath = multipath_aln.mutable_subpath(i);
                for (size_t j = 0; j < subpath->next_size(); j++) {
                    subpath->set_next(j, subpath->next(j) - removed_so_far[subpath->next(j)]);
                }
                
                for (size_t j = 0; j < subpath->connection_size(); j++) {
                    auto connection = subpath->mutable_connection(j);
                    connection->set_next(connection->next() - removed_so_far[connection->next()]);
                }
            }
            
            // update the indexes of the starts
            for (size_t i = 0; i < multipath_aln.start_size(); ++i) {
                multipath_aln.set_start(i, multipath_aln.start(i) - removed_so_far[multipath_aln.start(i)]);
            }
        }
    }
    
    vector<vector<int64_t>> connected_components(const multipath_alignment_t& multipath_aln) {
        
        int64_t comps = 0;
        
        vector<vector<int64_t>> reverse_edge_lists(multipath_aln.subpath_size());
        
        for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const subpath_t& subpath = multipath_aln.subpath(i);
            // collect edges and connections by their target
            for (size_t j = 0; j < subpath.next_size(); j++) {
                reverse_edge_lists[subpath.next(j)].push_back(i);
            }
            for (size_t j = 0; j < subpath.connection_size(); j++) {
                reverse_edge_lists[subpath.connection(j).next()].push_back(i);
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
                
                const subpath_t& subpath = multipath_aln.subpath(at);
                for (int64_t j = 0; j < subpath.next_size(); j++) {
                    int64_t idx = subpath.next(j);
                    if (!collected[idx]) {
                        collected[idx] = true;
                        stack.push_back(idx);
                    }
                }
                for (int64_t j = 0; j < subpath.connection_size(); j++) {
                    int64_t idx = subpath.connection(j).next();
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
    
    void extract_sub_multipath_alignment(const multipath_alignment_t& multipath_aln,
                                         const vector<int64_t>& subpath_indexes,
                                         multipath_alignment_t& sub_multipath_aln) {

        transfer_read_metadata(multipath_aln, sub_multipath_aln);
        
        // create subpaths for each of the ones we're retaining and record the translation
        unordered_map<int64_t, int64_t> new_index;
        for (int64_t i = 0; i < subpath_indexes.size(); i++) {
            int64_t old_idx = subpath_indexes[i];
            const subpath_t& old_subpath = multipath_aln.subpath(old_idx);
            
            subpath_t* subpath = sub_multipath_aln.add_subpath();
            *subpath->mutable_path() = old_subpath.path();
            subpath->set_score(old_subpath.score());
            
            new_index[old_idx] = i;
        }
        
        // add edges according to the translation
        for (int64_t i = 0; i < subpath_indexes.size(); i++) {
            const subpath_t& old_subpath = multipath_aln.subpath(subpath_indexes[i]);
            subpath_t* new_subpath = sub_multipath_aln.mutable_subpath(i);
            for (int64_t j = 0; j < old_subpath.next_size(); j++) {
                if (new_index.count(old_subpath.next(j))) {
                    new_subpath->add_next(new_index[old_subpath.next(j)]);
                }
            }
            for (int64_t j = 0; j < old_subpath.connection_size(); j++) {
                const connection_t& old_connection = old_subpath.connection(j);
                if (new_index.count(old_connection.next())) {
                    connection_t* new_connection = new_subpath->add_connection();
                    new_connection->set_next(new_index[old_connection.next()]);
                    new_connection->set_score(old_connection.score());
                }
            }
        }
        
        // assume that if we had starts labeled before, we want them again
        if (multipath_aln.start_size() > 0) {
            identify_start_subpaths(sub_multipath_aln);
        }
    }

    void append_multipath_alignment(multipath_alignment_t& multipath_aln,
                                    const multipath_alignment_t& to_append) {
        
        size_t original_size = multipath_aln.subpath().size();
        for (const subpath_t& appending_subpath : to_append.subpath()) {
            subpath_t* new_subpath = multipath_aln.add_subpath();
            new_subpath->set_score(appending_subpath.score());
            *new_subpath->mutable_path() = appending_subpath.path();
            new_subpath->mutable_next()->reserve(appending_subpath.next_size());
            for (auto n : appending_subpath.next()) {
                new_subpath->add_next(n + original_size);
            }
            for (const auto& c : appending_subpath.connection()) {
                auto new_connection = new_subpath->add_connection();
                new_connection->set_next(c.next() + original_size);
                new_connection->set_score(c.score());
            }
        }
        if (multipath_aln.start_size() != 0 && to_append.start_size() != 0) {
            for (auto s : to_append.start()) {
                multipath_aln.add_start(s + original_size);
            }
        }
        else if (multipath_aln.start_size() != 0) {
            identify_start_subpaths(multipath_aln);
        }
    }

    bool contains_connection(const multipath_alignment_t& multipath_aln) {
        bool no_connection = true;
        for (size_t i = 0; i < multipath_aln.subpath_size() && no_connection; ++i) {
            no_connection = multipath_aln.subpath(i).connection().empty();
        }
        return !no_connection;
    }

    vector<tuple<int64_t, int64_t, int64_t, int64_t>>
    search_multipath_alignment(const multipath_alignment_t& multipath_aln,
                               const pos_t& graph_pos, int64_t seq_pos) {
        
#ifdef debug_search
        cerr << "starting search for " << graph_pos << " at seq pos " << seq_pos << endl;
#endif
        
        vector<tuple<int64_t, int64_t, int64_t, int64_t>> return_val;
        
        vector<int64_t> subpath_seq_pos(multipath_aln.subpath_size(), 0);
        
        for (int64_t i = 0; i < multipath_aln.subpath_size(); ++i) {
#ifdef debug_search
            cerr << "subpath " << i << endl;
#endif
            const subpath_t& subpath = multipath_aln.subpath(i);
            const path_t& path = subpath.path();
            int64_t to_length_here = subpath_seq_pos[i];
            int64_t to_length_thru = to_length_here + path_to_length(subpath.path());
            // tell the next subpaths where they start
            for (auto j : subpath.next()) {
                subpath_seq_pos[j] = to_length_thru;
            }
            for (const auto& connection : subpath.connection()) {
                subpath_seq_pos[connection.next()] = to_length_thru;
            }
            
            if (to_length_here <= seq_pos && to_length_thru >= seq_pos) {
                // this is where we might expect to find the sequence position
                
#ifdef debug_search
                cerr << "interval " << to_length_here << " " << to_length_thru << " covers seq pos " << seq_pos << endl;
#endif
                
                for (int64_t j = 0; j < path.mapping_size(); ++j) {
                    
                    const auto& mapping = path.mapping(j);
                    const auto& pos = mapping.position();
                    
#ifdef debug_search
                    cerr << "mapping " << j << " at graph pos " << debug_string(pos) << endl;
#endif
                    
                    if (pos.node_id() == id(graph_pos) && pos.is_reverse() == is_rev(graph_pos)) {
                        int64_t offset_here = pos.offset();
                        
#ifdef debug_search
                        cerr << "position " << debug_string(pos) << " consistent with graph pos " << graph_pos << endl;
#endif
                        
                        // this mapping is on the right node to be a match
                        
                        for (int64_t k = 0; k < mapping.edit_size(); ++k) {
                            
                            const auto& edit = mapping.edit(k);
                            int64_t to_length_thru_edit = to_length_here + edit.to_length();
                            int64_t offset_thru_edit = offset_here + edit.from_length();
                            
#ifdef debug_search
                            cerr << "edit " << k << ", to length interval " << to_length_here << " " << to_length_thru_edit << ", offset interval " << offset_here << " " << offset_thru_edit << endl;
#endif
                            
                            // does this edit contain both the sequence and graph positions (allowing
                            // for a past-the-last position on the final edit)?
                            if (to_length_here <= seq_pos &&
                                (to_length_thru_edit > seq_pos || (to_length_thru_edit == seq_pos &&
                                                                   (k + 1 == mapping.edit_size() ||
                                                                    to_length_here == to_length_thru_edit))) &&
                                offset_here <= offset(graph_pos) &&
                                (offset_thru_edit > offset(graph_pos) || (offset_thru_edit == offset(graph_pos) &&
                                                                          (k + 1 == mapping.edit_size() ||
                                                                           offset_here == offset_thru_edit)))) {
                                
                                // are the offsets within the edit consistent with each other?
                                int64_t graph_l = offset(graph_pos) - offset_here;
                                int64_t seq_l = seq_pos - to_length_here;
                                bool consistent = (graph_l == seq_l ||
                                                   (graph_l == 0 && edit.from_length() == 0) ||
                                                   (seq_l == 0 && edit.to_length() == 0));
                                
#ifdef debug_search
                                cerr << "read interval " << to_length_here << " " << to_length_thru_edit << " covers seq pos " << seq_pos << ", offset interval " << offset_here << " " << offset_thru_edit << " covers offset "  << offset(graph_pos) << ", consistent? " << consistent << endl;
#endif
                                
                                // handle some special cases of the past-the-last position to make canonical results
                                // TODO: ugly
                                bool must_place_here = true;
                                if (consistent && k + 1 == mapping.edit_size() && to_length_thru_edit == seq_pos
                                    && offset_thru_edit == offset(graph_pos)) {
                                    // we're looking at locating this position at the past-the-last position on
                                    // a mapping, but it might also exist at the first position of the next mapping.
                                    // if so, we will canonicalize it to go there instead.
                                    
#ifdef debug_search
                                    cerr << "checking if must place past-the-last" << endl;
#endif
                                    if (j + 1 < path.mapping_size()) {
                                        // the next mapping is still on this subpath
                                        const auto& next_pos = path.mapping(j + 1).position();
                                        must_place_here &= (next_pos.node_id() != pos.node_id()
                                                            || next_pos.is_reverse() != pos.is_reverse()
                                                            || next_pos.offset() != offset_thru_edit);
                                    }
                                    else {
                                        // we have to check the next subpaths
                                        for (auto n : subpath.next()) {
                                            const auto& next_pos = multipath_aln.subpath(n).path().mapping(0).position();
                                            must_place_here &= (next_pos.node_id() != pos.node_id()
                                                                || next_pos.is_reverse() != pos.is_reverse()
                                                                || next_pos.offset() != offset_thru_edit);
                                        }
                                    }
                                }
                                
                                if (consistent && must_place_here) {
                                    // winner winner chicken dinner, record the match
                                    int64_t l = max(seq_l, graph_l);
#ifdef debug_search
                                    cerr << "recording match " << i << " " << j << " " << k << " " << l << endl;
#endif
                                    
                                    return_val.emplace_back(i, j, k, l);
                                }
                                
                            }
                            offset_here = offset_thru_edit;
                            to_length_here = to_length_thru_edit;
                        }
                    }
                    else {
                        to_length_here += mapping_to_length(mapping);
                    }
                }
            }
        }
        return return_val;
    }

    pair<tuple<int64_t, int64_t, int64_t>, vector<tuple<int64_t, int64_t, int64_t, int64_t>>>
    trace_path(const multipath_alignment_t& multipath_aln, const Path& path,
               int64_t subpath_idx, int64_t mapping_idx, int64_t edit_idx, int64_t base_idx,
               bool search_left, int64_t search_limit) {
        
#ifdef debug_trace
        cerr << "entering trace path algorithm, searching left? " << search_left << endl;
        cerr << "tracing path " << pb2json(path) << endl;
        cerr << "start coordinate: " << subpath_idx << ", " << mapping_idx << ", " << edit_idx << ", " << base_idx << endl;
        cerr << "search limit " << search_limit << endl;
#endif
        pair<tuple<int64_t, int64_t, int64_t>, vector<tuple<int64_t, int64_t, int64_t, int64_t>>> return_val;
        auto& pfarthest = return_val.first;
        auto& mfarthest = return_val.second;
        
        if (search_left) {
            // we like to index the base as if it's from the left even though it's from the right
            // to simplify some conditions later, so we have to reverse it now
            const auto& start_path = multipath_aln.subpath(subpath_idx).path();
            if (mapping_idx < start_path.mapping_size()) {
                const auto& start_mapping = start_path.mapping(mapping_idx);
                if (edit_idx < start_mapping.edit_size()) {
                    const auto& start_edit = start_mapping.edit(edit_idx);
                    base_idx = max(start_edit.from_length(), start_edit.to_length()) - base_idx;
#ifdef debug_trace
                    cerr << "flip leftward base index on edit " << debug_string(start_edit) << " to " << base_idx << endl;
#endif
                }
            }
        }
        
        // the farthest match along the path
        pfarthest = search_left ? tuple<int64_t, int64_t, int64_t>(path.mapping_size(), 0, 0)
                                : tuple<int64_t, int64_t, int64_t>(-1, 0, 0);
        
        // the position on the mp aln that corresponds to this match
        
        mfarthest.emplace_back(subpath_idx, mapping_idx, edit_idx, base_idx);
        
        // DFS stack
        vector<bool> stacked(multipath_aln.subpath_size(), false);
        vector<tuple<int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t>> stack;
        
        int64_t incr = search_left ? -1 : 1;
        
        // which end of the path are we starting on?
        int64_t p_start_mapping_idx = 0, p_start_edit_idx = 0;
        if (search_left) {
            p_start_mapping_idx = path.mapping_size() - 1;
            if (p_start_mapping_idx >= 0) {
                const auto& start_mapping = path.mapping(p_start_mapping_idx);
                p_start_edit_idx = start_mapping.edit_size() - 1;
            }
        }
        
        // we may need reverse adjacencies if we're searching leftwards
        vector<vector<int64_t>> reverse_adjacencies;
        if (search_left) {
            reverse_adjacencies.resize(multipath_aln.subpath_size());
            for (int64_t i = 0; i < multipath_aln.subpath_size(); ++i) {
                const auto& subpath = multipath_aln.subpath(i);
                for (auto n : subpath.next()) {
                    reverse_adjacencies[n].push_back(i);
                }
                for (const auto& c : subpath.connection()) {
                    reverse_adjacencies[c.next()].push_back(i);
                }
            }
        }
        
        // start at the indicated location on the mp aln and the beginning of the path
        // note: the logic is simpler if we treat base indexes from 0 regardless of which
        // direction w
        stack.emplace_back(subpath_idx, mapping_idx, edit_idx, base_idx,
                           p_start_mapping_idx, p_start_edit_idx, 0);
        stacked[subpath_idx] = true;
        bool first_iter = true;
        while (!stack.empty()) {
            // load up the indexes of where we're going to look for a match
            int64_t i, j, k, l, pj, pk, pl;
            tie(i, j, k, l, pj, pk, pl) = stack.back();
            stack.pop_back();
            
            
            // the indexes of the last non-empty match for each index
            int64_t ni, nj, nk, nl, npj, npk, npl;
            tie(ni, nj, nk, nl, npj, npk, npl) = tie(i, j, k, l, pj, pk, pl);
            bool any_new_matches = first_iter;
            first_iter = false;
                        
#ifdef debug_trace
            cerr << "destack (" << i << " " << j << " " << k << " " << l << ") (" << pj << " " << pk << " " << pl << ")" << endl;
#endif
            
            const subpath_t& subpath = multipath_aln.subpath(i);
            const path_t& mpath = subpath.path();
            bool reached_mismatch = false;
            while (j < mpath.mapping_size() && pj < path.mapping_size() && j >= 0 && pj >= 0 &&
                   (search_left || pj < search_limit) && (!search_left || pj >= search_limit) && !reached_mismatch) {
                
#ifdef debug_trace
                cerr << "mp mapping " << j << ", p mapping " << pj << endl;
#endif
                                
                const auto& mmapping = mpath.mapping(j);
                const auto& pmapping = path.mapping(pj);
                
                // TODO: these mappings can actually have positions that are inconsistent
                // but checking for consistency is tricky at mapping boundaries that
                // also correspond to node boundaries
                
                // skip over the subpath mapping if it's empty
                bool mnonempty = false;
                for (int64_t m = k, n = l; m < mmapping.edit_size() && m >= 0 && !mnonempty; m += incr) {
                    const auto& edit = mmapping.edit(m);
                    int64_t rem = max(edit.from_length(), edit.to_length()) - n;
                    if (rem) {
                        mnonempty = true;
                    }
                    n = 0;
                }
                if (!mnonempty) {
#ifdef debug_trace
                    cerr << "mp mapping is empty" << endl;
#endif
                    l = 0;
                    j += incr;
                    if (search_left && j >= 0) {
                        k = mpath.mapping(j).edit_size() - 1;
                    }
                    else {
                        k = 0;
                    }
                    continue;
                }
                // skip over the path mapping if it's empty
                bool pnonempty = false;
                for (int64_t m = pk, n = pl; m < pmapping.edit_size() && m >= 0 && !pnonempty; m += incr) {
                    const auto& edit = pmapping.edit(m);
                    int64_t rem = max(edit.from_length(), edit.to_length()) - n;
                    if (rem) {
                        pnonempty = true;
                    }
                    n = 0;
                }
                if (!pnonempty) {
#ifdef debug_trace
                    cerr << "p mapping is empty" << endl;
#endif
                    pl = 0;
                    pj += incr;
                    if (search_left && pj >= 0) {
                        pk = path.mapping(pj).edit_size() - 1;
                    }
                    else {
                        pk = 0;
                    }
                    continue;
                }
                
                // now that we know we're looking at non-empty mappings, the positions need to match
                
                // TODO: it would be nice if we didn't need to iterate over the entire mapping
                
                // find the graph position on the subpath
                const auto& mpos = mmapping.position();
                int64_t moffset = mpos.offset();
                for (int64_t m = 0; m < k; ++m) {
                    moffset += mmapping.edit(m).from_length();
                }
                if (search_left) {
                    moffset += mmapping.edit(k).from_length();
                }
                if (l > 0 && mmapping.edit(k).from_length() > 0) {
                    moffset += l * incr;
                }
                
                // find the graph position on the path
                const auto& ppos = pmapping.position();
                int64_t poffset = ppos.offset();
                for (int64_t m = 0; m < pk; ++m) {
                    poffset += pmapping.edit(m).from_length();
                }
                if (search_left) {
                    poffset += pmapping.edit(pk).from_length();
                }
                if (pl > 0 && pmapping.edit(pk).from_length() > 0) {
                    poffset += pl * incr;
                }
                
#ifdef debug_trace
                cerr << "mp pos " << mpos.node_id() << " " << mpos.is_reverse() << " " << moffset << ", p pos " << ppos.node_id() << " " << ppos.is_reverse() << " " << poffset << endl;
#endif
                
                if (mpos.node_id() != ppos.node_id() || mpos.is_reverse() != ppos.is_reverse()
                    || moffset != poffset) {
                    // these positions don't match
                    reached_mismatch = true;
                }
                
                
                // try to match edits
                while (k < mmapping.edit_size() && k >= 0 && pk < pmapping.edit_size() && pk >= 0 && !reached_mismatch) {
#ifdef debug_trace
                    cerr << "mp edit " << k << " " << l << ", p edit " << pk << " " << pl << endl;
#endif
                    const auto& medit = mmapping.edit(k);
                    const auto& pedit = pmapping.edit(pk);
                    if (medit.from_length() == 0 && medit.to_length() == 0) {
                        // skip over an empty edit
                        l = 0;
                        k += incr;
#ifdef debug_trace
                        cerr << "mp edit empty" << endl;
#endif
                    }
                    else if (pedit.from_length() == 0 && pedit.to_length() == 0) {
                        // skip over an empty edit
                        pl = 0;
                        pk += incr;
#ifdef debug_trace
                        cerr << "p edit empty" << endl;
#endif
                    }
                    else if ((medit.from_length() == 0) == (pedit.from_length() == 0) &&
                             (medit.to_length() == 0) == (pedit.to_length() == 0) &&
                             medit.sequence().empty() == pedit.sequence().empty()) {
                        
                        // the type of edit matches
                        
                        if ((medit.from_length() == 0 || medit.from_length() - l == pedit.from_length() - pl) &&
                            (medit.to_length() == 0 || medit.to_length() - l == pedit.to_length() - pl)) {
                            // the size of edit matches
                            l = 0;
                            k += incr;
                            pl = 0;
                            pk += incr;
                            
#ifdef debug_trace
                            cerr << "edits match" << endl;
#endif
                        }
                        else if ((medit.from_length() == 0 || medit.from_length() - l < pedit.from_length() - pl) &&
                                 (medit.to_length() == 0 || medit.to_length() - l < pedit.to_length() - pl)) {
                            
                            // subpath edit is a prefix of path edit
                            pl += max(medit.from_length(), medit.to_length()) - l;
                            if (pl == max(pedit.from_length(), pedit.to_length())) {
                                // TODO: won't this never happen because of the earlier condition?
                                pl = 0;
                                pk += incr;
                            }
                            l = 0;
                            k += incr;
#ifdef debug_trace
                            cerr << "mp edit is prefix" << endl;
#endif
                        }
                        else {
                            // path edit is a prefix of subpath edit
                            l += max(pedit.from_length(), pedit.to_length()) - pl;
                            if (l == max(medit.from_length(), medit.to_length())) {
                                // TODO: won't this never happen because of the earlier condition?
                                l = 0;
                                k += incr;
                            }
                            pl = 0;
                            pk += incr;
#ifdef debug_trace
                            cerr << "p edit is prefix" << endl;
#endif
                        }
                        // we made a non-empty match, update the non-empty index trackers
                        tie(ni, nj, nk, nl, npj, npk, npl) = tie(i, j, k, l, pj, pk, pl);
                        any_new_matches = true;
                    }
                    else {
                        // the edits do not match
                        reached_mismatch = true;
#ifdef debug_trace
                        cerr << "edits mismatch" << endl;
#endif
                    }
                }
                
                // did we finish off either mapping?
                if (k == mmapping.edit_size() || k < 0) {
#ifdef debug_trace
                    cerr << "finished mp mapping" << endl;
#endif
                    j += incr;
                    k = 0;
                    if (search_left && j >= 0) {
                        k = mpath.mapping(j).edit_size() - 1;
                    }
                    l = 0;
                }
                if (pk == pmapping.edit_size() || pk < 0) {
#ifdef debug_trace
                    cerr << "finished p mapping" << endl;
#endif
                    pj += incr;
                    pk = 0;
                    if (search_left && pj >= 0) {
                        pk = path.mapping(pj).edit_size() - 1;
                    }
                    pl = 0;
                }
            }
            // how far did we get along the path by walking this subpath (looking at non-empty matches only)?
            if (any_new_matches) {
                if ((search_left && (npj < get<0>(pfarthest) ||
                                     (npj == get<0>(pfarthest) && npk < get<1>(pfarthest)) ||
                                     (npj == get<0>(pfarthest) && npk == get<1>(pfarthest) && npl < get<2>(pfarthest)))) ||
                    (!search_left && (npj > get<0>(pfarthest) ||
                                      (npj == get<0>(pfarthest) && npk > get<1>(pfarthest)) ||
                                      (npj == get<0>(pfarthest) && npk == get<1>(pfarthest) && npl > get<2>(pfarthest))))) {
                    // we've traversed more of the path than on any previous subpath
#ifdef debug_trace
                    cerr << "new farthest at subpath index " << ni << ", " << nj << ", " << nk << ", " << nl << " and path index " << npj << ", " << npk << ", " << npl << endl;
#endif
                    pfarthest = make_tuple(npj, npk, npl);
                    mfarthest.clear();
                    mfarthest.emplace_back(ni, nj, nk, nl);
                }
                else if (npj == get<0>(pfarthest) && npk == get<1>(pfarthest) && npl == get<2>(pfarthest)) {
                    // we've tied the farthest we've gone along the path previously
#ifdef debug_trace
                    cerr << "tied existing farthest at subpath index " << ni << ", " << nj << ", " << nk << ", " << nl << endl;
#endif
                    mfarthest.emplace_back(ni, nj, nk, nl);
                }
            }
            
            if (pj == path.mapping_size() || pj < 0) {
                // we've found the farthest point possible
                break;
            }
            if ((j == mpath.mapping_size() || j < 0) && !reached_mismatch) {
                // we got to the end of the subpath without exhausting our match
                if (search_left) {
                    for (auto n : reverse_adjacencies[i]) {
                        if (!stacked[n]) {
                            stacked[n] = true;
                            const auto& next_path = multipath_aln.subpath(n).path();
                            int64_t next_mapping_idx = next_path.mapping_size() - 1;
                            int64_t next_edit_idx = 0;
                            if (next_mapping_idx >= 0) {
                                const auto& next_mapping = next_path.mapping(next_mapping_idx);
                                next_edit_idx = next_mapping.edit_size() - 1;
                            }
                            stack.emplace_back(n, next_mapping_idx, next_edit_idx, 0, pj, pk, pl);
#ifdef debug_trace
                            cerr << "stack up (" << n << " " << next_mapping_idx << " " << next_edit_idx << " " << 0 << ") (" << pj << " " << pk << " " << pl << ")" << endl;
#endif
                        }
                    }
                }
                else {
                    for (auto n : subpath.next()) {
                        if (!stacked[n]) {
                            stacked[n] = true;
                            stack.emplace_back(n, 0, 0, 0, pj, pk, pl);
#ifdef debug_trace
                            cerr << "stack up (" << n << " " << 0 << " " << 0 << " " << 0 << ") (" << pj << " " << pk << " " << pl << ")" << endl;
#endif
                        }
                    }
                    for (const auto& c : subpath.connection()) {
                        if (!stacked[c.next()]) {
                            stacked[c.next()] = true;
                            stack.emplace_back(c.next(), 0, 0, 0, pj, pk, pl);
#ifdef debug_trace
                            cerr << "stack up (" << c.next() << " " << 0 << " " << 0 << " " << 0 << ") (" << pj << " " << pk << " " << pl << ")" << endl;
#endif
                        }
                    }
                }
            }
        }
        
        if (search_left) {
            // we're set up to find past-the-first coordinates, but if we're going leftward
            // what we want is actually the final coordinates
            
#ifdef debug_trace
            cerr << "converting past-the-first coordinates to at-the-first" << endl;
            cerr << "p farthest (" << get<0>(pfarthest) << " " << get<1>(pfarthest) << " " << get<2>(pfarthest) << ") -> ";
#endif
            if (get<0>(pfarthest) < 0) {
                get<0>(pfarthest) = 0;
            }
            else if (get<1>(pfarthest) < 0) {
                get<1>(pfarthest) = 0;
                // even though we interpreted the base index of 0 differently before, it's already what we want here
            }
            else {
                // switch index of the base to from-the-end
                const auto& final_edit = path.mapping(get<0>(pfarthest)).edit(get<1>(pfarthest));
                get<2>(pfarthest) = max(final_edit.from_length(), final_edit.to_length()) - get<2>(pfarthest);
            }
#ifdef debug_trace
            cerr << "(" << get<0>(pfarthest) << " " << get<1>(pfarthest) << " " << get<2>(pfarthest) << ")" << endl;
#endif
            
            for (auto& coord : mfarthest) {
#ifdef debug_trace
                cerr << "m farthest (" << get<0>(coord) << " " << get<1>(coord) << " " << get<2>(coord) << " " << get<3>(coord) << ") -> ";
#endif
                if (get<1>(coord) < 0) {
                    get<1>(coord) = 0;
                }
                else if (get<2>(coord) < 0) {
                    get<2>(coord) = 0;
                    // even though we interpreted the base index of 0 differently before, it's already what we want here
                }
                else {
                    // switch index of the base to from-the-end
                    const auto& final_edit =  multipath_aln.subpath(get<0>(coord)).path().mapping(get<1>(coord)).edit(get<2>(coord));
                    get<3>(coord) = max(final_edit.from_length(), final_edit.to_length()) - get<3>(coord);
                }
#ifdef debug_trace
                cerr << "(" << get<0>(coord) << " " << get<1>(coord) << " " << get<2>(coord) << " " << get<3>(coord) << ")" << endl;
#endif
            }
        }
        
        return return_val;
    }

    bool contains_match(const multipath_alignment_t& multipath_aln, const pos_t& pos,
                        int64_t read_pos, int64_t match_length) {
        
#ifdef debug_find_match
        cerr << "starting search for match at graph pos " << pos << ", read pos " << read_pos << ", length " << match_length << endl;
        cerr << debug_string(multipath_aln) << endl;
#endif
        
        // to keep track of the read interval corresponding to each subpath
        vector<int64_t> to_length(multipath_aln.subpath_size(), 0);
        for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
            const subpath_t& subpath = multipath_aln.subpath(i);
            const path_t& path = subpath.path();
            
            int64_t to_length_here = to_length[i];
            
#ifdef debug_find_match
            cerr << "at subpath " << i << " at read pos " << to_length_here << endl;
#endif
            
            for (size_t j = 0; j < path.mapping_size(); ++j) {
                
                const path_mapping_t& mapping = path.mapping(j);
                const position_t& position = mapping.position();
                
#ifdef debug_find_match
                cerr << "at mapping " << j << ", graph pos " << make_pos_t(position) << ", read pos " << to_length_here << endl;
#endif
                
                if (id(pos) == position.node_id() && is_rev(pos) == position.is_reverse()) {
                    // we're on the node with the position we're looking for
                    int64_t offset_here = position.offset();
                    for (size_t k = 0; k < mapping.edit_size(); ++k) {
                        const edit_t& edit = mapping.edit(k);
#ifdef debug_find_match
                        cerr << "at edit " << k << ", offset " << offset_here << ", read pos " << to_length_here << endl;
#endif
                        if (offset_here <= offset(pos) && offset_here + edit.from_length() > offset(pos) &&
                            to_length_here <= read_pos && to_length_here + edit.to_length() > read_pos
                            && (offset(pos) - offset_here) == (read_pos - to_length_here)
                            && edit.sequence().empty()) {
                            // we're going to pass a place where we might find a match
                            
                            
#ifdef debug_find_match
                            cerr << "starting a DFS for the match on edit " << k << ", offset " << offset_here << ", read pos " << to_length_here << ", initial walked " << offset(pos) - offset_here << endl;
#endif
                            
                            // do DFS to find the match
                            // records of (remaining, subpath idx, mapping idx, edit idx, which next)
                            vector<tuple<int64_t, size_t, size_t, size_t>> stack;
                            stack.emplace_back(offset(pos) - offset_here, i, j, k);
                            while (!stack.empty()) {
                                
                                int64_t walked;
                                size_t di, dj, dk;
                                tie(walked, di, dj, dk) = stack.back();
                                stack.pop_back();

                                
                                const subpath_t& subpath_here = multipath_aln.subpath(di);
                                const path_t& path_here = subpath_here.path();
                                const path_mapping_t& mapping_here = path_here.mapping(dj);
                                const edit_t& edit_here = mapping_here.edit(dk);
                                
#ifdef debug_find_match
                                cerr << "DFS location " << di << ", " << dj << ", " << dk << ", walked " << walked << ", edit " << debug_string(edit_here) <<  endl;
#endif
                                
                                if (edit_here.to_length() && edit_here.from_length() && edit_here.sequence().empty()) {
                                    // this edit can continue the match
                                    walked += edit_here.to_length();
                                    if (walked >= match_length) {
                                        // we found the match
                                        return true;
                                    }
                                    // advance by one edit
                                    ++dk;
                                    if (dk == mapping_here.edit_size()) {
                                        // we're to the next mapping
                                        dk = 0;
                                        ++dj;
                                    }
                                    if (dj == path_here.mapping_size()) {
                                        // we're at the boundary to the next subpaths
                                        for (auto n : subpath_here.next()) {
                                            stack.emplace_back(walked, n, 0, 0);
#ifdef debug_find_match
                                            cerr << "queue up next subpath " << n << ", walked " << walked << endl;
#endif
                                        }
                                    }
                                    else {
                                        stack.emplace_back(walked, di, dj, dk);
#ifdef debug_find_match
                                        cerr << "queue up next edit " << di << ", " << dj << ", " << dk << ", walked " << walked << endl;
#endif
                                    }
                                }
                            }
                            
                        }
                        offset_here += edit.from_length();
                        to_length_here += edit.to_length();
                    }
                }
                else {
                    // we're not on the node we want, just record the read interval of this mapping
                    to_length_here += mapping_to_length(mapping);
                }
            }
            
            // record the read interval of the successors
            for (auto n : subpath.next()) {
                to_length[n] = to_length_here;
            }
            for (const auto& connection : subpath.connection()) {
                to_length[connection.next()] = to_length_here;
            }
        }
        
        // we never found the match
        return false;
    }

    // TODO: does it really make sense to split this algorithm into a separate surject
    // and CIGAR conversion? it seems like i'm replicating a lot of the same work
    vector<pair<int, char>> cigar_against_path(const multipath_alignment_t& multipath_aln, const string& path_name,
                                               bool rev, int64_t path_pos, const PathPositionHandleGraph& graph,
                                               int64_t min_splice_length) {
        
#ifdef debug_cigar
        cerr << "converting mp aln to CIGAR on path " << path_name << ", rev? " << rev << ", pos " << path_pos << endl;
        cerr << debug_string(multipath_aln) << endl;
#endif
        vector<pair<int, char>> cigar;
        if (path_pos < 0) {
            // read is unmapped
            return cigar;
        }
        
        // a graph of runs of alignments to the path
        // records of (subpath index, mapping index, num_mappings, final step, from connection, adj list of (index, distance))
        vector<tuple<size_t, size_t, size_t, step_handle_t, bool, vector<pair<size_t, size_t>>>> run_graph;
        
        path_handle_t path_handle = graph.get_path_handle(path_name);
        
        // the runs from the previous node that could be extended in the current iteration
        // second value in the pair is the number of bases remaining until the end of node
        unordered_map<step_handle_t, pair<size_t, size_t>> curr_runs;
        
        size_t num_mappings = 0;
        for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
            
            const auto& subpath = multipath_aln.subpath(i);
            const auto& path = subpath.path();
            
            if (subpath.next_size() + subpath.connection_size() > 1) {
                cerr << "error: cannot convert a multipath alignment to a CIGAR unless is consists of a single non-branching path" << endl;
            }
            
            for (size_t j = 0; j < path.mapping_size(); ++j, ++num_mappings) {
                
                const auto& mapping = path.mapping(j);
                const auto& pos = mapping.position();
                
#ifdef debug_cigar
                cerr << "look for runs on mapping " << i << " " << j << ": " << debug_string(mapping) << endl;
#endif
                
                if (i == 0 && j == 0 && !rev) {
                    // base case if the position is given for the beginning of the surjected alignment
                    // TODO: will this work if the position is between two nodes?
                    auto step = graph.get_step_at_position(path_handle, path_pos);
#ifdef debug_cigar
                    cerr << "got step " << graph.get_id(graph.get_handle_of_step(step)) << " " << graph.get_is_reverse(graph.get_handle_of_step(step)) << " as path pos " << graph.get_position_of_step(step) << endl;
#endif
                    if (graph.get_id(graph.get_handle_of_step(step)) != pos.node_id()
                        || graph.get_is_reverse(graph.get_handle_of_step(step)) != pos.is_reverse()
                        || path_pos - graph.get_position_of_step(step) != pos.offset()) {
                        // the step doesn't match our starting position, but this can sometimes happen when
                        // a position occurs right at a node boundary
                        auto prev_step = graph.get_previous_step(step);
                        
#ifdef debug_cigar
                        cerr << "didn't match, walk back to " << graph.get_id(graph.get_handle_of_step(prev_step)) << " " << graph.get_is_reverse(graph.get_handle_of_step(prev_step)) << " as path pos " << graph.get_position_of_step(prev_step) << endl;
#endif
                        if (prev_step != graph.path_front_end(path_handle)
                            && graph.get_id(graph.get_handle_of_step(prev_step)) == pos.node_id()
                            && graph.get_is_reverse(graph.get_handle_of_step(prev_step)) == pos.is_reverse()
                            && graph.get_length(graph.get_handle_of_step(prev_step)) == pos.offset()) {
                            step = prev_step;
                        }
                        else {
                            cerr << "error: couldn't find a matching subpath on " << path_name << " for read " << multipath_aln.sequence() << endl;
                            exit(1);
                        }
                    }
                    run_graph.emplace_back(0, 0, 1, step, false, vector<pair<size_t, size_t>>());
                    curr_runs[step] = pair<size_t, size_t>(0, graph.get_length(graph.get_handle_of_step(step))
                                                           - mapping_from_length(mapping) - pos.offset());
                    
#ifdef debug_cigar
                    cerr << "initializing on forward strand with step on " << graph.get_id(graph.get_handle_of_step(step)) << " " << graph.get_is_reverse(graph.get_handle_of_step(step)) << " at pos " << graph.get_position_of_step(step) << endl;
#endif
                    continue;
                }
                
                auto from_length = mapping_from_length(mapping);
                
                // get the next mapping's position, if there is one
                const position_t* next_pos = nullptr;
                if (j + 1 < path.mapping_size()) {
                    next_pos = &path.mapping(j + 1).position();
                }
                else if (i + 1 < multipath_aln.subpath_size()) {
                    next_pos = &multipath_aln.subpath(i + 1).path().mapping().front().position();
                }
                
                if (next_pos && next_pos->node_id() == pos.node_id() && next_pos->is_reverse() == pos.is_reverse()
                    && next_pos->offset() == pos.offset() + from_length) {
                    // we only care about transitions that are between nodes, this one is within a node
                    
#ifdef debug_cigar
                    cerr << "transition for " << i << " " << j << " does not exit node, stalling current runs" << endl;
#endif
                    
                    // but keep track of the number of mapping in each run
                    for (pair<const step_handle_t, pair<size_t, size_t>>& curr_run : curr_runs) {
                        ++get<2>(run_graph[curr_run.second.first]);
                        curr_run.second.second -= from_length;
                    }
                    continue;
                }
                
                // remember where previously created run nodes stop
                size_t num_runs_before_extend = run_graph.size();
                
                // is this step across a splice connection?
                bool across_connection = i > 0 ? j == 0 && !multipath_aln.subpath(i - 1).connection().empty() : false;
                
                // check whether steps on this handle extend previous runs or not
                unordered_map<step_handle_t, pair<size_t, size_t>> next_runs;
                next_runs.reserve(curr_runs.size());
                handle_t handle = graph.get_handle(pos.node_id(), pos.is_reverse());
#ifdef debug_cigar
                cerr << "iterating over steps on " << graph.get_id(handle) << " " << graph.get_is_reverse(handle) << endl;
                cerr << "curr runs to extend:" << endl;
                for (const auto& run : curr_runs) {
                    auto h = graph.get_handle_of_step(run.first);
                    cerr << "\trun " << run.second.first << ", node " << graph.get_id(h) << " " << graph.get_is_reverse(h) << ", rem " << run.second.second << ", pos " << graph.get_position_of_step(run.first) << endl;
                }
#endif
                graph.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
                    if (graph.get_path_handle_of_step(step) != path_handle ||
                        (graph.get_handle_of_step(step) != handle) != rev) {
                        // we're only concerned about one strand of this one path
                        return;
                    }
                    step_handle_t prev = rev ? graph.get_next_step(step) : graph.get_previous_step(step);

                    size_t remaining = (graph.get_length(graph.get_handle_of_step(step))
                                        - from_length - pos.offset());
                    auto it = curr_runs.find(prev);
                    if (it != curr_runs.end() && !across_connection && it->second.second == 0 &&
                        pos.offset() == 0) {
                        // this is the next step we would expect along a previous run, and it's
                        // not across a connection
                        next_runs[step] = make_pair(it->second.first, remaining);
                        auto& run_node = run_graph[it->second.first];
#ifdef debug_cigar
                        cerr << "extending run " << it->second.first << " from step at " << graph.get_position_of_step(get<3>(run_node)) << " to step at " << graph.get_position_of_step(step) << endl;
#endif
                        
                        ++get<2>(run_node);
                        get<3>(run_node) = step;
                        curr_runs.erase(it);
                    }
                    else {
                        // we're at the start of a new run, or we just crossed a connection, start
                        // a new run
                        next_runs[step] = make_pair(run_graph.size(), remaining);
                        run_graph.emplace_back(i, j, 1, step, across_connection, vector<pair<size_t, size_t>>());
#ifdef debug_cigar
                        cerr << "new run " << run_graph.size() - 1 << " for step at " << graph.get_position_of_step(step) << " and subpath indexes " << i << " " << j << endl;
#endif
                    }
                });
                
                // TODO: are there situations where i would need to split a run into multiple chunks
                // in order to find the full length?
                
                // check if any of the unextended runs can make a long-distance adjacency
                // to the new runs
                for (const auto& curr_run : curr_runs) {
                    // an unextended run from the previous iteration
                    for (size_t run_idx = num_runs_before_extend; run_idx < run_graph.size(); ++run_idx) {
                        // fresh new run that we just found
                        auto& run_node = run_graph[run_idx];
                        
                        int64_t dist;
                        if (rev) {
                            dist = (graph.get_position_of_step(curr_run.first)
                                    - graph.get_position_of_step(get<3>(run_node))
                                    + mapping.position().offset() + curr_run.second.second
                                    - graph.get_length(graph.get_handle_of_step(get<3>(run_node))));
                        }
                        else {
                            dist = (graph.get_position_of_step(get<3>(run_node))
                                    - graph.get_position_of_step(curr_run.first)
                                    + mapping.position().offset() + curr_run.second.second
                                    - graph.get_length(graph.get_handle_of_step(curr_run.first)));
                        }
                        if (dist >= 0) {
                            // they are in increasing order (relative to the strand)

                            // add an edge
                            get<5>(run_graph[curr_run.second.first]).emplace_back(run_idx, dist);
#ifdef debug_cigar
                            cerr << "adjacency of length " << dist << " from " << curr_run.second.first << " to " << run_idx << endl;
#endif
                        }
                    }
                }
                
                curr_runs = move(next_runs);
            }
        }
        
        // okay, now we did that whole business, it's time to find the path through the run graph
        // that corresponds to the surjected alignment
        
#ifdef debug_cigar
        cerr << "doing mapping length DP with a total number of mappings " << num_mappings << endl;
#endif
        
        // find the longest path, measured by number of mappings (should consume the whole alignment)
        vector<size_t> mapping_dp(run_graph.size(), 0);
        vector<bool> full_length(mapping_dp.size(), false);
        for (size_t i = 0; i < mapping_dp.size(); ++i) {
            
            const auto& run_node = run_graph[i];
            mapping_dp[i] += get<2>(run_node);
            for (const auto& edge : get<5>(run_node)) {
                mapping_dp[edge.first] = max(mapping_dp[i], mapping_dp[edge.first]);
            }
#ifdef debug_cigar
            cerr << "\t" << i << ": " << mapping_dp[i] << endl;
#endif
            
            // does it complete a full traversal of the alignment?
            full_length[i] = (mapping_dp[i] == num_mappings);
        }
        
#ifdef debug_cigar
        cerr << "identifying full length run combinations" << endl;
#endif
        
        // identify the run nodes that could be part of a full length alignment
        for (int64_t i = full_length.size() - 1; i >= 0; --i) {
            for (const auto& edge : get<5>(run_graph[i])) {
                full_length[i] = full_length[i] || full_length[edge.first];
            }
#ifdef debug_cigar
            cerr << "\t" << i << ": " << full_length[i] << endl;
#endif
        }
        
#ifdef debug_cigar
        cerr << "doing path distance DP" << endl;
#endif
        
        // identify the shortest path through the run graph based on total path length
        vector<size_t> path_dist_dp(run_graph.size(), numeric_limits<size_t>::max());
        vector<pair<int64_t, int64_t>> backpointer(path_dist_dp.size(), make_pair(-1, -1));
        int64_t best = -1;
        for (int64_t i = 0; i < path_dist_dp.size(); ++i) {
            if (!full_length[i]) {
                // this node can't be part of a full length alignment so we don't want
                // to find paths through it
                continue;
            }
            const auto& run_node = run_graph[i];
            if (get<0>(run_node) == 0 && get<1>(run_node) == 0) {
                // base case
                path_dist_dp[i] = 0;
            }
            for (const auto& edge : get<5>(run_node)) {
                size_t dist_thru = path_dist_dp[i] + edge.second;
                if (dist_thru < path_dist_dp[edge.first]) {
                    backpointer[edge.first] = pair<int64_t, int64_t>(i, edge.second);
                    path_dist_dp[edge.first] = dist_thru;
                }
            }
            if (mapping_dp[i] == num_mappings) {
#ifdef debug_cigar
                cerr << "\tpossible end " << i << ": is full length? " << full_length[i] << ", path dist " << path_dist_dp[i] << endl;
#endif
                if (rev) {
                    const auto& final_mapping = multipath_aln.subpath().back().path().mapping().back();
                    int64_t path_pos_here = (graph.get_position_of_step(get<3>(run_node))
                                             + graph.get_length(graph.get_handle_of_step(get<3>(run_node)))
                                             - final_mapping.position().offset()
                                             - mapping_from_length(final_mapping));
                    if (path_pos_here != path_pos) {
                        // this run doesn't end where it should based on our path position, it can't
                        // be part of the CIGAR
#ifdef debug_cigar
                        cerr << "path pos doesn't match expected " << path_pos << ", instead got " << path_pos_here << " from step pos " << graph.get_position_of_step(get<3>(run_node)) << ", node length " << graph.get_length(graph.get_handle_of_step(get<3>(run_node))) << " mapping offset " << final_mapping.position().offset() << " and mapping len " << mapping_from_length(final_mapping) << endl;
#endif
                        continue;
                    }
                }
                
                if (best == -1 || path_dist_dp[i] < path_dist_dp[best]) {
                    // this is the shortest full length run sequence we've seen
                    best = i;
                }
            }
        }
        
        if (best == -1) {
            cerr << "error: couldn't find a matching subpath on " << path_name << " for read " << multipath_aln.sequence() << endl;
            exit(1);
        }
        
#ifdef debug_cigar
        cerr << "backtracing" << endl;
#endif
        
        // compute the traceback
        vector<int64_t> traceback(1, best);
        while (backpointer[traceback.back()].first != -1) {
            traceback.push_back(backpointer[traceback.back()].first);
        }
        
#ifdef debug_cigar
        cerr << "forming CIGAR string" << endl;
#endif
        
        // now finally make the CIGAR
        for (int64_t i = traceback.size() - 1; i >= 0; --i) {
            
#ifdef debug_cigar
            cerr << "forward trace to " << traceback[i] << endl;
#endif
            
            auto& run_node = run_graph[traceback[i]];
            
            // handle the edge between these runs
            if (i != traceback.size() - 1) {
                int64_t dist = backpointer[traceback[i]].second;
#ifdef debug_cigar
                cerr << "handling adjacency of length " << dist << endl;
#endif
                if (dist >= min_splice_length || get<4>(run_node)) {
                    // let this be a splice either because it's long or because it's across a connection
                    cigar.emplace_back(dist, 'N');
                }
                else if (!cigar.empty() && cigar.back().second == 'D') {
                    cigar.back().first += dist;
                }
                else {
                    cigar.emplace_back(dist, 'D');
                }
            }
            
            // determine the bounds of the iteration over the mp aln that corresponds
            // to this run
            size_t j = get<0>(run_node);
            size_t k = get<1>(run_node);
            size_t j_end, k_end;
            if (i > 0) {
                auto& next_run_node = run_graph[traceback[i - 1]];
                j_end = get<0>(next_run_node);
                k_end = get<1>(next_run_node);
            }
            else {
                j_end = multipath_aln.subpath_size();
                k_end = 0;
            }
            
#ifdef debug_cigar
            cerr << "iteration bounds are " << j << " " << k << " to " << j_end << " " << k_end << endl;
#endif
            
            // convert this segment of the mp aln into CIGAR
            while (j != j_end || k != k_end) {
                const auto& path = multipath_aln.subpath(j).path();
                const auto& mapping = path.mapping(k);
#ifdef debug_cigar
                cerr << "on mapping " << debug_string(mapping) << endl;
#endif
                for (const auto& edit : mapping.edit()) {
                    char cigar_code;
                    int length;
                    if (edit.from_length() == edit.to_length()) {
                        cigar_code = 'M';
                        length = edit.from_length();
                    }
                    else if (edit.from_length() > 0 && edit.to_length() == 0) {
                        cigar_code = 'D';
                        length = edit.from_length();
                    }
                    else if (edit.to_length() > 0 && edit.from_length() == 0) {
                        cigar_code = 'I';
                        length = edit.to_length();
                    }
                    else {
                        throw std::runtime_error("Spliced CIGAR construction can only convert simple edits");
                    }
                    
                    if (!cigar.empty() && cigar_code == cigar.back().second) {
                        cigar.back().first += length;
                    }
                    else {
                        cigar.emplace_back(length, cigar_code);
                    }
                }
                k++;
                if (k == path.mapping_size()) {
                    ++j;
                    k = 0;
                }
            }
        }
        
        // change start/end insertions into softclips
        if (!cigar.empty()) {
            if (cigar.front().second == 'I') {
                cigar.front().second = 'S';
            }
            if (cigar.back().second == 'I') {
                cigar.back().second = 'S';
            }
        }
        
        if (rev) {
            // return the cigar relative to the forward strand
            for (size_t i = 0, end = cigar.size() / 2; i < end; ++i) {
                swap(cigar[i], cigar[cigar.size() - i - 1]);
            }
        }
        
#ifdef debug_cigar
        cerr << "got cigar: ";
        for (auto& cigar_record : cigar) {
            cerr << cigar_record.first << cigar_record.second;
        }
        cerr << endl;
        cerr << "coalescing runs of I/D..." << endl;
#endif

        consolidate_ID_runs(cigar);
        
#ifdef debug_cigar
        cerr << "final cigar: ";
        for (auto& cigar_record : cigar) {
            cerr << cigar_record.first << cigar_record.second;
        }
        cerr << endl;
#endif
        
        return cigar;
    }
    
    bool validate_multipath_alignment(const multipath_alignment_t& multipath_aln, const HandleGraph& handle_graph) {
        
        // are the subpaths in topological order?
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const subpath_t& subpath = multipath_aln.subpath(i);
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
                const subpath_t& subpath = multipath_aln.subpath(i);
                for (size_t j = 0; j < subpath.next_size(); j++) {
                    is_source[subpath.next(j)] = false;
                }
                for (const auto& connection : subpath.connection()) {
                    is_source[connection.next()] = false;
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
            
            const subpath_t& subpath = multipath_aln.subpath(i);
            int64_t subsequence_length = path_to_length(subpath.path());
            subpath_read_interval[i].second = subpath_read_interval[i].first + subsequence_length;
            
            if (subpath.next_size() == 0 && subpath.connection_size() == 0) {
                if (subpath_read_interval[i].second != multipath_aln.sequence().size()) {
#ifdef debug_verbose_validation
                    cerr << "validation failure on using complete read" << endl;
                    cerr << "subpath " <<  i << " ends on sequence index " << subpath_read_interval[i].second << " of " << multipath_aln.sequence().size() << endl;
                    for (size_t j = 0; j < multipath_aln.subpath_size(); j++) {
                        cerr << j << " (" << subpath_read_interval[j].first << ", " << subpath_read_interval[j].second << "): ";
                        for (size_t k = 0; k < multipath_aln.subpath(j).next_size(); k++) {
                            cerr << multipath_aln.subpath(j).next(k) << " ";
                        }
                        for (auto& connection : subpath.connection()) {
                            cerr << connection.next() << " ";
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
                            cerr << "validation failure on read contiguity from subpath " << i << " with read interval " << subpath_read_interval[i].first << ":" << subpath_read_interval[i].second << " to next " << subpath.next(j) << " with read interval " << subpath_read_interval[subpath.next(j)].first << ":" << subpath_read_interval[subpath.next(j)].second << endl;
#endif
                            return false;
                        }
                    }
                    else {
                        subpath_read_interval[subpath.next(j)].first = subpath_read_interval[i].second;
                    }
                }
                for (const auto& connection : subpath.connection()) {
                    if (subpath_read_interval[connection.next()].first >= 0) {
                        if (subpath_read_interval[connection.next()].first != subpath_read_interval[i].second) {
#ifdef debug_verbose_validation
                            cerr << "validation failure on read contiguity from subpath " << i << " with read interval " << subpath_read_interval[i].first << ":" << subpath_read_interval[i].second << " to connection " << connection.next() << " with read interval " << subpath_read_interval[connection.next()].first << ":" << subpath_read_interval[connection.next()].second << endl;
#endif
                            return false;
                        }
                    }
                    else {
                        subpath_read_interval[connection.next()].first = subpath_read_interval[i].second;
                    }
                }
            }
        }
        
        // are all of the subpaths nonempty?
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            if (multipath_aln.subpath(i).path().mapping_size() == 0) {
#ifdef debug_verbose_validation
                cerr << "validation failure on containing only nonempty paths" << endl;
                cerr << "subpath " << i << ": " << debug_string(multipath_aln.subpath(i)) << endl;
#endif
                return false;
            }
            for (size_t j = 0; j < multipath_aln.subpath(i).path().mapping_size(); j++) {
                if (multipath_aln.subpath(i).path().mapping(j).edit_size() == 0) {
#ifdef debug_verbose_validation
                    cerr << "validation failure on containing only nonempty mappings" << endl;
                    cerr << "subpath " << i << ": " << debug_string(multipath_aln.subpath(i)) << endl;
#endif
                    return false;
                }
            }
        }
        
        
        // are the subpaths contiguous within the graph?
        
        auto validate_adjacent_mappings = [&](const path_mapping_t& mapping_from, const path_mapping_t& mapping_to) {
            size_t mapping_from_end_offset = mapping_from.position().offset() + mapping_from_length(mapping_from);
            
            handle_t handle_from = handle_graph.get_handle(mapping_from.position().node_id(), mapping_from.position().is_reverse());
            handle_t handle_to = handle_graph.get_handle(mapping_to.position().node_id(), mapping_to.position().is_reverse());
            
            
            
            if (handle_from == handle_to) {
                if (!(mapping_to.position().offset() == 0 && mapping_from_end_offset == handle_graph.get_length(handle_from))) {
                    // We aren't going from the end of the handle back to its start (over an edge)
                
                    if (mapping_to.position().offset() != mapping_from_end_offset) {
                        // So then the mappings need to abut and they don't.
#ifdef debug_verbose_validation
                        cerr << "validation failure on within-node adjacency" << endl;
                        cerr << debug_string(mapping_from) << "->" << debug_string(mapping_to) << endl;
#endif
                        return false;
                    } else {
                        // No edge involved. We can succeed early.
                        return true;
                    }
                }
            }
            
            // If we get here, we must be crossing an edge.
            
            if (mapping_from_end_offset != handle_graph.get_length(handle_graph.get_handle(mapping_from.position().node_id()))) {
#ifdef debug_verbose_validation
                cerr << "validation failure on using edge at middle of node" << endl;
                cerr << debug_string(mapping_from) << "->" << debug_string(mapping_to) << endl;
#endif
                return false;
            }
            
            
            
            bool found_edge = false;
            function<bool(const handle_t&)> check_for_edge = [&](const handle_t& next_handle) {
                found_edge = (next_handle == handle_to);
                return !found_edge;
            };
            handle_graph.follow_edges(handle_from, false, check_for_edge);
            
            if (!found_edge) {
#ifdef debug_verbose_validation
                cerr << "validation failure on nodes not connected by an edge" << endl;
                cerr << debug_string(mapping_from) << "->" << debug_string(mapping_to) << endl;
#endif
                return false;
            }
            return true;
        };
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const subpath_t& subpath = multipath_aln.subpath(i);
            const path_t& path = subpath.path();
            for (size_t j = 1; j < path.mapping_size(); j++) {
                if (!validate_adjacent_mappings(path.mapping(j - 1), path.mapping(j))) {
                    return false;
                }
            }
            const path_mapping_t& final_mapping = path.mapping(path.mapping_size() - 1);
            for (size_t j = 0; j < subpath.next_size(); j++) {
                if (!validate_adjacent_mappings(final_mapping, multipath_aln.subpath(subpath.next(j)).path().mapping(0))) {
                    return false;
                }
            }
            // connections are not required to be contiguous, ignore them here
        }
        
        
        // do the paths represent valid alignments of the associated read string and graph path?
        
        auto validate_mapping_edits = [&](const path_mapping_t& mapping, const string& subseq) {
            handle_t handle = handle_graph.get_handle(mapping.position().node_id(), mapping.position().is_reverse());
            size_t node_idx = mapping.position().offset();
            size_t seq_idx = 0;
            for (size_t i = 0; i < mapping.edit_size(); i++) {
                const edit_t& edit = mapping.edit(i);
                if (edit.to_length() == edit.from_length() && edit.sequence().empty()) {
                    for (size_t j = 0; j < edit.from_length(); j++, node_idx++, seq_idx++) {
                        if (handle_graph.get_base(handle, node_idx) != subseq[seq_idx]) {
#ifdef debug_verbose_validation
                            cerr << "validation failure on match that does not match on node " << handle_graph.get_id(handle) << (handle_graph.get_is_reverse(handle) ? "-" : "+") << endl;
                            cerr << "node sequence: " << handle_graph.get_sequence(handle) << ", offset: " << node_idx << endl;
                            cerr << "read subsequence: " << subseq << ", offset: " << seq_idx << endl;
                            
                            cerr << debug_string(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                    }
                }
                else if (edit.to_length() == edit.from_length() && !edit.sequence().empty()) {
                    if (edit.sequence().size() != edit.to_length()) {
#ifdef debug_verbose_validation
                        cerr << "validation failure on mismatched sequence length and to length: " << debug_string(edit) << " in mapping " << debug_string(mapping) << endl;
#endif
                        return false;
                    }
                    
                    bool is_Ns = find_if(edit.sequence().begin(), edit.sequence().end(), [](char c) {return c != 'N';}) == edit.sequence().end();
                    for (size_t j = 0; j < edit.from_length(); j++, node_idx++, seq_idx++) {
                        // we will also let N's be marked as mismatches even if the node sequence is also Ns
                        if (handle_graph.get_base(handle, node_idx) == subseq[seq_idx] && !is_Ns) {
#ifdef debug_verbose_validation
                            cerr << "validation failure on mismatch that matches" << endl;
                            cerr << debug_string(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                        if (edit.sequence()[j] != subseq[seq_idx]) {
#ifdef debug_verbose_validation
                            cerr << "validation failure on substitution sequence that does not match read" << endl;
                            cerr << debug_string(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                    }
                }
                else if (edit.to_length() > 0 && edit.from_length() == 0) {
                    if (edit.sequence().size() != edit.to_length()) {
#ifdef debug_verbose_validation
                        cerr << "validation failure on mismatched sequence length and to length: " << debug_string(edit) << " on mapping " << debug_string(mapping) << endl;
#endif
                        return false;
                    }
                    
                    for (size_t j = 0; j < edit.to_length(); j++, seq_idx++) {
                        
                        if (edit.sequence()[j] != subseq[seq_idx]) {
#ifdef debug_verbose_validation
                            cerr << "validation failure on insertion sequence that does not match read" << endl;
                            cerr << debug_string(mapping) << ", " << subseq << endl;
#endif
                            return false;
                        }
                    }
                }
                else if (edit.from_length() > 0 && edit.to_length() == 0) {
                    node_idx += edit.from_length();
                }
                else {
#ifdef debug_verbose_validation
                    cerr << "validation failure on non-simple edit " << debug_string(edit) << endl;
#endif
                    return false;
                }
            }
            return true;
        };
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const subpath_t& subpath = multipath_aln.subpath(i);
            const path_t& path = subpath.path();
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
        //            const subpath_t& subpath = multipath_aln.subpath(i);
        //            Alignment& alignment;
        //            *alignment.mutable_sequence() = multipath_aln.sequence().substr(subpath_read_interval[i].first,
        //                                                                            subpath_read_interval[i].second - subpath_read_interval[i].first);
        //            *alignment.mutable_quality() = multipath_aln.quality().substr(subpath_read_interval[i].first,
        //                                                                          subpath_read_interval[i].second - subpath_read_interval[i].first);
        //            *alignment.mutable_path() = subpath.path();
        //        }
        
        
        return true;
    }
    
    void view_multipath_alignment(ostream& out, const multipath_alignment_t& multipath_aln, const HandleGraph& handle_graph) {
        
        size_t max_line_length = 128;
        
        vector<pair<int64_t, int64_t>> subpath_read_interval(multipath_aln.subpath_size(), pair<int64_t, int64_t>(0, 0));
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            
            const subpath_t& subpath = multipath_aln.subpath(i);
            subpath_read_interval[i].second = subpath_read_interval[i].first + path_to_length(subpath.path());
            
            for (int64_t j : subpath.next()) {
                subpath_read_interval[j].first = subpath_read_interval[i].second;
            }
        }
        
        auto format_position = [](const position_t& pos) {
            stringstream strm;
            strm << pos.node_id() << (pos.is_reverse() ? "-" : "+") << (pos.offset() ? (":" + to_string(pos.offset())) : "");
            return strm.str();
        };
        
        for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
            const subpath_t& subpath = multipath_aln.subpath(i);
            
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
                const path_mapping_t& mapping = subpath.path().mapping(j);
                string pos_string = format_position(mapping.position());
                
                
                stringstream mapping_read_strm;
                stringstream mapping_node_strm;
                mapping_read_strm << pos_string << " ";
                mapping_node_strm << string(pos_string.size() + 1, ' ');
                
                string node_seq = handle_graph.get_sequence(handle_graph.get_handle(mapping.position().node_id(),
                                                                                    mapping.position().is_reverse()));
                
                int64_t node_at = mapping.position().offset();
                for (const edit_t& edit : mapping.edit()) {
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
    
    void view_multipath_alignment_as_dot(ostream& out, const multipath_alignment_t& multipath_aln, bool show_graph) {
        out << "digraph graphname {" << endl;
        out << "rankdir=\"LR\";" << endl;
        
        // Track graph nodes so we get one node for each
        unordered_set<id_t> mentioned_nodes;
        // Similarly for graph edges
        unordered_set<pair<id_t, id_t>> mentioned_edges;
        
        // Do the start node
        out << "start[label=\"Start\" shape=circle];" << endl;
        for (size_t start_subpath : multipath_aln.start()) {
            // Hook it up to each start subpath
            out << "start -> s" << start_subpath << ";" << endl;
        }
        
        for (size_t i = 0; i < multipath_aln.subpath_size(); i++) {
            // For each subpath, say it with its score
            out << "s" << i << " [label=\"" << i << "\" shape=circle tooltip=\"" << multipath_aln.subpath(i).score() << "\"];" << endl;
            
            for (size_t next_subpath : multipath_aln.subpath(i).next()) {
                // For each edge from it, say where it goes
                out << "s" << i << " -> s" << next_subpath << ";" << endl;
            }
            
            if (show_graph) {
                auto& path = multipath_aln.subpath(i).path();
                for (size_t j = 0; j < path.mapping_size(); j++) {
                    // For each mapping in the path, show the vg node in the graph too
                    auto node_id = path.mapping(j).position().node_id();
                    
                    if (!mentioned_nodes.count(node_id)) {
                        // This graph node eneds to be made
                        mentioned_nodes.insert(node_id);
                        out << "g" << node_id << " [label=\"" << node_id << "\" shape=box];" << endl;
                    }
                    
                    // Attach the subpath to each involved graph node.
                    out << "s" << i << " -> g" << node_id << " [dir=none color=blue];" << endl;
                    
                    if (j != 0) {
                        // We have a previous node in this segment of path. What is it?
                        auto prev_id = path.mapping(j-1).position().node_id();
                        pair<id_t, id_t> edge_pair{prev_id, node_id};
                        
                        if (!mentioned_edges.count(edge_pair)) {
                            // This graph edge needs to be made
                            mentioned_edges.insert(edge_pair);
                            
                            out << "g" << prev_id << " -> g" << node_id << ";" << endl;
                        }
                    }
                }
            }
        }
        
        out << "}" << endl;
    }

    string debug_string(const connection_t& connection) {
        string to_return = "{next: " + to_string(connection.next()) + ", score: " + to_string(connection.score()) + "}";
        return to_return;
    }
    
    string debug_string(const subpath_t& subpath) {
        string to_return = "{path: " + debug_string(subpath.path());
        if (!subpath.next().empty()) {
            to_return += ", next: [";
            for (size_t i = 0; i < subpath.next_size(); ++i) {
                if (i > 0) {
                    to_return += ", ";
                }
                to_return += to_string(subpath.next(i));
            }
            to_return += "]";
        }
        if (!subpath.connection().empty()) {
            to_return += ", connection: [";
            for (size_t i = 0; i < subpath.connection_size(); ++i) {
                if (i > 0) {
                    to_return += ", ";
                }
                to_return += debug_string(subpath.connection(i));
            }
            to_return += "]";
        }
        to_return += ", score: " + to_string(subpath.score());
        to_return += "}";
        return to_return;
    }

    string debug_string(const multipath_alignment_t& multipath_aln) {
        string to_return = "{seq: " + multipath_aln.sequence();
        if (!multipath_aln.quality().empty()) {
            to_return += ", qual: " + string_quality_short_to_char(multipath_aln.quality());
        }
        if (!multipath_aln.subpath().empty()) {
            to_return += ", subpath: [";
            for (size_t i = 0; i < multipath_aln.subpath_size(); ++i) {
                if (i > 0) {
                    to_return += ", ";
                }
                to_return += debug_string(multipath_aln.subpath(i));
            }
            to_return += "]";
        }
        to_return += ", mapq: " + to_string(multipath_aln.mapping_quality());
        if (!multipath_aln.start().empty()) {
            to_return += ", start: [";
            for (size_t i = 0; i < multipath_aln.start_size(); ++i) {
                if (i > 0) {
                    to_return += ", ";
                }
                to_return += to_string(multipath_aln.start(i));
            }
            to_return += "]";
        }
        to_return += "}";
        return to_return;
    }

}







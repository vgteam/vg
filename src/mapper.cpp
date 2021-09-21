#include <unordered_set>
#include "mapper.hpp"
#include "graph.hpp"
#include "annotation.hpp"
#include "statistics.hpp"
#include "path.hpp"
#include "entropy.hpp"
#include "alignment.hpp"
#include "translator.hpp"
#include "algorithms/subgraph.hpp"
#include "algorithms/nearest_offsets_in_paths.hpp"
#include "algorithms/jump_along_path.hpp"
#include "algorithms/approx_path_distance.hpp"
#include "algorithms/path_string.hpp"
#include "algorithms/alignment_path_offsets.hpp"
#include "algorithms/extract_containing_graph.hpp"

//#define debug_mapper

//#define debug_strip_match

using namespace vg::io;

namespace vg {

// init the static memo
thread_local vector<size_t> BaseMapper::adaptive_reseed_length_memo;

BaseMapper::BaseMapper(PathPositionHandleGraph* xidex,
                       gcsa::GCSA* g,
                       gcsa::LCPArray* a,
                       haplo::ScoreProvider* haplo_score_provider) :
      AlignerClient(estimate_gc_content(g))
      , xindex(xidex)
      , gcsa(g)
      , lcp(a)
      , haplo_score_provider(haplo_score_provider)
      , min_mem_length(1)
      , mem_reseed_length(0)
      , fast_reseed(true)
      , fast_reseed_length_diff(0.45)
      , adaptive_reseed_diff(true)
      , adaptive_diff_exponent(0.065)
      , hit_max(0)
      , mapping_quality_method(Approx)
      , max_mapping_quality(60)
      , strip_bonuses(false)
      , assume_acyclic(false)
      , exclude_unaligned(false)
{

    if (xindex != nullptr) {
        avg_node_length = xindex->get_total_length() / xindex->get_node_count();
    }
    
    // TODO: removing these consistency checks because we seem to have violated them pretty wontonly in
    // the code base already by changing the members directly when they were still public
    
//    // allow a default constructor with no indexes
//    if (xidex || g || a) {
//        if(xidex == nullptr) {
//            // We need an PathPositionHandleGraph graph.
//            cerr << "error:[vg::Mapper] cannot create an xg-based Mapper with null xg index" << endl;
//            exit(1);
//        }
//        
//        if (g == nullptr) {
//            // We need a GCSA2 index too.
//            cerr << "error:[vg::Mapper] cannot create an xg-based Mapper with null GCSA2 index" << endl;
//            exit(1);
//        }
//        
//        if (a == nullptr) {
//            // And an LCP array
//            cerr << "error:[vg::Mapper] cannot create an xg-based Mapper with null LCP index" << endl;
//            exit(1);
//        }
//    }
}
   
BaseMapper::BaseMapper(void) : BaseMapper(nullptr, nullptr, nullptr) {
    // Nothing to do. Default constructed and can't really do anything.
}

// Use the GCSA2 index to find super-maximal exact matches.
vector<MaximalExactMatch>
BaseMapper::find_mems_simple(string::const_iterator seq_begin,
                             string::const_iterator seq_end,
                             int max_mem_length,
                             int min_mem_length,
                             int reseed_length) {
    
    if (!gcsa) {
        cerr << "error:[vg::Mapper] a GCSA2 index is required to query MEMs" << endl;
        exit(1);
    }
    
    vector<MaximalExactMatch> mems;
    
    // an empty sequence matches the entire bwt
    if (seq_begin == seq_end) {
        mems.emplace_back(MaximalExactMatch(seq_begin, seq_end,
                                            gcsa::range_type(0, gcsa->size() - 1)));
        return mems;
    }
    
    // find SMEMs using GCSA+LCP array
    // algorithm sketch:
    // set up a cursor pointing to the last position in the sequence
    // set up a structure to track our MEMs, and set it == "" and full range match
    // while our cursor is >= the beginning of the string
    //   try a step of backwards searching using LF mapping
    //   if our range goes to 0
    //       go back to the last non-empty range
    //       emit the MEM corresponding to this range
    //       start a new mem
    //           use the LCP array's parent function to cut off the end of the match
    //           (effectively, this steps up the suffix tree)
    //           and calculate the new end point using the LCP of the parent node
    // emit the final MEM, if we finished in a matching state
    
    // the temporary MEM we'll build up in this process
    auto full_range = gcsa::range_type(0, gcsa->size() - 1);
    string::const_iterator curr_end = seq_end;
    string::const_iterator cursor = seq_end - 1;
    // start off looking at the last character in the query
    gcsa::range_type range = accelerate_mem_query(seq_begin, cursor);
    while (cursor >= seq_begin) {
        // hold onto our previous range
        auto last_range = range;
        // execute one step of LF mapping
        range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
        if (gcsa::Range::empty(range)
            || (max_mem_length && curr_end - cursor > max_mem_length)
            || curr_end - cursor > gcsa->order()) {
            // break on N; which for DNA we assume is non-informative
            // this *will* match many places in assemblies; this isn't helpful
            if (*cursor == 'N' || last_range == full_range) {
                // we mismatched in a single character
                // there is no MEM here
                mems.emplace_back(cursor + 1, curr_end, last_range);
                curr_end = cursor;
                range = full_range;
                --cursor;
            } else {
                // we've exhausted our BWT range, so the last match range was maximal
                // or: we have exceeded the order of the graph (FPs if we go further)
                //     we have run over our parameter-defined MEM limit
                // record the last MEM
                mems.emplace_back(cursor + 1, curr_end, last_range);
                // set up the next MEM using the parent node range
                // length of last MEM, which we use to update our end pointer for the next MEM
                int64_t last_mem_length = (curr_end - cursor) - 1;
                // get the parent suffix tree node corresponding to the parent of the last MEM's STNode
                gcsa::STNode parent = lcp->parent(last_range);
                // change the end for the next mem to reflect our step size
                size_t step_size = last_mem_length - parent.lcp();
                curr_end = curr_end - step_size;
                // and set up the next MEM using the parent node range
                range = parent.range();
            }
        } else {
            // just step to the next position
            --cursor;
        }
    }
    // if we have a non-empty MEM at the end, record it
    if (curr_end > seq_begin) {
        mems.emplace_back(seq_begin, curr_end, range);
    }
    
    // find the SMEMs from the mostly-SMEM and some MEM list we've built
    // FIXME: un-hack this (it shouldn't be needed!)
    // the algorithm sometimes generates MEMs contained in SMEMs
    // with the pattern that they have the same beginning position
    map<string::const_iterator, string::const_iterator> smems_begin;
    for (auto& mem : mems) {
        auto x = smems_begin.find(mem.begin);
        if (x == smems_begin.end()) {
            smems_begin[mem.begin] = mem.end;
        } else {
            if (x->second < mem.end) {
                x->second = mem.end;
            }
        }
    }
    
    // remove zero-length entries and MEMs that aren't SMEMs
    // the zero-length ones are associated with single-base MEMs that tend to
    // match the entire index (typically Ns)
    // minor TODO: fix the above algorithm so they aren't introduced at all
    mems.erase(std::remove_if(mems.begin(), mems.end(),
                              [&smems_begin,
                               &min_mem_length](const MaximalExactMatch& m) {
                                  return ( m.end-m.begin == 0
                                          || m.length() < min_mem_length
                                          || smems_begin[m.begin] != m.end
                                          || m.count_Ns() > 0
                                          );
                              }),
               mems.end());
    // return the matches in natural order
    std::reverse(mems.begin(), mems.end());
    
    // fill the counts before deciding what to do
    for (auto& mem : mems) {
        if (mem.length() >= min_mem_length) {
            mem.match_count = gcsa->count(mem.range);
            if (hit_max) {
                gcsa->locate(mem.range, hit_max, mem.nodes);
            } else {
                gcsa->locate(mem.range, mem.nodes);
            }
        }
    }
    
    // reseed the long smems with shorter mems
    if (reseed_length) {
        // find if there are any mems that should be reseeded
        // iterate through MEMs
        vector<MaximalExactMatch> reseeded;
        for (auto& mem : mems) {
            // reseed if we have a long singular match
            if ((mem.length() >= reseed_length
                 && mem.match_count == 1)
                // or if we only have one mem for the entire read (even if it may have many matches)
                || mems.size() == 1) {
                // reseed at midway between here and the min mem length and at the min mem length
                int reseed_to = mem.length() / 2;
                int reseeds = 0;
                while (reseeds == 0 && reseed_to >= min_mem_length) {
#ifdef debug_mapper
                    cerr << "reseeding " << mem.sequence() << " with " << reseed_to << endl;
#endif
                    vector<MaximalExactMatch> remems = find_mems_simple(mem.begin,
                                                                        mem.end,
                                                                        reseed_to,
                                                                        min_mem_length,
                                                                        0);
                    reseed_to /= 2;
                    for (auto& rmem : remems) {
                        // keep if we have more than the match count of the parent
                        if (rmem.length() >= min_mem_length
                            && rmem.match_count > mem.match_count) {
                            ++reseeds;
                            reseeded.push_back(rmem);
                        }
                    }
                }
                // at least keep the original mem if needed
                if (reseeds == 0) {
                    reseeded.push_back(mem);
                }
            } else {
                reseeded.push_back(mem);
            }
        }
        mems = reseeded;
        // re-sort the MEMs by their start position
        std::sort(mems.begin(), mems.end(), [](const MaximalExactMatch& m1, const MaximalExactMatch& m2) { return m1.begin < m2.begin; });
    }
    return mems;
}

vector<MaximalExactMatch> BaseMapper::find_fanout_mems(string::const_iterator seq_begin,
                                                       string::const_iterator seq_end,
                                                       string::const_iterator qual_begin,
                                                       int max_fans_out,
                                                       char max_fanout_base_quality,
                                                       vector<deque<pair<string::const_iterator, char>>>* mem_fanout_breaks) {

#ifdef debug_mapper
    cerr << "find_fanout_mems: sequence " << string(seq_begin, seq_end) << ", max fan-out quality " << (int) max_fanout_base_quality << endl;
#endif
    
    vector<pair<uint8_t, string::const_iterator>> qual_pos;
    qual_pos.reserve(seq_end - seq_begin);
    for (auto seq_it = seq_begin, qual_it = qual_begin; seq_it != seq_end; ++seq_it, ++qual_it) {
        qual_pos.emplace_back(*(qual_it), seq_it);
    }
    make_heap(qual_pos.begin(), qual_pos.end(), greater<decltype(qual_pos)::value_type>());
    
    // find the k lowest base qualities and marking them as the fan out positions (breaking ties
    // by moving to start of read)
    vector<bool> do_fanout(seq_end - seq_begin, false);
    while (!qual_pos.empty() && qual_pos.front().first < max_fanout_base_quality
           && do_fanout.size() - qual_pos.size() < max_fans_out) {
        
        do_fanout[qual_pos.front().second - seq_begin] = true;
        pop_heap(qual_pos.begin(), qual_pos.end(), greater<decltype(qual_pos)::value_type>());
        qual_pos.pop_back();
    }
    
#ifdef debug_mapper
    cerr << "chose fan-out positions: ";
    for (size_t i = 0; i < do_fanout.size(); ++i) {
        if (do_fanout[i]) {
            cerr << i << " ";
        }
    }
    cerr << endl;
#endif

    // a struct to hold the information needed to continue a search from a fan-out
    struct FanOutSearch {
        // the SA range going into the  search
        gcsa::range_type prev_range;
        // where the search originally started
        string::const_iterator match_end;
        // where the search will continue from when this is dequeued
        string::const_iterator cursor;
        // the char we should use at the first search
        char fanout_char;
        // the where the MEM will need to be broken up
        deque<pair<string::const_iterator, char>> fanout_breaks;
        // record whether the iteration before the
        bool prev_iter_jumped_lcp;
        
        FanOutSearch(gcsa::range_type prev_range,
                     string::const_iterator seq_begin,
                     string::const_iterator match_end,
                     string::const_iterator cursor,
                     char fanout_char,
                     bool prev_iter_jumped_lcp,
                     const deque<pair<string::const_iterator, char>>& prev_fanout_breaks)
            : prev_range(prev_range), match_end(match_end), cursor(cursor), fanout_char(fanout_char),
              prev_iter_jumped_lcp(prev_iter_jumped_lcp), fanout_breaks(prev_fanout_breaks) {
            if (cursor >= seq_begin && *cursor != fanout_char) {
                fanout_breaks.emplace_front(cursor, fanout_char);
            }
        }
        FanOutSearch(gcsa::range_type prev_range,
                     string::const_iterator seq_begin,
                     string::const_iterator match_end,
                     string::const_iterator cursor,
                     char fanout_char)
            : prev_range(prev_range), match_end(match_end), cursor(cursor), fanout_char(fanout_char), prev_iter_jumped_lcp(false) {
            if (*cursor != fanout_char) {
                fanout_breaks.emplace_front(cursor, fanout_char);
            }
        }
    };
    
    gcsa::range_type full_range = gcsa::range_type(0, gcsa->size() - 1);
    
    // initialize the stack of search problems
    vector<FanOutSearch> stack;
    if (seq_begin < seq_end) {
        // walk backwards from end until finding a non-N base or a base where we want to fan out
        auto start_cursor = seq_end - 1;
        while (start_cursor >= seq_begin && *start_cursor == 'N' && !do_fanout[start_cursor - seq_begin]) {
            seq_end = start_cursor;
            --start_cursor;
        }
        if (start_cursor >= seq_begin) {
            if (do_fanout[start_cursor - seq_begin]) {
                // fan out at the first base
                for (char ch : {'A', 'C', 'G', 'T'}) {
                    stack.emplace_back(full_range, seq_begin, seq_end, start_cursor, ch);
                }
            }
            else {
                // don't fan out at the first base
                stack.emplace_back(full_range, seq_begin, seq_end, start_cursor, *start_cursor);
            }
        }
    }
    
    vector<tuple<gcsa::range_type, string::const_iterator, string::const_iterator, deque<pair<string::const_iterator, char>>>> search_results;
    
    while (!stack.empty()) {
                
        auto cursor = stack.back().cursor;
        auto match_end = stack.back().match_end;
        auto range = stack.back().prev_range;
        auto fanout_char = stack.back().fanout_char;
        auto fanout_breaks = move(stack.back().fanout_breaks);
        auto prev_iter_jumped_lcp = stack.back().prev_iter_jumped_lcp;
        stack.pop_back();
        
#ifdef debug_mapper
        cerr << "starting a fanout search with cursor at " << (cursor - seq_begin) << " and fanout char of " << fanout_char << endl;
        cerr << "current suffix: " << string(cursor, match_end) << endl;
        cerr << "breaks: " << endl;
        for (auto r : fanout_breaks) {
            cerr << "\t" << (r.first - seq_begin) << " -> " << r.second << endl;
        }
#endif
        
        
        bool use_fanout_char = true;
        
        while (cursor >= seq_begin) {
            
            // we only use the fan out character on the first iteration of a fan out search
            char ch = use_fanout_char ? fanout_char : *cursor;
            
#ifdef debug_mapper
            cerr << "LF iter at cursor " << (cursor - seq_begin) << " char " << ch << ", fanout? " << use_fanout_char << endl;
#endif
                        
            if (!use_fanout_char && do_fanout[cursor - seq_begin]) {
#ifdef debug_mapper
                cerr << "aborting to search to queue up fan-out searches" << endl;
#endif
                // fan out into all possiblities
                for (char nt : {'A', 'C', 'G', 'T'}) {
                    stack.emplace_back(range, seq_begin, match_end, cursor, nt, prev_iter_jumped_lcp, fanout_breaks);
                }
                // stop the current search
                break;
            }
            
            // hold onto our previous range
            auto prev_range = range;
            
            // execute one step of LF mapping
            range = gcsa->LF(range, gcsa->alpha.char2comp[ch]);
            
            if (gcsa::Range::empty(range) || match_end - cursor > gcsa->order() || ch == 'N') {
                
#ifdef debug_mapper
                cerr << "terminating search at end of a MEM" << endl;
#endif
                
                // we've exhausted our BWT range, so the last match range was maximal
                // or we have exceeded the order of the graph (FPs if we go further)
                
                if (cursor + 1 == match_end || ch == 'N') {
#ifdef debug_mapper
                    cerr << "skipping past a full index mismatch or N" << endl;
#endif
                    // avoid getting caught in infinite loop when a single character mismatches
                    // entire index (b/c then advancing the LCP doesn't move the search forward
                    // at all, need to move the cursor instead)
                    
                    if (!fanout_breaks.empty() && match_end - cursor - 1 <= fanout_length_threshold) {
                        // this fan-out search didn't find a long enough match for us to think
                        // that these are non-random matches
                        break;
                    }
                    
                    search_results.emplace_back(prev_range, cursor + 1, match_end, fanout_breaks);
                    
                    match_end = cursor;
                    range = full_range;
                    fanout_breaks.clear();
                    --cursor;
                    
                    prev_iter_jumped_lcp = false;
                }
                else {
#ifdef debug_mapper
                    cerr << "reached end of match" << endl;
#endif
                    auto match_begin = cursor + 1;
                    if (!fanout_breaks.empty() && match_end - match_begin <= fanout_length_threshold) {
                        // this fan-out search didn't find a long enough match for us to think
                        // that these are non-random matches
#ifdef debug_mapper
                        cerr << "match length of " << (match_end - match_begin) << " is below threshold " << fanout_length_threshold << ", aborting search" << endl;
#endif
                        break;
                    }

                    
                    // record the last MEM, but check to make sure were not actually still searching
                    // for the end of the next MEM
                    if (!prev_iter_jumped_lcp) {
#ifdef debug_mapper
                        cerr << "emitting a search result" << endl;
#endif
                        search_results.emplace_back(prev_range, match_begin, match_end, fanout_breaks);
                    }
                    
                    // init this outside of the if condition so that we can avoid calling the
                    // expensive LCP::parent function twice in the code path that checks max LCP
                    gcsa::STNode parent = lcp->parent(prev_range);
                    
                    // set the MEM to be the longest prefix that is shared with another MEM
                    match_end = match_begin + parent.lcp();
                    // and set up the next MEM using the parent node range
                    range = parent.range();
                    // forget about any fan outs that happened after the end of the new match
                    while (!fanout_breaks.empty() && fanout_breaks.back().first >= match_end) {
                        fanout_breaks.pop_back();
                    }
                    prev_iter_jumped_lcp = true;
                    
#ifdef debug_mapper
                    cerr << "after jumping LCP, suffix is " << string(cursor + 1, match_end) << endl;
#endif
                    
                    if (use_fanout_char) {
                        // this match ended on the fanout character, so we may need to remove
                        // that fan-out break from the search result
                        if (!get<3>(search_results.back()).empty() &&
                            get<3>(search_results.back()).front().first < get<1>(search_results.back())) {
                            get<3>(search_results.back()).pop_front();
                        }
                        // and skip the part of the iteration where we mark ourselves as being
                        // past the fan-out character
                        continue;
                    }
                    
                }
            }
            else {
                prev_iter_jumped_lcp = false;
                // just step to the next position
                --cursor;
            }
            
            use_fanout_char = false;
        }
        
        if (cursor < seq_begin) {
#ifdef debug_mapper
            cerr << "terminating search at start of the read" << endl;
#endif
            search_results.emplace_back(range, seq_begin, match_end, fanout_breaks);
        }
    }
    
    // filter out redundant MEMs, which are common in this algorithm (often when the
    // larger MEM includes a correct mismatch from a fan-out). locate() is the most
    // expensive part of this algorithm, so it's worth expending some effort to reduce
    // calls to it
    
    // first order by read interval start and then by decreasing match size and range
    // size, this ensures that a search result can only be fully redundant with another
    // one that is earlier in the vector
    sort(search_results.begin(), search_results.end(),
         [](const tuple<gcsa::range_type, string::const_iterator, string::const_iterator, deque<pair<string::const_iterator, char>>>& a,
            const tuple<gcsa::range_type, string::const_iterator, string::const_iterator, deque<pair<string::const_iterator, char>>>& b) {
        return (get<1>(a) < get<1>(b)
                || (get<1>(a) == get<1>(b) && get<2>(a) > get<2>(b))
                || (get<1>(a) == get<1>(b) && get<2>(a) == get<2>(b)
                    && gcsa::Range::length(get<0>(a)) > gcsa::Range::length(get<0>(b))));
    });
    
    size_t res_removed_so_far = 0;
    for (size_t i = 0, j = 1; j < search_results.size(); ++j) {
        // advance the starting index past the search results that cannot contain
        // this one
        while (get<2>(search_results[i]) < get<1>(search_results[j])) {
            ++i;
        }
        bool redundant = false;
        for (size_t k = i; k < j && !redundant; ++k) {
            // does the other result cover the entire read interval and all of its
            // locations in the graph?
            redundant = (get<2>(search_results[k]) >= get<2>(search_results[j])
                         && get<0>(search_results[k]).first <= get<0>(search_results[j]).first
                         && get<0>(search_results[k]).second >= get<0>(search_results[j]).second);
        }
        
        if (redundant) {
            ++res_removed_so_far;
        }
        else if (res_removed_so_far) {
            search_results[j - res_removed_so_far] = move(search_results[j]);
        }
    }
    
    if (res_removed_so_far) {
#ifdef debug_mapper
        cerr << "removing " << res_removed_so_far << " redundant search results" << endl;
#endif
        search_results.resize(search_results.size() - res_removed_so_far);
    }
    
    vector<MaximalExactMatch> mems;
    mems.reserve(search_results.size());
    if (mem_fanout_breaks) {
        mem_fanout_breaks->reserve(search_results.size());
    }
    for (auto it = search_results.begin(), end = search_results.end(); it != end; ++it) {
        
        const auto& search_result = *it;
        
        if (get<2>(search_result) - get<1>(search_result) < min_mem_length) {
            continue;
        }
        
        // there are no mismatches in the search, so the whole thing can become 1 MEM
        mems.emplace_back(get<1>(search_result), get<2>(search_result), get<0>(search_result),
                          gcsa->count(get<0>(search_result)));
        if (!hard_hit_max || mems.back().match_count < hard_hit_max) {
            if (hit_max) {
                gcsa->locate(mems.back().range, hit_max, mems.back().nodes);
                
            }
            else {
                gcsa->locate(mems.back().range, mems.back().nodes);
            }
        }
        
#ifdef debug_mapper
        cerr << "created MEM " << mems.back().sequence() << ", filled " << mems.back().nodes.size() << " of " << mems.back().match_count << " hits" << endl;
        for (auto n : mems.back().nodes) {
            cerr << "\t" << make_pos_t(n) << endl;
        }
#endif
        
        // keep track of the initial number of hits we query in case the nodes vector is
        // modified later (e.g. by prefiltering)
        mems.back().queried_count =  mems.back().nodes.size();
        
        if (mem_fanout_breaks) {
            // we're communicating fan-out breaks to the calling environment rather than
            // breaking them into exact matches right now
            
#ifdef debug_mapper
            cerr << "recording fanout breaks:" << endl;
            for (auto b : get<3>(search_result)) {
                cerr << "\t" << (b.first - seq_begin) << ": " << *b.first << " -> " << b.second << endl;
            }
#endif
            mem_fanout_breaks->emplace_back(move(get<3>(search_result)));
        }
        else if (!get<3>(search_result).empty()) {
#ifdef debug_mapper
            cerr << "need to break apart fanout MEM" << endl;
#endif
            
            // find the paths that each hit took
            vector<vector<pos_t>> paths;
            paths.reserve(mems.back().nodes.size());
            for (gcsa::node_type pos : mems.back().nodes) {
                paths.emplace_back(walk_fanout_path(get<1>(search_result),
                                                    get<2>(search_result),
                                                    get<3>(search_result),
                                                    pos));
            }
            
            // records of (pos index, bases past pos offset)
            vector<pair<size_t, size_t>> path_positions(paths.size(), pair<size_t, size_t>(0, 0));
            
            for (auto fanout_break : get<3>(search_result)) {
#ifdef debug_mapper
                cerr << "breaking fanout mems at " << (fanout_break.first - seq_begin) << " -> " << fanout_break.second << endl;
#endif
                
                // split up the match into two segments
                mems.back().end = fanout_break.first;
                if (fanout_break.first + 1 == seq_end) {
                    break;
                }
                mems.emplace_back(fanout_break.first + 1, get<2>(search_result),
                                  mems.back().range, mems.back().match_count);
                mems.back().queried_count = mems[mems.size() - 2].queried_count;
                
                // walk forward the specified length along each of the match paths
                size_t dist_to_walk = mems.back().begin - mems[mems.size() - 2].begin;
#ifdef debug_mapper
                cerr << "need to walk " << dist_to_walk << " from previous match positions" << endl;
#endif
                for (size_t i = 0; i < paths.size(); ++i) {
                    
                    pair<size_t, size_t>& path_pos = path_positions[i];
                    
                    size_t left_to_walk = dist_to_walk;
                    
                    pos_t path_step = paths[i][path_pos.first];
#ifdef debug_mapper
                    cerr << "starting a walk " << path_pos.second << " past " << path_step << endl;
#endif
                    size_t node_len = xindex->get_length(xindex->get_handle(id(path_step)));
                    size_t remain_len = node_len - offset(path_step) - path_pos.second;
                    while (remain_len < left_to_walk) {
                        left_to_walk -= remain_len;
                        ++path_pos.first;
                        path_pos.second = 0;
                        remain_len = xindex->get_length(xindex->get_handle(id(paths[i][path_pos.first])));
                    }
                    
                    path_pos.second += left_to_walk;
                    path_step = paths[i][path_pos.first];
                    
#ifdef debug_mapper
                    cerr << "\tadvance " << paths[i].front() << " to " << pos_t(id(path_step), is_rev(path_step), offset(path_step) + path_pos.second) << endl;
#endif
                    
                    // add the result as a position for this MEM
                    mems.back().nodes.emplace_back(gcsa::Node::encode(id(path_step),
                                                                      offset(path_step) + path_pos.second,
                                                                      is_rev(path_step)));
                    
                }
            }
        }
    }
    
    if (mems.size() == 1 && mems.front().length() >= mem_reseed_length && !mems.front().nodes.empty() &&
        (mem_fanout_breaks ? mem_fanout_breaks->front().empty() :
         find(mems.front().begin, mems.front().end, 'N') == mems.front().end)) {
        // try looking for reseed MEMs even if base qualities were high, just in case
        
        int min_sub_mem_length = max<int>(ceil(fast_reseed_length_diff * mems.front().length()), min_mem_length);
        
        vector<pair<MaximalExactMatch, vector<size_t>>> sub_mems;
        find_sub_mems_fast(mems, 0, 1, 0, mems.front().end, mems.front().begin, min_sub_mem_length, sub_mems);
        
        if (!sub_mems.empty()) {
            // we found a reseed sub-MEM
            
            // consolidate the sub MEMs into the MEM vector and record the parent relationships
            mems.reserve(mems.size() + sub_mems.size());
            vector<pair<int, vector<size_t>>> containment_graph(1);
            containment_graph.reserve(mems.size() + sub_mems.size());
            for (auto& sub_mem_and_parents : sub_mems) {
                mems.emplace_back(move(sub_mem_and_parents.first));
                if (mem_fanout_breaks) {
                    mem_fanout_breaks->emplace_back();
                }
                containment_graph.emplace_back(0, move(sub_mem_and_parents.second));
            }
            
            // try to remove redundant sub-MEMs
            if (prefilter_redundant_hits) {
                prefilter_redundant_sub_mems(mems, containment_graph);
            }
        }
    }
    
    auto cmp = [](const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
        if (m1.begin < m2.begin) {
            return true;
        }
        else if (m1.begin == m2.begin) {
            if (m1.end < m2.end) {
                return true;
            }
            else if (m1.end == m2.end) {
                return m1.nodes < m2.nodes;
            }
        }
        return false;
    };
    
    // put in lexicographic order
    if (!is_sorted(mems.begin(), mems.end(), cmp)) {
        vector<size_t> order(mems.size(), 0);
        for (size_t i = 1; i < order.size(); ++i) {
            order[i] = i;
        }
        stable_sort(order.begin(), order.end(), [&](size_t i, size_t j) { return cmp(mems[i], mems[j]); });
        vector<size_t> index(order.size());
        for (size_t i = 0; i < order.size(); ++i) {
            index[order[i]] = i;
        }
        for (size_t i = 0; i < index.size(); ++i) {
            while (index[i] != i) {
                swap(mems[i], mems[index[i]]);
                if (mem_fanout_breaks) {
                    swap((*mem_fanout_breaks)[i], (*mem_fanout_breaks)[index[i]]);
                }
                swap(index[i], index[index[i]]);
            }
        }
    }
    // remove null matches and duplicates
    size_t num_removed_so_far = 0;
    for (size_t i = 0; i < mems.size(); ++i) {
        if (mems[i].begin == mems[i].end ||
            (i >= num_removed_so_far + 1 && mems[i] == mems[i - num_removed_so_far - 1])) {
            // MEM is either empty or non-unique
            ++num_removed_so_far;
        }
        else if (num_removed_so_far) {
            // move the non-removed MEM past the removed ones
            mems[i - num_removed_so_far] = move(mems[i]);
            if (mem_fanout_breaks) {
                (*mem_fanout_breaks)[i - num_removed_so_far] = move((*mem_fanout_breaks)[i]);
            }
        }
    }
    
    if (num_removed_so_far) {
        mems.resize(mems.size() - num_removed_so_far);
        if (mem_fanout_breaks) {
            mem_fanout_breaks->resize(mem_fanout_breaks->size() - num_removed_so_far);
        }
    }
    
    return mems;
}

vector<pos_t> BaseMapper::walk_fanout_path(string::const_iterator begin,
                                           string::const_iterator end,
                                           const deque<pair<string::const_iterator, char>>& fanout_breaks,
                                           gcsa::node_type pos) {
#ifdef debug_mapper
    cerr << "beginning walk of fan-out path for sequence " << string(begin, end) << endl;
    cerr << "breaks:" << endl;
    for (auto b : fanout_breaks) {
        cerr << (b.first - begin) << " -> " << b.second << endl;
    }
#endif
    vector<pos_t> path;
    
    pos_t start_pos = make_pos_t(pos);
    handle_t start_handle = xindex->get_handle(id(start_pos), is_rev(start_pos));
    
    // records of (read pos, break pos, (node-offset records), next record to check)
    vector<tuple<string::const_iterator,
           deque<pair<string::const_iterator, char>>::const_iterator,
           vector<pair<handle_t, size_t>>,
           size_t>> stack;
    stack.emplace_back(begin,
                       fanout_breaks.begin(),
                       vector<pair<handle_t, size_t>>(1, make_pair(start_handle, offset(start_pos))),
                       0);
    
    while (!stack.empty()) {
        
        auto& stack_record = stack.back();
        
#ifdef debug_mapper
        cerr << "dequeueing stack record at relative idx " << (get<0>(stack_record) - begin) << ", option idx " << get<3>(stack_record) << " of " << get<2>(stack_record).size() << endl;
#endif
        
        if (get<3>(stack_record) == get<2>(stack_record).size()) {
#ifdef debug_mapper
            cerr << "popping stack record" << endl;
#endif
            
            // we've already traversed all of this nodes edges without finding a match
            stack.pop_back();
            
            continue;
        }
    
        
        // get the next position we're searching from
        handle_t handle;
        size_t off;
        tie(handle, off) = get<2>(stack_record)[get<3>(stack_record)++];
        
#ifdef debug_mapper
        cerr << "offset " << off << " on node " << xindex->get_id(handle) << ": " << xindex->get_sequence(handle) << endl;
#endif
        
        auto read_it = get<0>(stack_record);
        auto next_break = get<1>(stack_record);
        
        size_t node_len = xindex->get_length(handle);
        
        // check for a match along this node
        while (off < node_len && read_it != end) {
            // get either the read char or the fanout char
            char rch;
            if (next_break != fanout_breaks.end() && read_it == next_break->first) {
                rch = next_break->second;
                ++next_break;
            }
            else {
                rch = *read_it;
            }
#ifdef debug_mapper
            cerr << "\tlooking for match of " << rch << " to " << xindex->get_base(handle, off) << " at offset " << off << " on node length " << node_len << endl;
#endif
            
            if (rch != xindex->get_base(handle, off)) {
#ifdef debug_mapper
                cerr << "\tdoes not match" << endl;
#endif
                // this isn't a matching path
                break;
            }
            ++read_it;
            ++off;
        }
        
        if (read_it == end) {
            // we've finished walking out the match, the stack encodes
            // the path we took to get here
#ifdef debug_mapper
            cerr << "\tfound full length match" << endl;
#endif
            path.reserve(stack.size());
            for (const auto& search_record : stack) {
                handle_t h;
                size_t o;
                tie(h, o) = get<2>(search_record)[get<3>(search_record) - 1];
                path.emplace_back(xindex->get_id(h), xindex->get_is_reverse(h), o);
            }
            break;
        }
        else if (off == node_len) {
            // we walked to the end of the node, queue up the next nodes
#ifdef debug_mapper
            cerr << "\treached end of node" << endl;
#endif
            stack.emplace_back(read_it, next_break, vector<pair<handle_t, size_t>>(), 0);
            xindex->follow_edges(handle, false, [&](const handle_t& next) {
                get<2>(stack.back()).emplace_back(next, 0);
            });
        }
    }
    
#ifdef debug_mapper
    cerr << "walked path:";
    for (auto p : path) {
        cerr << " " << p;
    }
    cerr << endl;
#endif
    
    return path;
}

// Use the GCSA2 index to find super-maximal exact matches (and optionally sub-MEMs).
vector<MaximalExactMatch> BaseMapper::find_mems_deep(string::const_iterator seq_begin,
                                                     string::const_iterator seq_end,
                                                     double& longest_lcp,
                                                     double& fraction_filtered,
                                                     int max_mem_length,
                                                     int min_mem_length,
                                                     int reseed_length,
                                                     bool use_lcp_reseed_heuristic,
                                                     bool use_diff_based_fast_reseed,
                                                     bool include_parent_in_sub_mem_count,
                                                     bool record_max_lcp,
                                                     int reseed_below) {
#ifdef debug_mapper
    cerr << "find_mems: sequence " << string(seq_begin, seq_end) << ", max mem length " << max_mem_length << ", min mem length " << min_mem_length << ", reseed length " << reseed_length << endl;
#endif

    if (!gcsa) {
        cerr << "error:[vg::Mapper] a GCSA2 index is required to query MEMs" << endl;
        exit(1);
    }
    
    if (min_mem_length > reseed_length && reseed_length) {
        cerr << "error:[vg::Mapper] minimum reseed length for MEMs cannot be less than minimum MEM length" << endl;
        exit(1);
    }
    vector<MaximalExactMatch> mems;
    
    // an empty sequence matches the entire bwt
    if (seq_begin == seq_end) {
        mems.push_back(MaximalExactMatch(seq_begin, seq_end, gcsa::range_type(0, gcsa->size() - 1)));
    }
    
    // find SMEMs using GCSA+LCP array
    // algorithm sketch:
    // set up a cursor pointing to the last position in the sequence
    // set up a structure to track our MEMs, and set it == "" and full range match
    // while our cursor is >= the beginning of the string
    //   try a step of backwards searching using LF mapping
    //   if our range goes to 0
    //       go back to the last non-empty range
    //       emit the MEM corresponding to this range
    //       start a new mem
    //           use the LCP array's parent function to cut off the end of the match
    //           (effectively, this steps up the suffix tree)
    //           and calculate the new end point using the LCP of the parent node
    // emit the final MEM, if we finished in a matching state
    
    // next position we will extend matches to
    string::const_iterator cursor = seq_end - 1;
    // the end of the current MEM
    string::const_iterator curr_end = seq_end;
    // range of the current iteration
    gcsa::range_type range = accelerate_mem_query(seq_begin, cursor);
    
    // did we move the cursor or the end of the match last iteration?
    bool prev_iter_jumped_lcp = false;

    int filtered_mems = 0;
    int total_mems = 0;
    int max_lcp = 0;
    vector<int> lcp_maxima;
    

    // loop maintains invariant that match.range contains the hits for seq[cursor+1:match.end]
    while (cursor >= seq_begin) {
        
        // break the MEM on N; which for DNA we assume is non-informative
        // this *will* match many places in assemblies, but it isn't helpful
        if (*cursor == 'N') {
            auto curr_begin = cursor + 1;
            if (curr_end - curr_begin >= min_mem_length) {

                mems.emplace_back(curr_begin, curr_end, range);
                mems.back().match_count = gcsa->count(range);
                mems.back().primary = true;
                
                lcp_maxima.push_back(max_lcp);
                
#ifdef debug_mapper
                vector<gcsa::node_type> locations;
                if (hit_max) {
                    gcsa->locate(mems.back().range, hit_max, locations);
                } else {
                    gcsa->locate(mems.back().range, locations);
                }
                cerr << "adding MEM " << mems.back().sequence() << " after hitting N at positions ";
                for (auto nt : locations) {
                    cerr << make_pos_t(nt) << " ";
                }
                cerr << endl;
#endif
            }
            
            curr_end = cursor;
            --cursor;
            range = accelerate_mem_query(seq_begin, cursor);
            
            prev_iter_jumped_lcp = false;

            max_lcp = 0;

            // skip looking for matches since they are non-informative
            continue;
        }
        
        // hold onto our previous range
        auto last_range = range;
        
        // execute one step of LF mapping
        range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
        
        if (gcsa::Range::empty(range)
            || (max_mem_length && curr_end - cursor > max_mem_length)
            || curr_end - cursor > gcsa->order()) {
            
            // we've exhausted our BWT range, so the last match range was maximal
            // or: we have exceeded the order of the graph (FPs if we go further)
            // or: we have run over our parameter-defined MEM limit
            
            if (cursor + 1 == curr_end) {
                // avoid getting caught in infinite loop when a single character mismatches
                // entire index (b/c then advancing the LCP doesn't move the search forward
                // at all, need to move the cursor instead)
                auto curr_begin = cursor + 1;
                
                if (curr_end - curr_begin >= min_mem_length) {
                    mems.emplace_back(curr_begin, curr_end, last_range);
                    mems.back().match_count = gcsa->count(mems.back().range);
                    mems.back().primary = true;
                    lcp_maxima.push_back(max_lcp);
                }
                
                curr_end = cursor;
                --cursor;
                range = accelerate_mem_query(seq_begin, cursor);
                
                // don't reseed in empty MEMs
                prev_iter_jumped_lcp = false;
                max_lcp = 0;
            }
            else {
                auto curr_begin = cursor + 1;
                
                // record the last MEM, but check to make sure were not actually still searching
                // for the end of the next MEM
                bool add_mem = (curr_end - curr_begin >= min_mem_length && !prev_iter_jumped_lcp);
                if (add_mem) {
                    mems.emplace_back(curr_begin, curr_end, last_range);
                    mems.back().match_count = gcsa->count(mems.back().range);
                    mems.back().primary = true;
                    lcp_maxima.push_back(max_lcp);
                    
#ifdef debug_mapper
                    vector<gcsa::node_type> locations;
                    if (hit_max) {
                        gcsa->locate(mems.back().range, hit_max, locations);
                    } else {
                        gcsa->locate(mems.back().range, locations);
                    }
                    cerr << "adding MEM " << mems.back().sequence() << " after hitting an empty extension at positions ";
                    for (auto nt : locations) {
                        cerr << make_pos_t(nt) << " ";
                    }
                    cerr << endl;
#endif
                }
                
                // init this outside of the if condition so that we can avoid calling the
                // expensive LCP::parent function twice in the code path that checks max LCP
                gcsa::STNode parent;
                
                if (add_mem && use_greedy_mem_restarts
                    && curr_end - curr_begin >= greedy_restart_min_length
                    && mems.back().match_count <= greedy_restart_max_count) {
                    // the current match was relatively long and unique, so we think
                    // that we can get away with moving past it rather than looking
                    // for MEMs that overlap it on its left side
                    bool do_greedy_restart = true;
                    if (greedy_restart_max_lcp) {
                        // we also want to check whether this MEM has a long LCP (which
                        // would indicate that there is another match overlapping it)
                        parent = lcp->parent(last_range);
                        do_greedy_restart = (parent.lcp() <= greedy_restart_max_lcp);
#ifdef debug_mapper
                        cerr << "greedy restart aborted because of LCP of length " << parent.lcp() << endl;
#endif
                    }
                    if (do_greedy_restart) {
#ifdef debug_mapper
                        cerr << "doing greedy restart for next search iteration" << endl;
#endif
                        
                        // set up the search at the left end of the current MEM or one
                        // base further (which will create one fewer noise MEM if the current
                        // match was ended by a base substitution, but will make an artificially
                        // short MEM if it was ended by a deletion)
                        curr_end = greedy_restart_assume_substitution ? cursor : cursor + 1;
                        cursor = curr_end - 1;
                        range = accelerate_mem_query(seq_begin, cursor);
                        // we don't worry that we might still be jumping through prefixes
                        // of the current MEM because we've moved completely past it
                        prev_iter_jumped_lcp = false;
                        continue;
                    }
                }
                else {
                    // get the parent suffix tree node corresponding to the parent of the last MEM's STNode
                    parent = lcp->parent(last_range);
                }
                
                // set the match to be the longest prefix that is shared with another MEM
                curr_end = cursor + 1 + parent.lcp();
                // and set up the next MEM using the parent node range
                range = parent.range();
                // record our max lcp
                if (record_max_lcp) max_lcp = (int)parent.lcp();
                prev_iter_jumped_lcp = true;
            }
        }
        else {
            prev_iter_jumped_lcp = false;
            if (record_max_lcp) max_lcp = max(max_lcp, (int)lcp->parent(range).lcp());
            // just step to the next position
            --cursor;
        }
    }
    // TODO: is this where the bug with the duplicated MEMs is occurring? (when the prefix of a read
    // contains multiple non SMEM hits so that the iteration will loop through the LCP routine multiple
    // times before escaping out of the loop?
    
    // if we have a MEM at the beginning of the read, record it
    if (curr_end - seq_begin >= min_mem_length) {
        if (record_max_lcp) max_lcp = (int)lcp->parent(range).lcp();
        mems.emplace_back(seq_begin, curr_end, range);
        mems.back().match_count = gcsa->count(mems.back().range);
        mems.back().primary = true;
        lcp_maxima.push_back(max_lcp);
#ifdef debug_mapper
        vector<gcsa::node_type> locations;
        if (hit_max) {
            gcsa->locate(mems.back().range, hit_max, locations);
        } else {
            gcsa->locate(mems.back().range, locations);
        }
        cerr << "adding MEM " << mems.back().sequence() << " after hitting beginning of read at positions ";
        for (auto nt : locations) {
            cerr << make_pos_t(nt) << " ";
        }
        cerr << endl;
#endif
    }

    if (record_max_lcp) longest_lcp = lcp_maxima.empty() ? 0 : *max_element(lcp_maxima.begin(), lcp_maxima.end());

    assert(!record_max_lcp || lcp_maxima.size() == mems.size());
    
    // the graph of containment relationships between MEMs by index, also keeps track of what the sub-MEM min
    // length the children should be
    vector<pair<int, vector<size_t>>> sub_mem_containment_graph(mems.size());

    if (reseed_length) {
        // record what the minimum length of sub-MEMs they contain should be according to the parameters
        for (size_t i = 0; i < mems.size(); i++) {
            if (use_diff_based_fast_reseed && adaptive_reseed_diff) {
                sub_mem_containment_graph[i].first = max<int>(get_adaptive_min_reseed_length(mems[i].length()), min_mem_length);
            }
            else if (use_diff_based_fast_reseed) {
                sub_mem_containment_graph[i].first = max<int>(ceil(fast_reseed_length_diff * mems[i].length()), min_mem_length);
            }
            else {
                sub_mem_containment_graph[i].first = min_mem_length;
            }
        }
        
        // record where this layer of disjoint MEMs begins and ends in the return vector
        int layer_begin = 0, layer_end = mems.size();
        
        // apply sub-MEM algorithm to SMEMs and possibly then to the sub-MEMs
        for (int depth = 0; depth < max_sub_mem_recursion_depth; depth++) {
            // for all of the MEMs at this depth from the SMEMs in the containment graph
            for (int i = layer_begin; i < layer_end; i++) {
                
                MaximalExactMatch& mem = mems[i];
                // invalid mem
                // TODO: why are we ever producing these??
                if (mem.begin < seq_begin || mem.end > seq_end) continue;
                
                int min_sub_mem_length = sub_mem_containment_graph[i].first;
                
                // should we look for sub-MEMs from this MEM?
                // TODO: are the parentheses correct on the LCP heuristic
                if (mem.length() >= min_mem_length &&
                    mem.length() > min_sub_mem_length &&
                    (use_lcp_reseed_heuristic && record_max_lcp
                     && i < lcp_maxima.size() ? lcp_maxima[i] : mem.length()) >= reseed_length &&
                    (reseed_below == 0 || mem.match_count <= reseed_below)){
                    
                    // where should we start looking for sub-MEMs to ensure that we find independent hits
                    // from adjacent MEMs that might overlap this one?
                    string::const_iterator seed_boundary = i + 1 < layer_end ? mems[i + 1].end : mems[i].begin;
                    
                    // what's the leftmost index where a sub-MEM might be contained in the next MEM to the right?
                    string::const_iterator possible_containment_boundary = i - 1 >= layer_begin ? mems[i - 1].begin : mems[i].end;
                    
                    // reseed using the technique indicated by the mapper's parameters
                    vector<pair<MaximalExactMatch, vector<size_t>>> sub_mems;
                    if (fast_reseed) {
                        find_sub_mems_fast(mems, layer_begin, layer_end, i,
                                           possible_containment_boundary, seed_boundary,
                                           min_sub_mem_length, sub_mems);
                    }
                    else {
                        find_sub_mems(mems, layer_begin, layer_end, i, seed_boundary,
                                      min_sub_mem_length, sub_mems);
                    }
                    
                    for (pair<MaximalExactMatch, vector<size_t>>& sub_mem_and_parents : sub_mems) {
                        // move the MEM to the return vector and the parents to the containment graph
                        mems.emplace_back(move(sub_mem_and_parents.first));
                        sub_mem_containment_graph.emplace_back(0, move(sub_mem_and_parents.second));
                        
                        // set the minimum sub-MEM length to the maximum of its parents' minimum lengths
                        for (size_t j : sub_mem_containment_graph.back().second) {
                            sub_mem_containment_graph.back().first = max(sub_mem_containment_graph.back().first,
                                                                         sub_mem_containment_graph[j].first);
                        }
                    }
                }
            }
            
            // we're done with this containment layer, mark its boundaries for the next iteration
            layer_begin = layer_end;
            layer_end = mems.size();
        }
    }
    
    // determine the shortest length
    size_t min_mem_query_length = 0;
    if (filter_short_mems) {
        // we will filter out MEMs that are short relative to the longest MEMs
        size_t max_mem_length = 0;
        for (const auto& mem : mems) {
            max_mem_length = max<size_t>(max_mem_length, mem.length());
        }
        min_mem_query_length = max_mem_length * short_mem_filter_factor;
    }
    
    // query the locations of the hits
    // note: iterate in reverse so we remove the parent count from the children MEMs before decrementing
    // the parent count itself
    for (int64_t i = mems.size() - 1; i >= 0; i--) {
        
        MaximalExactMatch& mem = mems[i];
        
        if (!include_parent_in_sub_mem_count) {
            // remove the redundant hits with the parent from the total count
            for (size_t parent_idx : sub_mem_containment_graph[i].second) {
                mem.match_count -= mems[parent_idx].match_count;
            }
        }
        
        if (mem.match_count > 0 && mem.length() >= min_mem_query_length) {
            if (!hard_hit_max || mem.match_count < hard_hit_max) {
                if (hit_max) {
                    gcsa->locate(mem.range, hit_max, mem.nodes);
                    
                } else {
                    gcsa->locate(mem.range, mem.nodes);
                }
                // keep track of the initial number of hits we query in case the nodes vector is
                // modified later (e.g. by prefiltering)
                mem.queried_count = mem.nodes.size();
            }
        }
        filtered_mems += mem.match_count - mem.nodes.size();
        total_mems += mem.nodes.size();
#ifdef debug_mapper
        cerr << "MEM " << mem.sequence() << " has " << mem.match_count << " hits: ";
        for (auto nt : mem.nodes) {
            cerr << make_pos_t(nt) << ", ";
        }
        cerr << endl;
#endif
    }
    
    // try to identify which hits are redundant sub-MEMs and filter them out
    // note: do this before sorting the MEMs because the graph's adjacency list is ordered the same
    // as the return vector
    if (prefilter_redundant_hits && reseed_length) {
        prefilter_redundant_sub_mems(mems, sub_mem_containment_graph);
    }

    // TODO: uhhh... if we're actually producing these it's kind of a problem
    // remove strange MEMs
    mems.erase(std::remove_if(mems.begin(), mems.end(),
                              [&seq_begin, &seq_end](const MaximalExactMatch& mem)
                              { return mem.begin < seq_begin || mem.end > seq_end; }),
               mems.end());

    // return the MEMs in lexicographic order by the read interval
    std::sort(mems.begin(), mems.end(), [](const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
        return m1.begin < m2.begin || (m1.begin == m2.begin && m1.end < m2.end);
    });
    std::for_each(mems.begin(), mems.end(), [&total_mems](const MaximalExactMatch& mem) { total_mems += mem.nodes.size(); });
    fraction_filtered = (double) filtered_mems / (double) total_mems;
    
    // remove non-unique MEMs
    // TODO: I think I already fixed this
    mems.erase(unique(mems.begin(), mems.end()), mems.end());
    // remove MEMs that are overlapping positionally (they may be redundant)
    
    // sometimes we can merge hits that are only split up because we hit the order of the GCSA2 indexå
    if (precollapse_order_length_hits && gcsa->order() < seq_end - seq_begin) {
        precollapse_order_length_runs(seq_begin, mems);
    }

    return mems;
}

void BaseMapper::find_sub_mems(const vector<MaximalExactMatch>& mems,
                               int parent_layer_begin,
                               int parent_layer_end,
                               int mem_idx,
                               string::const_iterator next_mem_end,
                               int min_sub_mem_length,
                               vector<pair<MaximalExactMatch, vector<size_t>>>& sub_mems_out) {
    
    // get the most recently added MEM
    const MaximalExactMatch& mem = mems[mem_idx];
    
#ifdef debug_mapper
    cerr << "find_sub_mems: sequence ";
    for (auto iter = mem.begin; iter != mem.end; iter++) {
        cerr << *iter;
    }
    cerr << ", min mem length " << min_sub_mem_length << endl;
#endif
    
    // how many times does the parent MEM occur in the index?
    size_t parent_count = gcsa->count(mem.range);
    
    // next position where we will look for a match
    string::const_iterator cursor = mem.end - 1;
    
    // the righthand end of the sub-MEM we are building
    string::const_iterator sub_mem_end = mem.end;
    
    // the range that matches search_start:sub_mem_end
    gcsa::range_type range = gcsa::range_type(0, gcsa->size() - 1);
    
    // did we move the cursor or the end of the match last iteration?
    bool prev_iter_jumped_lcp = false;
    
    // look for matches that are contained in this MEM and not contained in the next MEM
    while (cursor >= mem.begin && sub_mem_end > next_mem_end) {
        // Note: there should be no need to handle N's or whole-index mismatches in this
        // routine (unlike the SMEM routine) since they should never make it into a parent
        // SMEM in the first place
        
        // hold onto our previous range
        gcsa::range_type last_range = range;
        // execute one step of LF mapping
        range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
        
        if (gcsa->count(range) <= parent_count) {
            // there are no more hits outside of parent MEM hits, record the previous
            // interval as a sub MEM
            string::const_iterator sub_mem_begin = cursor + 1;
            
            if (sub_mem_end - sub_mem_begin >= min_sub_mem_length && !prev_iter_jumped_lcp) {
                sub_mems_out.emplace_back(MaximalExactMatch(sub_mem_begin, sub_mem_end, last_range),
                                          vector<size_t>(1, mem_idx));
#ifdef debug_mapper
                vector<gcsa::node_type> locations;
                if (hit_max) {
                    gcsa->locate(last_range, hit_max, locations);
                } else {
                    gcsa->locate(last_range, locations);
                }
                cerr << "adding sub-MEM ";
                for (auto iter = sub_mem_begin; iter != sub_mem_end; iter++) {
                    cerr << *iter;
                }
                cerr << " at positions ";
                for (auto nt : locations) {
                    cerr << make_pos_t(nt) << " ";
                }
                cerr << endl;
#endif
                // identify all previous MEMs that also contain this sub-MEM
                for (int64_t i = mem_idx - 1; i >= parent_layer_begin; --i) {
                    if (sub_mem_begin >= mems[i].begin) {
                        // contined in next MEM, add its index to sub MEM's list of parents
                        sub_mems_out.back().second.push_back(i);
                    }
                    else {
                        // not contained in the next MEM, cannot be contained in earlier ones
                        break;
                    }
                }
            }
#ifdef debug_mapper
            else {
                cerr << "minimally more frequent MEM is too short ";
                for (auto iter = sub_mem_begin; iter != sub_mem_end; iter++) {
                    cerr << *iter;
                }
                cerr << endl;
            }
#endif
            
            // get the parent suffix tree node corresponding to the parent of the last MEM's STNode
            gcsa::STNode parent = lcp->parent(last_range);
            // set the MEM to be the longest prefix that is shared with another MEM
            sub_mem_end = sub_mem_begin + parent.lcp();
            // and get the next range as parent node range
            range = parent.range();
            
            prev_iter_jumped_lcp = true;
        }
        else {
            cursor--;
            prev_iter_jumped_lcp = false;
        }
    }
    
    // add a final sub MEM if there is one and it is not contained in the next parent MEM
    if (sub_mem_end > next_mem_end && sub_mem_end - mem.begin >= min_mem_length && !prev_iter_jumped_lcp) {
        sub_mems_out.emplace_back(MaximalExactMatch(mem.begin, sub_mem_end, range),
                                  vector<size_t>(1, mem_idx));
#ifdef debug_mapper
        cerr << "adding sub-MEM ";
        for (auto iter = mem.begin; iter != sub_mem_end; iter++) {
            cerr << *iter;
        }
        cerr << endl;
#endif
        // note: this sub MEM is at the far left side of the parent MEM, so we don't need to
        // check whether earlier MEMs contain it as well
    }
#ifdef debug_mapper
    else {
        cerr << "minimally more frequent MEM is too short ";
        for (auto iter = mem.begin; iter != sub_mem_end; iter++) {
            cerr << *iter;
        }
        cerr << endl;
    }
#endif
    
    
    // annotate the MEMs
    for (pair<MaximalExactMatch, vector<size_t>>& sub_mem_and_parents : sub_mems_out) {
        // count in entire range, including parents
        sub_mem_and_parents.first.match_count = gcsa->count(sub_mem_and_parents.first.range);
        // mark this is a sub-MEM
        sub_mem_and_parents.first.primary = false;
    }
}

// TODO: rewrite roughly as follows
// //// idea is to check all the possible sub mems starting from the same position so we obviate the need to call count for every LF step
// initally take the range of the MEM, then do parent operation to find the next-longest match with a long range
// get ranges for all suffixes of the MEM to reseed
// for each suffix range we keep the starting position fixed then use parent queries to find candidates (avoids excessive count calls)
// --- hope is that the first call to parent jumps us close to the limit for minimum length ...
// goal is to get to roughly a linear number of LFs an parent()s per MEM we reseed

void BaseMapper::find_sub_mems_fast(const vector<MaximalExactMatch>& mems,
                                    int parent_layer_begin,
                                    int parent_layer_end,
                                    int mem_idx,
                                    string::const_iterator leftmost_guaranteed_disjoint_bound,
                                    string::const_iterator leftmost_seeding_bound,
                                    int min_sub_mem_length,
                                    vector<pair<MaximalExactMatch, vector<size_t>>>& sub_mems_out) {
    
#ifdef debug_mapper
    cerr << "find_sub_mems_fast: mem ";
    for (auto iter = mems[mem_idx].begin; iter != mems[mem_idx].end; iter++) {
        cerr << *iter;
    }
    cerr << ", min_sub_mem_length " << min_sub_mem_length << endl;
#endif
    
    // how many times does the parent MEM occur in the index?
    size_t parent_range_count = use_approx_sub_mem_count ? gcsa::Range::length(mems[mem_idx].range) : mems[mem_idx].match_count;
    
    // the end of the leftmost substring that is at least the minimum length and not contained
    // in the next SMEM
    string::const_iterator probe_string_end = mems[mem_idx].begin + min_sub_mem_length;
    if (probe_string_end <= leftmost_seeding_bound) {
        probe_string_end = leftmost_seeding_bound + 1;
    }
    // the leftmost possible index of the next sub-MEM
    string::const_iterator leftmost_extension_bound = mems[mem_idx].begin;
    
    while (probe_string_end <= mems[mem_idx].end) {
        
        // locate the probe substring of length equal to the minimum length for a sub-MEM
        // that we are going to test to see if it's inside any sub-MEM
        string::const_iterator probe_string_begin = probe_string_end - min_sub_mem_length;
        
#ifdef debug_mapper
        cerr << "probe string is mem[" << probe_string_begin - mems[mem_idx].begin << ":" << probe_string_end - mems[mem_idx].begin << "] " << string(probe_string_begin, probe_string_end) << endl;
#endif
        
        // set up LF searching
        string::const_iterator cursor = probe_string_end - 1;
        gcsa::range_type range = accelerate_mem_query(probe_string_begin, cursor);
        
        // check if the probe substring is more frequent than the SMEM its contained in
        bool probe_string_more_frequent = true;
        while (cursor >= probe_string_begin) {
            
            range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
            
            // do a count operation if we've reached the beginning of the probe string or at invervals of the thinning parameter
            // past the burn-in parameter
            int64_t relative_idx = probe_string_end - cursor;
            if (cursor == probe_string_begin ||
                (relative_idx >= sub_mem_thinning_burn_in && (relative_idx - sub_mem_thinning_burn_in) % sub_mem_count_thinning == 0)) {
                
                if ((use_approx_sub_mem_count ? gcsa::Range::length(range) : gcsa->count(range)) <= parent_range_count) {
#ifdef debug_mapper
                    cerr << "partial probe only occurs in parent mem[" << cursor - mems[mem_idx].begin << ":" << probe_string_end - mems[mem_idx].begin << "] " << string(cursor, probe_string_end) << endl;
#endif
                    probe_string_more_frequent = false;
                    break;
                }
            }
            
            cursor--;
        }
        
        if (probe_string_more_frequent) {
            // this string is contained in a sub-MEM of length >= the minimum length, now we need to
            // find its ends
            
            if (probe_string_begin > leftmost_extension_bound) {
                // the probe string is contained in a sub-MEM, but that doesn't mean it's been
                // extended as far as possible to the left
                
                // extend match until beginning of SMEM or until the end of the independent hit
                while (cursor >= leftmost_extension_bound) {
                    gcsa::range_type last_range = range;
                    range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
                    
                    if ((use_approx_sub_mem_count ? gcsa::Range::length(range) : gcsa->count(range)) <= parent_range_count) {
                        range = last_range;
                        break;
                    }
                    
                    cursor--;
                }
                
                // mark this position as the beginning of the probe substring
                probe_string_begin = cursor + 1;
            }
            
            // now find the right end of the sub-MEM using binary search
            
            // inclusive interval that contains the past-the-last index of the sub-MEM
            string::const_iterator left_search_bound = probe_string_end;
            string::const_iterator right_search_bound = mems[mem_idx].end;
            
            // the match range of the longest prefix of the sub-MEM we've found so far (if we initialize
            // it here, the binary search is guaranteed to LF along the full sub-MEM in some iteration)
            gcsa::range_type sub_mem_range = range;
            
            // we'll keep track of the count of this sub-MEM during this loop and also remember it
            // for later
            size_t sub_mem_count = 0;
            
            // iterate until inteveral contains only one index
            while (right_search_bound > left_search_bound) {
                
                string::const_iterator middle = left_search_bound + (right_search_bound - left_search_bound + 1) / 2;
                
#ifdef debug_mapper
                cerr << "checking extension mem[" << probe_string_begin - mems[mem_idx].begin << ":" << middle - mems[mem_idx].begin << "] ";
                for (auto iter = probe_string_begin; iter != middle; iter++) {
                    cerr << *iter;
                }
                cerr << endl;
#endif
                
                // set up LF searching
                cursor = middle - 1;
                range = accelerate_mem_query(probe_string_begin, cursor);
                
                size_t extension_mem_count = 0;
                
                // check if there is an independent occurrence of this substring outside of the SMEM
                bool contained_in_independent_match = true;
                while (cursor >= probe_string_begin) {
                    
                    range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
                    
                    // do count operations on the final index and on intervals of the thinning parameter once we pass the
                    // burn in parameter
                    int64_t relative_idx = middle - cursor - 1;
                    if (cursor == probe_string_begin ||
                        (relative_idx >= sub_mem_thinning_burn_in
                         && (relative_idx - sub_mem_thinning_burn_in) % sub_mem_count_thinning == 0)) {
                        
                        extension_mem_count = use_approx_sub_mem_count ? gcsa::Range::length(range) : gcsa->count(range);
                        if (extension_mem_count <= parent_range_count) {
                            // this probe is too long and it no longer is contained in the indendent hit
                            // that we detected
                            contained_in_independent_match = false;
                            break;
                        }
                    }
                    
                    cursor--;
                }
                
                if (contained_in_independent_match) {
                    
                    // the end of the sub-MEM must be here or to the right
                    left_search_bound = middle;
                    
                    // update the range of matches (this is the longest match we've verified so far)
                    sub_mem_range = range;
                    
                    // this current count is the count of the longest match we've verified so far
                    sub_mem_count = extension_mem_count;
                }
                else {
                    // the end of the sub-MEM must be to the left
                    right_search_bound = middle - 1;
                }
            }
            
            bool possibly_redundant = (right_search_bound == mems[mem_idx].end && probe_string_begin >= leftmost_guaranteed_disjoint_bound);
            bool not_redundant = true;
            if (possibly_redundant) {
                // this sub-MEM is abutting the right end of the parent MEM and contained in the next MEM.
                // sometimes this means the sub-MEM is actually a truncated version of a sub-MEM of the next
                // MEM which was artificially cut short by the end of the current parent MEM, so we check that here
                
#ifdef debug_mapper
                cerr << "sub-MEM may be artificially truncated, checking extended sub-MEM" << endl;
                cerr << "current: " << MaximalExactMatch(probe_string_begin, right_search_bound, sub_mem_range).sequence() << endl;
#endif
                
                // get the GCSA range of the current sub-MEM extended one base past the end of the current parent MEM
                cursor = right_search_bound;
                range = accelerate_mem_query(probe_string_begin, cursor);
                size_t extended_count = numeric_limits<size_t>::max();
                bool contained_in_independent_match = true;
                while (cursor >= probe_string_begin) {
                    range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
                    int64_t relative_idx = right_search_bound - cursor + 1;
                    if (cursor == probe_string_begin ||
                        (relative_idx >= sub_mem_thinning_burn_in && (relative_idx - sub_mem_thinning_burn_in) % sub_mem_count_thinning == 0)) {
                        extended_count = use_approx_sub_mem_count ? gcsa::Range::length(range) : gcsa->count(range);
                        if (extended_count < sub_mem_count) {
                            // the count of the full lengthened sub-MEM will be fewer than the one we just found
                            break;
                        }
                    }
                    
                    cursor--;
                }
                
#ifdef debug_mapper
                cerr << "lengthened sub-MEM has count " << extended_count << ", compared to current count " << sub_mem_count << endl;
#endif
                
                // since the current sub-MEM is contained in the lengthened one, the lengthened sub-MEM can only have the same number
                // of hits or fewer. if it has the same number then the current sub-MEM is artificially truncated
                not_redundant = (extended_count < sub_mem_count);
                
            }
            
            if (not_redundant) {
#ifdef debug_mapper
                cerr << "final sub-MEM is mem[" << probe_string_begin - mems[mem_idx].begin << ":" << right_search_bound - mems[mem_idx].begin << "] ";
                for (auto iter = probe_string_begin; iter != right_search_bound; iter++) {
                    cerr << *iter;
                }
                cerr << endl;
#endif
                
                // record the sub-MEM
                sub_mems_out.emplace_back(MaximalExactMatch(probe_string_begin, right_search_bound, sub_mem_range),
                                          vector<size_t>(1, mem_idx));
                // annotate the sub-MEM
                // note: the count will not be accurate if we were using the approximate algorithm,
                // but we will handle that at the end of the function
                sub_mems_out.back().first.match_count = sub_mem_count;
                sub_mems_out.back().first.primary = false;
                
                // identify all previous MEMs that also contain this sub-MEM
                for (int64_t i = mem_idx - 1; i >= parent_layer_begin; --i) {
                    if (probe_string_begin >= mems[i].begin) {
                        // contained in next MEM, add its index to sub MEM's list of parents
                        sub_mems_out.back().second.push_back(i);
                    }
                    else {
                        // not contained in the next MEM, cannot be contained in earlier ones
                        break;
                    }
                }
                
                if (possibly_redundant && mem_idx + 1 != mems.size()) {
                    // the sub-MEM we might be redundant with is likely to be the last sub-MEM we identified (although not always).
                    // even if the current sub-MEM is non-redundant, it still is a child of that sub-MEM, so we should mark it as
                    // such to improve the efficiency of prefiltering. however, it's probably not worth looking too deeply into the
                    // MEM vector to find the parent if it is not the last one (avoids possible quadratic behavior)
                    if (mems.back().begin <= sub_mems_out.back().first.begin && mems.back().end >= sub_mems_out.back().first.end) {
                        sub_mems_out.back().second.push_back(mems.size() - 1);
                    }
                }
            }
            
            // the closest possible independent probe string will now occur one position
            // to the right of this sub-MEM
            probe_string_end = right_search_bound + 1;
            
            // any sub-MEMs will need to have their start position at least one position to the right
            // of this one
            leftmost_extension_bound = probe_string_begin + 1;
            
        }
        else {
            // we've found a suffix of the probe substring that is only contained inside the
            // parent SMEM, so we can move far enough to the right that the next probe substring
            // will not contain it
            
            probe_string_end = cursor + min_sub_mem_length + 1;
        }
    }
    
    if (use_approx_sub_mem_count) {
        // we didn't compute the exact count when looking for the sub-MEMs so now
        // we have to
        for (pair<MaximalExactMatch, vector<size_t>>& sub_mem_and_parents : sub_mems_out) {
            // count in entire range, including parents
            sub_mem_and_parents.first.match_count = gcsa->count(sub_mem_and_parents.first.range);
        }
    }
    
    
    // fast algorithm produces sub-MEMs left-to-right, switch the order so that they remain right-to-left
    // within the layer (as this algorithm expects in recursive calls)
    reverse(sub_mems_out.begin(), sub_mems_out.end());
}

vector<MaximalExactMatch> BaseMapper::find_stripped_matches(string::const_iterator seq_begin,
                                                            string::const_iterator seq_end,
                                                            size_t strip_length, size_t max_match_length,
                                                            size_t target_count) {
    if (!gcsa) {
        throw runtime_error("error:[vg::Mapper] a GCSA2 index is required to query matches");
    }
    if (strip_length <= 0) {
        throw runtime_error("error:[vg::Mapper] strip match length must be positive, set to " + to_string(strip_length));
    }
    if (target_count <= 0) {
        throw runtime_error("error:[vg::Mapper] target match count must be positive, set to " + to_string(target_count));
    }
    
#ifdef debug_strip_match
    cerr << "starting stripped match algorithm" << endl;
    cerr << "\tstrip length:" << strip_length << endl;
    cerr << "\tmax match length:" << max_match_length << endl;
    cerr << "\ttarget count:" << target_count << endl;
    cerr << "\tsequence:" << string(seq_begin, seq_end) << endl;
    
#endif
    
    // init the return value
    vector<MaximalExactMatch> matches;
    
    if (seq_end != seq_begin) {
        // we are not in the empty string
        int64_t seq_len = seq_end - seq_begin;
        int64_t num_strips = (seq_len - 1) / strip_length + 1;
        for (int64_t strip_num = 0; strip_num < num_strips; ++strip_num) {
            
#ifdef debug_strip_match
            cerr << "strip number " << strip_num << " of " << num_strips << endl;
#endif
            
            // the end of other strip match we will find
            auto strip_end = seq_end - strip_num * strip_length;
            // empty string starts matching entire index
            // note: because the stopping conditions are different here we can't
            // accelerate this MEM init query without changing results
            auto range = gcsa::range_type(0, gcsa->size() - 1);
            // a pointer to the next char we will match to
            auto cursor = strip_end - 1;
            
            while (cursor >= seq_begin &&
                   (!max_match_length || strip_end - cursor <= max_match_length)) {
                
                if (*cursor == 'N') {
                    // N matches are uninformative, so we don't want to match them
                    break;
                }
                
                // match one more char
                auto next_range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
                
#ifdef debug_strip_match
                cerr << "\tgot next range which is length " << gcsa::Range::length(next_range) << " and " << (gcsa::Range::empty(next_range) ? "" : "not ") << "empty" << endl;
#endif
                
                if (gcsa::Range::empty(next_range)) {
                    // we've gone too far, there are no more hits
                    break;
                }
                
                // the match was successful, advance to the range and move the cursor
                range = next_range;
                --cursor;
                
                if (target_count && gcsa::Range::length(next_range) <= target_count) {
                    // the (approximate) count is below the specified limit, so
                    // we've found reasonably unique sequence, stop looking for more
                    break;
                }
            }
            
            if (cursor + 1 == strip_end) {
                // edge case where one char mismatches the entire index, don't bother
                // with this
                continue;
            }
            
            if (!matches.empty()) {
                if (matches.back().begin <= cursor + 1 &&
                    matches.back().end >= strip_end) {
                    // this match is entirely contained within the larger match
                    // of the previous strip, so it's not likely to give us any
                    // new information
                    // also, filtering this helps us maintain reverse lexicographic
                    // order
                    continue;
                }
            }
            
#ifdef debug_strip_match
            cerr << "adding match of sequence " << string(cursor + 1, strip_end) << " and " << gcsa->count(range) << " hits" << endl;
#endif
            
            matches.emplace_back(cursor + 1, strip_end, range, gcsa->count(range));
            matches.back().primary = true;
            
            if (cursor < seq_begin) {
                // all further hits will be contained in ones we've already seen
                break;
            }
        }
    }
    
    // matches are queried in reverse lexicographic order, flip them around
    reverse(matches.begin(), matches.end());
    
    for (MaximalExactMatch& match : matches) {
        // figure out how many occurrences there are
        match.match_count = gcsa->count(match.range);
        if (!hard_hit_max || match.match_count < hard_hit_max) {
            // the total number of hits is low enough that we think it's at least
            // potentially worth querying hits
            if (hit_max) {
                // we may want to subsample
                gcsa->locate(match.range, hit_max, match.nodes);
                
            } else {
                // we won't subsample down to a prespecified maximum
                gcsa->locate(match.range, match.nodes);
            }
        }
        match.queried_count = match.nodes.size();
    }
    
    
    return matches;
}

gcsa::range_type BaseMapper::accelerate_mem_query(string::const_iterator begin,
                                                  string::const_iterator& cursor) const {
    
    if (!accelerator
        || cursor - begin < accelerator->length() - 1
        || find(cursor - accelerator->length() + 1, cursor + 1, 'N') <= cursor) {
        // we don't have an accelerator, the string we're querying is too short,
        // or there's an N (which we can't handle). decline to accelerate
        return gcsa::range_type(0, gcsa->size() - 1);
    }
    
    auto range = accelerator->memoized_LF(cursor);
    if (gcsa::Range::empty(range)) {
        // the acceleration lead to an empty range, so we might have been able to
        // get a non-empty one by actually doing each LF
        return gcsa::range_type(0, gcsa->size() - 1);
    }
    else {
        // success, update the cursor and return
        cursor -= accelerator->length();
        return range;
    }
}

void BaseMapper::prefilter_redundant_sub_mems(vector<MaximalExactMatch>& mems,
                                              vector<pair<int, vector<size_t>>>& sub_mem_containment_graph) {

    vector<size_t> sub_mem_range_end(mems.size());
    
    // for each MEM
    // note: reverse iteration is post order by construction
    for (int64_t i = mems.size() - 1; i >= 0; i--) {
        
        if (sub_mem_containment_graph[i].second.empty()) {
            // this MEM has no parents, so it cannot have redundant hits
            continue;
        }
        
#ifdef debug_mapper
        cerr << "prefiltering hits in MEM " << mems[i].sequence() << endl;
#endif
        // for all of its parents
        // get the supposed position this will be for each hit (if it is on the same node)
        unordered_set<pos_t> extended_hit_positions;
        for (size_t j : sub_mem_containment_graph[i].second) {
#ifdef debug_mapper
            cerr << "checking for extended hits from parent MEM " << mems[j].sequence() << endl;
#endif
            
            // how much farther should the sub-MEM be than the parent?
            int64_t relative_offset = mems[i].begin - mems[j].begin;
            for (gcsa::node_type hit : mems[j].nodes) {
                pos_t hit_pos = make_pos_t(hit);
                extended_hit_positions.emplace(id(hit_pos), is_rev(hit_pos), offset(hit_pos) + relative_offset);
#ifdef debug_mapper
                cerr << "\textending hit " << hit_pos << " to " << make_pos_t(id(hit_pos), is_rev(hit_pos), offset(hit_pos) + relative_offset) << endl;
#endif
            }
        }
        
        // remove any nodes whose hit location matches the extended position
        vector<gcsa::node_type>& nodes = mems[i].nodes;
        size_t removed_so_far = 0;
        for (size_t k = 0; k < nodes.size(); k++) {
            if (extended_hit_positions.count(make_pos_t(nodes[k]))) {
                // definitely redundant, skip this one
                removed_so_far++;
#ifdef debug_mapper
                cerr << "found redundant hit " << make_pos_t(nodes[k]) << ", removing" << endl;
#endif
            }
            else if (removed_so_far > 0) {
                // move this one up in the list to its new index accounting for skips
                nodes[k - removed_so_far] = nodes[k];
            }
        }
        
        // take off the end of the list, which only contains skipped hits and hits that have moved up
        if (removed_so_far) {
            nodes.resize(nodes.size() - removed_so_far);
        }
    }
}
    
void BaseMapper::precollapse_order_length_runs(string::const_iterator seq_begin,
                                               vector<MaximalExactMatch>& mems) {
    
    // find order length MEMs
    vector<size_t> order_length_mems;
    for (size_t i = 0; i < mems.size(); i++) {
        if (mems[i].length() == gcsa->order()) {
            order_length_mems.push_back(i);
        }
    }
    
    if (order_length_mems.empty()) {
        return;
    }
    
    vector<unordered_map<size_t, size_t>> collapsable_pairs(order_length_mems.size());
    for (size_t i = 0; i + 1 < order_length_mems.size(); i++) {
        
        MaximalExactMatch& extending_mem = mems[order_length_mems[i]];
        MaximalExactMatch& extend_to_mem = mems[order_length_mems[i + 1]];
        
        int64_t relative_offset = extend_to_mem.begin - extending_mem.begin;
        
        // find where hits should be if they are collapsable and on the same node
        unordered_map<pos_t, size_t> extended_hits;
        for (size_t j = 0; j < extending_mem.nodes.size(); j++) {
            pos_t pos = make_pos_t(extending_mem.nodes[j]);
            extended_hits[make_pos_t(id(pos), is_rev(pos), offset(pos) + relative_offset)] = j;
        }
        
        // check if we find any of these hits
        for (size_t j = 0; j < extend_to_mem.nodes.size(); j++) {
            pos_t pos = make_pos_t(extend_to_mem.nodes[j]);
            auto iter = extended_hits.find(pos);
            if (iter != extended_hits.end()) {
                collapsable_pairs[i].emplace(iter->second, j);
            }
        }
    }
    
    // did we find any order-length MEMs that we can collapse into longer MEMs a priori?
    if (!collapsable_pairs.empty()) {
        
        unordered_set<pair<size_t, size_t>> traversed;
        unordered_map<pair<int64_t, int64_t>, size_t> collapsed_mems;
        
        // make collapsed MEMs and identify their hits
        
        size_t num_uncollapsed_mems = mems.size();
        for (size_t i = 0; i + 1 < order_length_mems.size(); i++) {
            
            for (size_t j = 0; j < mems[order_length_mems[i]].nodes.size(); j++) {
                
                if (collapsable_pairs[i].count(j) && !traversed.count(make_pair(i, j))) {
                    // here's a collapsable hit we haven't yet traversed (i.e. the start of a new collapsed hit)
                    
#ifdef debug_mapper
                    cerr << "premerging MEM hits:" << endl;
                    cerr << "\t" << mems[order_length_mems[i]].sequence() << " " << make_pos_t(mems[order_length_mems[i]].nodes[j]) << endl;
#endif
                    // follow the forward links until we've traversed all the collapsable pairs
                    size_t at = j;
                    size_t col = i;
                    traversed.emplace(col, at);
                    pair<size_t, size_t> collapsable_run(at, at);
                    while (collapsable_pairs[col].count(at)) {
                        at = collapsable_pairs[col][at];
                        col++;
                        collapsable_run.second = at;
                        traversed.emplace(col, at);
#ifdef debug_mapper
                        cerr << "\t" << mems[order_length_mems[col]].sequence() << " " << make_pos_t(mems[order_length_mems[col]].nodes[at]) << endl;
#endif
                    }
                    
                    // find the interval of read indexes that the collapsed MEM will cover
                    pair<int64_t, int64_t> range(mems[order_length_mems[i]].begin - seq_begin,
                                                 mems[order_length_mems[col]].end - seq_begin);
                    
                    // make a new MEM for this read interval if we don't already have one
                    if (!collapsed_mems.count(range)) {
                        collapsed_mems[range] = mems.size();
                        mems.emplace_back(mems[order_length_mems[i]].begin,
                                          mems[order_length_mems[col]].end,
                                          gcsa::range_type());
                    }
                    
                    // get the MEM corresponding to this read interval
                    MaximalExactMatch& collapsed_mem = mems[collapsed_mems[range]];
                    
                    // add the hit to the collapsed MEM's hits
                    collapsed_mem.nodes.push_back(mems[order_length_mems[i]].nodes[j]);
                }
            }
        }
        
        // remove the hits that we collapsed from their MEMs
        
        for (size_t i = 0; i < order_length_mems.size(); i++) {
            
            size_t removed_so_far = 0;
            auto& mem = mems[order_length_mems[i]];
            auto& nodes = mem.nodes;
            for (size_t j = 0; j < nodes.size(); j++) {
                
                if (traversed.count(make_pair(i, j))) {
                    removed_so_far++;
                }
                else if (removed_so_far) {
                    nodes[j - removed_so_far] = nodes[j];
                }
            }
            
            if (removed_so_far) {
                nodes.resize(nodes.size() - removed_so_far);
            }
            
            mem.match_count -= removed_so_far;
            mem.queried_count -= removed_so_far;
            
        }
        
        // label the newly created MEMs with their hit count
        
        for (size_t i = num_uncollapsed_mems; i < mems.size(); i++) {
            auto& collapsed_mem = mems[i];
            collapsed_mem.match_count = collapsed_mem.nodes.size();
            collapsed_mem.queried_count = collapsed_mem.match_count;
            collapsed_mem.primary = true;
        }
        
        // resort the MEMs in lexicographic order by the read interval
        std::sort(mems.begin(), mems.end(), [](const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
            return m1.begin < m2.begin || (m1.begin == m2.begin && m1.end < m2.end);
        });
    }
}

void BaseMapper::rescue_high_count_order_length_mems(vector<MaximalExactMatch>& mems,
                                                     size_t max_rescue_hit_count) {
    
    vector<pair<size_t, size_t>> unfilled_mem_ranges;
    
    // identify the ranges of MEMs that are unfilled
    for (size_t i = 0; i < mems.size(); i++) {
        if (mems[i].nodes.empty()) {
            unfilled_mem_ranges.emplace_back(i, 0);
            while (i < mems.size() && mems[i].nodes.empty()) {
                i++;
            }
            unfilled_mem_ranges.back().second = i;
        }
    }
    
    for (pair<size_t, size_t>& mem_range : unfilled_mem_ranges) {
        // check that the range has a MEM that is maxed out at the order of the GCSA2
        // and also find the MEM in the tract with the minimum hit count
        bool has_order_length_mem = false;
        size_t min_hit_count = numeric_limits<size_t>::max();
        size_t min_hit_mem = numeric_limits<size_t>::max();
        for (size_t i = mem_range.first; i < mem_range.second; i++) {
            has_order_length_mem = has_order_length_mem || mems[i].length() == gcsa->order();
            if (mems[i].match_count < min_hit_count) {
                min_hit_count = mems[i].match_count;
                min_hit_mem = i;
            }
        }
        
        if (min_hit_count <= max_rescue_hit_count && has_order_length_mem) {
            // treat the minimum count MEM as a representative for this tract and fill it
            // to rescue the repeat mappings
            
#ifdef debug_mapper
            cerr << "found unfilled order length tract from MEM indexes " << mem_range.first << ":" << mem_range.second << ", filling with representative " << mems[min_hit_mem] << " with " << min_hit_count << " hits" << endl;
#endif

            if (hit_max) {
                gcsa->locate(mems[min_hit_mem].range, hit_max, mems[min_hit_mem].nodes);
            } else {
                gcsa->locate(mems[min_hit_mem].range, mems[min_hit_mem].nodes);
            }
        }
    }
}
    
size_t BaseMapper::get_adaptive_min_reseed_length(size_t parent_mem_length) {
    // extend memo until it contains this parent MEM length
    while (adaptive_reseed_length_memo.size() <= parent_mem_length) {
        size_t length = adaptive_reseed_length_memo.size();
        // a factor that decreases with longer MEMs but saturates
        double factor = (1.0 - fast_reseed_length_diff) * exp(-adaptive_diff_exponent * (length - mem_reseed_length + 1))
                         + fast_reseed_length_diff;
        adaptive_reseed_length_memo.push_back(length * factor);
    }
    return adaptive_reseed_length_memo[parent_mem_length];
}

bool PairedEndMapper::has_fixed_fragment_length_distr() {
    return fragment_length_distr.is_finalized();
}

void PairedEndMapper::force_fragment_length_distr(double mean, double stddev) {
    fragment_length_distr.force_parameters(mean, stddev);
}
    
void BaseMapper::apply_haplotype_consistency_scores(const vector<Alignment*>& alns) {
    if (haplo_score_provider == nullptr) {
        // There's no haplotype data available, so we can't add consistency scores.
        return;
    }
    
    if (xindex == nullptr) {
        // There's no database of haplotype names/counts available.
        // So we don't know how many haplotypes we should be looking for.
        return;
    }
    
    if (haplotype_consistency_exponent == 0) {
        // It won't matter either way
        return;
    }
   
    // Work out the population size. Try the score provider and then fall back to the xg.
    auto haplotype_count = haplo_score_provider->get_haplotype_count();
    if (haplotype_count == -1) {
        // The score provider doesn't have a haplotype count.
        haplotype_count = 0;
    }
   
    if (haplotype_count == 0 || haplotype_count == -1) {
        // We really should have a haplotype count
        throw runtime_error("Cannot score any haplotypes with a 0 or -1 haplotype count; are haplotypes available?");
    }
    
    // We don't look at strip_bonuses here, because we need these bonuses added
    // always in order to choose between alignments.
    
    // Build Yohei's recombination probability calculator. Feed it the haplotype
    // count from the PathPositionHandleGraph index that was generated alongside the GBWT.
    haplo::haploMath::RRMemo haplo_memo(recombination_penalty, haplotype_count);
    
    // This holds all the computed haplotype logprobs
    vector<double> haplotype_logprobs;
    haplotype_logprobs.reserve(alns.size());
    
    for (auto* aln : alns) {
        // On a first pass through, compute all the scores and make sure they all can be computed
        
        if (aln->path().mapping_size() == 0) {
            // Alignments with no actual mappings don't need scoring. But we
            // don't want to treat them as scoring failures, because we expect
            // some due to e.g. read pair mapping locations where one read maps
            // and the other needs rescue. We will skip them but continue on
            // with the rescoring, and also skip them when applying the scores.
            
            // Do a no-op adjustment
            haplotype_logprobs.push_back(0);
            
            continue;
        }
        
        // Score the path
        // This is a logprob (so, negative), and expresses the probability of the haplotype path being followed
        double haplotype_logprob;
        bool path_valid;
        std::tie(haplotype_logprob, path_valid) = haplo_score_provider->score(aln->path(), haplo_memo);
        
        if (std::isnan(haplotype_logprob) && path_valid) {
            // This shouldn't happen. Bail out on haplotype adjustment for this read and warn.
            cerr << "warning:[vg::Mapper]: NAN population score obtained for read with ostensibly successful query. Changing to failure." << endl;
            path_valid = false;
        }
        
        if (!path_valid) {
            // Our path does something the scorer doesn't like.
            // Bail out of applying haplotype scores.
            if (debug) {
                cerr << "Not applying haplotype consistency due to scoring failure" << endl;
            }
            return;
        }
        
        // Otherwise we haven't had a scoring failure yet, so keep going
        haplotype_logprobs.push_back(haplotype_logprob);
    }
    
    if (debug) {
        cerr << "Applying haplotype consistency to " << alns.size() << " alignment candidates" << endl;
    }
    
    for (size_t i = 0; i < alns.size(); i++) {
        // Get the aligner so we can convert from logprob to score points
        // TODO: This should always be the same aligner!
        auto* aligner = get_aligner(!alns[i]->quality().empty());
        assert(aligner->log_base != 0);
        
        if (alns[i]->path().mapping_size() != 0) {
            // We actually did rescore this one
            
            // This is a score "penalty" because it is usually negative. But positive = more score.
            // Convert to points, raise to haplotype consistency exponent power.
            double score_penalty = haplotype_consistency_exponent * (haplotype_logprobs[i] / aligner->log_base);
            
            // Apply "penalty"
            int64_t old_score = alns[i]->score();
            alns[i]->set_score(max((int64_t) 0, old_score + (int64_t) round(score_penalty)));
            // Note that we successfully corrected the score
            set_annotation(alns[i], "haplotype_score_used", true);
            // And save the score penalty/bonus
            set_annotation(alns[i], "haplotype_score", score_penalty);

            if (debug) {
                cerr << "Alignment statring at " << alns[i]->path().mapping(0).position().node_id()
                    << " got logprob " << haplotype_logprobs[i] << " vs " << haplotype_count
                    << " haplotypes, moving score by " << score_penalty
                    << " from " << old_score << " to " << alns[i]->score() << endl;
            }
        }
    }
}
    
double BaseMapper::estimate_gc_content(const gcsa::GCSA* gcsa) {
    
    uint64_t at = 0, gc = 0;
    
    if (gcsa) {
        at = gcsa::Range::length(gcsa->find(string("A"))) + gcsa::Range::length(gcsa->find(string("T")));
        gc = gcsa::Range::length(gcsa->find(string("G"))) + gcsa::Range::length(gcsa->find(string("C")));
    }
    
    if (at == 0 || gc == 0) {
        return default_gc_content;
    }
    
    return ((double) gc) / (at + gc);
}

int BaseMapper::random_match_length(double chance_random) {
    if (xindex) {
        // sum up the sequence length
        size_t length = 0;
        xindex->for_each_handle([&](const handle_t& handle) {
                length += xindex->get_length(handle);
            });
        return ceil(- (log(1.0 - pow(pow(1.0-chance_random, -1), (-1.0/length))) / log(4.0)));
    } else {
        return 0;
    }
}

void BaseMapper::set_alignment_scores(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend,
    int8_t full_length_bonus, double haplotype_consistency_exponent) {
    
    AlignerClient::set_alignment_scores(match, mismatch, gap_open, gap_extend, full_length_bonus);
    
    // Save the consistency exponent
    this->haplotype_consistency_exponent = haplotype_consistency_exponent;
}

void BaseMapper::set_alignment_scores(istream& matrix_stream, int8_t gap_open, int8_t gap_extend,
                                      int8_t full_length_bonus, double haplotype_consistency_exponent) {
    
    AlignerClient::set_alignment_scores(matrix_stream, gap_open, gap_extend, full_length_bonus);
    
    // Save the consistency exponent
    this->haplotype_consistency_exponent = haplotype_consistency_exponent;
}
    
void PairedEndMapper::set_fragment_length_distr_params(size_t maximum_sample_size, size_t reestimation_frequency,
                                                  double robust_estimation_fraction) {
    
    if (fragment_length_distr.is_finalized()) {
        cerr << "warning:[vg::PairedEndMapper] overwriting a fragment length distribution that has already been estimated" << endl;
    }
    
    fragment_length_distr = FragmentLengthDistribution(maximum_sample_size, reestimation_frequency,
                                                       robust_estimation_fraction);
}
    
Mapper::Mapper(PathPositionHandleGraph* xidex,
               gcsa::GCSA* g,
               gcsa::LCPArray* a,
               haplo::ScoreProvider* haplo_score_provider) :
    BaseMapper(xidex, g, a, haplo_score_provider)
    , thread_extension(10)
    , max_multimaps(1)
    , min_multimaps(4)
    , max_attempts(0)
    , min_cluster_length(0)
    , softclip_threshold(0)
    , max_softclip_iterations(10)
    , min_identity(0)
    , max_target_factor(128)
    , extra_multimaps(512)
    , band_multimaps(4)
    , always_rescue(false)
    , max_cluster_mapping_quality(1024)
    , use_cluster_mq(false)
    , simultaneous_pair_alignment(true)
    , drop_chain(0.2)
    , mq_overlap(0.2)
    , mate_rescues(0)
    , maybe_mq_threshold(0)
    , min_banded_mq(0)
    , max_band_jump(0)
    , patch_alignments(false)
    , identity_weight(2)
    , pair_rescue_hang_threshold(0.7)
    , pair_rescue_retry_threshold(0.5)
    , include_full_length_bonuses(true)
{
    // bench_init(bench[0]); bench_init(bench[1]); bench_init(bench[2]); bench_init(bench[3]);
    // counter[0] = 0; counter[1] = 0; counter[2] = 0; counter[3] = 0; counter[4] = 0; counter[5] = 0;
}

Mapper::Mapper(void) : BaseMapper() {
    // Nothing to do. Default constructed and can't really do anything.
    // bench_init(bench[0]); bench_init(bench[1]); bench_init(bench[2]); bench_init(bench[3]);
    // counter[0] = 0; counter[1] = 0; counter[2] = 0; counter[3] = 0; counter[4] = 0; counter[5] = 0;
}

Mapper::~Mapper(void) {
    //  Nothing to do, all managed memory handled in parent class
    /*
    fprintf(stderr, "mapper, time(%lu, %lu, %lu, %lu), cnt(%lu, %lu, %lu, %lu), counter(%lu, %lu, %lu, %lu)\n",
        bench_get(bench[0]) / 1000, bench_get(bench[1]) / 1000, bench_get(bench[2]) / 1000, bench_get(bench[3]) / 1000,
        bench_get_count(bench[0]), bench_get_count(bench[1]), bench_get_count(bench[2]), bench_get_count(bench[3]),
        counter[0], counter[1], counter[2], counter[3]);
    */
}


// todo add options for aligned global and pinned
Alignment Mapper::align_to_graph(const Alignment& aln,
                                 HandleGraph& graph,
                                 bool do_flip,
                                 bool traceback,
                                 bool pinned_alignment,
                                 bool pin_left,
                                 bool banded_global,
                                 bool keep_bonuses) {
    // do not use X-drop alignment when MEMs is not available
    vector<MaximalExactMatch> mems;
    return align_to_graph(aln, graph, mems, do_flip, traceback, pinned_alignment, pin_left, banded_global, false, keep_bonuses);
}

Alignment Mapper::align_to_graph(const Alignment& aln,
                                 HandleGraph& graph,
                                 const vector<MaximalExactMatch>& mems,
                                 bool do_flip,
                                 bool traceback,
                                 bool pinned_alignment,
                                 bool pin_left,
                                 bool banded_global,
                                 int xdrop_alignment,
                                 bool keep_bonuses) {

    // the longest path we could possibly align to (full gap and a full sequence)
    size_t target_length = aln.sequence().size() + get_aligner()->longest_detectable_gap(aln);

    // copy our alignment, which we'll then modify
    Alignment aligned = aln;

    // convert from bidirected to directed
    unordered_map<id_t, pair<id_t, bool> > node_trans;
    bdsg::HashGraph align_graph;
        
    // check if we can get away with using only one strand of the graph
    bool use_single_stranded = handlealgs::is_single_stranded(&graph);
    bool mem_strand = false;
    if (use_single_stranded) {
        if (mems.size()) {
            mem_strand = gcsa::Node::rc(mems[0].nodes.front());
            for (size_t i = 1; i < mems.size(); i++) {
                if (gcsa::Node::rc(mems[0].nodes.front()) != mem_strand) {
                    use_single_stranded = false;
                    break;
                }
            }
        } else if (do_flip) {
            mem_strand = true;
        }
    }
    bool flipped_alignment = false;
    if (use_single_stranded) {
        if (mem_strand && !xdrop_alignment) {
            aligned.set_sequence(reverse_complement(aligned.sequence()));
            if (!aligned.quality().empty()) {
                reverse(aligned.mutable_quality()->begin(),
                        aligned.mutable_quality()->end());
            }
            flipped_alignment = true;
        }
        if (mem_strand && xdrop_alignment) {
            // TODO -- investigate if reversing the mems is cheaper
            // xdrop requires that we reverse complement the mems or the graph
            handlealgs::reverse_complement_graph(&graph, &align_graph);
            node_trans.reserve(align_graph.get_node_count());
            align_graph.for_each_handle([&](const handle_t& handle) {
                node_trans[align_graph.get_id(handle)] = make_pair(align_graph.get_id(handle), true);
            });
        } else {
            // if we are using only the forward strand of the current graph, a make trivial node translation so
            // the later code's expectations are met
            // TODO: can we do this without the copy constructor?
            //align_graph = graph;
            node_trans.reserve(graph.get_node_count());
            graph.for_each_handle([&](const handle_t& handle) {
                    handle_t x = align_graph.create_handle(graph.get_sequence(handle), graph.get_id(handle));
                    assert(align_graph.get_id(x) == graph.get_id(handle));
                    node_trans[graph.get_id(handle)] = make_pair(graph.get_id(handle), false);
                    return true;
                });
            graph.for_each_edge([&](const edge_t& edge) {
                    align_graph.create_edge(graph.get_handle(graph.get_id(edge.first), graph.get_is_reverse(edge.first)),
                                            graph.get_handle(graph.get_id(edge.second), graph.get_is_reverse(edge.second)));
                });
        }
    }
    else {
        auto node_trans_tmp = handlealgs::split_strands(&graph, &align_graph);
        node_trans.reserve(node_trans_tmp.size());
        for (const auto& trans : node_trans_tmp) {
            node_trans[align_graph.get_id(trans.first)] = make_pair(graph.get_id(trans.second),
                                                                    graph.get_is_reverse(trans.second));
        }
    }

    // if necessary, convert from cyclic to acylic
    if (!handlealgs::is_directed_acyclic(&align_graph)) {
        // make a dagified graph and translation
        bdsg::HashGraph dagified;
        unordered_map<id_t,id_t> dagify_trans = handlealgs::dagify(&align_graph, &dagified, target_length);
        // replace the original with the dagified ones
        align_graph = move(dagified);
        node_trans = overlay_node_translations(dagify_trans, node_trans);
    }

    if (banded_global) {
        // the banded global alignment no longer constructs an internal representation of the graph
        // for topological sorting, etc., instead counting on the HandleGraph to do that itself. accordingly
        // we need a more serious/performant implementation of a graph here
        size_t max_span = aln.sequence().size();
        size_t band_padding_override = 0;
        bool permissive_banding = (band_padding_override == 0);
        size_t band_padding = permissive_banding ? max(max_span, (size_t) 1) : band_padding_override;
        get_aligner(!aln.quality().empty())->align_global_banded(aligned, align_graph, band_padding, false);
    } else if (pinned_alignment) {
        get_aligner(!aln.quality().empty())->align_pinned(aligned, align_graph, pin_left);
    } else if (xdrop_alignment) {
        get_aligner(!aln.quality().empty())->align_xdrop(aligned, align_graph,
                                                         translate_mems(mems, node_trans),
                                                         xdrop_alignment != 1,
                                                         max_xdrop_gap_length);
    } else {
        get_aligner(!aln.quality().empty())->align(aligned, align_graph, traceback);
    }
    if (traceback && !keep_bonuses && aligned.score()) {
        remove_full_length_bonuses(aligned);
    }
    // un-reverse complement the alignment
    if (flipped_alignment) {
        aligned = reverse_complement_alignment(
            aligned,
            (function<int64_t(int64_t)>) ([&](int64_t id) {
                    return align_graph.get_length(align_graph.get_handle(id));
                }));
    }
    if (node_trans.size()) {
        translate_oriented_node_ids(*aligned.mutable_path(), node_trans);
    }
    return aligned;
}

Alignment Mapper::align(const string& seq, int kmer_size, int stride, int max_mem_length, int band_width, int band_overlap, bool xdrop_alignment) {
    Alignment aln;
    aln.set_sequence(seq);
    return align(aln, kmer_size, stride, max_mem_length, band_width, band_overlap, xdrop_alignment);
}

vector<pos_t> Mapper::likely_mate_positions(const Alignment& aln, bool is_first_mate) {
    // when we don't have paths we can't properly estimate the mate position
    if (xindex->get_path_count() == 0) {
        return { };
    }
    unordered_map<path_handle_t, vector<pair<size_t, bool> > > offsets;
    for (auto& mapping : aln.path().mapping()) {
        auto pos_offs = algorithms::nearest_offsets_in_paths(xindex, make_pos_t(mapping.position()), aln.sequence().size());
        for (auto& p : pos_offs) {
            if (offsets.find(p.first)  == offsets.end()) {
                offsets[p.first] = p.second;
            }
        }
        if (offsets.size()) break; // find a single node that has a path position
    }
    // get our fragment model on the stack
    bool same_orientation = frag_stats.cached_fragment_orientation_same;
    bool forward_direction = frag_stats.cached_fragment_direction;
    int64_t delta = frag_stats.cached_fragment_length_mean;
    vector<pos_t> likely;
    for (auto& seq : offsets) {
        // find the likely position
        // then direction
        const std::string& seq_name = xindex->get_path_name(seq.first);
        for (auto& p : seq.second) { 
            size_t path_pos = p.first;
            bool on_reverse_path = p.second;
            int64_t mate_pos = 0;
            if (forward_direction) {
                if (is_first_mate) {
                    if (!on_reverse_path) {
                        mate_pos = path_pos + delta;
                    } else {
                        mate_pos = path_pos - delta;
                    }
                } else {
                    if (!on_reverse_path) {
                        mate_pos = path_pos + delta;
                    } else {
                        mate_pos = path_pos - delta;
                    }
                }
            } else {
                if (is_first_mate) {
                    if (!on_reverse_path) {
                        mate_pos = path_pos - delta;
                    } else {
                        mate_pos = path_pos + delta;
                    }
                } else {
                    if (!on_reverse_path) {
                        mate_pos = path_pos - delta;
                    } else {
                        mate_pos = path_pos + delta;
                    }
                }
            }
            path_handle_t path_handle = xindex->get_path_handle(seq_name);
            mate_pos = min(max((int64_t)0, mate_pos), (int64_t)xindex->get_path_length(path_handle)-1);
            //pos_t target = xindex->graph_pos_at_path_position(seq_name, mate_pos);
            step_handle_t step = xindex->get_step_at_position(path_handle, mate_pos);
            size_t pos_of_step = xindex->get_position_of_step(step);
            handle_t handle = xindex->get_handle_of_step(step);
            pos_t target = make_pos_t(xindex->get_id(handle), xindex->get_is_reverse(handle), mate_pos-pos_of_step);
            // what orientation should we use
            if ((same_orientation && on_reverse_path)
                || (!same_orientation && !on_reverse_path)) {
                target = reverse(target, get_node_length(id(target)));
            }
            likely.push_back(target);
        }
    }
    return likely;
}

pair<bool, bool> Mapper::pair_rescue(Alignment& mate1, Alignment& mate2,
                                     bool& tried1, bool& tried2,
                                     int match_score, int full_length_bonus, bool traceback, bool xdrop_alignment) {
    auto pair_sig = signature(mate1, mate2);
    // bail out if we can't figure out how far to go
    bool rescued1 = false;
    bool rescued2 = false;
    if (!frag_stats.fragment_size) return make_pair(false, false);
    double hang_threshold = pair_rescue_hang_threshold;
    double retry_threshold = pair_rescue_retry_threshold;
    double min_threshold = 0.5;
    double perfect_score = mate1.sequence().size() * match_score + full_length_bonus * 2;
    double accept_pval = 1e-6;
    // based on our statistics about the alignments
    // get the subgraph overlapping the likely candidate position of the second alignment
    bool rescue_off_first = false;
    bool rescue_off_second = false;
    double mate1_id = (double) mate1.score() / perfect_score;
    double mate2_id = (double) mate2.score() / perfect_score;
    //cerr << "---------------------------" << pb2json(mate1) << endl << pb2json(mate2) << endl;
    //cerr << "---------------------------" << endl;
#ifdef debug_rescue
    if (debug) cerr << "pair rescue: mate1 " << signature(mate1) << " " << mate1_id << " mate2 " << signature(mate2) << " " << mate2_id << endl;
    if (debug) cerr << "mate1: " << pb2json(mate1) << endl;
    if (debug) cerr << "mate2: " << pb2json(mate2) << endl;
#endif
    vector<pos_t> mate_positions;
    if (mate1_id > mate2_id && mate1_id > hang_threshold && mate2_id == 0) {
        // retry off mate1
#ifdef debug_rescue
        if (debug) cerr << "Rescue read 2 off of read 1" << endl;
#endif
        rescue_off_first = true;
        // record id and direction to second mate
        mate_positions = likely_mate_positions(mate1, true);
    } else if (mate2_id > mate1_id && mate2_id > hang_threshold && mate1_id == 0) {
        // retry off mate2
#ifdef debug_rescue
        if (debug) cerr << "Rescue read 1 off of read 2" << endl;
#endif
        rescue_off_second = true;
        // record id and direction to second mate
        mate_positions = likely_mate_positions(mate2, false);
    } else {
        return make_pair(false, false);
    }
    if (mate_positions.empty()) return make_pair(false, false); // can't rescue because the selected mate is unaligned
    bdsg::HashGraph graph;
    set<bool> orientations;
#ifdef debug_rescue
    if (debug) cerr << "got " << mate_positions.size() << " mate positions" << endl;
#endif
    for (auto& mate_pos : mate_positions) {
#ifdef debug_rescue
        if (debug) cerr << "aiming for " << mate_pos << endl;
#endif
        orientations.insert(is_rev(mate_pos));
        int get_at_least = (!frag_stats.cached_fragment_length_mean ? frag_stats.fragment_max
                            : min(frag_stats.fragment_max/2,
                                  (int64_t)max((double)frag_stats.cached_fragment_length_stdev * 10.0,
                                               mate1.sequence().size() * 3.0)));
        //cerr << "Getting at least " << get_at_least << endl;
        // TODO it may be possible to make this sugraph smaller with little penalty in terms of accuracy
        algorithms::extract_context(*xindex, graph, xindex->get_handle(id(mate_pos), is_rev(mate_pos)), offset(mate_pos), get_at_least);
        // if we're reversed, align the reverse sequence and flip it back
        // align against it
    }
    int max_mate1_score = mate1.score();
    int max_mate2_score = mate2.score();
    for (auto& orientation : orientations) {
        if (rescue_off_first) {
            Alignment aln2 = align_maybe_flip(mate2, graph, orientation, traceback, false, xdrop_alignment);
            tried2 = true;
            //vg::io::write_to_file(aln2, "rescue-" + h + ".gam");
#ifdef debug_rescue
            if (debug) cerr << "aln2 score/ident vs " << aln2.score() << "/" << aln2.identity()
                            << " vs " << mate2.score() << "/" << mate2.identity() << endl;
#endif
            if (aln2.score() > max_mate2_score && (double)aln2.score()/perfect_score > min_threshold && pair_consistent(mate1, aln2, accept_pval)) {
                if (!traceback) { // now get the traceback
                    aln2 = align_maybe_flip(mate2, graph, orientation, true, false, xdrop_alignment);
                }
#ifdef debug_rescue
                if (debug) cerr << "rescued aln2 " << pb2json(aln2) << endl;
#endif
                mate2 = aln2;
                mate2.clear_refpos();
                max_mate2_score = mate2.score();
                rescued2 = true;
            }
        } else if (rescue_off_second) {
            Alignment aln1 = align_maybe_flip(mate1, graph, orientation, traceback, false, xdrop_alignment);
            tried1 = true;
            //vg::io::write_to_file(aln1, "rescue-" + h + ".gam");
#ifdef debug_rescue
            if (debug) cerr << "aln1 score/ident vs " << aln1.score() << "/" << aln1.identity()
                            << " vs " << mate1.score() << "/" << mate1.identity() << endl;
#endif
            if (aln1.score() > max_mate1_score && (double)aln1.score()/perfect_score > min_threshold && pair_consistent(aln1, mate2, accept_pval)) {
                if (!traceback) { // now get the traceback
                    aln1 = align_maybe_flip(mate1, graph, orientation, true, false, xdrop_alignment);
                }
#ifdef debug_rescue
                if (debug) cerr << "rescued aln1 " << pb2json(aln1) << endl;
#endif
                mate1 = aln1;
                mate1.clear_refpos();
                max_mate1_score = mate1.score();
                rescued1 = true;
            }
        }
    }
    // if the new alignment is better
    // set the old alignment to it
    return make_pair(rescued1, rescued2);
}

bool Mapper::alignments_consistent(const map<string, double>& pos1,
                                   const map<string, double>& pos2,
                                   int fragment_size_bound) {
    set<string> comm_refs;
    for (auto& p : pos1) {
        auto& name = p.first;
        if (pos2.find(name) != pos2.end()) {
            comm_refs.insert(name);
        }
    }
    // Alignments are consistent if their median node id positions are within the fragment_size
    
    // get median mapping positions
    for (auto& ref : comm_refs) {
        // this is unsafe looking, but we know they have the same keys for these values
        auto mean1 = pos1.find(ref)->second;
        auto mean2 = pos2.find(ref)->second;
        if (abs(mean1 - mean2) < fragment_size_bound) {
            return true;
        }
    }
    return false;
}

bool Mapper::pair_consistent(Alignment& aln1,
                             Alignment& aln2,
                             double pval) {
    if (!(aln1.score() && aln2.score())) return false;
    bool length_ok = false;
    if (xindex->get_path_count() == 0) {
        // use the approximate distance
        //cerr << "using approx distance" << endl;
        int len = approx_fragment_length(aln1, aln2);
        if ((frag_stats.fragment_size && len > 0 && ((pval > 0 && frag_stats.fragment_length_pval(len) > pval)
                                                    || len < frag_stats.fragment_size))
            || (!frag_stats.fragment_size && len > 0 && len < frag_stats.fragment_max)) {
            length_ok = true;
        }
        bool aln1_is_rev = aln1.path().mapping(0).position().is_reverse();
        bool aln2_is_rev = aln2.path().mapping(0).position().is_reverse();
        bool same_orientation = frag_stats.cached_fragment_orientation_same;
        // XXX todo
        //bool direction_ok = frag_stats.cached_fragment_direction && 
        bool orientation_ok = (same_orientation && aln1_is_rev == aln2_is_rev)
            || (!same_orientation && aln1_is_rev != aln2_is_rev);
        return length_ok && orientation_ok;
    } else {
        // use the distance induced by the graph paths
        // won't recompute if we already have refpos
        algorithms::annotate_with_initial_path_positions(*xindex, aln1);
        algorithms::annotate_with_initial_path_positions(*xindex, aln2);
        map<string, vector<pair<size_t, bool> > > offsets1 = alignment_refpos_to_path_offsets(aln1);
        map<string, vector<pair<size_t, bool> > > offsets2 = alignment_refpos_to_path_offsets(aln2);
        auto pos_consistent =
            [&](const pair<size_t, bool>& p1, const pair<size_t, bool>& p2) {
            int64_t pos1 = p1.first;
            int64_t pos2 = p2.first;
            bool fwd1 = p1.second;
            bool fwd2 = p2.second;
            int64_t len = pos2 - pos1;
            if (frag_stats.fragment_size) {
                bool orientation_ok = (frag_stats.cached_fragment_orientation_same && fwd1 == fwd2) || fwd1 != fwd2;
                bool direction_ok = frag_stats.cached_fragment_direction && (!fwd1 && len >= 0 || fwd1 && len <= 0)
                    || (fwd1 && len >= 0 || !fwd1 && len <= 0);
                bool length_ok = frag_stats.fragment_length_pval(abs(len)) > pval;//|| pval == 0 && abs(len) < frag_stats.fragment_size;
                return orientation_ok && direction_ok && length_ok;
            } else {
                return abs(len) < frag_stats.fragment_max;
            }
        };
        for (auto& path : offsets1) {
            // see if the other alignment has it
            auto f = offsets2.find(path.first);
            if (f != offsets2.end()) {
                auto& pos1s = path.second;
                auto& pos2s = f->second;
                // in the cartesian product of the mapping positions is there one pair that would match our model?
                // TODO linearize this as it could get bad if we have lots of paths!
                for (auto& pos1 : pos1s) {
                    for (auto& pos2 : pos2s) {
                        if (pos_consistent(pos1, pos2)) return true;
                    }
                }
            }
        }
    }
    return false;
}

pair<vector<Alignment>, vector<Alignment>> Mapper::align_paired_multi(
    const Alignment& first_mate,
    const Alignment& second_mate,
    bool& queued_resolve_later,
    int max_mem_length,
    bool only_top_scoring_pair,
    bool retrying,
    bool xdrop_alignment) {

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    Alignment read1;
    read1.set_name(first_mate.name());
    read1.set_sequence(first_mate.sequence());
    read1.set_quality(first_mate.quality());
    Alignment read2;
    read2.set_name(second_mate.name());
    read2.set_sequence(second_mate.sequence());
    read2.set_quality(second_mate.quality());

    auto aligner = get_aligner(!read1.quality().empty());
    int8_t match = aligner->match;
    int8_t gap_extension = aligner->gap_extension;
    int8_t gap_open = aligner->gap_open;
    int8_t full_length_bonus = aligner->full_length_bonus;

    int total_multimaps = max(max_multimaps, extra_multimaps/2);
    double cluster_mq = 0;

    if(debug) {
        cerr << "align_paired_multi_simul "
             << "with " << read1.name() << " and "
             << read2.name() << endl
            //<< "read 1 " << read1.sequence() << endl
            //<< "read 2 " << read2.sequence() << endl
            //<< "read 1 " << pb2json(read1) << endl
            //<< "read 2 " << pb2json(read2) << endl
             << "fragment model " << frag_stats.fragment_max << ", "
             << frag_stats.fragment_size << ", "
             << frag_stats.cached_fragment_length_mean << ", "
             << frag_stats.cached_fragment_length_stdev << ", "
             << frag_stats.cached_fragment_orientation_same << ", "
             << frag_stats.cached_fragment_direction << ", "
             << frag_stats.since_last_fragment_length_estimate << ", " << endl;
    }

    pair<vector<Alignment>, vector<Alignment>> results;
    double longest_lcp1, longest_lcp2, fraction_filtered1, fraction_filtered2;
    // find the MEMs for the alignments
    vector<MaximalExactMatch> mems1 = find_mems_deep(read1.sequence().begin(),
                                                     read1.sequence().end(),
                                                     longest_lcp1,
                                                     fraction_filtered1,
                                                     max_mem_length,
                                                     min_mem_length,
                                                     mem_reseed_length,
                                                     false, true, true, true); // Make sure to actually fill in the longest LCP.

    vector<MaximalExactMatch> mems2 = find_mems_deep(read2.sequence().begin(),
                                                     read2.sequence().end(),
                                                     longest_lcp2,
                                                     fraction_filtered2,
                                                     max_mem_length,
                                                     min_mem_length,
                                                     mem_reseed_length,
                                                     false, true, true, true); // Make sure to actually fill in the longest LCP.

    double mq_cap1, mq_cap2;
    mq_cap1 = mq_cap2 = max_mapping_quality;

    int total_mem_length1 = 0;
    for (auto& mem : mems1) total_mem_length1 += mem.length(); // * mem.nodes.size();
    int total_mem_length2 = 0;
    for (auto& mem : mems2) total_mem_length2 += mem.length(); // * mem.nodes.size();
    double mem_read_ratio1 = min(1.0, (double)total_mem_length1 / (double)read1.sequence().size());
    double mem_read_ratio2 = min(1.0, (double)total_mem_length2 / (double)read2.sequence().size());

    int basis_length = min((int)read1.sequence().size(), (int)gcsa->order());
    double max_possible_mq = max_possible_mapping_quality(basis_length);

    int mem_max_length1 = 0;
    for (auto& mem : mems1) if (mem.primary && mem.match_count) mem_max_length1 = max(mem_max_length1, (int)mem.length());
    double maybe_mq1 = 0;
    if (mems1.size()) {
        maybe_mq1 = estimate_max_possible_mapping_quality(basis_length,
                                                          basis_length/max(1.0, (double)mem_max_length1),
                                                          basis_length/longest_lcp1);
    }
    int mem_max_length2 = 0;
    for (auto& mem : mems2) if (mem.primary && mem.match_count) mem_max_length2 = max(mem_max_length2, (int)mem.length());
    double maybe_mq2 = 0;
    if (mems2.size()) {
        maybe_mq2 = estimate_max_possible_mapping_quality(basis_length,
                                                          basis_length/max(1.0, (double)mem_max_length2),
                                                          basis_length/longest_lcp2);
    }

    // use the estimated mapping quality to avoid hard work when the results are likely noninformative
    double maybe_pair_mq = maybe_mq1+maybe_mq2;

    // if estimated mq is high scale difficulty using the estimated mapping quality
    /*
    if (maybe_pair_mq > max_mapping_quality) {
        total_multimaps = max(max(min_multimaps, max_multimaps), min(total_multimaps, (int)round(maybe_pair_mq)));
    }
    */

    if (debug) cerr << "maybe_mq1 " << read1.name() << " " << maybe_mq1 << " " << total_multimaps << " " << mem_max_length1 << " " << longest_lcp1 << " " << total_multimaps << " " << mem_read_ratio1 << " " << fraction_filtered1 << " " << max_possible_mq << " " << total_multimaps << endl;
    if (debug) cerr << "maybe_mq2 " << read2.name() << " " << maybe_mq2 << " " << total_multimaps << " " << mem_max_length2 << " " << longest_lcp2 << " " << total_multimaps << " " << mem_read_ratio2 << " " << fraction_filtered2 << " " << max_possible_mq << " " << total_multimaps << endl;
    
    if (debug) cerr << "mems for read 1 " << mems_to_json(mems1) << endl;
    if (debug) cerr << "mems for read 2 " << mems_to_json(mems2) << endl;

    auto transition_weight = [&](const MaximalExactMatch& m1, const MaximalExactMatch& m2) {

#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) cerr << "Compute distance from " << m1 << " to " << m2 << endl;
        }
#endif

        // set up positions for distance query
        pos_t m1_pos = make_pos_t(m1.nodes.front());
        pos_t m2_pos = make_pos_t(m2.nodes.front());

        // are the two mems in a different fragment?
        // we handle the distance metric differently in these cases
        if (m1.fragment < m2.fragment) {
            int64_t max_length = frag_stats.fragment_max;
            pair<int64_t, int64_t> d = mem_min_oriented_distances(m1, m2);
            // if we have a cached fragment orientation, use it to pick the min distance with the correct path relative orientation
            int64_t approx_dist = (!frag_stats.fragment_size ? min(d.first, d.second)
                                   : (frag_stats.cached_fragment_orientation_same ? d.first : d.second));
            if (approx_dist >= max_length) {
                // too far apart or wrong path/pair relative orientation
                return -std::numeric_limits<double>::max();
            } else if (frag_stats.fragment_size) {
                return frag_stats.fragment_length_pval(approx_dist) * m2.length();
            } else {
                return 1.0/approx_dist * (m1.length() + m2.length());
            }
        } else if (m1.fragment > m2.fragment) {
            // don't allow going backwards in the threads
            return -std::numeric_limits<double>::max();
        } else {
            pos_t m1_pos = make_pos_t(m1.nodes.front());
            pos_t m2_pos = make_pos_t(m2.nodes.front());
            int max_length = max(read1.sequence().size(), read2.sequence().size());
            double overlap_length = mems_overlap_length(m1, m2);
            pair<int64_t, int64_t> distances = mem_min_oriented_distances(m1, m2);
            int64_t dist_fwd = distances.first; // use only the forward orientation
            //int64_t dist_inv = distances.second;
            if (dist_fwd > max_length) {
                // too far
                return -std::numeric_limits<double>::max();
            } else {
                // accepted transition
                int read_dist = m2.begin - m1.begin;
                double jump_fwd = abs(read_dist - dist_fwd);
                //double jump_inv = abs(read_dist - dist_inv);
                return (double) -(gap_open*(jump_fwd>0) + jump_fwd*gap_extension) -overlap_length;
            }
        }
    };

    // build the paired-read MEM markov model
    vector<vector<MaximalExactMatch> > clusters;
    if (total_multimaps) {
        // We're going to run the chainer because we want to calculate alignments
        
        // What band width during the alignment should the chainer plan for?
        int64_t band_width = frag_stats.fragment_max;
        
#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) {
                cerr << "Invoking MEM chainer with band width " << band_width
                << ", fragment size " << frag_stats.fragment_size << ", frag_stats.fragment_max " 
                << frag_stats.fragment_max << endl;
            }
        }
#endif
    
        MEMChainModel chainer({ read1.sequence().size(), read2.sequence().size() },
                              { mems1, mems2 },
                              [&](pos_t n) -> int64_t {
                                  return approx_position(n);
                              },
                              [&](pos_t n) -> unordered_map<path_handle_t, vector<pair<size_t, bool> > > {
                                  return algorithms::nearest_offsets_in_paths(xindex, n, -1);
                              },
                              transition_weight,
                              band_width);
        clusters = chainer.traceback(total_multimaps, false, debug);
    }

    auto show_clusters = [&](void) {
        cerr << "clusters: " << endl;
        for (auto& cluster : clusters) {
            cerr << cluster.size() << " MEMs covering " << cluster_coverage(cluster) << " @ ";
            for (auto& mem : cluster) {
                size_t len = mem.begin - mem.end;
                for (auto& node : mem.nodes) {
                    id_t id = gcsa::Node::id(node);
                    size_t offset = gcsa::Node::offset(node);
                    bool is_rev = gcsa::Node::rc(node);
                    cerr << "|" << id << (is_rev ? "-" : "+") << ":" << offset << "," << mem.fragment << ",";
                }
                cerr << mem.sequence() << " ";
            }
            cerr << endl;
        }
    };

    if (debug) {
        cerr << "### clusters before filtering:" << endl;
        show_clusters();
    }

    vector<vector<MaximalExactMatch> > clusters1;
    vector<vector<MaximalExactMatch> > clusters2;
    for (auto& cluster : clusters) {
        // break the clusters into their fragments
        clusters1.emplace_back();
        auto& cluster1 = clusters1.back();
        clusters2.emplace_back();
        auto& cluster2 = clusters2.back();
        bool seen1=false, seen2=false;
        for (auto& mem : cluster) {
            if (!seen2 && mem.fragment == 1) {
                cluster1.push_back(mem);
                seen1 = true;
            } else if (mem.fragment == 2) {
                cluster2.push_back(mem);
                seen2 = true;
            } else {
                cerr << "vg map error misordered fragments in cluster" << endl;
                assert(false);
            }
        }
    }

    auto to_drop1 = clusters_to_drop(clusters1);
    auto to_drop2 = clusters_to_drop(clusters2);
    vector<pair<Alignment, Alignment> > alns;
    vector<pair<vector<MaximalExactMatch>*, vector<MaximalExactMatch>*> > cluster_ptrs;
    for (int i = 0; i < clusters1.size(); ++i) {
        auto& cluster1 = clusters1[i];
        auto& cluster2 = clusters2[i];
        cluster_ptrs.push_back(make_pair(&cluster1, &cluster2));
    }

    // sort the clusters by unique coverage in the read
    map<pair<vector<MaximalExactMatch>*, vector<MaximalExactMatch>*>, int> cluster_cov;
    for (auto& p : cluster_ptrs) {
        cluster_cov[p] = cluster_coverage(*p.first) + cluster_coverage(*p.second);
    }
    sort(cluster_ptrs.begin(),
         cluster_ptrs.end(), [&](const pair<vector<MaximalExactMatch>*, vector<MaximalExactMatch>*>& p1,
                                 const pair<vector<MaximalExactMatch>*, vector<MaximalExactMatch>*>& p2) {
             return cluster_cov[p1] > cluster_cov[p2];
         });

    auto show_paired_clusters = [&](void) {
        cerr << "clusters: " << endl;
        for (auto& cluster_ptr : cluster_ptrs) {
            auto& cluster1 = *cluster_ptr.first;
            auto& cluster2 = *cluster_ptr.second;
            cerr << cluster1.size() << " " << cluster2.size() << " MEMs covering " << cluster_coverage(cluster1) << " " << cluster_coverage(cluster2) << " @ ";
            cerr << " cluster1: ";
            for (auto& mem : cluster1) {
                size_t len = mem.begin - mem.end;
                for (auto& node : mem.nodes) {
                    id_t id = gcsa::Node::id(node);
                    size_t offset = gcsa::Node::offset(node);
                    bool is_rev = gcsa::Node::rc(node);
                    cerr << "|" << id << (is_rev ? "-" : "+") << ":" << offset << "," << mem.fragment << ",";
                }
                cerr << mem.sequence() << " ";
            }
            cerr << " cluster2: ";
            for (auto& mem : cluster2) {
                size_t len = mem.begin - mem.end;
                for (auto& node : mem.nodes) {
                    id_t id = gcsa::Node::id(node);
                    size_t offset = gcsa::Node::offset(node);
                    bool is_rev = gcsa::Node::rc(node);
                    cerr << "|" << id << (is_rev ? "-" : "+") << ":" << offset << "," << mem.fragment << ",";
                }
                cerr << mem.sequence() << " ";
            }
            cerr << endl;
        }
    };

    if (debug) {
        cerr << "### clusters after filtering:" << endl;
        show_paired_clusters();
    }

    set<pair<string, string> > seen_alignments;
    int filled1 = 0, filled2 = 0;
    for (auto& cluster_ptr : cluster_ptrs) {
        // break the cluster into two pieces
        auto& cluster1 = *cluster_ptr.first;
        auto& cluster2 = *cluster_ptr.second;
        alns.emplace_back();
        auto& p = alns.back();
        if (cluster1.size()
            && (filled1 < min_multimaps
                || !to_drop1.count(&cluster1))) {
            p.first = align_cluster(read1, cluster1, true, xdrop_alignment);
            ++filled1;
        } else {
            p.first = read1;
            p.first.clear_score();
            p.first.clear_identity();
            p.first.clear_path();
        }

        if (cluster2.size()
            && (filled2 < min_multimaps
                || !to_drop2.count(&cluster2))) {
            p.second = align_cluster(read2, cluster2, true, xdrop_alignment);
            ++filled2;
        } else {
            p.second = read2;
            p.second.clear_score();
            p.second.clear_identity();
            p.second.clear_path();
        }
        auto pair_sig = signature(p.first, p.second);
        if (seen_alignments.count(pair_sig)) {
            alns.pop_back();
            alns.emplace_back();
        } else {
            seen_alignments.insert(pair_sig);
        }
    }
    assert(cluster_ptrs.size() == alns.size());

    vector<pair<Alignment, Alignment>*> aln_ptrs;
    map<pair<Alignment, Alignment>*, int> aln_index;
    auto update_aln_ptrs = [&](void) {
        aln_ptrs.clear();
        aln_index.clear();
        int idx = 0;
        for (auto& alnp : alns) {
            aln_ptrs.push_back(&alnp);
            aln_index[&alnp] = idx++;
        }
    };

    update_aln_ptrs();

    auto show_alignments = [&](const string& arg) {
        if (debug) {
            for (auto& p : aln_ptrs) {
                cerr << arg << " ";
                auto& aln1 = p->first;
                cerr << "1:" << aln1.score();
                if (aln1.score()) cerr << "@" << aln1.path().mapping(0).position().node_id() << " ";
                auto& aln2 = p->second;
                cerr << " 2:" << aln2.score();
                if (aln2.score()) cerr << "@" << aln2.path().mapping(0).position().node_id() << " ";
                cerr << endl;
            }
        }
    };

    auto sort_and_dedup = [&](void) {
        // apply the fragment lengths for faster sorting
        for (auto& p : aln_ptrs) {
            auto& aln1 = p->first;
            auto& aln2 = p->second;
            auto approx_frag_lengths = min_pair_fragment_length(aln1, aln2);
            frag_stats.save_frag_lens_to_alns(aln1, aln2, approx_frag_lengths, xindex, pair_consistent(aln1, aln2, 1e-6));
        }
        // sort the aligned pairs by score
        sort_shuffling_ties(aln_ptrs.begin(), aln_ptrs.end(),
                            [&](pair<Alignment, Alignment>* pair1,
                                pair<Alignment, Alignment>* pair2) {
                                double weight1=0, weight2=0;
                                if (frag_stats.fragment_size) {
                                    weight1 = pair1->first.fragment_score();
                                    weight2 = pair2->first.fragment_score();
                                }
                                double score1 = (pair1->first.score() + pair1->second.score());
                                double score2 = (pair2->first.score() + pair2->second.score());
                                return score1 + weight1 > score2 + weight2;
                            });
        seen_alignments.clear();
        // remove duplicates (same score and same start position of both pairs)
        aln_ptrs.erase(
            std::remove_if(aln_ptrs.begin(), aln_ptrs.end(),
                           [&](pair<Alignment, Alignment>* p) {
                               auto pair_sig = signature(p->first, p->second);
                               if (seen_alignments.count(pair_sig)) {
                                   return true;
                               } else {
                                   seen_alignments.insert(pair_sig);
                                   return false;
                               }
                           }),
            aln_ptrs.end());
    };
    show_alignments("raw");

    sort_and_dedup();
    show_alignments("dedup");

    // now add back in single-ended versions of everything that we can use for rescue
    vector<pair<Alignment, Alignment> > se_alns;
    vector<pair<vector<MaximalExactMatch>*, vector<MaximalExactMatch>*> > se_cluster_ptrs;
    double max_possible_score = read1.sequence().size() * match + 2*full_length_bonus;
    double hang_threshold = max_possible_score * pair_rescue_hang_threshold;
    double retry_threshold = max_possible_score * pair_rescue_retry_threshold;
    for (auto& p : aln_ptrs) {
        auto& aln1 = p->first;
        auto& aln2 = p->second;
        auto cluster_ptr = cluster_ptrs[aln_index[p]];
        bool consistent = aln1.score() && aln2.score() && pair_consistent(aln1, aln2, 1e-6);
        // if both mates are aligned, add each single end into the mix
        if (aln1.score() > hang_threshold && (aln2.score() <= retry_threshold || !consistent)) {
            se_alns.emplace_back();
            auto& p = se_alns.back();
            p.first = aln1;
            p.second = read2;
            se_cluster_ptrs.push_back(make_pair(cluster_ptr.first, nullptr));
        }
        if (aln2.score() > hang_threshold && (aln1.score() <= retry_threshold || !consistent)) {
            se_alns.emplace_back();
            auto& q = se_alns.back();
            q.first = read1;
            q.second = aln2;
            se_cluster_ptrs.push_back(make_pair(nullptr, cluster_ptr.second));
        }
    }
    int k = 0;
    for (auto& se_aln : se_alns) {
        alns.push_back(se_aln);
        cluster_ptrs.push_back(se_cluster_ptrs[k]);
        ++k;
    }
    update_aln_ptrs();
    sort_and_dedup();
    show_alignments("singles");

    map<Alignment*, bool> rescued_aln;
    if (mate_rescues && frag_stats.fragment_size) {
        // go through the pairs and see if we need to rescue one side off the other
        bool rescued = false;
        int j = 0;
        for (auto& p : aln_ptrs) {
            if (j > mate_rescues) break;
            auto& aln1 = p->first;
            auto& aln2 = p->second;
            int score1 = aln1.score();
            int score2 = aln2.score();
            bool tried1 = false; bool tried2 = false;
            pair<bool, bool> rescues = pair_rescue(aln1, aln2, tried1, tried2, match, full_length_bonus, true, xdrop_alignment);
            if (tried1) ++j;
            if (tried2) ++j;
            rescued_aln[&aln1] = rescues.first;
            rescued_aln[&aln2] = rescues.second;
            rescued |= rescues.first || rescues.second;
        }
        if (rescued) {
            sort_and_dedup();
            show_alignments("rescue");
        }
    }

    vector<pair<Alignment, Alignment> > pe_alns;
    vector<pair<vector<MaximalExactMatch>*, vector<MaximalExactMatch>*> > pe_cluster_ptrs;
    double mix_threshold = max_possible_score * 0.1;
    int x = 0;
    for (auto& p1 : aln_ptrs) {
        //if (++x > total_multimaps) break;
        auto& aln1a = p1->first;
        auto& aln2a = p1->second;
        auto cluster_ptrA = cluster_ptrs[aln_index[p1]];
        int y = 0;
        for (auto& p2 : aln_ptrs) {
            //if (++y > total_multimaps) break;
            if (p1 == p2) continue;
            auto& aln1b = p2->first;
            auto& aln2b = p2->second;
            auto cluster_ptrB = cluster_ptrs[aln_index[p2]];
            bool consistent1a2b = aln1a.score() > mix_threshold && aln2b.score() > mix_threshold && pair_consistent(aln1a, aln2b, 1e-6);
            bool consistent1b2a = aln1b.score() > mix_threshold && aln2a.score() > mix_threshold && pair_consistent(aln1b, aln2a, 1e-6);
            if (consistent1a2b) {
                pe_alns.emplace_back();
                auto& p = pe_alns.back();
                p.first = aln1a;
                p.second = aln2b;
                pe_cluster_ptrs.push_back(make_pair(cluster_ptrA.first, cluster_ptrB.second));
            }
            if (consistent1b2a) {
                pe_alns.emplace_back();
                auto& p = pe_alns.back();
                p.first = aln1b;
                p.second = aln2a;
                pe_cluster_ptrs.push_back(make_pair(cluster_ptrB.first, cluster_ptrA.second));
            }
        }
    }
    k = 0;
    for (auto& pe_aln : pe_alns) {
        alns.push_back(pe_aln);
        cluster_ptrs.push_back(pe_cluster_ptrs[k]);
        ++k;
    }
    update_aln_ptrs();

    // Apply haplotype consistency scores if possible
    vector<Alignment*> flat_alns;
    flat_alns.reserve(aln_ptrs.size() * 2);
    for (auto& aln_pair : aln_ptrs) {
        flat_alns.push_back(&aln_pair->first);
        flat_alns.push_back(&aln_pair->second);
    }
    apply_haplotype_consistency_scores(flat_alns);

    sort_and_dedup();
    show_alignments("mixed");

    // apply pair consistency penalty
    if (frag_stats.fragment_size) {
        for (auto& p : aln_ptrs) {
            auto& aln1 = p->first;
            auto& aln2 = p->second;
            if (!pair_consistent(aln1, aln2, 0)) {
                aln1.set_score(max(0, aln1.score()-unpaired_penalty));
                aln2.set_score(max(0, aln2.score()-unpaired_penalty));
            }
        }
    }

    // calculate cluster mapping quality

    if (use_cluster_mq) {
        cluster_mq = compute_cluster_mapping_quality(clusters, read1.sequence().size() + read2.sequence().size());
#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) cerr << "cluster mq == " << cluster_mq << endl;
        }
#endif
    }

    show_alignments("end");
    
    int read1_max_score = 0;
    int read2_max_score = 0;
    int possible_pairs = 0;

    // build up the results
    int i = 0;
    for (auto& p : aln_ptrs) {
        read1_max_score = max(p->first.score(), read1_max_score);
        read2_max_score = max(p->second.score(), read2_max_score);
        results.first.push_back(p->first);
        results.second.push_back(p->second);
        possible_pairs += p->first.score() > 0 && p->second.score() > 0;
    }
    bool max_first = results.first.size() && (read1_max_score == results.first.front().score() && read2_max_score == results.second.front().score());

    // compute mapping quality
    compute_mapping_qualities(results, cluster_mq, maybe_mq1, maybe_mq2, mq_cap1, mq_cap2);

    // remove the extra pair used to compute mapping quality if necessary
    if (results.first.size() > max_multimaps) {
        results.first.resize(max_multimaps);
        results.second.resize(max_multimaps);
    }

    // mark primary and secondary
    for (int i = 0; i < results.first.size(); i++) {
        results.first[i].mutable_fragment_next()->set_name(read2.name());
        results.first[i].set_is_secondary(i > 0);
        results.second[i].mutable_fragment_prev()->set_name(read1.name());
        results.second[i].set_is_secondary(i > 0);
    }

    // optionally zap everything unless the primary alignments are individually top-scoring
    if (only_top_scoring_pair && results.first.size() &&
        (results.first[0].score() < read1_max_score ||
         results.second[0].score() < read2_max_score)) {
        results.first.clear();
        results.second.clear();
    }
    
    // Before we update the fragment length distribution, remember the
    // distribution used here.
    
    stringstream fragment_dist;
    fragment_dist << frag_stats.fragment_size
        << ':' << frag_stats.cached_fragment_length_mean
        << ':' << frag_stats.cached_fragment_length_stdev 
        << ':' << frag_stats.cached_fragment_orientation_same
        << ':' << frag_stats.cached_fragment_direction;
    
    // we tried to align
    // if we don't have a fragment_size yet determined
    // and we didn't get a perfect, unambiguous hit on both reads
    // we'll need to try it again later when we do have a fragment_size
    // so store it in a buffer local to this mapper

    // tag the results with their fragment lengths
    // record the lengths in a deque that we use to keep a running estimate of the fragment length distribution
    // we then set the fragment_size cutoff using the moments of the estimated distribution
    bool imperfect_pair = false;
    for (int i = 0; i < min(results.first.size(), results.second.size()); ++i) {
        auto& aln1 = results.first.at(i);
        auto& aln2 = results.second.at(i);
        //double ident1 = (double) aln1.score() / max_possible_score;
        //double ident2 = (double) aln2.score() / max_possible_score;
        int64_t length = std::numeric_limits<int64_t>::max();
        for (int j = 0; j < aln1.fragment_size(); ++j) {
            length = min(length, abs(aln1.fragment(j).length()));
        }
        bool consistent = ((frag_stats.fragment_size && length < frag_stats.fragment_size && pair_consistent(aln1, aln2, 1e-3))
                           || (!frag_stats.fragment_size && length < frag_stats.fragment_max));
        if (!retrying
            && consistent
            && results.first.size() == 1
            && results.second.size() == 1
            && results.first.front().identity() > frag_stats.perfect_pair_identity_threshold
            && results.second.front().identity() > frag_stats.perfect_pair_identity_threshold) { // hard cutoff
            //cerr << "aln\tperfect alignments" << endl;
            frag_stats.record_fragment_configuration(aln1, aln2, this);
        } else if (!retrying && !frag_stats.fragment_size) {
            // mark this pair to be buffered and remapped
            imperfect_pair = true;
        }
        
        if (consistent) {
            set_annotation(aln1, "proper_pair", true);
            set_annotation(aln2, "proper_pair", true);
        }
        else {
            set_annotation(aln1, "proper_pair", false);
            set_annotation(aln2, "proper_pair", false);
        }
    }

    if (!retrying && imperfect_pair && frag_stats.fragment_max) {
        imperfect_pairs_to_retry.push_back(make_pair(read1, read2));
        results.first.clear();
        results.second.clear();
        // we signal the fact that this isn't a perfect pair, so we don't write it out externally
        queued_resolve_later = true;
    }

    if(results.first.empty() && !exclude_unaligned) {
        results.first.push_back(read1);
        auto& aln = results.first.back();
        aln.clear_path();
        aln.clear_score();
        aln.clear_identity();
    }
    if(results.second.empty() && !exclude_unaligned) {
        results.second.push_back(read2);
        auto& aln = results.second.back();
        aln.clear_path();
        aln.clear_score();
        aln.clear_identity();
    }
    
    // Make sure to link up alignments even if they aren't mapped.
    for (auto& aln : results.first) {
        aln.set_name(read1.name());
        aln.mutable_fragment_next()->set_name(read2.name());
        aln.set_sequence(read1.sequence());
        aln.set_quality(read1.quality());
        aln.set_fragment_length_distribution(fragment_dist.str());
    }

    for (auto& aln : results.second) {
        aln.set_name(read2.name());
        aln.mutable_fragment_prev()->set_name(read1.name());
        aln.set_sequence(read2.sequence());
        aln.set_quality(read2.quality());
        aln.set_fragment_length_distribution(fragment_dist.str());
    }

    if (!results.first.empty() && !results.second.empty()) {
        // if we have references, annotate the alignments with their reference positions
        algorithms::annotate_with_initial_path_positions(*xindex, results.first);
        algorithms::annotate_with_initial_path_positions(*xindex, results.second);

        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        auto used_time = chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
        results.first.front().set_time_used(used_time);
        results.second.front().set_time_used(used_time);
    }

    return results;

}

double Mapper::compute_cluster_mapping_quality(const vector<vector<MaximalExactMatch> >& clusters,
                                               int read_length) {
    if (clusters.size() == 0) {
        return 0;
    }
    if (clusters.size() == 1) {
        return (double) max_cluster_mapping_quality;
    }
    vector<double> weights;
    for (auto& cluster : clusters) {
        weights.emplace_back();
        double& weight = weights.back();
        for (int i = 0; i < cluster.size(); ++i) {
            // for each mem, count half of its coverage with its neighbors towards this metric
            auto& mem = cluster[i];
            int shared_coverage = 0;
            if (i > 0) {
                auto& prev = cluster[i-1];
                if (prev.fragment == mem.fragment) {
                    shared_coverage += (prev.end <= mem.begin ? 0 : prev.end - mem.begin);
                }
            }
            if (i < cluster.size()-1) {
                auto& next = cluster[i+1];
                if (next.fragment == mem.fragment) {
                    shared_coverage += (mem.end <= next.begin ? 0 : mem.end - next.begin);
                }
            }
            weight +=
                (((double)mem.length() - (double)shared_coverage/2)
                 / read_length)
                / mem.match_count;
        }
        //cerr << "weight " << weight << endl;
    }
    // return the ratio between best and second best as quality
    std::sort(weights.begin(), weights.end(), std::greater<double>());
    // find how many maxes we have
    double max_weight = weights.front();
    int max_count = 0;
    while (max_weight == weights[max_count]) ++max_count;
    double best_chance = max_count > 1 ? prob_to_phred(1.0-(1.0/max_count)) : 0;
    if (weights[0] == 0) return 0;
    return min((double)max_cluster_mapping_quality,
               max(best_chance, prob_to_phred(weights[1]/weights[0])));
}

int sub_overlaps_of_first_aln(const vector<Alignment>& alns, float overlap_fraction) {
    // take the first
    // now look at the rest and measure overlap
    if (alns.empty()) return 0;
    auto& aln1 = alns.front();
    int seq_len = aln1.sequence().size();
    int overlaps = 0;
    for (int i = 1; i < alns.size(); ++i) {
        auto& aln2 = alns[i];
        // where the overlap is greater than overlap_fraction
        if ((double)query_overlap(aln1, aln2)/seq_len >= overlap_fraction) {
            ++overlaps;
        }
    }
    return overlaps;
}

set<const vector<MaximalExactMatch>* > Mapper::clusters_to_drop(const vector<vector<MaximalExactMatch> >& clusters) {
    set<const vector<MaximalExactMatch>* > to_drop;
    vector<double> covs; covs.resize(clusters.size());
    for (int i = 0; i < clusters.size(); ++i) {
        covs[i] = cluster_coverage(clusters[i]);
    }
    for (int i = 0; i < clusters.size(); ++i) {
        // establish overlaps with longer clusters for all clusters
        auto& this_cluster = clusters[i];
        double t = covs[i];
        if (t < min_cluster_length) {
            to_drop.insert(&this_cluster);
            continue;
        }
        int b = -1;
        int l = t;
        for (int j = i; j >= 0; --j) {
            if (j == i) continue;
            // are we overlapping?
            auto& other_cluster = clusters[j];
            if (to_drop.count(&other_cluster)) continue;
            // if the overlap length is more than our drop_chain fraction
            if (drop_chain > 0
                && clusters_overlap_in_read(this_cluster, other_cluster)) {
                double length_ratio = t/covs[j];
                if (length_ratio < drop_chain) {
                    to_drop.insert(&this_cluster);
                    break;
                }
            }
        }
    }
    return to_drop;
}

vector<Alignment>
Mapper::align_mem_multi(const Alignment& aln,
                        vector<MaximalExactMatch>& mems,
                        double& cluster_mq,
                        double longest_lcp,
                        double fraction_filtered,
                        int max_mem_length,
                        int keep_multimaps,
                        int additional_multimaps,
                        bool xdrop_alignment) {

    if (debug) cerr << "mems for read " << mems_to_json(mems) << endl;
    
    auto aligner = get_aligner(!aln.quality().empty());
    int8_t match = aligner->match;
    int8_t gap_extension = aligner->gap_extension;
    int8_t gap_open = aligner->gap_open;

    int total_multimaps = max(max_multimaps, additional_multimaps);
    double mq_cap = max_mapping_quality;

    int total_mem_length = 0;
    for (auto& mem : mems) total_mem_length += mem.length(); // * mem.nodes.size();
    double mem_read_ratio = min(1.0, (double)total_mem_length / (double)aln.sequence().size());
    int basis_length = min((int)aln.sequence().size(), (int)gcsa->order());
    double max_possible_mq = max_possible_mapping_quality(basis_length);

    // Estimate the maximum mapping quality we can get if the alignments based on the good MEMs are the best ones.
    int mem_max_length = 0;
    for (auto& mem : mems) if (mem.primary && mem.match_count) mem_max_length = max(mem_max_length, (int)mem.length());
    double maybe_mq = 0;
    if (mems.size()) {
        maybe_mq = estimate_max_possible_mapping_quality(basis_length,
                                                         basis_length/max(1.0, (double)mem_max_length),
                                                         basis_length/longest_lcp);
    }

    // scale difficulty using the estimated mapping quality
    /*
    if (maybe_mq > max_mapping_quality) {
        total_multimaps = max(max(min_multimaps*4, max_multimaps), min(total_multimaps, (int)round(maybe_mq)));
    }
    */

    if (debug) cerr << "maybe_mq " << aln.name() << " " << maybe_mq << " " << total_multimaps << " " << mem_max_length << " " << longest_lcp << " " << total_multimaps << " " << mem_read_ratio << " " << fraction_filtered << " " << max_possible_mq << " " << total_multimaps << endl;

    // go through the ordered single-hit MEMs
    // build the clustering model
    // find the alignments that are the best-scoring walks through it
    auto transition_weight = [&](const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
        pos_t m1_pos = make_pos_t(m1.nodes.front());
        pos_t m2_pos = make_pos_t(m2.nodes.front());
        int64_t max_length = aln.sequence().size();
        double overlap_length = mems_overlap_length(m1, m2);
        pair<int64_t, int64_t> distances = mem_min_oriented_distances(m1, m2);
        int64_t dist_fwd = distances.first; // use only the forward orientation
        //int64_t dist_inv = distances.second;
        if (dist_fwd > max_length) {
            // too far
            return -std::numeric_limits<double>::max();
        } else {
            // accepted transition
            int read_dist = m2.begin - m1.begin;
            double jump_fwd = abs(read_dist - dist_fwd);
            //double jump_inv = abs(read_dist - dist_inv);
            return (double) -(gap_open*(jump_fwd>0) + jump_fwd*gap_extension) -overlap_length;
        }
    };

    // establish the chains
    vector<vector<MaximalExactMatch> > clusters;
    if (total_multimaps) {
        MEMChainModel chainer({ aln.sequence().size() }, { mems },
                              [&](pos_t n) -> int64_t {
                                  return approx_position(n);
                              },
                              [&](pos_t n) -> unordered_map<path_handle_t, vector<pair<size_t, bool> > > {
                                  return algorithms::nearest_offsets_in_paths(xindex, n, -1);
                              },
                              transition_weight,
                              aln.sequence().size());
        clusters = chainer.traceback(total_multimaps, false, debug);
    }
    
    /*
    map<const vector<MaximalExactMatch>*, int> cluster_cov;
    for (auto& cluster : clusters) {
        cluster_cov[&cluster] = cluster_coverage(cluster);
    }
    */
    
    // don't attempt to align if we reach the maximum number of multimaps
    //if (clusters.size() == total_multimaps) clusters.clear();

    auto show_clusters = [&](void) {
        cerr << "clusters: " << endl;
        for (auto& cluster : clusters) {
            cerr << cluster.size() << " MEMs covering " << cluster_coverage(cluster) << " @ ";
            for (auto& mem : cluster) {
                size_t len = mem.begin - mem.end;
                for (auto& node : mem.nodes) {
                    id_t id = gcsa::Node::id(node);
                    size_t offset = gcsa::Node::offset(node);
                    bool is_rev = gcsa::Node::rc(node);
                    cerr << "|" << id << (is_rev ? "-" : "+") << ":" << offset << "," << mem.fragment << ",";
                }
                cerr << mem.sequence() << " ";
            }
            cerr << endl;
        }
    };

    if (debug) {
        cerr << "### clusters:" << endl;
        show_clusters();
    }

    if (use_cluster_mq) {
        cluster_mq = compute_cluster_mapping_quality(clusters, aln.sequence().size());
    }
#ifdef debug_mapper
#pragma omp critical
    {
        if (debug) {
            cerr << "cluster mapping quality " << cluster_mq << endl;
        }
    }
#endif
    auto to_drop = clusters_to_drop(clusters);

    // for up to our required number of multimaps
    // make the perfect-match alignment for the SMEM cluster
    // then fix it up with DP on the little bits between the alignments
    vector<Alignment> alns;
    vector<vector<MaximalExactMatch>*> used_clusters;
    set<string> seen_alignments;
    int multimaps = 0;
    int filled = 0;
    for (auto& cluster : clusters) {
        if (alns.size() >= total_multimaps) { break; }
        // skip if we've filtered the cluster
        if (to_drop.count(&cluster) && filled >= min_multimaps) {
            alns.push_back(aln);
            used_clusters.push_back(&cluster);
            continue;
        }
        ++filled;
        Alignment candidate = align_cluster(aln, cluster, true, xdrop_alignment);
        string sig = signature(candidate);

#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) {
                cerr << "Alignment with signature " << sig << " (seen: " << seen_alignments.count(sig) << ")" << endl;
                cerr << "\t" << pb2json(candidate) << endl;
            }
        }
#endif

        if (!seen_alignments.count(sig)) {
            alns.push_back(candidate);
            used_clusters.push_back(&cluster);
            seen_alignments.insert(sig);
        }
    }
    
    assert(alns.size() == used_clusters.size());
#ifdef debug_mapper
#pragma omp critical
    if (debug) {
        cerr << "alignments(" << alns.size() << ")" << endl;
        for (auto& aln : alns) {
            cerr << aln.score();
            if (aln.score()) cerr << " pos1 " << aln.path().mapping(0).position().node_id() << " ";
            cerr << endl;
        }
    }
#endif

    // Prepare a sortable vector of alignment pointers
    vector<Alignment*> aln_ptrs;
    map<Alignment*, int> aln_index;
    int idx = 0;
    for (auto& aln : alns) {
        aln_ptrs.push_back(&aln);
        aln_index[&aln] = idx++;
    }
    
    // Apply haplotype consistency scoring if possible
    apply_haplotype_consistency_scores(aln_ptrs);
    
    // sort alignments by score
    sort_shuffling_ties(aln_ptrs.begin(), aln_ptrs.end(),
                        [&](Alignment* a1, Alignment* a2) {
                            return a1->score() > a2->score();
                        });
    // remove likely perfect duplicates
    aln_ptrs.erase(
        std::unique(
            aln_ptrs.begin(), aln_ptrs.end(),
            [&](Alignment* aln1,
                Alignment* aln2) {
                return
                    aln1->score() == aln2->score()
                    && (aln1->score() == 0
                        || make_pos_t(aln1->path().mapping(0).position())
                        == make_pos_t(aln2->path().mapping(0).position()));
            }),
        aln_ptrs.end());
    if (aln_ptrs.size()) {
        vector<Alignment> best_alns;
        int i = 0;
        for ( ; i < min((int)aln_ptrs.size(), keep_multimaps); ++i) {
            Alignment* alnp = aln_ptrs.at(i);
            best_alns.push_back(*alnp);
        }
        for ( ; i < aln_ptrs.size(); ++i) {
            best_alns.push_back(*aln_ptrs[i]);
        }
        alns = score_sort_and_deduplicate_alignments(best_alns, aln);
    }
    // compute the mapping quality
    compute_mapping_qualities(alns, cluster_mq, maybe_mq, mq_cap);

    // final filter step
    filter_and_process_multimaps(alns, keep_multimaps);

    // if we didn't get anything, return an unaligned version of our input
    if (alns.empty() && !exclude_unaligned) {
        alns.push_back(aln);
        auto& unaligned = alns.back();
        unaligned.clear_path();
        unaligned.clear_score();
        unaligned.clear_identity();
    }

    return alns;
}

Alignment Mapper::align_maybe_flip(const Alignment& base, HandleGraph& graph, bool flip, bool traceback, bool banded_global, bool xdrop_alignment) {
    // do not use X-drop alignment when seed position is not available
    vector<MaximalExactMatch> mems;
    return(align_maybe_flip(base, graph, mems, flip, traceback, banded_global, xdrop_alignment));
}

Alignment Mapper::align_maybe_flip(const Alignment& base, HandleGraph& graph, const vector<MaximalExactMatch>& mems, bool flip, bool traceback, bool banded_global, bool xdrop_alignment) {
    Alignment aln = base;
    // do not modify aln.sequence() in X-drop alignment so that seed position is calculated by mem.begin - aln.sequence().begin()
    aln.set_sequence(base.sequence());
    if (!base.quality().empty()) {
        aln.set_quality(base.quality());
    }
    bool pinned_alignment = false;
    bool pinned_reverse = false;
    aln = align_to_graph(aln,
                         graph,
                         mems,
                         flip,
                         traceback,
                         pinned_alignment,
                         pinned_reverse,
                         banded_global,
                         xdrop_alignment ? 1 : 0,
                         include_full_length_bonuses);
    if (strip_bonuses && !banded_global && traceback) {
        // We want to remove the bonuses
        aln.set_score(get_aligner()->remove_bonuses(aln));
    }
    return aln;
}

double Mapper::compute_uniqueness(const Alignment& aln, const vector<MaximalExactMatch>& mems) {
    // compute the per-base copy number of the alignment based on the MEMs in the cluster
    vector<int> v; v.resize(aln.sequence().size());
    auto aln_begin = aln.sequence().begin();
    for (auto& mem : mems) {
        // from the beginning to the end of the mem in the read
        // add the count of hits to each position
        for (int i = mem.begin - aln_begin; i < mem.end - aln_begin; ++i) {
            v[i] += mem.match_count;
        }
    }
    // what is the fraction of the cluster that is in repetitive mems
    // we calculate the uniqueness metric as the average number of hits for bases in the read covered by any MEM
    double repeated = std::accumulate(v.begin(), v.end(), 0.0,
                                      [](const int& a, const int& b) { return b > 1 ? a + 1 : a; });
    return repeated / aln.sequence().length();
}

Alignment Mapper::align_cluster(const Alignment& aln, const vector<MaximalExactMatch>& mems, bool traceback, bool xdrop_alignment) {
    bdsg::HashGraph graph = cluster_subgraph_containing(*xindex, aln, mems, get_aligner());
    Alignment aligned = align_maybe_flip(aln, graph, mems, false, traceback, false, xdrop_alignment);
    return aligned;
}

// estimate the fragment length as the difference in mean positions of both alignments
unordered_map<path_handle_t, int64_t> Mapper::min_pair_fragment_length(const Alignment& aln1, const Alignment& aln2) {
    unordered_map<path_handle_t, int64_t> lengths;
    auto pos1 = algorithms::alignment_path_offsets(*xindex, aln1, true, false);
    auto pos2 = algorithms::alignment_path_offsets(*xindex, aln2, true, false);
    for (auto& p : pos1) {
        auto x = pos2.find(p.first);
        if (x != pos2.end()) {
            auto& d = lengths[p.first];
            int64_t l = std::numeric_limits<int64_t>::max();
            for (auto& pos1 : p.second) {
                for (auto& pos2 : x->second) {
                    int64_t f = (int64_t)pos2.first - (int64_t)pos1.first;
                    if (abs(f) < l) {
                        d = f; l = abs(f);
                    }
                }
            }
        }
    }
    //cerr << "got lengths "; for (auto& c : lengths) cerr << c.first << ":" << c.second << " "; cerr << endl;
    return lengths;
}

string FragmentLengthStatistics::fragment_model_str(void) {
    stringstream s;
    s << fragment_size << ":"
      << cached_fragment_length_mean << ":"
      << cached_fragment_length_stdev << ":"
      << cached_fragment_orientation_same << ":"
      << cached_fragment_direction;
    return s.str();
}

void FragmentLengthStatistics::save_frag_lens_to_alns(Alignment& aln1, Alignment& aln2,
                                                      const unordered_map<path_handle_t, int64_t>& approx_frag_lengths,
                                                      PathPositionHandleGraph* xindex,
                                                      bool is_consistent) {
    double max_score = 0;
    aln1.clear_fragment();
    aln2.clear_fragment();
    for (auto& j : approx_frag_lengths) {
        Path fragment;
        fragment.set_name(xindex->get_path_name(j.first));
        int length = j.second;
        fragment.set_length(length);
        *aln1.add_fragment() = fragment;
        *aln2.add_fragment() = fragment;
        if (fragment_size && is_consistent) {
            double score = 50 + max(1.0, log(min(100.0, prob_to_phred(1.0-fragment_length_pval(abs(length))))));
            max_score = max(max_score, score);
        }
    }
    aln1.set_fragment_score(max_score);
    aln2.set_fragment_score(max_score);
}

// XXX TODO this is busted because it doesn't respect the actual paths in the graph and assumes the alignment mapping position provides relative direction
void FragmentLengthStatistics::record_fragment_configuration(const Alignment& aln1, const Alignment& aln2, Mapper* mapper) {
    if (fixed_fragment_model) return;
    assert(aln1.path().mapping(0).has_position() && aln2.path().mapping(0).has_position());    
    unordered_map<path_handle_t, tuple<int64_t, bool, bool> > lengths;
    auto pos1 = algorithms::alignment_path_offsets(*mapper->xindex, aln1, true, false);
    auto pos2 = algorithms::alignment_path_offsets(*mapper->xindex, aln2, true, false);
    for (auto& p : pos1) {
        auto x = pos2.find(p.first);
        if (x != pos2.end()) {
            auto& d = lengths[p.first];
            int64_t l = std::numeric_limits<int64_t>::max();
            for (auto& pos1 : p.second) {
                for (auto& pos2 : x->second) {
                    int64_t f = (int64_t)pos2.first - (int64_t)pos1.first;
                    if (abs(f) < l) {
                        l = abs(f);
                        d = make_tuple(f, pos1.second, pos2.second);
                    }
                }
            }
        }
    }
    //cerr << "got lengths "; for (auto& c : lengths) cerr << c.first << ":" << c.second << " "; cerr << endl;
    //return lengths;
    for (auto& chr : lengths) {
        int64_t length = get<0>(chr.second);
        if (abs(length) > fragment_max // keep out super high values
            || (fragment_size && abs(length) > fragment_size)) continue;
        bool aln1_is_rev = get<1>(chr.second);
        bool aln2_is_rev = get<2>(chr.second);
        bool same_orientation = aln1_is_rev == aln2_is_rev;
        fragment_orientations.push_front(same_orientation);
        if (fragment_orientations.size() > fragment_length_cache_size) {
            fragment_orientations.pop_back();
        }
        bool same_direction = true;
        if (aln1_is_rev && length <= 0) {
            same_direction = true;
        } else if (!aln1_is_rev && length >= 0) {
            same_direction = true;
        } else if (aln1_is_rev && length >= 0) {
            same_direction = false;
        } else if (!aln1_is_rev && length <= 0) {
            same_direction = false;
        } else {
            assert(false);
        }
        fragment_directions.push_front(same_direction);
        if (fragment_directions.size() > fragment_length_cache_size) {
            fragment_directions.pop_back();
        }
        // assume we can record the fragment length
        fragment_lengths.push_front(abs(length));
        if (fragment_lengths.size() > fragment_length_cache_size) {
            auto last = fragment_lengths.back();
            fragment_lengths.pop_back();
        }
        if (++since_last_fragment_length_estimate > fragment_model_update_interval) {
            cached_fragment_length_mean = fragment_length_mean();
            cached_fragment_length_stdev = fragment_length_stdev();
            cached_fragment_orientation_same = fragment_orientation();
            cached_fragment_direction = fragment_direction();
            // set our fragment size cap to the cached mean + 10x the standard deviation
            fragment_size = cached_fragment_length_mean + fragment_sigma * cached_fragment_length_stdev;
            since_last_fragment_length_estimate = 1;
        }
    }
}

double FragmentLengthStatistics::fragment_length_stdev(void) {
    return stdev(fragment_lengths);
}

double FragmentLengthStatistics::fragment_length_mean(void) {
    double sum = std::accumulate(fragment_lengths.begin(), fragment_lengths.end(), 0.0);
    return sum / fragment_lengths.size();
}

double FragmentLengthStatistics::fragment_length_pdf(double length) {
    return normal_pdf(length, cached_fragment_length_mean, cached_fragment_length_stdev);
}

// that the value is at least as extreme as this one
double FragmentLengthStatistics::fragment_length_pval(double length) {
    double x = abs(length-cached_fragment_length_mean)/cached_fragment_length_stdev;
    return 1 - (Phi(x) - Phi(-x));
}

bool FragmentLengthStatistics::fragment_orientation(void) {
    int count_same = 0;
    int count_diff = 0;
    for (auto& same_strand : fragment_orientations) {
        if (same_strand) ++count_same;
        else ++count_diff;
    }
    return count_same > count_diff;
}

bool FragmentLengthStatistics::fragment_direction(void) {
    int count_fwd = 0;
    int count_rev = 0;
    for (auto& go_forward : fragment_directions) {
        if (go_forward) ++count_fwd;
        else ++count_rev;
    }
    return count_fwd > count_rev;
}

// We need a function to get the lengths of nodes, in case we need to
// reverse an Alignment, including all its Mappings and Positions.
int64_t Mapper::get_node_length(int64_t node_id) {
    // Grab the node sequence only from the PathPositionHandleGraph index and get its size.
    // Make sure to use the cache
    return xindex->get_length(xindex->get_handle(node_id));
}

bool Mapper::check_alignment(const Alignment& aln) {
    // use the graph to extract the sequence
    // assert that this == the alignment
    if (aln.path().mapping_size()) {
        // get the graph corresponding to the alignment path
        bdsg::HashGraph sub;
        for (int i = 0; i < aln.path().mapping_size(); ++ i) {
            auto& m = aln.path().mapping(i);
            if (m.has_position() && m.position().node_id()) {
                auto id = aln.path().mapping(i).position().node_id();
                algorithms::extract_id_range(*xindex, id, id, sub);
            }
        }
        auto seq = algorithms::path_string(sub, aln.path());
        //if (aln.sequence().find('N') == string::npos && seq != aln.sequence()) {
        if (aln.quality().size() && aln.quality().size() != aln.sequence().size()) {
            cerr << "alignment quality is not the same length as its sequence" << endl
                 << pb2json(aln) << endl;
            return false;
        }
        if (seq != aln.sequence()) {
            cerr << "alignment does not match graph " << endl
                 << pb2json(aln) << endl
                 << "expect:\t" << aln.sequence() << endl
                 << "got:\t" << seq << endl;
            // save alignment
            vg::io::write_to_file(aln, "fail-" + hash_alignment(aln) + ".gam");
            // save graph, bigger fragment
            algorithms::expand_subgraph_by_steps(*xindex, sub, 5);
            algorithms::add_subpaths_to_subgraph(*xindex, sub);
            std::ofstream out("fail-" + hash_alignment(aln) + ".hg");
            sub.serialize(out);
            return false;
        }
    }
    return true;
}

vector<Alignment> Mapper::make_bands(const Alignment& read, int band_width, int band_overlap, vector<pair<int, int>>& to_strip) {
    int segment_size = min((int)read.sequence().size(), band_width);
#ifdef debug_mapper
    if (debug) {
        cerr << "Segment size be " << segment_size << "/" << read.sequence().size() << endl;
    }
#endif
    vector<int> start_positions;

    // record the start positions
    int offset = 0;
    while (offset+segment_size < read.sequence().size()) {
        if (offset == 0) {
            start_positions.push_back(0);
            offset += band_width-band_overlap;
        } else {
            start_positions.push_back(offset);
            offset += band_width-band_overlap;
        }
    }
    // add in the last alignment if we need it
    if (offset < read.sequence().size()) {
        start_positions.push_back(read.sequence().size()-segment_size);
    }
    int m = band_overlap % 2;
    // set up the structures to hold onto the banded alignments
    int to_align = start_positions.size();
    to_strip.resize(to_align);
    vector<Alignment> bands; bands.resize(to_align);
    int i = 0;
    for (auto& p : start_positions) {
        if (&p == &start_positions.front()) {
            if (start_positions.size() > 1) {
                to_strip[i].second = band_overlap/2 + m;
            }
        } else if (&p == &start_positions.back()) {
            to_strip[i].first = (start_positions[i-1]+segment_size-band_overlap/2-m)-p;
        } else {
            to_strip[i].first = band_overlap/2;
            to_strip[i].second = band_overlap/2 + m;
        }
        //cerr << "position: " << p << " strip " << to_strip[i].first << " " << to_strip[i].second << endl;
        auto& aln = bands[i];
        aln.set_name(read.name());
        aln.set_sequence(read.sequence().substr(p, segment_size));
        if (!read.quality().empty()) aln.set_quality(read.quality().substr(p, segment_size));
        ++i;
    }
    return bands;
}

vector<Alignment> Mapper::align_banded(const Alignment& read, int kmer_size, int stride, int max_mem_length, int band_width, int band_overlap, bool xdrop_alignment) {
    // cerr << read.sequence() << endl;
    auto aligner = get_aligner(!read.quality().empty());
    int8_t match = aligner->match;
    int8_t gap_extension = aligner->gap_extension;
    int8_t gap_open = aligner->gap_open;

    //cerr << "top of align_banded " << pb2json(read) << endl;
    // split the alignment up into overlapping chunks of band_width size
    // force used bandwidth to be divisible by 4
    // round up so we have > band_width

#ifdef debug_mapper
#pragma omp critical
    {
        if (debug) {
            cerr << "trying band width " << band_width << endl;
        }
    }
#endif

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    // scan across the read choosing bands
    // these bands are hard coded to overlap by band_overlap
    // the second and next-to-last bands begin or end remander/2 into the read
    // the last band is guaranteed to be segment_size long
    // overlap scheme example
    // read: ----------------------
    //       --------
    //         --------
    //              --------
    //                   --------
    //                     --------
    // Afterwards, we align each fragment, trim the overlaps implied by the layout
    // and build up an AlignmentChain model to generate maximum likelihood chains through
    // the bands which define our best to Nth best alignment. Changes to the chain model cost
    // function can be used to enable direct detection of SVs and other large scale variations.
    vector<pair<int, int>> to_strip;
    vector<Alignment> bands = make_bands(read, band_width, band_overlap, to_strip);
    vector<vector<Alignment>> multi_alns;
    multi_alns.resize(bands.size());

    auto do_band = [&](int i) {
        // cerr << "aligning band " << i << endl;
        vector<Alignment>& malns = multi_alns[i];
        double cluster_mq = 0;
        malns = align_multi_internal(true, bands[i], kmer_size, stride, max_mem_length, bands[i].sequence().size(), 0, cluster_mq, band_multimaps, extra_multimaps, nullptr, xdrop_alignment);
        for (vector<Alignment>::iterator a = malns.begin(); a != malns.end(); ++a) {
            Alignment& aln = *a;
            int mapqual = aln.mapping_quality();
            //aln = simplify(aln);
            bool above_threshold = false;
            if (aln.score() > 0) {
                // strip overlaps and re-score the part of the alignment we keep
                aln = strip_from_start(aln, to_strip[i].first);
                aln = strip_from_end(aln, to_strip[i].second);
                aln.set_identity(identity(aln.path()));
                above_threshold = aln.identity() >= min_identity && mapqual >= min_banded_mq;
            }
            if (!above_threshold) {
                // treat as unmapped
                aln = bands[i];
                aln = strip_from_start(aln, to_strip[i].first);
                aln = strip_from_end(aln, to_strip[i].second);
            }
        }
    };

    // Always use OMP parallelism here; restrict threads by setting OMP thread count.
#pragma omp parallel for
    for (int i = 0; i < bands.size(); ++i) {
        do_band(i);
    }

    // cost function
    auto transition_weight = [&](const Alignment& aln1, const Alignment& aln2,
                                 const unordered_map<path_handle_t, vector<pair<size_t, bool> > >& pos1,
                                 const unordered_map<path_handle_t, vector<pair<size_t, bool> > >& pos2,
                                 int64_t band_distance) {
        if (aln1.has_path() && !aln2.has_path()) {
            // pay a lot to go into unaligned from aligned because then we risk dropping into a random place
            return -(double)(5*aln2.sequence().size() * gap_extension + gap_open);
        } else if (!aln1.has_path() && !aln2.has_path()) {
            // pay to continue unaligned
            return -(double)(2*aln2.sequence().size() * gap_extension);
        } else if (!aln1.has_path() && aln2.has_path()) {
            return 0.0;
        }
        auto aln1_end = make_pos_t(path_end_position(aln1.path()));
        auto aln2_begin = make_pos_t(path_start_position(aln2.path()));
        pair<int64_t, int64_t> distances = min_oriented_distances(pos1, pos2);
        // consider both the forward and inversion case
        // counter[3]++; bench_t b; bench_init(b); bench_start(b);
        int64_t dist_fwd = distances.first;
        dist_fwd -= band_distance;
        int64_t dist_inv = distances.second;
        dist_inv -= band_distance;
        // bench_end(b);

        double fwd_score = -((double)(gap_open * dist_fwd != 0) + (double)dist_fwd * (double)gap_extension);
        double inv_score = -2.0*((double)(gap_open * dist_fwd != 0) + (double)dist_inv * (double)gap_extension);
        return max(fwd_score, inv_score);
    };

    // counter[0] += multi_alns.size();
    // bench_start(bench[0]);
    AlignmentChainModel chainer(multi_alns, this, transition_weight, max_band_jump, 64, max_band_jump*2);
    // bench_end(bench[0]);
    if (debug) chainer.display(cerr);
    int total_multimaps = max(max_multimaps, extra_multimaps);
    vector<Alignment> alignments = chainer.traceback(read, total_multimaps, false, debug);
    if (patch_alignments) {
        for (auto& aln : alignments) {
            // patch the alignment to deal with short unaligned regions
            aln = patch_alignment(aln, band_width/2, true, false); // can't use xdrop_alignment here as it needs seeds to work correctly
        }
    }
    // sort the alignments by score
    std::sort(alignments.begin(), alignments.end(), [](const Alignment& aln1, const Alignment& aln2) { return aln1.score() > aln2.score(); });
    if (alignments.empty()) {
        alignments.push_back(read);
        auto& aln = alignments.back();
        aln.clear_score();
        aln.clear_path();
        aln.clear_identity();
        aln.clear_mapping_quality();
    } else if (alignments.size() == 1) {
        alignments.front().set_mapping_quality(max_mapping_quality);
    } else {
        compute_mapping_qualities(alignments, 0, max_mapping_quality, max_mapping_quality);
        filter_and_process_multimaps(alignments, max_multimaps);
    }
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    // set time for alignment
    alignments.front().set_time_used(chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count());
    return alignments;
}

bool Mapper::adjacent_positions(const Position& pos1, const Position& pos2) {
    // are they the same id, with offset differing by 1?
    if (pos1.node_id() == pos2.node_id()
        && pos1.is_reverse() == pos2.is_reverse()
        && pos1.offset() == pos2.offset()-1) {
        return true;
    }
    // otherwise, we're going to need to check via the index
    // look in the graph to figure out if we are adjacent
    handle_t handle1 = xindex->get_handle(pos1.node_id(), pos1.is_reverse());
    handle_t handle2 = xindex->get_handle(pos2.node_id(), pos2.is_reverse());
    bool adjacent = false;
    xindex->follow_edges(handle1, false, [&](const handle_t& next) {
            if (next == handle2) {
                adjacent = true;
                return false;
            }
            return true;
        });
    return adjacent;
}

void Mapper::compute_mapping_qualities(vector<Alignment>& alns, double cluster_mq, double mq_estimate, double mq_cap) {
    if (alns.empty()) return;
    double max_mq = min(mq_cap, (double)max_mapping_quality);
    const GSSWAligner* aligner = get_aligner();
    int sub_overlaps = sub_overlaps_of_first_aln(alns, mq_overlap);
    switch (mapping_quality_method) {
        case Approx:
            aligner->compute_mapping_quality(alns, max_mq, true, cluster_mq, use_cluster_mq, sub_overlaps, mq_estimate, maybe_mq_threshold, identity_weight);
            break;
        case Exact:
            aligner->compute_mapping_quality(alns, max_mq, false, cluster_mq, use_cluster_mq, sub_overlaps, mq_estimate, maybe_mq_threshold, identity_weight);
            break;
        default: // None
            break;
    }
}
    
void Mapper::compute_mapping_qualities(pair<vector<Alignment>, vector<Alignment>>& pair_alns, double cluster_mq, double mq_estimate1, double mq_estimate2, double mq_cap1, double mq_cap2) {
    if (pair_alns.first.empty() || pair_alns.second.empty()) return;
    double max_mq1 = min(mq_cap1, (double)max_mapping_quality);
    double max_mq2 = min(mq_cap2, (double)max_mapping_quality);
    const GSSWAligner* aligner = get_aligner();
    int sub_overlaps1 = sub_overlaps_of_first_aln(pair_alns.first, mq_overlap);
    int sub_overlaps2 = sub_overlaps_of_first_aln(pair_alns.second, mq_overlap);
    vector<double> frag_weights;
    for (int i = 0; i < pair_alns.first.size(); ++i) {
        auto& aln1 = pair_alns.first[i];
        frag_weights.push_back(aln1.fragment_score());
    }
    switch (mapping_quality_method) {
        case Approx:
            aligner->compute_paired_mapping_quality(pair_alns, frag_weights, max_mq1, max_mq2, true, cluster_mq, use_cluster_mq, sub_overlaps1, sub_overlaps2, mq_estimate1, mq_estimate2, maybe_mq_threshold, identity_weight);
            break;
        case Exact:
            aligner->compute_paired_mapping_quality(pair_alns, frag_weights, max_mq1, max_mq2, false, cluster_mq, use_cluster_mq, sub_overlaps1, sub_overlaps2, mq_estimate1, mq_estimate2, maybe_mq_threshold, identity_weight);
            break;
        default: // None
            break;
    }
}

double Mapper::estimate_max_possible_mapping_quality(int length, double min_diffs, double next_min_diffs) {
    return get_aligner()->estimate_max_possible_mapping_quality(length, min_diffs, next_min_diffs);
}

double Mapper::max_possible_mapping_quality(int length) {
    return get_aligner()->max_possible_mapping_quality(length);
}

vector<Alignment> Mapper::score_sort_and_deduplicate_alignments(vector<Alignment>& all_alns, const Alignment& original_alignment) {
    if (all_alns.size() == 0) {
        all_alns.emplace_back();
        Alignment& aln = all_alns.back();
        aln = original_alignment;
        aln.clear_path();
        aln.set_score(0);
        return all_alns;
    }

    // Hold all the alignments bucketed by score.
    // The buckets are vectors and not sets to avoild introducing a dependency on pointer value or unordered set order.
    map<int, vector<Alignment*> > alignment_by_score;
    for (auto& ta : all_alns) {
        Alignment* aln = &ta;
        alignment_by_score[aln->score()].push_back(aln);
    }
    // TODO: Filter down subject to a minimum score per base or something?
    // Collect all the unique alignments (to compute mapping quality) and order by score
    vector<Alignment> sorted_unique_alignments;
    for(auto it = alignment_by_score.rbegin(); it != alignment_by_score.rend(); ++it) {
        // Copy over all the alignments in descending score order (following the pointers into the "alignments" vector)
        // Iterating through a set keyed on ints backward is in descending order.
        
        // This is going to let us deduplicate our alignments with this score, by storing them serialized to strings in this set.
        set<string> serializedAlignmentsUsed;
        
        for(Alignment* pointer : (*it).second) {
            // We serialize the alignment to a string
            string serialized;
            pointer->SerializeToString(&serialized);
            
            if(!serializedAlignmentsUsed.count(serialized)) {
                // This alignment hasn't been produced yet. Produce it. The
                // order in the alignment vector doesn't matter for things with
                // the same score.
                sorted_unique_alignments.push_back(*pointer);
                
                // Save it so we can avoid putting it in the vector again
                serializedAlignmentsUsed.insert(serialized);
            }
        }
        
        if (it == alignment_by_score.rbegin()) {
            // We just did the top-scoring group of alignments.
            // Shuffle them to avoid bias.
            deterministic_shuffle(sorted_unique_alignments.begin(), sorted_unique_alignments.end());
        }
    
    }
    return sorted_unique_alignments;
}

// filters down to requested number of alignments and marks
void Mapper::filter_and_process_multimaps(vector<Alignment>& sorted_unique_alignments, int total_multimaps) {
    if (sorted_unique_alignments.size() > total_multimaps){
        sorted_unique_alignments.resize(total_multimaps);
    }
    
    // TODO log best alignment score?
    for(size_t i = 0; i < sorted_unique_alignments.size(); i++) {
        // Mark all but the first, best alignment as secondary
        sorted_unique_alignments[i].set_is_secondary(i > 0);
    }
}
    
vector<Alignment> Mapper::align_multi(const Alignment& aln, int kmer_size, int stride, int max_mem_length, int band_width, int band_overlap, bool xdrop_alignment) {
    double cluster_mq = 0;
    Alignment clean_aln;
    clean_aln.set_name(aln.name());
    clean_aln.set_sequence(aln.sequence());
    clean_aln.set_quality(aln.quality());
    clean_aln.clear_refpos();
    return align_multi_internal(true, clean_aln, kmer_size, stride, max_mem_length, band_width, band_overlap, cluster_mq, max_multimaps, extra_multimaps, nullptr, xdrop_alignment);
}
    
vector<Alignment> Mapper::align_multi_internal(bool compute_unpaired_quality,
                                               const Alignment& aln,
                                               int kmer_size, int stride,
                                               int max_mem_length,
                                               int band_width,
                                               int band_overlap,
                                               double& cluster_mq,
                                               int keep_multimaps,
                                               int additional_multimaps,
                                               vector<MaximalExactMatch>* restricted_mems,
                                               bool xdrop_alignment) {
    
    if(debug) {
#pragma omp critical
        cerr << "align_multi_internal("
            << compute_unpaired_quality << ", " 
            << aln.sequence() << ", " 
            << kmer_size << ", " 
            << stride << ", " 
            << band_width << ", "
            << keep_multimaps << ", "
            << additional_multimaps << ", "
            << restricted_mems << ")" 
            << endl;
        if (aln.has_path()) {
            // if we're realigning, show in the debugging output what we start with
            cerr << pb2json(aln) << endl;
        }
    }
    // make sure to respect the max multimaps if we haven't been givne a keep multimap count
    if (keep_multimaps == 0) keep_multimaps = max_multimaps;

    // trigger a banded alignment if we need to
    // note that this will in turn call align_multi_internal on fragments of the read
    if (aln.sequence().size() > band_width) {
        // TODO: banded alignment currently doesn't support mapping qualities because it only produces one alignment
#ifdef debug_mapper
#pragma omp critical
        if (debug) cerr << "switching to banded alignment" << endl;
#endif
        return vector<Alignment>{align_banded(aln, kmer_size, stride, max_mem_length, band_width, band_overlap, xdrop_alignment)};
    }

    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    // try to get at least 2 multimaps so that we can calculate mapping quality
    int additional_multimaps_for_quality;
    if (additional_multimaps == 0 && max_multimaps == 1 && mapping_quality_method != None) {
        additional_multimaps_for_quality = 1;
    }
    else {
        additional_multimaps_for_quality = additional_multimaps;
    }

    double longest_lcp, fraction_filtered;
    vector<Alignment> alignments;
    
    // use pre-restricted mems for paired mapping or find mems here
    if (restricted_mems != nullptr) {
        // mem hits will already have been queried
        alignments = align_mem_multi(aln, *restricted_mems, cluster_mq, longest_lcp, fraction_filtered, max_mem_length, keep_multimaps, additional_multimaps_for_quality, xdrop_alignment);
    }
    else {
        vector<MaximalExactMatch> mems = find_mems_deep(aln.sequence().begin(),
                                                        aln.sequence().end(),
                                                        longest_lcp,
                                                        fraction_filtered,
                                                        max_mem_length,
                                                        min_mem_length,
                                                        mem_reseed_length,
                                                        false, true, true, true); // Make sure to actually fill in the longest LCP.
        
        // query mem hits
        alignments = align_mem_multi(aln, mems, cluster_mq, longest_lcp, fraction_filtered, max_mem_length, keep_multimaps, additional_multimaps_for_quality, xdrop_alignment);
    }

#ifdef debug_mapper
    for (auto& aln : alignments) {
        // Make sure no alignments are wandering out of the graph
        for (size_t i = 0; i < aln.path().mapping_size(); i++) {
            // Look at each mapping
            auto& mapping = aln.path().mapping(i);
            
            if (mapping.position().node_id()) {
                // Get the size of its node from whatever index we have
                size_t node_size = get_node_length(mapping.position().node_id());
                
                // Make sure the mapping fits in the node
                assert(mapping.position().offset() + mapping_from_length(mapping) <= node_size);
            }
        }
    }
#endif
    
    algorithms::annotate_with_initial_path_positions(*xindex, alignments);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    // set time for alignment
    alignments.front().set_time_used(chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count());
    return alignments;
}

Alignment Mapper::align(const Alignment& aln, int kmer_size, int stride, int max_mem_length, int band_width, int band_overlap, bool xdrop_alignment) {
    // TODO computing mapping quality could be inefficient depending on the method chosen
    
    // Do the multi-mapping
    vector<Alignment> best = align_multi(aln, kmer_size, stride, max_mem_length, band_width, band_overlap, xdrop_alignment);

    if(best.size() == 0) {
        // Spit back an alignment that says we failed, but make sure it has the right sequence in it.
        Alignment failed = aln;
        failed.clear_path();
        failed.set_score(0);
        return failed;
    }

    // Otherwise, just report the best alignment, since we know one exists
    return best[0];
}

set<pos_t> gcsa_nodes_to_positions(const vector<gcsa::node_type>& nodes) {
    set<pos_t> positions;
    for(gcsa::node_type node : nodes) {
        positions.insert(make_pos_t(node));
    }
    return positions;    
}

int64_t Mapper::graph_mixed_distance_estimate(pos_t pos1, pos_t pos2, int64_t maximum) {
    int64_t path_dist = algorithms::min_approx_path_distance(xindex, pos1, pos2, maximum*2);
    int64_t approx_dist = abs(approx_distance(pos1, pos2));
    return min(path_dist, approx_dist);
}

int64_t Mapper::approx_position(pos_t pos) {
    // get nodes on the forward strand
    if (is_rev(pos)) {
        pos = reverse(pos, get_node_length(id(pos)));
    }
    return dynamic_cast<VectorizableHandleGraph*>(xindex)->node_vector_offset((nid_t)id(pos)) + offset(pos);
}

int64_t Mapper::approx_distance(pos_t pos1, pos_t pos2) {
    return approx_position(pos1) - approx_position(pos2);
}

/// returns approximate position of alignnment start in xindex
/// or -1.0 if alignment is unmapped
int64_t Mapper::approx_alignment_position(const Alignment& aln) {
    if (aln.path().mapping_size()) {
        for (int i = 0; i < aln.path().mapping_size(); ++i) {
            auto& mbeg = aln.path().mapping(i);
            if (mbeg.has_position()) {
                return approx_position(make_pos_t(mbeg.position()));
            }
        }
    }
    return -1.0;
}

/// returns approximate distance between alignment starts
/// or -1.0 if not possible to determine
int64_t Mapper::approx_fragment_length(const Alignment& aln1, const Alignment& aln2) {
    int64_t pos1 = approx_alignment_position(aln1);
    int64_t pos2 = approx_alignment_position(aln2);
    if (pos1 != -1 && pos2 != -1) {
        return abs(pos1 - pos2);
    } else {
        return -1;
    }
}

Position Mapper::alignment_end_position(const Alignment& aln) {
    if (!aln.has_path()) { Position pos; return pos; }
    Alignment b;
    *b.mutable_path()->add_mapping() = aln.path().mapping(aln.path().mapping_size()-1);
    b = reverse_complement_alignment(b,
                                     (function<int64_t(int64_t)>) ([&](int64_t id) {
                                             return (int64_t)get_node_length(id);
                                         }));
    return reverse(b.path().mapping(0).position(),
                   get_node_length(b.path().mapping(0).position().node_id()));
}

Alignment Mapper::patch_alignment(const Alignment& aln, int max_patch_length, bool trim_internal_deletions, bool xdrop_alignment) {
    //cerr << "top of patch_alignment" << endl;
    Alignment patched;
    // walk along the alignment and find the portions that are unaligned
    int read_pos = 0;
    double extend_fwd = 1.61803;
    double extend_rev = 1.0/1.61803;
    auto& path = aln.path();
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
        //cerr << "looking at mapping " << i << " " << pb2json(mapping) << endl;
        pos_t ref_pos = make_pos_t(mapping.position());
        for (int j = 0; j < mapping.edit_size(); ++j) {
            auto& edit = mapping.edit(j);
            //cerr << "looking at edit " << j << " " << pb2json(edit) << endl;
            if (edit_is_match(edit) || edit_is_sub(edit)) {
                patched.mutable_sequence()->append(aln.sequence().substr(read_pos, edit.to_length()));
                Mapping* new_mapping = patched.mutable_path()->add_mapping();
                *new_mapping->mutable_position() = make_position(ref_pos);
                *new_mapping->add_edit() = edit;
                get_offset(ref_pos) += edit.from_length();
            } else if (edit_is_deletion(edit)) {
                Mapping* new_mapping = patched.mutable_path()->add_mapping();
                *new_mapping->mutable_position() = make_position(ref_pos);
                *new_mapping->add_edit() = edit;
                get_offset(ref_pos) += edit.from_length();
            } else if (edit_is_insertion(edit)) {
                // let's try to patch this insertion into the gap between the neighboring bits
                //cerr << "patching " << pb2json(edit) << endl;
                Alignment patch;
                patch.set_sequence(edit.sequence());
                if (!aln.quality().empty()) {
                    patch.set_quality(aln.quality().substr(read_pos, edit.sequence().size()));
                }
                vector<Alignment> bands;
                vector<pair<int, int>> to_strip;
                if (edit.sequence().size() > max_patch_length) {
                    // do the banding thing
                    //cerr << "banding" << endl;
                    bands = make_bands(patch, max_patch_length, max_patch_length/4, to_strip);
                } else {
                    //cerr << "not banding" << endl;
                    to_strip.push_back(make_pair(0,0));
                    bands.push_back(patch);
                }
                //cerr << "got " << bands.size() << " bands" << endl;
                //pos_t band_ref_pos = ref_pos;
                set<pos_t> band_ref_pos;
                if (id(ref_pos)) band_ref_pos.insert(ref_pos);
                //cerr << "band ref pos size " << band_ref_pos.size() << endl;
                for (int k = 0; k < bands.size(); ++k) {
                    Alignment& band = bands[k];
                    //cerr << "on band " << k << " " << pb2json(band) << endl;
                    if (i == 0 && j == 0 && k == 0) {
                        //cerr << "soft clip at start " << band.sequence() << endl;
                        // reverse the position, we're going backwards to get the graph off the end of where we are
                        int max_score = -std::numeric_limits<int>::max();
                        for (auto& pos : band_ref_pos) {
                            //cerr << "trying position " << pos << endl;
                            bdsg::HashGraph graph;
                            algorithms::extract_context(*xindex, graph, xindex->get_handle(id(pos), is_rev(pos)), offset(pos),
                                                        band.sequence().size()*extend_fwd, true, false);
                            algorithms::extract_context(*xindex, graph, xindex->get_handle(id(pos), is_rev(pos)), offset(pos),
                                                        band.sequence().size()*extend_rev, false, true);
                            auto proposed_band = align_maybe_flip(band, graph, is_rev(pos), true, false, xdrop_alignment);
                            if (proposed_band.score() > max_score) { band = proposed_band; max_score = band.score(); }
                        }
                        // TODO
                        // ideal: edit it with the path fragment after this, to get the breakpoints right
                        // and trim off the part matching that bit
                    } else if (i == path.mapping_size()-1
                               && j == mapping.edit_size()-1
                               && k == bands.size()-1) {
                        //cerr << "soft clip at end " << band.sequence() << endl;
                        int max_score = -std::numeric_limits<int>::max();
                        for (auto& pos : band_ref_pos) {
                            //cerr << "trying position " << pos << endl;
                            bdsg::HashGraph graph;
                            algorithms::extract_context(*xindex, graph, xindex->get_handle(id(pos), is_rev(pos)), offset(pos),
                                                        band.sequence().size()*extend_fwd, true, false);
                            algorithms::extract_context(*xindex, graph, xindex->get_handle(id(pos), is_rev(pos)), offset(pos),
                                                        band.sequence().size()*extend_rev, false, true);
                            //cerr << "on graph " << pb2json(graph) << endl;
                            auto proposed_band = align_maybe_flip(band, graph, is_rev(pos), true, false, xdrop_alignment);
                            if (proposed_band.score() > max_score) { band = proposed_band; max_score = band.score(); }
                        }
                    } else {
                        //cerr << "internal addition" << endl;
                        int max_score = -std::numeric_limits<int>::max();
                        for (auto& pos : band_ref_pos) {
                            //cerr << "trying position " << pos << endl;
                            bdsg::HashGraph graph;
                            algorithms::extract_context(*xindex, graph, xindex->get_handle(id(pos), is_rev(pos)), offset(pos),
                                                        band.sequence().size()*extend_fwd, true, false);
                            algorithms::extract_context(*xindex, graph, xindex->get_handle(id(pos), is_rev(pos)), offset(pos),
                                                        band.sequence().size()*extend_rev, false, true);
                            //cerr << "on graph " << pb2json(graph) << endl;
                            auto proposed_band = align_maybe_flip(band, graph, is_rev(pos), true, false, xdrop_alignment);
                            if (proposed_band.score() > max_score) { band = proposed_band; max_score = band.score(); }
                        }
                    }
#ifdef debug_mapper
                    if (debug && !check_alignment(band)) {
                        cerr << "patching failure " << pb2json(band) << endl;
                        assert(false);
                    }
#endif
                    assert(band.sequence().size() > to_strip[k].first + to_strip[k].second);
                    if (band.path().mapping_size() == 0) { band.clear_path(); } // failed alignment
                    band = strip_from_start(band, to_strip[k].first);
                    band = strip_from_end(band, to_strip[k].second);
                    band = simplify(band);
                    band.set_identity(identity(band.path()));
                    // update the reference end position
                    if (band.has_path()) {
                        auto from_length = alignment_from_length(band);
                        if (from_length >= min_mem_length
                            && from_length >= min_cluster_length
                            && band.identity() > min_identity) {
                            band_ref_pos.clear();
                            // todo... step our position back just a little to match the banding
                            // right now we're relying on the chunkiness of the graph to get this for us
                            // strip back a little
                            if (to_strip.size() > k+1 && to_strip[k+1].first && band.sequence().size() > to_strip[k+1].first) {
                                auto band_chew = strip_from_end(band, to_strip[k+1].first);
                                pos_t pos = make_pos_t(alignment_end_position(band_chew));
                                if (id(pos)) band_ref_pos.insert(pos);
                            } else {
                                pos_t pos = make_pos_t(alignment_end_position(band));
                                if (id(pos)) band_ref_pos.insert(pos);
                            }
#ifdef debug_mapper
                            if (debug && !check_alignment(band)) {
                                cerr << "failed band " << pb2json(band) << endl;
                                assert(false);
                            }
#endif
                        } else {
                            //cerr << "clearing the path" << endl;
                            band.clear_path();
                            /*
                            Edit* e = band.mutable_path()->add_mapping()->add_edit();
                            e->set_to_length(band.sequence().size());
                            band.clear_score();
                            band.clear_identity();
                            */
                            // TODO try to align over a bigger chunk after this
                        }
                    }
                    if (band.path().mapping_size() == 0) {
                        //cerr << "the hack" << endl;
                        // hack to make merging work when we are unaligned
                        band.set_score(0);
                        band.clear_path();
                        Mapping* m = band.mutable_path()->add_mapping();
                        Edit* e = m->add_edit();
                        e->set_sequence(band.sequence());
                        e->set_to_length(band.sequence().size());
                        if (!band_ref_pos.empty()) {
                            *m->mutable_position() = make_position(*band_ref_pos.begin());
                        } // should we alternatively set ourselves back to the ref_pos?
                        // walk graph to estimate next position based on assumption we are ~ homologous to the graph
                        set<pos_t> next_pos;
                        for (auto& pos : band_ref_pos) {
                            for (auto& next : algorithms::jump_along_closest_path(xindex, pos, band.sequence().size(), band.sequence().size())) {
                                next_pos.insert(next);
                            }
                        }
                        // keep only 4 next positions
                        band_ref_pos.clear();
                        if (next_pos.size()) {
                            int to_keep = 4;
                            int p = 0;
                            for (auto& pos : next_pos) {
                                if (++p <= to_keep) {
                                    band_ref_pos.insert(pos);
                                }
                            }
                        }
                    }
                }

                /*
                if (debug) {
                    cerr << "done bands" << endl;
                    for (auto& band : bands) {
                        cerr << "band: " << pb2json(band) << endl;
                    }
                }
                */
                patch = simplify(merge_alignments(bands));
                if (patch.sequence() != edit.sequence()) {
                    cerr << "sequence mismatch" << endl;
                    cerr << "seq_expect: " << edit.sequence() << endl;
                    cerr << "seq_got:    " << patch.sequence() << endl;
                    assert(false);
                }
#ifdef debug_mapper
                if (debug && !check_alignment(patch)) {
                    cerr << "failed patch merge " << pb2json(patch) << endl;
                    assert(false);
                }
#endif
                if (patched.path().mapping_size()) {
                    extend_alignment(patched, patch, true);
                } else {
                    patched = patch;
                }
                /*
                if (!check_alignment(patched)) {
                    cerr << "failed patched extend " << pb2json(patched) << endl;
                    assert(false);
                }
                */
            }
            read_pos += edit.to_length();
        }
    }
    // finally, fix up the alignment score
    patched.set_name(aln.name());
    patched.set_sequence(aln.sequence());
    if (!aln.quality().empty()) {
        patched.set_quality(aln.quality());
    }
#ifdef debug_mapper
    if (debug && !check_alignment(patched)) {
        cerr << "failed final alignment " << pb2json(patched) << endl;
        assert(false);
    }
#endif
    // simplify the mapping representation
    auto clear_positions = [&](Alignment& aln) {
        for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
            Mapping* mapping = aln.mutable_path()->mutable_mapping(i);
            if (mapping->has_position() && from_length(*mapping) == 0) {
                mapping->clear_position();
            }
        }
    };
    clear_positions(patched);
    patched = simplify(patched, trim_internal_deletions);
    // set the identity
    patched.set_identity(identity(patched.path()));
    // recompute the score
    patched.set_score(score_alignment(patched, false));
    return patched;
}

void Mapper::remove_full_length_bonuses(Alignment& aln) {
    int32_t score = aln.score();
    int8_t bonus = get_aligner(!aln.quality().empty())->full_length_bonus;
    if (softclip_start(aln) == 0) score -= bonus;
    if (softclip_end(aln) == 0) score -= bonus;
    aln.set_score(score);
}

// generate a score from the alignment without realigning
// handles split alignments, where gaps of unknown length are
// by estimating length using the positional paths embedded in the graph
int32_t Mapper::score_alignment(const Alignment& aln, bool use_approx_distance) {

    // Find the right aligner to score with
    const GSSWAligner* aligner = get_aligner();
    
    if (use_approx_distance) {
        // Use an approximation
        return aligner->score_discontiguous_alignment(aln, [&](pos_t last, pos_t next, size_t max_search) {
            return approx_distance(last, next);
        }, strip_bonuses);
    } else {
        // Use the exact method, and if we hit the limit, fall back to the approximate method.
        return aligner->score_discontiguous_alignment(aln, [&](pos_t last, pos_t next, size_t max_search) {
                return graph_mixed_distance_estimate(last, next, min(32, (int)max_search));
        }, strip_bonuses);
    }
    
}

const int balanced_stride(int read_length, int kmer_size, int stride) {
    double r = read_length;
    double k = kmer_size;
    double j = stride;
    int i = (r > j) ? round((r-k)/round((r-k)/j)) : j;
    return max(1, i);
}

const vector<string> balanced_kmers(const string& seq, const int kmer_size, const int stride) {
    // choose the closest stride that will generate balanced kmers
    vector<string> kmers;
    int b = balanced_stride(seq.size(), kmer_size, stride);
    if (!seq.empty()) {
        for (int i = 0; i+kmer_size <= seq.size(); i+=b) {
            kmers.push_back(seq.substr(i,kmer_size));
        }
    }
    return kmers;
}

/////// banded long read alignment resolution

AlignmentChainModel::AlignmentChainModel(
    vector<vector<Alignment> >& bands,
    Mapper* mapper,
    const function<double(const Alignment&, const Alignment&, const unordered_map<path_handle_t, vector<pair<size_t, bool> > >&, const unordered_map<path_handle_t, vector<pair<size_t, bool> > >&, int64_t)>& transition_weight,
    int vertex_band_width,
    int position_depth,
    int max_connections) {
    // store the Alignments in the model
    int offset = 0;
    int idx = 0;
    // bench_start(mapper->bench[1]);
    for (auto& band : bands) {
        for (auto& aln : band) {
            if (!aln.has_path() || aln.score() == 0) continue;
            AlignmentChainModelVertex v;
            v.aln = &aln;
            v.band_begin = offset;
            v.band_idx = idx;
            v.weight = aln.score();
            v.prev = nullptr;
            v.score = 0;
            v.positions = algorithms::alignment_path_offsets(*mapper->xindex, aln, true, false);
            v.positions[handlegraph::as_path_handle(0)].push_back(make_pair(mapper->approx_alignment_position(aln), false));
            model.push_back(v);
            // mapper->counter[1]++;
        }
        assert(!band.empty());
        // save an unaligned band to fill in later
        Alignment aln = band.front();
        aln.clear_score(); aln.clear_path(); aln.clear_mapping_quality(); aln.clear_identity();
        unaligned_bands.push_back(aln);
        offset += band.front().sequence().length();
        ++idx;
    }
    for (vector<AlignmentChainModelVertex>::iterator v = model.begin(); v != model.end(); ++v) {
        for (auto u = v+1; u != model.end(); ++u) {
            if (v->band_idx != u->band_idx &&
                v->next_cost.size() < max_connections && u->prev_cost.size() < max_connections) {
                if (v->band_idx + vertex_band_width >= u->band_idx) {
                    // bench_start(mapper->bench[2]);
                    double weight = transition_weight(*v->aln, *u->aln, v->positions, u->positions,
                                                      u->band_begin - (v->band_begin+v->aln->sequence().size()));
                    // bench_end(mapper->bench[2]);
                    if (weight > -std::numeric_limits<double>::max()) {
                        v->next_cost.push_back(make_pair(&*u, weight));
                        u->prev_cost.push_back(make_pair(&*v, weight));
                    }
                    // mapper->counter[2]++;
                } else {
                    break;
                }
            }
        }
    }
    // bench_end(mapper->bench[1]);
}

void AlignmentChainModel::score(const unordered_set<AlignmentChainModelVertex*>& exclude) {
    // propagate the scores in the model
    for (auto& m : model) {
        // score is equal to the max inbound + mem.weight
        if (exclude.count(&m)) continue; // skip if vertex was whole cluster
        m.score = m.weight;
        for (auto& p : m.prev_cost) {
            if (p.first == nullptr) continue; // this transition is masked out
            double proposal = m.weight + p.second + p.first->score;
            if (proposal > m.score) {
                m.prev = p.first;
                m.score = proposal;
            }
        }
    }
}

AlignmentChainModelVertex* AlignmentChainModel::max_vertex(void) {
    AlignmentChainModelVertex* maxv = nullptr;
    for (auto& m : model) {
        if (maxv == nullptr || m.score > maxv->score) {
            maxv = &m;
        }
    }
    return maxv;
}

void AlignmentChainModel::clear_scores(void) {
    for (auto& m : model) {
        m.score = 0;
        m.prev = nullptr;
    }
}

vector<Alignment> AlignmentChainModel::traceback(const Alignment& read, int alt_alns, bool paired, bool debug) {
    debug = true;
    vector<vector<Alignment> > traces;
    traces.reserve(alt_alns); // avoid reallocs so we can refer to pointers to the traces
    unordered_set<AlignmentChainModelVertex*> exclude;
    for (auto& v : redundant_vertexes) exclude.insert(&*v);
    for (int i = 0; i < alt_alns; ++i) {
        // score the model, accounting for excluded traces
        clear_scores();
        score(exclude);
#ifdef debug_mapper
#pragma omp critical
        if (debug) {
            cerr << "AlignmentChainModel::traceback " << i << endl;
            display(cerr);
        }
#endif
        vector<AlignmentChainModelVertex*> vertex_trace;
        {
            // find the maximum score
            auto* vertex = max_vertex();
            // check if we've exhausted our Alignments
            if (vertex == nullptr || vertex->score == 0) break;
#ifdef debug_mapper
#pragma omp critical
            if (debug) cerr << "maximum score " << vertex->aln->sequence() << " " << vertex << ":" << vertex->score << endl;
#endif
            // make trace
            while (vertex != nullptr) {
                vertex_trace.push_back(vertex);
                if (vertex->prev != nullptr) {
                    vertex = vertex->prev;
                } else {
                    break;
                }
            }
        }
        // if we have a singular match or reads are not paired, record not to use it again
        if (paired && vertex_trace.size() == 1) {
            exclude.insert(vertex_trace.front());
        }
        // fill this out when we're paired to help mask out in-fragment transitions
        unordered_set<AlignmentChainModelVertex*> chain_members;
        if (paired) for (auto v : vertex_trace) chain_members.insert(v);
        traces.emplace_back();
        auto& aln_trace = traces.back();
        aln_trace = unaligned_bands;
        for (auto v = vertex_trace.rbegin(); v != vertex_trace.rend(); ++v) {
            auto& vertex = **v;
            if (!paired) exclude.insert(&vertex);
            if (v != vertex_trace.rbegin()) {
                auto y = v - 1;
                AlignmentChainModelVertex* prev = *y;
                // mask out used transitions
                for (auto& p : vertex.prev_cost) {
                    if (p.first == prev) {
                        p.first = nullptr;
                    } else if (paired && p.first != nullptr
                               && p.first->band_begin != vertex.band_begin
                               && chain_members.count(p.first)) {
                        p.first = nullptr;
                    }
                }
            }
            aln_trace[vertex.band_idx] = *vertex.aln;
        }
        //display_dot(cerr, vertex_trace);
    }
    vector<Alignment> alns;
    for (auto& trace : traces) {
        alns.emplace_back();
        Alignment& merged = alns.back();
        merged = simplify(merge_alignments(trace));
        merged.set_identity(identity(merged.path()));
        merged.set_quality(read.quality());
        merged.set_name(read.name());
    }
    return alns;
}

// show model
void AlignmentChainModel::display(ostream& out) {
    for (auto& vertex : model) {
        out << &vertex << ":" << vertex.band_begin << ":" << vertex.aln->sequence() << ":" << vertex.score << "@";
        out << "prev: ";
        for (auto& p : vertex.prev_cost) {
            auto& next = p.first;
            if (p.first == nullptr) continue;
            out << p.first << ":" << p.second << "@";
            out << " ; ";
        }
        out << " next: ";
        for (auto& p : vertex.next_cost) {
            auto& next = p.first;
            if (p.first == nullptr) continue;
            out << p.first << ":" << p.second << "@";
            out << " ; ";
        }
        out << endl;
    }
}

void AlignmentChainModel::display_dot(ostream& out, vector<AlignmentChainModelVertex*> vertex_trace) {
    map<AlignmentChainModelVertex*, int> vertex_ids;
    int i = 0;
    for (auto& vertex : model) {
        vertex_ids[&vertex] = ++i;
    }
    map<AlignmentChainModelVertex*,int> in_trace;
    i = 0;
    for (auto& v : vertex_trace) {
        in_trace[v] = ++i;
    }
    out << "digraph memchain {" << endl;
    out << "rankdir=LR;" << endl;
    for (auto& vertex : model) {
        out << vertex_ids[&vertex]
            << " [label=\""
            << vertex.aln->sequence()
            << "\nid:" << vertex_ids[&vertex]
            << " band:" << vertex.band_idx
            << " weight:" << vertex.weight
            << " score:" << vertex.score;
        out << "\npositions: ";
        for (auto& p : vertex.positions) {
            for (auto& q : p.second) {
                out << as_integer(p.first) << ":" << q.first << (q.second?"-":"+") << " ";
            }
        }
        pos_t p_start = make_pos_t(path_start_position(vertex.aln->path()));
        pos_t p_end = make_pos_t(path_end_position(vertex.aln->path()));
        out << "\nstart: " << id(p_start) << (is_rev(p_start)?"-":"+") << ":" << offset(p_start) << " ";
        out << "end: " << id(p_end) << (is_rev(p_end)?"-":"+") << ":" << offset(p_end) << " ";
        out << "\"";
        if (in_trace.find(&vertex) != in_trace.end()) {
            out << ",color=red";
        }
        out << ",shape=box];" << endl;
        for (auto& p : vertex.next_cost) {
            out << vertex_ids[&vertex] << " -> " << vertex_ids[p.first] << " [label=\"" << p.second << "\"";
            if (in_trace.find(&vertex) != in_trace.end() && in_trace.find(p.first) != in_trace.end() &&
                in_trace[&vertex] - 1 == in_trace[p.first]) {
                out << ",color=red,fontcolor=red";
            }
            out << "];" << endl;
        }
    }
    out << "}" << endl;
}

FragmentLengthDistribution::FragmentLengthDistribution(size_t maximum_sample_size,
                                                       size_t reestimation_frequency,
                                                       double robust_estimation_fraction) :
    maximum_sample_size(maximum_sample_size),
    reestimation_frequency(reestimation_frequency),
    robust_estimation_fraction(robust_estimation_fraction)
{
    assert(0.0 < robust_estimation_fraction && robust_estimation_fraction < 1.0);
}

FragmentLengthDistribution::FragmentLengthDistribution() : FragmentLengthDistribution(0, 1, 0.5)
{
    
}
    
FragmentLengthDistribution::~FragmentLengthDistribution() {
    
}

void FragmentLengthDistribution::force_parameters(double mean, double stddev) {
    mu = mean;
    sigma = stddev;
    is_fixed = true;
}

void FragmentLengthDistribution::register_fragment_length(int64_t length) {
    // allow this function to operate fully in parallel once the distribution is
    // fixed (and hence threadsafe)
    if (is_fixed) {
        return;
    }
#pragma omp critical
    {
        // in case the distribution became fixed while this thread was waiting
        // to execute the critical block
        if (!is_fixed) {
            lengths.insert((double) length);
            if (lengths.size() == maximum_sample_size) {
                // we've reached the maximum sample we wanted, so fix the estimation
                estimate_distribution();
                is_fixed = true;
            }
            else if (lengths.size() % reestimation_frequency == 0) {
                estimate_distribution();
            }
        }
    }
}
    
void FragmentLengthDistribution::estimate_distribution() {
    // remove the tails from the estimation
    size_t to_skip = (size_t) (lengths.size() * (1.0 - robust_estimation_fraction) * 0.5);
    auto begin = lengths.begin();
    auto end = lengths.end();
    for (size_t i = 0; i < to_skip; i++) {
        begin++;
        end--;
        
    }
    // compute cumulants
    double count = 0.0;
    double sum = 0.0;
    double sum_of_sqs = 0.0;
    for (auto iter = begin; iter != end; iter++) {
        count += 1.0;
        sum += *iter;
        sum_of_sqs += (*iter) * (*iter);
    }
    // use cumulants to compute moments
    mu = sum / count;
    double raw_var = sum_of_sqs / count - mu * mu;
    // apply method of moments estimation using the appropriate truncated normal distribution
    double a = Phi_inv(1.0 - 0.5 * (1.0 - robust_estimation_fraction));
    sigma = sqrt(raw_var / (1.0 - 2.0 * a * normal_pdf(a, 0.0, 1.0)));
}
    
double FragmentLengthDistribution::mean() const {
    return mu;
}

double FragmentLengthDistribution::std_dev() const {
    return sigma;
}

bool FragmentLengthDistribution::is_finalized() const {
    return is_fixed;
}
    
size_t FragmentLengthDistribution::max_sample_size() const {
    return maximum_sample_size;
}
    
size_t FragmentLengthDistribution::curr_sample_size() const {
    return lengths.size();
}
    
multiset<double>::const_iterator FragmentLengthDistribution::measurements_begin() const {
    return lengths.begin();
}

multiset<double>::const_iterator FragmentLengthDistribution::measurements_end() const {
    return lengths.end();
}
}

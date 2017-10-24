#include <unordered_set>
#include "mapper.hpp"
#include "algorithms/extract_containing_graph.hpp"

//#define debug_mapper

namespace vg {

BaseMapper::BaseMapper(xg::XG* xidex,
                       gcsa::GCSA* g,
                       gcsa::LCPArray* a) :
      xindex(xidex)
      , gcsa(g)
      , lcp(a)
      , min_mem_length(1)
      , mem_reseed_length(0)
      , fast_reseed(true)
      , fast_reseed_length_diff(0.75)
      , hit_max(0)
      , alignment_threads(1)
      , qual_adj_aligner(nullptr)
      , regular_aligner(nullptr)
      , adjust_alignments_for_base_quality(false)
      , mapping_quality_method(Approx)
      , max_mapping_quality(60)
      , strip_bonuses(false)
      , assume_acyclic(false)
{
    init_aligner(default_match, default_mismatch, default_gap_open,
                 default_gap_extension, default_full_length_bonus);
    
    // TODO: removing these consistency checks because we seem to have violated them pretty wontonly in
    // the code base already by changing the members directly when they were still public
    
//    // allow a default constructor with no indexes
//    if (xidex || g || a) {
//        if(xidex == nullptr) {
//            // We need an XG graph.
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

BaseMapper::~BaseMapper(void) {
    
    clear_aligners();
    
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
    
    string::const_iterator cursor = seq_end;
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
    MaximalExactMatch match(cursor, cursor, full_range);
    gcsa::range_type last_range = match.range;
    --cursor; // start off looking at the last character in the query
    while (cursor >= seq_begin) {
        // hold onto our previous range
        last_range = match.range;
        // execute one step of LF mapping
        match.range = gcsa->LF(match.range, gcsa->alpha.char2comp[*cursor]);
        if (gcsa::Range::empty(match.range)
            || max_mem_length && match.end-cursor > max_mem_length
            || match.end-cursor > gcsa->order()) {
            // break on N; which for DNA we assume is non-informative
            // this *will* match many places in assemblies; this isn't helpful
            if (*cursor == 'N' || last_range == full_range) {
                // we mismatched in a single character
                // there is no MEM here
                match.begin = cursor+1;
                match.range = last_range;
                mems.push_back(match);
                match.end = cursor;
                match.range = full_range;
                --cursor;
            } else {
                // we've exhausted our BWT range, so the last match range was maximal
                // or: we have exceeded the order of the graph (FPs if we go further)
                //     we have run over our parameter-defined MEM limit
                // record the last MEM
                match.begin = cursor+1;
                match.range = last_range;
                mems.push_back(match);
                // set up the next MEM using the parent node range
                // length of last MEM, which we use to update our end pointer for the next MEM
                size_t last_mem_length = match.end - match.begin;
                // get the parent suffix tree node corresponding to the parent of the last MEM's STNode
                gcsa::STNode parent = lcp->parent(last_range);
                // change the end for the next mem to reflect our step size
                size_t step_size = last_mem_length - parent.lcp();
                match.end = mems.back().end-step_size;
                // and set up the next MEM using the parent node range
                match.range = parent.range();
            }
        } else {
            // we are matching
            match.begin = cursor;
            // just step to the next position
            --cursor;
        }
    }
    // if we have a non-empty MEM at the end, record it
    if (match.end - match.begin > 0) mems.push_back(match);
    
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
            if (mem.match_count > 0 && (!hit_max || mem.match_count <= hit_max)) {
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
            if (mem.length() >= reseed_length
                && mem.match_count == 1
                // or if we only have one mem for the entire read (even if it may have many matches)
                || mems.size() == 1) {
                // reseed at midway between here and the min mem length and at the min mem length
                int reseed_to = mem.length() / 2;
                int reseeds = 0;
                while (reseeds == 0 && reseed_to >= min_mem_length) {
#ifdef debug_mapper
#pragma omp critical
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
    // print the matches
    /*
     for (auto& mem : mems) {
     cerr << mem << endl;
     }
     */
    // verify the matches (super costly at scale)
    /*
     #ifdef debug_mapper
     if (debug) { check_mems(mems); }
     #endif
     */
    return mems;
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
                                                     int reseed_below) {
#ifdef debug_mapper
#pragma omp critical
    {
        cerr << "find_mems: sequence ";
        for (auto iter = seq_begin; iter != seq_end; iter++) {
            cerr << *iter;
        }
        cerr << ", max mem length " << max_mem_length << ", min mem length " <<
        min_mem_length << ", reseed length " << reseed_length << endl;
    }
#endif

    if (!gcsa) {
        cerr << "error:[vg::Mapper] a GCSA2 index is required to query MEMs" << endl;
        exit(1);
    }
    
    if (min_mem_length > reseed_length && reseed_length) {
        cerr << "error:[vg::Mapper] minimimum reseed length for MEMs cannot be less than minimum MEM length" << endl;
        exit(1);
    }
    vector<MaximalExactMatch> mems;
    
    gcsa::range_type full_range = gcsa::range_type(0, gcsa->size() - 1);

    // an empty sequence matches the entire bwt
    if (seq_begin == seq_end) {
        mems.push_back(MaximalExactMatch(seq_begin, seq_end, full_range));
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
    
    // range of the last iteration
    gcsa::range_type last_range = full_range;
    
    // the temporary MEM we'll build up in this process
    MaximalExactMatch match(cursor, seq_end, full_range);
    
    // did we move the cursor or the end of the match last iteration?
    bool prev_iter_jumped_lcp = false;

    int filtered_mems = 0;
    int total_mems = 0;
    int max_lcp = 0;
    size_t mem_length = 0;
    vector<int> lcp_maxima;
    
    // loop maintains invariant that match.range contains the hits for seq[cursor+1:match.end]
    while (cursor >= seq_begin) {
        
        // break the MEM on N; which for DNA we assume is non-informative
        // this *will* match many places in assemblies, but it isn't helpful
        if (*cursor == 'N') {
            match.begin = cursor + 1;
            
            mem_length = match.length();
            
            if (mem_length >= min_mem_length) {

                mems.push_back(match);
                
#ifdef debug_mapper
#pragma omp critical
                {
                    vector<gcsa::node_type> locations;
                    gcsa->locate(match.range, locations);
                    cerr << "adding MEM " << match.sequence() << " at positions ";
                    for (auto nt : locations) {
                        cerr << make_pos_t(nt) << " ";
                    }
                    cerr << endl;
                }
#endif
            }
            
            match.end = cursor;
            match.range = full_range;
            --cursor;
            
            prev_iter_jumped_lcp = false;
            lcp_maxima.push_back(max_lcp);
            max_lcp = 0;
            // skip looking for matches since they are non-informative
            continue;
        }
        
        // hold onto our previous range
        last_range = match.range;
        
        // execute one step of LF mapping
        match.range = gcsa->LF(match.range, gcsa->alpha.char2comp[*cursor]);
        
        if (gcsa::Range::empty(match.range)
            || (max_mem_length && match.end - cursor > max_mem_length)
            || match.end-cursor > gcsa->order()) {
            
            // we've exhausted our BWT range, so the last match range was maximal
            // or: we have exceeded the order of the graph (FPs if we go further)
            // or: we have run over our parameter-defined MEM limit
            
            if (cursor + 1 == match.end) {
                // avoid getting caught in infinite loop when a single character mismatches
                // entire index (b/c then advancing the LCP doesn't move the search forward
                // at all, need to move the cursor instead)
                match.begin = cursor + 1;
                match.range = last_range;
                
                if (match.end - match.begin >= min_mem_length) {
                    mems.push_back(match);
                }
                
                match.end = cursor;
                match.range = full_range;
                --cursor;
                
                // don't reseed in empty MEMs
                prev_iter_jumped_lcp = false;
                lcp_maxima.push_back(max_lcp);
                max_lcp = 0;
            }
            else {
                match.begin = cursor + 1;
                match.range = last_range;
                mem_length = match.end - match.begin;
                // record the last MEM, but check to make sure were not actually still searching
                // for the end of the next MEM
                if (mem_length >= min_mem_length && !prev_iter_jumped_lcp) {
                    mems.push_back(match);
                    
#ifdef debug_mapper
#pragma omp critical
                    {
                        vector<gcsa::node_type> locations;
                        gcsa->locate(match.range, locations);
                        cerr << "adding MEM " << match.sequence() << " at positions ";
                        for (auto nt : locations) {
                            cerr << make_pos_t(nt) << " ";
                        }
                        cerr << endl;
                    }
#endif
                }
                
                // get the parent suffix tree node corresponding to the parent of the last MEM's STNode
                gcsa::STNode parent = lcp->parent(last_range);
                // set the MEM to be the longest prefix that is shared with another MEM
                match.end = match.begin + parent.lcp();
                // and set up the next MEM using the parent node range
                match.range = parent.range();
                // record our max lcp
                max_lcp = (int)parent.lcp();
                
                prev_iter_jumped_lcp = true;
                lcp_maxima.push_back(max_lcp);
                max_lcp = 0;
            }
        }
        else {
            prev_iter_jumped_lcp = false;
            max_lcp = (int)lcp->parent(match.range).lcp();
            ++mem_length;
            lcp_maxima.push_back(max_lcp);
            // just step to the next position
            --cursor;
        }
    }
    // TODO: is this where the bug with the duplicated MEMs is occurring? (when the prefix of a read
    // contains multiple non SMEM hits so that the iteration will loop through the LCP routine multiple
    // times before escaping out of the loop?
    
    // if we have a MEM at the beginning of the read, record it
    match.begin = seq_begin;
    mem_length = match.end - match.begin;
    if (mem_length >= min_mem_length) {
        max_lcp = (int)lcp->parent(match.range).lcp();
        mems.push_back(match);
        
#ifdef debug_mapper
#pragma omp critical
        {
            vector<gcsa::node_type> locations;
            gcsa->locate(match.range, locations);
            cerr << "adding MEM " << match.sequence() << " at positions ";
            for (auto nt : locations) {
                cerr << make_pos_t(nt) << " ";
            }
            cerr << endl;
        }
#endif
    }

    lcp_maxima.push_back(max_lcp);
    longest_lcp = *max_element(lcp_maxima.begin(), lcp_maxima.end());

    // filter weird MEMs ** that seem to occur when the input dBG to GCSA2 is made from path-only kmers
    mems.erase(std::remove_if(mems.begin(), mems.end(),
                              [&seq_begin, &seq_end](const MaximalExactMatch& mem)
                              { return mem.begin < seq_begin || mem.end > seq_end; }),
               mems.end());

    // fill the MEMs' node lists and indicate they are primary MEMs
    for (MaximalExactMatch& mem : mems) {
        mem.match_count = gcsa->count(mem.range);
        mem.primary = true;
        // if we aren't filtering on hit count, or if we have up to the max allowed hits
        if (mem.match_count > 0 && (!hit_max || mem.match_count <= hit_max)) {
            // extract the graph positions matching the range
            gcsa->locate(mem.range, mem.nodes);
            // it may be necessary to impose the cap on mem hits?
            if (hit_max > 0 && mem.nodes.size() > hit_max) {
                filtered_mems += mem.nodes.size();
                total_mems += mem.nodes.size();
                mem.nodes.clear();
            }
        }
    }

    if (reseed_length) {
        // get the sub_mem_and_parents
        vector<pair<MaximalExactMatch, vector<size_t> > > sub_mems;

        // a function that will implement the reseed according to various flags
        auto do_reseed = [&](int i) {
            if (fast_reseed) {
                if (use_diff_based_fast_reseed) {
                    find_sub_mems_fast(mems,
                                       i,
                                       (i+1 < mems.size() ? mems[i+1].begin : seq_begin),
                                       max<int>(ceil(fast_reseed_length_diff * mem_length),
                                                min_mem_length),
                                       sub_mems);
                }
                else {
                    find_sub_mems_fast(mems,
                                       i,
                                       (i+1 < mems.size() ? mems[i+1].begin : seq_begin),
                                       min_mem_length,
                                       sub_mems);
                }
            } else {
                find_sub_mems(mems,
                              i,
                              (i+1 < mems.size() ? mems[i+1].begin : seq_begin),
                              min_mem_length,
                              sub_mems);
            }
        };

        // run the reseeding
        for (int i = 0; i < mems.size(); ++i) {
            auto& mem = mems[i];
            // reseed when...
            if (mem.length() >= min_mem_length // our mem is greater than the min mem length (should be by default)
                && mem.length() >= reseed_length // is the right length to reseed
                && mem.nodes.size() // and wasn't filtered
                && (reseed_below == 0  // and has fewer hits than our threshold for reseeding, if there is a threshold
                    || mem.nodes.size() <= reseed_below)) {
                    do_reseed(i); // reseed the ith MEM
            }
        }

        // determine counts of matches
        for (pair<MaximalExactMatch, vector<size_t> >& sub_mem_and_parents : sub_mems) {
            // count in entire range, including parents
            sub_mem_and_parents.first.match_count = gcsa->count(sub_mem_and_parents.first.range);
            if (!include_parent_in_sub_mem_count) {
                // remove parents from count
                for (size_t parent_idx : sub_mem_and_parents.second) {
                    sub_mem_and_parents.first.match_count -= mems[parent_idx].match_count;
                }
            }
        }
        
        // fill MEMs with positions and set flag indicating they are submems
        for (auto& m : sub_mems) {
            auto& mem = m.first;
            mem.primary = false;
            if (mem.match_count > 0 && (!hit_max || mem.match_count <= hit_max)) {
                gcsa->locate(mem.range, mem.nodes);
                // it may be necessary to impose the cap on mem hits?
                if (hit_max > 0 && mem.nodes.size() > hit_max) {
                    filtered_mems += mem.nodes.size();
                    total_mems += mem.nodes.size();
                    mem.nodes.clear();
                }
            }
        }
        
        // combine the MEM and sub-MEM lists
        for (auto iter = sub_mems.begin(); iter != sub_mems.end(); iter++) {
            mems.push_back(std::move(iter->first));
        }
    }
    
    // return the MEMs in order along the read
    // TODO: there should actually be a linear time method to merge and order the sub-MEMs, since
    // they are ordered by the parent MEMs
    std::sort(mems.begin(), mems.end(), [](const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
        return m1.begin < m2.begin ? true : (m1.begin == m2.begin ? m1.end < m2.end : false);
    });

    std::for_each(mems.begin(), mems.end(), [&total_mems](const MaximalExactMatch& mem) { total_mems += mem.nodes.size(); });
    fraction_filtered = (double) filtered_mems / (double) total_mems;
    
    // remove non-unique MEMs
    // TODO: I think I already fixed this
    mems.erase(unique(mems.begin(), mems.end()), mems.end());
    // remove MEMs that are overlapping positionally (they may be redundant)
    return mems;
}

void BaseMapper::find_sub_mems(const vector<MaximalExactMatch>& mems,
                               int mem_idx,
                               string::const_iterator next_mem_end,
                               int min_mem_length,
                               vector<pair<MaximalExactMatch, vector<size_t>>>& sub_mems_out) {
    
    // get the most recently added MEM
    const MaximalExactMatch& mem = mems.back();
    
#ifdef debug_mapper
#pragma omp critical
    {
        cerr << "find_sub_mems: sequence ";
        for (auto iter = mem.begin; iter != mem.end; iter++) {
            cerr << *iter;
        }
        cerr << ", min mem length " << min_mem_length << endl;
    }
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
            
            if (sub_mem_end - sub_mem_begin >= min_mem_length && !prev_iter_jumped_lcp) {
                sub_mems_out.push_back(make_pair(MaximalExactMatch(sub_mem_begin, sub_mem_end, last_range),
                                                 vector<size_t>{(uint64_t)mem_idx}));
#ifdef debug_mapper
#pragma omp critical
                {
                    vector<gcsa::node_type> locations;
                    gcsa->locate(last_range, locations);
                    cerr << "adding sub-MEM ";
                    for (auto iter = sub_mem_begin; iter != sub_mem_end; iter++) {
                        cerr << *iter;
                    }
                    cerr << " at positions ";
                    for (auto nt : locations) {
                        cerr << make_pos_t(nt) << " ";
                    }
                    cerr << endl;
                    
                }
#endif
                // identify all previous MEMs that also contain this sub-MEM
                for (int64_t i = mem_idx - 1; i >= 0; --i) {
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
#pragma omp critical
                {
                    cerr << "minimally more frequent MEM is too short ";
                    for (auto iter = sub_mem_begin; iter != sub_mem_end; iter++) {
                        cerr << *iter;
                    }
                    cerr << endl;
                }
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
        sub_mems_out.push_back(make_pair(MaximalExactMatch(mem.begin, sub_mem_end, range),
                                         vector<size_t>{(uint64_t)mem_idx}));
#ifdef debug_mapper
#pragma omp critical
        {
            cerr << "adding sub-MEM ";
            for (auto iter = mem.begin; iter != sub_mem_end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
        }
#endif
        // note: this sub MEM is at the far left side of the parent MEM, so we don't need to
        // check whether earlier MEMs contain it as well
    }
#ifdef debug_mapper
    else {
#pragma omp critical
        {
            cerr << "minimally more frequent MEM is too short ";
            for (auto iter = mem.begin; iter != sub_mem_end; iter++) {
                cerr << *iter;
            }
            cerr << endl;
        }
    }
#endif
}

// TODO: rewrite roughly as follows
// //// idea is to check all the possible sub mems starting from the same position so we obviate the need to call count for every LF step
// initally take the range of the MEM, then do parent operation to find the next-longest match with a long range
// get ranges for all suffixes of the MEM to reseed
// for each suffix range we keep the starting position fixed then use parent queries to find candidates (avoids excessive count calls)
// --- hope is that the first call to parent jumps us close to the limit for minimum length ...
// goal is to get to roughly a linear number of LFs an parent()s per MEM we reseed

void BaseMapper::find_sub_mems_fast(const vector<MaximalExactMatch>& mems,
                                    int mem_idx,
                                    string::const_iterator next_mem_end,
                                    int min_sub_mem_length,
                                    vector<pair<MaximalExactMatch, vector<size_t>>>& sub_mems_out) {

    // get the most recently added MEM
    const MaximalExactMatch& mem = mems[mem_idx];
    
#ifdef debug_mapper
#pragma omp critical
    cerr << "find_sub_mems_fast: mem ";
    for (auto iter = mem.begin; iter != mem.end; iter++) {
        cerr << *iter;
    }
    cerr << ", min_sub_mem_length " << min_sub_mem_length << endl;
#endif
    
    // how many times does the parent MEM occur in the index?
    size_t parent_range_length = gcsa::Range::length(mem.range);
    
    // the end of the leftmost substring that is at least the minimum length and not contained
    // in the next SMEM
    string::const_iterator probe_string_end = mem.begin + min_sub_mem_length;
    if (probe_string_end <= next_mem_end) {
        probe_string_end = next_mem_end + 1;
    }
    // the leftmost possible index of the next sub-MEM
    string::const_iterator leftmost_bound = mem.begin;
    
    while (probe_string_end <= mem.end) {
        
        // locate the probe substring of length equal to the minimum length for a sub-MEM
        // that we are going to test to see if it's inside any sub-MEM
        string::const_iterator probe_string_begin = probe_string_end - min_sub_mem_length;
        
#ifdef debug_mapper
#pragma omp critical
        cerr << "probe string is mem[" << probe_string_begin - mem.begin << ":" << probe_string_end - mem.begin << "] ";
        for (auto iter = probe_string_begin; iter != probe_string_end; iter++) {
            cerr << *iter;
        }
        cerr << endl;
#endif
        
        // set up LF searching
        string::const_iterator cursor = probe_string_end - 1;
        gcsa::range_type range = gcsa::range_type(0, gcsa->size() - 1);
        
        // check if the probe substring is more frequent than the SMEM its contained in
        bool probe_string_more_frequent = true;
        while (cursor >= probe_string_begin) {
            
            range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
            
            if (gcsa::Range::length(range) <= parent_range_length) {
                probe_string_more_frequent = false;
                break;
            }
            
            cursor--;
        }
        
        if (probe_string_more_frequent) {
            // this is the prefix of a sub-MEM of length >= the minimum, now we need to
            // find its end using binary search
            
            if (probe_string_begin > leftmost_bound) {
                // the probe string is contained in a sub-MEM, but that doesn't mean it's been
                // extended as far as possible to the left
                
                // extend match until beginning of SMEM or until the end of the independent hit
                while (cursor >= leftmost_bound) {
                    gcsa::range_type last_range = range;
                    range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
                    
                    if (gcsa::Range::length(range) <= parent_range_length) {
                        range = last_range;
                        break;
                    }
                    
                    cursor--;
                }
                
                // mark this position as the beginning of the probe substring
                probe_string_begin = cursor + 1;
            }
            
            // inclusive interval that contains the past-the-last index of the sub-MEM
            string::const_iterator left_search_bound = probe_string_end;
            string::const_iterator right_search_bound = mem.end;
            
            // the match range of the longest prefix of the sub-MEM we've found so far (if we initialize
            // it here, the binary search is guaranteed to LF along the full sub-MEM in some iteration)
            gcsa::range_type sub_mem_range = range;
            
            // iterate until inteveral contains only one index
            while (right_search_bound > left_search_bound) {
                
                string::const_iterator middle = left_search_bound + (right_search_bound - left_search_bound + 1) / 2;
                
#ifdef debug_mapper
#pragma omp critical
                {
                    cerr << "checking extension mem[" << probe_string_begin - mem.begin << ":" << middle - mem.begin << "] ";
                    for (auto iter = probe_string_begin; iter != middle; iter++) {
                        cerr << *iter;
                    }
                    cerr << endl;
                }
#endif
                
                // set up LF searching
                cursor = middle - 1;
                range = gcsa::range_type(0, gcsa->size() - 1);
                
                // check if there is an independent occurrence of this substring outside of the SMEM
                bool contained_in_independent_match = true;
                while (cursor >= probe_string_begin) {

                    range = gcsa->LF(range, gcsa->alpha.char2comp[*cursor]);
                    
                    if (gcsa::Range::length(range) <= parent_range_length) {
                        // this probe is too long and it no longer is contained in the indendent hit
                        // that we detected
                        contained_in_independent_match = false;
                        break;
                    }
                    
                    cursor--;
                }
                
                if (contained_in_independent_match) {
                    
                    // the end of the sub-MEM must be here or to the right
                    left_search_bound = middle;
                    // update the range of matches (this is the longest match we've verified so far)
                    sub_mem_range = range;
                }
                else {
                    // the end of the sub-MEM must be to the left
                    right_search_bound = middle - 1;
                }
            }
            
#ifdef debug_mapper
#pragma omp critical
            cerr << "final sub-MEM is mem[" << probe_string_begin - mem.begin << ":" << right_search_bound - mem.begin << "] ";
            for (auto iter = probe_string_begin; iter != right_search_bound; iter++) {
                cerr << *iter;
            }
            cerr << endl;
#endif
            
            // record the sub-MEM
            sub_mems_out.push_back(make_pair(MaximalExactMatch(probe_string_begin, right_search_bound, sub_mem_range),
                                             vector<size_t>{(uint64_t)mem_idx}));
            
            // identify all previous MEMs that also contain this sub-MEM
            for (int64_t i = mem_idx - 1; i >= 0; --i) { //((int64_t) mems.size()) - 2; i >= 0; i--) {
                if (probe_string_begin >= mems[i].begin) {
                    // contained in next MEM, add its index to sub MEM's list of parents
                    sub_mems_out.back().second.push_back(i);
                }
                else {
                    // not contained in the next MEM, cannot be contained in earlier ones
                    break;
                }
            }
            
            // the closest possible independent probe string will now occur one position
            // to the right of this sub-MEM
            probe_string_end = right_search_bound + 1;
            
            // any sub-MEMs will need to have their start position at least one position to the right
            // of this one
            leftmost_bound = probe_string_begin + 1;
            
        }
        else {
            // we've found a suffix of the probe substring that is only contained inside the
            // parent SMEM, so we can move far enough to the right that the next probe substring
            // will not contain it
            
            probe_string_end = cursor + min_sub_mem_length + 1;
        }
    }
}

void BaseMapper::first_hit_positions_by_index(MaximalExactMatch& mem,
                                              vector<set<pos_t>>& positions_by_index_out) {
    // find the hit to the first index in the parent MEM's range
    vector<gcsa::node_type> all_first_hits;
    gcsa->locate(mem.range.first, all_first_hits, true, false);
    
    // find where in the graph the first hit of the parent MEM is at each index
    mem_positions_by_index(mem, make_pos_t(all_first_hits[0]), positions_by_index_out);
    
    // in case the first hit occurs in more than one place, accumulate all the hits
    if (all_first_hits.size() > 1) {
        for (size_t i = 1; i < all_first_hits.size(); i++) {
            vector<set<pos_t>> temp_positions_by_index;
            mem_positions_by_index(mem, make_pos_t(all_first_hits[i]),
                                   temp_positions_by_index);
            
            for (size_t i = 0; i < positions_by_index_out.size(); i++) {
                for (const pos_t& pos : temp_positions_by_index[i]) {
                    positions_by_index_out[i].insert(pos);
                }
            }
        }
    }
}

void BaseMapper::fill_nonredundant_sub_mem_nodes(vector<MaximalExactMatch>& parent_mems,
                                                 vector<pair<MaximalExactMatch, vector<size_t> > >::iterator sub_mem_records_begin,
                                                 vector<pair<MaximalExactMatch, vector<size_t> > >::iterator sub_mem_records_end) {
    
    
    // for each MEM, a vector of the positions that it touches at each index along the MEM
    vector<vector<set<pos_t>>> positions_by_index(parent_mems.size());
    
    for (auto iter = sub_mem_records_begin; iter != sub_mem_records_end; iter++) {
        
        pair<MaximalExactMatch, vector<size_t> >& sub_mem_and_parents = *iter;
        
        MaximalExactMatch& sub_mem = sub_mem_and_parents.first;
        vector<size_t>& parent_idxs = sub_mem_and_parents.second;
        
        // how many total hits does each parent MEM have?
        vector<size_t> num_parent_hits;
        // positions their first hits of the parent MEM takes at the start position of the sub-MEM
        vector<set<pos_t>*> first_parent_mem_hit_positions;
        for (size_t parent_idx : parent_idxs) {
            // get the parent MEM
            MaximalExactMatch& parent_mem = parent_mems[parent_idx];
            num_parent_hits.push_back(gcsa->count(parent_mem.range));
            
            if (positions_by_index[parent_idx].empty()) {
                // the parent MEM's positions by index haven't been calculated yet, so do it
                
                first_hit_positions_by_index(parent_mem, positions_by_index[parent_idx]);
                
            }
            
            // the index along the parent MEM that sub MEM starts
            size_t offset = sub_mem.begin - parent_mem.begin;
            first_parent_mem_hit_positions.push_back(&(positions_by_index[parent_idx][offset]));
        }
        
        for (gcsa::size_type i = sub_mem.range.first; i <= sub_mem.range.second; i++) {
            
            // add the locations of the hits, but do not remove duplicates yet
            vector<gcsa::node_type> hits;
            gcsa->locate(i, hits, true, false);
            
            // the number of subsequent hits (including these) that are inside a parent MEM
            size_t parent_hit_jump = 0;
            for (gcsa::node_type node : hits) {
                // look for the hit in each parent MEM
                for (size_t j = 0; j < first_parent_mem_hit_positions.size(); j++) {
                    if (first_parent_mem_hit_positions[j]->count(make_pos_t(node))) {
                        // this hit is also a node on a path of the first occurrence of the parent MEM
                        // that means that this is the first index of the sub-range that corresponds
                        // to the parent MEM's hits
                        
                        // calculate how many more positions to jump
                        parent_hit_jump = num_parent_hits[j];
                        break;
                    }
                }
            }
            
            if (parent_hit_jump > 0) {
                // we're at the start of an interval of parent hits, skip the rest of it
                i += (parent_hit_jump - 1);
            }
            else {
                // these are nonredundant sub MEM hits, add them
                for (gcsa::node_type node : hits) {
                    sub_mem.nodes.push_back(node);
                }
            }
        }
        
        // remove duplicates (copied this functionality from the gcsa locate function, but
        // I don't actually know what it's purpose is)
        gcsa::removeDuplicates(sub_mem.nodes, false);
    }
}

void BaseMapper::mem_positions_by_index(MaximalExactMatch& mem, pos_t hit_pos,
                                        vector<set<pos_t>>& positions_by_index_out) {
    
    // this is a specialized DFS that keeps track of both the distance along the MEM
    // and the position(s) in the graph in the stack by adding all of the next reachable
    // positions in a layer (i.e. vector) in the stack at the end of each iteration.
    // it also keeps track of whether a position in the graph matched to a position along
    // the MEM can potentially be extended to the full MEM to avoid combinatorially checking
    // all paths through bubbles
    
    size_t mem_length = std::distance(mem.begin, mem.end);
    
    // indicates a pairing of this graph position and this MEM index could be extended to a full match
    positions_by_index_out.clear();
    positions_by_index_out.resize(mem_length);
    
    // indicates a pairing of this graph position and this MEM index could not be extended to a full match
    vector<set<pos_t>> false_pos_by_mem_index(mem_length);
    
    // each record indicates the next edge index to traverse, the number of edges that
    // cannot reach a MEM end, and the positions along each edge out
    vector<pair<pair<size_t, size_t>, vector<pos_t> > > pos_stack;
    pos_stack.push_back(make_pair(make_pair((size_t) 0 , (size_t) 0), vector<pos_t>{hit_pos}));
    
    while (!pos_stack.empty()) {
        size_t mem_idx = pos_stack.size() - 1;
        
        // which edge are we going to search out of this node next?
        size_t next_idx = pos_stack.back().first.first;
        
        if (next_idx >= pos_stack.back().second.size()) {
            // we have traversed all of the edges out of this position
            
            size_t num_misses = pos_stack.back().first.second;
            bool no_full_matches_possible = (num_misses == pos_stack.back().second.size());
            
            // backtrack to previous node
            pos_stack.pop_back();
            
            // if necessary, mark the edge into this node as a miss
            if (no_full_matches_possible && !pos_stack.empty()) {
                // all of the edges out failed to reach the end of a MEM, this position is a dead end
                
                // get the position that traversed into the layer we just popped off
                pos_t prev_graph_pos = pos_stack.back().second[pos_stack.back().first.first - 1];
                
                // unlabel this node as a potential hit and instead mark it as a miss
                positions_by_index_out[mem_idx].erase(prev_graph_pos);
                false_pos_by_mem_index[mem_idx].insert(prev_graph_pos);
                
                // increase the count of misses in this layer
                pos_stack.back().first.second++;
            }
            
            // skip the forward search on this iteration
            continue;
        }
        
        // increment to the next edge
        pos_stack.back().first.first++;
        
        pos_t graph_pos = pos_stack.back().second[next_idx];
        
        
        // did we already find a MEM through this position?
        if (positions_by_index_out[mem_idx].count(graph_pos)) {
            // we don't need to check the same MEM suffix again
            continue;
        }
        
        // did we already determine that you can't reach a MEM through this position?
        if (false_pos_by_mem_index[mem_idx].count(graph_pos)) {
            // increase the count of misses in this layer
            pos_stack.back().first.second++;
            
            // we don't need to check the same MEM suffix again
            continue;
        }
        
        // does this graph position match the MEM?
        if (*(mem.begin + mem_idx) != xg_pos_char(graph_pos, xindex)) {
            // mark this node as a miss
            false_pos_by_mem_index[mem_idx].insert(graph_pos);
            
            // increase the count of misses in this layer
            pos_stack.back().first.second++;
        }
        else {
            // mark this node as a potential hit
            positions_by_index_out[mem_idx].insert(graph_pos);
            
            // are we finished with the MEM?
            if (mem_idx < mem_length - 1) {
                
                // add a layer onto the stack for all of the edges out
                pos_stack.push_back(make_pair(make_pair((size_t) 0 , (size_t) 0),
                                              vector<pos_t>()));
                
                // fill the layer with the next positions
                vector<pos_t>& nexts = pos_stack.back().second;
                for (const pos_t& next_graph_pos : positions_bp_from(graph_pos, 1, false)) {
                    nexts.push_back(next_graph_pos);
                }
            }
        }
    }
}
    
    
set<pos_t> BaseMapper::positions_bp_from(pos_t pos, int distance, bool rev) {
    return xg_positions_bp_from(pos, distance, rev, xindex);
}

void BaseMapper::check_mems(const vector<MaximalExactMatch>& mems) {
    for (auto mem : mems) {
#ifdef debug_mapper
#pragma omp critical
        cerr << "checking MEM: " << mem.sequence() << endl;
#endif
        // TODO: fix this for sub-MEMs
        if (sequence_positions(mem.sequence()) != gcsa_nodes_to_positions(mem.nodes)) {
            cerr << "SMEM failed! " << mem.sequence()
            << " expected " << sequence_positions(mem.sequence()).size() << " hits "
            << "but found " << gcsa_nodes_to_positions(mem.nodes).size()
            << "(aside: this consistency check is broken for sub-MEMs, oops)" << endl;
        }
    }
}
    
char BaseMapper::pos_char(pos_t pos) {
    return xg_pos_char(pos, xindex);
}

map<pos_t, char> BaseMapper::next_pos_chars(pos_t pos) {
    return xg_next_pos_chars(pos, xindex);
}
    
set<pos_t> BaseMapper::sequence_positions(const string& seq) {
    gcsa::range_type gcsa_range = gcsa->find(seq);
    std::vector<gcsa::node_type> gcsa_nodes;
    gcsa->locate(gcsa_range, gcsa_nodes);
    return gcsa_nodes_to_positions(gcsa_nodes);
}
    
void BaseMapper::set_alignment_threads(int new_thread_count) {
    alignment_threads = new_thread_count;
}

bool BaseMapper::has_fixed_fragment_length_distr() {
    return fragment_length_distr.is_finalized();
}

void BaseMapper::force_fragment_length_distr(double mean, double stddev) {
    fragment_length_distr.force_parameters(mean, stddev);
}
    
BaseAligner* BaseMapper::get_aligner(bool have_qualities) const {
    return (have_qualities && adjust_alignments_for_base_quality) ?
        (BaseAligner*) qual_adj_aligner :
        (BaseAligner*) regular_aligner;
}

QualAdjAligner* BaseMapper::get_qual_adj_aligner() const {
    assert(qual_adj_aligner != nullptr);
    return qual_adj_aligner;
}

Aligner* BaseMapper::get_regular_aligner() const {
    assert(regular_aligner != nullptr);
    return regular_aligner;
}

void BaseMapper::clear_aligners(void) {
    delete qual_adj_aligner;
    delete regular_aligner;
    qual_adj_aligner = nullptr;
    regular_aligner = nullptr;
}

void BaseMapper::init_aligner(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus) {
    // hacky, find max score so that scaling doesn't change score
    int8_t max_score = match;
    if (mismatch > max_score) max_score = mismatch;
    if (gap_open > max_score) max_score = gap_open;
    if (gap_extend > max_score) max_score = gap_extend;
    
    double gc_content = estimate_gc_content();
    
    qual_adj_aligner = new QualAdjAligner(match, mismatch, gap_open, gap_extend, full_length_bonus,
                                          max_score, 255, gc_content);
    regular_aligner = new Aligner(match, mismatch, gap_open, gap_extend, full_length_bonus);
}
    
double BaseMapper::estimate_gc_content(void) {
    
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
    if (gcsa) {
        size_t length = gcsa::Range::length(gcsa->find(string("")));
        return ceil(- (log(1.0 - pow(pow(1.0-chance_random, -1), (-1.0/length))) / log(4.0)));
    } else {
        return 0;
    }
}
    
void BaseMapper::set_alignment_scores(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus) {
    // clear the existing aligners and recreate them
    if (regular_aligner || qual_adj_aligner) {
        clear_aligners();
    }
    init_aligner(match, mismatch, gap_open, gap_extend, full_length_bonus);
}
    
void BaseMapper::set_fragment_length_distr_params(size_t maximum_sample_size, size_t reestimation_frequency,
                                                  double robust_estimation_fraction) {
    
    if (fragment_length_distr.is_finalized()) {
        cerr << "warning:[vg::Mapper] overwriting a fragment length distribution that has already been estimated" << endl;
    }
    
    fragment_length_distr = FragmentLengthDistribution(maximum_sample_size, reestimation_frequency,
                                                       robust_estimation_fraction);
}
    
Mapper::Mapper(xg::XG* xidex,
               gcsa::GCSA* g,
               gcsa::LCPArray* a) :
    BaseMapper(xidex, g, a)
    , thread_extension(10)
    , context_depth(1)
    , max_multimaps(1)
    , min_multimaps(4)
    , max_attempts(0)
    , min_cluster_length(0)
    , softclip_threshold(0)
    , max_softclip_iterations(10)
    , min_identity(0)
    , debug(false)
    , max_target_factor(128)
    , max_query_graph_ratio(128)
    , extra_multimaps(512)
    , band_multimaps(4)
    , always_rescue(false)
    , max_cluster_mapping_quality(1024)
    , use_cluster_mq(false)
    , simultaneous_pair_alignment(true)
    , drop_chain(0.2)
    , mq_overlap(0.2)
    , mate_rescues(0)
    , maybe_mq_threshold(10)
    , min_banded_mq(0)
    , max_band_jump(0)
    , identity_weight(2)
    , pair_rescue_hang_threshold(0.7)
    , pair_rescue_retry_threshold(0.5)
{
    
}

Mapper::Mapper(void) : BaseMapper() {
    // Nothing to do. Default constructed and can't really do anything.
}

Mapper::~Mapper(void) {
    //  Nothing to do, all managed memory handled in parent class
}

double Mapper::graph_entropy(void) {
    const size_t seq_bytes = xindex->sequence_bit_size() / 8;
    char* seq = (char*) xindex->sequence_data();
    return entropy(seq, seq_bytes);
}

// todo add options for aligned global and pinned
Alignment Mapper::align_to_graph(const Alignment& aln,
                                 Graph& graph,
                                 size_t max_query_graph_ratio,
                                 bool traceback,
                                 bool pinned_alignment,
                                 bool pin_left,
                                 bool banded_global) {
    // check if we need to make a vg graph to handle this graph
    if (!is_id_sortable(graph) || has_inversion(graph)) {
        VG vg;
        vg.extend(graph);
        if (aln.quality().empty() || !adjust_alignments_for_base_quality) {
            return vg.align(aln,
                            get_regular_aligner(),
                            traceback,
                            assume_acyclic,
                            max_query_graph_ratio,
                            pinned_alignment,
                            pin_left,
                            banded_global,
                            0, // band padding override
                            aln.sequence().size());
        } else {
            return vg.align_qual_adjusted(aln,
                                          get_qual_adj_aligner(),
                                          traceback,
                                          assume_acyclic,
                                          max_query_graph_ratio,
                                          pinned_alignment,
                                          pin_left,
                                          banded_global,
                                          0, // band padding override
                                          aln.sequence().size());
        }
    } else {
        // we've got an id-sortable graph and we can directly align with gssw
        Alignment aligned = aln;
        if (banded_global) {
            size_t max_span = aln.sequence().size();
            size_t band_padding_override = 0;
            bool permissive_banding = (band_padding_override == 0);
            size_t band_padding = permissive_banding ? max(max_span, (size_t) 1) : band_padding_override;
            get_aligner(!aln.quality().empty())->align_global_banded(aligned, graph, band_padding, false);
        } else if (pinned_alignment) {
            get_aligner(!aln.quality().empty())->align_pinned(aligned, graph, pin_left);
        } else {
            get_aligner(!aln.quality().empty())->align(aligned, graph, traceback, false);
        }
        return aligned;
    }
}

Alignment Mapper::align(const string& seq, int kmer_size, int stride, int max_mem_length, int band_width) {
    Alignment aln;
    aln.set_sequence(seq);
    return align(aln, kmer_size, stride, max_mem_length, band_width);
}

map<string, double> Mapper::alignment_mean_path_positions(const Alignment& aln, bool first_hit_only) {
    map<string, double> mean_pos;
    // Alignments are consistent if their median node id positions are within the fragment_size
    
    // We need the sets of nodes visited by each alignment
    set<id_t> ids;
    
    for(size_t i = 0; i < aln.path().mapping_size(); i++) {
        // Collect all the unique nodes visited by the first algnment
        ids.insert(aln.path().mapping(i).position().node_id());
    }
    map<string, map<int64_t, vector<id_t> > > node_positions;
    for(auto id : ids) {
        for (auto& ref : node_positions_in_paths(gcsa::Node::encode(id, 0))) {
            auto& name = ref.first;
            for (auto pos : ref.second) {
                node_positions[name][pos].push_back(id);
            }
        }
        // just get the first one
        if (first_hit_only && node_positions.size()) break;
    }
    // get mean mapping positions
    for (auto& ref : node_positions) {
        int idscount = 0;
        double idssum = 0;
        for (auto& p : ref.second) {
            for (auto& n : p.second) {
                auto pos = p.first + get_node_length(n)/2;
                if (ids.count(n)) {
                    idscount++;
                    idssum += pos;
                }
            }
        }
        mean_pos[ref.first] = idssum/idscount;
    }
    return mean_pos;
}
    
map<string, size_t> Mapper::alignment_initial_path_positions(const Alignment& aln) {
#ifdef debug_mapper
    cerr << "finding initial path positions for read " << aln.name() << endl;
#endif
    map<string, size_t> to_return;
    for (size_t i = 0; i < aln.path().mapping_size(); i++){
        const Position& pos = aln.path().mapping(i).position();
        map<string, vector<size_t>> path_positions = xindex->position_in_paths(pos.node_id(), pos.is_reverse(), pos.offset());
        for (const pair<string, vector<size_t>>& path_record : path_positions) {
            if (!to_return.count(path_record.first)) {
#ifdef debug_mapper
                cerr << "found first occurrence of path " << path_record.first << " on " << i << "-th mapping with position " << pb2json(aln.path().mapping(i).position()) << ", which occurs " << path_record.second.size() << " times on this path:" << endl;
                for (auto off : path_record.second) {
                    cerr << "\t" << off << endl;
                }
#endif
                to_return[path_record.first] = *min_element(path_record.second.begin(), path_record.second.end());
            }
        }
    }
    return to_return;
}

void Mapper::annotate_with_initial_path_positions(Alignment& aln) {
    map<string, size_t> init_path_positions = alignment_initial_path_positions(aln);
    for (const pair<string, size_t>& pos_record : init_path_positions) {
        Position* refpos = aln.add_refpos();
        refpos->set_name(pos_record.first);
        refpos->set_offset(pos_record.second);
    }
}

pos_t Mapper::likely_mate_position(const Alignment& aln, bool is_first_mate) {
    bool aln_is_rev = aln.path().mapping(0).position().is_reverse();
    int64_t aln_pos = approx_alignment_position(aln);
    //if (debug) cerr << "aln pos " << aln_pos << endl;
    // can't find the alignment position
    if (aln_pos < 0) return make_pos_t(0, false, 0);
    bool same_orientation = frag_stats.cached_fragment_orientation;
    bool forward_direction = frag_stats.cached_fragment_direction;
    int64_t delta = frag_stats.cached_fragment_length_mean;
    // which way is our delta?
    // we are on the forward strand
    id_t target;
    if (forward_direction) {
        if (is_first_mate) {
            if (!aln_is_rev) {
                target = node_approximately_at(aln_pos + delta);
            } else {
                target = node_approximately_at(aln_pos - delta);
            }
        } else {
            if (!aln_is_rev) {
                target = node_approximately_at(aln_pos + delta);
            } else {
                target = node_approximately_at(aln_pos - delta);
            }
        }
    } else {
        if (is_first_mate) {
            if (!aln_is_rev) {
                target = node_approximately_at(aln_pos - delta);
            } else {
                target = node_approximately_at(aln_pos + delta);
            }
        } else {
            if (!aln_is_rev) {
                target = node_approximately_at(aln_pos - delta);
            } else {
                target = node_approximately_at(aln_pos + delta);
            }
        }
    }
    if (same_orientation) {
        return make_pos_t(target, aln_is_rev, 0);
    } else {
        return make_pos_t(target, !aln_is_rev, 0);
    }
    /*
        && !aln_is_rev) {
    } else if (!same_direction && aln_is_rev) {
        target = (is_first_mate ? node_approximately_at(aln_pos + delta)
                  : node_approximately_at(aln_pos - delta));
    } else if (same_direction && aln_is_rev
               || !same_direction && !aln_is_rev) {
        target = (is_first_mate ? node_approximately_at(aln_pos - delta)
                  : node_approximately_at(aln_pos + delta));
    }
    */
    //bool target_is_rev = (same_orientation ? aln_is_rev : !aln_is_rev);
    //return make_pos_t(target, target_is_rev, 0);
}

pair<bool, bool> Mapper::pair_rescue(Alignment& mate1, Alignment& mate2, int match_score, int full_length_bonus, bool traceback) {
    auto pair_sig = signature(mate1, mate2);
    // bail out if we can't figure out how far to go
    bool rescued1 = false;
    bool rescued2 = false;
    if (!frag_stats.fragment_size) return make_pair(false, false);
    double hang_threshold = pair_rescue_hang_threshold;
    double retry_threshold = pair_rescue_retry_threshold;
    double perfect_score = mate1.sequence().size() * match_score + full_length_bonus;
    bool consistent = (mate1.score() > 0 && mate2.score() > 0 && pair_consistent(mate1, mate2, 0.0001));
    //double retry_threshold = mate1.sequence().size() * aligner->match * 0.3;
    // based on our statistics about the alignments
    // get the subgraph overlapping the likely candidate position of the second alignment
    bool rescue_off_first = false;
    bool rescue_off_second = false;
    double mate1_id = (double) mate1.score() / perfect_score;
    double mate2_id = (double) mate2.score() / perfect_score;
    pos_t mate_pos;
    //cerr << "---------------------------" << pb2json(mate1) << endl << pb2json(mate2) << endl;
    //if (debug) cerr << "pair rescue: mate1 " << signature(mate1) << " " << mate1_id << " mate2 " << signature(mate2) << " " << mate2_id << " consistent? " << consistent << endl;
    //cerr << "---------------------------" << endl;
    //if (debug) cerr << "mate1: " << pb2json(mate1) << endl;
    //if (debug) cerr << "mate2: " << pb2json(mate2) << endl;
    if (mate1_id > mate2_id && mate1_id > hang_threshold && mate2_id <= retry_threshold && !consistent) {
        // retry off mate1
#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) cerr << "Rescue read 2 off of read 1" << endl;
        }
#endif
        rescue_off_first = true;
        // record id and direction to second mate
        mate_pos = likely_mate_position(mate1, true);
    } else if (mate2_id > mate1_id && mate2_id > hang_threshold && mate1_id <= retry_threshold && !consistent) {
        // retry off mate2
#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) cerr << "Rescue read 1 off of read 2" << endl;
        }
#endif
        rescue_off_second = true;
        // record id and direction to second mate
        mate_pos = likely_mate_position(mate2, false);
    } else {
        return make_pair(false, false);
    }
#ifdef debug_mapper
#pragma omp critical
    {
        if (debug) cerr << "aiming for " << mate_pos << endl;
    }
#endif
    if (id(mate_pos) == 0) return make_pair(false, false); // can't rescue because the selected mate is unaligned
    int get_at_least = (!frag_stats.cached_fragment_length_mean ? frag_stats.fragment_max
                        : max((int)frag_stats.cached_fragment_length_stdev * 6 + mate1.sequence().size(),
                              mate1.sequence().size() * 4));
    Graph graph = xindex->graph_context_id(mate_pos, get_at_least/2);
    graph.MergeFrom(xindex->graph_context_id(reverse(mate_pos, get_node_length(id(mate_pos))), get_at_least/2));
    sort_by_id_dedup_and_clean(graph);
    //if (debug) cerr << "rescue got graph " << pb2json(graph.graph) << endl;
    // if we're reversed, align the reverse sequence and flip it back
    // align against it
    if (rescue_off_first) {
        Alignment aln2 = align_maybe_flip(mate2, graph, is_rev(mate_pos), traceback);
#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) cerr << "aln2 score/ident vs " << aln2.score() << "/" << aln2.identity()
                            << " vs " << mate2.score() << "/" << mate2.identity() << endl;
        }
#endif
        if (aln2.score() > mate2.score() && (double)aln2.score()/perfect_score > retry_threshold && pair_consistent(mate1, aln2, 0.0001)) {
            //cerr << "rescued aln2" << endl;
            mate2 = aln2;
            rescued2 = true;
        } else {
            return make_pair(false, false);
        }
    } else if (rescue_off_second) {
        Alignment aln1 = align_maybe_flip(mate1, graph, is_rev(mate_pos), traceback);
#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) cerr << "aln1 score/ident vs " << aln1.score() << "/" << aln1.identity()
                            << " vs " << mate1.score() << "/" << mate1.identity() << endl;
        }
#endif
        if (aln1.score() > mate1.score() && (double)aln1.score()/perfect_score > retry_threshold && pair_consistent(aln1, mate2, 0.0001)) {
            //cerr << "rescued aln1" << endl;
            mate1 = aln1;
            rescued1 = true;
        } else {
            return make_pair(false, false);
        }
    }
    // if the new alignment is better
    // set the old alignment to it
    return make_pair(rescued1, rescued2);
}

Alignment Mapper::realign_from_start_position(const Alignment& aln, int extra, int iteration) {
    if (!aln.path().mapping_size()) return aln;
    if (iteration > 3) return aln;
    int score = aln.score();
    pos_t pos = make_pos_t(aln.path().mapping(0).position());
    int get_at_least = 1.61803 * aln.sequence().size() + extra;
    Graph graph = xindex->graph_context_id(pos, get_at_least/1.61803);
    graph.MergeFrom(xindex->graph_context_id(reverse(pos, get_node_length(id(pos))), get_at_least*(1-1/1.61803)));
    sort_by_id_dedup_and_clean(graph);
    Alignment result = align_maybe_flip(aln, graph, is_rev(pos), true);
    if (result.score() >= score) {
        return result;
    } else {
        return realign_from_start_position(aln, 2*extra, ++iteration);
    }
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

bool Mapper::pair_consistent(const Alignment& aln1,
                             const Alignment& aln2,
                             double pval) {
    if (!(aln1.score() && aln2.score())) return false;
    bool length_ok = false;
    if (aln1.fragment_size() == 0 || aln2.fragment_size() == 0 || aln1.fragment_size() != aln2.fragment_size()) {
        // use the approximate distance
        int len = approx_fragment_length(aln1, aln2);
        if (frag_stats.fragment_size && len > 0 && (pval > 0 && frag_stats.fragment_length_pval(len) > pval
                                                    || len < frag_stats.fragment_size)
            || !frag_stats.fragment_size && len > 0 && len < frag_stats.fragment_max) {
            length_ok = true;
        }
    } else {
        // use the distance induced by the graph paths
        assert(aln1.fragment_size() == aln2.fragment_size());
        for (size_t i = 0; i < aln1.fragment_size(); ++i) {
            int len = abs(aln1.fragment(i).length());
            if (frag_stats.fragment_size && len > 0 && (pval > 0 && frag_stats.fragment_length_pval(len) > pval
                                             || len < frag_stats.fragment_size)
                || !frag_stats.fragment_size && len > 0 && len < frag_stats.fragment_max) {
                length_ok = true;
                break;
            }
        }
    }
    bool aln1_is_rev = aln1.path().mapping(0).position().is_reverse();
    bool aln2_is_rev = aln2.path().mapping(0).position().is_reverse();
    bool same_orientation = frag_stats.cached_fragment_orientation;
    bool orientation_ok = same_orientation && aln1_is_rev == aln2_is_rev
        || !same_orientation && aln1_is_rev != aln2_is_rev;
    return length_ok && orientation_ok;
}

pair<vector<Alignment>, vector<Alignment>> Mapper::align_paired_multi(
    const Alignment& first_mate,
    const Alignment& second_mate,
    bool& queued_resolve_later,
    int max_mem_length,
    bool only_top_scoring_pair,
    bool retrying) {

    Alignment read1;
    read1.set_name(first_mate.name());
    read1.set_sequence(first_mate.sequence());
    read1.set_quality(first_mate.quality());
    Alignment read2;
    read2.set_name(second_mate.name());
    read2.set_sequence(second_mate.sequence());
    read2.set_quality(second_mate.quality());

    double avg_node_len = average_node_length();

    auto aligner = get_aligner(!read1.quality().empty());
    int8_t match = aligner->match;
    int8_t gap_extension = aligner->gap_extension;
    int8_t gap_open = aligner->gap_open;
    int8_t full_length_bonus = aligner->full_length_bonus;

    int total_multimaps = max(max_multimaps, extra_multimaps);
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
             << frag_stats.cached_fragment_orientation << ", "
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
                                                     false, false, false, 0);
    vector<MaximalExactMatch> mems2 = find_mems_deep(read2.sequence().begin(),
                                                     read2.sequence().end(),
                                                     longest_lcp2,
                                                     fraction_filtered2,
                                                     max_mem_length,
                                                     min_mem_length,
                                                     mem_reseed_length,
                                                     false, false, false, 0);

    double mq_cap1, mq_cap2;
    mq_cap1 = mq_cap2 = max_mapping_quality;

    int total_mem_length1 = 0;
    for (auto& mem : mems1) total_mem_length1 += mem.length() * mem.nodes.size();
    int total_mem_length2 = 0;
    for (auto& mem : mems2) total_mem_length2 += mem.length() * mem.nodes.size();

    int mem_max_length1 = 0;
    for (auto& mem : mems1) if (mem.primary && mem.match_count) mem_max_length1 = max(mem_max_length1, (int)mem.length());
    double maybe_mq1 = estimate_max_possible_mapping_quality(read1.sequence().size(),
                                                             read1.sequence().size()/max(1.0, (double)mem_max_length1),
                                                             read1.sequence().size()/longest_lcp1);
    int mem_max_length2 = 0;
    for (auto& mem : mems2) if (mem.primary && mem.match_count) mem_max_length2 = max(mem_max_length2, (int)mem.length());
    double maybe_mq2 = estimate_max_possible_mapping_quality(read2.sequence().size(),
                                                             read2.sequence().size()/max(1.0, (double)mem_max_length2),
                                                             read2.sequence().size()/longest_lcp2);
    // use the estimated mapping quality to avoid hard work when the results are likely noninformative
    double maybe_min = min(maybe_mq1, maybe_mq2);
    if (maybe_min < maybe_mq_threshold) {
        mq_cap1 = maybe_min;
        mq_cap2 = maybe_min;
    }
    total_multimaps = max(min_multimaps, (int)round(total_multimaps/min(maybe_mq1,maybe_mq2)));
    if (debug) cerr << "maybe_mq1 " << read1.name() << " " << maybe_mq1 << " " << total_multimaps << " " << mem_max_length1 << " " << longest_lcp1 << endl;
    if (debug) cerr << "maybe_mq2 " << read2.name() << " " << maybe_mq2 << " " << total_multimaps << " " << mem_max_length2 << " " << longest_lcp2 << endl;

//#ifdef debug_mapper
#pragma omp critical
    {
        if (debug) cerr << "mems for read 1 " << mems_to_json(mems1) << endl;
    }
//#endif
//#ifdef debug_mapper
#pragma omp critical
    {
        if (debug) cerr << "mems for read 2 " << mems_to_json(mems2) << endl;
    }
//#endif


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
            int64_t approx_dist = mem_min_distance(m1, m2); //graph_mixed_distance_estimate(m1_pos, m2_pos, 0); // use graph/path based estimate here
#ifdef debug_mapper
#pragma omp critical
            {
                if (debug) cerr << "between fragment approx distance " << approx_dist << endl;
            }
#endif
            if (approx_dist >= max_length) {
                // Seem to be too far appart
#ifdef debug_mapper
#pragma omp critical
                {
                    if (debug) cerr << "seem too far apart by approx_dist" << endl;
                }
#endif
                return -std::numeric_limits<double>::max();
            } else {

                if (approx_dist >= max_length) {
#ifdef debug_mapper
#pragma omp critical
                    {
                        if (debug) cerr << "too far apart by path distance" << endl;
                    }
#endif
                    return -std::numeric_limits<double>::max();
                } else if (frag_stats.fragment_size) {
                    // exclude cases that don't match our model
                    if (!frag_stats.cached_fragment_orientation
                        && is_rev(m1_pos) == is_rev(m2_pos)
                        || frag_stats.cached_fragment_orientation
                        && is_rev(m1_pos) != is_rev(m2_pos)
                        || approx_dist > frag_stats.fragment_size) {
#ifdef debug_mapper
#pragma omp critical
                        {
                            if (debug) cerr << "bad orientations or dist of " << approx_dist
                                << " beyond fragment_size of " << frag_stats.fragment_size << endl;
                        }
#endif
                        return -std::numeric_limits<double>::max();
                    } else {
#ifdef debug_mapper
#pragma omp critical
                        {
                            if (debug) cerr << "OK with known fragment size" << endl;
                        }
#endif
                        return frag_stats.fragment_length_pval(approx_dist) * (m1.length() + m2.length());
                    }
                } else {
#ifdef debug_mapper
#pragma omp critical
                    {
                        if (debug) cerr << "OK with no fragment size" << endl;
                    }
#endif
                    return 1.0/approx_dist * (m1.length() + m2.length());
                }
            }
        } else if (m1.fragment > m2.fragment) {
            // don't allow going backwards in the threads
#ifdef debug_mapper
#pragma omp critical
            {
                if (debug) cerr << "can't go backward" << endl;
            }
#endif
            return -std::numeric_limits<double>::max();
        } else {
            //int max_length = (m1.length() + m2.length());
            int max_length = max(read1.sequence().size(), read2.sequence().size());
            // find the difference in m1.end and m2.begin
            // find the positional difference in the graph between m1.end and m2.begin
            int duplicate_coverage = mems_overlap_length(m1, m2);
            int64_t approx_dist = mem_min_distance(m1, m2); //graph_mixed_distance_estimate(m1_pos, m2_pos, m1.length() + m2.length());
#ifdef debug_mapper
#pragma omp critical
            {
                if (debug) cerr << "in fragment approx distance " << approx_dist << endl;
            }
#endif
            if (approx_dist > max_length) {
                // too far
                
#ifdef debug_mapper
#pragma omp critical
                {
                    if (debug) cerr << "too far apart on same node by approx_dist" << endl;
                }
#endif
                
                return -std::numeric_limits<double>::max();
            } else {
                // we may want to switch back to exact measurement, although the approximate metric is simpler and more reliable despite being less precise
                //int64_t distance = min(approx_dist, graph_distance(m1_pos, m2_pos, m1.length())); // enable for exact distance calculation
                int64_t distance = approx_dist;
#ifdef debug_mapper
#pragma omp critical
                {
                    if (debug) cerr << "---> true distance " << distance << endl;
                }
#endif
                if (distance >= max_length) {
#ifdef debug_mapper
#pragma omp critical
                    {
                        if (debug) cerr << "exactly too far apart" << endl;
                    }
#endif
                    return -std::numeric_limits<double>::max();
                }
                if (is_rev(m1_pos) != is_rev(m2_pos)) {
                    // disable inversions
                    // TODO: shouldn't we be using cached_fragment_orientation like in the two-node case???
                    
#ifdef debug_mapper
#pragma omp critical
                    {
                        if (debug) cerr << "no inversions allowed on the same node" << endl;
                    }
#endif
                    return -std::numeric_limits<double>::max();
                } else {
                    // accepted transition
                    double jump = abs((m2.begin - m1.begin) - distance);
                    if (jump) {
#ifdef debug_mapper
#pragma omp critical
                        {
                            if (debug) cerr << "accept distance with jump" << endl;
                        }
#endif
                        return (double) -duplicate_coverage * match - (gap_open + jump * gap_extension);
                    } else {
#ifdef debug_mapper
#pragma omp critical
                        {
                            if (debug) cerr << "accept distance without jump" << endl;
                        }
#endif
                        return (double) -duplicate_coverage * match;
                    }
                }
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
                              [&](pos_t n) -> map<string, vector<size_t> > {
                                  return xindex->position_in_paths(id(n), is_rev(n), offset(n));
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

#pragma omp critical
    {
        if (debug) {
            cerr << "### clusters before filtering:" << endl;
            show_clusters();
        }
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
        
#pragma omp critical
    {
        if (debug) {
            cerr << "### clusters after filtering:" << endl;
            show_paired_clusters();
        }
    }
    
    set<pair<string, string> > seen_alignments;
    for (auto& cluster_ptr : cluster_ptrs) {
        // break the cluster into two pieces
        auto& cluster1 = *cluster_ptr.first;
        auto& cluster2 = *cluster_ptr.second;
        alns.emplace_back();
        if ((to_drop1.count(&cluster1) || to_drop2.count(&cluster2)) && alns.size() >= min_multimaps) {
            continue;
        }
        auto& p = alns.back();
        if (cluster1.size()) {
            p.first = align_cluster(read1, cluster1, false);
        } else {
            p.first = read1;
            p.first.clear_score();
            p.first.clear_identity();
            p.first.clear_path();
        }
        if (cluster2.size()) {
            p.second = align_cluster(read2, cluster2, false);
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
            if (aln1.fragment_size() == 0 || aln2.fragment_size() == 0) {
                auto approx_frag_lengths = approx_pair_fragment_length(aln1, aln2);
                frag_stats.save_frag_lens_to_alns(aln1, aln2, approx_frag_lengths, pair_consistent(aln1, aln2, 0.0001));
            }
        }
        // sort the aligned pairs by score
        std::sort(aln_ptrs.begin(), aln_ptrs.end(),
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

    double max_possible_score = read1.sequence().size() * match + 2*full_length_bonus;

    // now add back in single-ended versions of everything
    vector<pair<Alignment, Alignment> > se_alns;
    vector<pair<vector<MaximalExactMatch>*, vector<MaximalExactMatch>*> > se_cluster_ptrs;
    for (auto& p : aln_ptrs) {
        auto& aln1 = p->first;
        auto& aln2 = p->second;
        // if both mates are aligned, add each single end into the mix
        if (aln1.score() && aln2.score()) {
            auto cluster_ptr = cluster_ptrs[aln_index[p]];
            // if these can be used for rescue
            if (aln1.score() > max_possible_score * pair_rescue_hang_threshold) {
                se_alns.emplace_back();
                auto& p = se_alns.back();
                p.first = aln1;
                p.second = read2;
                se_cluster_ptrs.push_back(make_pair(cluster_ptr.first, nullptr));
            }
            if (aln2.score() > max_possible_score * pair_rescue_hang_threshold) {
                se_alns.emplace_back();
                auto& q = se_alns.back();
                q.first = read1;
                q.second = aln2;
                se_cluster_ptrs.push_back(make_pair(nullptr, cluster_ptr.second));
            }
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
            if (++j > mate_rescues) break;
            auto& aln1 = p->first;
            auto& aln2 = p->second;
            int score1 = aln1.score();
            int score2 = aln2.score();
            pair<bool, bool> rescues = pair_rescue(aln1, aln2, match, full_length_bonus, false);
            rescued_aln[&aln1] = rescues.first;
            rescued_aln[&aln2] = rescues.second;
            rescued |= rescues.first || rescues.second;
        }
        if (rescued) {
            sort_and_dedup();
            show_alignments("rescue");
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

    // realign to generate traceback if needed
    int i = 0;
    bool rescored = false;
    bool rescued = false;
    
    auto traceback_alns = [&](void) {
        rescored = false;
        rescued = false;
        for (auto& p : aln_ptrs) {
            auto& aln1 = p->first;
            auto& aln2 = p->second;
            if (++i > max_multimaps) {
                // we don't need to realign because we are not emitting this multimap
                //aln1.clear_fragment();
                //aln2.clear_fragment();
            } else {
                // we give exactly the same cluster to the alignments to get the traceback
                assert(aln_index.find(p) != aln_index.end());
                auto cluster_ptr = cluster_ptrs[aln_index[p]];
                auto s1 = aln1.score();
                auto s2 = aln2.score();
                if (rescued_aln[&aln1]) {
                    assert(!alignment_to_length(aln1) && aln1.path().mapping(0).has_position());
                    // realign based on alignment end position
                    aln1 = realign_from_start_position(aln1, aln1.sequence().size()/1.61803, 0);
                } else if (cluster_ptr.first != nullptr && cluster_ptr.first->size()) {
                    if (!alignment_to_length(aln1)) { // traceback needs to be generated
                        aln1 = align_cluster(read1, *cluster_ptr.first, true);
                    }
                }
                if (rescued_aln[&aln2]) {
                    assert(!alignment_to_length(aln2) && aln2.path().mapping(0).has_position());
                    // realign based on alignment end position
                    aln2 = realign_from_start_position(aln2, aln2.sequence().size()/1.61803, 0);
                } else if (cluster_ptr.second != nullptr && cluster_ptr.second->size()) {
                    if (!alignment_to_length(aln2)) { // traceback needs to be generated
                        aln2 = align_cluster(read2, *cluster_ptr.second, true);
                    }
                }
                //assert(aln1.score() >= s1);
                //assert(aln2.score() >= s2);
                if (aln1.score() > s1 || aln2.score() > s2) rescored = true;
                // we can reassign based on paths to get a more accurate fragment estimate
                aln1.clear_fragment();
                aln2.clear_fragment();
                auto approx_frag_lengths = approx_pair_fragment_length(aln1, aln2);
                frag_stats.save_frag_lens_to_alns(aln1, aln2, approx_frag_lengths, pair_consistent(aln1, aln2, 0.0001));
            }
        }
    };
    
    // traceback alignments
    traceback_alns();
    if (rescored) {
        // sync if we need to
        sort_and_dedup();
        traceback_alns();
    }
    show_alignments("end");
    
    int read1_max_score = 0;
    int read2_max_score = 0;
    // we must be considering more than one pair to do the paired-end mq calculation
    int possible_pairs = 0;

    // build up the results
    i = 0;
    for (auto& p : aln_ptrs) {
        read1_max_score = max(p->first.score(), read1_max_score);
        read2_max_score = max(p->second.score(), read2_max_score);
        results.first.push_back(p->first);
        results.second.push_back(p->second);
        possible_pairs += p->first.score() > 0 && p->second.score() > 0;
    }
    bool max_first = results.first.size() && (read1_max_score == results.first.front().score() && read2_max_score == results.second.front().score());

    double mem_read_ratio1 = min(1.0, (double)total_mem_length1 / (double)read1.sequence().size());
    double mem_read_ratio2 = min(1.0, (double)total_mem_length2 / (double)read2.sequence().size());
    double mqmax1 = max_mapping_quality;
    double mqmax2 = max_mapping_quality;
    // calculate paired end quality if the model assumptions are not obviously violated
    if (results.first.size() && results.second.size()
        && (fraction_filtered1 < 0.1 && fraction_filtered2 < 0.1 && maybe_mq1 > 1 && maybe_mq2 > 1 && max_first && (mem_read_ratio1 > 0.5 || mem_read_ratio2 > 0.5) || possible_pairs > 1) // may help in human context
        && pair_consistent(results.first.front(), results.second.front(), 0.0001)) {
        compute_mapping_qualities(results, cluster_mq, mq_cap1, mq_cap2, mqmax1, mqmax2);
    } else {
        // through filtering of candidate hits we've ended up with only one possible pair
        if (results.first.size() < 2 && results.second.size() < 2) {
            mqmax1 = maybe_mq1;
            mqmax2 = maybe_mq2;
        }
        mqmax1 = mem_read_ratio1 > 0.5 ? mqmax1 : mem_read_ratio1 * mqmax1;
        mqmax2 = mem_read_ratio2 > 0.5 ? mqmax2 : mem_read_ratio2 * mqmax2;
        // compute mq independently
        compute_mapping_qualities(results.first, cluster_mq, mq_cap1, mqmax1);
        compute_mapping_qualities(results.second, cluster_mq, mq_cap2, mqmax2);
    }

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
        << ':' << frag_stats.cached_fragment_orientation 
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
        if (retrying) break;
        auto& aln1 = results.first.at(i);
        auto& aln2 = results.second.at(i);
        for (int j = 0; j < aln1.fragment_size(); ++j) {
            int length = aln1.fragment(j).length();
            // if we have a perfect mapping, and we're under our hard fragment length cutoff
            // push the result into our deque of fragment lengths
            if (results.first.size() == 1
                && results.second.size() == 1
                && results.first.front().identity() > frag_stats.perfect_pair_identity_threshold
                && results.second.front().identity() > frag_stats.perfect_pair_identity_threshold
                && (frag_stats.fragment_size && abs(length) < frag_stats.fragment_size
                    || !frag_stats.fragment_size && abs(length) < frag_stats.fragment_max)) { // hard cutoff
                //cerr << "aln\tperfect alignments" << endl;
                frag_stats.record_fragment_configuration(length, aln1, aln2);
            } else if (!frag_stats.fragment_size) {
                imperfect_pair = true;
                break;
            }
        }
    }

    if (!retrying && imperfect_pair && frag_stats.fragment_max) {
        imperfect_pairs_to_retry.push_back(make_pair(read1, read2));
        results.first.clear();
        results.second.clear();
        // we signal the fact that this isn't a perfect pair, so we don't write it out externally?
        queued_resolve_later = true;
    }

    if(results.first.empty()) {
        results.first.push_back(read1);
        auto& aln = results.first.back();
        aln.clear_path();
        aln.clear_score();
        aln.clear_identity();
    }
    if(results.second.empty()) {
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

    // if we have references, annotate the alignments with their reference positions
    annotate_with_mean_path_positions(results.first);
    annotate_with_mean_path_positions(results.second);

    return results;

}

void Mapper::annotate_with_mean_path_positions(vector<Alignment>& alns) {
    for (auto& aln : alns) {
        for (auto& ref : alignment_mean_path_positions(aln)) {
            Position* refpos = aln.add_refpos();
            refpos->set_name(ref.first);
            refpos->set_offset(round(ref.second));
        }
    }
}
double Mapper::compute_cluster_mapping_quality(const vector<vector<MaximalExactMatch> >& clusters,
                                               int read_length) {
    if (clusters.size() == 0) {
        return 0;
    }
    if (clusters.size() == 1) {
        return { (double)max_cluster_mapping_quality };
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

double
Mapper::average_node_length(void) {
    return (double) xindex->seq_length / (double) xindex->node_count;
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
    map<const vector<MaximalExactMatch>*, int> cluster_cov;
    for (auto& cluster : clusters) {
        cluster_cov[&cluster] = cluster_coverage(cluster);
    }
    for (int i = 0; i < clusters.size(); ++i) {
        // establish overlaps with longer clusters for all clusters
        auto& this_cluster = clusters[i];
        int t = cluster_cov[&this_cluster];
        int b = -1;
        int l = t;
        for (int j = i; j >= 0; --j) {
            if (j == i) continue;
            // are we overlapping?
            auto& other_cluster = clusters[j];
            //if (to_drop.count(&other_cluster)) continue;
            if (clusters_overlap_in_graph(this_cluster, other_cluster)) {
                to_drop.insert(&this_cluster);
                break;
            }
            if (clusters_overlap_in_read(this_cluster, other_cluster)) {
                int c = cluster_cov[&other_cluster];
                if (c > l) {
                    l = c;
                    b = j;
                }
                if (b >= 0 && (float) t / (float) l < drop_chain) {
                    to_drop.insert(&this_cluster);
                    break;
                }
            }
        }
        if (b >= 0 && (float) t / (float) l < drop_chain) {
            to_drop.insert(&this_cluster);
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
                        int additional_multimaps) {

//#ifdef debug_mapper
#pragma omp critical
    {
        if (debug) cerr << "mems for read " << mems_to_json(mems) << endl;
    }
//#endif
    
    auto aligner = get_aligner(!aln.quality().empty());
    int8_t match = aligner->match;
    int8_t gap_extension = aligner->gap_extension;
    int8_t gap_open = aligner->gap_open;

    int total_multimaps = max(max_multimaps, additional_multimaps);
    double mq_cap = max_mapping_quality;

    int total_mem_length = 0;
    for (auto& mem : mems) total_mem_length += mem.length() * mem.nodes.size();

    // Estimate the maximum mapping quality we can get if the alignments based on the good MEMs are the best ones.
    int mem_max_length = 0;
    for (auto& mem : mems) if (mem.primary && mem.match_count) mem_max_length = max(mem_max_length, (int)mem.length());
    double maybe_mq = estimate_max_possible_mapping_quality(aln.sequence().size(),
                                                            aln.sequence().size()/max(1.0, (double)mem_max_length),
                                                            aln.sequence().size()/longest_lcp);
    // use the estimated mapping quality to avoid hard work when the outcome will be noninformative
    if (maybe_mq < maybe_mq_threshold) {
        mq_cap = maybe_mq;
    }
    // TODO: why should we limit our number of MEM chains to examine to max_multimaps / max estimated mapping quality?
    total_multimaps = max(min_multimaps, (int)round(total_multimaps/maybe_mq));
    if (debug) cerr << "maybe_mq " << aln.name() << " max estimate: " << maybe_mq << " estimated multimap limit: " << total_multimaps << " max mem length: " << mem_max_length << " min mem length: " << min_mem_length << " longest LCP: " << longest_lcp << endl;

    double avg_node_len = average_node_length();
    // go through the ordered single-hit MEMs
    // build the clustering model
    // find the alignments that are the best-scoring walks through it
    auto transition_weight = [&](const MaximalExactMatch& m1, const MaximalExactMatch& m2) {

        int duplicate_coverage = mems_overlap_length(m1, m2);
        pos_t m1_pos = make_pos_t(m1.nodes.front());
        pos_t m2_pos = make_pos_t(m2.nodes.front());
        int64_t max_length = aln.sequence().size();
        int64_t approx_dist = mem_min_distance(m1, m2);
        
#ifdef debug_mapper
#pragma omp critical
        {
            if (debug) cerr << "mems " << &m1 << ":" << m1 << " -> " << &m2 << ":" << m2 << " approx_dist " << approx_dist << " duplicate_coverage " << duplicate_coverage << endl;
        }
#endif
        if (approx_dist > max_length) {
            // too far
            return -std::numeric_limits<double>::max();
        } else {
            if (is_rev(m1_pos) != is_rev(m2_pos)) {
                // disable inversions
                return -std::numeric_limits<double>::max();
            } else {
                //int64_t distance = min(approx_dist, graph_distance(m1_pos, m2_pos, m1.length())); // enable for exact distance calculation
                int64_t distance = approx_dist;
                // accepted transition
                double jump = abs((m2.begin - m1.begin) - distance);
                if (jump) {
                    return (double) -duplicate_coverage * match - (gap_open + jump * gap_extension);
                } else {
                    return (double) -duplicate_coverage * match;
                }
            }
        }
    };

    // establish the chains
    vector<vector<MaximalExactMatch> > clusters;
    if (total_multimaps) {
        MEMChainModel chainer({ aln.sequence().size() }, { mems },
                              [&](pos_t n) {
                                  return approx_position(n);
                              },
                              [&](pos_t n) -> map<string, vector<size_t> > {
                                  return xindex->position_in_paths(id(n), is_rev(n), offset(n));
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
                    /*
                    for (auto& ref : node_positions_in_paths(gcsa::Node::encode(id, 0, is_rev))) {
                        auto& name = ref.first;
                        for (auto pos : ref.second) {
                            //cerr << name << (is_rev?"-":"+") << pos + offset;
                            cerr << "|" << id << (is_rev ? "-" : "+") << ":" << offset << ",";
                        }
                    }
                    */
                }
                cerr << mem.sequence() << " ";
            }
            cerr << endl;
        }
    };

//#ifdef debug_mapper
#pragma omp critical
    {
        if (debug) {
            cerr << "### clusters:" << endl;
            show_clusters();
        }
    }
//#endif

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
    for (auto& cluster : clusters) {
        if (alns.size() >= total_multimaps) { break; }
        // skip if we've filtered the cluster
        if (to_drop.count(&cluster) && alns.size() >= min_multimaps) {
            alns.push_back(aln);
            used_clusters.push_back(&cluster);
            continue;
        }
        // skip if we've got enough multimaps to get MQ and we're under the min cluster length
        if (min_cluster_length && cluster_coverage(cluster) < min_cluster_length && alns.size() > 1) {
            alns.emplace_back(aln);
            used_clusters.push_back(&cluster);
            continue;
        }
        Alignment candidate = align_cluster(aln, cluster, false);
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

#pragma omp critical
    if (debug) {
        cerr << "alignments" << endl;
        for (auto& aln : alns) {
            cerr << aln.score();
            if (aln.score()) cerr << " pos1 " << aln.path().mapping(0).position().node_id() << " ";
            cerr << endl;
        }
    }

    vector<Alignment*> aln_ptrs;
    map<Alignment*, int> aln_index;
    int idx = 0;
    for (auto& aln : alns) {
        aln_ptrs.push_back(&aln);
        aln_index[&aln] = idx++;
    }
    // sort alignments by score
    std::sort(aln_ptrs.begin(), aln_ptrs.end(),
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
            // only realign if we haven't yet
            if (!alignment_to_length(*alnp)) {
                auto& cluster = *used_clusters.at(aln_index.at(alnp));
                Alignment candidate = align_cluster(aln, cluster, true);
                best_alns.push_back(candidate);
            }
        }
        for ( ; i < aln_ptrs.size(); ++i) {
            best_alns.push_back(*aln_ptrs[i]);
        }
        alns = score_sort_and_deduplicate_alignments(best_alns, aln);
    }
    // compute the mapping quality
    compute_mapping_qualities(alns, cluster_mq, mq_cap, max_mapping_quality);

    // final filter step
    filter_and_process_multimaps(alns, keep_multimaps);

    // if we didn't get anything, return an unaligned version of our input
    if (alns.empty()) {
        alns.push_back(aln);
        auto& unaligned = alns.back();
        unaligned.clear_path();
        unaligned.clear_score();
        unaligned.clear_identity();
    }

    return alns;
}

Alignment Mapper::align_maybe_flip(const Alignment& base, Graph& graph, bool flip, bool traceback, bool banded_global) {
    Alignment aln = base;
    map<id_t, int64_t> node_length;
    if (flip) {
        for (auto& node : graph.node()) {
            node_length[node.id()] = node.sequence().size();
        }
        aln.set_sequence(reverse_complement(base.sequence()));
        if (!base.quality().empty()) {
            aln.set_quality(base.quality());
            reverse(aln.mutable_quality()->begin(),
                    aln.mutable_quality()->end());
        }
    } else {
        aln.set_sequence(base.sequence());
        if (!base.quality().empty()) {
            aln.set_quality(base.quality());
        }
    }
    bool pinned_alignment = false;
    bool pinned_reverse = false;

    aln = align_to_graph(aln,
                         graph,
                         max_query_graph_ratio,
                         traceback,
                         pinned_alignment,
                         pinned_reverse,
                         banded_global);

    if (strip_bonuses && !banded_global && traceback) {
        // We want to remove the bonuses
        aln.set_score(get_aligner()->remove_bonuses(aln));
    }
    if (flip) {
        aln = reverse_complement_alignment(
            aln,
            (function<int64_t(int64_t)>) ([&](int64_t id) {
                    return node_length[id];
                }));
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

Alignment Mapper::align_cluster(const Alignment& aln, const vector<MaximalExactMatch>& mems, bool traceback) {
    // poll the mems to see if we should flip
    int count_fwd = 0, count_rev = 0;
    for (auto& mem : mems) {
        bool is_rev = gcsa::Node::rc(mem.nodes.front());
        if (is_rev) {
            ++count_rev;
        } else {
            ++count_fwd;
        }
    }
    // get the graph with cluster.hpp's cluster_subgraph
    Graph graph = cluster_subgraph(*xindex, aln, mems);
    // and test each direction for which we have MEM hits
    Alignment aln_fwd;
    Alignment aln_rev;
    if (count_fwd) {
        aln_fwd = align_maybe_flip(aln, graph, false, traceback);
    }
    if (count_rev) {
        aln_rev = align_maybe_flip(aln, graph, true, traceback);
    }
    // TODO check if we have soft clipping on the end of the graph and if so try to expand the context
    if (aln_fwd.score() + aln_rev.score() == 0) {
        // abject failure, nothing aligned with score > 0
        Alignment result = aln;
        result.clear_path();
        result.clear_score();
        return result;
    } else if (aln_rev.score() > aln_fwd.score()) {
        // reverse won
        return aln_rev;
    } else {
        // forward won
        return aln_fwd;
    }
}

VG Mapper::cluster_subgraph_strict(const Alignment& aln, const vector<MaximalExactMatch>& mems) {
#ifdef debug_mapper
#pragma omp critical
    {
        if (debug) {
            cerr << "Getting a cluster graph for " << mems.size() << " MEMs" << endl;
        }
    }
#endif

    // As in the multipath aligner, we work out how far we can get from a MEM
    // with gaps and use that for how much graph to grab.
    vector<pos_t> positions;
    vector<size_t> forward_max_dist;
    vector<size_t> backward_max_dist;
    
    positions.reserve(mems.size());
    forward_max_dist.reserve(mems.size());
    backward_max_dist.reserve(mems.size());
    
    // What aligner are we using?
    BaseAligner* aligner = get_aligner();
    
    for (const auto& mem : mems) {
        // get the start position of the MEM
        assert(!mem.nodes.empty());
        positions.push_back(make_pos_t(mem.nodes.front()));
        
        // search far enough away to get any hit detectable without soft clipping
        forward_max_dist.push_back(aligner->longest_detectable_gap(aln, mem.end)
                                   + (aln.sequence().end() - mem.begin));
        backward_max_dist.push_back(aligner->longest_detectable_gap(aln, mem.begin)
                                    + (mem.begin - aln.sequence().begin()));
    }
    
    
    // extract the protobuf Graph
    Graph proto_graph;
    algorithms::extract_containing_graph(xindex, proto_graph, positions, forward_max_dist, backward_max_dist);
                                         
    // Wrap it in a vg
    VG graph;
    graph.extend(proto_graph);
    
    graph.remove_orphan_edges();
    
#ifdef debug_mapper
#pragma omp critical
    {
        if (debug) {
            cerr << "\tFound " << graph.node_count() << " nodes " << graph.min_node_id() << " - " << graph.max_node_id()
                << " and " << graph.edge_count() << " edges" << endl;
        }
    }
#endif
    return graph;
}

VG Mapper::alignment_subgraph(const Alignment& aln, int context_size) {
    set<id_t> nodes;
    auto& path = aln.path();
    for (int i = 0; i < path.mapping_size(); ++i) {
        nodes.insert(path.mapping(i).position().node_id());
    }
    VG graph;
    for (auto& node : nodes) {
        *graph.graph.add_node() = xindex->node(node);
    }
    xindex->expand_context(graph.graph, max(1, context_size), false); // get connected edges
    graph.rebuild_indexes();
    return graph;
}

// estimate the fragment length as the difference in mean positions of both alignments
map<string, int64_t> Mapper::approx_pair_fragment_length(const Alignment& aln1, const Alignment& aln2) {
    map<string, int64_t> lengths;
    auto pos1 = alignment_mean_path_positions(aln1);
    auto pos2 = alignment_mean_path_positions(aln2);
    for (auto& p : pos1) {
        auto x = pos2.find(p.first);
        if (x != pos2.end()) {
            lengths[p.first] = x->second - p.second;
        }
    }
    return lengths;
}

string FragmentLengthStatistics::fragment_model_str(void) {
    stringstream s;
    s << fragment_size << ":"
      << cached_fragment_length_mean << ":"
      << cached_fragment_length_stdev << ":"
      << cached_fragment_orientation << ":"
      << cached_fragment_direction;
    return s.str();
}

int64_t Mapper::first_approx_pair_fragment_length(const Alignment& aln1, const Alignment& aln2) {
    auto pos1 = alignment_mean_path_positions(aln1);
    auto pos2 = alignment_mean_path_positions(aln2);
    for (auto& p : pos1) {
        auto x = pos2.find(p.first);
        if (x != pos2.end()) {
            return x->second - p.second;
        }
    }
    return -1;
}

void FragmentLengthStatistics::save_frag_lens_to_alns(Alignment& aln1, Alignment& aln2,
    const map<string, int64_t>& approx_frag_lengths, bool is_consistent) {
    
    for (auto& j : approx_frag_lengths) {
        Path fragment;
        fragment.set_name(j.first);
        int length = j.second;
        fragment.set_length(length);
        *aln1.add_fragment() = fragment;
        *aln2.add_fragment() = fragment;
        if (fragment_size && is_consistent) {
            double pval = fragment_length_pval(abs(length));
            double score = pval > 0.01 ? 10 + pval : 0;
            aln1.set_fragment_score(score);
            aln2.set_fragment_score(score);
        } else if (length < fragment_max) {
            aln1.set_fragment_score(0);
            aln2.set_fragment_score(0);
        } else {
            aln1.set_fragment_score(0);
            aln2.set_fragment_score(0);
        }
    }
}

void FragmentLengthStatistics::record_fragment_configuration(int length, const Alignment& aln1, const Alignment& aln2) {
    if (fixed_fragment_model) return;
    // record the relative orientations
    assert(aln1.path().mapping(0).has_position() && aln2.path().mapping(0).has_position());
    bool aln1_is_rev = aln1.path().mapping(0).position().is_reverse();
    bool aln2_is_rev = aln2.path().mapping(0).position().is_reverse();
    bool same_orientation = aln1_is_rev == aln2_is_rev;
    fragment_orientations.push_front(same_orientation);
    if (fragment_orientations.size() > fragment_length_cache_size) {
        fragment_orientations.pop_back();
    }
    // assuming a dag-like graph
    // which direction do we go relative to the orientation of our first mate to find the second?
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
        cached_fragment_orientation = fragment_orientation();
        cached_fragment_direction = fragment_direction();
        // set our fragment size cap to the cached mean + 10x the standard deviation
        fragment_size = cached_fragment_length_mean + fragment_sigma * cached_fragment_length_stdev;
        since_last_fragment_length_estimate = 1;
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
    return 1 - phi(-x,x);
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

set<MaximalExactMatch*> Mapper::resolve_paired_mems(vector<MaximalExactMatch>& mems1,
                                                    vector<MaximalExactMatch>& mems2) {
    // find the MEMs that are within estimated_fragment_size of each other

    set<MaximalExactMatch*> pairable;

    // do a wide clustering and then do all pairs within each cluster
    // we will use these to determine the alignment strand
    //map<id_t, StrandCounts> node_strands;
    // records a mapping of id->MEMs, for cluster ranking
    map<id_t, vector<MaximalExactMatch*> > id_to_mems;
    // for clustering
    set<id_t> ids1, ids2;
    vector<id_t> ids;

    // run through the mems
    for (auto& mem : mems1) {
        for (auto& node : mem.nodes) {
            id_t id = gcsa::Node::id(node);
            id_to_mems[id].push_back(&mem);
            ids1.insert(id);
            ids.push_back(id);
        }
    }
    for (auto& mem : mems2) {
        for (auto& node : mem.nodes) {
            id_t id = gcsa::Node::id(node);
            id_to_mems[id].push_back(&mem);
            ids2.insert(id);
            ids.push_back(id);
        }
    }
    // remove duplicates
    //std::sort(ids.begin(), ids.end());
    //ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

    // get each hit's path-relative position
    map<string, map<int, vector<id_t> > > node_positions;
    for (auto& id : ids) {
        for (auto& ref : node_positions_in_paths(gcsa::Node::encode(id, 0))) {
            auto& name = ref.first;
            for (auto pos : ref.second) {
                node_positions[name][pos].push_back(id);
            }
        }
    }

    vector<vector<id_t> > clusters;
    for (auto& g : node_positions) {
        //if (g.second.empty()) continue; // should be impossible
        //cerr << g.first << endl;
        clusters.emplace_back();
        int prev = -1;
        for (auto& x : g.second) {
            auto cluster = &clusters.back();
            //auto& prev = clusters.back().back();
            auto curr = x.first;
            if(debug) {
                cerr << "p/c " << prev << " " << curr << endl;
            }
            if (prev != -1) {
                if (curr - prev <= frag_stats.fragment_size) {
                    // in cluster
#ifdef debug_mapper
#pragma omp critical
                    {
                        if (debug) {
                            cerr << "in cluster" << endl;
                        }
                    }
#endif
                } else {
                    // It's a new cluster
                    clusters.emplace_back();
                    cluster = &clusters.back();
                }
            }
            //cerr << " " << x.first << endl;
            for (auto& y : x.second) {
                //cerr << "  " << y << endl;
                cluster->push_back(y);
            }
            prev = curr;
        }
    }

    for (auto& cluster : clusters) {
        // for each pair of ids in the cluster
        // which are not from the same read
        // estimate the distance between them
        // we're roughly in the expected range
        bool has_first = false;
        bool has_second = false;
        for (auto& id : cluster) {
            has_first |= ids1.count(id);
            has_second |= ids2.count(id);
        }
        if (!has_first || !has_second) continue;
        for (auto& id : cluster) {
            for (auto& memptr : id_to_mems[id]) {
                pairable.insert(memptr);
            }
        }
    }

    return pairable;
}

// We need a function to get the lengths of nodes, in case we need to
// reverse an Alignment, including all its Mappings and Positions.
int64_t Mapper::get_node_length(int64_t node_id) {
    // Grab the node sequence only from the XG index and get its size.
    // Make sure to use the cache
    return xg_node_length(node_id, xindex);
}

bool Mapper::check_alignment(const Alignment& aln) {
    // use the graph to extract the sequence
    // assert that this == the alignment
    if (aln.path().mapping_size()) {
        // get the graph corresponding to the alignment path
        Graph sub;
        for (int i = 0; i < aln.path().mapping_size(); ++ i) {
            auto& m = aln.path().mapping(i);
            if (m.has_position() && m.position().node_id()) {
                auto id = aln.path().mapping(i).position().node_id();
                // XXXXXX this is single-threaded!
                xindex->neighborhood(id, 2, sub);
            }
        }
        VG g; g.extend(sub);
        auto seq = g.path_string(aln.path());
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
            write_alignment_to_file(aln, "fail-" + hash_alignment(aln) + ".gam");
            // save graph, bigger fragment
            xindex->expand_context(sub, 5, true);
            VG gn; gn.extend(sub);
            gn.serialize_to_file("fail-" + gn.hash() + ".vg");
            return false;
        }
    }
    return true;
}

vector<Alignment> Mapper::make_bands(const Alignment& read, int band_width, vector<pair<int, int>>& to_strip) {
    if (band_width % 4) {
        band_width -= band_width % 4; band_width += 4;
    }
    int div = 2;
    while (read.sequence().size()/div > band_width) {
        ++div;
    }
    int segment_size = read.sequence().size()/div;
    // use segment sizes divisible by 4, as this simplifies math
    // round up as well
    // we'll divide the overlap by 2 and 2 and again when stripping from the start
    // and end of sub-reads
    if (segment_size % 4) {
        segment_size -= segment_size % 4; segment_size += 4;
    }
#ifdef debug_mapper
    if (debug) {
        cerr << "Segment size be " << segment_size << "/" << read.sequence().size() << endl;
    }
#endif
    // and overlap them too
    //int to_align = div * 2 - 1; // number of alignments we'll do
    //vector<pair<size_t, size_t>> to_strip; to_strip.resize(to_align);
    //vector<Alignment> bands; bands.resize(to_align);

    int remainder = (int)read.sequence().size() % segment_size;
    //cerr << "remainder " << remainder << endl;
    if (remainder % 2) {
        remainder -= remainder % 2; remainder += 2; // make divisible by 2
    }
    //cerr << "remainder adj " << remainder << endl;

    vector<int> start_positions;

    // record the start positions
    int offset = 0;
    while (offset+segment_size < read.sequence().size()) {
        if (offset == 0) {
            start_positions.push_back(0);
            if (remainder) {
                offset += remainder/2;
            } else {
                offset += segment_size/2;
            }
        } else {
            start_positions.push_back(offset);
            offset += segment_size/2;
        }
    }
    // add in the last alignment
    start_positions.push_back(read.sequence().size()-segment_size);
    // set up the structures to hold onto the banded alignments
    int to_align = start_positions.size();
    to_strip.resize(to_align);
    vector<Alignment> bands; bands.resize(to_align);
    int i = 0;
    int q = segment_size/4;
    for (auto& p : start_positions) {
        // so that we tend to obtain the path component derived from the middle of alignments
        // with the exception of the first and last alignments
        // we will remove half the overlap length from each end
        // the overlap length is 1/2 the bandwidth, so we remove 1/4 length from each end
        if (&p == &start_positions.front()) {
            to_strip[i].second = segment_size - (start_positions[i+1] + q);
        } else if (&p == &start_positions.back()) {
            to_strip[i].first = start_positions[i-1]+segment_size-q - p;
            // if we only have two bands, handle the potential non-divisibility by 2
            if (start_positions.size() == 2) {
                to_strip[i].first -= (int)read.sequence().size() % 2;
            }
        } else {
            to_strip[i].first = q;
            to_strip[i].second = q;
        }
        //cerr << "position: " << p << " strip " << to_strip[i].first << " " << to_strip[i].second << endl;
        auto& aln = bands[i];
        aln.set_sequence(read.sequence().substr(p, segment_size));
        if (!read.quality().empty()) aln.set_quality(read.quality().substr(p, segment_size));
        ++i;
    }
    return bands;
}

vector<Alignment> Mapper::align_banded(const Alignment& read, int kmer_size, int stride, int max_mem_length, int band_width) {

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

    // scan across the read choosing bands
    // these bands are hard coded to overlap by 50%
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
    vector<Alignment> bands = make_bands(read, band_width, to_strip);
    vector<vector<Alignment>> multi_alns;
    multi_alns.resize(bands.size());

    auto do_band = [&](int i) {
        //cerr << "aligning band " << i << endl;
        vector<Alignment>& malns = multi_alns[i];
        double cluster_mq = 0;
        malns = align_multi_internal(true, bands[i], kmer_size, stride, max_mem_length, bands[i].sequence().size(), cluster_mq, band_multimaps, extra_multimaps, nullptr);
        // always include an unaligned mapping
        malns.push_back(bands[i]);
        for (vector<Alignment>::iterator a = malns.begin(); a != malns.end(); ++a) {
            Alignment& aln = *a;
            aln = simplify(aln);
            bool above_threshold = false;
            if (aln.score() > 0) {
                // strip overlaps and re-score the part of the alignment we keep
                aln = strip_from_start(aln, to_strip[i].first);
                aln = strip_from_end(aln, to_strip[i].second);
                aln.set_identity(identity(aln.path()));
                above_threshold = aln.identity() >= min_identity && aln.mapping_quality() >= min_banded_mq;
            }
            if (!above_threshold) {
                // treat as unmapped
                aln = bands[i];
                aln = strip_from_start(aln, to_strip[i].first);
                aln = strip_from_end(aln, to_strip[i].second);
            }
        }
    };

    if (alignment_threads > 1) {
#pragma omp parallel for
        for (int i = 0; i < bands.size(); ++i) {
            do_band(i);
        }
    } else {
        for (int i = 0; i < bands.size(); ++i) {
            do_band(i);
        }
    }

    // cost function
    auto transition_weight = [&](const Alignment& aln1, const Alignment& aln2,
                                 const map<string, double>& pos1, const map<string, double>& pos2) {
        // scoring scheme for unaligned reads
        if (!aln1.has_path() || !aln2.has_path()) {
            return -(double)(10*(aln1.sequence().size() + aln2.sequence().size()) * gap_extension + gap_open);
        }
        auto aln1_end = make_pos_t(path_end(aln1.path()));
        auto aln2_begin = make_pos_t(path_start(aln2.path()));
        //auto dist = graph_mixed_distance_estimate(aln1_end, aln2_begin, 32);
        int64_t graph_dist = graph_distance(aln1_end, aln2_begin, 32);
        int64_t dist = std::numeric_limits<int64_t>::max();
        dist = min(graph_dist, dist);
        for (auto& p : pos1) {
            auto f = pos2.find(p.first);
            if (f != pos2.end()) {
                dist = min((int64_t)round(abs(f->second - p.second)), dist);
            }
        }
        if (debug) cerr << "dist " << dist << endl;
        return -((double)gap_open + (double)dist * (double)gap_extension);
    };

    AlignmentChainModel chainer(multi_alns, this, transition_weight, 1);
    if (debug) chainer.display(cerr);
    vector<Alignment> alignments = chainer.traceback(read, max_multimaps, false, debug);
    for (auto& aln : alignments) {
        // patch the alignment to deal with short unaligned regions
        aln = patch_alignment(aln, band_width);
    }
    // sort the alignments by score
    std::sort(alignments.begin(), alignments.end(), [](const Alignment& aln1, const Alignment& aln2) { return aln1.score() > aln2.score(); });
    if (alignments.size() == 1) {
        alignments.front().set_mapping_quality(max_mapping_quality);
    } else {
        compute_mapping_qualities(alignments, 0, max_mapping_quality, max_mapping_quality);
    }
    filter_and_process_multimaps(alignments, max_multimaps);
    //cerr << "got alignment " << pb2json(alignments.front()) << endl;
    return alignments;
}

bool Mapper::adjacent_positions(const Position& pos1, const Position& pos2) {
    // are they the same id, with offset differing by 1?
    if (pos1.node_id() == pos2.node_id()
        && pos1.offset() == pos2.offset()-1) {
        return true;
    }
    // otherwise, we're going to need to check via the index
    VG graph;
    // pick up a graph that's just the neighborhood of the start and end positions
    int64_t id1 = pos1.node_id();
    int64_t id2 = pos2.node_id();
    if(xindex) {
        // Grab the node sequence only from the XG index and get its size.
        xindex->get_id_range(id1, id1, graph.graph);
        xindex->get_id_range(id2, id2, graph.graph);
        xindex->expand_context(graph.graph, 1, false);
        graph.rebuild_indexes();
    } else {
        throw runtime_error("No index to get nodes from.");
    }
    // now look in the graph to figure out if we are adjacent
    return graph.adjacent(pos1, pos2);
}

void Mapper::compute_mapping_qualities(vector<Alignment>& alns, double cluster_mq, double mq_estimate, double mq_cap) {
    if (alns.empty()) return;
    double max_mq = min(mq_cap, (double)max_mapping_quality);
    BaseAligner* aligner = get_aligner();
    int sub_overlaps = 0; //sub_overlaps_of_first_aln(alns, mq_overlap);
    switch (mapping_quality_method) {
        case Approx:
            aligner->compute_mapping_quality(alns, max_mq, true, cluster_mq, use_cluster_mq, sub_overlaps, mq_estimate, identity_weight);
            break;
        case Exact:
            aligner->compute_mapping_quality(alns, max_mq, false, cluster_mq, use_cluster_mq, sub_overlaps, mq_estimate, identity_weight);
            break;
        default: // None
            break;
    }
}
    
void Mapper::compute_mapping_qualities(pair<vector<Alignment>, vector<Alignment>>& pair_alns, double cluster_mq, double mq_estimate1, double mq_estimate2, double mq_cap1, double mq_cap2) {
    if (pair_alns.first.empty() || pair_alns.second.empty()) return;
    double max_mq1 = min(mq_cap1, (double)max_mapping_quality);
    double max_mq2 = min(mq_cap2, (double)max_mapping_quality);
    BaseAligner* aligner = get_aligner();
    int sub_overlaps1 = 0; //sub_overlaps_of_first_aln(pair_alns.first, mq_overlap);
    int sub_overlaps2 = 0; //sub_overlaps_of_first_aln(pair_alns.second, mq_overlap);
    vector<double> frag_weights;
    for (int i = 0; i < pair_alns.first.size(); ++i) {
        auto& aln1 = pair_alns.first[i];
        frag_weights.push_back(aln1.fragment_score());
    }
    switch (mapping_quality_method) {
        case Approx:
            aligner->compute_paired_mapping_quality(pair_alns, frag_weights, max_mq1, max_mq2, true, cluster_mq, use_cluster_mq, sub_overlaps1, sub_overlaps2, mq_estimate1, mq_estimate2, identity_weight);
            break;
        case Exact:
            aligner->compute_paired_mapping_quality(pair_alns, frag_weights, max_mq1, max_mq2, false, cluster_mq, use_cluster_mq, sub_overlaps1, sub_overlaps2, mq_estimate1, mq_estimate2, identity_weight);
            break;
        default: // None
            break;
    }
}

double Mapper::estimate_max_possible_mapping_quality(int length, double min_diffs, double next_min_diffs) {
    return get_aligner()->estimate_max_possible_mapping_quality(length, min_diffs, next_min_diffs);
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
    
    map<int, set<Alignment*> > alignment_by_score;
    for (auto& ta : all_alns) {
        Alignment* aln = &ta;
        alignment_by_score[aln->score()].insert(aln);
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
    
vector<Alignment> Mapper::align_multi(const Alignment& aln, int kmer_size, int stride, int max_mem_length, int band_width) {
    double cluster_mq = 0;
    Alignment clean_aln;
    clean_aln.set_name(aln.name());
    clean_aln.set_sequence(aln.sequence());
    clean_aln.set_quality(aln.quality());
    return align_multi_internal(true, clean_aln, kmer_size, stride, max_mem_length, band_width, cluster_mq, max_multimaps, extra_multimaps, nullptr);
}
    
vector<Alignment> Mapper::align_multi_internal(bool compute_unpaired_quality,
                                               const Alignment& aln,
                                               int kmer_size, int stride,
                                               int max_mem_length,
                                               int band_width,
                                               double& cluster_mq,
                                               int keep_multimaps,
                                               int additional_multimaps,
                                               vector<MaximalExactMatch>* restricted_mems) {
    
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
        return vector<Alignment>{align_banded(aln, kmer_size, stride, max_mem_length, band_width)};
    }
    
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
        alignments = align_mem_multi(aln, *restricted_mems, cluster_mq, longest_lcp, fraction_filtered, max_mem_length, keep_multimaps, additional_multimaps_for_quality);
    }
    else {
        vector<MaximalExactMatch> mems = find_mems_deep(aln.sequence().begin(),
                                                        aln.sequence().end(),
                                                        longest_lcp,
                                                        fraction_filtered,
                                                        max_mem_length,
                                                        min_mem_length,
                                                        mem_reseed_length,
                                                        false, false, false, 0);
        // query mem hits
        alignments = align_mem_multi(aln, mems, cluster_mq, longest_lcp, fraction_filtered, max_mem_length, keep_multimaps, additional_multimaps_for_quality);
    }

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
    
    annotate_with_mean_path_positions(alignments);

    return alignments;
}

Alignment Mapper::align(const Alignment& aln, int kmer_size, int stride, int max_mem_length, int band_width) {
    // TODO computing mapping quality could be inefficient depending on the method chosen
    
    // Do the multi-mapping
    vector<Alignment> best = align_multi(aln, kmer_size, stride, max_mem_length, band_width);

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


int64_t Mapper::graph_distance(pos_t pos1, pos_t pos2, int64_t maximum) {
    return xg_distance(pos1, pos2, maximum, xindex);
}

int64_t Mapper::graph_mixed_distance_estimate(pos_t pos1, pos_t pos2, int64_t maximum) {
    if (maximum) {
        int64_t graph_dist = graph_distance(pos1, pos2, maximum);
        if (graph_dist < maximum) return graph_dist;
    }
    int64_t path_dist = xindex->min_approx_path_distance(id(pos1), id(pos2));
    int64_t approx_dist = abs(approx_distance(pos1, pos2));
    return min(path_dist, approx_dist);
}

int64_t Mapper::approx_position(pos_t pos) {
    // get nodes on the forward strand
    if (is_rev(pos)) {
        pos = reverse(pos, xg_node_length(id(pos), xindex));
    }
    return (int64_t)xg_node_start(id(pos), xindex) + (int64_t)offset(pos);
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

id_t Mapper::node_approximately_at(int64_t approx_pos) {
    return xindex->node_at_seq_pos(
        min(xindex->seq_length,
            (size_t)max(approx_pos, (int64_t)1)));
}

// use LRU caching to get the most-recent node positions
map<string, vector<size_t> > Mapper::node_positions_in_paths(gcsa::node_type node) {
    return xindex->position_in_paths(gcsa::Node::id(node), gcsa::Node::rc(node), gcsa::Node::offset(node));
}

Alignment Mapper::walk_match(const string& seq, pos_t pos) {
    //cerr << "in walk match with " << seq << " " << seq.size() << " " << pos << endl;
    Alignment aln;
    aln.set_sequence(seq);
    auto alns = walk_match(aln, seq, pos);
    if (!alns.size()) {
        //cerr << "no alignments returned from walk match with " << seq << " " << seq.size() << " " << pos << endl;
        //assert(false);
        return aln;
    }
    aln = alns.front(); // take the first one we found
    //assert(alignment_to_length(aln) == alignment_from_length(aln));
    if (alignment_to_length(aln) != alignment_from_length(aln)
        || alignment_to_length(aln) != seq.size()) {
        //cerr << alignment_to_length(aln) << " is not " << seq.size() << endl;
        //cerr << pb2json(aln) << endl;
        //assert(false);
        aln.clear_path();
    }
#ifdef debug_mapper
    if (debug) {
        cerr << "walk_match result " << pb2json(aln) << endl;
        if (!check_alignment(aln)) {
            cerr << "aln is invalid!" << endl;
            exit(1);
        }
    }
#endif
    return aln;
}

vector<Alignment> Mapper::walk_match(const Alignment& base, const string& seq, pos_t pos) {
    //cerr << "in walk_match " << seq << " from " << pos << " with base " << pb2json(base) << endl;
    // go to the position in the xg index
    // and step in the direction given
    // until we exhaust our sequence
    // or hit another node
    vector<Alignment> alns;
    Alignment aln = base;
    Path& path = *aln.mutable_path();
    Mapping* mapping = path.add_mapping();
    *mapping->mutable_position() = make_position(pos);
#ifdef debug_mapper
#pragma omp critical
    if (debug) cerr << "walking match for seq " << seq << " at position " << pb2json(*mapping) << endl;
#endif
    // get the first node we match
    int total = 0;
    size_t match_len = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        char c = seq[i];
        //cerr << string(base.path().mapping_size(), ' ') << pos << " @ " << i << " on " << c << endl;
        auto nexts = next_pos_chars(pos);
        // we can have a match on the current node
        if (nexts.size() == 1 && id(nexts.begin()->first) == id(pos)) {
            pos_t npos = nexts.begin()->first;
            // check that the next position would match
            if (i+1 < seq.size()) {
                // we can't step, so we break
                //cerr << "Checking if " << pos_char(npos) << " != " << seq[i+1] << endl;
                if (pos_char(npos) != seq[i+1]) {
#ifdef debug_mapper
#pragma omp critical
                    if (debug) cerr << "MEM does not match position, returning without creating alignment" << endl;
#endif
                    return alns;
                }
            }
            // otherwise we step our counters
            ++match_len;
            ++get_offset(pos);
        } else { // or we go into the next node
            // we must be going into another node
            // emit the mapping for this node
            //cerr << "we are going into a new node" << endl;
            // finish the last node
            {
                // we must have matched / we already checked
                ++match_len;
                Edit* edit = mapping->add_edit();
                edit->set_from_length(match_len);
                edit->set_to_length(match_len);
                // reset our counter
                match_len = 0;
            }
            // find the next node that matches our MEM
            bool got_match = false;
            if (i+1 < seq.size()) {
                //cerr << "nexts @ " << i << " " << nexts.size() << endl;
                for (auto& p : nexts) {
                    //cerr << " next : " << p.first << " " << p.second << " (looking for " << seq[i+1] << ")" << endl;
                    if (p.second == seq[i+1]) {
                        if (!got_match) {
                            pos = p.first;
                            got_match = true;
                        } else {
                            auto v = walk_match(aln, seq.substr(i+1), p.first);
                            if (v.size()) {
                                alns.reserve(alns.size() + distance(v.begin(), v.end()));
                                alns.insert(alns.end(), v.begin(), v.end());
                            }
                        }
                    }
                }
                if (!got_match) {
                    // this matching ends here
                    // and we haven't finished matching
                    // thus this path doesn't contain the match
                    //cerr << "got no match" << endl;
                    return alns;
                }

                // set up a new mapping
                mapping = path.add_mapping();
                *mapping->mutable_position() = make_position(pos);
            } else {
                //cerr << "done!" << endl;
            }
        }
    }
    if (match_len) {
        Edit* edit = mapping->add_edit();
        edit->set_from_length(match_len);
        edit->set_to_length(match_len);
    }
    alns.push_back(aln);
#ifdef debug_mapper
#pragma omp critical
    if (debug) {
        cerr << "walked alignment(s):" << endl;
        for (auto& aln : alns) {
            cerr << pb2json(aln) << endl;
        }
    }
#endif
    //cerr << "returning " << alns.size() << endl;
    return alns;
}

// convert one mem into a set of alignments, one for each exact match
vector<Alignment> Mapper::mem_to_alignments(MaximalExactMatch& mem) {
    vector<Alignment> alns;
    const string seq = mem.sequence();
    for (auto& node : mem.nodes) {
        pos_t pos = make_pos_t(node);
        alns.emplace_back(walk_match(seq, pos));
    }
    return alns;
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

Alignment Mapper::patch_alignment(const Alignment& aln, int max_patch_length) {
    //cerr << "top of patch_alignment" << endl;
    Alignment patched;
    // walk along the alignment and find the portions that are unaligned
    int read_pos = 0;
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
                    bands = make_bands(patch, max_patch_length, to_strip);
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
                            pos_t pos_rev = reverse(pos, xg_node_length(id(pos), xindex));
                            Graph graph = xindex->graph_context_id(pos_rev, band.sequence().size()*4);
                            graph.MergeFrom(xindex->graph_context_id(pos, band.sequence().size()*2));
                            sort_by_id_dedup_and_clean(graph);
                            auto proposed_band = align_maybe_flip(band, graph, is_rev(pos), true);
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
                            Graph graph = xindex->graph_context_id(pos, band.sequence().size()*4);
                            pos_t pos_rev = reverse(pos, xg_node_length(id(pos), xindex));
                            graph.MergeFrom(xindex->graph_context_id(pos_rev, band.sequence().size()*2));
                            sort_by_id_dedup_and_clean(graph);
                            //cerr << "on graph " << pb2json(graph) << endl;
                            auto proposed_band = align_maybe_flip(band, graph, is_rev(pos), true);
                            if (proposed_band.score() > max_score) { band = proposed_band; max_score = band.score(); }
                        }
                    } else {
                        //cerr << "internal addition" << endl;
                        int max_score = -std::numeric_limits<int>::max();
                        for (auto& pos : band_ref_pos) {
                            //cerr << "trying position " << pos << endl;
                            Graph graph = xindex->graph_context_id(pos, band.sequence().size()*4);
                            pos_t pos_rev = reverse(pos, xg_node_length(id(pos), xindex));
                            graph.MergeFrom(xindex->graph_context_id(pos_rev, band.sequence().size()*2));
                            sort_by_id_dedup_and_clean(graph);
                            //cerr << "on graph " << pb2json(graph) << endl;
                            auto proposed_band = align_maybe_flip(band, graph, is_rev(pos), true);
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
                            //cerr << "thing worked " << pb2json(band) << endl;
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
                            for (auto& next : positions_bp_from(pos, band.sequence().size(), false)) {
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
                cerr << "done bands" << endl;
                for (auto& band : bands) {
                    cerr << "band: " << pb2json(band) << endl;
                }
                */
                patch = merge_alignments(bands);
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
                //cerr << "adding " << pb2json(patch) << endl;
                if (patched.path().mapping_size()) {
                    extend_alignment(patched, patch);
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
    patched = simplify(patched);
    // set the identity
    patched.set_identity(identity(patched.path()));
    // recompute the score
    patched.set_score(score_alignment(patched, false));
    return patched;
}

// generate a score from the alignment without realigning
// handles split alignments, where gaps of unknown length are
// by estimating length using the positional paths embedded in the graph
int32_t Mapper::score_alignment(const Alignment& aln, bool use_approx_distance) {
    
    // Find the right aligner to score with
    BaseAligner* aligner = get_aligner();
    
    if (use_approx_distance) {
        // Use an approximation
        return aligner->score_gappy_alignment(aln, [&](pos_t last, pos_t next, size_t max_search) {
            return approx_distance(last, next);
        }, strip_bonuses);
    } else {
        // Use the exact method, and if we hit the limit, fall back to the approximate method.
        return aligner->score_gappy_alignment(aln, [&](pos_t last, pos_t next, size_t max_search) {
                return graph_mixed_distance_estimate(last, next, min(32, (int)max_search));
        }, strip_bonuses);
    }
    
}

// make a perfect-match alignment out of a vector of MEMs which each have only one recorded hit
// use the base alignment sequence (which the SMEMs relate to) to fill in the gaps
Alignment Mapper::mems_to_alignment(const Alignment& aln, vector<MaximalExactMatch>& mems) {
    // base case--- empty alignment
    if (mems.empty()) {
        Alignment aln; return aln;
    }
    vector<Alignment> alns;
    // get reference to the start and end of the sequences
    string::const_iterator seq_begin = aln.sequence().begin();
    string::const_iterator seq_end = aln.sequence().end();
    // we use this to track where we need to add sequence
    string::const_iterator last_end = seq_begin;
    for (int i = 0; i < mems.size(); ++i) {
        auto& mem = mems.at(i);
        //cerr << "looking at " << mem.sequence() << endl;
        // this mem is contained in the last
        if (mem.end <= last_end) {
            continue;
        }
        // handle unaligned portion between here and the last SMEM or start of read
        if (mem.begin > last_end) {
            alns.emplace_back();
            alns.back().set_sequence(aln.sequence().substr(last_end - seq_begin, mem.begin - last_end));
        }
        Alignment aln = mem_to_alignment(mem);
        // find and trim overlap with previous
        if (i > 0) {
            // use the end of the last mem we touched (we may have skipped several)
            int overlap = last_end - mem.begin;
            if (overlap > 0) {
                aln = strip_from_start(aln, overlap);
            }
        }
        alns.push_back(aln);
        last_end = mem.end;
    }
    // handle unaligned portion at end of read
    int start = last_end - seq_begin;
    int length = seq_end - (seq_begin + start);
    
    alns.emplace_back();
    alns.back().set_sequence(aln.sequence().substr(start, length));

    auto alnm = merge_alignments(alns);
    *alnm.mutable_quality() = aln.quality();
    return alnm;
}

// convert one mem into an alignment; validates that only one node is given
Alignment Mapper::mem_to_alignment(MaximalExactMatch& mem) {
    const string seq = mem.sequence();
    if (mem.nodes.size() > 1) {
        cerr << "[vg::Mapper] warning: generating first alignment from MEM with multiple recorded hits" << endl;
    }
    auto& node = mem.nodes.front();
    pos_t pos = make_pos_t(node);
    return walk_match(seq, pos);
}

// transform the path into a path relative to another path (defined by path_name)
// source -> surjection (in path_name coordinate space)
// the product is equivalent to a pairwise alignment between this path and the other

// new approach
// get path sequence
// get graph component overlapping path
// removing elements which aren't in the path of interest
// realign to this graph
// cross fingers

Alignment Mapper::surject_alignment(const Alignment& source,
                                    set<string>& path_names,
                                    string& path_name,
                                    int64_t& path_pos,
                                    bool& path_reverse,
                                    int window) {

    Alignment surjection = source;
    // Leave the original mapping quality in place (because that's the quality
    // on the placement of this read in this region at all)
    surjection.clear_score();
    surjection.clear_identity();
    surjection.clear_path();

    // get start and end nodes in path
    // get range between +/- window
    if (!source.has_path() || source.path().mapping_size() == 0) {
#ifdef debug_mapper

#pragma omp critical (cerr)
        cerr << "Alignment " << source.name() << " is unmapped and cannot be surjected" << endl;

#endif
        return surjection;
    }

    set<id_t> nodes;
    for (int i = 0; i < source.path().mapping_size(); ++ i) {
        nodes.insert(source.path().mapping(i).position().node_id());
    }
    VG graph;
    for (auto& node : nodes) {
        *graph.graph.add_node() = xindex->node(node);
    }
    xindex->expand_context(graph.graph, context_depth, true); // get connected edges and path
    graph.paths.append(graph.graph);
    graph.rebuild_indexes();

    set<string> kept_paths;
    graph.keep_paths(path_names, kept_paths);

    // We need this for inverting mappings to the correct strand
    function<int64_t(id_t)> node_length = [&graph](id_t node) {
        return graph.get_node(node)->sequence().size();
    };
    
    // What is our alignment to surject spelled the other way around? We can't
    // just use the normal alignment RC function because the mappings reference
    // nonexistent nodes.
    // Make sure to copy all the things about the alignment (name, etc.)

    Alignment surjection_rc = surjection;
    surjection_rc.set_sequence(reverse_complement(surjection.sequence()));
    
    // Align the old alignment to the graph in both orientations. Apparently
    // align only does a single oriantation, and we have no idea, even looking
    // at the mappings, which of the orientations will correspond to the one the
    // alignment is actually in.

    Graph subgraph = graph.graph;
    sort_by_id_dedup_and_clean(subgraph);
    auto surjection_forward = align_to_graph(surjection, subgraph, max_query_graph_ratio, true);
    auto surjection_reverse = align_to_graph(surjection_rc, subgraph, max_query_graph_ratio, true);

#ifdef debug_mapper
#pragma omp critical (cerr)
    cerr << surjection.name() << " " << surjection_forward.score() << " forward score, " << surjection_reverse.score() << " reverse score" << endl;
#endif
    
    if(surjection_reverse.score() > surjection_forward.score()) {
        // Even if we have to surject backwards, we have to send the same string out as we got in.
        surjection = reverse_complement_alignment(surjection_reverse, node_length);
    } else {
        surjection = surjection_forward;
    }
    
    
#ifdef debug_mapper

#pragma omp critical (cerr)
        cerr << surjection.path().mapping_size() << " mappings, " << kept_paths.size() << " paths" << endl;

#endif

    if (surjection.path().mapping_size() > 0 && kept_paths.size() == 1) {
        // determine the paths of the node we mapped into
        //  ... get the id of the first node, get the paths of it
        assert(kept_paths.size() == 1);
        path_name = *kept_paths.begin();

        int64_t path_id = xindex->path_rank(path_name);
        auto& first_pos = surjection.path().mapping(0).position();
        int64_t hit_id = surjection.path().mapping(0).position().node_id();
        bool hit_backward = surjection.path().mapping(0).position().is_reverse();
        // we pick up positional information using the index

        auto path_posns = xindex->position_in_path(hit_id, path_name);
        if (path_posns.size() > 1) {
            cerr << "[vg map] surject_alignment: warning, multiple positions for node " << hit_id << " in " << path_name << " but will use only first: " << path_posns.front() << endl;
        } else if (path_posns.size() == 0) {
            cerr << "[vg map] surject_alignment: error, no positions for alignment " << source.name() << endl;
            exit(1);
        }

        // if we are reversed
        path_pos = path_posns.front();
        bool reversed_path = xindex->mapping_at_path_position(path_name, path_pos).position().is_reverse();
        if (reversed_path) {
            // if we got the start of the node position relative to the path
            // we need to offset to make thinsg right
            // but which direction
            if (hit_backward) {
                path_pos = path_posns.front() + first_pos.offset();
            } else {
                auto pos = reverse_complement_alignment(surjection, node_length).path().mapping(0).position();
                path_pos = xindex->position_in_path(pos.node_id(), path_name).front() + pos.offset();
            }
            path_reverse = !hit_backward;
        } else {
            if (!hit_backward) {
                path_pos = path_posns.front() + first_pos.offset();
            } else {
                auto pos = reverse_complement_alignment(surjection, node_length).path().mapping(0).position();
                path_pos = xindex->position_in_path(pos.node_id(), path_name).front() + pos.offset();
            }
            path_reverse = hit_backward;
        }

    } else {

        surjection = source;
#ifdef debug_mapper

#pragma omp critical (cerr)
        cerr << "Alignment " << source.name() << " did not align to the surjection subgraph" << endl;

#endif

    }

#ifdef debug_mapper
    
#pragma omp critical (cerr)
    cerr << "Surjection on reverse strand? " << path_reverse << endl;
    cerr << "Surjected alignment: " << pb2json(surjection) << endl;
    
#endif
    
    return surjection;
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
    const function<double(const Alignment&, const Alignment&, const map<string, double>&, const map<string, double>&)>& transition_weight,
    int vertex_band_width,
    int position_depth,
    int max_connections) {
    // store the Alignments in the model
    int offset = 0;
    int idx = 0;
    for (auto& band : bands) {
        for (auto& aln : band) {
            AlignmentChainModelVertex v;
            v.aln = &aln;
            v.band_begin = offset;
            v.band_idx = idx;
            v.weight = aln.sequence().size() + aln.score() + aln.mapping_quality();
            v.prev = nullptr;
            v.score = 0;
            v.positions = mapper->alignment_mean_path_positions(aln);
            v.positions[""] = mapper->approx_alignment_position(aln);
            model.push_back(v);
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
            if (v->next_cost.size() < max_connections && u->prev_cost.size() < max_connections) {
                if (v->band_idx + vertex_band_width <= u->band_idx) {
                    double weight = transition_weight(*v->aln, *u->aln, v->positions, u->positions);
                    if (weight > -std::numeric_limits<double>::max()) {
                        v->next_cost.push_back(make_pair(&*u, weight));
                        u->prev_cost.push_back(make_pair(&*v, weight));
                    }
                } else if (u->band_idx + vertex_band_width <= v->band_idx) {
                    double weight = transition_weight(*u->aln, *v->aln, u->positions, v->positions);
                    if (weight > -std::numeric_limits<double>::max()) {
                        u->next_cost.push_back(make_pair(&*v, weight));
                        v->prev_cost.push_back(make_pair(&*u, weight));
                    }
                }
            }            
        }
    }
}

void AlignmentChainModel::score(const set<AlignmentChainModelVertex*>& exclude) {
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
    set<AlignmentChainModelVertex*> exclude;
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
        set<AlignmentChainModelVertex*> chain_members;
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
    }
    vector<Alignment> alns;
    for (auto& trace : traces) {
        alns.emplace_back();
        Alignment& merged = alns.back();
        merged = merge_alignments(trace);
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
    double a = normal_inverse_cdf(1.0 - 0.5 * (1.0 - robust_estimation_fraction));
    sigma = sqrt(raw_var / (1.0 - 2.0 * a * normal_pdf(a, 0.0, 1.0)));
}
    
double FragmentLengthDistribution::mean() const {
    return mu;
}

double FragmentLengthDistribution::stdev() const {
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
